// #include <cilk/cilk.h>
// #include <cilk/cilk_api.h>

#include "globals-new.hpp"
#include "grammar.hpp"

/***
 * V-centric (forward, temporal ptrs, grammar driven) - the real one
 * Bi directional traversing of the Graph (incoming and outgoing edges)
 * One buffer list instead of three for OLD, NEW, and FUTURE edges
 * Buffer struct (adjcency list is made with this struct) holds the pointes for the OLD, NEW, and FUTURE edges in the single buffer list
 * Check the Future edge flag whenever an edge is created
 * Grammar index
 */

void getPeakMemoryUsage()
{   
    std::string line;
    std::ifstream statusFile("/proc/self/status");
    
    while (getline(statusFile, line))
    {   
        if (line.substr(0, 7) == "VmPeak:")
        {   
            long memoryKb;
            std::istringstream iss(line.substr(7));
            iss >> memoryKb;
            double memGB = memoryKb / (1024.0 * 1024.0);
            std::cout << "Peak Virtual Memory Usage: " << memoryKb << " KB" << std::endl;
            std::cout << "Peak Virtual Memory Usage (in GB): " << memGB << " GB" << std::endl;
        }
        
        if (line.substr(0, 6) == "VmHWM:")
        {   
            long memoryKb;
            std::istringstream iss(line.substr(6));
            iss >> memoryKb;
            double memGB = memoryKb / (1024.0 * 1024.0); 
            std::cout << "Peak Physical  Memory Usage: " << memoryKb << " KB" << std::endl;
            std::cout << "Peak Physical Memory Usage (in GB): " << memGB << " GB" << std::endl;
        }
        
        if (line.substr(0, 6) == "VmRSS:")
        {   
            long memoryKb;
            std::istringstream iss(line.substr(6));
            iss >> memoryKb;
            double memGB = memoryKb / (1024.0 * 1024.0);
            std::cout << "VmRSS  Memory Usage: " << memoryKb << " KB" << std::endl;
            std::cout << "VmRSS Memory Usage (in GB): " << memGB << " GB" << std::endl;
        }
    }
}


int main(int argc, char **argv)
{

    bool debug = false;
    // Get the graph file path and grammar file path from command line argument
    if (argc == 1)
    {
        cout << "Please provide the graph file path and grammar file path. For example, ./topo-driven <graph_file> <grammar_file>" << endl;
        return 0;
    }
    else if (argc == 2)
    {
        cout << "Please provide the grammar file path. For example, ./topo-driven <graph_file> <grammar_file>" << endl;
        return 0;
    }

    // read Graph and Grammar
    const string inputGraph = argv[1];
    Grammar grammar(argv[2]);

    //        grammar.printGrammarIndex3_2();
    //      exit(1);

    uint num_nodes = 0;
    uint num_edges = 0;

    ifstream infile(inputGraph);

    vector<EdgeForReading> edges;
    unordered_set<uint> nodes;
    EdgeForReading newEdge;
    uint from, to;
    string label;
    while (infile >> newEdge.from)
    {
        infile >> newEdge.to;
        infile >> label;
        // if the label does not exist in the Grammar, discard
        if (grammar.hashSym.find(label) == grammar.hashSym.end())
        {
            continue;
        }
        newEdge.label = grammar.hashSym[label];

        edges.push_back(newEdge);

        num_edges++;
        num_nodes = max(num_nodes, max(newEdge.from + 1, newEdge.to + 1));
        // If we insert the nodes into the set, the duplicate nodes will automatically get
        // discarded by the set.
        nodes.insert(newEdge.from);
        nodes.insert(newEdge.to);
    }

    infile.close();

    // level-1: vertex ID, level 2:  grammar labels level-3: NEW, OLD, FUTURE pointers and incoming edges
    // vector<vector<Buffer>> inEdgeVecs(num_nodes, vector<Buffer>(grammar.labelSize));
    // level-1: vertex ID, level-2: NEW, OLD, FUTURE pointers and outgoing edges
    vector<vector<Buffer>> inEdgeVecs(grammar.labelSize, vector<Buffer>(num_nodes));

    vector<vector<unordered_set<ull>>> inHashset(num_nodes, vector<unordered_set<ull>>(grammar.labelSize, unordered_set<ull>()));
    // hashset for incoming edges
    // unordered_set<ull> *inHashset = new unordered_set<ull>[num_nodes];

    cout << "#nodes " << num_nodes << endl;
    cout << "SF::#nodes " << nodes.size() << endl;
    cout << "#edges " << num_edges << endl;

    cout << "Start making sets and the hash ...\n";

    for (uint i = 0; i < num_edges; i++)
    {
        inEdgeVecs[edges[i].label][edges[i].to].vertexList.push_back(edges[i].from);
        // inEdgeVecs[edges[i].to][edges[i].label].vertexList.push_back(edges[i].from);

        // update the buffer pointers
        inEdgeVecs[edges[i].label][edges[i].to].NEW_END++;
        // inEdgeVecs[edges[i].to][edges[i].label].NEW_END++;

        // insert edge into the hashset
        inHashset[edges[i].to][edges[i].label].insert(edges[i].from);
    }

    edges.clear();

    // count the total no of initial unique edge
    uint initialEdgeCount = countEdge(inHashset, num_nodes, grammar.labelSize);

    cout << "Done!\n";

    // currently exclude the initialization time
    std::chrono::time_point<std::chrono::steady_clock> start, finish;
    std::chrono::time_point<std::chrono::steady_clock> startTransfer, finishTransfer;
    std::chrono::time_point<std::chrono::steady_clock> startItr, finishItr;
    start = std::chrono::steady_clock::now();

    uint newEdgeCnt = initialEdgeCount;
    uint oldEdgeCnt = 0;
    uint futureEdgeCnt = 0;

    atomic<int> newEdgeCounter(0);
    bool finished;
    int itr = 0;
    ull calcCnt = 0;
    double elapsed_seconds_comp = 0.0;
    double elapsed_seconds_transfer = 0.0;

    // handle epsilon rules: add an edge to itself
    // grammar1 is for epsilon rules A --> e
    // grammar2 is for one symbol on RHS A --> B
    // grammar3 is for two symbols on RHS A --> BC
    for (uint l = 0; l < grammar.grammar1.size(); l++)
    {
        for (uint i = 0; i < num_nodes; i++)
        {
            calcCnt++;
            // check if the new edge based on an epsilon grammar rule exists or not. l: grammar ID, 0: LHS
            if (inHashset[i][grammar.grammar1[l][0]].find(i) == inHashset[i][grammar.grammar1[l][0]].end())
            {
                // insert edge into hashset
                inHashset[i][grammar.grammar1[l][0]].insert(i);
                // insert edge into the graph
                inEdgeVecs[grammar.grammar1[l][0]][i].vertexList.push_back(i);
                // inEdgeVecs[i][grammar.grammar1[l][0]].vertexList.push_back(i);

                // update the sliding/temporal pointers
                inEdgeVecs[grammar.grammar1[l][0]][i].NEW_END++;
                // inEdgeVecs[i][grammar.grammar1[l][0]].NEW_END++;

                // newEdgeCounter++;
		        newEdgeCnt++;
            }
        }
    }

    cout << "********************\n";

    do
    {
    	//startItr = std::chrono::steady_clock::now();

        finished = true;
        itr++;

        // for each grammar rule like A --> B
        for (uint g = 0; g < grammar.labelSize; g++)
        {
            for (int i = num_nodes - 1; i >= 0; i--)
            {
                uint inNbr1;
                //  the valid index range is [START_NEW, END_NEW-1]
                uint START_NEW = inEdgeVecs[g][i].OLD_END;
                uint END_NEW = inEdgeVecs[g][i].NEW_END;

                // for each new edge
                for (uint j = START_NEW; j < END_NEW; j++)
                {
                    inNbr1 = inEdgeVecs[g][i].vertexList[j];
                    // rule: A = B
                    for (uint m = 0; m < grammar.grammar2index[g].size(); m++)
                    {
                        // calcCnt++;
                        uint A = grammar.grammar2index[g][m];

                        if (inHashset[i][A].find(inNbr1) == inHashset[i][A].end())
                        {
                            finished = false;
                            inHashset[i][A].insert(inNbr1);
                            inEdgeVecs[A][i].vertexList.push_back(inNbr1);

			    futureEdgeCnt++;
                        }
                    }

                    // rule: A = BC
                    // In this case the C is read in first hop, then B.
                    // all the OLD, and NEW incoming edges of the first edge
                    for (uint m = 0; m < grammar.grammar3indexRight[g].size(); m++)
                    {
                        // C = g
                        uint B = grammar.grammar3indexRight[g][m].first;
                        uint A = grammar.grammar3indexRight[g][m].second;

                        uint START_OLD_OUT = 0;
                        uint END_NEW_OUT = inEdgeVecs[B][inNbr1].NEW_END;

                        for (uint h = START_OLD_OUT; h < END_NEW_OUT; h++)
                        {
                            uint inNbr2 = inEdgeVecs[B][inNbr1].vertexList[h];

                            // calcCnt++;
                            if (inHashset[i][A].find(inNbr2) == inHashset[i][A].end())
                            {
                                finished = false;
                                inHashset[i][A].insert(inNbr2);
                                inEdgeVecs[A][i].vertexList.push_back(inNbr2);

				futureEdgeCnt++;
                            }
                        }
                    }
                }

                uint START_OLD = 0;
                uint END_OLD = inEdgeVecs[g][i].OLD_END;

                // for each old edge
                for (uint j = START_OLD; j < END_OLD; j++)
                {
                    inNbr1 = inEdgeVecs[g][i].vertexList[j];

                    // rule: A = BC
                    // all the NEW incoming edges of the first edge
                    for (uint m = 0; m < grammar.grammar3indexRight[g].size(); m++)
                    {
                        // C = g
                        uint B = grammar.grammar3indexRight[g][m].first;
                        uint A = grammar.grammar3indexRight[g][m].second;

                        uint START_NEW_OUT = inEdgeVecs[B][inNbr1].OLD_END;
                        uint END_NEW_OUT = inEdgeVecs[B][inNbr1].NEW_END;

                        for (uint h = START_NEW_OUT; h < END_NEW_OUT; h++)
                        {
                            uint inNbr2 = inEdgeVecs[B][inNbr1].vertexList[h];

                            // calcCnt++;
                            if (inHashset[i][A].find(inNbr2) == inHashset[i][A].end())
                            {
                                finished = false;
                                inHashset[i][A].insert(inNbr2);
                                inEdgeVecs[A][i].vertexList.push_back(inNbr2);

				futureEdgeCnt++;
                            }
                        }
                    }
                }
            }
        }

    	//finishItr = std::chrono::steady_clock::now();

        cout << "---------------------------------" << itr << endl;
        cout << "Iteration number " << itr << endl;
        cout << "---------------------------------" << itr << endl;


        startTransfer = std::chrono::steady_clock::now();


        for (uint g = 0; g < grammar.labelSize; g++)
        {
            // update the sliding/temporal pointers
            for (int i = 0; i < num_nodes; i++)
            {
                inEdgeVecs[g][i].OLD_END = inEdgeVecs[g][i].NEW_END;
                inEdgeVecs[g][i].NEW_END = inEdgeVecs[g][i].vertexList.size();
            }
        }

        //finishTransfer = std::chrono::steady_clock::now();
       	//std::chrono::duration<double> elapsedSecondsTransferCurr = finishTransfer - startTransfer;
        //elapsed_seconds_transfer += elapsedSecondsTransferCurr.count();

       	//std::chrono::duration<double> elapsed_secondsItr = finishItr - startItr;
       	//double elapsed_seconds_itr = elapsed_secondsItr.count();

	    //cout <<"TIME:\t" << elapsed_seconds_itr << endl;
	    //cout <<"FUTURE EDGES:\t" << futureEdgeCnt << endl; 
	    //cout << "NEW EDGE:\t" << newEdgeCnt << endl;
        //cout << "OLD EDGE:\t" << oldEdgeCnt << endl;

	    oldEdgeCnt += newEdgeCnt;
        newEdgeCnt = futureEdgeCnt;
        futureEdgeCnt = 0;

    } while (!finished);

    finish = std::chrono::steady_clock::now();
    std::chrono::duration<double> elapsed_seconds = finish - start;
    // std::time_t finish_time = std::chrono::steady_clock::to_time_t(finish);

    uint totalNewEdgeCount = countEdge(inHashset, num_nodes, grammar.labelSize) - initialEdgeCount;

    cout << "**************************" << endl;
    std::cout << "# Total time = " << elapsed_seconds.count() << std::endl;
    std::cout << "# Total computation time = " << elapsed_seconds_comp << std::endl;
    std::cout << "# Total data transfer time = " << elapsed_seconds_transfer << std::endl;

    cout << "SF:: # Number of new edges: " << totalNewEdgeCount << endl;
    cout << "AM:: # Number of new edges: " << newEdgeCounter << endl;

    cout << "# Iterations: " << itr << endl;
    cout << "# Total Calculations: " << calcCnt << endl;

    getPeakMemoryUsage();
}
