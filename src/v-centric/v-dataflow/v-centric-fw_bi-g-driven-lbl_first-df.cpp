#include "../globals-new.hpp"
#include "../grammar.hpp"

/***
 * V-centric (bi, temporal ptrs, grammar driven) - the real one
 * adjacency list: first level: grammar label, second level: vertex
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
    }
}

bool calcFirst(Grammar &grammar, vector<vector<Buffer>> &edgeVecs, vector<vector<Buffer>> &inEdgeVecs, vector<vector<unordered_set<ull>>> &hashset, uint num_nodes, ull &calcCnt)
{
    bool finished = true;

    // for each grammar rule like A --> B
    for (uint g = 0; g < grammar.labelSize; g++)
    {
        for (uint i = 0; i < num_nodes; i++)
        {
            uint nbr;
            //  the valid index range is [START_NEW, END_NEW-1]
            uint START_NEW = edgeVecs[g][i].OLD_END;
            uint END_NEW = edgeVecs[g][i].NEW_END;

            // for each new edge
            for (uint j = START_NEW; j < END_NEW; j++)
            {
                nbr = edgeVecs[g][i].vertexList[j];
                // rule: A = BC
                // new * new
                uint g3SizeLeft = grammar.grammar3indexLeft[g].size();
                for (uint m = 0; m < g3SizeLeft; m++)
                {
                    uint C = grammar.grammar3indexLeft[g][m].first;
                    uint A = grammar.grammar3indexLeft[g][m].second;

                    uint START_NEW_OUT = edgeVecs[C][nbr].OLD_END;
                    uint END_NEW_OUT = edgeVecs[C][nbr].NEW_END;

                    for (uint h = START_NEW_OUT; h < END_NEW_OUT; h++)
                    {
                        uint outNbr = edgeVecs[C][nbr].vertexList[h];

                        calcCnt++;
                        if (hashset[i][A].find(outNbr) == hashset[i][A].end())
                        {
                            finished = false;
                            hashset[i][A].insert(outNbr);

                            edgeVecs[A][i].vertexList.push_back(outNbr);
                            inEdgeVecs[A][outNbr].vertexList.push_back(i);
                        }
                    }
                }
            }
        }
    }

    for (uint g = 0; g < grammar.labelSize; g++)
    {
        // update the sliding/temporal pointers
        for (int i = 0; i < num_nodes; i++)
        {
            edgeVecs[g][i].OLD_END = edgeVecs[g][i].NEW_END;
            inEdgeVecs[g][i].OLD_END = inEdgeVecs[g][i].NEW_END;

            edgeVecs[g][i].NEW_END = edgeVecs[g][i].vertexList.size();
            inEdgeVecs[g][i].NEW_END = inEdgeVecs[g][i].vertexList.size();
        }
    }

    return finished;
}

bool calcDf(Grammar &grammar, vector<vector<Buffer>> &edgeVecs, vector<vector<Buffer>> &inEdgeVecs, vector<vector<unordered_set<ull>>> &hashset, uint num_nodes, ull &calcCnt)
{
    bool finished = true;

    // for each grammar rule like A --> B
    for (uint g = 0; g < grammar.labelSize; g++)
    {
        for (uint i = 0; i < num_nodes; i++)
        {
            uint nbr;
            //  the valid index range is [START_NEW, END_NEW-1]
            uint START_NEW = edgeVecs[g][i].OLD_END;
            uint END_NEW = edgeVecs[g][i].NEW_END;

            // for each new edge
            for (uint j = START_NEW; j < END_NEW; j++)
            {
                nbr = edgeVecs[g][i].vertexList[j];

                // rule: A = BC
                // new * old
                uint g3SizeLeft = grammar.grammar3indexLeft[g].size();
                for (uint m = 0; m < g3SizeLeft; m++)
                {
                    uint C = grammar.grammar3indexLeft[g][m].first;
                    uint A = grammar.grammar3indexLeft[g][m].second;

                    uint START_OLD_OUT = 0;
                    uint END_OLD_OUT = edgeVecs[C][nbr].OLD_END;

                    for (uint h = START_OLD_OUT; h < END_OLD_OUT; h++)
                    {
                        uint outNbr = edgeVecs[C][nbr].vertexList[h];

                        calcCnt++;
                        if (hashset[i][A].find(outNbr) == hashset[i][A].end())
                        {
                            finished = false;
                            hashset[i][A].insert(outNbr);

                            edgeVecs[A][i].vertexList.push_back(outNbr);
                            inEdgeVecs[A][outNbr].vertexList.push_back(i);
                        }
                    }
                }
            }
        }
    }

    for (uint g = 0; g < grammar.labelSize; g++)
    {
        // update the sliding/temporal pointers
        for (int i = 0; i < num_nodes; i++)
        {
            edgeVecs[g][i].OLD_END = edgeVecs[g][i].NEW_END;
            inEdgeVecs[g][i].OLD_END = inEdgeVecs[g][i].NEW_END;

            edgeVecs[g][i].NEW_END = edgeVecs[g][i].vertexList.size();
            inEdgeVecs[g][i].NEW_END = inEdgeVecs[g][i].vertexList.size();
        }
    }

    return finished;
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
    vector<vector<Buffer>> inEdgeVecs(grammar.labelSize, vector<Buffer>(num_nodes));
    // level-1: vertex ID, level-2: NEW, OLD, FUTURE pointers and outgoing edges
    vector<vector<Buffer>> edgeVecs(grammar.labelSize, vector<Buffer>(num_nodes));

    vector<vector<unordered_set<ull>>> hashset(num_nodes, vector<unordered_set<ull>>(grammar.labelSize, unordered_set<ull>()));
    // hashset for incoming edges
    // unordered_set<ull> *inHashset = new unordered_set<ull>[num_nodes];

    cout << "#nodes " << num_nodes << endl;
    cout << "SF::#nodes " << nodes.size() << endl;
    cout << "#edges " << num_edges << endl;

    cout << "Start making sets and the hash ...\n";

    for (uint i = 0; i < num_edges; i++)
    {
        edgeVecs[edges[i].label][edges[i].from].vertexList.push_back(edges[i].to);
        inEdgeVecs[edges[i].label][edges[i].to].vertexList.push_back(edges[i].from);

        // update the buffer pointers
        edgeVecs[edges[i].label][edges[i].from].NEW_END++;
        inEdgeVecs[edges[i].label][edges[i].to].NEW_END++;

        // insert edge into the hashset
        hashset[edges[i].from][edges[i].label].insert(edges[i].to);
    }

    edges.clear();

    // count the total no of initial unique edge
    uint initialEdgeCount = countEdge(hashset, num_nodes, grammar.labelSize);

    cout << "Done!\n";

    // currently exclude the initialization time
    std::chrono::time_point<std::chrono::steady_clock> start, finish;
    start = std::chrono::steady_clock::now();

    atomic<int> newEdgeCounter(0);
    bool finished;
    int itr = 0;
    ull calcCnt = 0;
    double elapsed_seconds_comp = 0.0;
    double elapsed_seconds_transfer = 0.0;

    std::chrono::time_point<std::chrono::steady_clock> start_comp, finish_comp, start_transfer, end_transfer;

    std::cout << "********************\n";

    finished = calcFirst(grammar, edgeVecs, inEdgeVecs, hashset, num_nodes, calcCnt);
    std::cout << "Iteration number " << itr << endl;
    itr++;

    while (!finished)
    {
        itr++;

        // start_comp = std::chrono::steady_clock::now();
        finished = calcDf(grammar, edgeVecs, inEdgeVecs, hashset, num_nodes, calcCnt);

        std::cout << "Iteration number " << itr << endl;

        // start_transfer = std::chrono::steady_clock::now();
        // std::chrono::duration<double> interval_tranfer = (std::chrono::steady_clock::now() - start_transfer);
        // elapsed_seconds_transfer += interval_tranfer.count();
    }

    finish = std::chrono::steady_clock::now();
    std::chrono::duration<double> elapsed_seconds = finish - start;
    // std::time_t finish_time = std::chrono::steady_clock::to_time_t(finish);

    uint totalNewEdgeCount = countEdge(hashset, num_nodes, grammar.labelSize) - initialEdgeCount;

    std::cout << "**************************" << endl;
    std::cout << "# Total time = " << elapsed_seconds.count() << std::endl;
    std::cout << "# Total computation time = " << elapsed_seconds_comp << std::endl;
    std::cout << "# Total data transfer time = " << elapsed_seconds_transfer << std::endl;
    std::cout << "# Total sum of comp + transfer = " << elapsed_seconds_comp + elapsed_seconds_transfer << std::endl;

    std::cout << "SF:: # Number of new edges: " << totalNewEdgeCount << endl;
    std::cout << "AM:: # Number of new edges: " << newEdgeCounter << endl;

    std::cout << "# Iterations: " << itr << endl;
    std::cout << "# Total Calculations: " << calcCnt << endl;

    getPeakMemoryUsage();
}
