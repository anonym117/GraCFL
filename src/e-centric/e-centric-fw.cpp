#include "globals.hpp"
#include "grammar.hpp"

/***
 * E-centric (forward, sliding ptrs, grammar driven) 
 * Grammar driven
 * One buffer list instead of three for OLD, NEW, and FUTURE edges
 * Buffer struct (adjcency list is made with this struct) holds the pointes for the OLD, NEW, and FUTURE edges in the single buffer list
 * Check the Future edge flag whenever an edge is created
 */

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
    vector<vector<Buffer>> edgeVecs(grammar.labelSize, vector<Buffer>(num_nodes));

    vector<vector<unordered_set<ull>>> hashset(num_nodes, vector<unordered_set<ull>>(grammar.labelSize, unordered_set<ull>()));

    queue<EdgeForReading> newWorklist;
    queue<EdgeForReading> futureNewWorklist;
    vector<EdgeForReading> oldWorklist;
    // queue<EdgeForReading> futureOldWorklist;
    // hashset for incoming edges
    // unordered_set<ull> *inHashset = new unordered_set<ull>[num_nodes];

    std::cout << "#nodes " << num_nodes << endl;
    std::cout << "SF::#nodes " << nodes.size() << endl;
    std::cout << "#edges " << num_edges << endl;

    std::cout << "Start making sets and the hash ...\n";

    for (uint i = 0; i < num_edges; i++)
    {
        edgeVecs[edges[i].label][edges[i].from].vertexList.push_back(edges[i].to);
        // inEdgeVecs[edges[i].to][edges[i].label].vertexList.push_back(edges[i].from);

        // update the buffer pointers
        edgeVecs[edges[i].label][edges[i].from].NEW_END++;
        // inEdgeVecs[edges[i].to][edges[i].label].NEW_END++;

        // insert edge into the hashset
        hashset[edges[i].from][edges[i].label].insert(edges[i].to);

        newWorklist.push(edges[i]);
    }

    edges.clear();

    // count the total no of initial unique edge
    uint initialEdgeCount = countEdge(hashset, num_nodes, grammar.labelSize);

    cout << "Done!\n";

    // currently exclude the initialization time
    std::chrono::time_point<std::chrono::steady_clock> start, finish;
    std::chrono::time_point<std::chrono::steady_clock> startTransfer, finishTransfer;
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
            if (hashset[i][grammar.grammar1[l][0]].find(i) == hashset[i][grammar.grammar1[l][0]].end())
            {
                // insert edge into hashset
                hashset[i][grammar.grammar1[l][0]].insert(i);
                // insert edge into the graph
                edgeVecs[grammar.grammar1[l][0]][i].vertexList.push_back(i);
                // inEdgeVecs[i][grammar.grammar1[l][0]].vertexList.push_back(i);

                // update the sliding/temporal pointers
                edgeVecs[grammar.grammar1[l][0]][i].NEW_END++;
                // inEdgeVecs[i][grammar.grammar1[l][0]].NEW_END++;

                newWorklist.push(EdgeForReading(i, i, grammar.grammar1[l][0]));

                // newEdgeCounter++;
                newEdgeCnt++;
            }
        }
    }

    cout << "********************\n";

    do
    {

        finished = true;
        itr++;

        uint CURR_OLD_SIZE = oldWorklist.size();
        std::cout << "CURR_OLD_SIZE: " << CURR_OLD_SIZE << endl;

        while (!newWorklist.empty())
        {
            EdgeForReading currEdge = newWorklist.front();
            newWorklist.pop();

            // for each grammar rule, A --> B
            for (uint g = 0; g < grammar.grammar2index[currEdge.label].size(); g++)
            {
                calcCnt++;
                // calcCnt2++;
                uint leftLabel = grammar.grammar2index[currEdge.label][g];
                // if the edge to the neighbor is labeled with B
                if (hashset[currEdge.from][leftLabel].find(currEdge.to) == hashset[currEdge.from][leftLabel].end())
                {
                    // finished = false;
                    hashset[currEdge.from][leftLabel].insert(currEdge.to);

                    // the new edge is added to FUTURE list, due to potential conflicts in parallel version
                    // if topo-driven is run in serial, only OLD and NEW are needed.
                    // Graspan did not mention FUTURE list, but it puts the new edges into a seperate list, which
                    // is similar to FUTURE list's purpose.
                    // for data-driven, we only need one list (no NEW OLD FUTURE)
                    edgeVecs[leftLabel][currEdge.from].vertexList.push_back(currEdge.to);
                    futureNewWorklist.push(EdgeForReading(currEdge.from, currEdge.to, leftLabel));

                    // newEdgeCounter++;
                }
            }

            // A = BC
            for (uint g = 0; g < grammar.grammar3indexLeft[currEdge.label].size(); g++)
            {
                uint B = currEdge.label;
                uint A = grammar.grammar3indexLeft[currEdge.label][g].second;
                uint C = grammar.grammar3indexLeft[currEdge.label][g].first;

                uint OLD_START = 0;
                uint NEW_END = edgeVecs[C][currEdge.to].NEW_END;

                for (uint j = OLD_START; j < NEW_END; j++)
                {
                    uint nbr;

                    // A = BC, B = grammar.grammar3[g][1] = currEdge.label
                    // NEW outgoing edges of currEdge.to
                    nbr = edgeVecs[C][currEdge.to].vertexList[j];
                    // A = BC, C = grammar.grammar3[g][2]
                    calcCnt++;
                    // calcCnt3_left++;

                    if (hashset[currEdge.from][A].find(nbr) == hashset[currEdge.from][A].end())
                    {
                        // finished = false;
                        hashset[currEdge.from][A].insert(nbr);
                        futureNewWorklist.push(EdgeForReading(currEdge.from, nbr, A));
                        edgeVecs[A][currEdge.from].vertexList.push_back(nbr);
                        // newEdgeCounter++;
                    }
                }
            }

            oldWorklist.push_back(currEdge);
        }

        for (uint i = 0; i < CURR_OLD_SIZE; i++)
        {
            EdgeForReading currEdge = oldWorklist[i];
            // A = BC
            for (uint g = 0; g < grammar.grammar3indexLeft[currEdge.label].size(); g++)
            {
                uint B = currEdge.label;
                uint A = grammar.grammar3indexLeft[currEdge.label][g].second;
                uint C = grammar.grammar3indexLeft[currEdge.label][g].first;

                uint NEW_START = edgeVecs[C][currEdge.to].OLD_END;
                uint NEW_END = edgeVecs[C][currEdge.to].NEW_END;

                for (uint j = NEW_START; j < NEW_END; j++)
                {
                    uint nbr;

                    // A = BC, B = grammar.grammar3[g][1] = currEdge.label
                    // NEW outgoing edges of currEdge.to
                    nbr = edgeVecs[C][currEdge.to].vertexList[j];
                    // A = BC, C = grammar.grammar3[g][2]
                    calcCnt++;
                    // calcCnt3_left++;

                    if (hashset[currEdge.from][A].find(nbr) == hashset[currEdge.from][A].end())
                    {
                        // finished = false;
                        hashset[currEdge.from][A].insert(nbr);
                        futureNewWorklist.push(EdgeForReading(currEdge.from, nbr, A));
                        edgeVecs[A][currEdge.from].vertexList.push_back(nbr);
                        // newEdgeCounter++;
                    }
                }
            }
        }

        startTransfer = std::chrono::steady_clock::now();
        for (uint g = 0; g < grammar.labelSize; g++)
        {
            // update the sliding/temporal pointers
            for (int i = 0; i < num_nodes; i++)
            {
                edgeVecs[g][i].OLD_END = edgeVecs[g][i].NEW_END;
                edgeVecs[g][i].NEW_END = edgeVecs[g][i].vertexList.size();
            }
        }

        finishTransfer = std::chrono::steady_clock::now();

        std::chrono::duration<double> intervalTransfer = finishTransfer - startTransfer;
        elapsed_seconds_transfer += intervalTransfer.count();

        if (!futureNewWorklist.empty())
        {
            finished = false;
        }

        swap(newWorklist, futureNewWorklist);

        cout << "# Iteration " << itr << endl;
        cout << "futureEdges: " << newWorklist.size() << endl;

    } while (!finished);

    finish = std::chrono::steady_clock::now();
    std::chrono::duration<double> elapsed_seconds = finish - start;
    // std::time_t finish_time = std::chrono::steady_clock::to_time_t(finish);

    uint totalNewEdgeCount = countEdge(hashset, num_nodes, grammar.labelSize) - initialEdgeCount;

    std::cout << "**************************" << endl;
    std::cout << "# Total time = " << elapsed_seconds.count() << std::endl;
    std::cout << "# Total computation time = " << elapsed_seconds_comp << std::endl;
    std::cout << "# Total data transfer time = " << elapsed_seconds_transfer << std::endl;

    std::cout << "SF:: # Number of new edges: " << totalNewEdgeCount << endl;
    std::cout << "AM:: # Number of new edges: " << newEdgeCounter << endl;

    std::cout << "# Iterations: " << itr << endl;
    std::cout << "# Total Calculations: " << calcCnt << endl;
}
