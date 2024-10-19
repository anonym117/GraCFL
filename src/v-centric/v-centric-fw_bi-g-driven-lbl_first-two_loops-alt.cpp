#include "globals-new.hpp"
#include "grammar.hpp"
#include <unordered_map>

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

    vector<vector<unordered_map<ull, ull>>> newEdgeCntMap(grammar.labelSize, vector<unordered_map<ull, ull>>(num_nodes, unordered_map<ull, ull>()));

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
    double stepTime1 = 0.0;
    double stepTime2 = 0.0;

    std::chrono::time_point<std::chrono::steady_clock> start_comp, finish_comp, start_transfer, end_transfer, stepStart1, stepEnd1, stepStart2, stepEnd2;
    start_comp = std::chrono::steady_clock::now();

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
                inEdgeVecs[grammar.grammar1[l][0]][i].vertexList.push_back(i);

                // update the sliding/temporal pointers
                edgeVecs[grammar.grammar1[l][0]][i].NEW_END++;
                inEdgeVecs[grammar.grammar1[l][0]][i].NEW_END++;

                newEdgeCounter++;
            }
        }
    }

    std::chrono::duration<double> interval_comp = (std::chrono::steady_clock::now() - start_comp);
    elapsed_seconds_comp += interval_comp.count();

    cout << "********************\n";

    do
    {
        finished = true;
        itr++;

        start_comp = std::chrono::steady_clock::now();

        // for each grammar rule like A --> B
        stepStart1 = std::chrono::steady_clock::now();

        for (uint g = 0; g < grammar.labelSize; g++)
        {
            double newEdgeCurrItr = 0.0;
            for (uint i = 0; i < num_nodes; i++)
            {
                uint nbr;
                //  the valid index range is [START_NEW, END_NEW-1]
                uint START_NEW = edgeVecs[g][i].OLD_END;
                uint END_NEW = edgeVecs[g][i].NEW_END;

                newEdgeCurrItr += (END_NEW - START_NEW);
                auto mapItr = newEdgeCntMap[g][i].find(itr);

                if (mapItr == newEdgeCntMap[g][i].end())
                {
                    newEdgeCntMap[g][i][itr] = END_NEW - START_NEW;
                }
                else
                {
                    newEdgeCntMap[g][i][itr] += (END_NEW - START_NEW);
                }

                // for each new edge
                for (uint j = START_NEW; j < END_NEW; j++)
                {
                    nbr = edgeVecs[g][i].vertexList[j];
                    // rule: A = B
                    uint g2Size = grammar.grammar2index[g].size();
                    for (uint m = 0; m < g2Size; m++)
                    {
                        calcCnt++;
                        uint A = grammar.grammar2index[g][m];

                        if (hashset[i][A].find(nbr) == hashset[i][A].end())
                        {
                            finished = false;
                            hashset[i][A].insert(nbr);

                            edgeVecs[A][i].vertexList.push_back(nbr);
                            inEdgeVecs[A][nbr].vertexList.push_back(i);

                            // No need to update the pointers. Because the FUTURE_START starts from the
                            // NEW_END, and NEW_END is already updated .
                        }
                    }

                    // rule: A = BC
                    // C (old + new) * B (new)
                    // all the OLD, and NEW outgoing edges of the first edge
                    uint g3SizeLeft = grammar.grammar3indexLeft[g].size();
                    for (uint m = 0; m < g3SizeLeft; m++)
                    {
                        uint C = grammar.grammar3indexLeft[g][m].first;
                        uint A = grammar.grammar3indexLeft[g][m].second;

                        uint START_OLD_OUT = 0;
                        uint END_NEW_OUT = edgeVecs[C][nbr].NEW_END;

                        for (uint h = START_OLD_OUT; h < END_NEW_OUT; h++)
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

            // double newEdgeCounterD = newEdgeCounter;
            // double num_nodesD = num_nodes;
            double avgNewEdge = newEdgeCurrItr / num_nodes;
            // std::cout << std::fixed << std::setprecision(4) << "Itr = " << itr << ", Label = " << g << ", Total NEW edge: "
            //           << newEdgeCurrItr << ", Avg Edge Per Vtx = " << avgNewEdge << endl;
        }

        stepEnd1 = std::chrono::steady_clock::now();
        std::chrono::duration<double> stepInterval1 = stepEnd1 - stepStart1;
        stepTime1 += stepInterval1.count();

        for (uint g = 0; g < grammar.labelSize; g++)
        {
            for (uint i = 0; i < num_nodes; i++)
            {
                uint nbr;
                //  the valid index range is [START_NEW, END_NEW-1]
                uint START_NEW = edgeVecs[g][i].OLD_END;
                uint END_NEW = edgeVecs[g][i].NEW_END;

                for (uint j = START_NEW; j < END_NEW; j++)
                {
                    nbr = edgeVecs[g][i].vertexList[j];
                    uint g3SizeRight = grammar.grammar3indexRight[g].size();
                    for (uint m = 0; m < g3SizeRight; m++)
                    {
                        uint C = grammar.grammar3indexRight[g][m].first;
                        uint A = grammar.grammar3indexRight[g][m].second;

                        uint START_OLD_IN = 0;
                        uint END_OLD_IN = inEdgeVecs[C][i].OLD_END;

                        for (uint h = START_OLD_IN; h < END_OLD_IN; h++)
                        {
                            uint inNbr = inEdgeVecs[C][i].vertexList[h];

                            calcCnt++;
                            if (hashset[inNbr][A].find(nbr) == hashset[inNbr][A].end())
                            {
                                finished = false;
                                hashset[inNbr][A].insert(nbr);

                                edgeVecs[A][inNbr].vertexList.push_back(nbr);
                                inEdgeVecs[A][nbr].vertexList.push_back(inNbr);
                            }
                        }
                    }
                }
            }
        }

        stepEnd2 = std::chrono::steady_clock::now();
        std::chrono::duration<double> stepInterval2 = stepEnd2 - stepEnd1;
        stepTime2 += stepInterval2.count();

        std::chrono::duration<double> interval_comp = (std::chrono::steady_clock::now() - start_comp);
        elapsed_seconds_comp += interval_comp.count();

        std::cout << "Iteration number " << itr << endl;

        start_transfer = std::chrono::steady_clock::now();

        for (uint g = 0; g < grammar.labelSize; g++)
        {
            // update the sliding/temporal pointers
            for (int i = 0; i < num_nodes; i++)
            {
                edgeVecs[g][i].OLD_END = edgeVecs[g][i].NEW_END;
                inEdgeVecs[g][i].OLD_END = inEdgeVecs[g][i].NEW_END;

                // check if new edges are created in this iteration or not.
                // if (edgeVecs[g][i].NEW_END < edgeVecs[g][i].vertexList.size())
                // {
                //     finished = false;
                // }

                edgeVecs[g][i].NEW_END = edgeVecs[g][i].vertexList.size();
                inEdgeVecs[g][i].NEW_END = inEdgeVecs[g][i].vertexList.size();
            }
        }

        std::chrono::duration<double> interval_tranfer = (std::chrono::steady_clock::now() - start_transfer);
        elapsed_seconds_transfer += interval_tranfer.count();

    } while (!finished);

    finish = std::chrono::steady_clock::now();
    std::chrono::duration<double> elapsed_seconds = finish - start;
    // std::time_t finish_time = std::chrono::steady_clock::to_time_t(finish);

    vector<vector<unordered_map<uint, uint>>> revNewEdgeCntMap(grammar.labelSize, vector<unordered_map<uint, uint>>(itr + 1, unordered_map<uint, uint>()));

    for (uint g = 0; g < grammar.labelSize; g++)
    {
        for (uint i = 0; i < num_nodes; i++)
        {
            for (uint j = 1; j <= itr; j++)
            {
                uint newEdgeCnt = newEdgeCntMap[g][i][j];
                if (revNewEdgeCntMap[g][j].find(newEdgeCnt) == revNewEdgeCntMap[g][j].end())
                {
                    revNewEdgeCntMap[g][j][newEdgeCnt] = 1;
                }
                else
                {
                    revNewEdgeCntMap[g][j][newEdgeCnt]++;
                }
            }
        }
    }

    cout << "---------Edge Distribution--------------" << endl;

    double globalSumAvgEdgeVtxP2 = 0.0;
    double globalSumAvgEdgeTotalP2 = 0.0;
    double globalSumAvgEdgeVtxP1 = 0.0;
    double globalSumAvgEdgeTotalP1 = 0.0;

    for (uint j = 1; j <= itr; j++)
    {
        for (uint g = 0; g < grammar.labelSize; g++)
        {
            ull cnt = 0;
            ull edgeVtxTotal = 0;
            ull edgeVtx1 = 0;
            ull edgeVtx2 = 0;

            for (auto itrRev = revNewEdgeCntMap[g][j].begin(); itrRev != revNewEdgeCntMap[g][j].end(); itrRev++)
            {
                if (itrRev->first == 0)
                {
                    // double p = (double)itrRev->second / num_nodes;
                    // std::cout << "Itr = " << j << ", Label = " << g << ", (New edge #" << itrRev->first << ": vtx " << itrRev->second << "), percentage (vtx/totalVtx): " << p * 100 << endl;
                }
                else if (itrRev->first == 1)
                {
                    edgeVtxTotal += itrRev->second;
                    edgeVtx1 = itrRev->second;
                    // double p = (double)itrRev->second / num_nodes;
                    // std::cout << "Itr = " << j << ", Label = " << g << ", (New edge #" << itrRev->first << ": vtx " << itrRev->second << "), percentage (vtx/totalVtx): " << p * 100 << endl;
                }
                else
                {
                    edgeVtxTotal += itrRev->second;
                    edgeVtx2 += itrRev->second;
                }
            }
            double avgTotalEdgeP2 = (double)edgeVtx2 / num_nodes;
            double avgEdgeVtxP2 = (double)edgeVtx2 / edgeVtxTotal;
            double avgEdgeVtxP1 = 1 - avgEdgeVtxP2;

            globalSumAvgEdgeTotalP2 += avgTotalEdgeP2;
            globalSumAvgEdgeTotalP1 += (1 - avgTotalEdgeP2);

            globalSumAvgEdgeVtxP2 += avgEdgeVtxP2;
            globalSumAvgEdgeVtxP1 += avgEdgeVtxP1;

            // std::cout << "Itr = " << j << ", Label = " << g << ", (New edge # >= 2 : vtx " << cnt << "), percentage (vtx/totalVtx): " << p * 100 << endl;
        }
    }

    double globalAvgEdgeTotalP2 = globalSumAvgEdgeTotalP2 / (itr * grammar.labelSize);
    double globalAvgEdgeTotalP1 = globalSumAvgEdgeTotalP1 / (itr * grammar.labelSize);

    double globalAvgEdgeVtxP2 = globalSumAvgEdgeVtxP2 / (itr * grammar.labelSize);
    double globalAvgEdgeVtxP1 = globalSumAvgEdgeVtxP1 / (itr * grammar.labelSize);

    uint totalNewEdgeCount = countEdge(hashset, num_nodes, grammar.labelSize) - initialEdgeCount;

    cout << "**************************" << endl;
    std::cout << "# Total time = " << elapsed_seconds.count() << std::endl;
    std::cout << "# Total computation time = " << elapsed_seconds_comp << std::endl;
    std::cout << "# Total data transfer time = " << elapsed_seconds_transfer << std::endl;
    std::cout << "# Total sum of comp + transfer = " << elapsed_seconds_comp + elapsed_seconds_transfer << std::endl;
    std::cout << "# Total step 1 (FW) time = " << stepTime1 << std::endl;
    std::cout << "# Total step 2 (BI) time = " << stepTime2 << std::endl;

    cout << "SF:: # Number of new edges: " << totalNewEdgeCount << endl;
    cout << "AM:: # Number of new edges: " << newEdgeCounter << endl;

    cout << "# Iterations: " << itr << endl;
    cout << "# Total Calculations: " << calcCnt << endl;

    cout << "globalAvgEdgeTotalP2 (P2Edges / num_nodes) % = " << globalAvgEdgeTotalP2 << endl;
    cout << "globalAvgEdgeTotalP1 (P1Edges / num_nodes) % = " << globalAvgEdgeTotalP1 << endl;
    cout << "No edges % = " << (1 - (globalAvgEdgeTotalP2 + globalAvgEdgeTotalP1)) << endl;
    cout << "globalAvgEdgeVtxP2 (P2Edges / EdgeVtx) % = " << globalAvgEdgeVtxP2 << endl;
    cout << "globalAvgEdgeVtxP1 (P1Edges / EdgeVtx) % = " << globalAvgEdgeVtxP1 << endl;

    getPeakMemoryUsage();
}
