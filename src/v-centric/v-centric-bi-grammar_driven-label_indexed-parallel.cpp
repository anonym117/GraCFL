#include <omp.h>
#include "globals.hpp"
#include "globals-par.hpp"
#include "grammar.hpp"

/***
 * V-centric (bi, sliding ptrs, grammar driven)
 * adjacency list: first level: grammar label, second level: vertex
 * Bi directional traversing of the Graph (incoming and outgoing edges)
 * One buffer list instead of three for OLD, NEW, and FUTURE edges
 * ConcurrBuff struct (adjcency list is made with this struct) holds the pointes for the OLD, NEW, and FUTURE edges in the single buffer list
 * Check the Future edge flag whenever an edge is created
 */

int main(int argc, char **argv)
{

    // Get the graph file path and grammar file path from command line argument: output file path is optional
	if (argc == 1)
	{
		std::cout << "Please provide the graph file path and grammar file path. For example, ./exc.out <graph_file> <grammar_file>" << std::endl;
		return 0;
	}
	else if (argc == 2)
	{
		std::cout << "Please provide the grammar file path. For example, ./exc.out <graph_file> <grammar_file>" << std::endl;
		return 0;
	}

	std::cout << "-----------START----------" << std::endl;
	std::cout << "--------------------------" << std::endl;

	const std::string inputGraph = argv[1];
	std::cout << "GraphFile:\t" << inputGraph << std::endl;

	std::string grammarFilePath = argv[2];
	std::cout << "GrammarFile:\t" << grammarFilePath << endl;
	uint TOTAL_THREADS = 32;

	if (argc == 4)
	{
		TOTAL_THREADS = stoi(argv[3]);
	}

	std::cout << "No of Threads:\t" << TOTAL_THREADS << std::endl; 
	std::cout << "--------------------------" << std::endl;
    
	Grammar grammar(grammarFilePath); // Read grammar

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

    std::cout << "# Vertex Count:\t" << num_nodes << std::endl;
	std::cout << "# Initial Edge Count:\t" << num_edges << std::endl;
	std::cout << "Start initializing the lists, hashset and worklists ..." << std::endl;

    // level-1: vertex ID, level-2: NEW, OLD, FUTURE pointers and outgoing edges
    ConcurrBuff **edgeVecs = new ConcurrBuff *[grammar.labelSize];
    ConcurrBuff **inEdgeVecs = new ConcurrBuff *[grammar.labelSize];

    for (uint i = 0; i < grammar.labelSize; i++)
    {
        edgeVecs[i] = new ConcurrBuff[num_nodes];
        inEdgeVecs[i] = new ConcurrBuff[num_nodes];
    }

    tbb::concurrent_unordered_set<ull> **hashset = new tbb::concurrent_unordered_set<ull> *[num_nodes];
    for (uint i = 0; i < num_nodes; i++)
    {
        hashset[i] = new tbb::concurrent_unordered_set<ull>[grammar.labelSize];
    }

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

    std::cout << "Initialization Done!" << std::endl;

    bool finished; // fixed-point iteration flag
	int itr = 0; // Iteration counter for fixed-point iteration

	std::cout << "Start Calculations...\n";
    std::chrono::time_point<std::chrono::steady_clock> start, finish;
    start = std::chrono::steady_clock::now();
    
    // handle epsilon rules: add an edge to itself
    for (uint l = 0; l < grammar.grammar1.size(); l++)
    {
#pragma omp parallel for schedule(static) num_threads(TOTAL_THREADS)
        for (uint i = 0; i < num_nodes; i++)
        {
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
            }
        }
    }

    do
    {
        finished = true;
        itr++;
        std::cout << "Iteration " << itr << std::endl;

        // for each grammar rule like A --> B
        for (uint g = 0; g < grammar.labelSize; g++)
        {
#pragma omp parallel for schedule(static) num_threads(TOTAL_THREADS)
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
                    // rule: A = B
                    uint g2Size = grammar.grammar2index[g].size();
                    for (uint m = 0; m < g2Size; m++)
                    {
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

                    // rule: A = CB
                    // C (old + new) * B (new)
                    // all the OLD, and NEW outgoing edges of the first edge
                    uint g3SizeRight = grammar.grammar3indexRight[g].size();
                    for (uint m = 0; m < g3SizeRight; m++)
                    {
                        uint C = grammar.grammar3indexRight[g][m].first;
                        uint A = grammar.grammar3indexRight[g][m].second;

                        uint START_OLD_OUT = 0;
                        uint END_NEW_OUT = inEdgeVecs[C][i].NEW_END;

                        for (uint h = START_OLD_OUT; h < END_NEW_OUT; h++)
                        {
                            uint inNbr = inEdgeVecs[C][i].vertexList[h];

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

                // Rule: A = BC
                uint START_NEW_IN = inEdgeVecs[g][i].OLD_END;
                uint END_NEW_IN = inEdgeVecs[g][i].NEW_END;

                for (uint j = START_NEW_IN; j < END_NEW_IN; j++)
                {
                    uint inNbr = inEdgeVecs[g][i].vertexList[j];

                    uint g3SizeLeft = grammar.grammar3indexLeft[g].size();
                    for (uint m = 0; m < g3SizeLeft; m++)
                    {
                        uint C = grammar.grammar3indexLeft[g][m].first;
                        uint A = grammar.grammar3indexLeft[g][m].second;

                        uint START_OLD_OUT = 0;
                        uint END_OLD_OUT = edgeVecs[C][i].OLD_END;

                        for (uint h = START_OLD_OUT; h < END_OLD_OUT; h++)
                        {
                            uint nbr = edgeVecs[C][i].vertexList[h];

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

        for (uint g = 0; g < grammar.labelSize; g++)
        {
            // update the sliding/temporal pointers
#pragma omp parallel for schedule(static) num_threads(TOTAL_THREADS)
            for (int i = 0; i < num_nodes; i++)
            {
                edgeVecs[g][i].OLD_END = edgeVecs[g][i].NEW_END;
                inEdgeVecs[g][i].OLD_END = inEdgeVecs[g][i].NEW_END;

                edgeVecs[g][i].NEW_END = edgeVecs[g][i].vertexList.size();
                inEdgeVecs[g][i].NEW_END = inEdgeVecs[g][i].vertexList.size();
            }
        }
    } while (!finished);

    finish = std::chrono::steady_clock::now();
    std::chrono::duration<double> elapsed_seconds = finish - start;
    std::cout << "Calculation Done!" << std::endl;

    uint totalNewEdgeCount = countEdge(hashset, num_nodes, grammar.labelSize) - initialEdgeCount;

    std::cout << "----------RESULTS----------" << std::endl;
	std::cout << "Graph File: " << inputGraph << std::endl;
	std::cout << "--------------------------" << std::endl;
	std::cout << "# Total Calculation Time =\t" << elapsed_seconds.count() << std::endl;
	std::cout << "# Total NEW Edge Created =\t" << totalNewEdgeCount << std::endl;
	std::cout << "# Total Iterations =\t" << itr << std::endl;

	// Get memory usage: comment this if not needed
	getPeakMemoryUsage();
	std::cout << "--------------------------" << std::endl;
	std::cout << "------------END-----------" << std::endl;
    
    // Delete dynamically alloted memories
    for (uint i = 0; i < grammar.labelSize; i++)
    {
        delete[] edgeVecs[i];
        delete[] inEdgeVecs[i];
    }

    for (uint i = 0; i < num_nodes; i++)
    {
        delete[] hashset[i];
    }

    delete[] edgeVecs;
    delete[] inEdgeVecs;
    delete[] hashset;
}
