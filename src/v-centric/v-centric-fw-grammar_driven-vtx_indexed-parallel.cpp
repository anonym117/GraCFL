#include <omp.h>
#include "globals.hpp"
#include "grammar.hpp"

#define TOTAL_THREADS 16
/***
 * V-centric (forward, sliding ptrs, grammar driven)
 * Parallel
 * First Level Index: Vertex, Second Level Index: Edge Labels
 * Forward directional traversing of the Graph (outgoing edges)
 * One BufferEdge list instead of three for OLD, NEW, and FUTURE edges
 * BufferEdge struct (adjcency list is made with this struct) holds the pointes for the OLD, NEW, and FUTURE edges in the single BufferEdge list
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
    vector<vector<BufferEdge>> edgeVecs(num_nodes, vector<BufferEdge>(grammar.labelSize));

    vector<vector<unordered_set<ull>>> hashset(num_nodes, vector<unordered_set<ull>>(grammar.labelSize, unordered_set<ull>()));
    // hashset for incoming edges
    // unordered_set<ull> *inHashset = new unordered_set<ull>[num_nodes];

    for (uint i = 0; i < num_edges; i++)
    {
        edgeVecs[edges[i].from][edges[i].label].vertexList.push_back(edges[i].to);
        // update the BufferEdge pointers
        edgeVecs[edges[i].from][edges[i].label].NEW_END++;
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
	std::chrono::time_point<std::chrono::system_clock> start, finish;
	start = std::chrono::system_clock::now();

// handle epsilon rules: add an edge to itself
#pragma omp parallel for schedule(static) num_threads(TOTAL_THREADS)
    for (uint i = 0; i < num_nodes; i++)
    {
        for (uint l = 0; l < grammar.grammar1.size(); l++)
        {
            // check if the new edge based on an epsilon grammar rule exists or not. l: grammar ID, 0: LHS
            if (hashset[i][grammar.grammar1[l][0]].find(i) == hashset[i][grammar.grammar1[l][0]].end())
            {
                // insert edge into hashset
                hashset[i][grammar.grammar1[l][0]].insert(i);
                // insert edge into the graph
                edgeVecs[i][grammar.grammar1[l][0]].vertexList.push_back(i);
                // update the sliding/temporal pointers
                edgeVecs[i][grammar.grammar1[l][0]].NEW_END++;
            }
        }
    }

    do
    {
        finished = true;
        itr++;
        std::cout << "Iteration number " << itr << std::endl;

// for each grammar rule like A --> B
#pragma omp parallel for schedule(static) num_threads(TOTAL_THREADS)
        for (uint i = 0; i < num_nodes; i++)
        {
            for (uint g = 0; g < grammar.labelSize; g++)
            {
                uint nbr;
                //  the valid index range is [START_NEW, END_NEW-1]
                uint START_NEW = edgeVecs[i][g].OLD_END;
                uint END_NEW = edgeVecs[i][g].NEW_END;

                // for each new edge
                for (uint j = START_NEW; j < END_NEW; j++)
                {
                    nbr = edgeVecs[i][g].vertexList[j];
                    // rule: A = B
                    for (uint m = 0; m < grammar.grammar2index[g].size(); m++)
                    {
                        uint A = grammar.grammar2index[g][m];

                        if (hashset[i][A].find(nbr) == hashset[i][A].end())
                        {
                            finished = false;
                            hashset[i][A].insert(nbr);

                            edgeVecs[i][A].vertexList.push_back(nbr);
                            // No need to update the pointers. Because the FUTURE_START starts from the
                            // NEW_END, and NEW_END is already updated .
                        }
                    }

                    // rule: A = BC
                    // all the OLD, and NEW outgoing edges of the first edge
                    for (uint m = 0; m < grammar.grammar3indexLeft[g].size(); m++)
                    {
                        uint C = grammar.grammar3indexLeft[g][m].first;
                        uint A = grammar.grammar3indexLeft[g][m].second;

                        uint START_OLD_OUT = 0;
                        uint END_NEW_OUT = edgeVecs[nbr][C].NEW_END;

                        for (uint h = START_OLD_OUT; h < END_NEW_OUT; h++)
                        {
                            uint outNbr = edgeVecs[nbr][C].vertexList[h];

                            if (hashset[i][A].find(outNbr) == hashset[i][A].end())
                            {
                                finished = false;
                                hashset[i][A].insert(outNbr);

                                edgeVecs[i][A].vertexList.push_back(outNbr);
                            }
                        }
                    }
                }

                uint START_OLD = 0;
                uint END_OLD = edgeVecs[i][g].OLD_END;

                // for each old edge
                for (uint j = START_OLD; j < END_OLD; j++)
                {
                    nbr = edgeVecs[i][g].vertexList[j];

                    // rule: A = BC
                    // all the NEW outgoing edges of the first edge
                    for (uint m = 0; m < grammar.grammar3indexLeft[g].size(); m++)
                    {
                        uint C = grammar.grammar3indexLeft[g][m].first;
                        uint A = grammar.grammar3indexLeft[g][m].second;

                        uint START_NEW_OUT = edgeVecs[nbr][C].OLD_END;
                        uint END_NEW_OUT = edgeVecs[nbr][C].NEW_END;

                        for (uint h = START_NEW_OUT; h < END_NEW_OUT; h++)
                        {
                            uint outNbr = edgeVecs[nbr][C].vertexList[h];

                            if (hashset[i][A].find(outNbr) == hashset[i][A].end())
                            {
                                finished = false;
                                hashset[i][A].insert(outNbr);

                                edgeVecs[i][A].vertexList.push_back(outNbr);
                            }
                        }
                    }
                }
            }
        }

#pragma omp parallel for schedule(static) num_threads(TOTAL_THREADS)
        for (int i = 0; i < num_nodes; i++)
        {
            // update the sliding/temporal pointers
            for (uint g = 0; g < grammar.labelSize; g++)
            {
                edgeVecs[i][g].OLD_END = edgeVecs[i][g].NEW_END;
                edgeVecs[i][g].NEW_END = edgeVecs[i][g].vertexList.size();
            }
        }

    } while (!finished);

    finish = std::chrono::system_clock::now();
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
}