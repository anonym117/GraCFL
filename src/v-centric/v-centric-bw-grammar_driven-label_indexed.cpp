#include "globals.hpp"
#include "grammar.hpp"

/***
 * V-centric (backward, sliding ptrs, grammar driven)
 * Backward directional traversing of the Graph (incoming edges)
 * One BufferEdge list instead of three for OLD, NEW, and FUTURE edges
 * BufferEdge struct (adjcency list is made with this struct) holds the pointes for the OLD, NEW, and FUTURE edges in the single BufferEdge list
 * Check the Future edge flag whenever an edge is created
 * Grammar index
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
    vector<vector<BufferEdge>> inEdgeVecs(grammar.labelSize, vector<BufferEdge>(num_nodes));
    vector<vector<unordered_set<ull>>> inHashset(num_nodes, vector<unordered_set<ull>>(grammar.labelSize, unordered_set<ull>()));

    for (uint i = 0; i < num_edges; i++)
    {
        inEdgeVecs[edges[i].label][edges[i].to].vertexList.push_back(edges[i].from);
        // update the BufferEdge pointers
        inEdgeVecs[edges[i].label][edges[i].to].NEW_END++;
        // insert edge into the hashset
        inHashset[edges[i].to][edges[i].label].insert(edges[i].from);
    }

    edges.clear();

    // count the total no of initial unique edge
    uint initialEdgeCount = countEdge(inHashset, num_nodes, grammar.labelSize);

    std::cout << "Initialization Done!" << std::endl;

	bool finished; // fixed-point iteration flag
	int itr = 0; // Iteration counter for fixed-point iteration

	std::cout << "Start Calculations...\n";

    // currently exclude the initialization time
    std::chrono::time_point<std::chrono::steady_clock> start, finish;
    start = std::chrono::steady_clock::now();

    // handle epsilon rules: add an edge to itself
    for (uint l = 0; l < grammar.grammar1.size(); l++)
    {
        for (uint i = 0; i < num_nodes; i++)
        {
            // check if the new edge based on an epsilon grammar rule exists or not. l: grammar ID, 0: LHS
            if (inHashset[i][grammar.grammar1[l][0]].find(i) == inHashset[i][grammar.grammar1[l][0]].end())
            {
                // insert edge into hashset
                inHashset[i][grammar.grammar1[l][0]].insert(i);
                // insert edge into the graph
                inEdgeVecs[grammar.grammar1[l][0]][i].vertexList.push_back(i);
                // update the sliding/temporal pointers
                inEdgeVecs[grammar.grammar1[l][0]][i].NEW_END++;
            }
        }
    }

    do
    {
        finished = true;
        itr++;
        std::cout << "Iteration number " << itr << std::endl;

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
                        uint A = grammar.grammar2index[g][m];

                        if (inHashset[i][A].find(inNbr1) == inHashset[i][A].end())
                        {
                            finished = false;
                            inHashset[i][A].insert(inNbr1);
                            inEdgeVecs[A][i].vertexList.push_back(inNbr1);
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

                            if (inHashset[i][A].find(inNbr2) == inHashset[i][A].end())
                            {
                                finished = false;
                                inHashset[i][A].insert(inNbr2);
                                inEdgeVecs[A][i].vertexList.push_back(inNbr2);
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

                            if (inHashset[i][A].find(inNbr2) == inHashset[i][A].end())
                            {
                                finished = false;
                                inHashset[i][A].insert(inNbr2);
                                inEdgeVecs[A][i].vertexList.push_back(inNbr2);
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
                inEdgeVecs[g][i].OLD_END = inEdgeVecs[g][i].NEW_END;
                inEdgeVecs[g][i].NEW_END = inEdgeVecs[g][i].vertexList.size();
            }
        }
    } while (!finished);

    finish = std::chrono::steady_clock::now();
    std::chrono::duration<double> elapsed_seconds = finish - start;
    std::cout << "Calculation Done!" << std::endl;

    uint totalNewEdgeCount = countEdge(inHashset, num_nodes, grammar.labelSize) - initialEdgeCount;

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
