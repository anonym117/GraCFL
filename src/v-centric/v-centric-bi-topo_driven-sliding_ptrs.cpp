#include "globals.hpp"
#include "grammar.hpp"

/***
 * Bi directional traversing of the Graph (incoming and outgoing edges)
 * One buffer list instead of three for OLD, NEW, and FUTURE edges
 * Buffer struct (adjcency list is made with this struct) holds the pointes for the OLD, NEW, and FUTURE edges in the single buffer list
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

	// level-1: vertex ID, level-2: NEW, OLD, FUTURE pointers and incoming edges
	vector<Buffer> inEdgeVecs(num_nodes);
	// level-1: vertex ID, level-2: NEW, OLD, FUTURE pointers and outgoing edges
	vector<Buffer> edgeVecs(num_nodes);

	// check if an edge exist or not
	unordered_set<ull> *hashset = new unordered_set<ull>[num_nodes];
	// hashset for incoming edges
	// unordered_set<ull> *inHashset = new unordered_set<ull>[num_nodes];

	for (uint i = 0; i < num_edges; i++)
	{
		edgeVecs[edges[i].from].vertexList.push_back(OutEdge(edges[i].to, edges[i].label));
		inEdgeVecs[edges[i].to].vertexList.push_back(OutEdge(edges[i].from, edges[i].label));

		// update the buffer pointers
		edgeVecs[edges[i].from].NEW_END++;
		inEdgeVecs[edges[i].to].NEW_END++;

		hashset[edges[i].from].insert(COMBINE(edges[i].to, edges[i].label));
		// inHashset[edges[i].to].insert(COMBINE(edges[i].from, edges[i].label));
	}

	edges.clear();

	// count the total no of initial unique edge
	uint initialEdgeCount = countEdge(hashset, num_nodes);

	std::cout << "Initialization Done!" << std::endl;

	bool finished; // fixed-point iteration flag
	int itr = 0; // Iteration counter for fixed-point iteration

	std::cout << "Start Calculations...\n";
	std::chrono::time_point<std::chrono::system_clock> start, finish;
	start = std::chrono::system_clock::now();

	// handle epsilon rules: add an edge to itself
	for (uint l = 0; l < grammar.grammar1.size(); l++)
	{
		for (uint i = 0; i < num_nodes; i++)
		{
			// check if the new edge based on an epsilon grammar rule exists or not. l: grammar ID, 0: LHS
			if (hashset[i].find(COMBINE(i, grammar.grammar1[l][0])) == hashset[i].end())
			{
				// insert into hashset
				hashset[i].insert(COMBINE(i, grammar.grammar1[l][0]));
				// inHashset[i].insert(COMBINE(i, grammar.grammar1[l][0]));
				// insert into the graph
				edgeVecs[i].vertexList.push_back(OutEdge(i, grammar.grammar1[l][0]));
				inEdgeVecs[i].vertexList.push_back(OutEdge(i, grammar.grammar1[l][0]));

				// update the pointers
				edgeVecs[i].NEW_END++;
				inEdgeVecs[i].NEW_END++;
			}
		}
	}

	do
	{
		finished = true;
		itr++;
		std::cout << "Iteration number " << itr << std::endl;

		// for each grammar rule like A --> B
		for (uint i = 0; i < num_nodes; i++)
		{
			OutEdge inNbr;
			//  the valid index range is [START_NEW, END_NEW-1]
			uint START_NEW = inEdgeVecs[i].OLD_END;
			uint END_NEW = inEdgeVecs[i].NEW_END;

			// for each new edge
			for (uint j = START_NEW; j < END_NEW; j++)
			{
				inNbr = inEdgeVecs[i].vertexList[j];
				// if the edge to the neighbor is labeled with B

				for (uint g = 0; g < grammar.grammar2.size(); g++)
				{
					if (inNbr.label == grammar.grammar2[g][1])
					{
						if (hashset[inNbr.end].find(COMBINE(i, grammar.grammar2[g][0])) == hashset[inNbr.end].end())
						{
							finished = false;
							hashset[inNbr.end].insert(COMBINE(i, grammar.grammar2[g][0]));
							
							edgeVecs[inNbr.end].vertexList.push_back(OutEdge(i, grammar.grammar2[g][0]));
							inEdgeVecs[i].vertexList.push_back(OutEdge(inNbr.end, grammar.grammar2[g][0]));

							// No need to update the pointers. Because the FUTURE_START starts from the
							// NEW_END, and NEW_END is already updated.
						}
					}
				}

				// all the OLD, and NEW outgoing edges of the first edge
				uint START_OLD_OUT = 0;
				uint END_NEW_OUT = edgeVecs[i].NEW_END;
				// grammar3 calc
				for (uint h = START_OLD_OUT; h < END_NEW_OUT; h++)
				{
					OutEdge outNbr = edgeVecs[i].vertexList[h];
					for (uint g = 0; g < grammar.grammar3.size(); g++)
					{
						if (inNbr.label == grammar.grammar3[g][1] && outNbr.label == grammar.grammar3[g][2])
						{
							if (hashset[inNbr.end].find(COMBINE(outNbr.end, grammar.grammar3[g][0])) == hashset[inNbr.end].end())
							{
								finished = false;
								hashset[inNbr.end].insert(COMBINE(outNbr.end, grammar.grammar3[g][0]));
								
								edgeVecs[inNbr.end].vertexList.push_back(OutEdge(outNbr.end, grammar.grammar3[g][0]));
								inEdgeVecs[outNbr.end].vertexList.push_back(OutEdge(inNbr.end, grammar.grammar3[g][0]));
								// No need to update the pointers. Because the FUTURE_START starts from the
								// NEW_END, and NEW_END is already updated.
							}
						}
					}
				}
			}

			uint START_NEW_OUT = edgeVecs[i].OLD_END;
			uint END_NEW_OUT = edgeVecs[i].NEW_END;

			for (uint j = START_NEW_OUT; j < END_NEW_OUT; j++)
			{

				OutEdge outNbr = edgeVecs[i].vertexList[j];
				// OLD incoming neighbors of the first edge
				uint START_OLD_IN = 0;
				uint END_OLD_IN = inEdgeVecs[i].OLD_END;

				for (uint h = START_OLD_IN; h < END_OLD_IN; h++)
				{
					OutEdge inNbr = inEdgeVecs[i].vertexList[h];
					for (uint g = 0; g < grammar.grammar3.size(); g++)
					{
						if (inNbr.label == grammar.grammar3[g][1] && outNbr.label == grammar.grammar3[g][2])
						{
							if (hashset[inNbr.end].find(COMBINE(outNbr.end, grammar.grammar3[g][0])) == hashset[inNbr.end].end())
							{
								finished = false;
								hashset[inNbr.end].insert(COMBINE(outNbr.end, grammar.grammar3[g][0]));
								
								edgeVecs[inNbr.end].vertexList.push_back(OutEdge(outNbr.end, grammar.grammar3[g][0]));
								inEdgeVecs[outNbr.end].vertexList.push_back(OutEdge(inNbr.end, grammar.grammar3[g][0]));
								// No need to update the pointers. Because the FUTURE_START starts from the
								// NEW_END, and NEW_END is already updated.
							}
						}
					}
				}
			}
		}

		// update the sliding pointers
		for (int i = 0; i < num_nodes; i++)
		{
			edgeVecs[i].OLD_END = edgeVecs[i].NEW_END;
			inEdgeVecs[i].OLD_END = inEdgeVecs[i].NEW_END;

			edgeVecs[i].NEW_END = edgeVecs[i].vertexList.size();
			inEdgeVecs[i].NEW_END = inEdgeVecs[i].vertexList.size();
		}
	} while (!finished);

	finish = std::chrono::system_clock::now();
	std::chrono::duration<double> elapsed_seconds = finish - start;
	std::cout << "Calculation Done!" << std::endl;

	uint totalNewEdgeCount = countEdge(hashset, num_nodes) - initialEdgeCount;

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
