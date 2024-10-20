#include "globals.hpp"
#include "grammar.hpp"

/***
 * V-centric (bi, temporal ptrs, rule-lookup table)
 * Bi directional traversing of the Graph (incoming and outgoing edges) only
 * One buffer list instead of three for OLD, NEW, and FUTURE edges
 * Buffer struct (adjcency list is made with this struct) holds the pointes for the OLD, NEW, and FUTURE edges in the single buffer list
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

	// level-1: vertex ID, level-2: NEW, OLD, FUTURE pointers and incoming edges
	vector<Buffer> inEdgeVecs(num_nodes);
	// level-1: vertex ID, level-2: NEW, OLD, FUTURE pointers and outgoing edges
	vector<Buffer> edgeVecs(num_nodes);

	// check if an edge exist or not
	// unordered_set<ull> *hashset = new unordered_set<ull>[num_nodes];
	vector<vector<unordered_set<ull>>> hashset(num_nodes, vector<unordered_set<ull>>(grammar.labelSize, unordered_set<ull>()));

	for (uint i = 0; i < num_edges; i++)
	{
		edgeVecs[edges[i].from].vertexList.push_back(OutEdge(edges[i].to, edges[i].label));
		inEdgeVecs[edges[i].to].vertexList.push_back(OutEdge(edges[i].from, edges[i].label));

		// update the temporal pointers
		edgeVecs[edges[i].from].NEW_END++;
		inEdgeVecs[edges[i].to].NEW_END++;

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
	for (uint l = 0; l < grammar.grammar1.size(); l++)
	{
		for (uint i = 0; i < num_nodes; i++)
		{
			// check if the new edge based on an epsilon grammar rule exists or not. l: grammar ID, 0: LHS
			if (hashset[i][grammar.grammar1[l][0]].find(i) == hashset[i][grammar.grammar1[l][0]].end())
			{
				// insert into hashset
				hashset[i][grammar.grammar1[l][0]].insert(i);
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
			uint START_NEW_IN = inEdgeVecs[i].OLD_END;
			uint END_NEW_IN = inEdgeVecs[i].NEW_END;

			// for each new edge
			for (uint j = START_NEW_IN; j < END_NEW_IN; j++)
			{
				inNbr = inEdgeVecs[i].vertexList[j];
				// if the edge to the neighbor is labeled with B
				vector<uint> leftLabels = grammar.rule2(inNbr.label);

				for (uint g = 0; g < leftLabels.size(); g++)
				{
					if (hashset[inNbr.end][leftLabels[g]].find(i) == hashset[inNbr.end][leftLabels[g]].end())
					{
						finished = false;
						hashset[inNbr.end][leftLabels[g]].insert(i);

						edgeVecs[inNbr.end].vertexList.push_back(OutEdge(i, leftLabels[g]));
						inEdgeVecs[i].vertexList.push_back(OutEdge(inNbr.end, leftLabels[g]));
					}
				}

				// OLD incoming neighbors of the first edge
				uint START_OLD = 0;
				uint END_NEW = edgeVecs[i].NEW_END;

				for (uint h = START_OLD; h < END_NEW; h++)
				{
					OutEdge outNbr = edgeVecs[i].vertexList[h];
					vector<uint> leftLabels = grammar.rule3(inNbr.label, outNbr.label);

					for (uint g = 0; g < leftLabels.size(); g++)
					{
						if (hashset[inNbr.end][leftLabels[g]].find(outNbr.end) == hashset[inNbr.end][leftLabels[g]].end())
						{
							finished = false;
							hashset[inNbr.end][leftLabels[g]].insert(outNbr.end);
							edgeVecs[inNbr.end].vertexList.push_back(OutEdge(outNbr.end, leftLabels[g]));
							inEdgeVecs[outNbr.end].vertexList.push_back(OutEdge(inNbr.end, leftLabels[g]));
						}
					}
				}
			}

			uint NEW_START_OUT = edgeVecs[i].OLD_END;
			uint NEW_END_OUT = edgeVecs[i].NEW_END;

			for (uint j = NEW_START_OUT; j < NEW_END_OUT; j++)
			{
				OutEdge outNbr = edgeVecs[i].vertexList[j];

				uint OLD_START_IN = 0;
				uint OLD_END_IN = inEdgeVecs[i].OLD_END;

				for (uint h = OLD_START_IN; h < OLD_END_IN; h++)
				{
					OutEdge inNbr = inEdgeVecs[i].vertexList[h];
					vector<uint> leftLabels = grammar.rule3(inNbr.label, outNbr.label);

					for (uint g = 0; g < leftLabels.size(); g++)
					{
						if (hashset[inNbr.end][leftLabels[g]].find(outNbr.end) == hashset[inNbr.end][leftLabels[g]].end())
						{
							finished = false;
							hashset[inNbr.end][leftLabels[g]].insert(outNbr.end);
							edgeVecs[inNbr.end].vertexList.push_back(OutEdge(outNbr.end, leftLabels[g]));
							inEdgeVecs[outNbr.end].vertexList.push_back(OutEdge(inNbr.end, leftLabels[g]));
						}
					}
				}
			}
		}

		for (int i = 0; i < num_nodes; i++)
		{
			// update the pointers
			edgeVecs[i].OLD_END = edgeVecs[i].NEW_END;
			inEdgeVecs[i].OLD_END = inEdgeVecs[i].NEW_END;

			edgeVecs[i].NEW_END = edgeVecs[i].vertexList.size();
			inEdgeVecs[i].NEW_END = inEdgeVecs[i].vertexList.size();
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
