#include "globals.hpp"
#include "grammar.hpp"

/***
 * Bi directional traversing of the Graph (incoming and outgoing edges)
 * Topology Driven
 * Three separate buffer lists for OLD, NEW, and FUTURE edges
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

	vector<vector<vector<OutEdge>>> inEdgeVecs(num_nodes, vector<vector<OutEdge>>(3));
	// level-1: vertex ID, level-2: NEW, OLD, FUTURE, level-3: outgoing edges
	vector<vector<vector<OutEdge>>> edgeVecs(num_nodes, vector<vector<OutEdge>>(3));
	// check if an edge exist or not
	unordered_set<ull> *hashset = new unordered_set<ull>[num_nodes];
	// hashset for incoming edges

	for (uint i = 0; i < num_edges; i++)
	{
		edgeVecs[edges[i].from][NEW].push_back(OutEdge(edges[i].to, edges[i].label));
		inEdgeVecs[edges[i].to][NEW].push_back(OutEdge(edges[i].from, edges[i].label));

		// edgeVecs[edges[i].from][OLD].push_back(OutEdge(edges[i].to, edges[i].label));
		// inEdgeVecs[edges[i].to][OLD].push_back(OutEdge(edges[i].from, edges[i].label));

		hashset[edges[i].from].insert(COMBINE(edges[i].to, edges[i].label));
		// inHashset[edges[i].to].insert(COMBINE(edges[i].from, edges[i].label));
	}

	edges.clear();

	// SF:: count the total no of initial unique edge
	uint initialEdgeCount = countEdge(hashset, num_nodes);

	std::cout << "Initialization Done!" << std::endl;

	bool finished; // fixed-point iteration flag
	int itr = 0; // Iteration counter for fixed-point iteration

	std::cout << "Start Calculations...\n";

	// currently exclude the initialization time
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
				// insert into the graph
				edgeVecs[i][NEW].push_back(OutEdge(i, grammar.grammar1[l][0]));
				inEdgeVecs[i][NEW].push_back(OutEdge(i, grammar.grammar1[l][0]));
			}
		}
	}

	do
	{
		itr++;
		std::cout << "Iteration number " << itr << std::endl;

		// for each grammar rule like A --> B
		for (uint i = 0; i < num_nodes; i++)
		{
			for (uint j = 0; j < inEdgeVecs[i][NEW].size(); j++)
			{
				OutEdge inNbr = inEdgeVecs[i][NEW][j];

				for (uint g = 0; g < grammar.grammar2.size(); g++)
				{
					if (inNbr.label == grammar.grammar2[g][1])
					{
						if (hashset[inNbr.end].find(COMBINE(i, grammar.grammar2[g][0])) == hashset[inNbr.end].end())
						{
							hashset[inNbr.end].insert(COMBINE(i, grammar.grammar2[g][0]));
							edgeVecs[inNbr.end][FUTURE].push_back(OutEdge(i, grammar.grammar2[g][0]));
							inEdgeVecs[i][FUTURE].push_back(OutEdge(inNbr.end, grammar.grammar2[g][0]));
						}
					}
				}

				// grammar3 calc
				// A= BC; B = first edge, C= outgoing neighbor B
				for (uint h = 0; h < edgeVecs[i][NEW].size(); h++)
				{
					OutEdge outNbr = edgeVecs[i][NEW][h];
					for (uint g = 0; g < grammar.grammar3.size(); g++)
					{
						if (inNbr.label == grammar.grammar3[g][1] && outNbr.label == grammar.grammar3[g][2])
						{
							if (hashset[inNbr.end].find(COMBINE(outNbr.end, grammar.grammar3[g][0])) == hashset[inNbr.end].end())
							{
								hashset[inNbr.end].insert(COMBINE(outNbr.end, grammar.grammar3[g][0]));
								// inHashset[outNbr.end].insert(COMBINE(i, grammar.grammar3[g][0]));
								edgeVecs[inNbr.end][FUTURE].push_back(OutEdge(outNbr.end, grammar.grammar3[g][0]));
								inEdgeVecs[outNbr.end][FUTURE].push_back(OutEdge(inNbr.end, grammar.grammar3[g][0]));
							}
						}
					}
				}

				for (uint h = 0; h < edgeVecs[i][OLD].size(); h++)
				{
					OutEdge outNbr = edgeVecs[i][OLD][h];
					for (uint g = 0; g < grammar.grammar3.size(); g++)
					{
						if (inNbr.label == grammar.grammar3[g][1] && outNbr.label == grammar.grammar3[g][2])
						{
							if (hashset[inNbr.end].find(COMBINE(outNbr.end, grammar.grammar3[g][0])) == hashset[inNbr.end].end())
							{
								hashset[inNbr.end].insert(COMBINE(outNbr.end, grammar.grammar3[g][0]));
								// inHashset[outNbr.end].insert(COMBINE(i, grammar.grammar3[g][0]));
								edgeVecs[inNbr.end][FUTURE].push_back(OutEdge(outNbr.end, grammar.grammar3[g][0]));
								inEdgeVecs[outNbr.end][FUTURE].push_back(OutEdge(inNbr.end, grammar.grammar3[g][0]));
							}
						}
					}
				}
			}

			for (uint j = 0; j < edgeVecs[i][NEW].size(); j++)
			{
				OutEdge outNbr = edgeVecs[i][NEW][j];

				for (uint h = 0; h < inEdgeVecs[i][OLD].size(); h++)
				{
					OutEdge inNbr = inEdgeVecs[i][OLD][h];
					for (uint g = 0; g < grammar.grammar3.size(); g++)
					{
						if (inNbr.label == grammar.grammar3[g][1] && outNbr.label == grammar.grammar3[g][2])
						{
							if (hashset[inNbr.end].find(COMBINE(outNbr.end, grammar.grammar3[g][0])) == hashset[inNbr.end].end())
							{
								hashset[inNbr.end].insert(COMBINE(outNbr.end, grammar.grammar3[g][0]));
								// inHashset[outNbr.end].insert(COMBINE(i, grammar.grammar3[g][0]));
								edgeVecs[inNbr.end][FUTURE].push_back(OutEdge(outNbr.end, grammar.grammar3[g][0]));
								inEdgeVecs[outNbr.end][FUTURE].push_back(OutEdge(inNbr.end, grammar.grammar3[g][0]));
							}
						}
					}
				}

				// edgeVecs[i][OLD].push_back(OutEdge(nbr.end, nbr.label));
				// inEdgeVecs[nbr.end][OLD].push_back(OutEdge(i, nbr.label));
			}
		}

		finished = true;

		// add NEW to OLD, then copy FUTURE to NEW
		for (int i = 0; i < num_nodes; i++)
		{
			edgeVecs[i][OLD].insert(edgeVecs[i][OLD].begin(), edgeVecs[i][NEW].begin(), edgeVecs[i][NEW].end());
			inEdgeVecs[i][OLD].insert(inEdgeVecs[i][OLD].begin(), inEdgeVecs[i][NEW].begin(), inEdgeVecs[i][NEW].end());

			edgeVecs[i][NEW] = edgeVecs[i][FUTURE];
			inEdgeVecs[i][NEW] = inEdgeVecs[i][FUTURE];
			if (edgeVecs[i][NEW].size() > 0 || inEdgeVecs[i][NEW].size() > 0)
				finished = false;
			edgeVecs[i][FUTURE].clear();
			inEdgeVecs[i][FUTURE].clear();
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
