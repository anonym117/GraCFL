#include "globals.hpp"
#include "grammar.hpp"

/***
 * V-centric (Forward)
 * Grammar Driven
 * Separate OLD, NEW, and FUTURE lists for edges
 * Vector-copying
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
			// cout << "Unknown edge label!\n";
			// exit(0);
		}
		newEdge.label = grammar.hashSym[label];

		edges.push_back(newEdge);

		num_edges++;
		num_nodes = max(num_nodes, max(newEdge.from + 1, newEdge.to + 1));
		// SF:: If we insert the nodes into the set, the duplicate nodes will automatically get
		// discarded by the set.
		nodes.insert(newEdge.from);
		nodes.insert(newEdge.to);
	}

	infile.close();

	std::cout << "# Vertex Count:\t" << num_nodes << std::endl;
	std::cout << "# Initial Edge Count:\t" << num_edges << std::endl;
	std::cout << "Start initializing the lists, hashset and worklists ..." << std::endl;

	// level-1: vertex ID, level-2: NEW, OLD, FUTURE, level-3: outgoing edges
	vector<vector<vector<OutEdge>>> edgeVecs(num_nodes, vector<vector<OutEdge>>(3));
	// check if an edge exist or not
	unordered_set<ull> *hashset = new unordered_set<ull>[num_nodes];

	for (uint i = 0; i < num_edges; i++)
	{
		edgeVecs[edges[i].from][NEW].push_back(OutEdge(edges[i].to, edges[i].label));
	}

	edges.clear();

	for (uint i = 0; i < num_nodes; i++)
	{
		for (uint j = 0; j < edgeVecs[i][NEW].size(); j++)
		{
			hashset[i].insert(COMBINE(edgeVecs[i][NEW][j].end, edgeVecs[i][NEW][j].label));
		}
	}

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
				// insert into the graph
				edgeVecs[i][NEW].push_back(OutEdge(i, grammar.grammar1[l][0]));
			}
		}
	}

	do
	{
		itr++;
		std::cout << "Iteration number " << itr << std::endl;

		for (uint i = 0; i < num_nodes; i++)
		{
			for (uint j = 0; j < edgeVecs[i][NEW].size(); j++)
			{
				OutEdge nbr;
				nbr = edgeVecs[i][NEW][j];

				for (uint g = 0; g < grammar.grammar2.size(); g++)
				{
					if (nbr.label == grammar.grammar2[g][1])
					{
						if (hashset[i].find(COMBINE(nbr.end, grammar.grammar2[g][0])) == hashset[i].end())
						{
							hashset[i].insert(COMBINE(nbr.end, grammar.grammar2[g][0]));
							edgeVecs[i][FUTURE].push_back(OutEdge(nbr.end, grammar.grammar2[g][0]));
						}
					}
				}
			}
		}

		// for each grammar like A --> BC
		for (uint i = 0; i < num_nodes; i++)
		{
			OutEdge nbr;
			OutEdge nbrOfNbr;

			// new x new , new x old
			// for each new edge to a neighbor
			for (uint j = 0; j < edgeVecs[i][NEW].size(); j++)
			{
				nbr = edgeVecs[i][NEW][j];

				for (uint k = 0; k < edgeVecs[nbr.end][NEW].size(); k++)
				{
					nbrOfNbr = edgeVecs[nbr.end][NEW][k];

					for (uint g = 0; g < grammar.grammar3.size(); g++)
					{
						// if the neighbor's label is B
						if (nbr.label == grammar.grammar3[g][1] && nbrOfNbr.label == grammar.grammar3[g][2])
						{
							// new x new
							// for each new neighbor of this neighbor
							// if this neighbor's neighbor is labeled with C
							if (hashset[i].find(COMBINE(nbrOfNbr.end, grammar.grammar3[g][0])) == hashset[i].end())
							{
								hashset[i].insert(COMBINE(nbrOfNbr.end, grammar.grammar3[g][0]));
								edgeVecs[i][FUTURE].push_back(OutEdge(nbrOfNbr.end, grammar.grammar3[g][0]));
							}
						}
					}
				}
				// new x old
				// for each old neighbor of this neighbor
				for (uint k = 0; k < edgeVecs[nbr.end][OLD].size(); k++)
				{

					nbrOfNbr = edgeVecs[nbr.end][OLD][k];
					for (uint g = 0; g < grammar.grammar3.size(); g++)
					{
						if (nbr.label == grammar.grammar3[g][1] && nbrOfNbr.label == grammar.grammar3[g][2])
						{
							if (hashset[i].find(COMBINE(nbrOfNbr.end, grammar.grammar3[g][0])) == hashset[i].end())
							{
								hashset[i].insert(COMBINE(nbrOfNbr.end, grammar.grammar3[g][0]));
								edgeVecs[i][FUTURE].push_back(OutEdge(nbrOfNbr.end, grammar.grammar3[g][0]));
							}
						}
					}
				}
			}

			// old x new
			for (int j = 0; j < edgeVecs[i][OLD].size(); j++)
			{
				nbr = edgeVecs[i][OLD][j];

				for (int k = 0; k < edgeVecs[nbr.end][NEW].size(); k++)
				{
					nbrOfNbr = edgeVecs[nbr.end][NEW][k];
					for (uint g = 0; g < grammar.grammar3.size(); g++)
					{
						if (nbr.label == grammar.grammar3[g][1] && nbrOfNbr.label == grammar.grammar3[g][2])
						{
							if (hashset[i].find(COMBINE(nbrOfNbr.end, grammar.grammar3[g][0])) == hashset[i].end())
							{
								hashset[i].insert(COMBINE(nbrOfNbr.end, grammar.grammar3[g][0]));
								edgeVecs[i][FUTURE].push_back(OutEdge(nbrOfNbr.end, grammar.grammar3[g][0]));
							}
						}
					}
				}
			}
		}

		finished = true;

		// add NEW to OLD, then copy FUTURE to NEW
		for (int i = 0; i < num_nodes; i++)
		{
			edgeVecs[i][OLD].insert(edgeVecs[i][OLD].begin(), edgeVecs[i][NEW].begin(), edgeVecs[i][NEW].end());
			// edgeVecs[i][OLD] = edgeVecs[i][OLD_UPD];
			edgeVecs[i][NEW] = edgeVecs[i][FUTURE];
			if (edgeVecs[i][NEW].size() > 0)
				finished = false;
			edgeVecs[i][FUTURE].clear();
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
