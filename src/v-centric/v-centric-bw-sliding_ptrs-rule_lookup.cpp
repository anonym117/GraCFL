// #include <cilk/cilk.h>
// #include <cilk/cilk_api.h>

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

	// level-1: vertex ID, level-2: NEW, OLD, FUTURE pointers and incoming edges
	//vector<Buffer> inEdgeVecs(num_nodes);
	// level-1: vertex ID, level-2: NEW, OLD, FUTURE pointers and outgoing edges
	vector<Buffer> inEdgeVecs(num_nodes);

	// check if an edge exist or not
	// unordered_set<ull> *hashset = new unordered_set<ull>[num_nodes];
	vector<vector<unordered_set<ull>>> inHashset(num_nodes, vector<unordered_set<ull>>(grammar.labelSize, unordered_set<ull>()));

	cout << "#nodes " << num_nodes << endl;
	cout << "SF::#nodes " << nodes.size() << endl;
	cout << "#edges " << num_edges << endl;

	cout << "Start making sets and the hash ...\n";
	for (uint i = 0; i < num_edges; i++)
	{
		//edgeVecs[edges[i].from].vertexList.push_back(OutEdge(edges[i].to, edges[i].label));
		inEdgeVecs[edges[i].to].vertexList.push_back(OutEdge(edges[i].from, edges[i].label));

		// update the temporal pointers
		//edgeVecs[edges[i].from].NEW_END++;
		inEdgeVecs[edges[i].to].NEW_END++;

		inHashset[edges[i].to][edges[i].label].insert(edges[i].from);
	}

	edges.clear();

	// count the total no of initial unique edge
	uint initialEdgeCount = countEdge(inHashset, num_nodes, grammar.labelSize);

	cout << "Done!\n";

	// currently exclude the initialization time
	std::chrono::time_point<std::chrono::system_clock> start, finish;
	start = std::chrono::system_clock::now();

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
			if (inHashset[i][grammar.grammar1[l][0]].find(i) == inHashset[i][grammar.grammar1[l][0]].end())
			{
				// insert into hashset
				inHashset[i][grammar.grammar1[l][0]].insert(i);
				// insert into the graph
				//edgeVecs[i].vertexList.push_back(OutEdge(i, grammar.grammar1[l][0]));
				inEdgeVecs[i].vertexList.push_back(OutEdge(i, grammar.grammar1[l][0]));
				// update the pointers
				//edgeVecs[i].NEW_END++;
				inEdgeVecs[i].NEW_END++;
			}
		}
	}

	cout << "********************\n";

	do
	{
		finished = true;
		itr++;

		std::chrono::time_point<std::chrono::system_clock> startC, finishC;
		startC = std::chrono::system_clock::now();

		// for each grammar rule like A --> B
		for (uint i = 0; i < num_nodes; i++)
		{
			OutEdge nbr;
			//  the valid index range is [START_NEW, END_NEW-1]
			uint START_NEW_OUT = inEdgeVecs[i].OLD_END;
			uint END_NEW_OUT = inEdgeVecs[i].NEW_END;

			// for each new edge
			for (uint j = START_NEW_OUT; j < END_NEW_OUT; j++)
			{
				nbr = inEdgeVecs[i].vertexList[j];
				// if the edge to the neighbor is labeled with B
				vector<uint> leftLabels = grammar.rule2(nbr.label);

				for (uint g = 0; g < leftLabels.size(); g++)
				{
					calcCnt++;
					if (inHashset[i][leftLabels[g]].find(nbr.end) == inHashset[i][leftLabels[g]].end())
					{
						finished = false;
						//hashset[][leftLabels[g]].insert(nbr.end);
						inHashset[i][leftLabels[g]].insert(nbr.end);

						//edgeVecs[i].vertexList.push_back(OutEdge(nbr.end, leftLabels[g]));
						inEdgeVecs[i].vertexList.push_back(OutEdge(nbr.end, leftLabels[g]));
					}
				}

				// NEW + OLD  neighbors of the first edge
				uint START_OLD = 0;
				uint END_NEW = inEdgeVecs[nbr.end].NEW_END;

				for (uint h = START_OLD; h < END_NEW; h++)
				{
					OutEdge outInNbr = inEdgeVecs[nbr.end].vertexList[h];
					vector<uint> leftLabels = grammar.rule3(outInNbr.label, nbr.label);

					for (uint g = 0; g < leftLabels.size(); g++)
					{
						calcCnt++;
						if (inHashset[i][leftLabels[g]].find(outInNbr.end) == inHashset[i][leftLabels[g]].end())
						{
							finished = false;
							inHashset[i][leftLabels[g]].insert(outInNbr.end);
							//edgeVecs[i].vertexList.push_back(OutEdge(outNbr.end, leftLabels[g]));
							inEdgeVecs[i].vertexList.push_back(OutEdge(outInNbr.end, leftLabels[g]));
						}
					}
				}
			}

			uint OLD_START_IN = 0;
			uint OLD_END_IN = inEdgeVecs[i].OLD_END;

			for (uint j = OLD_START_IN; j < OLD_END_IN; j++)
			{
				OutEdge nbr = inEdgeVecs[i].vertexList[j];

				uint NEW_START_IN = inEdgeVecs[nbr.end].OLD_END;
				uint NEW_END_IN = inEdgeVecs[nbr.end].NEW_END;

				for (uint h = NEW_START_IN; h < NEW_END_IN; h++)
				{
					OutEdge outInNbr = inEdgeVecs[nbr.end].vertexList[h];
					vector<uint> leftLabels = grammar.rule3(outInNbr.label, nbr.label);

					for (uint g = 0; g < leftLabels.size(); g++)
					{
						calcCnt++;
						if (inHashset[i][leftLabels[g]].find(outInNbr.end) == inHashset[i][leftLabels[g]].end())
						{
							finished = false;
							inHashset[i][leftLabels[g]].insert(outInNbr.end);
							//edgeVecs[i].vertexList.push_back(OutEdge(outNbr.end, leftLabels[g]));
							inEdgeVecs[i].vertexList.push_back(OutEdge(outInNbr.end, leftLabels[g]));
						}
					}
				}
			}
		}

		cout << "Iteration number " << itr << endl;

		for (int i = 0; i < num_nodes; i++)
		{
			// update the pointers
			inEdgeVecs[i].OLD_END = inEdgeVecs[i].NEW_END;
			//inEdgeVecs[i].OLD_END = inEdgeVecs[i].NEW_END;

			inEdgeVecs[i].NEW_END = inEdgeVecs[i].vertexList.size();
			//inEdgeVecs[i].NEW_END = inEdgeVecs[i].vertexList.size();
		}
	} while (!finished);

	finish = std::chrono::system_clock::now();
	std::chrono::duration<double> elapsed_seconds = finish - start;
	std::time_t finish_time = std::chrono::system_clock::to_time_t(finish);

	uint totalNewEdgeCount = countEdge(inHashset, num_nodes, grammar.labelSize) - initialEdgeCount;

	cout << "**************************" << endl;
	std::cout << "# Total time = " << elapsed_seconds.count() << std::endl;
	std::cout << "# Total computation time = " << elapsed_seconds_comp << std::endl;
	std::cout << "# Total data transfer time = " << elapsed_seconds_transfer << std::endl;

	cout << "SF:: # Number of new edges: " << totalNewEdgeCount << endl;
	cout << "AM:: # Number of new edges: " << newEdgeCounter << endl;

	cout << "# Iterations: " << itr << endl;
	cout << "# Total Calculations: " << calcCnt << endl;

    getPeakMemoryUsage();
}
