#include <cilk/cilk.h>
#include <cilk/cilk_api.h>

#include "globals.hpp"
#include "grammar.hpp"

/***
 * V-centric (bi, temporal ptrs, rule-lookup table)
 * Bi directional traversing of the Graph (incoming and outgoing edges)
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
	vector<Buffer> inEdgeVecs(num_nodes);
	// level-1: vertex ID, level-2: NEW, OLD, FUTURE pointers and outgoing edges
	vector<Buffer> edgeVecs(num_nodes);

	// check if an edge exist or not
	unordered_set<ull> *hashset = new unordered_set<ull>[num_nodes];
	// hashset for incoming edges
	unordered_set<ull> *inHashset = new unordered_set<ull>[num_nodes];

	cout << "#nodes " << num_nodes << endl;
	cout << "SF::#nodes " << nodes.size() << endl;
	cout << "#edges " << num_edges << endl;

	cout << "Start making sets and the hash ...\n";
	for (uint i = 0; i < num_edges; i++)
	{
		edgeVecs[edges[i].from].vertexList.push_back(OutEdge(edges[i].to, edges[i].label));
		inEdgeVecs[edges[i].to].vertexList.push_back(OutEdge(edges[i].from, edges[i].label));

		// update the buffer pointers
		edgeVecs[edges[i].from].NEW_END++;
		inEdgeVecs[edges[i].to].NEW_END++;

		hashset[edges[i].from].insert(COMBINE(edges[i].to, edges[i].label));
		inHashset[edges[i].to].insert(COMBINE(edges[i].from, edges[i].label));
	}

	edges.clear();

	// count the total no of initial unique edge
	uint initialEdgeCount = countEdge(hashset, num_nodes);

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

	if (debug)
	{
		cout << "[DEBUG]: grammar1 size: " << grammar.grammar1.size() << endl;
	}
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
			if (hashset[i].find(COMBINE(i, grammar.grammar1[l][0])) == hashset[i].end())
			{
				// insert into hashset
				hashset[i].insert(COMBINE(i, grammar.grammar1[l][0]));
				inHashset[i].insert(COMBINE(i, grammar.grammar1[l][0]));
				// insert into the graph
				edgeVecs[i].vertexList.push_back(OutEdge(i, grammar.grammar1[l][0]));
				inEdgeVecs[i].vertexList.push_back(OutEdge(i, grammar.grammar1[l][0]));

				// update the pointers
				edgeVecs[i].NEW_END++;
				inEdgeVecs[i].NEW_END++;

				newEdgeCounter++;
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

			// for each vertex
			// cilk_for(uint i=0; i<num_nodes; i++)

			OutEdge nbr;
			// uint START_NEW = pointerEdgeVecs[i][OLD_CNT];
			// uint END_NEW = START_NEW + pointerEdgeVecs[i][NEW_CNT];
			//  the valid index range is [START_NEW, END_NEW-1]
			uint START_NEW = edgeVecs[i].OLD_END;
			uint END_NEW = edgeVecs[i].NEW_END;

			// for each new edge
			for (uint j = START_NEW; j < END_NEW; j++)
			{
				nbr = edgeVecs[i].vertexList[j];
				// if the edge to the neighbor is labeled with B

				vector<uint> leftLabels = grammar.rule2(nbr.label);

				for (uint g = 0; g < leftLabels.size(); g++)
				{
					calcCnt++;
					if (hashset[i].find(COMBINE(nbr.end, leftLabels[g])) == hashset[i].end())
					{
						finished = false;
						hashset[i].insert(COMBINE(nbr.end, leftLabels[g]));
						inHashset[nbr.end].insert(COMBINE(i, leftLabels[g]));

						edgeVecs[i].vertexList.push_back(OutEdge(nbr.end, leftLabels[g]));
						inEdgeVecs[nbr.end].vertexList.push_back(OutEdge(i, leftLabels[g]));

						// No need to update the pointers. Because the FUTURE_START starts from the
						// NEW_END, and NEW_END is already updated.

						newEdgeCounter++;
					}
				}

				// all the OLD, and NEW outgoing edges of the first edge
				uint START_OLD_OUT = 0;
				uint END_NEW_OUT = edgeVecs[nbr.end].NEW_END;
				// grammar3 calc
				for (uint h = START_OLD_OUT; h < END_NEW_OUT; h++)
				{
					OutEdge outNbr = edgeVecs[nbr.end].vertexList[h];
					vector<uint> leftLabels = grammar.rule3(nbr.label, outNbr.label);

					for (uint g = 0; g < leftLabels.size(); g++)
					{
						calcCnt++;
						if (hashset[i].find(COMBINE(outNbr.end, leftLabels[g])) == hashset[i].end())
						{
							finished = false;
							hashset[i].insert(COMBINE(outNbr.end, leftLabels[g]));
							inHashset[outNbr.end].insert(COMBINE(i, leftLabels[g]));
							edgeVecs[i].vertexList.push_back(OutEdge(outNbr.end, leftLabels[g]));
							inEdgeVecs[outNbr.end].vertexList.push_back(OutEdge(i, leftLabels[g]));
							// No need to update the pointers. Because the FUTURE_START starts from the
							// NEW_END, and NEW_END is already updated.
						}
					}
				}

				// OLD incoming neighbors of the first edge
				uint START_OLD_IN = 0;
				uint END_NEW_IN = inEdgeVecs[i].OLD_END;

				for (uint h = START_OLD_IN; h < END_NEW_IN; h++)
				{
					OutEdge inNbr = inEdgeVecs[i].vertexList[h];
					vector<uint> leftLabels = grammar.rule3(inNbr.label, nbr.label);

					for (uint g = 0; g < leftLabels.size(); g++)
					{
						calcCnt++;
						if (hashset[inNbr.end].find(COMBINE(nbr.end, leftLabels[g])) == hashset[inNbr.end].end())
						{
							finished = false;
							hashset[inNbr.end].insert(COMBINE(nbr.end, leftLabels[g]));
							inHashset[nbr.end].insert(COMBINE(inNbr.end, leftLabels[g]));
							edgeVecs[inNbr.end].vertexList.push_back(OutEdge(nbr.end, leftLabels[g]));
							inEdgeVecs[nbr.end].vertexList.push_back(OutEdge(inNbr.end, leftLabels[g]));
							// No need to update the pointers. Because the FUTURE_START starts from the
							// NEW_END, and NEW_END is already updated.
						}
					}
				}
			}
		}

		cout << "Iteration number " << itr << endl;

		if (debug)
		{
			cout << "Number of new edges so far: " << newEdgeCounter << endl;
		}

		finishC = std::chrono::system_clock::now();
		std::chrono::duration<double> elapsed_secondsC = finishC - startC;
		std::time_t finish_timeC = std::chrono::system_clock::to_time_t(finishC);
		elapsed_seconds_comp += elapsed_secondsC.count();

		if (debug)
		{
			cout << "This Iteration Computation Time = " << elapsed_secondsC.count() << endl;
		}

		std::chrono::time_point<std::chrono::system_clock> startT, finishT;
		startT = std::chrono::system_clock::now();

		for (int i = 0; i < num_nodes; i++)
		{
			// update the pointers
			edgeVecs[i].OLD_END = edgeVecs[i].NEW_END;
			inEdgeVecs[i].OLD_END = inEdgeVecs[i].NEW_END;

			edgeVecs[i].NEW_END = edgeVecs[i].vertexList.size();
			inEdgeVecs[i].NEW_END = inEdgeVecs[i].vertexList.size();
		}

		finishT = std::chrono::system_clock::now();
		std::chrono::duration<double> elapsed_secondsT = finishT - startT;
		std::time_t finish_timeT = std::chrono::system_clock::to_time_t(finishT);
		elapsed_seconds_transfer += elapsed_secondsT.count();

		if (debug)
		{
			std::cout << "This Iteration Transformation Time = " << elapsed_secondsT.count() << std::endl;
		}

	} while (!finished);

	finish = std::chrono::system_clock::now();
	std::chrono::duration<double> elapsed_seconds = finish - start;
	std::time_t finish_time = std::chrono::system_clock::to_time_t(finish);

	uint totalNewEdgeCount = countEdge(hashset, num_nodes) - initialEdgeCount;

	cout << "**************************" << endl;
	std::cout << "# Total time = " << elapsed_seconds.count() << std::endl;
	std::cout << "# Total computation time = " << elapsed_seconds_comp << std::endl;
	std::cout << "# Total data transfer time = " << elapsed_seconds_transfer << std::endl;

	cout << "SF:: # Number of new edges: " << totalNewEdgeCount << endl;
	cout << "AM:: # Number of new edges: " << newEdgeCounter << endl;

	cout << "# Iterations: " << itr << endl;
	cout << "# Total Calculations: " << calcCnt << endl;
}
