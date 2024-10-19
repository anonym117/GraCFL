#include "globals.hpp"
#include "grammar.hpp"

/***
 * Bi directional traversing of the Graph (incoming and outgoing edges)
 * Three separate buffer lists for OLD, NEW, and FUTURE edges
 */

int main(int argc, char **argv)
{

	bool debug = true;
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

	// *** read grammer
	// *** read graph
	const string inputGraph = argv[1];
	Grammar grammar(argv[2]);

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

	vector<vector<vector<OutEdge>>> inEdgeVecs(num_nodes, vector<vector<OutEdge>>(3));
	// level-1: vertex ID, level-2: NEW, OLD, FUTURE, level-3: outgoing edges
	vector<vector<vector<OutEdge>>> edgeVecs(num_nodes, vector<vector<OutEdge>>(3));
	// check if an edge exist or not
	unordered_set<ull> *hashset = new unordered_set<ull>[num_nodes];
	// hashset for incoming edges
	// unordered_set<ull> *inHashset = new unordered_set<ull>[num_nodes];

	cout << "#nodes " << num_nodes << endl;
	cout << "SF::#nodes " << nodes.size() << endl;
	cout << "#edges " << num_edges << endl;

	cout << "Start making sets and the hash ...\n";
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
				// inHashset[i].insert(COMBINE(i, grammar.grammar1[l][0]));
				// insert into the graph
				edgeVecs[i][NEW].push_back(OutEdge(i, grammar.grammar1[l][0]));
				inEdgeVecs[i][NEW].push_back(OutEdge(i, grammar.grammar1[l][0]));

				// edgeVecs[i][OLD].push_back(OutEdge(i, grammar.grammar1[l][0]));
				// inEdgeVecs[i][OLD].push_back(OutEdge(i, grammar.grammar1[l][0]));

				newEdgeCounter++;
			}
		}
	}

	cout << "********************\n";

	do
	{
		itr++;

		// std::chrono::time_point<std::chrono::system_clock> startC, finishC;
		// startC = std::chrono::system_clock::now();

		// for each grammar rule like A --> B
		for (uint i = 0; i < num_nodes; i++)
		{
			for (uint j = 0; j < inEdgeVecs[i][NEW].size(); j++)
			{
				OutEdge inNbr = inEdgeVecs[i][NEW][j];

				for (uint g = 0; g < grammar.grammar2.size(); g++)
				{
					calcCnt++;
					if (inNbr.label == grammar.grammar2[g][1])
					{
						if (hashset[inNbr.end].find(COMBINE(i, grammar.grammar2[g][0])) == hashset[inNbr.end].end())
						{
							hashset[inNbr.end].insert(COMBINE(i, grammar.grammar2[g][0]));
							// inHashset[nbr.end].insert(COMBINE(i, grammar.grammar2[g][0]));

							edgeVecs[inNbr.end][FUTURE].push_back(OutEdge(i, grammar.grammar2[g][0]));
							inEdgeVecs[i][FUTURE].push_back(OutEdge(inNbr.end, grammar.grammar2[g][0]));

							newEdgeCounter++;
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
						calcCnt++;
						if (inNbr.label == grammar.grammar3[g][1] && outNbr.label == grammar.grammar3[g][2])
						{
							if (hashset[inNbr.end].find(COMBINE(outNbr.end, grammar.grammar3[g][0])) == hashset[inNbr.end].end())
							{
								hashset[inNbr.end].insert(COMBINE(outNbr.end, grammar.grammar3[g][0]));
								// inHashset[outNbr.end].insert(COMBINE(i, grammar.grammar3[g][0]));
								edgeVecs[inNbr.end][FUTURE].push_back(OutEdge(outNbr.end, grammar.grammar3[g][0]));
								inEdgeVecs[outNbr.end][FUTURE].push_back(OutEdge(inNbr.end, grammar.grammar3[g][0]));

								newEdgeCounter++;
							}
						}
					}
				}

				for (uint h = 0; h < edgeVecs[i][OLD].size(); h++)
				{
					OutEdge outNbr = edgeVecs[i][OLD][h];
					for (uint g = 0; g < grammar.grammar3.size(); g++)
					{
						calcCnt++;
						if (inNbr.label == grammar.grammar3[g][1] && outNbr.label == grammar.grammar3[g][2])
						{
							if (hashset[inNbr.end].find(COMBINE(outNbr.end, grammar.grammar3[g][0])) == hashset[inNbr.end].end())
							{
								hashset[inNbr.end].insert(COMBINE(outNbr.end, grammar.grammar3[g][0]));
								// inHashset[outNbr.end].insert(COMBINE(i, grammar.grammar3[g][0]));
								edgeVecs[inNbr.end][FUTURE].push_back(OutEdge(outNbr.end, grammar.grammar3[g][0]));
								inEdgeVecs[outNbr.end][FUTURE].push_back(OutEdge(inNbr.end, grammar.grammar3[g][0]));

								newEdgeCounter++;
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
						calcCnt++;
						if (inNbr.label == grammar.grammar3[g][1] && outNbr.label == grammar.grammar3[g][2])
						{
							if (hashset[inNbr.end].find(COMBINE(outNbr.end, grammar.grammar3[g][0])) == hashset[inNbr.end].end())
							{
								hashset[inNbr.end].insert(COMBINE(outNbr.end, grammar.grammar3[g][0]));
								// inHashset[outNbr.end].insert(COMBINE(i, grammar.grammar3[g][0]));
								edgeVecs[inNbr.end][FUTURE].push_back(OutEdge(outNbr.end, grammar.grammar3[g][0]));
								inEdgeVecs[outNbr.end][FUTURE].push_back(OutEdge(inNbr.end, grammar.grammar3[g][0]));

								newEdgeCounter++;
							}
						}
					}
				}

				// edgeVecs[i][OLD].push_back(OutEdge(nbr.end, nbr.label));
				// inEdgeVecs[nbr.end][OLD].push_back(OutEdge(i, nbr.label));
			}
		}

		cout << "Iteration number " << itr << endl;

		// if (debug)
		// {
		// 	cout << "Number of new edges so far: " << newEdgeCounter << endl;
		// }

		// finishC = std::chrono::system_clock::now();
		// std::chrono::duration<double> elapsed_secondsC = finishC - startC;
		// std::time_t finish_timeC = std::chrono::system_clock::to_time_t(finishC);
		// elapsed_seconds_comp += elapsed_secondsC.count();

		// if (debug)
		// {
		// 	cout << "This Iteration Computation Time = " << elapsed_secondsC.count() << endl;
		// }

		std::chrono::time_point<std::chrono::system_clock> startT, finishT;
		startT = std::chrono::system_clock::now();

		finished = true;

		// // add NEW to OLD, then copy FUTURE to NEW
		// for (int i = 0; i < num_nodes; i++)
		// {
		// 	// a.insert(a.end(), b.begin(), b.end());
		// 	edgeVecs[i][OLD].insert(edgeVecs[i][OLD].begin(), edgeVecs[i][NEW].begin(), edgeVecs[i][NEW].end());
		// 	// edgeVecs[i][OLD] = edgeVecs[i][OLD_UPD];
		// 	edgeVecs[i][NEW] = edgeVecs[i][FUTURE];
		// 	if (edgeVecs[i][FUTURE].size() > 0)
		// 		finished = false;
		// 	edgeVecs[i][FUTURE].clear();
		// }

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

		finishT = std::chrono::system_clock::now();
		std::chrono::duration<double> elapsed_secondsT = finishT - startT;
		// std::time_t finish_timeT = std::chrono::system_clock::to_time_t(finishT);
		elapsed_seconds_transfer += elapsed_secondsT.count();

		// if (debug)
		// {
		// 	std::cout << "This Iteration Transformation Time = " << elapsed_secondsT.count() << std::endl;
		// }

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
