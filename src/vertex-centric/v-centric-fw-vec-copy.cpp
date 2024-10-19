#include "globals.hpp"
#include "grammar.hpp"

/***
 * Basic vertex-centric forward direction
 *
 */

int main(int argc, char **argv)
{

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

	// level-1: vertex ID, level-2: NEW, OLD, FUTURE, level-3: outgoing edges
	vector<vector<vector<OutEdge>>> edgeVecs(num_nodes, vector<vector<OutEdge>>(3));
	// check if an edge exist or not
	unordered_set<ull> *hashset = new unordered_set<ull>[num_nodes];

	cout << "#nodes " << num_nodes << endl;
	cout << "SF::#nodes " << nodes.size() << endl;
	cout << "#edges " << num_edges << endl;

	cout << "Start making sets and the hash ...\n";
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

	cout << "Done!\n";

	// currently exclude the initialization time
	std::chrono::time_point<std::chrono::system_clock> start, finish;
	start = std::chrono::system_clock::now();

	atomic<int> newEdgeCounter(0);
	bool finished;
	int itr = 0;
	ull calcCnt = 0;

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
				// insert into the graph
				edgeVecs[i][NEW].push_back(OutEdge(i, grammar.grammar1[l][0]));
				newEdgeCounter++;
			}
		}
	}

	cout << "********************\n";

	do
	{
		itr++;

		std::chrono::time_point<std::chrono::system_clock> startC, finishC;
		startC = std::chrono::system_clock::now();

		for (uint i = 0; i < num_nodes; i++)
		{
			for (uint j = 0; j < edgeVecs[i][NEW].size(); j++)
			{
				OutEdge nbr;
				nbr = edgeVecs[i][NEW][j];

				for (uint g = 0; g < grammar.grammar2.size(); g++)
				{

					calcCnt++;
					if (nbr.label == grammar.grammar2[g][1])
					{
						if (hashset[i].find(COMBINE(nbr.end, grammar.grammar2[g][0])) == hashset[i].end())
						{
							hashset[i].insert(COMBINE(nbr.end, grammar.grammar2[g][0]));
							edgeVecs[i][FUTURE].push_back(OutEdge(nbr.end, grammar.grammar2[g][0]));
							newEdgeCounter++;
						}
					}
				}

				// edgeVecs[i][OLD_UPD].push_back(nbr);
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
						calcCnt++;
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
								newEdgeCounter++;
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
						calcCnt++;
						if (nbr.label == grammar.grammar3[g][1] && nbrOfNbr.label == grammar.grammar3[g][2])
						{
							if (hashset[i].find(COMBINE(nbrOfNbr.end, grammar.grammar3[g][0])) == hashset[i].end())
							{
								hashset[i].insert(COMBINE(nbrOfNbr.end, grammar.grammar3[g][0]));
								edgeVecs[i][FUTURE].push_back(OutEdge(nbrOfNbr.end, grammar.grammar3[g][0]));
								newEdgeCounter++;
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
						calcCnt++;
						if (nbr.label == grammar.grammar3[g][1] && nbrOfNbr.label == grammar.grammar3[g][2])
						{
							if (hashset[i].find(COMBINE(nbrOfNbr.end, grammar.grammar3[g][0])) == hashset[i].end())
							{
								hashset[i].insert(COMBINE(nbrOfNbr.end, grammar.grammar3[g][0]));
								edgeVecs[i][FUTURE].push_back(OutEdge(nbrOfNbr.end, grammar.grammar3[g][0]));
								newEdgeCounter++;
							}
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

		if (debug)
		{
			cout << "This Iteration Computation Time = " << elapsed_secondsC.count() << endl;
		}

		std::chrono::time_point<std::chrono::system_clock> startT, finishT;
		startT = std::chrono::system_clock::now();

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

		finishT = std::chrono::system_clock::now();
		std::chrono::duration<double> elapsed_secondsT = finishT - startT;
		std::time_t finish_timeT = std::chrono::system_clock::to_time_t(finishT);

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

	cout << "SF:: # Number of new edges: " << totalNewEdgeCount << endl;
	cout << "AM:: # Number of new edges: " << newEdgeCounter << endl;

	cout << "# Iterations: " << itr << endl;
	cout << "# Total Calculations: " << calcCnt << endl;
}
