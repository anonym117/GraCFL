#include "globals.hpp"
#include "grammar.hpp"

/*
 * E-centric w/ rule-idx and label-idx-adj_list
 * new edge worklist, iterative
 * supports pocr graphs
 */

int main(int argc, char **argv)
{

	// bool debug = true;
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
	std::string label;

	std::string line;
	while (getline(infile, line))
	{
		EdgeForReading newEdge;
		std::stringstream ss;
		ss << line;
		// cout << line << endl;
		std::string symbol;
		vector<std::string> symbolsVec;
		while (ss >> symbol)
		{
			symbolsVec.push_back(symbol);
		}
		uint size = symbolsVec.size();
		// cout << "size: " << size << endl;
		if (size == 0)
		{
			continue;
		}
		else if (size >= 3)
		{
			newEdge.from = stoi(symbolsVec[0]);
			newEdge.to = stoi(symbolsVec[1]);
			label = symbolsVec[2];
			if (grammar.hashSym.find(label) == grammar.hashSym.end())
			{
				// Unknown edge label
				continue;
			}
			newEdge.label = grammar.hashSym[label];

			if (size == 3)
			{
				// if there is no contextID, then put 0
				newEdge.contextID = 0;
			}
			else
			{
				// Add 1 with the contextID, as 0 means no contextID
				// and there can be contextID 0, so we will have to make 0 to 1
				// to remove the ambiguity
				newEdge.contextID = stoi(symbolsVec[3]);
				newEdge.contextID++;
			}
		}

		// cout << "newEdge:: " << newEdge.from << "	" << newEdge.to << "	" << newEdge.label << "	" << newEdge.contextID << endl;
		edges.push_back(newEdge);

		num_edges++;
		num_nodes = max(num_nodes, max(newEdge.from + 1, newEdge.to + 1));
		// SF:: If we insert the nodes into the set, the duplicate nodes will automatically get
		// discarded by the set.
		nodes.insert(newEdge.from);
		nodes.insert(newEdge.to);
	}

	infile.close();

	vector<pair<uint, uint>> **inEdgeVecs = new vector<pair<uint, uint>> *[grammar.labelSize];
	// level-1: vertex ID, level-2: NEW, OLD, FUTURE, level-3: outgoing edges
	vector<pair<uint, uint>> **edgeVecs = new vector<pair<uint, uint>> *[grammar.labelSize];

	// check if an edge exist or not
	// unordered_set<ull> **hashset = new unordered_set<ull> *[num_nodes];
	vector<vector<unordered_set<ull>>> hashset(num_nodes, vector<unordered_set<ull>>(grammar.labelSize, unordered_set<ull>()));
	// hashset for incoming edges
	// unordered_set<ull> *inHashset = new unordered_set<ull>[num_nodes];

	queue<EdgeForReading> activeQueue;
	queue<EdgeForReading> futureQueue;

	cout << "#nodes " << num_nodes << endl;
	cout << "SF::#nodes " << nodes.size() << endl;
	cout << "#edges " << num_edges << endl;

	cout << "Start making sets and the hash ...\n";

	for (uint i = 0; i < grammar.labelSize; i++)
	{
		edgeVecs[i] = new vector<pair<uint, uint>>[num_nodes];
		inEdgeVecs[i] = new vector<pair<uint, uint>>[num_nodes];
	}

	// for (uint i = 0; i < num_nodes; i++)
	// {
	// 	hashset[i] = new unordered_set<ull>[grammar.labelSize];
	// }

	for (uint i = 0; i < num_edges; i++)
	{
		edgeVecs[edges[i].label][edges[i].from].push_back(make_pair(edges[i].to, edges[i].contextID));
		inEdgeVecs[edges[i].label][edges[i].to].push_back(make_pair(edges[i].from, edges[i].contextID));

		hashset[edges[i].from][edges[i].label].insert(COMBINE(edges[i].to, edges[i].contextID));
		// inHashset[edges[i].to].insert(COMBINE(edges[i].from, edges[i].label));

		activeQueue.push(edges[i]);
	}

	edges.clear();

	// count the total no of initial unique edge
	uint initialEdgeCount = countEdge(hashset, num_nodes, grammar.labelSize);

	cout << "Done!\n";

	// currently exclude the initialization time

	atomic<int> newEdgeCounter(0);
	bool finished;
	int itr = 0;
	ull calcCnt = 0;
	ull calcCnt1 = 0;
	ull calcCnt2 = 0;
	ull calcCnt3_left = 0;
	ull calcCnt3_right = 0;

	double elapsed_seconds_comp = 0.0;
	double elapsed_seconds_transfer = 0.0;

	std::chrono::time_point<std::chrono::system_clock> start, finish;
	start = std::chrono::system_clock::now();

	// if (debug)
	// {
	// 	cout << "[DEBUG]: grammar1 size: " << grammar.grammar1.size() << endl;
	// }
	// handle epsilon rules: add an edge to itself
	// grammar1 is for epsilon rules A --> e
	// grammar2 is for one symbol on RHS A --> B
	// grammar3 is for two symbols on RHS A --> BC
	for (uint l = 0; l < grammar.grammar1.size(); l++)
	{
		for (uint i = 0; i < num_nodes; i++)
		{
			// calcCnt++;
			// calcCnt1++;

			// check if the new edge based on an epsilon grammar rule exists or not. l: grammar ID, 0: LHS
			if (hashset[i][grammar.grammar1[l][0]].find(i) == hashset[i][grammar.grammar1[l][0]].end())
			{
				// insert into hashset
				hashset[i][grammar.grammar1[l][0]].insert(COMBINE(i, 0));
				// inHashset[i].insert(COMBINE(i, grammar.grammar1[l][0]));
				// insert into the graph
				edgeVecs[grammar.grammar1[l][0]][i].push_back(make_pair(i, 0));
				inEdgeVecs[grammar.grammar1[l][0]][i].push_back(make_pair(i, 0));

				activeQueue.push(EdgeForReading(i, i, grammar.grammar1[l][0], 0));

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

		while (!activeQueue.empty())
		{
			EdgeForReading currEdge = activeQueue.front();
			activeQueue.pop();

			// for each grammar rule like A --> B
			for (uint g = 0; g < grammar.grammar2index[currEdge.label].size(); g++)
			{
				// calcCnt++;
				// calcCnt2++;
				uint leftLabel = grammar.grammar2index[currEdge.label][g];
				// if the edge to the neighbor is labeled with B
				if (hashset[currEdge.from][leftLabel].find(COMBINE(currEdge.to, currEdge.contextID)) == hashset[currEdge.from][leftLabel].end())
				{
					hashset[currEdge.from][leftLabel].insert(COMBINE(currEdge.to, currEdge.contextID));
					// inHashset[currEdge.to].insert(COMBINE(currEdge.from, grammar.grammar2index[currEdge.label][g]));

					// the new edge is added to FUTURE list, due to potential conflicts in parallel version
					// if topo-driven is run in serial, only OLD and NEW are needed.
					// Graspan did not mention FUTURE list, but it puts the new edges into a seperate list, which
					// is similar to FUTURE list's purpose.
					// for data-driven, we only need one list (no NEW OLD FUTURE)

					futureQueue.push(EdgeForReading(currEdge.from, currEdge.to, grammar.grammar2index[currEdge.label][g], currEdge.contextID));

					// newEdgeCounter++;
				}
			}

			// SF: start
			// A = BC
			for (uint g = 0; g < grammar.grammar3indexLeft[currEdge.label].size(); g++)
			{
				uint B = currEdge.label;
				uint A = grammar.grammar3indexLeft[currEdge.label][g].second;
				uint C = grammar.grammar3indexLeft[currEdge.label][g].first;

				for (uint j = 0; j < edgeVecs[C][currEdge.to].size(); j++)
				{
					// uint nbr;

					// A = BC, B = grammar.grammar3[g][1] = currEdge.label
					// NEW outgoing edges of currEdge.to
					pair<uint, uint> nbr = edgeVecs[C][currEdge.to][j];
					// A = BC, C = grammar.grammar3[g][2]
					// calcCnt++;
					// calcCnt3_left++;

					if (currEdge.contextID == nbr.second)
					{
						uint cxId = 0;
						if (grammar.contextLabels[A])
						{
							cxId = currEdge.contextID;
						}

						if (hashset[currEdge.from][A].find(COMBINE(nbr.first, cxId)) == hashset[currEdge.from][A].end())
						{
							hashset[currEdge.from][A].insert(COMBINE(nbr.first, cxId));
							// inHashset[nbr].insert(COMBINE(currEdge.from, A));

							futureQueue.push(EdgeForReading(currEdge.from, nbr.first, A, cxId));
							// newEdgeCounter++;
						}
					}
					else
					{
						if (currEdge.contextID == 0)
						{
							uint cxId = nbr.second;
							if (hashset[currEdge.from][A].find(COMBINE(nbr.first, cxId)) == hashset[currEdge.from][A].end())
							{
								hashset[currEdge.from][A].insert(COMBINE(nbr.first, cxId));
								// inHashset[nbr].insert(COMBINE(currEdge.from, A));

								futureQueue.push(EdgeForReading(currEdge.from, nbr.first, A, cxId));
								// newEdgeCounter++;
							}
						}
						else if (nbr.second == 0)
						{
							uint cxId = currEdge.contextID;
							if (hashset[currEdge.from][A].find(COMBINE(nbr.first, cxId)) == hashset[currEdge.from][A].end())
							{
								hashset[currEdge.from][A].insert(COMBINE(nbr.first, cxId));
								// inHashset[nbr].insert(COMBINE(currEdge.from, A));

								futureQueue.push(EdgeForReading(currEdge.from, nbr.first, A, cxId));
								// newEdgeCounter++;
							}
						}
					}
				}
			}

			// A = CB, B = grammar.grammar3[g][2] = currEdge.label
			for (uint g = 0; g < grammar.grammar3indexRight[currEdge.label].size(); g++)
			{
				uint B = currEdge.label;
				uint A = grammar.grammar3indexRight[currEdge.label][g].second;
				uint C = grammar.grammar3indexRight[currEdge.label][g].first;

				for (uint h = 0; h < inEdgeVecs[C][currEdge.from].size(); h++)
				{

					pair<uint, uint> inNbr = inEdgeVecs[C][currEdge.from][h];
					// A = CB, C = grammar.grammar3[g][1]
					// calcCnt++;
					// calcCnt3_right++;

					if (currEdge.contextID == inNbr.second)
					{
						uint cxId = 0;
						if (grammar.contextLabels[A])
						{
							cxId = currEdge.contextID;
						}

						if (hashset[inNbr.first][A].find(COMBINE(currEdge.to, cxId)) == hashset[inNbr.first][A].end())
						{
							hashset[inNbr.first][A].insert(COMBINE(currEdge.to, cxId));
							// inHashset[currEdge.to].insert(COMBINE(inNbr, A));

							futureQueue.push(EdgeForReading(inNbr.first, currEdge.to, A, cxId));
							// newEdgeCounter++;
						}
					}
					else
					{
						if (currEdge.contextID == 0)
						{
							uint cxId = inNbr.second;

							if (hashset[inNbr.first][A].find(COMBINE(currEdge.to, cxId)) == hashset[inNbr.first][A].end())
							{
								hashset[inNbr.first][A].insert(COMBINE(currEdge.to, cxId));
								// inHashset[currEdge.to].insert(COMBINE(inNbr, A));

								futureQueue.push(EdgeForReading(inNbr.first, currEdge.to, A, cxId));
								// newEdgeCounter++;
							}
						}

						if (inNbr.second == 0)
						{
							uint cxId = currEdge.contextID;

							if (hashset[inNbr.first][A].find(COMBINE(currEdge.to, cxId)) == hashset[inNbr.first][A].end())
							{
								hashset[inNbr.first][A].insert(COMBINE(currEdge.to, cxId));
								// inHashset[currEdge.to].insert(COMBINE(inNbr, A));

								futureQueue.push(EdgeForReading(inNbr.first, currEdge.to, A, cxId));
								// newEdgeCounter++;
							}
						}
					}
				}
			}
		}

		cout << "Iteration number " << itr << endl;

		finished = true;

		queue<EdgeForReading> tempQueue = futureQueue;

		while (!tempQueue.empty())
		{
			EdgeForReading edge = tempQueue.front();
			tempQueue.pop();

			edgeVecs[edge.label][edge.from].push_back(make_pair(edge.to, edge.contextID));
			inEdgeVecs[edge.label][edge.to].push_back(make_pair(edge.from, edge.contextID));
		}

		if (!futureQueue.empty())
			finished = false;
		swap(activeQueue, futureQueue);

	} while (!finished);

	finish = std::chrono::system_clock::now();
	std::chrono::duration<double> elapsed_seconds = finish - start;
	std::time_t finish_time = std::chrono::system_clock::to_time_t(finish);

	uint totalNewEdgeCount = countEdge(hashset, num_nodes, grammar.labelSize) - initialEdgeCount;

	std::cout << "**************************" << endl;
	std::cout << "# Total time = " << elapsed_seconds.count() << std::endl;
	std::cout << "# Total iterations = " << itr << endl;
	// std::cout << "# Total computation time = " << elapsed_seconds_comp << std::endl;
	// std::cout << "# Total data transfer time = " << elapsed_seconds_transfer << std::endl;

	std::cout << "SF:: # Number of new edges: " << totalNewEdgeCount << endl;
	std::cout << "AM:: # Number of new edges: " << newEdgeCounter << endl;

	std::cout << "# Iterations: " << itr << endl;
	std::cout << "# Total Calculations: " << calcCnt << endl;

	long peakMemoryUsage = getPeakMemoryUsage();
	double memInGB = peakMemoryUsage / (1024.0 * 1024.0);
	std::cout << "Peak Memory Usage: " << peakMemoryUsage << " KB" << std::endl;
	std::cout << "#Peak Memory Usage (in GB): " << memInGB << " GB" << std::endl;

	// extended::
	cout << "Print total edges (old+new) for each grammar labels:: " << endl;
	// countEdge(hashset, num_nodes, grammar.labelSize, grammar.hashSymRev);

	for (uint i = 0; i < grammar.labelSize; i++)
	{
		delete[] edgeVecs[i];
		delete[] inEdgeVecs[i];
	}

	delete[] edgeVecs;
	delete[] inEdgeVecs;
}
