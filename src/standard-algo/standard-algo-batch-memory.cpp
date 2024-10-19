#include "globals.hpp"
#include "grammar.hpp"

/*
 * E-centric w/ rule-idx and label-idx-adj_list
 * new edge worklist, iterative
 */

long getPeakMemoryUsage()
{
	std::string line;
	std::ifstream statusFile("/proc/self/status");

	while (getline(statusFile, line))
	{
		if (line.substr(0, 7) == "VmPeak:")
		{
			long memoryKb;
			std::istringstream iss(line.substr(7));
			iss >> memoryKb;
			return memoryKb; // Peak memory usage in kilobytes
		}
	}

	return -1; // Return -1 if VmPeak is not found
}

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

	vector<EdgeForReading2> edges;
	unordered_set<uint> nodes;
	EdgeForReading2 newEdge;
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

	vector<uint> **inEdgeVecs = new vector<uint> *[grammar.labelSize];
	// level-1: vertex ID, level-2: NEW, OLD, FUTURE, level-3: outgoing edges
	vector<uint> **edgeVecs = new vector<uint> *[grammar.labelSize];

	// check if an edge exist or not
	// unordered_set<ull> **hashset = new unordered_set<ull> *[num_nodes];
	vector<vector<unordered_set<ull>>> hashset(num_nodes, vector<unordered_set<ull>>(grammar.labelSize, unordered_set<ull>()));

	queue<EdgeForReading2> activeQueue;
	queue<EdgeForReading2> futureQueue;

	cout << "#nodes " << num_nodes << endl;
	cout << "SF::#nodes " << nodes.size() << endl;
	cout << "#edges " << num_edges << endl;

	cout << "Start making sets and the hash ...\n";

	for (uint i = 0; i < grammar.labelSize; i++)
	{
		edgeVecs[i] = new vector<uint>[num_nodes];
		inEdgeVecs[i] = new vector<uint>[num_nodes];
	}

	for (uint i = 0; i < num_edges; i++)
	{
		edgeVecs[edges[i].label][edges[i].from].push_back(edges[i].to);
		inEdgeVecs[edges[i].label][edges[i].to].push_back(edges[i].from);

		hashset[edges[i].from][edges[i].label].insert(edges[i].to);

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
			// check if the new edge based on an epsilon grammar rule exists or not. l: grammar ID, 0: LHS
			if (hashset[i][grammar.grammar1[l][0]].find(i) == hashset[i][grammar.grammar1[l][0]].end())
			{
				// insert into hashset
				hashset[i][grammar.grammar1[l][0]].insert(i);
				// insert into the graph
				edgeVecs[grammar.grammar1[l][0]][i].push_back(i);
				inEdgeVecs[grammar.grammar1[l][0]][i].push_back(i);

				activeQueue.push(EdgeForReading2(i, i, grammar.grammar1[l][0]));

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
			EdgeForReading2 currEdge = activeQueue.front();
			activeQueue.pop();

			// for each grammar rule like A --> B
			for (uint g = 0; g < grammar.grammar2index[currEdge.label].size(); g++)
			{
				// calcCnt++;
				// calcCnt2++;
				uint leftLabel = grammar.grammar2index[currEdge.label][g];
				// if the edge to the neighbor is labeled with B
				if (hashset[currEdge.from][leftLabel].find(currEdge.to) == hashset[currEdge.from][leftLabel].end())
				{
					hashset[currEdge.from][leftLabel].insert(currEdge.to);

					// the new edge is added to FUTURE list, due to potential conflicts in parallel version
					// if topo-driven is run in serial, only OLD and NEW are needed.
					// Graspan did not mention FUTURE list, but it puts the new edges into a seperate list, which
					// is similar to FUTURE list's purpose.
					// for data-driven, we only need one list (no NEW OLD FUTURE)

					futureQueue.push(EdgeForReading2(currEdge.from, currEdge.to, grammar.grammar2index[currEdge.label][g]));

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
					uint nbr;

					// A = BC, B = grammar.grammar3[g][1] = currEdge.label
					// NEW outgoing edges of currEdge.to
					nbr = edgeVecs[C][currEdge.to][j];
					// A = BC, C = grammar.grammar3[g][2]
					// calcCnt++;
					// calcCnt3_left++;

					if (hashset[currEdge.from][A].find(nbr) == hashset[currEdge.from][A].end())
					{
						hashset[currEdge.from][A].insert(nbr);
						futureQueue.push(EdgeForReading2(currEdge.from, nbr, A));
						// newEdgeCounter++;
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
					uint inNbr;
					inNbr = inEdgeVecs[C][currEdge.from][h];
					// A = CB, C = grammar.grammar3[g][1]
					// calcCnt++;
					// calcCnt3_right++;

					if (hashset[inNbr][A].find(currEdge.to) == hashset[inNbr][A].end())
					{
						hashset[inNbr][A].insert(currEdge.to);
						futureQueue.push(EdgeForReading2(inNbr, currEdge.to, A));
						// newEdgeCounter++;
					}
				}
			}
		}

		finished = true;

		queue<EdgeForReading2> tempQueue = futureQueue;

		while (!tempQueue.empty())
		{
			EdgeForReading2 edge = tempQueue.front();
			tempQueue.pop();

			edgeVecs[edge.label][edge.from].push_back(edge.to);
			inEdgeVecs[edge.label][edge.to].push_back(edge.from);
		}

		if (!futureQueue.empty())
			finished = false;
		swap(activeQueue, futureQueue);
	} while (!finished);

	finish = std::chrono::system_clock::now();
	std::chrono::duration<double> elapsed_seconds = finish - start;
	std::time_t finish_time = std::chrono::system_clock::to_time_t(finish);

	uint totalNewEdgeCount = countEdge(hashset, num_nodes, grammar.labelSize) - initialEdgeCount;

	cout << "**************************" << endl;
	std::cout << "# Total time = " << elapsed_seconds.count() << std::endl;
	std::cout << "# Total iterations = " << itr << endl;

	cout << "SF:: # Number of new edges: " << totalNewEdgeCount << endl;
	cout << "AM:: # Number of new edges: " << newEdgeCounter << endl;

	cout << "# Iterations: " << itr << endl;
	cout << "# Total Calculations: " << calcCnt << endl;
	long peakMemoryUsage = getPeakMemoryUsage();
	double memInGB = peakMemoryUsage / (1024.0 * 1024.0);
	std::cout << "Peak Memory Usage: " << peakMemoryUsage << " KB" << std::endl;
	std::cout << "#Peak Memory Usage (in GB): " << memInGB << " GB" << std::endl;
}