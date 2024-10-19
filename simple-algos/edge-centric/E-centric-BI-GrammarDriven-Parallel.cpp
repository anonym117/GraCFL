#include <omp.h>
#include "globals.hpp"
#include "grammar.hpp"

#define TOTAL_THREADS 16

/*
 * E-centric w/ rule-idx and label-idx-adj_list
 * new edge worklist, iterative, standard
 * Parallel, concurrent vector, unordered_set used from
 * Intel tbb library
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

	tbb::concurrent_vector<uint> **edgeVecs = new tbb::concurrent_vector<uint> *[grammar.labelSize];
	tbb::concurrent_vector<uint> **inEdgeVecs = new tbb::concurrent_vector<uint> *[grammar.labelSize];

	for (uint i = 0; i < grammar.labelSize; i++)
	{
		edgeVecs[i] = new tbb::concurrent_vector<uint>[num_nodes];
		inEdgeVecs[i] = new tbb::concurrent_vector<uint>[num_nodes];
	}

	// check if an edge exist or not
	tbb::concurrent_unordered_set<ull> **hashset = new tbb::concurrent_unordered_set<ull> *[num_nodes];
	for (uint i = 0; i < num_nodes; i++)
	{
		hashset[i] = new tbb::concurrent_unordered_set<ull>[grammar.labelSize];
	}

	tbb::concurrent_vector<EdgeForReading2> activeQueue;
	tbb::concurrent_vector<EdgeForReading2> futureQueue;

	cout << "#nodes " << num_nodes << endl;
	cout << "SF::#nodes " << nodes.size() << endl;
	cout << "#edges " << num_edges << endl;

	cout << "Start making sets and the hash ...\n";

	for (uint i = 0; i < num_edges; i++)
	{
		edgeVecs[edges[i].label][edges[i].from].push_back(edges[i].to);
		inEdgeVecs[edges[i].label][edges[i].to].push_back(edges[i].from);

		hashset[edges[i].from][edges[i].label].insert(edges[i].to);

		activeQueue.push_back(edges[i]);
	}

	edges.clear();

	// count the total no of initial unique edge
	uint initialEdgeCount = countEdge(hashset, num_nodes, grammar.labelSize);

	cout << "Done!\n";

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

	// handle epsilon rules: add an edge to itself
	// grammar1 is for epsilon rules A --> e
	// grammar2 is for one symbol on RHS A --> B
	// grammar3 is for two symbols on RHS A --> BC

	for (uint l = 0; l < grammar.grammar1.size(); l++)
	{
#pragma omp parallel for schedule(static) num_threads(TOTAL_THREADS)
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

				activeQueue.push_back(EdgeForReading2(i, i, grammar.grammar1[l][0]));

				newEdgeCounter++;
			}
		}
	}

	cout << "********************\n";

	do
	{
		itr++;

#pragma omp parallel for schedule(static) num_threads(TOTAL_THREADS)
		for (uint z = 0; z < activeQueue.size(); z++)
		{
			EdgeForReading2 currEdge = activeQueue[z];

			// for each grammar rule like A --> B
			for (uint g = 0; g < grammar.grammar2index[currEdge.label].size(); g++)
			{
				uint leftLabel = grammar.grammar2index[currEdge.label][g];
				// if the edge to the neighbor is labeled with B
				if (hashset[currEdge.from][leftLabel].find(currEdge.to) == hashset[currEdge.from][leftLabel].end())
				{
					hashset[currEdge.from][leftLabel].insert(currEdge.to);
					futureQueue.push_back(EdgeForReading2(currEdge.from, currEdge.to, grammar.grammar2index[currEdge.label][g]));
				}
			}

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

					if (hashset[currEdge.from][A].find(nbr) == hashset[currEdge.from][A].end())
					{
						hashset[currEdge.from][A].insert(nbr);
						futureQueue.push_back(EdgeForReading2(currEdge.from, nbr, A));
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
					if (hashset[inNbr][A].find(currEdge.to) == hashset[inNbr][A].end())
					{
						hashset[inNbr][A].insert(currEdge.to);
						futureQueue.push_back(EdgeForReading2(inNbr, currEdge.to, A));
						// newEdgeCounter++;
					}
				}
			}
		}

		cout << "Iteration number " << itr << endl;

		finished = true;
		if (futureQueue.size() > 0)
		{
			finished = false;
		}

#pragma omp parallel for schedule(static) num_threads(TOTAL_THREADS)
		for (uint z = 0; z < futureQueue.size(); z++)
		{
			EdgeForReading2 edge = futureQueue[z];

			edgeVecs[edge.label][edge.from].push_back(edge.to);
			inEdgeVecs[edge.label][edge.to].push_back(edge.from);
		}

		activeQueue.swap(futureQueue);
		tbb::concurrent_vector<EdgeForReading2>().swap(futureQueue);
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

	for (uint i = 0; i < grammar.labelSize; i++)
	{
		delete[] edgeVecs[i];
		delete[] inEdgeVecs[i];
	}

	for (uint i = 0; i < num_nodes; i++)
	{
		delete[] hashset[i];
	}

	delete[] edgeVecs;
	delete[] inEdgeVecs;
	delete[] hashset;
}
