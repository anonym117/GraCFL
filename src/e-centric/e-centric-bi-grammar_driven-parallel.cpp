#include <omp.h>
#include "globals.hpp"
#include "globals-par.hpp"
#include "grammar.hpp"

#define TOTAL_THREADS 16

/*
 * E-centric-BI Parallel
 * Grammar Driven
 * Label and vertex indexed adjacency lists
 * Uses concurrent datastructures from Intel tbb library
 */

int main(int argc, char **argv)
{

	// Get the graph file path and grammar file path from command line argument: output file path is optional
	if (argc == 1)
	{
		std::cout << "Please provide the graph file path and grammar file path. For example, ./exc.out <graph_file> <grammar_file>" << endl;
		return 0;
	}
	else if (argc == 2)
	{
		std::cout << "Please provide the grammar file path. For example, ./exc.out <graph_file> <grammar_file>" << endl;
		return 0;
	}

	std::cout << "-----------START----------" << std::endl;
	std::cout << "--------------------------" << std::endl;

	const std::string inputGraph = argv[1];
	std::cout << "GraphFile:\t" << inputGraph << std::endl;

	std::string grammarFilePath = argv[2];
	std::cout << "GrammarFile:\t" << grammarFilePath << std::endl;
	std::cout << "--------------------------" << std::endl;

	Grammar grammar(grammarFilePath); // Read grammar

	uint num_nodes = 0;
	uint num_edges = 0;

	std::ifstream infile(inputGraph);

	std::vector<EdgeForReading> edges;
	std::unordered_set<uint> nodes;
	EdgeForReading newEdge;
	uint from, to;
	std::string label;
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

	tbb::concurrent_vector<EdgeForReading> activeQueue;
	tbb::concurrent_vector<EdgeForReading> futureQueue;

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
	// std::cout << "#Intital Unique Edge Count:\t" << initialEdgeCount << std::endl

	std::cout << "Initialization Done!" << std::endl;

	bool finished; // fixed-point iteration flag
	int itr = 0;   // Iteration counter for fixed-point iteration

	std::cout << "Start Calculations...\n";
	std::chrono::time_point<std::chrono::system_clock> start, finish;
	start = std::chrono::system_clock::now();

	// handle epsilon rules: add an edge to itself
	for (uint l = 0; l < grammar.grammar1.size(); l++)
	{
#pragma omp parallel for schedule(static) num_threads(TOTAL_THREADS)
		for (uint i = 0; i < num_nodes; i++)
		{
			// check if the new edge based on an epsilon grammar rule exists or not. l: grammar ID, 0: LHS
			if (hashset[i][grammar.grammar1[l][0]].find(i) == hashset[i][grammar.grammar1[l][0]].end())
			{
				hashset[i][grammar.grammar1[l][0]].insert(i);

				edgeVecs[grammar.grammar1[l][0]][i].push_back(i);
				inEdgeVecs[grammar.grammar1[l][0]][i].push_back(i);

				activeQueue.push_back(EdgeForReading(i, i, grammar.grammar1[l][0]));
			}
		}
	}

	do
	{
		itr++;
		std::cout << "Iteration " << itr << std::endl;

#pragma omp parallel for schedule(static) num_threads(TOTAL_THREADS)
		for (uint z = 0; z < activeQueue.size(); z++)
		{
			EdgeForReading currEdge = activeQueue[z];

			// for each grammar rule like, A = B
			for (uint g = 0; g < grammar.grammar2index[currEdge.label].size(); g++)
			{
				uint leftLabel = grammar.grammar2index[currEdge.label][g];
				// if the edge to the neighbor is labeled with B
				if (hashset[currEdge.from][leftLabel].find(currEdge.to) == hashset[currEdge.from][leftLabel].end())
				{
					hashset[currEdge.from][leftLabel].insert(currEdge.to);
					futureQueue.push_back(EdgeForReading(currEdge.from, currEdge.to, grammar.grammar2index[currEdge.label][g]));
				}
			}

			// for each grammar rule like, A = BC
			for (uint g = 0; g < grammar.grammar3indexLeft[currEdge.label].size(); g++)
			{
				uint B = currEdge.label;
				uint A = grammar.grammar3indexLeft[currEdge.label][g].second;
				uint C = grammar.grammar3indexLeft[currEdge.label][g].first;

				for (uint j = 0; j < edgeVecs[C][currEdge.to].size(); j++)
				{
					uint nbr;
					nbr = edgeVecs[C][currEdge.to][j];

					if (hashset[currEdge.from][A].find(nbr) == hashset[currEdge.from][A].end())
					{
						hashset[currEdge.from][A].insert(nbr);
						futureQueue.push_back(EdgeForReading(currEdge.from, nbr, A));
					}
				}
			}

			// for each grammar rule like, A = CB
			for (uint g = 0; g < grammar.grammar3indexRight[currEdge.label].size(); g++)
			{
				uint B = currEdge.label;
				uint A = grammar.grammar3indexRight[currEdge.label][g].second;
				uint C = grammar.grammar3indexRight[currEdge.label][g].first;

				for (uint h = 0; h < inEdgeVecs[C][currEdge.from].size(); h++)
				{
					uint inNbr;
					inNbr = inEdgeVecs[C][currEdge.from][h];
					
					if (hashset[inNbr][A].find(currEdge.to) == hashset[inNbr][A].end())
					{
						hashset[inNbr][A].insert(currEdge.to);
						futureQueue.push_back(EdgeForReading(inNbr, currEdge.to, A));
					}
				}
			}
		}

		// std::cout << "# Generated edge in this itreation:\t" << futureQueue.size() << std::endl;

		// Check if fixed-point is reached
		finished = true;
		if (futureQueue.size() > 0)
		{
			finished = false;
		}

		// Add the FUTURE edges to the adjacency lists
#pragma omp parallel for schedule(static) num_threads(TOTAL_THREADS)
		for (uint z = 0; z < futureQueue.size(); z++)
		{
			EdgeForReading edge = futureQueue[z];
			edgeVecs[edge.label][edge.from].push_back(edge.to);
			inEdgeVecs[edge.label][edge.to].push_back(edge.from);
		}

		activeQueue.swap(futureQueue);
		tbb::concurrent_vector<EdgeForReading>().swap(futureQueue);
	} while (!finished);

	// Get the calculation time
	finish = std::chrono::system_clock::now();
	std::chrono::duration<double> elapsed_seconds = finish - start;
	std::time_t finish_time = std::chrono::system_clock::to_time_t(finish);

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
