//#include <cilk/cilk.h>
//#include <cilk/cilk_api.h>

#include "globals.hpp"
#include "grammar.hpp"

/***
 * e-centric (bi, rule-lookup table, new edge worklist)
 * does not have temporal pointers for NEW/OLD/FUTURE edges
 */

int main(int argc, char **argv)
{

	// bool debug = true;
	//  Get the graph file path and grammar file path from command line argument
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
		// SF:: If we insert the nodes into the set, the duplicate nodes will automatically get
		// discarded by the set.
		nodes.insert(newEdge.from);
		nodes.insert(newEdge.to);
	}

	infile.close();

	// incoming edge adjacency list
	vector<vector<OutEdge>> inEdgeVecs(num_nodes);
	// level-1: vertex ID, level-2: NEW, OLD, FUTURE, level-3: outgoing edges
	vector<vector<OutEdge>> edgeVecs(num_nodes);
	// hashset for outgoing edges check if an edge exist or not
	unordered_set<ull> *hashset = new unordered_set<ull>[num_nodes];
	// hashset for incoming edges
	// unordered_set<ull> *inHashset = new unordered_set<ull>[num_nodes];
	// current iteration edge worklist
	queue<EdgeForReading> activeQueue;
	// next iteration edge worklist
	queue<EdgeForReading> futureQueue;

	cout << "#nodes " << num_nodes << endl;
	cout << "#edges " << num_edges << endl;

	cout << "Start making sets and the hash ...\n";
	for (uint i = 0; i < num_edges; i++)
	{
		edgeVecs[edges[i].from].push_back(OutEdge(edges[i].to, edges[i].label));
		inEdgeVecs[edges[i].to].push_back(OutEdge(edges[i].from, edges[i].label));

		hashset[edges[i].from].insert(COMBINE(edges[i].to, edges[i].label));
		// inHashset[edges[i].to].insert(COMBINE(edges[i].from, edges[i].label));

		// populate the current iteration worklist with initial edges
		activeQueue.push(edges[i]);
	}

	edges.clear();

	// count the total no of initial unique edges
	ull initialEdgeCount = countEdge(hashset, num_nodes);

	cout << "Done!\n";

	// atomic<int> newEdgeCounter(0);
	bool finished;
	int itr = 0;
	ull calcCnt = 0;
	double elapsed_seconds_comp = 0.0;
	double elapsed_seconds_transfer = 0.0;

	// start the timer
	// currently excludes the preprocessing time
	std::chrono::time_point<std::chrono::system_clock> start, finish;
	start = std::chrono::system_clock::now();

	// handle epsilon rules: add an edge to itself
	// grammar1 is for epsilon rules A --> e
	// grammar2 is for one symbol on RHS A --> B
	// grammar3 is for two symbols on RHS A --> BC
	for (uint l = 0; l < grammar.grammar1.size(); l++)
	{
		for (uint i = 0; i < num_nodes; i++)
		{
			// calcCnt++;

			// check if the new edge based on an epsilon grammar rule exists or not.
			if (hashset[i].find(COMBINE(i, grammar.grammar1[l][0])) == hashset[i].end())
			{
				// insert into hashset
				hashset[i].insert(COMBINE(i, grammar.grammar1[l][0]));
				// inHashset[i].insert(COMBINE(i, grammar.grammar1[l][0]));
				// insert into the graph
				edgeVecs[i].push_back(OutEdge(i, grammar.grammar1[l][0]));
				inEdgeVecs[i].push_back(OutEdge(i, grammar.grammar1[l][0]));

				activeQueue.push(EdgeForReading(i, i, grammar.grammar1[l][0]));

				// newEdgeCounter++;
			}
		}
	}

	cout << "********************\n";

	do
	{
		itr++;

		while (!activeQueue.empty())
		{
			EdgeForReading currEdge = activeQueue.front();
			activeQueue.pop();

			// for each grammar rule, A --> B
			// take all the left labels, A for currEdge label B
			vector<uint> leftLabels = grammar.rule2(currEdge.label);
			uint leftLabelSize = leftLabels.size();

			for (uint g = 0; g < leftLabelSize; g++)
			{
				// calcCnt++;

				// if the edge to the neighbor is labeled with B
				if (hashset[currEdge.from].find(COMBINE(currEdge.to, leftLabels[g])) == hashset[currEdge.from].end())
				{
					hashset[currEdge.from].insert(COMBINE(currEdge.to, leftLabels[g]));
					// inHashset[currEdge.to].insert(COMBINE(currEdge.from, leftLabels[g]));

					futureQueue.push(EdgeForReading(currEdge.from, currEdge.to, leftLabels[g]));
					// newEdgeCounter++;
				}
			}

			// A = BC, B = grammar.grammar3[g][1] = currEdge.label
			for (uint j = 0; j < edgeVecs[currEdge.to].size(); j++)
			{
				OutEdge nbr;

				// NEW outgoing edges of currEdge.to
				nbr = edgeVecs[currEdge.to][j];
				// get the left labels associated with rightlabel1, and rightlabel2
				leftLabels = grammar.rule3(currEdge.label, nbr.label);
				leftLabelSize = leftLabels.size();

				// A = BC, C = grammar.grammar3[g][2]
				for (uint g = 0; g < leftLabelSize; g++)
				{
					// calcCnt++;

					if (hashset[currEdge.from].find(COMBINE(nbr.end, leftLabels[g])) == hashset[currEdge.from].end())
					{
						hashset[currEdge.from].insert(COMBINE(nbr.end, leftLabels[g]));
						// inHashset[nbr.end].insert(COMBINE(currEdge.from, leftLabels[g]));

						futureQueue.push(EdgeForReading(currEdge.from, nbr.end, leftLabels[g]));
						// newEdgeCounter++;
					}
				}
			}

			// A = CB, B = grammar.grammar3[g][2] = currEdge.label
			for (uint h = 0; h < inEdgeVecs[currEdge.from].size(); h++)
			{
				OutEdge inNbr;
				inNbr = inEdgeVecs[currEdge.from][h];
				// get the left labels associated with rightlabel1, and rightlabel2
				leftLabels = grammar.rule3(inNbr.label, currEdge.label);
				leftLabelSize = leftLabels.size();

				for (uint g = 0; g < leftLabelSize; g++)
				{
					// calcCnt++;

					if (hashset[inNbr.end].find(COMBINE(currEdge.to, leftLabels[g])) == hashset[inNbr.end].end())
					{
						hashset[inNbr.end].insert(COMBINE(currEdge.to, leftLabels[g]));
						// inHashset[currEdge.to].insert(COMBINE(inNbr.end, leftLabels[g]));

						futureQueue.push(EdgeForReading(inNbr.end, currEdge.to, leftLabels[g]));
						// newEdgeCounter++;
					}
				}
			}
		}

		cout << "Iteration number " << itr << endl;

		// add the future edges in the adjacency list and swap the two worklist edges
		queue<EdgeForReading> tempQueue = futureQueue;
		while (!tempQueue.empty())
		{
			EdgeForReading currEdge = tempQueue.front();
			tempQueue.pop();
			edgeVecs[currEdge.from].push_back(OutEdge(currEdge.to, currEdge.label));
			inEdgeVecs[currEdge.to].push_back(OutEdge(currEdge.from, currEdge.label));
		}

		finished = true;

		if (!futureQueue.empty())
			finished = false;
		swap(activeQueue, futureQueue);

	} while (!finished);

	finish = std::chrono::system_clock::now();
	std::chrono::duration<double> elapsed_seconds = finish - start;
	std::time_t finish_time = std::chrono::system_clock::to_time_t(finish);

	ull totalNewEdgeCount = countEdge(hashset, num_nodes) - initialEdgeCount;

	cout << "**************************" << endl;
	std::cout << "# Total time = " << elapsed_seconds.count() << std::endl;
	// std::cout << "# Total computation time = " << elapsed_seconds_comp << std::endl;
	// std::cout << "# Total data transfer time = " << elapsed_seconds_transfer << std::endl;

	cout << "SF:: # Number of new edges: " << totalNewEdgeCount << endl;
	// cout << "AM:: # Number of new edges: " << newEdgeCounter << endl;

	cout << "# Iterations: " << itr << endl;
	cout << "# Total Calculations: " << calcCnt << endl;
}