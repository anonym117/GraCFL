#include "globals.hpp"
#include "grammar.hpp"

/***
 * E-centric (Bi, rule-lookup table, new edge worklist)
 * does not have temporal pointers for NEW/OLD/FUTURE edges
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

	const std::string executableName = argv[0];
	std::cout << "ExecutableName:\t" << executableName << endl;

	const std::string inputGraph = argv[1];
	std::cout << "GraphFile:\t" << inputGraph << std::endl;

	std::string grammarFilePath = argv[2];
	std::cout << "GrammarFile:\t" << grammarFilePath << endl;
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
		// SF:: If we insert the nodes into the set, the duplicate nodes will automatically get
		// discarded by the set.
		nodes.insert(newEdge.from);
		nodes.insert(newEdge.to);
	}

	infile.close();

	std::cout << "# Vertex Count:\t" << num_nodes << std::endl;
	std::cout << "# Initial Edge Count:\t" << num_edges << std::endl;
	std::cout << "Start initializing the lists, hashset and worklists ..." << std::endl;

	// incoming edge adjacency list
	std::vector<std::vector<OutEdge>> inEdgeVecs(num_nodes);
	// level-1: vertex ID, level-2: NEW, OLD, FUTURE, level-3: outgoing edges
	std::vector<std::vector<OutEdge>> edgeVecs(num_nodes);
	// hashset for outgoing edges check if an edge exist or not
	std::unordered_set<ull> *hashset = new unordered_set<ull>[num_nodes];
	// hashset for incoming edges
	// unordered_set<ull> *inHashset = new unordered_set<ull>[num_nodes];
	// current iteration edge worklist
	std::queue<EdgeForReading> activeQueue;
	// next iteration edge worklist
	std::queue<EdgeForReading> futureQueue;

	for (uint i = 0; i < num_edges; i++)
	{
		edgeVecs[edges[i].from].push_back(OutEdge(edges[i].to, edges[i].label));
		inEdgeVecs[edges[i].to].push_back(OutEdge(edges[i].from, edges[i].label));

		hashset[edges[i].from].insert(COMBINE(edges[i].to, edges[i].label));

		// populate the current iteration worklist with initial edges
		activeQueue.push(edges[i]);
	}

	edges.clear();

	// count the total no of initial unique edges
	ull initialEdgeCount = countEdge(hashset, num_nodes);

	std::cout << "Initialization Done!" << std::endl;

	bool finished;
	int itr = 0;

	std::cout << "Start Calculations...\n";
	std::chrono::time_point<std::chrono::system_clock> start, finish;
	start = std::chrono::system_clock::now();

	// handle epsilon rules: add an edge to itself
	for (uint l = 0; l < grammar.grammar1.size(); l++)
	{
		for (uint i = 0; i < num_nodes; i++)
		{
			// check if the new edge based on an epsilon grammar rule exists or not.
			if (hashset[i].find(COMBINE(i, grammar.grammar1[l][0])) == hashset[i].end())
			{
				hashset[i].insert(COMBINE(i, grammar.grammar1[l][0]));
				
				edgeVecs[i].push_back(OutEdge(i, grammar.grammar1[l][0]));
				inEdgeVecs[i].push_back(OutEdge(i, grammar.grammar1[l][0]));

				activeQueue.push(EdgeForReading(i, i, grammar.grammar1[l][0]));
			}
		}
	}

	do
	{
		itr++;
		std::cout << "Iteration " << itr << std::endl;

		while (!activeQueue.empty())
		{
			EdgeForReading currEdge = activeQueue.front();
			activeQueue.pop();

			// for each grammar rule like, A = B
			vector<uint> leftLabels = grammar.rule2(currEdge.label);
			uint leftLabelSize = leftLabels.size();

			for (uint g = 0; g < leftLabelSize; g++)
			{
				// if the edge to the neighbor is labeled with B
				if (hashset[currEdge.from].find(COMBINE(currEdge.to, leftLabels[g])) == hashset[currEdge.from].end())
				{
					hashset[currEdge.from].insert(COMBINE(currEdge.to, leftLabels[g]));
					futureQueue.push(EdgeForReading(currEdge.from, currEdge.to, leftLabels[g]));
				}
			}

			
			// for each grammar rule like, A = BC, B = grammar.grammar3[g][1] = currEdge.label
			for (uint j = 0; j < edgeVecs[currEdge.to].size(); j++)
			{
				OutEdge nbr;

				// NEW outgoing edges of currEdge.to
				nbr = edgeVecs[currEdge.to][j];
				// get the left labels associated with rightlabel1, and rightlabel2
				leftLabels = grammar.rule3(currEdge.label, nbr.label);
				leftLabelSize = leftLabels.size();

				for (uint g = 0; g < leftLabelSize; g++)
				{
					if (hashset[currEdge.from].find(COMBINE(nbr.end, leftLabels[g])) == hashset[currEdge.from].end())
					{
						hashset[currEdge.from].insert(COMBINE(nbr.end, leftLabels[g]));
						futureQueue.push(EdgeForReading(currEdge.from, nbr.end, leftLabels[g]));
					}
				}
			}

			// for each grammar rule like, A = CB, B = grammar.grammar3[g][2] = currEdge.label
			for (uint h = 0; h < inEdgeVecs[currEdge.from].size(); h++)
			{
				OutEdge inNbr;
				inNbr = inEdgeVecs[currEdge.from][h];
				// get the left labels associated with rightlabel1, and rightlabel2
				leftLabels = grammar.rule3(inNbr.label, currEdge.label);
				leftLabelSize = leftLabels.size();

				for (uint g = 0; g < leftLabelSize; g++)
				{
					if (hashset[inNbr.end].find(COMBINE(currEdge.to, leftLabels[g])) == hashset[inNbr.end].end())
					{
						hashset[inNbr.end].insert(COMBINE(currEdge.to, leftLabels[g]));
						futureQueue.push(EdgeForReading(inNbr.end, currEdge.to, leftLabels[g]));
					}
				}
			}
		}

		// Add the future edges in the adjacency list and swap the two worklist edges
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

	std::cout << "Calculation Done!" << std::endl;
	ull totalNewEdgeCount = countEdge(hashset, num_nodes) - initialEdgeCount;

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
