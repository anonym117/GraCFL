#include "globals.hpp"
#include "grammar.hpp"

/***
 * V-centric (bi, sliding ptrs, grammar driven)
 * adjacency list: first level: grammar label, second level: vertex
 * Bi directional traversing of the Graph (incoming and outgoing edges)
 * Three BufferEdge lists for OLD, NEW, and FUTURE edges
 * Check the Future edge flag whenever an edge is created
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

    std::cout << "# Vertex Count:\t" << num_nodes << std::endl;
	std::cout << "# Initial Edge Count:\t" << num_edges << std::endl;
	std::cout << "Start initializing the lists, hashset and worklists ..." << std::endl;

    // level-1: labelID, level 2: nodeID, level-3: NEW, OLD, FUTURE vecs, level-4: and incoming edges
    vector<vector<vector<vector<uint>>>> inEdgeVecs(grammar.labelSize, vector<vector<vector<uint>>>(num_nodes, vector<vector<uint>>(3, vector<uint>())));
    // level-1: labelID, level 2: nodeID, level-3: NEW, OLD, FUTURE vecs, level-4: and outgoing edges
    vector<vector<vector<vector<uint>>>> edgeVecs(grammar.labelSize, vector<vector<vector<uint>>>(num_nodes, vector<vector<uint>>(3, vector<uint>())));
    vector<vector<unordered_set<ull>>> hashset(num_nodes, vector<unordered_set<ull>>(grammar.labelSize, unordered_set<ull>()));

    for (uint i = 0; i < num_edges; i++)
    {
        edgeVecs[edges[i].label][edges[i].from][NEW].push_back(edges[i].to);
        inEdgeVecs[edges[i].label][edges[i].to][NEW].push_back(edges[i].from);
        // insert edge into the hashset
        hashset[edges[i].from][edges[i].label].insert(edges[i].to);
    }

    edges.clear();

    // count the total no of initial unique edge
    uint initialEdgeCount = countEdge(hashset, num_nodes, grammar.labelSize);

    std::cout << "Initialization Done!" << std::endl;

    bool finished; // fixed-point iteration flag
	int itr = 0; // Iteration counter for fixed-point iteration

	std::cout << "Start Calculations...\n";
    std::chrono::time_point<std::chrono::steady_clock> start, finish;
    start = std::chrono::steady_clock::now();

    // handle epsilon rules: add an edge to itself
    for (uint l = 0; l < grammar.grammar1.size(); l++)
    {
        for (uint i = 0; i < num_nodes; i++)
        {
            // check if the new edge based on an epsilon grammar rule exists or not. l: grammar ID, 0: LHS
            if (hashset[i][grammar.grammar1[l][0]].find(i) == hashset[i][grammar.grammar1[l][0]].end())
            {
                // insert edge into hashset
                hashset[i][grammar.grammar1[l][0]].insert(i);
                // insert edge into the graph
                edgeVecs[grammar.grammar1[l][0]][i][NEW].push_back(i);
                inEdgeVecs[grammar.grammar1[l][0]][i][NEW].push_back(i);
            }
        }
    }

    do
    {
        finished = true;
        itr++;

        std::cout << "Iteration number " << itr << std::endl;

        // for each grammar rule like A --> B
        for (uint g = 0; g < grammar.labelSize; g++)
        {
            for (uint i = 0; i < num_nodes; i++)
            {
                uint nbr;
                //  the valid index range is [START_NEW, END_NEW-1]
                uint START_NEW = 0;
                uint END_NEW = inEdgeVecs[g][i][NEW].size();

                // for each new edge
                for (uint j = START_NEW; j < END_NEW; j++)
                {
                    uint inNbr = inEdgeVecs[g][i][NEW][j];
                    // rule: A = B
                    uint g2Size = grammar.grammar2index[g].size();
                    for (uint m = 0; m < g2Size; m++)
                    {
                        // calcCnt++;
                        uint A = grammar.grammar2index[g][m];

                        if (hashset[inNbr][A].find(i) == hashset[inNbr][A].end())
                        {
                            finished = false;
                            hashset[inNbr][A].insert(i);

                            edgeVecs[A][inNbr][FUTURE].push_back(i);
                            inEdgeVecs[A][i][FUTURE].push_back(inNbr);

                            // No need to update the pointers. Because the FUTURE_START starts from the
                            // NEW_END, and NEW_END is already updated
                        }
                    }

                    // Rule: A = BC
                    uint g3SizeLeft = grammar.grammar3indexLeft[g].size();
                    for (uint m = 0; m < g3SizeLeft; m++)
                    {
                        uint C = grammar.grammar3indexLeft[g][m].first;
                        uint A = grammar.grammar3indexLeft[g][m].second;

                        uint START_NEW_OUT = 0;
                        uint END_NEW_OUT = edgeVecs[C][i][NEW].size();

                        // new * new
                        for (uint h = START_NEW_OUT; h < END_NEW_OUT; h++)
                        {
                            uint nbr = edgeVecs[C][i][NEW][h];
                            if (hashset[inNbr][A].find(nbr) == hashset[inNbr][A].end())
                            {
                                finished = false;
                                hashset[inNbr][A].insert(nbr);
                                edgeVecs[A][inNbr][FUTURE].push_back(nbr);
                                inEdgeVecs[A][nbr][FUTURE].push_back(inNbr);
                            }
                        }

                        // new * old
                        uint START_OLD_OUT = 0;
                        uint END_OLD_OUT = edgeVecs[C][i][OLD].size();

                        for (uint h = START_OLD_OUT; h < END_OLD_OUT; h++)
                        {
                            uint nbr = edgeVecs[C][i][OLD][h];
                            if (hashset[inNbr][A].find(nbr) == hashset[inNbr][A].end())
                            {
                                finished = false;
                                hashset[inNbr][A].insert(nbr);
                                edgeVecs[A][inNbr][FUTURE].push_back(nbr);
                                inEdgeVecs[A][nbr][FUTURE].push_back(inNbr);
                            }
                        }

                    }
                }

                // rule: A = CB
                // C (old) * B (new)
                // all the OLD, and NEW outgoing edges of the first edge
                uint START_NEW_OUT = 0;
                uint END_NEW_OUT = edgeVecs[g][i][NEW].size();

                for (uint j = START_NEW_OUT; j < END_NEW_OUT; j++)
                {   
                    uint nbr = edgeVecs[g][i][NEW][j];
                    uint g3SizeRight = grammar.grammar3indexRight[g].size();

                    for (uint m = 0; m < g3SizeRight; m++)
                    {
                        uint C = grammar.grammar3indexRight[g][m].first;
                        uint A = grammar.grammar3indexRight[g][m].second;

                        uint START_OLD_IN = 0;
                        uint END_OLD_IN = inEdgeVecs[C][i][OLD].size();

                        for (uint h = START_OLD_IN; h < END_OLD_IN; h++)
                        {
                            uint inNbr = inEdgeVecs[C][i][OLD][h];

                            // calcCnt++;
                            if (hashset[inNbr][A].find(nbr) == hashset[inNbr][A].end())
                            {
                                finished = false;
                                hashset[inNbr][A].insert(nbr);

                                edgeVecs[A][inNbr][FUTURE].push_back(nbr);
                                inEdgeVecs[A][nbr][FUTURE].push_back(inNbr);
                            }
                        }
                    }
                    
                }
            }
        }

        for (uint g = 0; g < grammar.labelSize; g++)
        {
            // update the sliding/temporal pointers
            for (int i = 0; i < num_nodes; i++)
            {
                edgeVecs[g][i][OLD].insert(edgeVecs[g][i][OLD].begin(), edgeVecs[g][i][NEW].begin(), edgeVecs[g][i][NEW].end());
                inEdgeVecs[g][i][OLD].insert(inEdgeVecs[g][i][OLD].begin(), inEdgeVecs[g][i][NEW].begin(), inEdgeVecs[g][i][NEW].end());

                edgeVecs[g][i][NEW] = edgeVecs[g][i][FUTURE];
                inEdgeVecs[g][i][NEW] = inEdgeVecs[g][i][FUTURE];

                edgeVecs[g][i][FUTURE].clear();
                inEdgeVecs[g][i][FUTURE].clear();
            }
        }
    } while (!finished);

    finish = std::chrono::steady_clock::now();
    std::chrono::duration<double> elapsed_seconds = finish - start;
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
}
