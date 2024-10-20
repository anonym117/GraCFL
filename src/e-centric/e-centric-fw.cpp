#include "globals.hpp"
#include "grammar.hpp"

/***
 * E-centric (forward, sliding ptrs, grammar driven)
 * Grammar driven
 * One BufferEdge list instead of three for OLD, NEW, and FUTURE edges
 * BufferEdge struct (adjcency list is made with this struct) holds the pointes for the OLD, NEW, and FUTURE edges in the single BufferEdge list
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

    // level-1: vertex ID, level-2: NEW, OLD, FUTURE pointers and outgoing edges
    std::vector<std::vector<BufferEdge>> edgeVecs(grammar.labelSize, std::vector<BufferEdge>(num_nodes));
    std::vector<std::vector<std::unordered_set<ull>>> hashset(num_nodes, std::vector<std::unordered_set<ull>>(grammar.labelSize, unordered_set<ull>()));
    std::queue<EdgeForReading> newWorklist;
    std::queue<EdgeForReading> futureNewWorklist;
    std::vector<EdgeForReading> oldWorklist;

    for (uint i = 0; i < num_edges; i++)
    {
        edgeVecs[edges[i].label][edges[i].from].vertexList.push_back(edges[i].to);
        // update the BufferEdge pointers
        edgeVecs[edges[i].label][edges[i].from].NEW_END++;
        // insert edge into the hashset
        hashset[edges[i].from][edges[i].label].insert(edges[i].to);
        newWorklist.push(edges[i]);
    }

    edges.clear();

    // count the total no of initial unique edge
    uint initialEdgeCount = countEdge(hashset, num_nodes, grammar.labelSize);

    std::cout << "Initialization Done!" << std::endl;

    std::chrono::time_point<std::chrono::steady_clock> start, finish;
    std::chrono::time_point<std::chrono::steady_clock> startTransfer, finishTransfer;
    start = std::chrono::steady_clock::now();

    bool finished; // fixed-point iteration flag
	int itr = 0; // Iteration counter for fixed-point iteration

    double elapsed_seconds_comp = 0.0;
    double elapsed_seconds_transfer = 0.0;

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
                edgeVecs[grammar.grammar1[l][0]][i].vertexList.push_back(i);
                // update the sliding/temporal pointers
                edgeVecs[grammar.grammar1[l][0]][i].NEW_END++;
                newWorklist.push(EdgeForReading(i, i, grammar.grammar1[l][0]));
            }
        }
    }

    do
    {

        finished = true;
        itr++;
        std::cout << "Iteration " << itr << std::endl;

        uint CURR_OLD_SIZE = oldWorklist.size();

        while (!newWorklist.empty())
        {
            EdgeForReading currEdge = newWorklist.front();
            newWorklist.pop();

            // for each grammar rule, A --> B
            for (uint g = 0; g < grammar.grammar2index[currEdge.label].size(); g++)
            {
                uint leftLabel = grammar.grammar2index[currEdge.label][g];
                // if the edge to the neighbor is labeled with B
                if (hashset[currEdge.from][leftLabel].find(currEdge.to) == hashset[currEdge.from][leftLabel].end())
                {
                    // finished = false;
                    hashset[currEdge.from][leftLabel].insert(currEdge.to);

                    // the new edge is added to FUTURE list, due to potential conflicts in parallel version
                    // if topo-driven is run in serial, only OLD and NEW are needed.
                    // Graspan did not mention FUTURE list, but it puts the new edges into a seperate list, which
                    // is similar to FUTURE list's purpose.
                    // for data-driven, we only need one list (no NEW OLD FUTURE)
                    edgeVecs[leftLabel][currEdge.from].vertexList.push_back(currEdge.to);
                    futureNewWorklist.push(EdgeForReading(currEdge.from, currEdge.to, leftLabel));
                }
            }

            // A = BC
            for (uint g = 0; g < grammar.grammar3indexLeft[currEdge.label].size(); g++)
            {
                uint B = currEdge.label;
                uint A = grammar.grammar3indexLeft[currEdge.label][g].second;
                uint C = grammar.grammar3indexLeft[currEdge.label][g].first;

                uint OLD_START = 0;
                uint NEW_END = edgeVecs[C][currEdge.to].NEW_END;

                for (uint j = OLD_START; j < NEW_END; j++)
                {
                    uint nbr;
                    // A = BC, B = grammar.grammar3[g][1] = currEdge.label
                    // NEW outgoing edges of currEdge.to
                    nbr = edgeVecs[C][currEdge.to].vertexList[j];
                    // A = BC, C = grammar.grammar3[g][2]

                    if (hashset[currEdge.from][A].find(nbr) == hashset[currEdge.from][A].end())
                    {
                        // finished = false;
                        hashset[currEdge.from][A].insert(nbr);
                        futureNewWorklist.push(EdgeForReading(currEdge.from, nbr, A));
                        edgeVecs[A][currEdge.from].vertexList.push_back(nbr);
                    }
                }
            }

            oldWorklist.push_back(currEdge);
        }

        for (uint i = 0; i < CURR_OLD_SIZE; i++)
        {
            EdgeForReading currEdge = oldWorklist[i];
            // A = BC
            for (uint g = 0; g < grammar.grammar3indexLeft[currEdge.label].size(); g++)
            {
                uint B = currEdge.label;
                uint A = grammar.grammar3indexLeft[currEdge.label][g].second;
                uint C = grammar.grammar3indexLeft[currEdge.label][g].first;

                uint NEW_START = edgeVecs[C][currEdge.to].OLD_END;
                uint NEW_END = edgeVecs[C][currEdge.to].NEW_END;

                for (uint j = NEW_START; j < NEW_END; j++)
                {
                    uint nbr;
                    // A = BC, B = grammar.grammar3[g][1] = currEdge.label
                    // NEW outgoing edges of currEdge.to
                    nbr = edgeVecs[C][currEdge.to].vertexList[j];

                    if (hashset[currEdge.from][A].find(nbr) == hashset[currEdge.from][A].end())
                    {
                        hashset[currEdge.from][A].insert(nbr);
                        futureNewWorklist.push(EdgeForReading(currEdge.from, nbr, A));
                        edgeVecs[A][currEdge.from].vertexList.push_back(nbr);
                    }
                }
            }
        }

        startTransfer = std::chrono::steady_clock::now();
        for (uint g = 0; g < grammar.labelSize; g++)
        {
            // update the sliding/temporal pointers
            for (int i = 0; i < num_nodes; i++)
            {
                edgeVecs[g][i].OLD_END = edgeVecs[g][i].NEW_END;
                edgeVecs[g][i].NEW_END = edgeVecs[g][i].vertexList.size();
            }
        }

        if (!futureNewWorklist.empty())
        {
            finished = false;
        }

        swap(newWorklist, futureNewWorklist);
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
