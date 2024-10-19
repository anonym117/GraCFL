#ifndef GLOBALS_HPP
#define GLOBALS_HPP

#include <iostream>
#include <stdlib.h>
#include <ctime>
#include <fstream>
#include <string>
#include <ctime>
#include <random>
#include <stdio.h>
#include <iomanip>
#include <locale>
#include <sstream>
#include <string>
#include <cstring>
#include <vector>
#include <cstdlib>
#include <math.h>
#include <chrono>
#include <queue>
#include <stdexcept>
#include <atomic>
#include <mutex>
#include <thread>
#include <unordered_set>
#include <unordered_map>
#include <algorithm>

using namespace std;

typedef unsigned int uint;
typedef unsigned long long ull;

// SF: Keeps the destination vertex and label of an edge
struct OutEdge
{
	uint end;
	uint label;

	OutEdge() {}

	OutEdge(uint end, uint label)
	{
		this->end = end;
		this->label = label;
	}
};

// SF: Keeps the source, destination vertex, and label of an edge
struct Edge
{
	uint source;
	OutEdge out;

	Edge() {}

	Edge(uint source, OutEdge out)
	{
		this->source = source;
		this->out = out;
	}
};

// SF:: Keeps the source, destination vertex and label of an edge
// For label the original label string is not saved. Instead the value from
// hashSym map is saved, hashMap[label] = value;
struct EdgeForReading
{
	uint from;
	uint to;
	uint label;

	EdgeForReading() {}

	EdgeForReading(uint from, uint to, uint label)
	{
		this->from = from;
		this->to = to;
		this->label = label;
	}
};

// SF: Keeps the source and destination vertex of an edge
struct UnweightedEdge
{
	uint source;
	uint end;

	UnweightedEdge() {}

	UnweightedEdge(uint source, uint end)
	{
		this->source = source;
		this->end = end;
	}
};

#define COMBINE(i, j) ((ull)i << 32 | j & 0xFFFFFFFFULL)
#define LEFT(i) ((int)(i >> 32))
#define RIGHT(i) ((int)i)

// prints meta-data if true
bool debug = false;

// SF: OLD: old edge, NEW: new edge, FUTURE: future edge
enum
{
	OLD,
	NEW,
	FUTURE,
	OLD_UPD
};

// SF:: counts the total no of unique edges in a 2D hash set
ull countEdge(unordered_set<ull> *hashset, uint hashsetSize)
{
	ull size = 0;
	for (uint i = 0; i < hashsetSize; i++)
	{
		size += hashset[i].size();
	}
	return size;
}

// holds the OLD, NEW, FUTURE edges
struct Buffer
{
	// OLD_START = 0;
	// NEW_START = OLD_END;
	// FUTURE_START = NEW_END;
	// FUTURE_END = vertexList.size()
	uint OLD_END = 0;
	uint NEW_END = 0;
	vector<uint> vertexList;
	// this srcV fields keeps the srcV, it works as indexToSrcMap
	int srcV = -1;

	Buffer() {}
};

// holds the OLD, NEW, FUTURE edges
struct EdgeListBuffer
{
	// This vector's index is the index of the srcV and it holds the destination vertices
	// along with the NEW/OLD/FUTURE pointers
	// actual srcV is saved inside the Buffer struct
	vector<Buffer> edgeList;
	// This vector keeps the srcIndexes that have new edges in it
	vector<int> newEdgeIndexList;

	EdgeListBuffer() {}
};

struct hashFunction
{
	int operator()(const vector<uint> &V) const
	{
		uint hash = V.size();
		for (auto &i : V)
		{
			hash ^= i + 0x9e3779b9 + (hash << 6) + (hash >> 2);
		}
		return hash;
	}
};

enum Color
{
	WHITE,
	GREY,
	BLACK
};

#endif
