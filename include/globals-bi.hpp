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
#include <stdint.h>

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
struct EdgeForReading2
{
  uint from;
  uint to;
  uint label;
  uint contextID;

  EdgeForReading2() {}

  EdgeForReading2(uint from, uint to, uint label, int contextID)
  {
    this->from = from;
    this->to = to;
    this->label = label;
    this->contextID = contextID;
  }
};

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
#define COMBINE_INT(i, j, k) ((int64_t)(((int64_t)i << 32) | ((int64_t)j << 8) | (int64_t)k))
#define CANTOR2(a, b) ((a + b + 1) * (a + b) / 2 + b)
#define CANTOR3(a, b, c) (CANTOR2(a, CANTOR2(b, c)))
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
uint countEdge(unordered_set<ull> *hashset, uint hashsetSize)
{
  uint size = 0;
  for (uint i = 0; i < hashsetSize; i++)
  {
    size += hashset[i].size();
  }
  return size;
}

ull countEdge(vector<vector<unordered_set<ull>>> &hashset, uint hashsetSize, uint grammarLabelSize)
{
  ull size = 0;
  for (uint i = 0; i < hashsetSize; i++)
  {
    for (uint j = 0; j < grammarLabelSize; j++)
    {
      size += hashset[i][j].size();
    }
  }
  return size;
}

// holds the OLD, NEW, FUTURE edges
struct Buffer2
{
  // OLD_START = 0;
  // NEW_START = OLD_END;
  // FUTURE_START = NEW_END;
  // FUTURE_END = vertexList.size()
  uint OLD_END = 0;
  uint NEW_END = 0;
  // first: destV, second: contextID
  vector<pair<int, int>> vertexList;
  // this srcV fields keeps the srcV, it works as indexToSrcMap
  int srcV = -1;

  Buffer2() {}
};

struct Buffer
{
  // OLD_START = 0;
  // NEW_START = OLD_END;
  // FUTURE_START = NEW_END;
  // FUTURE_END = vertexList.size()
  uint OLD_END = 0;
  uint NEW_END = 0;

  vector<int> vertexList;
  // this srcV fields keeps the srcV, it works as indexToSrcMap
  int srcV = -1;

  Buffer() {}
};

struct BufferOpt
{
  // OLD_START = 0;
  // NEW_START = OLD_END;
  // FUTURE_START = NEW_END;
  // FUTURE_END = vertexList.size()
  uint OLD_END = 0;
  uint NEW_END = 0;
  uint IN_OLD_END = 0;
  uint IN_NEW_END = 0;
  vector<uint> outVertexList;
  vector<uint> inVertexList;
  // this srcV fields keeps the srcV, it works as indexToSrcMap
  int srcV = -1;

  BufferOpt() {}
};

struct BufferOrder
{

  // pair.first = OLD_END, pair.second = NEW_END
  // index of the vector is the rule ID in the dependency grammar graph
  vector<pair<uint, uint>> temporalPointers;
  // first: destV, second: contextID
  vector<int> vertexList;
  // this srcV fields keeps the srcV, it works as indexToSrcMap
  int srcV = -1;

  BufferOrder() {}
};

// holds the OLD, NEW, FUTURE edges
struct EdgeListBuffer
{
  // This vector's index is the index of the srcV and it holds the destination vertices
  // along with the NEW/OLD/FUTURE pointers
  // actual srcV is saved inside the Buffer struct
  vector<BufferOpt> edgeList;
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

bool hasEnding(std::string const &baseStr, std::string const &endStr)
{
  if (baseStr.length() >= endStr.length())
  {
    return (0 == baseStr.compare(baseStr.length() - endStr.length(), endStr.length(), endStr));
  }
  else
  {
    return false;
  }
}

void getPeakMemoryUsage()
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
      double memGB = memoryKb / (1024.0 * 1024.0);
      std::cout << "# Peak Virtual Memory Usage =\t" << memoryKb << " KB" << std::endl;
      std::cout << "# Peak Virtual Memory Usage (in GB) =\t" << memGB << " GB" << std::endl;
    }

    if (line.substr(0, 6) == "VmHWM:")
    {
      long memoryKb;
      std::istringstream iss(line.substr(6));
      iss >> memoryKb;
      double memGB = memoryKb / (1024.0 * 1024.0);
      std::cout << "# Peak Physical  Memory Usage =\t" << memoryKb << " KB" << std::endl;
      std::cout << "# Peak Physical Memory Usage (in GB) =\t" << memGB << " GB" << std::endl;
    }

    if (line.substr(0, 6) == "VmRSS:")
    {
      long memoryKb;
      std::istringstream iss(line.substr(6));
      iss >> memoryKb;
      double memGB = memoryKb / (1024.0 * 1024.0);
      std::cout << "# VmRSS  Memory Usage =\t" << memoryKb << " KB" << std::endl;
      std::cout << "# VmRSS Memory Usage (in GB) =\t" << memGB << " GB" << std::endl;
    }
  }
}

#endif
