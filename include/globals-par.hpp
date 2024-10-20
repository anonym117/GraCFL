#ifndef GLOBALSPAR_HPP
#define GLOBALSPAR_HPP

#include "tbb/concurrent_vector.h"
#include "tbb/concurrent_queue.h"
#include "tbb/concurrent_unordered_set.h"

ull countEdge(tbb::concurrent_unordered_set<ull> **hashset, uint hashsetSize, uint grammarLabelSize)
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

#endif
