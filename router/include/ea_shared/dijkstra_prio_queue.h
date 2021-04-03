#pragma once
#include <algorithm>
#include <iostream>
#include <list>
#include <utility>
#include <vector>

#include "core/data.h"

class dijkPQ {
 private:
  std::vector<std::list<int>> queue;
  int pointer = 0;
  int max_candidate = 0;

 public:
  explicit dijkPQ(int initialSize) { queue.resize(initialSize, std::list<int>()); }

  auto top() -> std::pair<int, int> {
    while (pointer <= max_candidate) {
      if (!queue[pointer].empty()) {
        return {pointer, queue[pointer].front()};
      }
      pointer++;
    }
    return {-1, -1};
  }

  void pop() {
    while (pointer <= max_candidate) {
      if (!queue[pointer].empty()) {
        queue[pointer].pop_front();
        return;
      }
      pointer++;
    }
  }

  auto empty() -> bool {
    while (pointer <= max_candidate) {
      if (!queue[pointer].empty()) {
        return false;
      }
      pointer++;
    }
    return true;
  }

  auto getSize() -> int { return static_cast<int>(queue.size()); }

  void push(std::pair<int, int> p) {
    // NOLINTNEXTLINE(readability-implicit-bool-conversion, google-runtime-int)
    if (__builtin_expect(p.first >= static_cast<long>(queue.size()), 0)) {
      // std::cout << "resizing dijkstra queue array from " << queue.size() << " to " <<
      // std::max(p.first+1, static_cast<int>(queue.size()*2)) << std::endl;
      queue.resize(std::max(p.first + 1, static_cast<int>(queue.size() * 2)), std::list<int>());
    }
    queue[p.first].push_back(p.second);
    max_candidate = std::max(max_candidate, p.first);
  }
};