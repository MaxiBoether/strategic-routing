#pragma once
#include <unordered_map>
#include <vector>

#include "core/data.h"
#include "core/globals.h"

struct usage {
  std::vector<double> routeFlow;
  std::unordered_map<link*, double> edgeFlow;
};

class individual {
 public:
  std::vector<route> routes;
  usage u;
  usage oldUsage;
  double score = -1;
  explicit individual(int n) : routes(n) {}
  bool operator<(const individual& str) const { return (score < str.score); }
};

class island {
  public:
    std::vector<individual> parents;
    std::vector<individual> offspring;
    std::vector<individual> crossover_offspring;
    std::vector<int> lastExchanges;
    int iterationsWithoutChange;
    int lastMigration;
    explicit island(int indivs, int routes) : parents(indivs, individual(routes)), offspring(indivs, individual(routes)), lastExchanges(indivs), iterationsWithoutChange(0), lastMigration(0) {}
};