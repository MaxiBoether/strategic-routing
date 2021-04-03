#pragma once

#include <tuple>
#include <unordered_map>

#include "core/data.h"
#include "core/globals.h"
#include "ea/ea_constants.h"
#include "ea/ea_data.h"

#include "json.hpp"

using json = nlohmann::json;

class Scorer {
 public:
  virtual std::tuple<double, json> scoreIndividual(int k, individual& indiv) = 0;
};

class FrankWolfeScorer : public Scorer {
 private:
  const bool partanEnabled = false;
  int hardlimit;
  void calcInitialUsage(int k, individual& indiv);
  void calcEdgeFlow(individual& indiv);
  void captureOldUsage(individual& indiv);

  double lineSearch(unordered_map<link*, double>& x, unordered_map<link*, double>& y, double min,
                    double max);
  double findGamma(int chosenRouteIndex, int k, individual& indiv);
  double findMu(unordered_map<link*, double>& x, unordered_map<link*, double>& y, int k);

 public:
  explicit FrankWolfeScorer() : hardlimit(iterations_hardlimit) {}
  explicit FrankWolfeScorer(int _hardlimit) : hardlimit(_hardlimit) {}
  std::tuple<double, json> scoreIndividual(int k, individual& indiv) override;
};

class SimulationScorer : public Scorer {
 private:
  const int driver_approximation_factor = 1;
  void calculateUsage(int k, individual& indiv);

 public:
  std::tuple<double, json> scoreIndividual (int k, individual& indiv) override;
};