#include "ea/ea_scoring.h"

#include <algorithm>
#include <iostream>
#include <unordered_map>
#include <vector>

#include <gsl/gsl_poly.h>

#include "core/data.h"
#include "ea/ea_data.h"
#include "ea/ea_util.h"

#include "json.hpp"

using json = nlohmann::json;

void FrankWolfeScorer::calcInitialUsage(int k, individual& indiv) {
  indiv.u.edgeFlow.clear();
  indiv.u.routeFlow.clear();
  indiv.u.routeFlow.resize(indiv.routes.size(), 0);

  indiv.oldUsage.edgeFlow.clear();
  indiv.oldUsage.routeFlow.clear();
  indiv.oldUsage.routeFlow.resize(indiv.routes.size(), 0);

  for (auto& r : indiv.routes) {
    for (link* l : r.links) {
      indiv.u.edgeFlow[l] = 0;
    }
  }

  std::vector<double> latencies;
  for (auto& r : indiv.routes) {
    latencies.push_back(route_latency(indiv.u.edgeFlow, r));
  }
  int chosenRouteIndex = std::min_element(latencies.begin(), latencies.end()) - latencies.begin();

  for (link* l : indiv.routes[chosenRouteIndex].links) {
    indiv.u.edgeFlow[l] = k;
  }

  indiv.u.routeFlow[chosenRouteIndex] = k;
}

void FrankWolfeScorer::calcEdgeFlow(individual& indiv) {
  indiv.u.edgeFlow.clear();
  for (unsigned int i = 0; i < indiv.routes.size(); i++) {
    for (link* l : indiv.routes[i].links) {
      indiv.u.edgeFlow[l] += indiv.u.routeFlow[i];
    }
  }
}

void FrankWolfeScorer::captureOldUsage(individual& indiv) {
  for (auto& p : indiv.u.edgeFlow) {
    indiv.oldUsage.edgeFlow[p.first] = p.second;
  }

  for (unsigned int i = 0; i < indiv.routes.size(); i++) {
    indiv.oldUsage.routeFlow[i] = indiv.u.routeFlow[i];
  }
}

double FrankWolfeScorer::lineSearch(unordered_map<link*, double>& x,
                                    unordered_map<link*, double>& y, double min, double max) {
  // this currently assumes quadratic cost functions
  // see paper/bachelor thesis for what this solves

  double A = 0, B = 0, C = 0;
  for (auto& p : x) {
    double a = p.first->a();
    double b = p.first->b();

    double y_min_x = y[p.first] - x[p.first];
    A += a * y_min_x * y_min_x * y_min_x;
    B += 2 * a * x[p.first] * y_min_x * y_min_x;
    C += (a * x[p.first] * x[p.first] + b) * y_min_x;
  }

  double x0, x1;
  gsl_poly_solve_quadratic(A, B, C, &x0, &x1);

  // check costs on the edge (of the world)
  double zpp_min = 0, zpp_max = 0;
  for (auto& p : x) {
    double a = p.first->a();
    double b = p.first->b();

    double min_u = (1 - min) * p.second + min * y[p.first];
    double max_u = (1 - max) * p.second + max * y[p.first];

    zpp_min += a / 3.0f * min_u * min_u * min_u + b * min_u;
    zpp_max += a / 3.0f * max_u * max_u * max_u + b * max_u;
  }

  if (x0 >= min && x0 <= max)
    return x0;
  if (x1 >= min && x1 <= max)
    return x1;

  return (zpp_min < zpp_max) ? min : max;
}

double FrankWolfeScorer::findGamma(int chosenRouteIndex, int k, individual& indiv) {
  unordered_map<link*, double> y;
  for (auto& p : indiv.u.edgeFlow) {
    y[p.first] = 0;  // just necessary to have a map entry for *all* used edges, even if the new
                     // route does not use them we can subtract in the line search
  }
  for (link* l : indiv.routes[chosenRouteIndex].links) {
    y[l] = k;
  }

  return lineSearch(indiv.u.edgeFlow, y, 0, 1);
}

double FrankWolfeScorer::findMu(unordered_map<link*, double>& x, unordered_map<link*, double>& y,
                                int k) {
  // probably wrong currently
  double max_mu = -1;
  for (auto& p : x) {
    double diff = x[p.first] - y[p.first];
    double mu = 0;
    if (diff == 0)
      continue;
    else if (diff < 0) {
      // when do we hit 0?
      mu = -x[p.first] / diff;
    } else {
      // when do we hit k?
      mu = (static_cast<double>(k) - x[p.first]) / diff;
    }

    if (max_mu == -1 || mu < max_mu)
      max_mu = mu;
  }

  double lambda = lineSearch(x, y, -max_mu, 0);
  return -lambda;
}

std::tuple<double, json> FrankWolfeScorer::scoreIndividual(int k, individual& indiv) {
  // calculate starting point x0 for Frank-Wolfe
  // we set the starting point to sending all drivers over the dijkstra route
  // std::cout << "Called a FrankWolfeScorer!" << std::endl;
  calcInitialUsage(k, indiv);

  int iterationCounter = 0;
  int lastRSME = -1;
  double rsmeEqualCounter = 0;
  int rsmeBelowTenCounter = 0;

  json _log = json::object();
  std::vector<double> latencies2;
  while (rsmeBelowTenCounter < iterations_RSMEBelowTen &&
         rsmeEqualCounter < iterations_RSMENotChanging && iterationCounter < hardlimit) {
    iterationCounter++;

    if (partanEnabled) {
      captureOldUsage(indiv);
    }

    // find point x1 to interpolate to - this is the fastest route assuming the current flow
    // distribution why is this the fastest route? read our paper!

    std::vector<double> latencies;
    for (auto& r : indiv.routes) {
      double l = route_latency(indiv.u.edgeFlow, r);
      latencies.push_back(l);
    }
    unsigned int chosenRouteIndex =
        std::min_element(latencies.begin(), latencies.end()) - latencies.begin();

    // distribute drivers_at_once drivers over the chosen route
    // interpolate

    // double gamma = 2.0f / ((double)cc + 2.0f); for very simple gamma calculation
    double gamma = findGamma(chosenRouteIndex, k, indiv);

    for (unsigned int i = 0; i < indiv.routes.size(); i++) {
      if (i == chosenRouteIndex) {
        indiv.u.routeFlow[i] = (1 - gamma) * indiv.u.routeFlow[i] + gamma * k;
      } else {
        indiv.u.routeFlow[i] = (1 - gamma) * indiv.u.routeFlow[i];
      }
    }

    calcEdgeFlow(indiv);

    if (partanEnabled) {
      double mu = findMu(indiv.u.edgeFlow, indiv.oldUsage.edgeFlow, k);
      for (unsigned int i = 0; i < indiv.routes.size(); i++) {
        indiv.u.routeFlow[i] =
            indiv.u.routeFlow[i] + mu * (indiv.u.routeFlow[i] - indiv.oldUsage.routeFlow[i]);
      }

      calcEdgeFlow(indiv);
    }

    // store latencies again - just to be able to calculate the rsme...

    latencies2.clear();
    for (auto& r : indiv.routes) {
      double l = route_latency(indiv.u.edgeFlow, r);
      latencies2.push_back(l);
    }

    double mean = std::accumulate(latencies2.begin(), latencies2.end(), 0.0F) /
                  static_cast<float>(latencies2.size());

    std::vector<double> errors(latencies2.size());
    std::transform(latencies2.begin(), latencies2.end(), errors.begin(),
                   [&mean](double lat) { return std::pow(lat - mean, 2); });

    if (iterationCounter > 50000) {
      double mu = 1337;
      std::cout << "mean: " << mean << "gamma: " << gamma << " mu: " << mu << std::endl;

      for (unsigned int i = 0; i < latencies2.size(); i++) {
        std::cout << " latency: " << latencies2[i] << " error: " << errors[i]
                  << " flow: " << indiv.u.routeFlow[i];
      }
    }
    double rsme = std::sqrt(std::accumulate(errors.begin(), errors.end(), 0.0F) /
                            static_cast<float>(errors.size()));

    if (rsme < 10) {
      rsmeBelowTenCounter++;
    } else {
      rsmeBelowTenCounter = 0;
    }

    if (iterationCounter > 50000) {
      std::cout.precision(std::numeric_limits<double>::max_digits10);
      auto Compare = [](float a, float b, float epsilon = std::numeric_limits<float>::epsilon()) {
        return (std::fabs(a - b) <= epsilon);
      };
      std::cout << " rmse after iteration " << iterationCounter << " :" << rsme
                << " last rsme: " << lastRSME << " equal counter: " << rsmeEqualCounter
                << " are they equal? " << (lastRSME == rsme) << " are they ccequal lambda? "
                << Compare(rsme, lastRSME) << std::endl;
    }

    if (static_cast<int>(std::round(rsme)) == lastRSME) {
      rsmeEqualCounter++;
    } else {
      rsmeEqualCounter = 0;
      lastRSME = static_cast<int>(std::round(rsme));
    }
  }

  // output FW iteration count convergeance statistic here
  // also output quality of UE to compare with SimulationScorer
  _log["hardlimit_reached"] = 0;
  if (iterationCounter >= iterations_hardlimit) {
    std::cout << "Warning! Reached Hardlimit for Scoring." << std::endl;
    _log["hardlimit_reached"] =
        1;  // todo: maybe if this happens too often then do simulation scoring
  }
  _log["fw_iterations"] = iterationCounter;
  _log["latencies"] = json::array();
  for (double l : latencies2) {
    _log["latencies"].push_back(l);
  }

  double _score = 0;
  for (auto& p : indiv.u.edgeFlow) {
    _score += p.second * p.first->latency(p.second);
  }

  indiv.score = _score;  // set class variable
  return std::make_tuple(_score, _log);
}

void SimulationScorer::calculateUsage(int k, individual& indiv) {
  // initialize edge and route flow

  indiv.u.edgeFlow.clear();
  indiv.u.routeFlow.clear();
  indiv.u.routeFlow.resize(indiv.routes.size(), 0);

  for (auto& r : indiv.routes) {
    for (link* l : r.links) {
      indiv.u.edgeFlow[l] = 0;
    }
  }

  // calculate flow using user equilibrium as psychological model
  for (int i = 0; i < k; i += driver_approximation_factor) {
    // for each drivers_at_once, choose the fastet route given the current edge flows
    std::vector<double> latencies;
    for (auto& r : indiv.routes) {
      latencies.push_back(route_latency(indiv.u.edgeFlow, r));
    }
    int chosenRouteIndex = std::min_element(latencies.begin(), latencies.end()) - latencies.begin();

    // distribute drivers_at_once drivers over the chosen route
    for (link* l : indiv.routes[chosenRouteIndex].links) {
      indiv.u.edgeFlow[l] += driver_approximation_factor;
    }
    indiv.u.routeFlow[chosenRouteIndex] += driver_approximation_factor;
  }
}

std::tuple<double, json> SimulationScorer::scoreIndividual(int k, individual& indiv) {
  calculateUsage(k, indiv);
  json _log = json::object();
  _log["latencies"] = json::array();

  double _score = 0;
  for (auto& p : indiv.u.edgeFlow) {
    _score += p.second * p.first->latency(p.second);
  }

  indiv.score = _score;  // set class variable

  // output equilibium quality
  for (auto& r : indiv.routes) {
    double l = route_latency(indiv.u.edgeFlow, r);
    _log["latencies"].push_back(l);
  }

  return std::make_tuple(_score, _log);
}