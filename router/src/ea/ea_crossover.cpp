#include <chrono>
#include <cmath>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <tuple>
#include <vector>

#include "ea/ea_data.h"
#include "ea/ea_globals.h"
#include "ea/ea_logging.h"
#include "ea/ea_random.h"
#include "ea/ea_scoring.h"
#include "ea/ea_util.h"

#include "json.hpp"

using json = nlohmann::json;

individual makeDiverseBaby(individual& parent1, individual& parent2, int k) {
  // auto start = std::chrono::steady_clock::now();
  individual baby(n);
  FrankWolfeScorer fastScorer(50);
  std::string bitmask(n, 1);  // K leading 1's
  bitmask.resize(2 * n, 0);   // N-K trailing 0's

  std::vector<route> currRoutes(n);
  int curScr = std::numeric_limits<int>::max();  // want to minimize this
  int permut = 0;
  do {
    int next = 0;
    for (int i = 0; i < 2 * n; ++i)  // [0..N-1] integers
    {
      if (bitmask[i]) {
        individual* _src = nullptr;
        int _id = i;
        if (i < n) {
          _src = &parent1;
        } else {
          _src = &parent2;
          _id -= n;
        }

        // std::cout << "Setting route " << (next + 1) << " to indiv " << _id << " of" << (_src ==
        // &parent1 ? " parent 1" : " parent 2") << std::endl;
        currRoutes[next] = _src->routes[_id];
        next++;
      }
    }
    individual _tmp(n);
    _tmp.routes = currRoutes;
    int scr = static_cast<int>(std::round(std::get<0>(fastScorer.scoreIndividual(k, _tmp))));

    if (scr < curScr) {
      // if (permut > 0) std::cout << "New Score better in permutation " << permut << "Old: " <<
      // curScr << " New: " << scr << std::endl;
      baby.routes = currRoutes;
      curScr = scr;
    }

    currRoutes.clear();
    currRoutes.resize(n);
    permut++;
  } while (std::prev_permutation(bitmask.begin(), bitmask.end()));

  return baby;
}

individual makeDiverseBaby2(individual& parent1, individual& parent2, int k) {
  (void)k;
  // auto start = std::chrono::steady_clock::now();
  individual baby(n);
  std::string bitmask(n, 1);  // K leading 1's
  bitmask.resize(2 * n, 0);   // N-K trailing 0's

  std::vector<route> currRoutes(n);
  double currentDistScore = std::numeric_limits<double>::max();  // want to minimize this
  int permut = 0;
  do {
    int next = 0;
    for (int i = 0; i < 2 * n; ++i)  // [0..N-1] integers
    {
      if (bitmask[i]) {
        individual* _src = nullptr;
        int _id = i;
        if (i < n) {
          _src = &parent1;
        } else {
          _src = &parent2;
          _id -= n;
        }

        // std::cout << "Setting route " << (next + 1) << " to indiv " << _id << " of" << (_src ==
        // &parent1 ? " parent 1" : " parent 2") << std::endl;
        currRoutes[next] = _src->routes[_id];
        next++;
      }
    }
    double _score = distScoreRoutes(currRoutes);
    individual _tmp(n);
    _tmp.routes = currRoutes;

    if (_score < currentDistScore) {
      // if (permut > 0) std::cout << "New Score better in permutation " << permut << "Old: " <<
      // currentDistScore << " New: " << _score << std::endl;
      baby.routes = currRoutes;
      currentDistScore = _score;
    }

    currRoutes.clear();
    currRoutes.resize(n);
    permut++;
  } while (std::prev_permutation(bitmask.begin(), bitmask.end()));

  return baby;
}

individual makeHeuristicBaby(individual& parent1, individual& parent2, int k) {
  (void)k;

  individual baby(n);
  std::vector<double> weights(2 * n, 1);
  for (unsigned int i = 0; i < parent1.routes.size(); i++) {
    weights[i] = 1.0F / static_cast<double>(parent1.u.routeFlow[i]);
  }
  for (unsigned int i = 0; i < parent2.routes.size(); i++) {
    weights[i + n] = 1.0F / static_cast<double>(parent2.u.routeFlow[i]);
  }
  std::discrete_distribution<> d(weights.begin(), weights.end());
  int start_route = d(getGenerator());

  if (start_route < n) {
    baby.routes[0] = parent1.routes[start_route];
  } else {
    baby.routes[0] = parent2.routes[start_route - n];
  }

  for (unsigned int i = 1; i < baby.routes.size(); i++) {
    double best_dist = std::numeric_limits<double>::max();
    int best = -1;
    individual* _src = nullptr;
    for (unsigned int j = 0; j < parent1.routes.size(); j++) {
      std::vector<route> _tmproutes(baby.routes);
      _tmproutes.push_back(parent1.routes[j]);
      double score = distScoreRoutes(_tmproutes);

      if (score < best_dist) {
        best_dist = score;
        best = j;
        _src = &parent1;
      }
    }

    for (unsigned int j = 0; j < parent2.routes.size(); j++) {
      std::vector<route> _tmproutes(baby.routes);
      _tmproutes.push_back(parent2.routes[j]);
      double score = distScoreRoutes(_tmproutes);

      if (score < best_dist) {
        best_dist = score;
        best = j;
        _src = &parent2;
      }
    }

    if (best == -1) {
      std::cout << "holy shit the house is on fire, couldnt determine best individual greedily."
                << std::endl;
      baby = parent1;
      return baby;
    }
    baby.routes[i] = _src->routes[best];
  }
  return baby;
}

individual makeHeuristicBaby2(individual& parent1, individual& parent2, int k) {
  (void)k;

  individual baby(n);
  std::vector<double> weights(2 * n, 1);
  for (unsigned int i = 0; i < parent1.routes.size(); i++) {
    weights[i] = 1.0F / static_cast<double>(parent1.u.routeFlow[i]);
  }
  for (unsigned int i = 0; i < parent2.routes.size(); i++) {
    weights[i + n] = 1.0F / static_cast<double>(parent2.u.routeFlow[i]);
  }
  std::discrete_distribution<> d(weights.begin(), weights.end());
  int start_route = d(getGenerator());

  if (start_route < n) {
    baby.routes[0] = parent1.routes[start_route];
  } else {
    baby.routes[0] = parent2.routes[start_route - n];
  }

  for (unsigned int i = 1; i < baby.routes.size(); i++) {
    std::vector<double> scores;
    scores.reserve(2 * n);
    for (unsigned int j = 0; j < parent1.routes.size(); j++) {
      std::vector<route> _tmproutes(baby.routes);
      _tmproutes.push_back(parent1.routes[j]);
      double score = distScoreRoutes(_tmproutes);

      scores.push_back(1.0F / score);
    }

    for (unsigned int j = 0; j < parent2.routes.size(); j++) {
      std::vector<route> _tmproutes(baby.routes);
      _tmproutes.push_back(parent2.routes[j]);
      double score = distScoreRoutes(_tmproutes);

      scores.push_back(1.0F / score);
    }

    std::discrete_distribution<> nrouted(scores.begin(), scores.end());
    int new_route = nrouted(getGenerator());
    if (new_route < n) {
      baby.routes[i] = parent1.routes[new_route];
    } else {
      baby.routes[i] = parent2.routes[new_route - n];
    }
  }
  return baby;
}

int binomialCoeff(int _n, int _k) {
  int res = 1;

  // Since C(n, k) = C(n, n-k)
  if (_k > _n - _k)
    _k = _n - _k;

  // Calculate value of
  // [n * (n-1) *---* (n-k+1)] / [k * (k-1) *----* 1]
  for (int i = 0; i < _k; ++i) {
    res *= (_n - i);
    res /= (i + 1);
  }

  return res;
}

void handle_innerisland_crossover(island& _island, json& iteration_json, int k) {
  auto start = std::chrono::steady_clock::now();
  json crossover_log = json::object();
  crossover_log["type"] = "innerisland";
  _island.crossover_offspring.clear();
  int offspring_amount = static_cast<int>(
      std::round(std::sqrt(binomialCoeff(static_cast<int>(_island.parents.size()), 2))));
  crossover_log["offspring_amount"] = offspring_amount;
  _island.crossover_offspring.resize(offspring_amount, individual(n));

  std::vector<int> weights(static_cast<int>(_island.parents.size()), 1);
  std::discrete_distribution<> d(weights.begin(), weights.end());

  auto cstart = std::chrono::steady_clock::now();
  crossover_log["crossovers"] = json::array();
#pragma omp parallel for default(none) \
    shared(offspring_amount, d, _island, k, n, scorer, crossover_func, crossover_log, std::wcout)
  for (int i = 0; i < offspring_amount; i++) {
    auto it_start = std::chrono::steady_clock::now();
    int parent1 = d(getGenerator());
    int parent2 = -1;

    do {
      parent2 = d(getGenerator());
    } while (parent1 == parent2);

    json tmplog = {{"parent1", parent1}, {"parent2", parent2}};

    _island.crossover_offspring[i] =
        crossover_func(_island.parents[parent1], _island.parents[parent2], k);
    int score = static_cast<int>(
        std::round(std::get<0>(scorer->scoreIndividual(k, _island.crossover_offspring[i]))));
    auto it_end = std::chrono::steady_clock::now();
    tmplog["duration"] =
        std::chrono::duration_cast<std::chrono::microseconds>(it_end - it_start).count();
    tmplog["score"] = score;

#pragma omp critical(crossover_log_access)
    {
      try {
        crossover_log["crossovers"].push_back(tmplog);
      } catch (std::exception const& e) {
        std::wstringstream ss;
        ss << L"An error of type '" << typeid(e).name()
           << L"' occurred while pushing tmplog into crossover_log" << std::endl;
        std::wcout << ss.str() << std::endl;
      }
    }
  }
  auto cend = std::chrono::steady_clock::now();

  crossover_log["duration_crossovers"] =
      std::chrono::duration_cast<std::chrono::microseconds>(cend - cstart).count();
  crossover_log["duration_total"] =
      std::chrono::duration_cast<std::chrono::microseconds>(cend - start).count();
  crossover_log["duration_initializing"] =
      std::chrono::duration_cast<std::chrono::microseconds>(cstart - start).count();

  inform_log_about_crossover_json(iteration_json, crossover_log);
}

void handle_crossisland_crossover(std::vector<island>& islands, json& iteration_json, int k) {
  auto start = std::chrono::steady_clock::now();
  json crossover_log = json::object();
  crossover_log["type"] = "crossisland";
  unsigned int wifes_per_island = static_cast<int>(std::round(islands.size() / 2.0F));
  crossover_log["wifes_per_island"] = wifes_per_island;
  crossover_log["offspring_per_wife"] = offspring_per_wife;

  std::vector<int> weights(static_cast<int>(islands.size()), 1);
  std::discrete_distribution<> d(weights.begin(), weights.end());
  std::vector<std::vector<int>> wifes;
  wifes.resize(islands.size(), std::vector<int>(wifes_per_island));
  auto wife_start = std::chrono::steady_clock::now();
#pragma omp parallel for default(none) shared(islands, wifes, d, wifes_per_island, n)
  for (unsigned int i = 0; i < islands.size(); i++) {
    unsigned int currWife = 0;
    while (currWife < wifes_per_island) {
      unsigned int _wife = d(getGenerator());
      if (_wife != i) {
        wifes[i][currWife] = _wife;
        currWife++;
      }
    }

    islands[i].crossover_offspring.clear();
    islands[i].crossover_offspring.resize(wifes_per_island * offspring_per_wife, individual(n));
  }

  std::vector<int> weights2(static_cast<int>(islands[0].parents.size()), 1);
  std::discrete_distribution<> d2(weights2.begin(), weights2.end());

  crossover_log["crossovers"] = json::array();
  auto cstart = std::chrono::steady_clock::now();

#pragma omp parallel for default(none)                                                            \
    shared(islands, wifes, wifes_per_island, k, n, offspring_per_wife, crossover_func, std::cout, \
           d2, scorer, std::wcout, crossover_log) collapse(3)
  for (unsigned int i = 0; i < islands.size(); i++) {
    for (unsigned int j = 0; j < wifes_per_island; j++) {
      for (int l = 0; l < offspring_per_wife; l++) {
        auto it_start = std::chrono::steady_clock::now();

        int parent1 =
            d2(getGenerator());  // from islands[i].parents[parent1] nehmen wir das eine Elternteil
        int parent2 = -1;

        do {
          parent2 = d2(getGenerator());
        } while (parent1 == parent2);
        json tmplog = {{"island", i},
                       {"wife", wifes[i][j]},
                       {"parent1", parent1},
                       {"parent2", parent2},
                       {"child", j * offspring_per_wife + l}};

        // from islands[wifes[i][j]].parents[parent2] nehmen wir das zweite Elternteil
        // und damit machen wir jetzt ein Kind (und nennen es nach toten Personen)
        /*
        #pragma omp critical
                std::cout << std::endl
                          << "Cross-Island-Crossver. Putting into Island " << i << " the child of
        Island "
                          << i << ", Individual " << parent1 << " and Island " << wifes[i][j]
                          << ", Individual " << parent2 << std::endl;
        */
        islands[i].crossover_offspring.at(j * offspring_per_wife + l) = crossover_func(
            islands[i].parents[parent1], islands[wifes[i][j]].parents[parent2], k);
        int score = static_cast<int>(std::round(std::get<0>(scorer->scoreIndividual(
            k, islands[i].crossover_offspring[j * offspring_per_wife + l]))));
        auto it_end = std::chrono::steady_clock::now();
        tmplog["duration"] =
            std::chrono::duration_cast<std::chrono::microseconds>(it_end - it_start).count();
        tmplog["score"] = score;

#pragma omp critical(crossover_log_access)
        {
          try {
            crossover_log["crossovers"].push_back(tmplog);
          } catch (std::exception const& e) {
            std::wstringstream ss;
            ss << L"An error of type '" << typeid(e).name()
               << L"' occurred while pushing tmplog into crossover_log" << std::endl;
            std::wcout << ss.str() << std::endl;
          }
        }
      }
    }
  }
  auto cend = std::chrono::steady_clock::now();

  crossover_log["duration_crossovers"] =
      std::chrono::duration_cast<std::chrono::microseconds>(cend - cstart).count();
  crossover_log["duration_total"] =
      std::chrono::duration_cast<std::chrono::microseconds>(cend - start).count();
  crossover_log["duration_initializing"] =
      std::chrono::duration_cast<std::chrono::microseconds>(wife_start - start).count();
  crossover_log["duration_wife_calculation"] =
      std::chrono::duration_cast<std::chrono::microseconds>(cstart - wife_start).count();

  inform_log_about_crossover_json(iteration_json, crossover_log);
}