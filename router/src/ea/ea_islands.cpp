#include "ea/ea_islands.h"

#include <omp.h>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <iostream>
#include <random>
#include <tuple>
#include <utility>
#include <vector>

#include "ea/ea_data.h"
#include "ea/ea_globals.h"
#include "ea/ea_logging.h"
#include "ea/ea_random.h"
#include "ea/ea_util.h"

#include "json.hpp"

using json = nlohmann::json;

void score_islands(std::vector<island>& islands, int k) {
#pragma omp parallel for default(none) collapse(2) shared(islands, k, scorer, inds, islandsCount)
  for (int j = 0; j < islandsCount; j++) {  // NOLINT(modernize-loop-convert)
    for (int i = 0; i < inds; i++) {        // NOLINT(modernize-loop-convert)
      scorer->scoreIndividual(k, islands[j].parents[i]);
    }
  }
}

void sort_islands(std::vector<island>& islands) {
#pragma omp parallel for default(none) shared(islands)
  for (unsigned int j = 0; j < islands.size(); j++) {  // NOLINT(modernize-loop-convert)
    std::sort(islands[j].parents.begin(), islands[j].parents.end());
  }
}

void initialize_individuals(std::vector<individual>& individuals, int k, int origin,
                            int destination) {
#pragma omp taskloop collapse(2) default(none) shared(k, origin, destination, individuals, inds, n)
  for (int i = 0; i < inds; i++) {
    for (int j = 0; j < n; j++) {
      individuals[i].routes[j] = randDijkstra(k / static_cast<double>(n), origin, destination);
    }
  }
}

void initialize_islands(std::vector<island>& islands, int k, int origin, int destination) {
#pragma omp parallel for default(none) shared(k, origin, destination, inds, n, islands)
  for (unsigned int m = 0; m < islands.size(); m++) {  // NOLINT(modernize-loop-convert)
    initialize_individuals(islands[m].parents, k, origin, destination);
  }
}

void do_migrations(std::vector<island>& islands, int iteration, int& iterationsWithoutChange,
                   json& iteration_json) {
  int topIndivCount = std::min(3, (islandsCount - 1) * inds);

  for (unsigned int i = 0; i < islands.size(); i++) {
    island& _island = islands[i];
    if ((_island.iterationsWithoutChange >= iterations_for_migration) &&
        ((iteration - _island.lastMigration) >= iterations_for_migration)) {
      // finde die min(3, islands * population size) besten Individuuen aus den anderen Inseln
      auto start = std::chrono::steady_clock::now();

      std::vector<std::tuple<double, int, int>> individuals(inds * (islandsCount - 1));
      int idx = 0;
      for (unsigned int k = 0; k < islands.size(); k++) {
        if (k != i) {
          for (int l = 0; l < inds; l++) {
            individual& _indiv = islands[k].parents[l];
            individuals.at(idx) = std::make_tuple(_indiv.score, k, l);
            idx++;
          }
        }
      }
      int scoreBefore = static_cast<int>(std::round(_island.parents[0].score));
      std::sort(individuals.begin(), individuals.end());
      std::vector<double> scores;

      scores.reserve(topIndivCount);
      for (int k = 0; k < topIndivCount; k++) {
        scores.push_back(std::get<0>(individuals.at(k)));
      }
      auto selection = std::chrono::steady_clock::now();

      // wähle gewichtet nach ihrer fitness ein Individuum aus
      std::discrete_distribution<> indivd(scores.begin(), scores.begin() + topIndivCount);
      int indiv = indivd(getGenerator());
      int island_id = std::get<1>(individuals[indiv]);
      int indiv_id = std::get<2>(individuals[indiv]);
      // ersetze das schlechteste individuum aus _island damit
      // ODER (TODO): Erstelle Crossover aus dem besten von _island und dem gewählten und schau ob
      // das besser ist als ??
      auto selection2 = std::chrono::steady_clock::now();
      _island.parents[inds - 1] = islands[island_id].parents[indiv_id];

      // restore sorted invariant and set important variables
      _island.lastMigration = iteration;
      std::sort(_island.parents.begin(), _island.parents.end());
      std::string improvement_str = "no";
      if (scoreBefore == static_cast<int>(std::round(_island.parents[0].score))) {
        _island.iterationsWithoutChange++;
      } else {
        std::cout << std::endl;
        iterationsWithoutChange = 0;
        _island.iterationsWithoutChange = 0;
        improvement_str = "yes";
      }
      auto end = std::chrono::steady_clock::now();
      inform_log_about_migration(iteration_json, i, island_id, indiv_id, improvement_str, start,
                                 selection, selection2, end);
    }
  }
}

// returns pointer to the best individual
auto update_parents(std::vector<individual>& parents, std::vector<individual>& offspring,
                    std::vector<individual>* crossover = nullptr) -> individual* {
  // as a first step, sort the two vectors (in parallel)
  // todo parallelisieren über taskloop
  std::sort(parents.begin(), parents.end());
  std::sort(offspring.begin(), offspring.end());

  if (crossover != nullptr) {
    std::sort((*crossover).begin(), (*crossover).end());
  }

  // now merge the two vectors
  std::vector<individual> concatenated(2 * inds, individual(n));

  if (crossover != nullptr) {
    int cr_score = static_cast<int>(std::round((*crossover)[0].score));
    int bp_score = static_cast<int>(std::round(parents[0].score));
    int bp2_score = static_cast<int>(std::round(parents[0].score));
    bool improved_against_parents = false;
    /*for (auto& indiv : parents) {
      if (best_crossover_indiv.score < indiv.score) {
        improved_against_parents = true;
      }
    } */

    if (cr_score < bp_score || cr_score < bp2_score) {
      improved_against_parents = true;
    }
    //std::cout << "Best Parent: " << bp_score << " Best CR: " << cr_score << std::endl;
    if (improved_against_parents) {
      //std::cout << "Crossover made improvement, taking it!" << std::endl;
      // just if we improved against parents, we include the crossover
      concatenated.resize(2 * inds + static_cast<int>((*crossover).size()), individual(n));
      std::vector<individual> concatenated2(2 * inds, individual(n));
      std::merge(parents.begin(), parents.end(), offspring.begin(), offspring.end(),
                 concatenated2.begin());
      std::merge(concatenated2.begin(), concatenated2.end(), (*crossover).begin(),
                 (*crossover).end(), concatenated.begin());
    } else {
      //std::cout << "not taking crossover!" << std::endl;

      std::merge(parents.begin(), parents.end(), offspring.begin(), offspring.end(),
                 concatenated.begin());
    }

  } else {
    std::merge(parents.begin(), parents.end(), offspring.begin(), offspring.end(),
               concatenated.begin());
  }

  // move the first inds individuals from concatenated to parents
  std::move(concatenated.begin(), concatenated.begin() + inds, parents.begin());

  return &parents[0];
}

auto update_islands(std::vector<island>& islands, int& iterationsWithoutChange)
    -> std::pair<individual*, int> {
  int change = 0;
#pragma omp parallel for default(none) \
    shared(islands, islandsCount, change, crossover_func, std::cout)
  for (int i = 0; i < islandsCount; i++) {  // NOLINT(modernize-loop-convert)
    int scoreBefore = static_cast<int>(std::round(islands[i].parents[0].score));
    if (crossover_func == nullptr) {
      update_parents(islands[i].parents, islands[i].offspring);
    } else {
      update_parents(islands[i].parents, islands[i].offspring, &islands[i].crossover_offspring);
    }
    if (scoreBefore == static_cast<int>(std::round(islands[i].parents[0].score))) {
      islands[i].iterationsWithoutChange++;
    } else {
#pragma omp critical(change_update_island)
      change = 1;
      islands[i].iterationsWithoutChange = 0;
    }
  }

  if (change != 0) {
    iterationsWithoutChange = 0;
  } else {
    iterationsWithoutChange++;
  }

  individual* best = getBestIndividual(islands);

  return std::make_pair(best, change);
}