#include <omp.h>

#include <cmath>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include "core/data.h"
#include "core/globals.h"
#include "core/io.h"
#include "ea/ea_constants.h"
#include "ea/ea_crossover.h"
#include "ea/ea_data.h"
#include "ea/ea_globals.h"
#include "ea/ea_io.h"
#include "ea/ea_islands.h"
#include "ea/ea_logging.h"
#include "ea/ea_mutations.h"
#include "ea/ea_random.h"
#include "ea/ea_scoring.h"
#include "ea/ea_util.h"
#include "ea_shared/random.h"

#include "json.hpp"

using json = nlohmann::json;

int queue_size;
Scorer* scorer;

void handle_mutations(std::vector<island>& islands, int k, int origin, int destination,
                      int iteration, json& iteration_json) {
#pragma omp parallel for collapse(2) default(none)                                      \
    shared(k, origin, destination, inds, islandsCount, islands, scorer, iteration_json, \
           std::wcout, iteration, std::cout)
  for (int i = 0; i < islandsCount; i++) {
    for (int j = 0; j < inds; j++) {
      island& _island = islands[i];
      std::vector<individual>& offspring = _island.offspring;
      std::vector<individual>& parents = _island.parents;
      offspring[j] = parents[j];

      try {
        mutate(offspring[j], k, origin, destination, j, iteration_json, iteration,
               _island.iterationsWithoutChange, _island.lastExchanges, i);
      } catch (std::exception const& e) {
        std::wstringstream ss;
        ss << L"An error of type '" << typeid(e).name() << L"' occurred while mutating individual "
           << j << " in island " << i << std::endl;
        std::wcout << ss.str() << std::endl;
        throw;
      }
      auto mutate_time = std::chrono::steady_clock::now();
      auto score_result = scorer->scoreIndividual(k, offspring[j]);
      auto scoring_time = std::chrono::steady_clock::now();
      inform_log_about_mutations(iteration_json, score_result, scoring_time, mutate_time);
    }
  }
}

void do_ea(int origin, int destination, const std::vector<int>& pids) {
  int k = pids.size();
  SCORING_METHOD _scorer;
  scorer = &_scorer;

  parse_env();
  setup_logs();

  std::vector<island> islands(islandsCount, island(inds, n));
  initialize_islands(islands, k, origin, destination);

  queue_size = determineQueueSize(islands[0].parents[0].routes, k);
  // do_debug_stuff();

  score_islands(islands, k);
  sort_islands(islands);

  individual* best = getBestIndividual(islands);
  int best_score = static_cast<int>(std::round(best->score));
  assign_routes(best->routes, best->u);
  outputPlansToNewFile(get_intermediate_filename(0, best_score).c_str());

  inform_about_settings(k, origin, destination, best_score);

  int iteration = 1;
  int overallIterationsWithoutChange = 0;
  log_fitness(islands, 0);

  while (overallIterationsWithoutChange <= iterations_needed_for_convergeance && !interrupted) {
    print_iteration_infos(iteration, islands);
    json iteration_json = setup_iteration_json(iteration, best_score);

    if (islandsCount > 1 && inds > 1) {
      do_migrations(islands, iteration, overallIterationsWithoutChange, iteration_json);
    }

    if (inds > 1 && n > 1 && crossover_func != nullptr) {
      if (crossisland_iterations > 0 && iteration % crossisland_iterations == 0 && islandsCount > 1) {
        handle_crossisland_crossover(islands, iteration_json, k);
      } else {
#pragma omp parallel for default(none) \
    shared(islands, iteration_json, islandsCount, k, n, crossover_func)
        for (int i = 0; i < islandsCount; i++) {
          handle_innerisland_crossover(islands[i], iteration_json, k);
        }
      }
    }

    auto mut = std::chrono::steady_clock::now();
    handle_mutations(islands, k, origin, destination, iteration, iteration_json);
    inform_log_about_offspring_scores(iteration_json, islands);
    auto update = std::chrono::steady_clock::now();

    auto update_islands_result = update_islands(islands, overallIterationsWithoutChange);
    best = update_islands_result.first;

    if (update_islands_result.second != 0) {
      std::cout << std::endl;
    }
    auto end_up = std::chrono::steady_clock::now();

    log_fitness(islands, iteration);
    inform_log_about_newparent_scores(iteration_json, islands);

    int new_score = static_cast<int>(std::round(best->score));
    if (best_score != new_score) {
      assign_routes(best->routes, best->u);
      outputPlansToNewFile(get_intermediate_filename(iteration, new_score).c_str());
    }

    best_score = new_score;
    finalize_iteration_log(iteration_json, end_up, mut, update);

    iteration++;
  }

  std::cout << std::endl
            << "Converged after " << iteration - overallIterationsWithoutChange << " Iterations!"
            << std::endl;
  std::cout << "Assigning routes" << std::endl;
  assign_routes(best->routes, best->u);

  std::cout << "Dumping logs" << std::endl;
  close_logs(iteration, overallIterationsWithoutChange);
  std::cout << "Done." << std::endl;
}

// NOLINTNEXTLINE(cppcoreguidelines-avoid-c-arrays,hicpp-avoid-c-arrays,modernize-avoid-c-arrays)
void do_routing(int argc, char* argv[]) {
  (void)argc;
  (void)argv;
  std::cout.precision(std::numeric_limits<double>::max_digits10);
  setup_signal_handling();
  queue_size = minimum_queue_size;

  std::map<std::pair<std::pair<int, int>, std::string>, std::vector<int>>
      c;  // cluster map -> it maps OD pairs and timestring s to person ids that go from O to D at
          // s
  for (unsigned int pid = 0; pid < persons.size(); pid++) {
    auto& p = persons[pid];
    c[{{p.origin, p.destination}, p.timestr}].push_back(pid);
  }

  // do EA for each OD-time-triple
  for (auto& [sdts, pv] : c) {
    do_ea(sdts.first.first, sdts.first.second, pv);
  }

  clearGSLRandom();
}
