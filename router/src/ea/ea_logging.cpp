#include "ea/ea_logging.h"

#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <utility>

#include "ea/ea_globals.h"

#include "json.hpp"

using json = nlohmann::json;

std::ofstream fitness_log;
std::ofstream iteration_log;
json iteration_log_json;

void setup_logs() {
  std::stringstream fitnesslog_filename;
  fitnesslog_filename << intermediate_folder << "/fitnesslog.csv";
  fitness_log.open(fitnesslog_filename.str());
  fitness_log << "Iteration,";
  for (int i = 0; i < islandsCount; i++) {
    for (int j = 0; j < inds; j++) {
      fitness_log << "(" << i << "," << j << ")"
                  << ",";
    }
  }

  fitness_log << "\n";

  std::stringstream iterationlog_filename;
  iterationlog_filename << intermediate_folder << "/iterationlog.json";
  iteration_log.open(iterationlog_filename.str());
  iteration_log_json = json::array();
}

void close_logs(int iteration, int iterationsWithoutChange) {
  if (iterationsWithoutChange > 0) {
    iteration_log_json.erase(iteration_log_json.begin() + (iteration - iterationsWithoutChange),
                             iteration_log_json.end());
  }

  iteration_log << iteration_log_json.dump(4) << std::endl;
  fitness_log.close();
  iteration_log.close();
}

void log_fitness(std::vector<island>& islands, int iteration) {
  fitness_log << iteration << ",";
  for (auto& _island : islands) {
    for (auto& indiv : _island.parents) {
      fitness_log << std::setprecision(std::numeric_limits<double>::max_digits10)
                  << std::round(indiv.score) << ",";
    }
  }
  fitness_log << "\n";
}

auto setup_iteration_json(int iteration, int score) -> json {
  json iteration_json = {
      {"iteration", iteration},
      {"score", score},
  };

  iteration_json["migrations"] = json::array();
  iteration_json["mutations"] = json::array();
  iteration_json["times"] = json::object();
  iteration_json["times"]["scoring"] = json::array();

  return iteration_json;
}

void inform_log_about_migration(json& iteration_json, int stuck_island, int source_island,
                                int source_individual, const std::string& improvement_str,
                                std::chrono::time_point<std::chrono::steady_clock> start,
                                std::chrono::time_point<std::chrono::steady_clock> selection,
                                std::chrono::time_point<std::chrono::steady_clock> selection2,
                                std::chrono::time_point<std::chrono::steady_clock> end) {
  json migration_json = {
      {"stuck_island", stuck_island},
      {"source_island", source_island},
      {"source_individual", source_individual},
      {"improvement", improvement_str},
      {"times",
       {{"finding_best",
         std::chrono::duration_cast<std::chrono::microseconds>(selection - start).count()},
        {"selecting_random",
         std::chrono::duration_cast<std::chrono::microseconds>(selection2 - selection).count()},
        {"replacing",
         std::chrono::duration_cast<std::chrono::microseconds>(end - selection2).count()}}}};

#pragma omp critical(iteration_json)
  {
    try {
      iteration_json["migrations"].push_back(migration_json);
    } catch (std::exception const& e) {
      // Print the error to the console.
      std::wstringstream ss;
      ss << L"An error of type '" << typeid(e).name()
         << L"' occurred while pushing migration_json into interation_json['migrations']."
         << std::endl;
      std::wcout << ss.str() << std::endl;
    }
  }
}

void inform_log_about_mutations(json& iteration_json, std::tuple<double, json>& score_result,
                                std::chrono::time_point<std::chrono::steady_clock> scoring_time,
                                std::chrono::time_point<std::chrono::steady_clock> mutate_time) {
#pragma omp critical(iteration_json)
  {
    try {
      iteration_json["scoring_info"].push_back(std::get<1>(score_result));
      iteration_json["times"]["scoring"].push_back(
          std::chrono::duration_cast<std::chrono::microseconds>(scoring_time - mutate_time)
              .count());
    } catch (std::exception const& e) {
      std::wstringstream ss;
      ss << L"An error of type '" << typeid(e).name()
         << L"' occurred while pushing scoring time into iteration_json" << std::endl;
      std::wcout << ss.str() << std::endl;
    }
  }
}

void inform_log_about_offspring_scores(json& iteration_json, std::vector<island>& islands) {
#pragma omp critical(iteration_json)
  {
    try {
      for (int i = 0; i < islandsCount; i++) {
        for (int j = 0; j < inds; j++) {
          iteration_json["offspring_scores"]["island" + std::to_string(i)].push_back(
              islands[i].offspring[j].score);
        }
      }
    } catch (std::exception const& e) {
      std::wstringstream ss;
      ss << L"An error of type '" << typeid(e).name()
         << L"' occurred while pushing offspring_scores into iteration_json" << std::endl;
      std::wcout << ss.str() << std::endl;
    }
  }
}

void inform_log_about_newparent_scores(json& iteration_json, std::vector<island>& islands) {
#pragma omp critical(iteration_json)
  {
    try {
      for (int i = 0; i < islandsCount; i++) {
        for (int j = 0; j < inds; j++) {
          iteration_json["new_parents_scores"]["island" + std::to_string(i)].push_back(
              islands[i].parents[j].score);
        }
      }
    } catch (std::exception const& e) {
      std::wstringstream ss;
      ss << L"An error of type '" << typeid(e).name()
         << L"' occurred while pushing newparent_scores into iteration_json" << std::endl;
      std::wcout << ss.str() << std::endl;
    }
  }
}

void finalize_iteration_log(json& iteration_json,
                            std::chrono::time_point<std::chrono::steady_clock> end_up,
                            std::chrono::time_point<std::chrono::steady_clock> mut,
                            std::chrono::time_point<std::chrono::steady_clock> update) {
#pragma omp critical(iteration_json)
  {
    try {
      iteration_json["times"]["total"] =
          std::chrono::duration_cast<std::chrono::microseconds>(end_up - mut).count();
      iteration_json["times"]["individuals_loop"] =
          std::chrono::duration_cast<std::chrono::microseconds>(update - mut).count();
      iteration_json["times"]["parent_update"] =
          std::chrono::duration_cast<std::chrono::microseconds>(end_up - update).count();
    } catch (std::exception const& e) {
      std::wstringstream ss;
      ss << L"An error of type '" << typeid(e).name()
         << L"' occurred while pushing final times into iteration_json" << std::endl;
      std::wcout << ss.str() << std::endl;
    }
  }

#pragma omp critical(iteration_log_json)
  {
    try {
      iteration_log_json.push_back(iteration_json);
    } catch (std::exception const& e) {
      std::wstringstream ss;
      ss << L"An error of type '" << typeid(e).name()
         << L"' occurred while pushing iteration_json into iteration_log_json" << std::endl;
      std::wcout << ss.str() << std::endl;
    }
  }
}

void inform_log_about_mutations_json(json& iteration_json, json& mutations_json) {
#pragma omp critical(iteration_json)
  {
    try {
      iteration_json["mutations"].push_back(mutations_json);
    } catch (std::exception const& e) {
      // Print the error to the console.
      std::wstringstream ss;
      ss << L"An error of type '" << typeid(e).name()
         << L"' occurred while pushing mutations_json into interation_json['mutations']."
         << std::endl;
      std::wcout << ss.str() << std::endl;
    }
  }
}

void inform_log_about_crossover_json(json& iteration_json, json& crossover_json) {
#pragma omp critical(iteration_json)
  {
    try {
      iteration_json["crossover"] = crossover_json;
    } catch (std::exception const& e) {
      // Print the error to the console.
      std::wstringstream ss;
      ss << L"An error of type '" << typeid(e).name()
         << L"' occurred while pushing crossover_json into iteration_json['crossover']."
         << std::endl;
      std::wcout << ss.str() << std::endl;
    }
  }
}