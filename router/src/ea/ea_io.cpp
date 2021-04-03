#include "ea/ea_io.h"

#include <csignal>
#include <filesystem>
#include <iostream>
#include <sstream>

#include "ea/ea_crossover.h"
#include "ea/ea_defaults.h"
#include "ea/ea_globals.h"

#define STRINGIFY(x) #x           // NOLINT(cppcoreguidelines-macro-usage)
#define TOSTRING(x) STRINGIFY(x)  // NOLINT(cppcoreguidelines-macro-usage)

int inds;
int n;
int islandsCount;
int iterations_needed_for_convergeance;
int iterations_for_migration;
int crossisland_iterations;
std::string intermediate_results_prefix;
std::string intermediate_folder;
individual (*crossover_func)(individual&, individual&, int);
std::string crossover_func_str;

bool interrupted = false;
void signal_handler(int signal_num) {
  std::cout << std::endl << "Received interrupt signal (" << signal_num << ")." << std::endl;
  interrupted = true;
}

void setup_signal_handling() {
  // std::signal(SIGABRT, signal_handler);
  // std::signal(SIGTERM, signal_handler);
  std::signal(SIGINT, signal_handler);
}

void parse_env() {
  std::cout << "Parsing the environment." << std::endl;
  char* irp_env = getenv("EA_INTERMEDIATE_RESULTS_PREFIX");
  if (irp_env == nullptr) {
    std::cout << "You may specify the env variable EA_INTERMEDIATE_RESULTS_PREFIX" << std::endl;
    intermediate_results_prefix = "andreaspolze";
  } else {
    intermediate_results_prefix = std::string(irp_env);
  }

  char* imf_env = getenv("EA_INTERMEDIATE_RESULTS_FOLDER");
  if (imf_env == nullptr) {
    std::cout << "You may specify the env variable EA_INTERMEDIATE_RESULTS_FOLDER" << std::endl;
    intermediate_folder = std::filesystem::current_path().string();
  } else {
    intermediate_folder = std::string(imf_env);
  }

  char* indiv_env = getenv("EA_INDIVIDUALS");
  if (indiv_env == nullptr) {
    std::cout << "You may specify the env variable EA_INDIVIDUALS" << std::endl;
    inds = default_inds;
  } else {
    inds = static_cast<int>(std::strtol(indiv_env, nullptr, 0));
  }

  char* islands_env = getenv("EA_ISLANDS");
  if (islands_env == nullptr) {
    std::cout << "You may specify the env variable EA_ISLANDS" << std::endl;
    islandsCount = default_islands;
  } else {
    islandsCount = static_cast<int>(std::strtol(islands_env, nullptr, 0));
  }

  char* routes_env = getenv("EA_ROUTECOUNT");
  if (routes_env == nullptr) {
    std::cout << "You may specify the env variable EA_ROUTECOUNT" << std::endl;
    n = default_n;
  } else {
    n = static_cast<int>(std::strtol(routes_env, nullptr, 0));
  }

  char* conv_env = getenv("EA_CONVERGEANCE_COUNT");
  if (conv_env == nullptr) {
    std::cout << "You may specify the env variable EA_CONVERGEANCE_COUNT" << std::endl;
    iterations_needed_for_convergeance = default_iterations_needed_for_convergeance;
  } else {
    iterations_needed_for_convergeance = static_cast<int>(std::strtol(conv_env, nullptr, 0));
  }

  char* migration_env = getenv("EA_MIGRATION_ITERATIONS");
  if (migration_env == nullptr) {
    std::cout << "You may specify the env variable EA_MIGRATION_ITERATIONS" << std::endl;
    iterations_for_migration = default_iterations_for_migration;
  } else {
    iterations_for_migration = static_cast<int>(std::strtol(migration_env, nullptr, 0));
  }

  char* crossisland_env = getenv("EA_CROSSISLAND_ITERATIONS");
  if (crossisland_env == nullptr) {
    std::cout << "You may specify the env variable EA_CROSSISLAND_ITERATIONS (-1 for disabling it)" << std::endl;
    crossisland_iterations = default_crossisland_iterations;
  } else {
    crossisland_iterations = static_cast<int>(std::strtol(crossisland_env, nullptr, 0));
  }

  char* crossover_env = getenv("EA_CROSSOVER_FUNC");
  if (crossover_env == nullptr) {
    std::cout << "You may specify the env variable EA_CROSSOVER_FUNC (crossover thus disabled)"
              << std::endl;
    crossover_func = nullptr;
    crossover_func_str = "no crossover";
  } else {
    std::string _crossover_func_str(crossover_env);
    crossover_func_str = "Unknown Crossover Method specified!";
    if (_crossover_func_str == "fw") {
      crossover_func_str = "FrankWolfe-Permutation";
      crossover_func = &makeDiverseBaby;
    }
    if (_crossover_func_str == "heur-all") {
      crossover_func_str = "Heuristic-Permutation";
      crossover_func = &makeDiverseBaby2;
    }
    if (_crossover_func_str == "heur-greed") {
      crossover_func_str = "Heuristic-Greedy";
      crossover_func = &makeHeuristicBaby;
    }
    if (_crossover_func_str == "heur-greed-rand") {
      crossover_func_str = "Heuristic-RandomGreedy";
      crossover_func = &makeHeuristicBaby2;
    }
  }
}

auto get_intermediate_filename(int iteration, int score) -> std::string {
  std::stringstream _filename;
  _filename << intermediate_folder << "/" << intermediate_results_prefix << "_" << iteration << "_"
            << score << ".xml";

  return _filename.str();
}

void inform_about_settings(int k, int origin, int destination, int score) {
  std::cout << "You will find the intermediate output in " << intermediate_folder << "/"
            << intermediate_results_prefix << "*.xml" << std::endl;
  std::cout << "Islands: " << islandsCount << std::endl;
  std::cout << "Individuals per Island: " << inds << std::endl;
  std::cout << "Routes: " << n << std::endl;
  std::cout << "Iterations needed for convergeance: " << iterations_needed_for_convergeance
            << std::endl;
  std::cout << "Drivers (k): " << k << std::endl;
  std::cout << "Routing from " << origin << " to " << destination << std::endl;
  std::cout << "Queue Size: " << queue_size << std::endl;
  std::cout << "Initial best score: " << score << std::endl;
  std::cout << "ScoringMethod: " << TOSTRING(SCORING_METHOD) << std::endl;
  std::cout << "CrossoverMethod: " << crossover_func_str << std::endl;
  std::cout << "Stuck Iterations for Migration: " << iterations_for_migration << std::endl;
  std::cout << "Cross-Island Crossover Iterations: " << crossisland_iterations << std::endl;
}

void print_iteration_infos(int iteration, std::vector<island>& islands) {
  std::cout << "\33[2K\r" << std::flush;
  std::cout << "Iteration " << iteration << ": ";
  for (int i = 0; i < islandsCount - 1; i++) {
    std::cout << islands[i].parents[0].score << ", ";
  }
  std::cout << islands[islandsCount - 1].parents[0].score;
  std::cout << std::flush;
}