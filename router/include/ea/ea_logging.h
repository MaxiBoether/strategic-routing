#pragma once
#include <chrono>
#include <string>
#include <tuple>
#include <vector>

#include "ea/ea_data.h"

#include "json.hpp"

using json = nlohmann::json;

extern void setup_logs();
extern void close_logs(int iteration, int iterationsWithoutChange);
extern void log_fitness(std::vector<island>& islands, int iteration);

extern json setup_iteration_json(int iteration, int score);

extern void inform_log_about_migration(
    json& iteration_json, int stuck_island, int source_island, int source_individual,
    const std::string& improvement_str, std::chrono::time_point<std::chrono::steady_clock> start,
    std::chrono::time_point<std::chrono::steady_clock> selection,
    std::chrono::time_point<std::chrono::steady_clock> selection2,
    std::chrono::time_point<std::chrono::steady_clock> end);

extern void inform_log_about_mutations(
    json& iteration_json, std::tuple<double, json>& score_result,
    std::chrono::time_point<std::chrono::steady_clock> scoring_time,
    std::chrono::time_point<std::chrono::steady_clock> mutate_time);

extern void inform_log_about_offspring_scores(json& iteration_json, std::vector<island>& islands);
extern void inform_log_about_newparent_scores(json& iteration_json, std::vector<island>& islands);

extern void finalize_iteration_log(json& iteration_json,
                                   std::chrono::time_point<std::chrono::steady_clock> end_up,
                                   std::chrono::time_point<std::chrono::steady_clock> mut,
                                   std::chrono::time_point<std::chrono::steady_clock> update);

extern void inform_log_about_mutations_json(json& iteration_json, json& mutations_json);

extern void inform_log_about_crossover_json(json& iteration_json, json& crossover_json);