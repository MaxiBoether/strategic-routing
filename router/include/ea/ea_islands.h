#pragma once
#include <utility>
#include <vector>

#include "ea/ea_data.h"

#include "json.hpp"

using json = nlohmann::json;

extern void initialize_islands(std::vector<island>& islands, int k, int origin, int destination);

extern void score_islands(std::vector<island>& islands, int k);
extern void sort_islands(std::vector<island>& islands);

extern void do_migrations(std::vector<island>& islands, int iteration, int& iterationsWithoutChange,
                          json& iteration_json);

extern auto update_islands(std::vector<island>& islands, int& iterationsWithoutChange)
    -> std::pair<individual*, int>;