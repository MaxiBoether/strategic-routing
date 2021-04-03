#pragma once
#include <vector>

#include "ea/ea_data.h"

#include "json.hpp"

using json = nlohmann::json;

extern void mutate(individual& indiv, int k, int origin, int destination, int indiv_id,
                   json& iteration_json, int iteration, int iterationsWithoutChange,
                   std::vector<int>& lastExchanges, int island_id);