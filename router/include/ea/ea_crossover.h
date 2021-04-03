#pragma once
#include "ea/ea_data.h"

#include "json.hpp"

using json = nlohmann::json;

individual makeDiverseBaby(individual& parent1, individual& parent2, int k);
individual makeDiverseBaby2(individual& parent1, individual& parent2, int k);
individual makeHeuristicBaby(individual& parent1, individual& parent2, int k);
individual makeHeuristicBaby2(individual& parent1, individual& parent2, int k);

void handle_innerisland_crossover(island& _island, json& iteration_json, int k);
void handle_crossisland_crossover(std::vector<island>& islands, json& iteration_json, int k);