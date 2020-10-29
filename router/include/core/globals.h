#pragma once

#include <unordered_map>
#include <vector>

#include "core/data.h"
#include "core/psychmod.h"

// variables
extern std::vector<person> persons;
extern std::vector<node*> nodes;
extern std::vector<std::vector<link*>> adj;

extern int number_agents;

extern PSYCH_MODEL_CLASS psychological_model;