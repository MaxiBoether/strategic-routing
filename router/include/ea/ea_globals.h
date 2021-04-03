#pragma once
#include "ea/ea_scoring.h"

#ifndef SCORING_METHOD
#define SCORING_METHOD FrankWolfeScorer  // NOLINT(cppcoreguidelines-macro-usage)
#endif

extern Scorer* scorer;
extern int inds;
extern int n;
extern int islandsCount;
extern int iterations_needed_for_convergeance;
extern int iterations_for_migration;
extern int crossisland_iterations;
extern int queue_size;
extern individual (*crossover_func)(individual&, individual&, int);
extern std::string crossover_func_str;
extern std::string intermediate_results_prefix;
extern std::string intermediate_folder;