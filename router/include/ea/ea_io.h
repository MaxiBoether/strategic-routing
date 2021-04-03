#pragma once
#include <string>
#include <vector>

#include "ea/ea_data.h"

extern bool interrupted;
extern void setup_signal_handling();
extern void parse_env();
extern std::string get_intermediate_filename(int iteration, int score);
extern void inform_about_settings(int k, int origin, int destination, int score);
extern void print_iteration_infos(int iteration, std::vector<island>& islands);