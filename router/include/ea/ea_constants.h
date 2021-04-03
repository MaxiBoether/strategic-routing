#pragma once
static const double stddev_factor_randdijk = 0.8;
static const double min_weight_randdijk = 1;

static const double average_factor_mutate = 0.25;
static const double stddev_factor_mutate = 0.5;
static const double average_mutations = 1.5;

static const double stddev_mutation_capacity_start = 10;
static const double stddev_mutation_capacity_advance = 10;

static const double prob_mutate_randpair = 10;
static const double prob_mutate_weightpair = 30;
static const double prob_mutate_newroute = 30;
static const double prob_mutate_capacity = 50;
static const double prob_mutate_exchange = 15;

static constexpr int gap_between_exchanges = 4;
static constexpr int max_prob_mutate_exchange = 40;
static constexpr float max_prob_mutate_reached_after =
    0.2F;  // 20% of iterations_needed_for_convergeance without change for max prob

static const double iterations_after_newrout_prob_one = 500;
static const int iterations_constant_newroute_prob = 10;

static constexpr int minimum_queue_size = 65565;
static const double queue_size_factor = 0.3;

static constexpr int iterations_RSMEBelowTen = 5;
static constexpr int iterations_RSMENotChanging = 10;
static constexpr int iterations_hardlimit = 10000;

static constexpr int offspring_per_wife = 2;
