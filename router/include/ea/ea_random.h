#pragma once
#include <random>

#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>

extern auto getGenerator() -> std::mt19937&;
auto getGSLRng() -> gsl_rng*;

/* This file is a little fix to ensure compatability of random.h with Leon's SSOTD EA. In a
 * nutshell, we implemented stuff in random.h instead of just doing declarations there. Now, after
 * the S(0+n)OTD EA is split up in several files, they cannot all import random.h without creating
 * multiple definitions in their object files, which is why the linker cannot link the final binary.
 * The solution to this is importing random.h in ea.cpp to craete the real bindings there and use
 * there declarions here everywhere else. K thx goodbye. */