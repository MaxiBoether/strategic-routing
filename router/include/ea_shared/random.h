#pragma once
#include <sys/time.h>

#include <algorithm>
#include <random>
#include <vector>

#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>

// thanks
// to https://codereview.stackexchange.com/questions/109260/seed-stdmt19937-from-stdrandom-device
auto ProperlySeededRandomEngine() -> std::mt19937 {
  std::vector<uint32_t> random_data(std::mt19937::state_size);
  std::random_device source;
  std::generate(random_data.begin(), random_data.end(), std::ref(source));
  std::seed_seq seeds(random_data.begin(), random_data.end());
  std::mt19937 engine(seeds);
  return engine;
}

auto getGenerator() -> std::mt19937& {
  static std::mt19937 gen =
      ProperlySeededRandomEngine();  // obtain a random engine once, reuse whenever we want to. this
                                     // avoid having to re-seed every time
  return gen;
}

auto createGSL() -> gsl_rng* {
  const gsl_rng_type* T;
  gsl_rng* rand;
  gsl_rng_env_setup();
  struct timeval tv{};  // Seed generation based on time
  gettimeofday(&tv, nullptr);
  unsigned long mySeed = tv.tv_sec + tv.tv_usec;
  T = gsl_rng_taus2;
  rand = gsl_rng_alloc(T);
  gsl_rng_set(rand, mySeed);
  std::cout << "seed: " << mySeed << std::endl;
  return rand;
}

auto getGSLRng() -> gsl_rng* {
  static gsl_rng* rand = createGSL();
  return rand;
}

void clearGSLRandom() {
  gsl_rng_free(getGSLRng());
}
