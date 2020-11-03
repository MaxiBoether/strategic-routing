#pragma once
#include <vector>

#include "core/data.h"

using namespace std;

class psychmod {
 public:
  virtual double b(link& l) = 0;
  virtual double a(link& l) = 0;
  virtual double latency(double a, double b, double x) = 0;
  virtual vector<double> calc_usage(double ap, double bp, double aq, double bq, double apnq,
                                    double bpnq, int k) = 0; //apnq and bpnq refer tp the parameters a and b of the shared edges of p and q while the others refer to the whole paths p and q
  virtual bool dominating(shared_ptr<ParetoElement> par1, shared_ptr<ParetoElement> par2) = 0;
  virtual bool strongly_dominating(shared_ptr<ParetoElement> par1, shared_ptr<ParetoElement> par2) = 0;
  pair<double, int> score_route(route& p, route& q, int k);
  pair<double, int> score_route(double ap, double bp, double aq, double bq, double apnq,
                                    double bpnq, int k);
                                    
  pair<double, double> score_routes_individually(route& p, route& q, int k);
};

class linear_simple_example_model_2r : public psychmod {
 public:
  virtual double b(link& l);
  virtual double a(link& l);
  virtual double latency(double a, double b, double x);
  virtual vector<double> calc_usage(double ap, double bp, double aq, double bq, double apnq,
                                    double bpnq, int k);
  virtual bool dominating(shared_ptr<ParetoElement> par1, shared_ptr<ParetoElement> par2);
  virtual bool strongly_dominating(shared_ptr<ParetoElement> par1, shared_ptr<ParetoElement> par2);
};

class user_equilibrium_2r : public linear_simple_example_model_2r {
 public:
  virtual vector<double> calc_usage(double ap, double bp, double aq, double bq, double apnq,
                                    double bpnq, int k);
};

class system_optimum_2r : public linear_simple_example_model_2r {
 public:
  virtual vector<double> calc_usage(double ap, double bp, double aq, double bq, double apnq,
                                    double bpnq, int k);
};
