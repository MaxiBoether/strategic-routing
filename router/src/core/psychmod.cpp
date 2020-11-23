
#include <cmath>
#include <iostream>
#include <unordered_map>

#include <gsl/gsl_poly.h>

#include "core/data.h"
#include "core/psychmod.h"

#ifndef ALPHA
#define ALPHA 1.0f
#endif

using namespace std;


pair<double, double> psychmod::score_routes_individually(route& p, route& q, int k) {
  double ap = 0, bp = 0, aq = 0, bq = 0, as = 0, bs = 0;
  unordered_map<link*, bool> on_q;
  for (link* l : q.links) {
    on_q[l] = true;
    aq += l->a();
    bq += l->b();
  }

  for (link* l : p.links) {
    ap += l->a();
    bp += l->b();
    if (on_q[l]) {
      as += l->a();
      bs += l->b();
    }
  }

  auto usage = score_route(ap, bp, aq, bq, as, bs, k).second;

  double time_shared = latency(as, bs, k);
  double time_o = latency(aq - as, bq - bs, k - usage) + time_shared;
  double time_a = latency(ap - as, bp - bs, usage) + time_shared;
  return {time_o, time_a};
}

pair<double, int> psychmod::score_route(route& p, route& q, int k) {
  double ap = 0, bp = 0, aq = 0, bq = 0, as = 0, bs = 0;
  unordered_map<link*, bool> on_q;
  for (link* l : q.links) {
    on_q[l] = true;
    aq += l->a();
    bq += l->b();
  }

  for (link* l : p.links) {
    ap += l->a();
    bp += l->b();
    if (on_q[l]) {
      as += l->a();
      bs += l->b();
    }
  }
  return score_route(ap, bp, aq, bq, as, bs, k);
}

pair<double, int> psychmod::score_route(double ap, double bp, double aq, double bq, double apnq,
                                    double bpnq, int k) {
  double overlapping_latency = latency(apnq, bpnq, k);

  vector<double> usages =
      calc_usage(ap, bp, aq, bq, apnq, bpnq, k);

  double score = numeric_limits<double>::max();
  int usage = 9;
  for (double u : usages) {
    double l_p = latency(ap - apnq, bp - bpnq, u);
    double l_q = latency(aq - apnq, bq - bpnq, k - u);
    int _score = u * l_p + (k - u) * l_q + k * overlapping_latency;
    if (_score < 0) //Overflow!
      _score = numeric_limits<double>::max();
    if (_score < score) {
      score = _score;
      usage = u;
    }
  }
  return {score, usage};
}

double linear_simple_example_model_2r::latency(double a, double b, double x) {
  return a * x * x + b;
}

double linear_simple_example_model_2r::a(link& l) {
  return (0.15f * l.length) / (l.freespeed * pow(l.capacity, 2));
}

double linear_simple_example_model_2r::b(link& l) { return l.length / l.freespeed; }

bool linear_simple_example_model_2r::dominating(shared_ptr<ParetoElement> par1,
                                                         shared_ptr<ParetoElement> par2) {
   return par1->b() <= par2->b() && par1->taud() <= par2->taud();
}

bool linear_simple_example_model_2r::strongly_dominating(shared_ptr<ParetoElement> par1,
                                                         shared_ptr<ParetoElement> par2) {
   return par1->b() <= par2->b() && par1->taud() <= par2->taud() && par1->shared_a() <= par2->shared_a();
}

vector<double> linear_simple_example_model_2r::calc_usage(double ap, double bp, double aq,
                                                          double bq, double apnq, double bpnq,
                                                          int k) {
  vector<double> in_bounds;
  if (ap == apnq && bp == bpnq && (ap < aq || bp < bq)) {
    in_bounds.push_back(k);
    return in_bounds;
  }
  double alpha = ALPHA;
  double overlapping_latency = latency(apnq, bpnq, k);
  double a = aq - apnq;
  double b = bq - bpnq + overlapping_latency;
  double c = ap - apnq;
  double d = bp - bpnq + overlapping_latency;
  double k2 = k * k;
  double aa = alpha * a;
  double A = -(aa * k) / c;
  double B = (d + 2 * k2 * aa) / c;
  double C = -(aa * k2 * k + alpha * k * b) / c;
  double x0, x1, x2;
  gsl_poly_solve_cubic(A, B, C, &x0, &x1, &x2);
  if (0 <= x0 && x0 <= k)
    in_bounds.push_back(x0);
  else if (x0 > k)
    in_bounds.push_back(k);
  else if (x0 < 0)
    in_bounds.push_back(0);
  if (0 <= x1 && x1 <= k)
    in_bounds.push_back(x1);
  else if (x1 > k)
    in_bounds.push_back(k);
  else if (x1 < 0)
    in_bounds.push_back(0);
  if (0 <= x2 && x2 <= k)
    in_bounds.push_back(x2);
  else if (x2 > k)
    in_bounds.push_back(k);
  else if (x2 < 0)
    in_bounds.push_back(0);
  if (in_bounds.empty()) {
    cerr << "WARNING!" << endl;
    cerr << "could not determine in-bound score for a pareto route" << endl;
    return {0.0f};
  }
  return in_bounds;
}

vector<double> user_equilibrium_2r::calc_usage(double a_p, double b_p, double a_q, double b_q,
                                               double apnq, double bpnq, int k) {
  double overlapping_latency = latency(apnq, bpnq, k);
  double aq = a_q - apnq;
  double bq = b_q - bpnq + overlapping_latency;
  double ap = a_p - apnq;
  double bp = b_p - bpnq + overlapping_latency;

  double A = ap - aq;
  double B = 2 * aq * k;
  double C = bp - bq - aq * k * k;

  double x0 = -1, x1 = -1;
  gsl_poly_solve_quadratic(A, B, C, &x0, &x1);
  vector<double> in_bounds;
  if (0 <= x0 && x0 <= k) {
    in_bounds.push_back(x0);
    return in_bounds;
  }
  if (0 <= x1 && x1 <= k) {
    in_bounds.push_back(x1);
    return in_bounds;
  }

  if (latency(ap, bp, k) <= latency(aq, bq, 0)) {
    in_bounds.push_back(k);
    return in_bounds;
  }
  in_bounds.push_back(0);
  return in_bounds;
}

  vector<double> system_optimum_2r::calc_usage(double ap, double bp, double aq, double bq,
                                               double apnq, double bpnq, int k) {
    double qa = aq - apnq;
    double qb = bq - bpnq;
    double pa = ap - apnq;
    double pb = bp - bpnq;

    double A = (3 * pa - 3 * qa);  // x^2,
    double B = (6 * qa * k);       // x,
    double C = (pb - qb - 3 * qa * k * k);
    double x0, x1;
    gsl_poly_solve_quadratic(A, B, C, &x0, &x1);
    vector<double> in_bounds;
    in_bounds.push_back(k);
    in_bounds.push_back(0);
  if (0 <= x0 && x0 <= k)
    in_bounds.push_back(x0);
  if (0 <= x1 && x1 <= k)
    in_bounds.push_back(x1);
  if (in_bounds.empty()) {
    cerr << "WARNING!" << endl;
    cerr << "could not determine in-bound score for a pareto route" << endl;
    return {0.0f};
  }
  return in_bounds;
}