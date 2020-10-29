#pragma once
#include <any>
#include <list>
#include <memory>
#include <unordered_map>

#include "core/data.h"

using namespace std;
void fill_best_pars_dijkstra(int to, unordered_map<int, bool> inactive=unordered_map<int, bool>());

void fill_best_pars_dijkstra_forward(int from, unordered_map<int, bool> inactive=unordered_map<int, bool>());

extern vector<double> bestAs, bestBs, bestAsForward, bestBsForward; // For dijkstra & airline local opt

extern vector<double> origTt, origPartA, origPartB;  // original route prefix sums

extern unordered_map<int, int> nodes_original_route;  // maps a node id to its index in orig route

shared_ptr<route> dijkstra(int a, int b, shared_ptr<route> original_route = nullptr);

// pair<double, int> route_score_disjoint(shared_ptr<route> original_route,
//                                           double a, double b, int k);

// pair<double, int> route_score_disjoint(shared_ptr<route> original_route,
//                                           shared_ptr<ParetoElement>& pareto_route, int k);

// void pareto_dijkstra_local_opt(int a, int from, int to, vector<list<shared_ptr<ParetoElement>>>& pareto,
//                                shared_ptr<route> original_route, int k, double qot,
//                                unordered_map<int, bool> inactive,
//                                pair<double, double> (*lower_bound_score)(shared_ptr<ParetoElement>, int, int, int, int,
//                                                            shared_ptr<route>, int));

bool standard_prio(pair<shared_ptr<ParetoElement>, int> left, pair<shared_ptr<ParetoElement>, int> right);
 
bool astar_prio_dijkstra(pair<shared_ptr<ParetoElement>, int> left, pair<shared_ptr<ParetoElement>, int> right);

bool astar_prio_airline(pair<shared_ptr<ParetoElement>, int> left, pair<shared_ptr<ParetoElement>, int> right);

pair<double, long long> pareto_dijkstra_local_opt(int a, int from, int to, vector<vector<shared_ptr<ParetoElement>>>& pareto,
                               shared_ptr<route> original_route, int k, double qot,
                               unordered_map<int, bool> inactive = unordered_map<int, bool>(),
                               pair<double, double> (*lower_bound_score)(shared_ptr<ParetoElement>, int, int, int, int,
                                                           shared_ptr<route>, int) = nullptr,
       bool (*prio)(pair<shared_ptr<ParetoElement>, int> left, pair<shared_ptr<ParetoElement>, int> right)=&standard_prio);

// void pareto_dijsktra(int a, vector<list<shared_ptr<ParetoElement>>>& pareto, unordered_map<int, bool> inactive);

long long pareto_dijsktra(int a, int b, vector<vector<shared_ptr<ParetoElement>>>& pareto, shared_ptr<route> original_route, unordered_map<int, bool> inactive=unordered_map<int, bool>(),
       bool (*prio)(pair<shared_ptr<ParetoElement>, int> left, pair<shared_ptr<ParetoElement>, int> right)=&standard_prio);

long long pareto_dijsktra_4d(int a, int b, vector<vector<shared_ptr<ParetoElement>>>& pareto, shared_ptr<route> original_route, unordered_map<int, bool> inactive=unordered_map<int, bool>(),
       bool (*prio)(pair<shared_ptr<ParetoElement>, int> left, pair<shared_ptr<ParetoElement>, int> right)=&standard_prio);

void pareto_dijsktra_4d_1D(int a, int b, vector<vector<shared_ptr<ParetoElement>>>& pareto, shared_ptr<route> original_route, unordered_map<int, bool> inactive=unordered_map<int, bool>(),
       bool (*prio)(pair<shared_ptr<ParetoElement>, int> left, pair<shared_ptr<ParetoElement>, int> right)=&standard_prio);

pair<double, long long> pareto_dijkstra_local_opt_4d(int a, int from, int to, vector<vector<shared_ptr<ParetoElement>>>& pareto,
                               shared_ptr<route> original_route, int k, double qot, unordered_map<int, bool> is_orig_edge,
                               pair<double, double> (*lower_bound_score)(shared_ptr<ParetoElement>, int, int, int, int,
                                                           shared_ptr<route>, int),
                     bool (*prio)(pair<shared_ptr<ParetoElement>, int>, pair<shared_ptr<ParetoElement>, int>));

pair<double, long long> pareto_dijkstra_local_opt_4d_1D(int a, int from, int to, vector<vector<shared_ptr<ParetoElement>>>& pareto,
              shared_ptr<route> original_route, int k, double qot, unordered_map<int, bool> is_orig_edge,
              pair<double, double> (*lower_bound_score)(shared_ptr<ParetoElement>, int, int, int, int,
                                          shared_ptr<route>, int),
bool (*prio)(pair<shared_ptr<ParetoElement>, int>, pair<shared_ptr<ParetoElement>, int>));


shared_ptr<vector<double>> dijkstra_for_opt(int v, bool doA, unordered_map<int, bool> inactive, vector<vector<link*>> adj, bool forward);

bool insert_and_dominate(list<shared_ptr<ParetoElement>>& A, shared_ptr<ParetoElement>& frag);

double airline_dist(int from, int to);

link* artificial_link(int from, int to, double freespeed, double capacity);

double score_for_relax(int idc, int idv, shared_ptr<ParetoElement> par, double k);

int index_in_original(int v);

int is_orig_node(int node, shared_ptr<route> orig);

void prepare_original_route(shared_ptr<route> original_route, unordered_map<int, bool>& inactive);

pair<double, double> get_lazy_airval(int node, int to);

pair<double, double> get_lazy_airval_forward(int node, int from);

void check_route_sanity(route& r, string routeName);