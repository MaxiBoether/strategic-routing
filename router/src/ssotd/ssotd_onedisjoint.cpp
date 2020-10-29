// vim: et ts=4 sw=4

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <functional>
#include <iostream>
#include <iterator>
#include <limits>
#include <map>
#include <memory>
#include <numeric>
#include <set>
#include <sstream>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include "core/data.h"
#include "core/globals.h"
#include "ssotd/ssotd_core.h"
using namespace std;
using ll = long long;

#define MAX_ORIGINAL_ROUTE_NODES 1024

pair<double, double> lower_bound_score_airline(shared_ptr<ParetoElement> par, int from, int to,
                                               int c, int v, shared_ptr<route> original_route,
                                               int k) {
  (void) from;
  int idc = index_in_original(c);
  auto [a, b] = get_lazy_airval(v, to);
  auto score = psychological_model.score_route(par->a() + a + origPartA[idc], par->b() + b + origPartB[idc],
                                               original_route->a(), original_route->b(), origPartA[idc], origPartB[idc], k);
  if (score.second > 0)
    return make_pair(score.first, score_for_relax(idc, index_in_original(v), par, k));
  return make_pair(HUGE_VAL, -1);
}

pair<double, double> lower_bound_score_zero(shared_ptr<ParetoElement> par, int from, int to, int c,
                                            int v, shared_ptr<route> original_route, int k) {
  (void) from;  (void) to;
  int idc = index_in_original(c);
  double newA = par->a() + origPartA.at(idc);
  double newB = par->b() + origPartB.at(idc);
  auto score = psychological_model.score_route(newA, newB, original_route->a(),
                                               original_route->b(), origPartA.at(idc), origPartB.at(idc), k);
  if (score.second > 0)
    return make_pair(score.first, score_for_relax(idc, index_in_original(v), par, k));
  return make_pair(HUGE_VAL, -1);
}

pair<double, double> lower_bound_score_dijkstra(shared_ptr<ParetoElement> par, int from, int to,
                                                int c, int v, shared_ptr<route> original_route,
                                                int k) {
  (void) from; (void) to;
  int idc = index_in_original(c);
  double newA = par->a() + origPartA.at(idc) + bestAs[v];
  double newB = par->b() + origPartB.at(idc) + bestBs[v];
  auto score = psychological_model.score_route(newA, newB, original_route->a(),
                                               original_route->b(), origPartA.at(idc), origPartB.at(idc), k);
  if (score.second > 0)
    return make_pair(score.first, score_for_relax(idc, index_in_original(v), par, k));
  return make_pair(HUGE_VAL, -1);
}

void sanity_check_1D(unordered_map<int, bool> original_edges, shared_ptr<route> alternative, shared_ptr<route> original, double score, int usage, int k) {

  int crosses = 0;
  int splits = 0;
  bool splitted = false;
  for (size_t i=0; i<alternative->links.size(); i++) {
    auto l = alternative->links[i];
    if (original_edges[l->id]) {
      splitted = false;
    } else {
      if (!splitted)
        splits++;
      splitted = true;

      if (i>0 && !original_edges[alternative->links[i-1]->id] && is_orig_node(l->from, original))
        crosses++;
    }
  }
  cout << "splits: " << splits << endl;
  
  cout << "Crosses: " << crosses << endl;

  double a = 0, b = 0, sa = 0, sb = 0;
  for (auto l : alternative->links) {
    a += l->a();
    b += l->b();
    if (original_edges[l->id]) {
      sa += l->a();
      sb += l->b();
    }
  }

  double othera = alternative->a();
  double otherb = alternative->b();
  if (a != othera || b != otherb)
    cout << "Holy shit! Route parameters gone wrong: " << a << " " << b << " != " << othera << " " << otherb << endl;
  auto [actual_score, actual_usage] = psychological_model.score_route(a,b, original->a(), original->b(), sa, sb, k);
  if (actual_score != score || actual_usage != usage)
    cout << "Holy shit! Scoring gone wrong: " << actual_score << " (" << actual_usage << ") != " << score << " (" << usage << ")" << endl;
}

pair<shared_ptr<route>, double> ssotd_route(int a, int b, shared_ptr<route> original_route, int k,
                                            string optimization) {
  map<pair<int, int>, vector<shared_ptr<ParetoElement>>> paretoFronts;
  unordered_map<int, bool> inactive;
  prepare_original_route(original_route, inactive);
  double qot = k * psychological_model.latency(original_route->a(), original_route->b(), k);
  cout << "DIJKSTRA OT: " << qot << std::endl;
  cout << "Calculating pareto fronts." << endl;
  function<pair<double, int>(int, vector<vector<shared_ptr<ParetoElement>>>*)> pareto_dijk;

  if ("none" == optimization || "hop" == optimization) {
    cout << "Doing no optimization" << endl;
    pareto_dijk = [&inactive, &original_route, qot, b](int c, vector<vector<shared_ptr<ParetoElement>>>* pareto) {
      pareto_dijsktra(c, b, *pareto, original_route, inactive);
      return make_pair(qot,-1);
    };
  } else if ("simple_local_opt" == optimization) {
    cout << "Doing simple optimization" << endl;
    pareto_dijk = [a, b, &original_route, k, &qot, &inactive](
                      int c, vector<vector<shared_ptr<ParetoElement>>>* pareto) {
      return pareto_dijkstra_local_opt(c, a, b, *pareto, original_route, k, qot, inactive, &lower_bound_score_zero);
    };
  } else if ("dijkstra_local_opt" == optimization) {
    cout << "Doing dijkstra optimization" << endl;
    auto start = chrono::steady_clock::now();
    fill_best_pars_dijkstra(b);
    auto end = chrono::steady_clock::now();
    cout << "Route specific precalculation time: "
         << chrono::duration_cast<chrono::microseconds>(end - start).count() << endl;
    pareto_dijk = [a, b, &original_route, k, &qot, &inactive](
                      int c, vector<vector<shared_ptr<ParetoElement>>>* pareto) {
      return pareto_dijkstra_local_opt(c, a, b, *pareto, original_route, k, qot, inactive, &lower_bound_score_dijkstra);
    };
  } else if ("airline_local_opt" == optimization) {
    cout << "Doing airline optimization" << endl;
    
    bestAs = vector<double>(nodes.size(),0);
    bestBs = vector<double>(nodes.size(),0);
    cerr << "Assuming best freeflowspeed = 130 and best capacity = 2288 (sane-berlin.xml)!!!!"
         << endl;
    bestAs = vector<double>(nodes.size());
    bestBs = vector<double>(nodes.size());
    pareto_dijk = [a, b, &original_route, k, &qot, &inactive](
                      int c, vector<vector<shared_ptr<ParetoElement>>>* pareto) {
      return pareto_dijkstra_local_opt(c, a, b, *pareto, original_route, k, qot, inactive, &lower_bound_score_airline);
    };
  } else if ("airline_astar_opt" == optimization) {
    cout << "Doing airline astar optimization" << endl;
    cerr << "Assuming best freeflowspeed = 130 and best capacity = 2288 (sane-berlin.xml)!!!!"
         << endl;
    bestAs = vector<double>(nodes.size());
    bestBs = vector<double>(nodes.size());
    pareto_dijk = [a, b, &original_route, k, &qot, &inactive](
                      int c, vector<vector<shared_ptr<ParetoElement>>>* pareto) {
      return pareto_dijkstra_local_opt(c, a, b, *pareto, original_route, k, qot, inactive, &lower_bound_score_airline, &astar_prio_airline);
    };
  } else if ("dijkstra_astar_opt" == optimization) {
    cout << "Doing dijkstra astar optimization" << endl;
    bestAs = vector<double>(nodes.size(),0);
    bestBs = vector<double>(nodes.size(),0);
    auto start = chrono::steady_clock::now();
    fill_best_pars_dijkstra(b);
    auto end = chrono::steady_clock::now();
    cout << "Route specific precalculation time: "
         << chrono::duration_cast<chrono::microseconds>(end - start).count() << endl;
    pareto_dijk = [a, b, &original_route, k, &qot, &inactive](
                      int c, vector<vector<shared_ptr<ParetoElement>>>* pareto) {
      return pareto_dijkstra_local_opt(c, a, b, *pareto, original_route, k, qot, inactive, &lower_bound_score_dijkstra, &astar_prio_dijkstra);
    };
  } else
    cerr << "Invalid parameter for optimization given: " << optimization << endl;


  auto start = chrono::steady_clock::now();
  double upperBound = qot;
  ll visits = 0;
#pragma omp parallel for default(none) shared(paretoFronts, original_route, adj, inactive, upperBound, k, \
                                              b, pareto_dijk, visits) schedule(dynamic, 2) num_threads(8)

  for (unsigned int lid = 0; lid < original_route->links.size(); lid++) {
    // iterate over all vertices of the original route except the last
    int v = original_route->links[lid]->from;
    vector<vector<shared_ptr<ParetoElement>>> pareto(adj.size() + 1);
    // vector<list<shared_ptr<ParetoElement>>> pareto(adj.size() +1);
    auto [newUpperBound, new_visits] = pareto_dijk(v, &pareto);
    visits += new_visits;
    if (newUpperBound < upperBound)
      upperBound = newUpperBound;

#pragma omp critical
    for (unsigned int _lid = lid + 1; _lid < original_route->links.size(); _lid++) {
      paretoFronts[{lid, _lid}] = pareto[original_route->links[_lid]->from];
    }
    paretoFronts[{lid, original_route->links.size()}] = pareto[original_route->links.back()->to];
  }
  auto end = chrono::steady_clock::now();
  cout << "Pareto-dijkstra time: " << chrono::duration_cast<chrono::microseconds>(end - start).count() << endl;
  cout << "Node visits: " << visits << endl;
  cout << "start evaluation" << endl;
 

  start = chrono::steady_clock::now();

  double best_ot = numeric_limits<double>::max();
  double best_usage = 0.0;
  shared_ptr<ParetoElement> best;
  int bestI = 0, bestJ = 0;
  double shared_a = 0, shared_b = 0;
  int n = (original_route->links.size()+1);
  auto pareto_sizes = vector<int>(n*n/2-n);  //n^2/2 - n Pareto-Fronts

  for (unsigned int i = 0; i < original_route->links.size(); i++) {
    for (unsigned int j = i+1; j < original_route->links.size() + 1; j++) {
      pareto_sizes.push_back(paretoFronts[{i,j}].size());
      for (auto par : paretoFronts[{i,j}]) {
        shared_a = origPartA[i] + origPartA.back() - origPartA[j];
        shared_b = origPartB[i] + origPartB.back() - origPartB[j];
        auto [ot, usage] = psychological_model.score_route(par->a() + shared_a, par->b() + shared_b, origPartA.back(),
                                      origPartB.back(), shared_a, shared_b, k);
         if (ot < best_ot) {
          best_ot = ot;
          best_usage = usage;
          best = par;
          bestI = i;
          bestJ = j;
        }
      }
    }
  }
  end = chrono::steady_clock::now();
  cout << "Evaluation time: "
        << chrono::duration_cast<chrono::microseconds>(end - start).count() << endl;

  shared_a = origPartA[bestI] + origPartA.back() - origPartA[bestJ];
  shared_b = origPartB[bestI] + origPartB.back() - origPartB[bestJ];
  // cout << "Latency pareto-route: " << psychological_model.latency(shared_a, shared_b, k) + psychological_model.latency(best->a(), best->b(), best_usage) << endl;
  // cout << "Latency original-route: " << psychological_model.latency(shared_a, shared_b, k) + psychological_model.latency(origPartA[bestJ] - origPartA[bestI], origPartB[bestJ] - origPartB[bestI], k- best_usage) << endl;
  // cout << "Latency if all on original-route: " << psychological_model.latency(origPartA.back(), origPartB.back(),k) << endl;
  
  auto total_candidates = accumulate(pareto_sizes.begin(), pareto_sizes.end(), 0);
  nth_element(pareto_sizes.begin(), pareto_sizes.begin() + pareto_sizes.size() / 2, pareto_sizes.end());
  auto mean_pareto_set_size = pareto_sizes[pareto_sizes.size()/2];
  cout << "Found " << total_candidates << " pareto-optimal routes" << endl;
  cout << "Mean Pareto-set size: " << mean_pareto_set_size << endl;
  cout << "Sum Pareto-set size: " << total_candidates << endl;
  if (total_candidates > 0) {
  cout << "\nBEST PARETO OT: " << best_ot << endl;
  cout << "BEST PARETO SHARES " << bestI + original_route->links.size() - bestJ << " of " << original_route->links.size() << " edges of the original route (leaving at " << bestI 
       << " and reuniting at " << bestJ << ")" <<endl;
  }
  if (best_ot > qot)
    return {original_route, 0.0};

  vector<link*> parLinks = *(best->collectLinks());
  vector<link*> routeLinks;
  routeLinks.reserve(original_route->links.size() - bestJ + bestI + parLinks.size());
  copy(original_route->links.begin(), original_route->links.begin() + bestI, back_inserter(routeLinks));
  copy(parLinks.begin(), parLinks.end(), back_inserter(routeLinks));
  copy(original_route->links.begin() + bestJ, original_route->links.end(), back_inserter(routeLinks));
  auto res = make_shared<route>(routeLinks);
  // sanity_check_1D(inactive, res, original_route, best_ot, best_usage, k);
  // check_route_sanity(*res, "alternative route");
  cout << "Collected SSOTD route" << endl;

  return {res, best_usage};
}

void ssotd(int source, int destination, vector<int> pids, string optimization) {
  shared_ptr<route> original_route = dijkstra(source, destination);
    cout << "Length original: " << original_route->links.size() << endl;
    cout << "K: " << pids.size() << endl;
  if (original_route->links.size() > MAX_ORIGINAL_ROUTE_NODES) {
    cerr << "original route is too big (" << original_route->links.size() << " nodes )" << endl;
    exit(1);
  }
  auto start = chrono::steady_clock::now();
  pair<shared_ptr<route>, double> ssotd_res =
      ssotd_route(source, destination, original_route, pids.size(), optimization);
  auto end = chrono::steady_clock::now();
  cout << "time used: " << chrono::duration_cast<chrono::microseconds>(end - start).count() << endl;
  double usage = ssotd_res.second / static_cast<double>(pids.size());
  cout << "normalized usage of the pareto route: " << usage << endl;
  // check_route_sanity(*(ssotd_res.first), "alternative route");
  // check_route_sanity(*original_route, "original route");
  for (int pid : pids) {
    if ((rand() % (1 << 16)) / static_cast<double>(1 << 16) < usage)
      persons[pid].r = ssotd_res.first;
    else
      persons[pid].r = original_route;
  }
  cout << "SSOTD assignment completed." << endl;
}

void do_routing(int argc, char* argv[]) {
  string optimization;
  if (argc > 0) {
    optimization = argv[0];
    int pos1 = optimization.find_first_not_of("\t\n\v\f\r ");
    int pos2 = optimization.find_last_not_of("\t\n\v\f\r ");
    optimization = optimization.substr(pos1, pos2 - pos1 + 1);
  } else
    optimization = "none";

  map<pair<pair<int, int>, string>, vector<int>> c;
  for (unsigned int pid = 0; pid < persons.size(); pid++) {
    auto& p = persons[pid];
    string s = p.timestr;  // maybe
    c[{{p.origin, p.destination}, s}].push_back(pid);
  }
  // do ssotd for all
  for (auto& [sdts, pv] : c) {
    number_agents = pv.size();
    ssotd(sdts.first.first, sdts.first.second, pv, optimization);
  }
  cout << "entire SSOTD routing complete" << endl;
}
