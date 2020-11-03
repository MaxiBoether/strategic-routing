// vim: et ts=4 sw=4

#include <algorithm>
#include <bitset>
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
#include <queue>
#include <random>
#include <set>
#include <sstream>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include "argh.h"
#include "core/data.h"
#include "core/globals.h"
#include "core/routing.h"
#include "ssotd/ssotd_core.h"
using namespace std;
using ll = long long;
template <class T>
using minq = priority_queue<T, vector<T>, greater<>>;

//This file refers to the SAP algorithm


pair<double, double> lower_bound_score_dijkstra(shared_ptr<ParetoElement> par, int from, int to,
                                                int c, int v, shared_ptr<route> original_route,
                                                int k) {
  (void) from;  (void) to;
  double newA = par->a() + bestAs[v];
  double newB = par->b() + bestBs[v];
  auto score = psychological_model.score_route(newA, newB, original_route->a(), original_route->b(),
                                               par->shared_a(), par->shared_b(), k);
  if (score.second > 0)
    return make_pair(score.first, score_for_relax(index_in_original(c), index_in_original(v), par, k));
  return make_pair(HUGE_VAL, -1);
}


pair<shared_ptr<route>, double> ssotd_route(int a, int b, shared_ptr<route> original_route, int k, string optimization) {
  vector<vector<shared_ptr<ParetoElement>>> paretoFront(adj.size());
  unordered_map<int, bool> is_orig_edge;
  prepare_original_route(original_route, is_orig_edge);
  double qot = k * psychological_model.latency(original_route->a(), original_route->b(), k);

  cout << "DIJKSTRA OT: " << qot << endl;
  cout << "Calculating pareto fronts." << endl;
  double upperBound = qot;
  auto start = chrono::steady_clock::now();
  auto end = start;
  ll visits;

 
  cout << "Doing dijkstra-astar optimization" << endl;
  start = chrono::steady_clock::now();
  fill_best_pars_dijkstra(b);
  end = chrono::steady_clock::now();
  cout << "Route specific precalculation time: "
        << chrono::duration_cast<chrono::microseconds>(end - start).count() << endl;
  start = chrono::steady_clock::now();
  visits = pareto_dijkstra_local_opt_4d(a, a, b, paretoFront, original_route, k, upperBound, is_orig_edge, &lower_bound_score_dijkstra, &astar_prio_dijkstra).second;
    
   end = chrono::steady_clock::now();
   cout << "Pareto-dijkstra time: " << chrono::duration_cast<chrono::microseconds>(end - start).count() << endl;
   cout << "Node visits: " << visits << endl;

   start = chrono::steady_clock::now();

    double best_ot = numeric_limits<double>::max();
    double best_usage = 0.0;
    shared_ptr<ParetoElement> best;

    if (paretoFront[b].size() > 0) {

      for (shared_ptr<ParetoElement>& par : paretoFront[b]) {
        auto [ot, usage] =
          psychological_model.score_route(par->a(), par->b(), original_route->a(),
                                          original_route->b(), par->shared_a(), par->shared_b(), k);
        if (ot < best_ot) {
          best_ot = ot;
          best_usage = usage;
          best = par;
        }
      }

        end = chrono::steady_clock::now();
        cout << "Evaluation time: " << chrono::duration_cast<chrono::microseconds>(end - start).count() << endl;
        cout << "Mean Pareto-set size: " << paretoFront[b].size() << endl;
        cout << "Sum Pareto-set size: " << paretoFront[b].size() << endl;
        cout << "Found " << paretoFront[b].size() << " pareto-optimal routes" << endl; 
        cout << "\nBEST PARETO OT: " << best_ot << endl;
        cout << "\nSELECTED ALTERNATIVE: a=" << best->a() << " b=" << best->b() << " sa=" << best->shared_a() << " sb=" << best->shared_b() << endl;
        cout << "b/a=" << best->b() / best->a() << endl;
      } else {
        end = chrono::steady_clock::now();
        cout << "Evaluation time: " << chrono::duration_cast<chrono::microseconds>(end - start).count() << endl;
        cout << "Mean Pareto-set size: 0" << endl;
        cout << "Sum Pareto-set size: 0" << endl;
        cout << "Found 0 pareto-optimal routes" << endl; 
        cout << "\nBEST PARETO OT: " << qot << endl;
        return {original_route, 0.0};
      }

  if (best_ot > qot)
    return {original_route, 0.0};

  auto res = shared_ptr<route>(best->collectRoute());
  return {res, best_usage};
  cout << "Collected SSOTD route" << endl;

  return {res, best_usage};
}

void ssotd(int source, int destination, vector<int> pids, string optimization) {
  shared_ptr<route> original_route = dijkstra(source, destination);
    cout << "Length original: " << original_route->links.size() << endl;
    cout << "K: " << pids.size() << endl;
  auto start = chrono::steady_clock::now();
  pair<shared_ptr<route>, double> ssotd_res =
      ssotd_route(source, destination, original_route, pids.size(), optimization);
  auto end = chrono::steady_clock::now();
  cout << "time used: " << chrono::duration_cast<chrono::microseconds>(end - start).count() << endl;
  double usage = ssotd_res.second / static_cast<double>(pids.size());
  cout << "normalized usage of the pareto route: " << usage << endl;
 
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

