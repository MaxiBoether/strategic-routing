#include <algorithm>
#include <any>
#include <chrono>
#include <cmath>
#include <iostream>
#include <list>
#include <map>
#include <string>

#include <cstdlib>
#include <fstream>
#include <iostream>
#include <iterator>
#include <limits>
#include <memory>
#include <tuple>
#include <utility>
#include <vector>

#include "core/data.h"
#include "core/globals.h"
#include "core/routing.h"
#include "ssotd/ssotd_core.h"

using namespace std;

//This file refers to the D-SAP algorithm

pair<double, double> lower_bound_score_dijkstra(shared_ptr<ParetoElement> par, int from, int to, int c, int v,
                                  shared_ptr<route> original_route, int k) {
  (void) from; (void) c;
  auto score = psychological_model.score_route(par->a() + bestAs[v], par->b() + bestBs[v], original_route->a(), original_route->b(), 0 , 0, k);
  if (score.second > 0)
    return make_pair(score.first, to == v ? score.first : -1);
  return make_pair(HUGE_VAL, -1);
}

pair<shared_ptr<route>, double> ssotd_route(int a, int b, shared_ptr<route> original_route, int k,
                                            string optimization) {
  (void)optimization;
  vector<vector<shared_ptr<ParetoElement>>> pareto(adj.size());
  unordered_map<int, bool> inactive;
  for_each(original_route->links.begin(), original_route->links.end(),
           [&inactive](link* l) { inactive[l->id] = true; });
  pareto.resize(adj.size());

  double qot = k * psychological_model.latency(original_route->a(), original_route->b(), k);
  std::cout << "DIJKSTRA OT: " << qot << std::endl;
  cout << "Doing dijkstra-astar optimization" << endl;
  auto start = chrono::steady_clock::now();
  fill_best_pars_dijkstra(b);
  auto end = chrono::steady_clock::now();
  cout << "Route specific precalculation time: "
        << chrono::duration_cast<chrono::microseconds>(end - start).count() << endl;
  start = chrono::steady_clock::now();
  auto [bound, visits] = pareto_dijkstra_local_opt(a, a, b, pareto, original_route, k, qot, inactive, &lower_bound_score_dijkstra, &astar_prio_dijkstra);
  (void) bound;
  end = chrono::steady_clock::now();
  cout << "Node visits: " << visits << endl;
  cout << "Pareto-dijkstra time: " << chrono::duration_cast<chrono::microseconds>(end - start).count() << endl;
  
  if (pareto[b].empty()) {
    cout << "Found no useful pareto-routes." << endl;
    return {original_route, 0.0};
  }
  int mycounter = 0;
  for (size_t i=0; i < nodes.size(); i++) {
    if (pareto[i].size() > 0)
      mycounter++;
  }
  cout << "Visited " << mycounter << " nodes" << endl;
  cout << "Found " << pareto[b].size() << " pareto-optimal routes" << endl;
  cout << "Mean Pareto-set size: " << pareto[b].size() << endl;
  cout << "Sum Pareto-set size: " << pareto[b].size() << endl;
    

  auto start = chrono::steady_clock::now();
  
  pair<double, int> score, best_score = {HUGE_VAL, 0};
  auto best_elem = pareto[b].begin();
  for (auto current_elem = pareto[b].begin(); current_elem != pareto[b].end(); current_elem++) {
    score = psychological_model.score_route((*current_elem)->a(), (*current_elem)->b(), original_route->a(), original_route->b(), 0 , 0, k);
    if (best_score.first > score.first) {
      best_score = score;
      best_elem = current_elem;
    }
  }
    auto end = chrono::steady_clock::now();
    cout << "Evaluation time: "
         << chrono::duration_cast<chrono::microseconds>(end - start).count() << endl;

  cout << "BEST PARETO OT: " << best_score.first << endl;
  if (best_score.first > qot)
    return {original_route, 0.0};

  cout << "a: " << (*best_elem)->a() << "  b: " << (*best_elem)->b() <<endl;
  auto res = shared_ptr<route>((*best_elem)->collectRoute());
  return {res, best_score.second};
}

void ssotd(int source, int destination, vector<int> pids, string optimization) {
  shared_ptr<route> original_route = dijkstra(source, destination);
    cout << "Length original: " << original_route->links.size() << endl;
    cout << "d: " << pids.size() << endl;
  cout << "doing another dijkstra" << endl;
  auto route_dijk = dijkstra(source, destination, original_route);  // checkup
  auto score = psychological_model.score_route(route_dijk->a(), route_dijk->b(), original_route->a(), original_route->b(), 0 , 0, pids.size());
  cout << "Score other dijkstra: " << score.first << " (" << score.second << ")" << endl;

  auto start = chrono::steady_clock::now();
  pair<shared_ptr<route>, double> ssotd_res =
      ssotd_route(source, destination, original_route, pids.size(), optimization);
  auto end = chrono::steady_clock::now();
  cout << "time used: " << chrono::duration_cast<chrono::microseconds>(end - start).count() << endl;
  double usage = ssotd_res.second / static_cast<double>(pids.size());
  cout << "normalized usage of the pareto route: " << usage << endl;
  for (int pid : pids)
    if ((rand() % (1 << 16)) / static_cast<double>(1 << 16) < usage)
      persons[pid].r = ssotd_res.first;
    else
      persons[pid].r = original_route;
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
}
