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

#define MAX_ORIGINAL_ROUTE_NODES 1024



class RouteFragment {
 public:
  virtual void add_to(vector<link*>&) const = 0;
  virtual double a() = 0;
  virtual double b() = 0;
  virtual double taud() = 0;
  virtual double shared_a() { return a(); }  // by default everything shared
  virtual double shared_b() { return b(); }
  virtual double shared_taud() { return taud(); }
  virtual double latency(int k) { return psychological_model.latency(a(), b(), k); }
  virtual double slatency(int k) { return psychological_model.latency(shared_a(), shared_b(), k); }
  virtual void add_shared_links(vector<link*>& links) { (void)links; }
  shared_ptr<ParetoElement> to_par_elem() {
  return make_shared<ParetoElement>(a(), b(), taud(), shared_a(), shared_b(), shared_taud());
}
  bool strongly_dominating(shared_ptr<RouteFragment>& other) {
    return psychological_model.strongly_dominating(this->to_par_elem(), other->to_par_elem());
  }
};

class EmptyRouteFragment : public RouteFragment {
 public:
  EmptyRouteFragment() = default;
  virtual void add_to(vector<link*>&) const override {}
  virtual double a() override { return 0.0; }
  virtual double b() override { return 0.0; }
  virtual double taud() override { return 0.0; }
};

class CompositeRouteFragment : public RouteFragment {
 protected:
  double _a = -1.0;
  double _b = -1.0;
  double _taud = -1.0;
  double _sa = -1.0;
  double _sb = -1.0;
  double _staud = -1.0;
  vector<shared_ptr<RouteFragment>>
      components;  // TODO try to replace with 2 explicit components (variables)
                   // maybe this vector is just too much overhead

 public:
  CompositeRouteFragment(vector<shared_ptr<RouteFragment>>& fragments) {
    copy(fragments.begin(), fragments.end(), back_inserter(components));
  }
  CompositeRouteFragment(shared_ptr<RouteFragment> head, shared_ptr<RouteFragment> appendix) {
    components.push_back(head);
    components.push_back(appendix);
    _a = head->a() + appendix->a();
    _b = head->b() + appendix->b();
    _taud = head->taud() + appendix->taud();
    _sa = head->shared_a() + appendix->shared_a();
    _sb = head->shared_b() + appendix->shared_b();
    _staud = head->shared_taud() + appendix->shared_taud();
  }
  virtual void add_to(vector<link*>& links) const override {
    for (const shared_ptr<RouteFragment> rf : components)
      rf->add_to(links);
  }
  virtual double a() override {
    if (_a == -1.0)
      _a = accumulate(components.begin(), components.end(), 0.0,
                      [](double su, const shared_ptr<RouteFragment> rf) { return su + rf->a(); });
    return _a;
  }
  virtual double b() override {
    if (_b == -1.0)
      _b = accumulate(components.begin(), components.end(), 0.0,
                      [](double su, const shared_ptr<RouteFragment> rf) { return su + rf->b(); });
    return _b;
  }
  virtual double taud() override {
    if (_taud == -1.0)
      _taud = accumulate(components.begin(), components.end(), 0.0,
                      [](double su, const shared_ptr<RouteFragment> rf) { return su + rf->taud(); });
    return _taud;
  }
  virtual double shared_a() override {
    if (_sa == -1.0)
      _sa = accumulate(
          components.begin(), components.end(), 0.0,
          [](double su, const shared_ptr<RouteFragment> rf) { return su + rf->shared_a(); });
    return _sa;
  }
  virtual double shared_b() override {
    if (_sb == -1.0)
      _sb = accumulate(
          components.begin(), components.end(), 0.0,
          [](double su, const shared_ptr<RouteFragment> rf) { return su + rf->shared_b(); });
    return _sb;
  }
  virtual double shared_taud() override {
    if (_staud == -1.0)
      _staud = accumulate(components.begin(), components.end(), 0.0,
                      [](double su, const shared_ptr<RouteFragment> rf) { return su + rf->shared_taud(); });
    return _staud;
  }
};

class LinkFragment : public RouteFragment {
 protected:
  link* l;
  int origIdx;

 public:
  LinkFragment(link* mLink, int mOrigIdx) : l(mLink), origIdx(mOrigIdx) {}
  virtual void add_to(vector<link*>& links) const override { links.push_back(l); }
  virtual double a() override { return l->a(); }
  virtual double b() override { return l->b(); }
  virtual double taud() override { return l->taud(); }
  virtual double shared_a() override { return l->a(); }
  virtual double shared_b() override { return l->b(); }
  virtual double shared_taud() override { return l->taud(); }
};

class ParetoElementFragment : public RouteFragment {
 protected:
  shared_ptr<ParetoElement> p;

 public:
  ParetoElementFragment(shared_ptr<ParetoElement> pe) : p(pe) {}
  virtual void add_to(vector<link*>& links) const override {
    auto pLinks = p->collectLinks();
    for (link* l : *pLinks)
      links.push_back(l);
  }
  virtual double a() override { return p->a(); }
  virtual double b() override { return p->b(); }
  virtual double taud() override { return p->taud(); }
  virtual double shared_a() override { return 0.0; }
  virtual double shared_b() override { return 0.0; }
  virtual double shared_taud() override { return 0.0; }
};

pair<double, double> lower_bound_score_zero(shared_ptr<ParetoElement> par, int from, int to, int c,
                                            int v, shared_ptr<route> original_route, int k) {
  (void) from;  (void) to;
  auto score = psychological_model.score_route(par->a(), par->b(), original_route->a(),
                                               original_route->b(), 0, 0, k);
  if (score.second > 0)
    return make_pair(score.first, score_for_relax(index_in_original(c), index_in_original(v), par, k));
  return make_pair(HUGE_VAL, -1);
}

pair<double, double> lower_bound_score_airline(shared_ptr<ParetoElement> par, int from, int to,
                                               int c, int v, shared_ptr<route> original_route, int k) {
  auto [a, b] = get_lazy_airval(v, to);
  auto [af, bf] = get_lazy_airval_forward(c, from);
  
  auto score = psychological_model.score_route(par->a() + a + af, par->b() + b + bf,
                                               original_route->a(), original_route->b(), 0, 0, k);
                               
  if (score.second > 0)
    return make_pair(score.first, score_for_relax(index_in_original(c), index_in_original(v), par, k));
  return make_pair(HUGE_VAL, -1);
}

pair<double, double> lower_bound_score_dijkstra(shared_ptr<ParetoElement> par, int from, int to,
                                                int c, int v, shared_ptr<route> original_route,
                                                int k) {
  (void) from;  (void) to;
  double newA = par->a() + bestAsForward[c] + bestAs[v];
  double newB = par->b() + bestBsForward[c] + bestBs[v];
  auto score = psychological_model.score_route(newA, newB, original_route->a(), original_route->b(),
                                               0, 0, k);
  if (score.second > 0)
    return make_pair(score.first, score_for_relax(index_in_original(c), index_in_original(v), par, k));
  return make_pair(HUGE_VAL, -1);
}

bool insert_and_dominate(list<shared_ptr<RouteFragment>>& A, shared_ptr<RouteFragment> frag) {
  bool appended = false;
  auto toBeDeleted = A.end();
  for (auto current_elem = A.begin(); current_elem != A.end(); current_elem++) {
    if (toBeDeleted != A.end()) {
      A.erase(toBeDeleted);
      toBeDeleted = A.end();
    }
    if ((*current_elem)->strongly_dominating(frag)) {
      A.push_front(*current_elem);
      A.erase(current_elem);
      return false;
    } else if (frag->strongly_dominating(*current_elem)) {
      if (appended) {
        toBeDeleted = current_elem;
      } else {
        *current_elem = frag;
        appended = true;
      }
    }
  }
  if (!appended)
    A.push_back(frag);
  return true;
}

pair<shared_ptr<route>, double> ssotd_route(int a, int b, shared_ptr<route> original_route, int k,
                                            string optimization) {
  map<pair<int, int>, vector<shared_ptr<ParetoElement>>> paretoFronts;
  unordered_map<int, bool> inactive;
  prepare_original_route(original_route, inactive);
  double qot = k * psychological_model.latency(original_route->a(), original_route->b(), k);

  cout << "DIJKSTRA OT: " << qot << endl;
  cout << "Calculating pareto fronts." << endl;
  function<pair<double,int>(int, vector<vector<shared_ptr<ParetoElement>>>*)> pareto_dijk;
  double upperBound = qot;

  if ("none" == optimization || "hop" == optimization) {
    cout << "Doing no optimization" << endl;
    pareto_dijk = [&inactive, &original_route, upperBound, b](int c, vector<vector<shared_ptr<ParetoElement>>>* pareto) {
      auto vis = pareto_dijsktra(c, b, *pareto, original_route, inactive);
      return make_pair(upperBound, vis);
    };
  } else if ("simple_local_opt" == optimization) {
    cout << "Doing simple optimization" << endl;
    pareto_dijk = [a, b, &original_route, k, &upperBound, &inactive](
                      int c, vector<vector<shared_ptr<ParetoElement>>>* pareto) {
      return pareto_dijkstra_local_opt(c, a, b, *pareto, original_route, k, upperBound, inactive, &lower_bound_score_zero);
    };
  } else if ("dijkstra_local_opt" == optimization) {
    cout << "Doing dijkstra optimization" << endl;
    auto start = chrono::steady_clock::now();
    fill_best_pars_dijkstra(b);
    fill_best_pars_dijkstra_forward(a);
    auto end = chrono::steady_clock::now();
    cout << "Route specific precalculation time: "
         << chrono::duration_cast<chrono::microseconds>(end - start).count() << endl;
    pareto_dijk = [a, b, &original_route, k, &upperBound, &inactive](
                      int c, vector<vector<shared_ptr<ParetoElement>>>* pareto) {
      return pareto_dijkstra_local_opt(c, a, b, *pareto, original_route, k, upperBound, inactive, &lower_bound_score_dijkstra);
    };
  } else if ("airline_local_opt" == optimization) {
    cout << "Doing airline optimization" << endl;
    cerr << "Assuming best freeflowspeed = 130 and best capacity = 2288 (sane-berlin.xml)!!!!"
         << endl;
    bestAs = vector<double>(nodes.size());
    bestBs = vector<double>(nodes.size());
    bestAsForward = vector<double>(nodes.size());
    bestBsForward = vector<double>(nodes.size());
    pareto_dijk = [a, b, &original_route, k, &upperBound, &inactive](
                      int c, vector<vector<shared_ptr<ParetoElement>>>* pareto) {
      return pareto_dijkstra_local_opt(c, a, b, *pareto, original_route, k, upperBound, inactive, &lower_bound_score_airline);
    };
   } else if ("airline_astar_opt" == optimization) {
    cout << "Doing airline astar optimization" << endl;
    cerr << "Assuming best freeflowspeed = 130 and best capacity = 2288 (sane-berlin.xml)!!!!"
         << endl;
    bestAs = vector<double>(nodes.size());
    bestBs = vector<double>(nodes.size());
    bestAsForward = vector<double>(nodes.size());
    bestBsForward = vector<double>(nodes.size());
    pareto_dijk = [a, b, &original_route, k, &upperBound, &inactive](
                      int c, vector<vector<shared_ptr<ParetoElement>>>* pareto) {
      return pareto_dijkstra_local_opt(c, a, b, *pareto, original_route, k, upperBound, inactive, &lower_bound_score_airline, &astar_prio_airline);
    };
  } else if ("dijkstra_astar_opt" == optimization) {
     cout << "Doing dijkstra astar optimization" << endl;
    auto start = chrono::steady_clock::now();
    fill_best_pars_dijkstra(b);
    fill_best_pars_dijkstra_forward(a);
    auto end = chrono::steady_clock::now();
    cout << "Route specific precalculation time: "
         << chrono::duration_cast<chrono::microseconds>(end - start).count() << endl;
    pareto_dijk = [a, b, &original_route, k, &upperBound, &inactive](
                      int c, vector<vector<shared_ptr<ParetoElement>>>* pareto) {
      return pareto_dijkstra_local_opt(c, a, b, *pareto, original_route, k, upperBound, inactive, &lower_bound_score_dijkstra, &astar_prio_dijkstra);
    };
  } else
    cerr << "Invalid parameter for optimization given: " << optimization << endl;

  long long visits = 0;
  auto start = chrono::steady_clock::now();
#pragma omp parallel for default(none) shared(paretoFronts, original_route, adj, inactive, upperBound, k, \
                                              b, pareto_dijk, visits) schedule(dynamic, 2) num_threads(8)
  for (unsigned int lid = 0; lid < original_route->links.size(); lid++) {
    // iterate over all vertices of the original route except the last
    int v = original_route->links[lid]->from;
    vector<vector<shared_ptr<ParetoElement>>> pareto(adj.size() + 1);
    // vector<list<shared_ptr<ParetoElement>>> pareto(adj.size() +1);
    auto [newUpperbound, new_visits] = pareto_dijk(v, &pareto);
    visits += new_visits;
    if (newUpperbound < upperBound)
      upperBound = newUpperbound;


#pragma omp critical
    for (unsigned int _lid = lid + 1; _lid < original_route->links.size(); _lid++) {
      paretoFronts[{lid, _lid}] = pareto[original_route->links[_lid]->from];
    }
    paretoFronts[{lid, original_route->links.size()}] = pareto[original_route->links.back()->to];
  }
  auto end = chrono::steady_clock::now();
  cout << "Node visits: " << visits << endl;
  cout << "Pareto-dijkstra time: " << chrono::duration_cast<chrono::microseconds>(end - start).count() << endl;

  // DP
  cout << "starting arbitrary disjoint dp" << endl;


  stringstream statsCsv;

  start = chrono::steady_clock::now();

  int counter = 0;
  int n = original_route->links.size() + 1;
  auto pareto_sizes = vector<int>(n*n/2-n);  //n^2/2 - n Pareto-Fronts
  vector<list<shared_ptr<RouteFragment>>> A(original_route->links.size() + 1);
  A[0].push_back(make_shared<EmptyRouteFragment>());
  for (size_t i = 1; i <= original_route->links.size(); i++) {
    auto appendix = make_shared<LinkFragment>(original_route->links[i - 1], i - 1);
    for (auto& frag : A[i - 1]) {
      counter++;
      auto newFrag = make_shared<CompositeRouteFragment>(frag, appendix);
      insert_and_dominate(A[i], newFrag);
    }
    for (size_t j = 0; j < i; j++) {
      pareto_sizes.push_back(paretoFronts[{j, i}].size());
      for (auto& frag : A[j]) {
        for (auto& bridge : paretoFronts[{j, i}]) {
          counter++;
          auto bridgeFragment = make_shared<ParetoElementFragment>(bridge);
          auto newFrag = make_shared<CompositeRouteFragment>(frag, bridgeFragment);  // create copy
          insert_and_dominate(A[i], newFrag);
        }
      }
    }

    // // statistics recording
    // for (auto& frag : A[i]) {
    //   auto [ot, usage] = psychological_model.score_route(
    //       frag->a(), frag->b(), origPartA[i], origPartB[i], frag->shared_a(), frag->shared_b(), k);
    //   (void)usage;
    //   // scale down ot by length
    //   ot /= origTt[i];
    //   statsCsv << i << "," << ot << "," << frag->b() / frag->a() << ","
    //             << frag->shared_b() / frag->shared_a() << "\n";
  // }
  }
  double best_ot = numeric_limits<double>::max();
  double best_usage = 0.0;
  shared_ptr<RouteFragment> best;

  for (shared_ptr<RouteFragment>& frag : A[original_route->links.size()]) {
    auto [ot, usage] =
      psychological_model.score_route(frag->a(), frag->b(), original_route->a(),
                                      original_route->b(), frag->shared_a(), frag->shared_b(), k);
    if (ot < best_ot) {
      best_ot = ot;
      best_usage = usage;
      best = frag;
    }
  } 
  end = chrono::steady_clock::now();
  cout << "Evaluation time: "
        << chrono::duration_cast<chrono::microseconds>(end - start).count() << endl;
    
    auto total_paretosizes = accumulate(pareto_sizes.begin(), pareto_sizes.end(), 0);
    nth_element(pareto_sizes.begin(), pareto_sizes.begin() + pareto_sizes.size() / 2, pareto_sizes.end());
    auto mean_pareto_set_size = pareto_sizes[pareto_sizes.size()/2];
    auto A_sizes = vector<int>(A.size());
    for (auto elem : A)
      A_sizes.push_back(elem.size());
    auto total_dp = accumulate(A_sizes.begin(), A_sizes.end(), 0);
    nth_element(A_sizes.begin(), A_sizes.begin() + A_sizes.size() / 2, A_sizes.end());
    auto mean_dp_set_size = A_sizes[A_sizes.size()/2];
    cout << "Mean Pareto-set size: " << mean_pareto_set_size << endl;
    cout << "Sum Pareto-set size: " << total_paretosizes << endl;
    cout << "Mean DP-set size: " << mean_dp_set_size << endl;
    cout << "Sum DP-set size: " << total_dp << endl;

    cout << "Found " << A[original_route->links.size()].size() << " pareto-optimal routes" << endl;
    cout << "evaluated " << counter << " pareto parts" << endl;

  cout << "\nBEST PARETO OT: " << best_ot << endl;
  cout << "\nSELECTED ALTERNATIVE: a=" << best->a() << " b=" << best->b()
       << " sa=" << best->shared_a() << " sb=" << best->shared_b() << endl;
  cout << "b/a=" << best->b() / best->a() << endl;
  if (best_ot > qot)
    return {original_route, 0.0};

  vector<link*> altLinks;
  best->add_to(altLinks);
  auto res = make_shared<route>(altLinks);
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
