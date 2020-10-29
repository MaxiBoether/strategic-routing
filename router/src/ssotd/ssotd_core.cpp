#include "ssotd/ssotd_core.h"

#include <algorithm>
#include <any>
#include <chrono>
#include <cmath>
#include <functional>
#include <iostream>
#include <list>
#include <memory>
#include <numeric>
#include <queue>
#include <vector>

#include "core/data.h"
#include "core/globals.h"

template <class T>
using minq = priority_queue<T, vector<T>, greater<>>;
using ll = long long;

using namespace std;

int number_agents;

static double mean_a = 0.0026179770605695984;
static double mean_b = 39.811862162760164;


vector<double> bestAs, bestBs, bestAsForward, bestBsForward; // For dijkstra & airline local opt

vector<double> origTt, origPartA, origPartB;  // original route prefix sums

unordered_map<int, int> nodes_original_route;  // maps a node id to its index in orig route

int to_node;
double max_sharedA;
double mean_taud;
shared_ptr<route> orig_path;

shared_ptr<route> dijkstra(int a, int b, shared_ptr<route> original_route) {
  unordered_map<int, bool> inactive;
  if (original_route)
    for_each(original_route->links.begin(), original_route->links.end(),
             [&inactive](link* l) { inactive[l->id] = true; });
  vector<double> dist(adj.size(), HUGE_VAL);
  vector<pair<int, link*>> prec(adj.size(), {-1, nullptr});
  minq<pair<double, int>> q;
  dist[a] = 0.0f;
  q.push({0.0f, a});
  while (!q.empty()) {
    auto [d, cur] = q.top();
    q.pop();
    if (cur == b)
      break;
    if (d > dist[cur])
      continue;
    for (link* l : adj[cur]) {
      if (original_route && inactive[l->id])
        continue;
      double newDist = d + l->b();
      if (newDist < dist[l->to]) {
        dist[l->to] = newDist;
        q.push({newDist, l->to});
        prec[l->to] = {cur, l};
      }
    }
  }
  if (prec[b].first == -1) {
    cerr << "HOLY SH*T THE HOUSE IS ON FIRE" << endl;
    cerr << "(djikstra) could not find any route from " << a << " to " << b << endl;
    exit(1);
  }

  vector<link*> newRt;
  int cur = b;
  while (prec[cur].first != -1) {
    newRt.push_back(prec[cur].second);
    cur = prec[cur].first;
  }
  reverse(newRt.begin(), newRt.end());
  shared_ptr<route> r = make_shared<route>(newRt);
  return r;
}


shared_ptr<vector<double>> dijkstra_for_opt(int v, bool doA, unordered_map<int, bool> inactive, vector<vector<link*>> adj2, bool forward) {
  auto dist = make_shared<vector<double>>(adj2.size(), HUGE_VAL);
  minq<pair<double, int>> q;
  (*dist)[v] = 0.0f;
  q.push({0.0f, v});
  while (!q.empty()) {
    auto [d, cur] = q.top();
    q.pop();
    if (d > (*dist)[cur])
      continue;
    for (link* l : adj2[cur]) {
      if (inactive[l->id])
        continue;
      double newDist = doA ? d + l->a() : d + l->b();
      if (newDist < (*dist)[forward ? l->to : l->from]) {
        (*dist)[forward ? l->to : l->from] = newDist;
        q.push({newDist, forward ? l->to : l->from});
      }
    }
  }
  return dist;
}

bool standard_prio(pair<shared_ptr<ParetoElement>, int> left, pair<shared_ptr<ParetoElement>, int> right) {
  return right.first->k() < left.first->k();
}

bool astar_prio_dijkstra2(pair<shared_ptr<ParetoElement>, int> left, pair<shared_ptr<ParetoElement>, int> right) {
  double leftTaud = (left.first->taud() + psychological_model.latency(bestAs[left.second], bestBs[left.second], number_agents)) / mean_taud;
  double leftB = (left.first->b() + bestBs[left.second]) / mean_b;
  double rightTaud = (right.first->taud() + psychological_model.latency(bestAs[right.second], bestBs[right.second], number_agents)) / mean_taud;
  double rightB = (right.first->b() + bestBs[right.second]) / mean_b;
  return leftTaud + leftB + left.first->shared_a() * 2 > rightTaud + rightB + right.first->shared_a() * 2;
}

bool astar_prio_dijkstra(pair<shared_ptr<ParetoElement>, int> left, pair<shared_ptr<ParetoElement>, int> right) {
  double leftA = left.first->a() + bestAs[left.second];
  double leftB = left.first->b() + bestBs[left.second];
  double rightA = right.first->a() + bestAs[right.second];
  double rightB = right.first->b() + bestBs[right.second];
  return psychological_model.score_route(leftA, leftB, orig_path->a(), orig_path->b(), left.first->shared_a(), left.first->shared_b(), number_agents) > 
  psychological_model.score_route(rightA, rightB, orig_path->a(), orig_path->b(), right.first->shared_a(), right.first->shared_b(), number_agents);
}

bool astar_prio_airline2(pair<shared_ptr<ParetoElement>, int> left, pair<shared_ptr<ParetoElement>, int> right) {
  (void) left;
  (void) right;
 return false;
}

bool astar_prio_airline(pair<shared_ptr<ParetoElement>, int> left, pair<shared_ptr<ParetoElement>, int> right) {
  auto [leftA, leftB] = get_lazy_airval(left.second, to_node);
  leftA += left.first->a();
  leftB += left.first->b();
  auto [rightA, rightB] = get_lazy_airval(right.second, to_node);
  rightA += right.first->a();
  rightB += right.first->b();
  return psychological_model.score_route(leftA, leftB, orig_path->a(), orig_path->b(), left.first->shared_a(), left.first->shared_b(), number_agents) > 
  psychological_model.score_route(rightA, rightB, orig_path->a(), orig_path->b(), right.first->shared_a(), right.first->shared_b(), number_agents);
}

bool can_ignore(shared_ptr<ParetoElement> par, link* l) {
  return par->myLink && par->myLink->from == l->to;
}

double min_score(int u, int to, shared_ptr<ParetoElement> par, int k) {
  auto [a, b] = get_lazy_airval(u, to);
  a += par->a();
  b += par->b();
  return psychological_model.score_route(a, b, orig_path->a(), orig_path->b(), par->shared_a(), par->shared_b(), k).second;
}

pair<double, ll> pareto_dijkstra_local_opt(int a, int from, int to, vector<vector<shared_ptr<ParetoElement>>>& pareto,
                               shared_ptr<route> original_route, int k, double qot,
                               unordered_map<int, bool> inactive,
                               pair<double, double> (*lower_bound_score)(shared_ptr<ParetoElement>, int, int, int, int,
                                                           shared_ptr<route>, int),
                     bool (*prio)(pair<shared_ptr<ParetoElement>, int>, pair<shared_ptr<ParetoElement>, int>)) {
  cout << "Finding pareto routes for " << a << "  using qot  " << qot << endl;
  to_node = to;
  max_sharedA = original_route->a();
  mean_taud = psychological_model.latency(mean_a, mean_b, number_agents);
  orig_path = original_route;
  auto cmp = [prio](pair<shared_ptr<ParetoElement>, int> left,
                pair<shared_ptr<ParetoElement>, int> right) {
    return (*prio)(left, right);
  };
  priority_queue<pair<shared_ptr<ParetoElement>, int>,
                 std::vector<pair<shared_ptr<ParetoElement>, int>>, decltype(cmp)>
      q(cmp);
  ll visits = 0;
  auto zero_el = make_shared<ParetoElement>();
  q.push({zero_el, a});
  while (!q.empty()) {
    visits++;
    auto [par, u] = q.top();
    q.pop();
    // if (qot < min_score(u, to, par, k))  //Does not work!  + Remove if not used with A*-prio ordered by score_route value !!!!!!!!!!
    //   break;
    for (link* l : adj[u]) {
      if (can_ignore(par,l) || inactive[l->id])
        continue;
      int v = l->to;
      auto newPar = make_shared<ParetoElement>(par, l);
      pair<double, double> ot = lower_bound_score(newPar, from, to, a, v, original_route, k);

      if (ot.first > qot + 100) {
        continue;
      }
      if (ot.second > 0 && ot.second < qot) {
        cout << "relaxed ot cap" << endl;
        qot = ot.second;
      }
      if (!any_of(pareto[v].begin(), pareto[v].end(), [&newPar](shared_ptr<ParetoElement>& vpar) {
            return psychological_model.dominating(vpar, newPar);
          })) {
        pareto[v].erase(remove_if(pareto[v].begin(), pareto[v].end(),
                                  [&newPar](shared_ptr<ParetoElement>& vpar) {
                                    return psychological_model.dominating(newPar, vpar);
                                  }),
                        pareto[v].end());
        pareto[v].push_back(newPar);
        q.push({newPar, v});
      }
    }
  }
  return make_pair(qot, visits);
}

pair<double, ll> pareto_dijkstra_local_opt_4d(int a, int from, int to, vector<vector<shared_ptr<ParetoElement>>>& pareto,
                               shared_ptr<route> original_route, int k, double qot, unordered_map<int, bool> is_orig_edge,
                               pair<double, double> (*lower_bound_score)(shared_ptr<ParetoElement>, int, int, int, int,
                                                           shared_ptr<route>, int),
                     bool (*prio)(pair<shared_ptr<ParetoElement>, int>, pair<shared_ptr<ParetoElement>, int>)) {
  cout << "Finding pareto routes for " << a << "  using qot  " << qot << endl;
  to_node = to;
  max_sharedA = original_route->a();
mean_taud = psychological_model.latency(mean_a, mean_b, number_agents);
orig_path = original_route;
  auto cmp = [prio](pair<shared_ptr<ParetoElement>, int> left,
                pair<shared_ptr<ParetoElement>, int> right) {
    return (*prio)(left, right);
  };
  priority_queue<pair<shared_ptr<ParetoElement>, int>,
                 std::vector<pair<shared_ptr<ParetoElement>, int>>, decltype(cmp)>
      q(cmp);
  ll visits = 0;
  auto zero_el = make_shared<ParetoElement>();
  q.push({zero_el, a});
  while (!q.empty()) {
    visits++;
    auto [par, u] = q.top();
    q.pop();
    // if (qot < min_score(u, to, par, k))  //Does not work!  + Remove if not used with A*-prio ordered by score_route value !!!!!!!!!!
    //   break;
    for (link* l : adj[u]) {
      if (can_ignore(par,l))
	  continue;
      int v = l->to;
      auto newPar = make_shared<ParetoElement>(par, l, is_orig_edge[l->id]);
      pair<double, double> ot = lower_bound_score(newPar, from, to, a, v, original_route, k);

      if (ot.first > qot + 100) {
        continue;
      }
      if (ot.second > 0 && ot.second < qot) {
        cout << "relaxed ot cap" << endl;
        qot = ot.second;
      }
      if (!any_of(pareto[v].begin(), pareto[v].end(), [&newPar](shared_ptr<ParetoElement>& vpar) {
            return psychological_model.strongly_dominating(vpar, newPar);
          })) {
        pareto[v].erase(remove_if(pareto[v].begin(), pareto[v].end(),
                                  [&newPar](shared_ptr<ParetoElement>& vpar) {
                                    return psychological_model.strongly_dominating(newPar, vpar);
                                  }),
                        pareto[v].end());
        pareto[v].push_back(newPar);
        q.push({newPar, v});
      }
    }
  }
  return make_pair(qot, visits);
}


pair<double, ll> pareto_dijkstra_local_opt_4d_1D(int a, int from, int to, vector<vector<shared_ptr<ParetoElement>>>& pareto,
                               shared_ptr<route> original_route, int k, double qot, unordered_map<int, bool> is_orig_edge,
                               pair<double, double> (*lower_bound_score)(shared_ptr<ParetoElement>, int, int, int, int,
                                                           shared_ptr<route>, int),
                     bool (*prio)(pair<shared_ptr<ParetoElement>, int>, pair<shared_ptr<ParetoElement>, int>)) {
  cout << "Finding pareto routes for " << a << "  using qot  " << qot << endl;
  to_node = to;
  max_sharedA = original_route->a();
mean_taud = psychological_model.latency(mean_a, mean_b, number_agents);
orig_path = original_route;

  auto cmp = [prio](pair<shared_ptr<ParetoElement>, int> left,
                pair<shared_ptr<ParetoElement>, int> right) {
    return (*prio)(left, right);
  };
  priority_queue<pair<shared_ptr<ParetoElement>, int>,
                 std::vector<pair<shared_ptr<ParetoElement>, int>>, decltype(cmp)>
      q(cmp);
  ll visits = 0;
  auto zero_el = make_shared<ParetoElement>();
  q.push({zero_el, a});
  while (!q.empty()) {
    visits++;
    auto [par, u] = q.top();
    q.pop();
    // if (qot < min_score(u, to, par, k))  //Does not work!  + Remove if not used with A*-prio ordered by score_route value !!!!!!!!!!
    //   break;
    for (link* l : adj[u]) {
      if (can_ignore(par,l) || (par->hasSplit && !is_orig_edge[l->id]))
        continue;
      int v = l->to;
      auto newPar = make_shared<ParetoElement>(par, l, is_orig_edge[l->id]);
 
      if (par->hasSplit || (is_orig_edge[l->id] && (par->myLink && !is_orig_edge[par->myLink->id])))
        newPar->hasSplit = true;
      pair<double, double> ot = lower_bound_score(newPar, from, to, a, v, original_route, k);

      if (ot.first > qot + 100) {
        continue;
      }
      if (ot.second > 0 && ot.second < qot) {
        cout << "relaxed ot cap" << endl;
        qot = ot.second;
      }
      if (!any_of(pareto[v].begin(), pareto[v].end(), [&newPar](shared_ptr<ParetoElement>& vpar) {
            return psychological_model.strongly_dominating(vpar, newPar);
          })) {
        pareto[v].erase(remove_if(pareto[v].begin(), pareto[v].end(),
                                  [&newPar](shared_ptr<ParetoElement>& vpar) {
                                   return psychological_model.strongly_dominating(newPar, vpar);
                                  }),
                        pareto[v].end());
        pareto[v].push_back(newPar);
        q.push({newPar, v});
      }
    }
  }
  return make_pair(qot, visits);
}

ll pareto_dijsktra(int a, int b, vector<vector<shared_ptr<ParetoElement>>>& pareto, shared_ptr<route> original_route, 
                     unordered_map<int, bool> inactive,
                     bool (*prio)(pair<shared_ptr<ParetoElement>, int> left, pair<shared_ptr<ParetoElement>, int> right)) {
  cout << "Finding pareto routes for " << a << endl;
  to_node = b;
  max_sharedA = original_route->a();
  mean_taud = psychological_model.latency(mean_a, mean_b, number_agents);
orig_path = original_route;
  auto cmp = [prio](pair<shared_ptr<ParetoElement>, int> left,
                pair<shared_ptr<ParetoElement>, int> right) {
    return (*prio)(left, right);
  };
  priority_queue<pair<shared_ptr<ParetoElement>, int>,
                 std::vector<pair<shared_ptr<ParetoElement>, int>>, decltype(cmp)>
      q(cmp);
  auto zero_el = make_shared<ParetoElement>();
  q.push({zero_el, a});
  ll visits = 0;
  while (!q.empty()) {
    visits++;
    auto [par, u] = q.top();
    q.pop();
    for (link* l : adj[u]) {
      if (can_ignore(par,l) || inactive[l->id])
        continue;
      int v = l->to;
      auto newPar = make_shared<ParetoElement>(par, l);
      if (!any_of(pareto[v].begin(), pareto[v].end(), [&newPar](shared_ptr<ParetoElement>& vpar) {
            return psychological_model.dominating(vpar, newPar);
          })) {
        pareto[v].erase(remove_if(pareto[v].begin(), pareto[v].end(),
                                  [&newPar](shared_ptr<ParetoElement>& vpar) {
                                    return psychological_model.dominating(newPar, vpar);
                                  }),
                        pareto[v].end());
        pareto[v].push_back(newPar);
        q.push({newPar, v});
      }
    }
  }
  return visits;
}

ll pareto_dijsktra_4d(int a, int b, vector<vector<shared_ptr<ParetoElement>>>& pareto, shared_ptr<route> original_route,
                     unordered_map<int, bool> is_orig_edge,
                     bool (*prio)(pair<shared_ptr<ParetoElement>, int> left, pair<shared_ptr<ParetoElement>, int> right)) {
  cout << "Finding pareto routes for " << a << endl;
  to_node = b;
  max_sharedA = original_route->a();
mean_taud = psychological_model.latency(mean_a, mean_b, number_agents);
orig_path = original_route;
  auto cmp = [prio](pair<shared_ptr<ParetoElement>, int> left,
                pair<shared_ptr<ParetoElement>, int> right) {
    return (*prio)(left, right);
  };
  priority_queue<pair<shared_ptr<ParetoElement>, int>,
                 std::vector<pair<shared_ptr<ParetoElement>, int>>, decltype(cmp)>
      q(cmp);
  auto zero_el = make_shared<ParetoElement>();
  q.push({zero_el, a});
  ll visits = 0;
  while (!q.empty()) {
    visits++;
    auto [par, u] = q.top();
    q.pop();
    for (link* l : adj[u]) {
      if (can_ignore(par,l))
	  continue;
      int v = l->to;
      auto newPar = make_shared<ParetoElement>(par, l, is_orig_edge[l->id]);
      if (!any_of(pareto[v].begin(), pareto[v].end(), [&newPar](shared_ptr<ParetoElement>& vpar) {
            return psychological_model.strongly_dominating(vpar, newPar);
          })) {
        pareto[v].erase(remove_if(pareto[v].begin(), pareto[v].end(),
                                  [&newPar](shared_ptr<ParetoElement>& vpar) {
                                    return psychological_model.strongly_dominating(newPar, vpar);
                                  }),
                        pareto[v].end());
        pareto[v].push_back(newPar);
        q.push({newPar, v});
      }
    }
  }
  return visits;
}


void pareto_dijsktra_4d_1D(int a, int b, vector<vector<shared_ptr<ParetoElement>>>& pareto, shared_ptr<route> original_route, 
                     unordered_map<int, bool> is_orig_edge,
                     bool (*prio)(pair<shared_ptr<ParetoElement>, int> left, pair<shared_ptr<ParetoElement>, int> right)) {
  cout << "Finding pareto routes for " << a << endl;
  to_node = b;
  max_sharedA = original_route->a();
mean_taud = psychological_model.latency(mean_a, mean_b, number_agents);
orig_path = original_route;
  auto cmp = [prio](pair<shared_ptr<ParetoElement>, int> left,
                pair<shared_ptr<ParetoElement>, int> right) {
    return (*prio)(left, right);
  };
  priority_queue<pair<shared_ptr<ParetoElement>, int>,
                 std::vector<pair<shared_ptr<ParetoElement>, int>>, decltype(cmp)>
      q(cmp);
  auto zero_el = make_shared<ParetoElement>();
  q.push({zero_el, a});
  while (!q.empty()) {
    auto [par, u] = q.top();
    q.pop();
    for (link* l : adj[u]) {
      if (can_ignore(par, l) || (par->hasSplit && !is_orig_edge[l->id]))
        continue;
      int v = l->to;
      auto newPar = make_shared<ParetoElement>(par, l, is_orig_edge[l->id]);
 
      if (par->hasSplit || (is_orig_edge[l->id] && (par->myLink && !is_orig_edge[par->myLink->id])))
        newPar->hasSplit = true;
      
      if (!any_of(pareto[v].begin(), pareto[v].end(), [&newPar](shared_ptr<ParetoElement>& vpar) {
            return psychological_model.strongly_dominating(vpar, newPar);;
          })) {
        pareto[v].erase(remove_if(pareto[v].begin(), pareto[v].end(),
                                  [&newPar](shared_ptr<ParetoElement>& vpar) {
                                    return psychological_model.strongly_dominating(newPar, vpar);
                                  }),
                        pareto[v].end());
        pareto[v].push_back(newPar);
        q.push({newPar, v});
      }
    }
  }
}

double airline_dist(int from, int to) {
  double fromX = atof(nodes[from]->x);
  double fromY = atof(nodes[from]->y);
  double toX = atof(nodes[to]->x);
  double toY = atof(nodes[to]->y);
  return sqrt((fromX - toX) * (fromX - toX) + (fromY - toY) * (fromY - toY));
}

link* artificial_link(int from, int to, double freespeed, double capacity) {
  link* l = new link;
  l->id = -77;
  l->from = from;
  l->to = to;
  l->length = airline_dist(from, to);
  l->capacity = capacity;
  l->freespeed = freespeed;
  return l;
}

pair<double, double> get_lazy_airval(int node, int to) {
  double b = bestBs[node];
  
  if (b == 0) {
    link* art_link = artificial_link(node, to, 130, 2288);
    bestAs[node] = art_link->a();
    bestBs[node] = art_link->b();
  }
  return make_pair(bestAs[node], bestBs[node]);
}

pair<double, double> get_lazy_airval_forward(int node, int from) {
  double b = bestBsForward[node];
  
  if (b == 0) {
    link* art_link = artificial_link(from, node, 130, 2288);
    bestAsForward[node] = art_link->a();
    bestBsForward[node] = art_link->b();
  }
  return make_pair(bestAsForward[node], bestBsForward[node]);
}

void fill_best_pars_dijkstra(int to, unordered_map<int, bool> inactive) {
  std::vector<std::vector<link*>> adj_inv(adj.size());
  for (unsigned int i=0; i<adj.size(); i++) {
      for (auto link : adj[i])
        adj_inv[link->to].push_back(link);
    }
    bestAs = *dijkstra_for_opt(to, true, inactive, adj_inv, false);
    bestBs = *dijkstra_for_opt(to, false, inactive, adj_inv, false);
}

void fill_best_pars_dijkstra_forward(int from, unordered_map<int, bool> inactive) {
  bestAsForward = *dijkstra_for_opt(from, true, inactive, adj, true);
  bestBsForward = *dijkstra_for_opt(from, false, inactive, adj, true);
}

int index_in_original(int v) {
  if (auto val = nodes_original_route.find(v); val != nodes_original_route.end()) {
    return val->second;
  }
  return -1;
}

int is_orig_node(int node, shared_ptr<route> orig) {  
if (orig->links[0]->from == node)
	return 0;
 for (unsigned int i=0;i<orig->links.size(); i++) {
	if (orig->links[i]->to ==node)
		return i+1;
 }
return -1;
}

void prepare_original_route(shared_ptr<route> original_route, unordered_map<int, bool>& inactive) {
  int count = 0;
  nodes_original_route[original_route->links[0]->from] = count++;
  for_each(original_route->links.begin(), original_route->links.end(),
           [&inactive, &count](link* l) {
             inactive[l->id] = true;
             nodes_original_route[l->to] = count++;
           });

  // original route prefix sums
  origTt = *new vector<double>(original_route->links.size() + 1, 0.0);
  origPartA = *new vector<double>(original_route->links.size() + 1, 0.0);
  origPartB = *new vector<double>(original_route->links.size() + 1, 0.0);
  for (size_t i = 1; i < original_route->links.size() + 1; i++) {
    origTt[i] = origTt[i - 1] + original_route->links[i - 1]->length;
    origPartA[i] = origPartA[i - 1] + original_route->links[i - 1]->a();
    origPartB[i] = origPartB[i - 1] + original_route->links[i - 1]->b();
  }
}

double score_for_relax(int idc, int idv, shared_ptr<ParetoElement> par, double k) {
  if (idv < 0)
    return -1;
  double shared_a = origPartA.at(idc) + origPartA.back() - origPartA.at(idv);
  double shared_b = origPartB.at(idc) + origPartB.back() - origPartB.at(idv);
  auto [score, usage] =
      psychological_model.score_route(par->a() + shared_a, par->b() + shared_b, origPartA.back(),
                                      origPartB.back(), shared_a + par->shared_a(), shared_b + par->shared_b(), k);
  return usage > 0 ? score + 10 : -1;
}


void check_route_sanity(route& r, string routeName) {
  int last_node = r.links[0]->from;
  for (link* l : r.links) {
    if (l->from != last_node) {
      std::cout << "invalid route " << routeName << ". HIER KOMMT DIE ROUTENPOLIZEI! "
                << l->from << " " << last_node << std::endl;
    }
    last_node = l->to;
  }
}
