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
#include <queue>

#include "core/data.h"
#include "core/globals.h"
#include "core/routing.h"
#include "ssotd/ssotd_core.h"

using namespace std;

template <class T>
using minq = priority_queue<T, vector<T>, greater<>>;


shared_ptr<route> dijkstra_all(int a, int b, int k) {
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
      double newDist = d + l->latency(k);
      if (newDist < dist[l->to]) {
        dist[l->to] = newDist;
        q.push({newDist, l->to});
        prec[l->to] = {cur, l};
      }
    }
  }
  if (prec[b].first == -1) {
    cerr << "WARNING!" << endl;
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

void route(int source, int destination, vector<int> pids) {
  auto start = chrono::steady_clock::now();
  auto original_route = dijkstra_all(source, destination, pids.size());
  auto end = chrono::steady_clock::now();
  cout << "time used: " << chrono::duration_cast<chrono::microseconds>(end - start).count() << endl;
  cout << "Length original: " << original_route->links.size() << endl;
    cout << "K: " << pids.size() << endl;
  auto score = pids.size() * psychological_model.latency(original_route->a(), original_route->b(), pids.size());

  cout << "normalized usage of the pareto route: 0" << endl;
  cout << "\nBEST PARETO OT: " << score << endl;
  for (int pid : pids)
      persons[pid].r = original_route;
}

void do_routing(int argc, char* argv[]) {
    (void) argc;
    (void) argv;

  map<pair<pair<int, int>, string>, vector<int>> c;
  for (unsigned int pid = 0; pid < persons.size(); pid++) {
    auto& p = persons[pid];
    string s = p.timestr;  // maybe
    c[{{p.origin, p.destination}, s}].push_back(pid);
  }
  // do ssotd for all
  for (auto& [sdts, pv] : c) {
    number_agents = pv.size();
    route(sdts.first.first, sdts.first.second, pv);
  }
}
