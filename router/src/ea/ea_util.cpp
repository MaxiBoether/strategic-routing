#include "ea/ea_util.h"

#include <algorithm>
#include <random>
#include <unordered_map>

#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>

#include "core/data.h"
#include "core/globals.h"
#include "ea/ea_constants.h"
#include "ea/ea_globals.h"
#include "ea/ea_random.h"
#include "ea_shared/dijkstra_prio_queue.h"

auto randDijkstra(double usage, int origin, int destination, route* notPreferredLinks, int k)
    -> route {
  // auto start = std::chrono::steady_clock::now();

  dijkPQ pq(queue_size);
  std::vector<int> dist(nodes.size(), -1);
  std::vector<link*> a(nodes.size(), nullptr);
  std::vector<bool> visited(nodes.size(), false);

  std::unordered_map<link*, bool> forbidden_links;
  if (notPreferredLinks != nullptr) {
    for (link* l : notPreferredLinks->links) {
      forbidden_links[l] = true;
    }
  }

  pq.push({0, origin});
  dist[origin] = 0;
  // auto dijk_start = std::chrono::steady_clock::now();
  // std::vector<int> rand_times;
  while (!pq.empty()) {
    int u = pq.top().second;
    auto w = pq.top().first;
    visited[u] = true;
    pq.pop();

    if (__builtin_expect(dist[u] != w, 0)) {  // NOLINT(readability-implicit-bool-conversion)
      continue;
    }
    if (__builtin_expect(u == destination, 0)) {  // NOLINT(readability-implicit-bool-conversion)
      break;
    }

    for (link* l : adj[u]) {
      int v = l->to;
      if (visited[v]) {
        continue;
      }

      // auto rand_start = std::chrono::steady_clock::now();
      double lat = l->latency(usage);
      int weight = -1;

      if (forbidden_links[l]) {
        weight = static_cast<int>(l->latency(k));
      } else {
        weight = static_cast<int>(
            std::max(min_weight_randdijk,
                     gsl_ran_gaussian_ziggurat(getGSLRng(), stddev_factor_randdijk * lat) + lat));
      }

      // int weight = static_cast<int>(lat);

      /* auto rand_end = std::chrono::steady_clock::now();
      auto rand_dur =
          std::chrono::duration_cast<std::chrono::microseconds>(rand_end - rand_start).count();
      rand_times.push_back(rand_dur); */

      if (dist[v] > dist[u] + weight || dist[v] == -1) {
        a[v] = l;
        dist[v] = dist[u] + weight;
        // std::cout << "dist to " << u << ": " << dist[u] << " weight: " << weight << " new dist["
        // << v << "]: " << dist[v] << std::endl;
        pq.push({dist[v], v});
      }
    }
  }

  // auto route_start = std::chrono::steady_clock::now();
  route r;
  int currentNode = destination;
  while (currentNode != origin) {
    r.links.push_back(a[currentNode]);
    currentNode = a[currentNode]->from;
  }

  // auto reverse_start = std::chrono::steady_clock::now();
  std::reverse(r.links.begin(), r.links.end());
  // auto end = std::chrono::steady_clock::now();

  /*std::cout << "random dijkstra total took "
            << std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() << " us"
            << std::endl;
    std::cout
        << "random dijkstra from " << origin << " to " << destination << " queue size "
        << pq.getSize() << " usage " << usage << " total took "
        << std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() << " us"
        << "(init took "
        << std::chrono::duration_cast<std::chrono::microseconds>(dijk_start - start).count()
        << "us,"
        << "real dijkstra took "
        << std::chrono::duration_cast<std::chrono::microseconds>(route_start - dijk_start).count()
        << "us,"
        << "the randomness took " << std::accumulate(rand_times.begin(), rand_times.end(), 0)
        << "us,"
        << "builing the route took "
        << std::chrono::duration_cast<std::chrono::microseconds>(reverse_start - route_start)
               .count()
        << "us,"
        << "reversing it took "
        << std::chrono::duration_cast<std::chrono::microseconds>(end - reverse_start).count()
        << "us"
        << ")" << std::endl;*/
  r.calculate_params();
  return r;
}

auto determineQueueSize(std::vector<route>& routes, int k) -> int {
  std::vector<int> latencies;
  std::transform(routes.begin(), routes.end(), std::back_inserter(latencies),
                 [k](route r) -> int { return static_cast<int>(r.latency(k)); });

  return std::max(
      static_cast<int>(*max_element(latencies.begin(), latencies.end()) * queue_size_factor),
      minimum_queue_size);
}

auto getBestIndividual(std::vector<island>& islands) -> individual* {
  individual* best = &islands[0].parents[0];

  for (auto& island : islands) {
    if (best->score > island.parents[0].score) {
      best = &island.parents[0];
    }
  }

  return best;
}

void assign_routes(std::vector<route>& routes, const usage& u) {
  std::discrete_distribution<> d(u.routeFlow.begin(), u.routeFlow.end());
  for (auto& p : persons) {
    p.r = std::make_shared<route>(routes[d(getGenerator())]);
  }
}

auto findInVector(const std::vector<int>& vecOfElements, const int& element)
    -> std::pair<bool, int> {
  std::pair<bool, int> result;

  // Find given element in vector
  auto it = std::find(vecOfElements.begin(), vecOfElements.end(), element);

  if (it != vecOfElements.end()) {
    result.second = distance(vecOfElements.begin(), it);
    result.first = true;
  } else {
    result.first = false;
    result.second = -1;
  }

  return result;
}

auto intersection(std::vector<int> v1, std::vector<int> v2) -> std::vector<int> {
  std::vector<int> v3;

  std::sort(v1.begin(), v1.end());
  std::sort(v2.begin(), v2.end());

  std::set_intersection(v1.begin(), v1.end(), v2.begin(), v2.end(), back_inserter(v3));
  return v3;
}

void delete_circle(route& r) {
  std::unordered_map<int, int> last_occ;
  for (size_t i = 0; i < r.links.size(); i++) {
    auto p = last_occ.find(r.links[i]->from);
    if (p == last_occ.end()) {
      last_occ[r.links[i]->from] = i;
    } else {
      r.links.erase(r.links.begin() + p->second, r.links.begin() + i);
      i = 0;
      last_occ.clear();
    }
  }
}

// calculates latency for a route in a set of multiple routes, where edgeFlow defines the usage of
// each edge over all routes
auto route_latency(std::unordered_map<link*, double>& edgeFlow, route& r) -> double {
  double latency = 0;
  for (link* l : r.links) {
    latency += l->latency(edgeFlow[l]);
  }
  return latency;
}

double distScoreRoutes(std::vector<route>& routes) {
  std::unordered_map<link*, int> occurences;
  for (auto& r : routes) {
    for (link* l : r.links) {
      occurences[l]++;
    }
  }
  int individuals = 0;
  int shared = 0;
  for (auto const& [l, occ] : occurences) {
    if (occ == 1) {
      individuals++;
    } else {
      shared += std::pow(occ, 2);
    }
  }

  individuals = (individuals == 0) ? 1 : individuals;
  return static_cast<double>(shared) / static_cast<double>(individuals);
}