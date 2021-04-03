#include "ea/ea_mutations.h"

#include <algorithm>
#include <chrono>
#include <fstream>
#include <iostream>
#include <random>
#include <set>
#include <sstream>
#include <string>
#include <vector>

#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>

#include "ea/ea_constants.h"
#include "ea/ea_data.h"
#include "ea/ea_globals.h"
#include "ea/ea_logging.h"
#include "ea/ea_random.h"
#include "ea/ea_util.h"

#include "json.hpp"

using json = nlohmann::json;

bool debugEnabled = false;

auto mutation_capacity(individual& indiv, int k, int overall_origin, int overall_destination)
    -> json {
  (void)overall_origin;
  (void)overall_destination;

  std::vector<double> weights(indiv.routes.size(), 1);
  for (unsigned int i = 0; i < indiv.routes.size(); i++) {
    weights[i] = 1.0F / static_cast<double>(indiv.u.routeFlow[i]);
  }
  std::discrete_distribution<> d(weights.begin(), weights.end());
  // now, we have a discrete distribution over all routes depending on the flow on the routes

  json result = {{"modifications", json::array()}};
  int route_count = std::min(static_cast<int>(indiv.routes.size()),
                             std::max(1, static_cast<int>(gsl_ran_poisson(getGSLRng(), 1))));

  std::set<int> routes;
  for (int i = 0; i < route_count; i++) {
    routes.insert(d(getGenerator()));
  }

  for (unsigned int j = 0; j < routes.size(); j++) {
    int i = *std::next(routes.begin(), j);
    auto start = std::chrono::steady_clock::now();
    route& r = indiv.routes[i];

    std::vector<double> link_weights(r.links.size(), 0);
    std::vector<double> link_weights2(r.links.size(), 0);
    for (unsigned int currLink = 0; currLink < r.links.size(); currLink++) {
      link_weights[currLink] += std::pow(r.links[currLink]->a(), 3);
      for (link* l : adj[r.links[currLink]->from]) {
        if (l != r.links[currLink]) {
          link_weights2[currLink] += 1.0F / l->a();
        }
      }
    }

    std::discrete_distribution<> dis(link_weights.begin(), link_weights.end() - 1);

    int start_index =
        std::max(0, dis(getGenerator()) + static_cast<int>(gsl_ran_gaussian_ziggurat(
                                              getGSLRng(), stddev_mutation_capacity_start)));
    start_index = std::min(start_index, static_cast<int>(r.links.size()) - 4);

    std::discrete_distribution<> end_dis(link_weights2.begin() + start_index + 1,
                                         link_weights2.end());

    int origin = r.links[start_index]->from;

    int advance = std::max(
        end_dis(getGenerator()) + static_cast<int>(gsl_ran_gaussian_ziggurat(
                                      getGSLRng(), stddev_mutation_capacity_advance)),
        4);  // this is a bit fuzzy with end_dis but as we add noise anyways we don't really care
    int end_index = std::min(advance + start_index, static_cast<int>(r.links.size()) - 1);

    // std::cout << "start_index " << start_index << " end_index " << end_index << std::endl;
    if (end_index >= static_cast<int>(r.links.size()) || end_index <= start_index) {
      std::cout << "HOLY SHIT HOUSE UNDSO" << std::endl;
      std::cout << "end: " << end_index << " start: " << start_index
                << " size: " << static_cast<int>(r.links.size()) << std::endl;
      continue;
    }

    int destination = r.links[end_index]->from;
    end_index--;  // decrease end_index because if we choose edge e from s' to t', s' is the node
                  // we mean, thus we have to replace until the edge before that

    auto selection_time = std::chrono::steady_clock::now();
    route new_r = randDijkstra(indiv.u.routeFlow[i], origin, destination, &r, k);
    auto rand_dijk = std::chrono::steady_clock::now();
    if (new_r.links.size() > 1) {
      r.links.erase(r.links.begin() + start_index, r.links.begin() + end_index + 1);
      r.links.insert(r.links.begin() + start_index, new_r.links.begin(), new_r.links.end());
      auto insertion = std::chrono::steady_clock::now();
      delete_circle(r);
      auto end = std::chrono::steady_clock::now();
      r.calculate_params();
      json under_result = {
          {"route", i},
          {"origin", origin},
          {"destination", destination},
          {"times",
           {{"total", std::chrono::duration_cast<std::chrono::microseconds>(end - start).count()},
            {"init",
             std::chrono::duration_cast<std::chrono::microseconds>(selection_time - start).count()},
            {"random_dijkstra",
             std::chrono::duration_cast<std::chrono::microseconds>(rand_dijk - selection_time)
                 .count()},
            {"update",
             std::chrono::duration_cast<std::chrono::microseconds>(insertion - rand_dijk).count()},
            {"circle_deletion",
             std::chrono::duration_cast<std::chrono::microseconds>(end - insertion).count()}}}};

      result["modifications"].push_back(under_result);
    } else {
      /*std::cout << "holy shit the house is on fire in capacity mutation" << std::endl;
      std::cout << "end: " << end_index << " start: " << start_index
                << " size: " << static_cast<int>(r.links.size()) << " advance: " << advance <<
      std::endl; */
    }
  }

  return result;
}

auto mutation_weightpair(individual& indiv, int k, int overall_origin, int overall_destination)
    -> json {
  (void)overall_origin;
  (void)overall_destination;

  std::vector<double> weights(indiv.routes.size(), 1);
  for (unsigned int i = 0; i < indiv.routes.size(); i++) {
    weights[i] = 1.0F / static_cast<double>(indiv.u.routeFlow[i]);
  }
  std::discrete_distribution<> d(weights.begin(), weights.end());
  // now, we have a discrete distribution over all routes depending on the flow on the routes

  json result = {{"modifications", json::array()}};

  int route_count = std::min(static_cast<int>(indiv.routes.size()),
                             std::max(1, static_cast<int>(gsl_ran_poisson(getGSLRng(), 1))));

  std::set<int> routes;
  for (int i = 0; i < route_count; i++) {
    routes.insert(d(getGenerator()));
  }

  for (unsigned int j = 0; j < routes.size(); j++) {
    int i = *std::next(routes.begin(), j);
    auto start = std::chrono::steady_clock::now();
    route& r = indiv.routes[i];

    std::vector<double> link_weights(r.links.size(), 0);
    for (unsigned int currLink = 0; currLink < r.links.size(); currLink++) {
      for (link* l : adj[r.links[currLink]->from]) {
        if (l != r.links[currLink]) {
          link_weights[currLink] += l->capacity;
          // link_weights[i] += 1;
        }
      }
    }

    std::discrete_distribution<> dis(
        link_weights.begin(),
        link_weights.end() - 2);  // -1 ist neu - bisher konnte es vorkommen, dass bei einer size
                                  // von 144 der startindex 143 war - wieso!?
    int start_index = dis(getGenerator());

    std::discrete_distribution<> end_dis(link_weights.begin() + start_index + 2,
                                         link_weights.end());

    int origin = r.links[start_index]->from;

    int end_index = end_dis(getGenerator()) + start_index + 2;

    if (end_index >= static_cast<int>(r.links.size()) || end_index <= start_index) {
      std::cout << "HOLY SHIT HOUSE UNDSO" << std::endl;
      std::cout << "end: " << end_index << " start: " << start_index
                << " size: " << static_cast<int>(r.links.size()) << std::endl;
      continue;
    }

    int destination = r.links[end_index]->from;
    end_index--;  // decrease end_index because if we choose edge e from s' to t', s' is the node
                  // we mean, thus we have to replace until the edge before that

    auto selection_time = std::chrono::steady_clock::now();
    route new_r = randDijkstra(indiv.u.routeFlow[i], origin, destination, &r, k);
    auto rand_dijk = std::chrono::steady_clock::now();
    if (new_r.links.size() > 1) {
      r.links.erase(r.links.begin() + start_index, r.links.begin() + end_index + 1);
      r.links.insert(r.links.begin() + start_index, new_r.links.begin(), new_r.links.end());
      auto insertion = std::chrono::steady_clock::now();
      delete_circle(r);
      auto end = std::chrono::steady_clock::now();
      r.calculate_params();
      json under_result = {
          {"route", i},
          {"origin", origin},
          {"destination", destination},
          {"times",
           {{"total", std::chrono::duration_cast<std::chrono::microseconds>(end - start).count()},
            {"init",
             std::chrono::duration_cast<std::chrono::microseconds>(selection_time - start).count()},
            {"random_dijkstra",
             std::chrono::duration_cast<std::chrono::microseconds>(rand_dijk - selection_time)
                 .count()},
            {"update",
             std::chrono::duration_cast<std::chrono::microseconds>(insertion - rand_dijk).count()},
            {"circle_deletion",
             std::chrono::duration_cast<std::chrono::microseconds>(end - insertion).count()}}}};

      result["modifications"].push_back(under_result);
    } else {
      /*std::cout << "holy shit the house is on fire in weightpair mutation" << std::endl; */
    }
  }

  return result;
}

auto mutation_randpair(individual& indiv, int k, int overall_origin, int overall_destination)
    -> json {
  (void)overall_origin;
  (void)overall_destination;

  std::vector<double> weights(indiv.routes.size(), 1);
  for (unsigned int i = 0; i < indiv.routes.size(); i++) {
    weights[i] = 1.0F / static_cast<double>(indiv.u.routeFlow[i]);
  }
  std::discrete_distribution<> d(weights.begin(), weights.end());
  // now, we have a discrete distribution over all routes depending on the flow on the routes

  json result = {{"modifications", json::array()}};

  int route_count = std::min(static_cast<int>(indiv.routes.size()),
                             std::max(1, static_cast<int>(gsl_ran_poisson(getGSLRng(), 1))));

  std::set<int> routes;
  for (int i = 0; i < route_count; i++) {
    routes.insert(d(getGenerator()));
  }

  for (unsigned int j = 0; j < routes.size(); j++) {
    int i = *std::next(routes.begin(), j);
    auto start = std::chrono::steady_clock::now();
    route& r = indiv.routes[i];

    std::uniform_int_distribution<> dis(0, r.links.size() - 2);
    std::normal_distribution<> lendis{average_factor_mutate * r.links.size(),
                                      stddev_factor_mutate * r.links.size()};
    int start_index = std::min(dis(getGenerator()), static_cast<int>(r.links.size()) - 3);
    int origin = r.links[start_index]->from;
    int offset = std::max(2, static_cast<int>(std::round(lendis(getGenerator()))));
    int end_index = std::min(start_index + offset, static_cast<int>(r.links.size()) - 1);
    if (end_index <= start_index) {
      std::cout << "HOLY SHIT THE HOUSE IS ON FIRE IN RANDPAIR MUTATION." << std::endl;
      continue;
    }
    int destination = r.links[end_index]->to;
    auto selection_time = std::chrono::steady_clock::now();
    route new_r = randDijkstra(indiv.u.routeFlow[i], origin, destination, &r, k);
    auto rand_dijk = std::chrono::steady_clock::now();
    if (new_r.links.size() > 1) {
      r.links.erase(r.links.begin() + start_index, r.links.begin() + end_index + 1);
      r.links.insert(r.links.begin() + start_index, new_r.links.begin(), new_r.links.end());
      auto insertion = std::chrono::steady_clock::now();
      delete_circle(r);
      auto end = std::chrono::steady_clock::now();
      r.calculate_params();
      json under_result = {
          {"route", i},
          {"origin", origin},
          {"destination", destination},
          {"times",
           {{"total", std::chrono::duration_cast<std::chrono::microseconds>(end - start).count()},
            {"init",
             std::chrono::duration_cast<std::chrono::microseconds>(selection_time - start).count()},
            {"random_dijkstra",
             std::chrono::duration_cast<std::chrono::microseconds>(rand_dijk - selection_time)
                 .count()},
            {"update",
             std::chrono::duration_cast<std::chrono::microseconds>(insertion - rand_dijk).count()},
            {"circle_deletion",
             std::chrono::duration_cast<std::chrono::microseconds>(end - insertion).count()}}}};

      result["modifications"].push_back(under_result);
    } else {
      /*std::cout << "holy shit the house is on fire in randpair mutation" << std::endl;
      std::cout << "end: " << end_index << " start: " << start_index
                << " size: " << static_cast<int>(r.links.size()) << std::endl; */
    }
  }

  return result;
}

auto mutation_exchangepart(individual& indiv, int k, int origin, int destination) -> json {
  (void)k;
  int overall_orig = origin;
  int overall_dest = destination;
  (void)origin;
  (void)destination;

  auto start = std::chrono::steady_clock::now();
  int route1id = -1;
  int route2id = -1;
  std::vector<int> divergence_points;
  int dp = -1;
  int goto_node = -1;
  // 1. select two random routes we want to crossover
  if (indiv.routes.size() > 1) {
    std::uniform_int_distribution<> dis(0, indiv.routes.size() - 1);
    route1id = dis(getGenerator());
    route2id = dis(getGenerator());

    while (route1id == route2id) {
      route2id = dis(getGenerator());
    }

    // 2. find divergence points of these routes
    // divergence points are points that are shared, but then go on to different vertices
    delete_circle(indiv.routes[route1id]);
    delete_circle(indiv.routes[route2id]);

    std::vector<int> route1 = indiv.routes[route1id].to_node_vec();
    std::vector<int> route2 = indiv.routes[route2id].to_node_vec();

    std::vector<int> shared_points = intersection(route1, route2);
    if (debugEnabled) {
      /*std::cout << "route1: " << std::endl;

      for (int p : route1) {
        std::cout << p << " ";
      }
      std::cout << std::endl;

      std::cout << "route2: " << std::endl;

      for (int p : route2) {
        std::cout << p << " ";
      }
      std::cout << std::endl; */

      std::cout << "shared points: " << std::endl;

      for (int p : shared_points) {
        std::cout << p << " ";
      }
      std::cout << std::endl;
    }

    for (int node : shared_points) {
      int idx1 = findInVector(route1, node).second;
      int idx2 = findInVector(route2, node).second;

      if (idx1 + 1 < static_cast<int>(route1.size()) &&
          idx2 + 1 < static_cast<int>(route2.size())) {
        if (route1[idx1 + 1] != route2[idx2 + 1]) {
          divergence_points.push_back(node);
        }
      }
    }

    if (debugEnabled) {
      std::cout << "divergence points: " << std::endl;

      for (int p : divergence_points) {
        std::cout << p << " ";
      }
      std::cout << std::endl;
    }
    if (!divergence_points.empty()) {
      std::uniform_int_distribution<> dis2(0, divergence_points.size() - 1);
      int dp_index = dis2(getGenerator());
      dp = divergence_points[dp_index];
      int dp_route_index = findInVector(route1, dp).second;
      int dp_route_index2 = findInVector(route2, dp).second;

      std::vector<int> gotoPoints;

      for (int node : shared_points) {
        if (node != dp) {
          if (findInVector(route1, node).second > dp_route_index && findInVector(route2, node).second > dp_route_index2) {
            gotoPoints.push_back(node);
          }
        }
      }
      if (!gotoPoints.empty()) {
        std::uniform_int_distribution<> dis3(0, gotoPoints.size() - 1);

        int goto_index = dis3(getGenerator());
        goto_node = gotoPoints[goto_index];

        if (debugEnabled) {
          std::cout << "exchanging from " << dp << " to " << goto_node << std::endl;
          std::cout << "goto points: " << std::endl;

          for (int p : gotoPoints) {
            std::cout << p << " ";
          }
          std::cout << std::endl;
        }
        route& r1 = indiv.routes[route1id];
        route& r2 = indiv.routes[route2id];

        int l1 = -1;
        int l2 = -1;
        int l3 = -1;
        int l4 = -1;
        // find link that has FROM dp TO something in route1
        for (int i = 0; i < static_cast<int>(r1.links.size()); i++) {
          link* l = r1.links[i];
          if (l->from == dp) {
            l1 = i;
            break;
          }
        }

        // find link that has FROM dp TO something in route2
        for (int i = 0; i < static_cast<int>(r2.links.size()); i++) {
          link* l = r2.links[i];
          if (l->from == dp) {
            l2 = i;
            break;
          }
        }

        // find link that has FROM something' TO goto_node in route1
        for (int i = 0; i < static_cast<int>(r1.links.size()); i++) {
          link* l = r1.links[i];
          if (l->to == goto_node) {
            l3 = i;
            break;
          }
        }

        // find link that has FROM something' TO goto_node in route2
        for (int i = 0; i < static_cast<int>(r2.links.size()); i++) {
          link* l = r2.links[i];
          if (l->to == goto_node) {
            l4 = i;
            break;
          }
        }

        if (l1 < 0 || l2 < 0 || l3 < 0 || l4 < 0) {
          std::uniform_int_distribution<> dis4(0, 5000);
          std::ofstream outfile("error_rand" + std::to_string(dis4(getGenerator())) + ".log");
          outfile << "Error in mutation_exchangepart. " << std::endl
                  << "Overall Origin: " << overall_orig << " and Destination: " << overall_dest
                  << std::endl
                  << "Route1ID: " << route1id << " Route2ID: " << route2id << std::endl;

          std::cout << "HOLY SH*T THE HOUSE IS ON FIRE!" << std::endl;
          std::cout << "couldn't find links for exchange" << std::endl;
          std::cout << "exchanging from " << dp << " to " << goto_node << std::endl;
          std::cout << "l1: " << l1 << " l2: " << l2 << " l3: " << l3 << " l4: " << l4 << std::endl;
          std::cout << "route1: " << std::endl;

          outfile << "HOLY SH*T THE HOUSE IS ON FIRE!" << std::endl;
          outfile << "couldn't find links for exchange" << std::endl;
          outfile << "exchanging from " << dp << " to " << goto_node << std::endl;
          outfile << "l1: " << l1 << " l2: " << l2 << " l3: " << l3 << " l4: " << l4 << std::endl;
          outfile << "route1: " << std::endl;

          for (int p : route1) {
            std::cout << p << " ";
            outfile << p << " ";
          }
          std::cout << std::endl;
          outfile << std::endl;

          std::cout << "route2: " << std::endl;
          outfile << "route2: " << std::endl;

          for (int p : route2) {
            std::cout << p << " ";
            outfile << p << " ";
          }
          std::cout << std::endl;
          outfile << std::endl;

          std::cout << "shared points: " << std::endl;
          outfile << "shared points: " << std::endl;

          for (int p : shared_points) {
            std::cout << p << " ";
            outfile << p << " ";
          }
          std::cout << std::endl;
          outfile << std::endl;

          std::cout << "divergence points: " << std::endl;
          outfile << "divergence points: " << std::endl;

          for (int p : divergence_points) {
            std::cout << p << " ";
            outfile << p << " ";
          }

          std::cout << std::endl;
          outfile << std::endl;

          std::cout << "goto points: " << std::endl;
          outfile << "goto points: " << std::endl;

          for (int p : gotoPoints) {
            outfile << p << " ";

            std::cout << p << " ";
          }
          std::cout << std::endl;
          outfile << std::endl;
          outfile << std::endl;
          outfile << std::endl;
          outfile << "route1 with tostring method (nodes): " << std::endl;
          outfile << indiv.routes[route1id].to_string() << std::endl;
          outfile << "route1 with to linkid string method (links): " << std::endl;
          outfile << indiv.routes[route1id].to_linkid_string() << std::endl << std::endl;
          outfile << "route2 with tostring method (nodes): " << std::endl;
          outfile << indiv.routes[route2id].to_string() << std::endl;
          outfile << "route2 with to linkid string method (links): " << std::endl;
          outfile << indiv.routes[route2id].to_linkid_string() << std::endl
                  << std::endl
                  << std::endl;

          json result = {{"exchange_fail", 1},
                         {"route1id", route1id},
                         {"route2id", route2id},
                         {"div_point", dp},
                         {"goto_node", goto_node},
                         {"divergence_points", divergence_points},
                         {"route1", route1},
                         {"route2", route2},
                         {"shared_points", shared_points},
                         {"goto_points", gotoPoints}};

          outfile.close();
          return result;
        }

        std::vector<link*> _route1part = {r1.links.begin() + l1, r1.links.begin() + l3 + 1};
        // std::cout << "found link was " << indiv.routes[route1id].links[l3]->id << ", last is " <<
        // _route1part[_route1part.size() - 1]->id;

        /*for (link* l : _route1part) {
          std::cout << l->id << " ";
        }
        std::cout << std::endl;
        for (int i = l2; i < l4 + 1; i++) {
          std::cout << r2.links[i]->id << " ";
        }
        std::cout << std::endl; */

        if (r1.links[l3]->id != _route1part[_route1part.size() - 1]->id) {
          std::cout << "holy shit, again fire!" << std::endl;
          std::cout << "ids sind nicht gleich" << std::endl;
        }
        // std::cout << "removing " << indiv.routes[route1id].links[l1]->id << " to " <<
        // indiv.routes[route1id].links[l3]->id << std::endl;

        // std::cout << indiv.routes[route1id].to_id_string() << std::endl;
        // std::cout << std::endl;
        // std::cout << indiv.routes[route2id].to_id_string() << std::endl;

        r1.links.erase(r1.links.begin() + l1, r1.links.begin() + l3 + 1);

        // r.links.insert(r.links.begin() + start_index, new_r.links.begin(), new_r.links.end());
        try {
          r1.links.insert(r1.links.begin() + l1, r2.links.begin() + l2, r2.links.begin() + l4 + 1);
        } catch (std::exception const& e) {
          std::uniform_int_distribution<> dis4(0, 5000);
          std::ofstream outfile("error_rand" + std::to_string(dis4(getGenerator())) + ".log");
          outfile << "Error in mutation_exchangepart That we were looking for!!. " << std::endl
                  << "Overall Origin: " << overall_orig << " and Destination: " << overall_dest
                  << std::endl
                  << "Route1ID: " << route1id << " Route2ID: " << route2id << std::endl;

          std::cout << "HOLY SH*T THE HOUSE IS ON FIRE!" << std::endl;
          std::cout << "couldn't find links for exchange" << std::endl;
          std::cout << "exchanging from " << dp << " to " << goto_node << std::endl;
          std::cout << "l1: " << l1 << " l2: " << l2 << " l3: " << l3 << " l4: " << l4 << std::endl;
          std::cout << "route1: " << std::endl;

          outfile << "HOLY SH*T THE HOUSE IS ON FIRE!" << std::endl;
          outfile << "couldn't find links for exchange" << std::endl;
          outfile << "exchanging from " << dp << " to " << goto_node << std::endl;
          outfile << "l1: " << l1 << " l2: " << l2 << " l3: " << l3 << " l4: " << l4 << std::endl;
          outfile << "route1: " << std::endl;

          for (int p : route1) {
            std::cout << p << " ";
            outfile << p << " ";
          }
          std::cout << std::endl;
          outfile << std::endl;

          std::cout << "route2: " << std::endl;
          outfile << "route2: " << std::endl;

          for (int p : route2) {
            std::cout << p << " ";
            outfile << p << " ";
          }
          std::cout << std::endl;
          outfile << std::endl;

          std::cout << "shared points: " << std::endl;
          outfile << "shared points: " << std::endl;

          for (int p : shared_points) {
            std::cout << p << " ";
            outfile << p << " ";
          }
          std::cout << std::endl;
          outfile << std::endl;

          std::cout << "divergence points: " << std::endl;
          outfile << "divergence points: " << std::endl;

          for (int p : divergence_points) {
            std::cout << p << " ";
            outfile << p << " ";
          }

          std::cout << std::endl;
          outfile << std::endl;

          std::cout << "goto points: " << std::endl;
          outfile << "goto points: " << std::endl;

          for (int p : gotoPoints) {
            outfile << p << " ";

            std::cout << p << " ";
          }
          std::cout << std::endl;
          outfile << std::endl;
          outfile << std::endl;
          outfile << std::endl;
          outfile << "route1 with tostring method (nodes): " << std::endl;
          outfile << indiv.routes[route1id].to_string() << std::endl;
          outfile << "route1 with to linkid string method (links): " << std::endl;
          outfile << indiv.routes[route1id].to_linkid_string() << std::endl << std::endl;
          outfile << "route2 with tostring method (nodes): " << std::endl;
          outfile << indiv.routes[route2id].to_string() << std::endl;
          outfile << "route2 with to linkid string method (links): " << std::endl;
          outfile << indiv.routes[route2id].to_linkid_string() << std::endl
                  << std::endl
                  << std::endl;

          delete_circle(indiv.routes[route1id]);
          delete_circle(indiv.routes[route2id]);
          outfile << "after deleting circles:" << std::endl;
          outfile << "route1 with tostring method (nodes): " << std::endl;
          outfile << indiv.routes[route1id].to_string() << std::endl;
          outfile << "route1 with to linkid string method (links): " << std::endl;
          outfile << indiv.routes[route1id].to_linkid_string() << std::endl << std::endl;
          outfile << "route2 with tostring method (nodes): " << std::endl;
          outfile << indiv.routes[route2id].to_string() << std::endl;
          outfile << "route2 with to linkid string method (links): " << std::endl;
          outfile << indiv.routes[route2id].to_linkid_string() << std::endl
                  << std::endl
                  << std::endl;

          json result = {{"exchange_fail", 1},
                         {"route1id", route1id},
                         {"route2id", route2id},
                         {"div_point", dp},
                         {"goto_node", goto_node},
                         {"divergence_points", divergence_points},
                         {"route1", route1},
                         {"route2", route2},
                         {"shared_points", shared_points},
                         {"goto_points", gotoPoints}};

          outfile.close();
          throw;
        }
        // std::cout << "after" << std::endl;
        // std::cout << indiv.routes[route1id].to_id_string() << std::endl;
        r2.links.erase(r2.links.begin() + l2, r2.links.begin() + l4 + 1);
        r2.links.insert(r2.links.begin() + l2, _route1part.begin(), _route1part.end());

        r1.calculate_params();
        r2.calculate_params();

        // std::cout << "done" << std::endl;
        // std::cout << indiv.routes[route1id].to_id_string() << std::endl;
        // std::cout << std::endl;
        // std::cout << indiv.routes[route2id].to_id_string() << std::endl;
      }
    }
  }
  auto end = std::chrono::steady_clock::now();

  json result = {
      {"route1id", route1id},
      {"route2id", route2id},
      {"div_point", dp},
      {"goto_node", goto_node},
      {"divergence_points", divergence_points},
      {"times",
       {{"total", std::chrono::duration_cast<std::chrono::microseconds>(end - start).count()}}}};
  return result;
}

auto mutation_newroute(individual& indiv, int k, int origin, int destination) -> json {
  auto start = std::chrono::steady_clock::now();
  std::vector<double> weights(indiv.routes.size(), 1);
  for (unsigned int i = 0; i < indiv.routes.size(); i++) {
    weights[i] = 1.0F / static_cast<double>(indiv.u.routeFlow[i]);
  }
  std::discrete_distribution<> d(weights.begin(), weights.end());
  int routeToBeRemoved = d(getGenerator());

  auto rand_dijk = std::chrono::steady_clock::now();
  indiv.routes[routeToBeRemoved] =
      randDijkstra(static_cast<float>(k) / static_cast<float>(indiv.routes.size()), origin,
                   destination, &indiv.routes[routeToBeRemoved], k);

  indiv.routes[routeToBeRemoved].calculate_params();
  auto end = std::chrono::steady_clock::now();

  json result = {
      {"replaced_route", routeToBeRemoved},
      {"times",
       {{"total", std::chrono::duration_cast<std::chrono::microseconds>(end - start).count()},
        {"init", std::chrono::duration_cast<std::chrono::microseconds>(rand_dijk - start).count()},
        {"random_dijkstra",
         std::chrono::duration_cast<std::chrono::microseconds>(end - rand_dijk).count()}}}};

  return result;
}

std::vector<json (*)(individual&, int, int, int)> mutation_operators{
    &mutation_randpair, &mutation_weightpair, &mutation_newroute, &mutation_capacity,
    &mutation_exchangepart};

std::vector<std::string> mutation_operator_names{"randpair", "weightpair", "newroute", "capacity",
                                                 "exchange"};

std::vector<double> update_probs(int iteration, std::vector<int>& lastExchanges, int indiv_id,
                                 int iterationsWithoutChange) {
  std::vector<double> result{prob_mutate_randpair, prob_mutate_weightpair, prob_mutate_newroute,
                             prob_mutate_capacity, prob_mutate_exchange};
  if (iteration > iterations_constant_newroute_prob) {
    if (prob_mutate_newroute > 0) {
      double newroute_weight =
          std::max(((1 - prob_mutate_newroute) / iterations_after_newrout_prob_one) *
                          static_cast<double>(iteration) +
                      prob_mutate_newroute,
                  static_cast<double>(1));
      result[2] = newroute_weight;
    } else {
      result[2] = 0;
    }
  }

  if (iteration - lastExchanges[indiv_id] > gap_between_exchanges && prob_mutate_exchange > 0) {
    result[4] = static_cast<double>(std::min(
        max_prob_mutate_exchange,
        static_cast<int>((static_cast<float>(max_prob_mutate_exchange - prob_mutate_exchange) /
                          (max_prob_mutate_reached_after *
                           static_cast<float>(iterations_needed_for_convergeance))) *
                             static_cast<float>(iterationsWithoutChange) +
                         prob_mutate_exchange)));
  } else {
    result[4] = 0;
  }

  return result;
}

void check_route_sanity(route& r, int op) {
  int last_node = r.links[0]->from;
  for (link* l : r.links) {
    if (l->from != last_node) {
      std::cout << "invalid route after operator " << op << ". HIER KOMMT DIE ROUTENPOLIZEI! "
                << l->from << " " << last_node << std::endl;
      exit(1337);
    }
    last_node = l->to;
  }
}

void mutate(individual& indiv, int k, int origin, int destination, int indiv_id,
            json& iteration_json, int iteration, int iterationsWithoutChange,
            std::vector<int>& lastExchanges, int island_id) {
  int mutations = std::max(1, static_cast<int>(gsl_ran_poisson(getGSLRng(), average_mutations)));

  std::vector<double> local_probs =
      update_probs(iteration, lastExchanges, indiv_id, iterationsWithoutChange);

  std::discrete_distribution<> probd(local_probs.begin(), local_probs.end());

  json mutations_json = {{"island", island_id},
                         {"individual", indiv_id},
                         {"mutation_count", mutations},
                         {"probs", local_probs}};

  std::vector<int> mutation_vector;
  for (int i = 0; i < mutations; i++) {
    int op = probd(getGenerator());
    if (op != 4) {
      mutation_vector.push_back(op);
    } else {
      mutation_vector.clear();
      mutation_vector.push_back(op);
      lastExchanges[indiv_id] = iteration;
      break;
    }
  }

  for (int op : mutation_vector) {
    json result = {{"type", mutation_operator_names[op]}};
    result.update(mutation_operators[op](indiv, k, origin, destination));
    mutations_json["details"].push_back(result);

    for (auto& r : indiv.routes) {
      check_route_sanity(r, op);
    }
  }

  inform_log_about_mutations_json(iteration_json, mutations_json);
}