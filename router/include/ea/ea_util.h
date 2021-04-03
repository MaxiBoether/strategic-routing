#pragma once
#include <vector>

#include "core/data.h"
#include "ea/ea_data.h"

extern auto randDijkstra(double usage, int origin, int destination,
                         route* notPreferredLinks = nullptr, int k = 1) -> route;

extern int determineQueueSize(std::vector<route>& routes, int k);

extern auto getBestIndividual(std::vector<island>& islands) -> individual*;

extern void assign_routes(std::vector<route>& routes, const usage& u);

extern auto findInVector(const std::vector<int>& vecOfElements, const int& element)
    -> std::pair<bool, int>;

extern auto intersection(std::vector<int> v1, std::vector<int> v2) -> std::vector<int>;

extern void delete_circle(route& r);

extern double route_latency(std::unordered_map<link*, double>& edgeFlow, route& r);

extern double distScoreRoutes(std::vector<route>& routes);