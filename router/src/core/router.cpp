#include <iostream>
#include <vector>
#include <unordered_map>
#include <cstdlib>
#include <cstring>
#include <tinyxml2.h>

#include "core/io.h"
#include "core/data.h"
#include "core/globals.h"
#include "core/routing.h"

std::vector<person> persons = {};
std::vector<node*> nodes = {};
std::vector<std::vector<link*>> adj = {};
PSYCH_MODEL_CLASS psychological_model;

int main(int argc, char *argv[]) {
    if (argc < 4) {
        std::cerr << "<graph> <plans> <output>" << std::endl;
        return 0;
    }
    // load graph
    loadGraph(argv[1]);

    // load od stuff
    loadPlans(argv[2]);

    std::cout << "Done loading Graph and Plans" << std::endl;

    // routing
    char *routing_args[argc-4];
    for (int i=0; i<argc-4;i++) {
        routing_args[i] = argv[i+4];
    }
    do_routing(argc - 4, routing_args);
    std::cout << "Done routing" << std::endl;

    // output
    outputPlans(argv[3]);

    for (node * n : nodes)
        delete n;
    for (auto& n : adj) {
        for (link* l: n)
            delete l;
    }

    return 0;
}
