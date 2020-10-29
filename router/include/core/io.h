#pragma once

#include <tinyxml2.h>

void readPersons(const tinyxml2::XMLDocument&);
void loadGraph(const char* graphFile);
void loadPlans(const char* plansFile);
void outputPlans(const char* plansFile);
void outputPlansToNewFile(const char* plansFile);

extern tinyxml2::XMLDocument plansXml;
