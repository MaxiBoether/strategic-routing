
#include "core/data.h"

#include <tinyxml2.h>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <iterator>
#include <numeric>
#include <sstream>
#include <string>
#include <vector>

#include "core/globals.h"
#include "core/io.h"
// To pull best routes from the queue first, adjust this value such that it fits your input data.
// Let mean_a and mean_b refer to the mean of the parameters over all edges in your graph. 
// Then set WEIGHT_A to mean_b/mean_a
#define WEIGHT_A 20000
 
using namespace tinyxml2;

person::person(int _origin, int _destination, const char* _timestr, XMLElement* _element)
    : origin(_origin), destination(_destination), timestr(_timestr), element(_element) {}

void person::write() { write(plansXml, element); }

void person::write(tinyxml2::XMLDocument& xml, tinyxml2::XMLElement* _element) {
  std::string routeText = r->to_string();

  _element->DeleteAttribute("origin_node");
  _element->DeleteAttribute("destination_node");
  _element->DeleteAttribute("time");

  XMLElement* planChild = xml.NewElement("plan");
  planChild->SetAttribute("selected", "yes");

  XMLElement* activityOrigin = xml.NewElement("act");
  activityOrigin->SetAttribute("type", "dummy");
  activityOrigin->SetAttribute("x", nodes[origin]->x);
  activityOrigin->SetAttribute("y", nodes[origin]->y);
  activityOrigin->SetAttribute("link", r->links.front()->id);
  activityOrigin->SetAttribute("end_time", timestr.c_str());
  planChild->InsertEndChild(activityOrigin);

  XMLElement* leg = xml.NewElement("leg");
  leg->SetAttribute("mode", "car");
  XMLElement* route = xml.NewElement("route");

  XMLText* routeDesc = xml.NewText(routeText.c_str());
  route->InsertEndChild(routeDesc);
  leg->InsertEndChild(route);
  planChild->InsertEndChild(leg);

  XMLElement* activityDestination = xml.NewElement("act");
  activityDestination->SetAttribute("type", "dummy");
  activityDestination->SetAttribute("x", nodes[destination]->x);
  activityDestination->SetAttribute("y", nodes[destination]->y);
  activityDestination->SetAttribute("link", adj[r->links.back()->to].front()->id);
  planChild->InsertEndChild(activityDestination);

  _element->InsertEndChild(planChild);
}

double link::b() { return _b ? _b : _b = psychological_model.b(*this); }
double link::a() { return _a ? _a : _a = psychological_model.a(*this); }
double link::taud() { return _taud ? _taud : _taud = latency(number_agents); }
double link::latency(double x) {
  return psychological_model.latency(a(), b(), x);  // a()*x*x + b();
}

double route::a() const { return _a; }

double route::b() const { return _b; }

double route::length() {
  double length = 0;
  for (link* l : links)
    length += l->length;
  return length;
}

route route::operator+(const route& b) {
  std::vector<link*> new_links;
  std::copy(links.begin(), links.end(), back_inserter(new_links));
  std::copy(b.links.begin(), b.links.end(), back_inserter(new_links));
  auto r = new route();
  r->_a = _a + b._a;
  r->_b = _b + b._b;
  r->links = new_links;
  return *r;
}

void route::calculate_params() {
  _a = std::accumulate(links.begin(), links.end(), 0.0f,
                       [](double su, link* l) { return su + l->a(); });
  _b = std::accumulate(links.begin(), links.end(), 0.0f,
                       [](double su, link* l) { return su + l->b(); });
}

route::route(route& source, link* l) {
  _a = source.a() + l->a();
  _b = source.b() + l->b();
  std::copy(source.links.begin(), source.links.end(), back_inserter(links));
  links.push_back(l);
}

route::route(std::string nodeString) {
  std::stringstream ss(nodeString);
  std::istream_iterator<std::string> begin(ss);
  std::istream_iterator<std::string> end;
  std::vector<std::string> vstrings(begin, end);
  std::vector<int> nodeVec;
  std::transform(vstrings.begin(), vstrings.end(), std::back_inserter(nodeVec),
                 [](std::string& str) { return std::stoi(str); });

  initByNodeVec(nodeVec);
}

route::route(std::vector<int>& nodeVec) { initByNodeVec(nodeVec); }

void route::initByNodeVec(std::vector<int>& nodeVec) {
  for (unsigned int i = 0; i < nodeVec.size() - 1; i++) {
    int node = nodeVec[i];
    int nextNode = nodeVec[i + 1];
    int found = 0;
    for (link* l : adj[node]) {
      if (l->to == nextNode) {
        if (!found) {
          links.push_back(l);
          found = 1;
        } else {
          std::cout << "Warning, there's another option from " << node << " to " << nextNode
                    << " while constructing route." << std::endl;
          links[links.size() - 1] = l;
        }
      }
    }
    if (!found) {
      std::cout << "WARNING!" << std::endl;
      std::cout << "Could not find outgoing edge from " << node << " to " << nextNode
                << " while constructing route." << std::endl;
    }
  }

  calculate_params();
}

route::route(std::vector<link*>& linkVec) {
  std::copy(linkVec.begin(), linkVec.end(), back_inserter(links));
  calculate_params();
}

route::route(double a, double b) {
  _a = a;
  _b = b;
}

std::string route::to_string() {
  return std::accumulate(
      links.begin(), links.end(), std::to_string(links.front()->from),
      [](std::string a, link* l) { return std::move(a) + " " + std::to_string(l->to); });
}

std::string route::to_linkid_string() {
  return std::accumulate(
      std::next(links.begin()), links.end(), std::to_string(links.front()->id),
      [](std::string a, link* l) { return std::move(a) + " " + std::to_string(l->id); });
}

std::vector<int> route::to_node_vec() {
  std::vector<int> _nodes;
  _nodes.push_back(links.front()->from);
  for (link* l : links) {
    _nodes.push_back(l->to);
  }

  return _nodes;
}

double route::compareTo(route& other) {
  std::vector<link*> r1(links);
  std::vector<link*> r2(other.links);
#pragma omp parallel sections
  {
#pragma omp section
    std::sort(r1.begin(), r1.end(),
              [](const link* a, const link* b) -> bool { return a->id > b->id; });
#pragma omp section
    std::sort(r2.begin(), r2.end(),
              [](const link* a, const link* b) -> bool { return a->id > b->id; });
  }

  std::vector<link*> symmetric_difference;
  std::set_symmetric_difference(r1.begin(), r1.end(), r2.begin(), r2.end(),
                                back_inserter(symmetric_difference),
                                [](const link* a, const link* b) -> bool { return a->id > b->id; });

  double diff =
      static_cast<double>(symmetric_difference.size()) / static_cast<double>(r1.size() + r2.size());

  return diff;
}

double ParetoElement::k() const {
  return a() * WEIGHT_A + b();
}
ParetoElement::ParetoElement(shared_ptr<ParetoElement> par, link* l) {
  parent = par;
  myLink = l;
  if (l->to == 0 && l->from == 0) {
    _a = parent->a();
    _b = parent->b();
    _taud = parent->taud();
  } else {
    _a = parent->a() + l->a();
    _b = parent->b() + l->b();
    _taud = parent->taud() + l->taud();
  }
}

ParetoElement::ParetoElement(shared_ptr<ParetoElement> par, link* l, bool shared) {
  parent = par;
  myLink = l;
  _a = parent->a();
  _b = parent->b();
  _taud = parent->taud();
  _shared_a = parent->shared_a();
  _shared_b = parent->shared_b();
  _shared_taud = parent->shared_taud();
  if (!(l->to == 0 && l->from == 0)) {
    _a += l->a();
    _b += l->b();
    _taud += l->taud();
    if (shared) {
      _shared_a += l->a();
      _shared_b += l->b();
      _shared_taud += l->taud();
    }
  }
}

ParetoElement::ParetoElement(shared_ptr<ParetoElement> par) {  // recursive memory copy
  if (par->parent != nullptr)
    parent = make_shared<ParetoElement>(par->parent);
  myLink = par->myLink;
  _a = par->a();
  _b = par->b();
  _taud = par->shared_taud();
  _shared_a = par->shared_a();
  _shared_b = par->shared_b();
  _shared_taud = par->shared_taud();
}

ParetoElement::ParetoElement(double a, double b, double taud, double sa, double sb, double staud) {
  _a = a;
  _b = b;
  _taud = taud;
  _shared_a = sa;
  _shared_b = sb;
  _shared_taud = staud;
}
double ParetoElement::a() const { return _a; }
double ParetoElement::b() const { return _b; }
double ParetoElement::taud() const { return _taud; }
double ParetoElement::shared_a() const { return _shared_a; }
double ParetoElement::shared_b() const { return _shared_b; }
double ParetoElement::shared_taud() const { return _shared_taud; }
bool ParetoElement::operator<(const ParetoElement& other) const { return k() < other.k(); }

unique_ptr<vector<link*>> ParetoElement::collectLinks() {
  if (parent != nullptr) {
    auto links = parent->collectLinks();
    links->push_back(myLink);
    return links;
  }
  return make_unique<vector<link*>>();
}
shared_ptr<route> ParetoElement::collectRoute() {
  auto links = collectLinks();
  route* r = new route(*links);
  return shared_ptr<route>(r);
}

