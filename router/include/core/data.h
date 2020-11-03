#pragma once

#include <vector>
#include <string>
#include <memory>
#include <tinyxml2.h>

struct node {
    char *x, *y;
    ~node()  {
        free(x);
        free(y);
    }
};

class link {
    public:
        double _a = 0, _b = 0, _taud = 0;
        int id;
        int from, to;
        double length, capacity, freespeed;
        double b();
        double a();
        double taud(); //The time spend on the link per agent if all d agents use the link.
        double latency(double x); //The time spend on the link per agent if there are x * d agents on the link.
};

class route {
private:
    double _a, _b;
    void initByNodeVec(std::vector<int>& nodeVec);
public:
    std::vector<link*> links; // TODO prevent unauthorized modification
    double a() const;
    double b() const;
    route() = default;
    double length();
    void calculate_params();
    route operator+(const route& b);
    explicit route(route& source, link* l);
    explicit route(std::vector<link*>& linkVec);
    explicit route(std::vector<int>& nodeVec);
    explicit route(std::string nodeString);

    explicit route(double a, double b);
    std::string to_string();
    std::string to_linkid_string();
    std::vector<int> to_node_vec();
    double compareTo(route& other);
};

struct person {
    int origin, destination;
    std::string timestr;
    std::shared_ptr<route> r = nullptr;
    tinyxml2::XMLElement *element;
    person(int _origin, int _destination, const char* _timestr, tinyxml2::XMLElement* _element);
    void write(tinyxml2::XMLDocument& xml, tinyxml2::XMLElement* _element);
    void write();
};

class ParetoElement {
 private:
  double _a = 0.0, _b = 0.0, _shared_a = 0.0, _shared_b = 0.0, _taud = 0.0, _shared_taud = 0.0;
  std::shared_ptr<ParetoElement> parent = nullptr;

 public:
  bool hasSplit = false;
  link* myLink = nullptr;
  ParetoElement() = default;
  ParetoElement(std::shared_ptr<ParetoElement> par, link* l);
  ParetoElement(std::shared_ptr<ParetoElement> par, link* l, bool shared);
  ParetoElement(double a, double b, double taud, double sa, double sb, double staud);
  double a() const;
  double b() const;
  double taud() const;
  double shared_a() const;
  double shared_b() const;
  double shared_taud() const;
  double k() const;
  bool operator<(const ParetoElement& other) const;
  std::unique_ptr<std::vector<link*>> collectLinks();
  std::shared_ptr<route> collectRoute();
  ParetoElement(std::shared_ptr<ParetoElement> par);  // recursive memory copy
  void set_parent(const std::shared_ptr<ParetoElement> par);
};
