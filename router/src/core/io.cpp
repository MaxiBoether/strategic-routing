#include "core/io.h"
#include "core/data.h"
#include "core/globals.h"
#include "core/data.h"

#include <algorithm>
#include <cstdlib>
#include <cstring>
#include <vector>
#include <numeric>
#include <iterator>
#include <tinyxml2.h>

using namespace tinyxml2;

XMLDocument plansXml = XMLDocument();


class LinkCollector : public XMLVisitor {
    void handleLink(const XMLElement& e) {
        link* l = new link;// psychological_model.newLink();
        l->id = e.IntAttribute("id");
        l->from = e.IntAttribute("from");
        l->to = e.IntAttribute("to");
        l->length = std::atof(e.Attribute("length"));
        l->capacity = std::atof(e.Attribute("capacity"));
        l->freespeed = std::atof(e.Attribute("freespeed"));
        adj.at(l->from).push_back(l);
    }
    void handleNode(const XMLElement& e) {
        node* n = new node;
        n->x = strdup(e.Attribute("x"));
        n->y = strdup(e.Attribute("y"));
        nodes[e.IntAttribute("id")] = n;
    }
    virtual bool VisitEnter(const XMLElement& e, const XMLAttribute *attrs) override {
        (void)attrs;
        if (std::strcmp(e.Name(), "node") == 0)
            handleNode(e);
        if (std::strcmp(e.Name(), "link") == 0)
            handleLink(e);
        return true;
    };
};

class NodeCounter : public XMLVisitor {
    virtual bool VisitEnter(const XMLElement& e, const XMLAttribute *attrs) override {
        (void)attrs;
        if (std::strcmp(e.Name(), "node") == 0)
            maxId = std::max(maxId, e.IntAttribute("id"));
        return true;
    };
public:
    int maxId = 0;
};

void loadGraph(const char* graphFile) {
    XMLDocument input;
    input.LoadFile(graphFile);
    NodeCounter cnt;
    input.Accept(&cnt);
    nodes.resize(cnt.maxId+1);
    adj.resize(cnt.maxId+1);
    LinkCollector linkCollector;
    input.Accept(&linkCollector);
}

void _handlePersonElement(const XMLElement *e)
{
    if (std::strcmp(e->Name(), "person") == 0) {
        if (e->Attribute("origin_node")) {
            persons.emplace_back(
                    e->IntAttribute("origin_node"),
                    e->IntAttribute("destination_node"),
                    e->Attribute("time"),
                    const_cast<XMLElement*>(e)
                    );
        }
    }
    if (e->NextSiblingElement())
        _handlePersonElement(e->NextSiblingElement());
}

void readPersons(const XMLDocument& doc) {
    const XMLElement *root = doc.RootElement();
    if (root->FirstChildElement())
        _handlePersonElement(root->FirstChildElement());
}

void loadPlans(const char* plansFile) {
    plansXml.LoadFile(plansFile);
    readPersons(plansXml);
}

void outputPlansToNewFile(const char* plansFile) {
    XMLDocument xml;

    // declaration
    XMLDeclaration* decl = xml.NewDeclaration("xml version=\"1.0\" ");
    xml.InsertEndChild(decl);

    // doctype
    XMLUnknown* doctype = xml.NewUnknown("DOCTYPE plans SYSTEM \"http://www.matsim.org/files/dtd/plans_v4.dtd\"");
    xml.InsertEndChild(doctype);

    // root element
    XMLElement* root = xml.NewElement("plans");
    xml.InsertEndChild(root);

    int i = 0;
    for (auto p : persons) {
        XMLElement* pelem = xml.NewElement("person");
        pelem->SetAttribute("id", i);
        root->InsertEndChild(pelem);
        p.write(xml, pelem);
        i++;
    }
    xml.SaveFile(plansFile);
}

void outputPlans(const char* plansFile) {
    for (auto p : persons) {
        p.write();
    }
    plansXml.SaveFile(plansFile);
}

