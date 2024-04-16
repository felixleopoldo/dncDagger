#include <cassert>
#include <chrono>
#include <iostream>
#include <cmath>
#include <RInside.h>
#include <Rcpp.h>
#include "OrderScoring.h"
#include "auxiliary.h"
#include "dnc.h"

using namespace std;
using namespace boost;

//typedef adjacency_list< vecS, vecS, directedS > Graph;


/*  Divide and conquer algorithm.
    Returns the optimal DAG.
*/
vector<vector<bool>> dnc(OrderScoring &scoring,
                        vector<vector<bool>> &h_min_adj,
                        vector<vector<bool>> &h_max_adj)
{
    size_t p = h_min_adj.size();
    Graph G_H_min(p);
    Graph G_H_max(p);
    // read h_min_adj and h_max_adj into boost graphs
    for (size_t i = 0; i < h_min_adj.size(); i++)
    {
        for (size_t j = 0; j < h_min_adj[i].size(); j++)
        {
            if (h_min_adj[i][j])
            {
                add_edge(i, j, G_H_min);
            }
            if (h_max_adj[i][j])
            {
                add_edge(i, j, G_H_max);
            }
        }
    }

    IsoComps iso_comps = structure_components(G_H_min, G_H_max);

    // iterate over isolated components and print their content
    for (size_t i = 0; i < iso_comps.iso_comps.size(); i++)
    {
        IsoComp & iso_comp = iso_comps.iso_comps[i];
        cout << "Isolated component " << i << " has " << iso_comp.nodes.size() << " nodes and " << iso_comp.subcomps.size() << " subcomponents." << endl;
        for (size_t j = 0; j < iso_comp.subcomps.size(); j++)
        {
            SubComp & subcomp = iso_comp.subcomps[j];
            cout << "Subcomponent " << j << " has " << subcomp.nodes.size() << " nodes." << endl;
        }
    }

    return(h_min_adj);
}

IsoComps structure_components(Graph & G_H_min, Graph & G_H_max) {
    // get components of G_H_min and G_H_max
    vector<int> component_sub(num_vertices(G_H_min));
    int num_sub = connected_components(G_H_min, &component_sub[0]);

    vector<int> component_max(num_vertices(G_H_max));
    int num_max = connected_components(G_H_max, &component_max[0]);
    // get isolated components
    IsoComps iso_comps;
    iso_comps.iso_comps = vector<IsoComp>(num_max);

    for (int i = 0; i < num_max; i++)
    {
        // Find all the nodes in this isloated component
        vector<int> subnodes;
        for (size_t j = 0; j < num_vertices(G_H_max); j++)
        {
            if (component_max[j] == i)
            {
                subnodes.push_back(j);
            }
        }

        IsoComp iso_comp;
        iso_comp.nodes = subnodes;
        iso_comp.subcompids = vector<int>();
        iso_comp.scores = vector<double>();
        iso_comp.subcomps = vector<SubComp>();

        iso_comps.iso_comps[i] = iso_comp;
    }

    
    // add sub components to the right isolated components
    // iterate over subcomp ids
    for (int comp_id = 0; comp_id < num_sub; comp_id++)
    {
        // store all the nodes in sub component i in a vector called subnodes
        vector<int> subnodes;
        for (size_t j = 0; j < num_vertices(G_H_min); j++)
        {
            if (component_sub[j] == comp_id)
            {
                subnodes.push_back(j);
            }
        }
        
        SubComp subcomp;
        subcomp.nodes = subnodes;

        // Now get the Isolated component that this is a subcomponent of
        // Just take the first node in the subcomponent and get the isolated component
        IsoComp & iso_comp = iso_comps.iso_comps[component_max[subnodes[0]]];

        // Add the subcomponent to the isolated component
        iso_comp.subcomps.push_back(subcomp);
        iso_comp.subcompids.push_back(comp_id);

    }

    return(iso_comps);
}