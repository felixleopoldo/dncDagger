#include <cassert>
#include <chrono>
#include <iostream>
#include <cmath>
#include <RInside.h>
#include <Rcpp.h>
#include "OrderScoring.h"
#include "opruner_right.h"
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

    //printIsoComps(iso_comps);

    // iterate over isolated components references and find the component dependencies
    for (size_t i = 0; i < iso_comps.iso_comps.size(); i++)
    {
        IsoComp & iso_comp = iso_comps.iso_comps[i];
        component_dependence(iso_comp, scoring);
    }

    printIsoComps(iso_comps);

    return(h_min_adj);
}

void component_dependence(IsoComp & iso_comp, OrderScoring & scoring){

    size_t p = scoring.numparents.size();
    // create a matrix, comp_dep, of component dependencies
    iso_comp.comp_dep = vector<vector<int>>(iso_comp.subcomps.size(), vector<int>(iso_comp.subcomps.size(), 0));

    // iterate over the subcomponents
    for (size_t i = 0; i < iso_comp.subcomps.size(); i++)
    {
        SubComp & subcomp = iso_comp.subcomps[i];

        // If this subcomponent hasnt been touched before
        if (subcomp.score == 0) {

            // find the rest of the nodes, not in subcomp.nodes
            vector<int> restnodes;
            for (size_t j = 0; j < p; j++)
            {
                // add the node to restnodes if it is not in subcomp.nodes
                if (find(subcomp.nodes.begin(), subcomp.nodes.end(), j) == subcomp.nodes.end())
                {
                    restnodes.push_back(j);
                }
            }

            // run opruner_right conditioned on restnodes (initial nodes)
            const auto &[order, log_score, suborder, suborder_log_score, node_scores, max_n_particles, tot_n_particles] = opruner_right(scoring, restnodes);
            // find the correspodning DAG
            vector<vector<int>> adjmat = order_to_dag(order, scoring);

            // updated the isocom and subcomp with the new information
            subcomp.suborder = suborder; // is this the suborder?
            subcomp.score = suborder_log_score;
            subcomp.tot_n_particles = tot_n_particles;
            subcomp.max_n_particles = max_n_particles;
            subcomp.opt_adjmat = adjmat;

            // update the isolated component
            iso_comp.tot_n_particles += tot_n_particles;
            iso_comp.max_n_particles = max(iso_comp.max_n_particles, max_n_particles);
        }

        // // print optimal adjmat
        // cout << "Optimal adjmat for subcomponent " << i << endl;
        // for (size_t j = 0; j < subcomp.opt_adjmat.size(); j++)
        // {
        //     for (size_t k = 0; k < subcomp.opt_adjmat.size(); k++)
        //     {
        //         cout << subcomp.opt_adjmat[j][k] << " ";
        //     }
        //     cout << endl;
        // }

        // iterate over the subcomponents
        for (size_t j = 0; j < iso_comp.subcomps.size(); j++)
        {
            SubComp & subcomp2 = iso_comp.subcomps[j];

            if (iso_comp.connected) break;
            if (i == j) continue; // same subcomponent

            // check if the are any edges between nodes in the subcomponents 
            // with id i and j in subcomp.opt_adjmat
            
            for (size_t k = 0; k < subcomp.nodes.size(); k++)
            {
                for (size_t l = 0; l < subcomp2.nodes.size(); l++)
                {
                    if (subcomp.opt_adjmat[subcomp2.nodes[l]][subcomp.nodes[k]] == 1)
                    {
                        iso_comp.comp_dep[j][i] = 1;
                        continue;
                    }
                }
            }
        }
           
    }
}

void printIsoComps(IsoComps & iso_comps) {
    // iterate over isolated components and print their content
    for (size_t i = 0; i < iso_comps.iso_comps.size(); i++)
    {
        IsoComp & iso_comp = iso_comps.iso_comps[i];
        cout << "Isolated component " << i << " has " << iso_comp.nodes.size() << " nodes and " << iso_comp.subcomps.size() << " subcomponents." << endl;
        for (size_t j = 0; j < iso_comp.subcomps.size(); j++)
        {
            SubComp & subcomp = iso_comp.subcomps[j];
            cout << "Subcomponent " << j << " has " << subcomp.nodes.size() << " nodes." << endl;
            // pring the nodes
            cout << "Nodes: ";
            for (size_t k = 0; k < subcomp.nodes.size(); k++)
            {
                cout << subcomp.nodes[k] << " ";
            }
            cout << endl;
            // print all the other info of the subcomponent
            cout << "Suborder: ";
            for (size_t k = 0; k < subcomp.suborder.size(); k++)
            {
                cout << subcomp.suborder[k] << " ";
            }
            cout << endl;
            cout << "Scores: ";
            for (size_t k = 0; k < subcomp.scores.size(); k++)
            {
                cout << subcomp.scores[k] << " ";
            }
        }
        // print comp_dep matrix
        cout << "Component dependence matrix: " << endl;
        for (size_t j = 0; j < iso_comp.comp_dep.size(); j++)
        {
            for (size_t k = 0; k < iso_comp.comp_dep.size(); k++)
            {
                cout << iso_comp.comp_dep[j][k] << " ";
            }
            cout << endl;
        }
        cout << endl;
    }

}

IsoComps structure_components(Graph & G_H_min, Graph & G_H_max) {
    // get components of G_H_min and G_H_max
    vector<int> component_sub(num_vertices(G_H_min));
    int num_sub = connected_components(G_H_min, &component_sub[0]);
    // print the compnent vector
    cout << "Number of subcomponents: " << num_sub << endl; // "Number of components:
    cout << "SubComponent vector: ";
    for (size_t i = 0; i < component_sub.size(); i++)
    {
        cout << component_sub[i] << " ";
    }
    cout << endl;
    

    vector<int> component_max(num_vertices(G_H_max));
    int num_max = connected_components(G_H_max, &component_max[0]);
    // get isolated components
    IsoComps iso_comps;
    iso_comps.iso_comps = vector<IsoComp>(num_max);

    for (int i = 0; i < num_max; i++)
    {
        // Find all the nodes in this isoloated component
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

        // print the nodes in the subcomponent
        cout << "Subcomponent " << comp_id << " has " << subnodes.size() << " nodes." << endl;
        cout << "Nodes: ";
        for (size_t j = 0; j < subnodes.size(); j++)
        {
            cout << subnodes[j] << " ";
        }
        cout << endl;
        // print nodes in isocomp
        cout << "Nodes in isocomp: ";
        for (size_t j = 0; j < iso_comp.nodes.size(); j++)
        {
            cout << iso_comp.nodes[j] << " ";
        }
        cout << endl;

        // Add the subcomponent to the isolated component
        iso_comp.subcomps.push_back(subcomp);
        iso_comp.subcompids.push_back(comp_id);

    }

    return(iso_comps);
}