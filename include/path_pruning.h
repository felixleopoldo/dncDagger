#ifndef PATH_PRUNER_H
#define PATH_PRUNER_H

#include <cassert>
#include <chrono>
#include <iostream>
#include "auxiliary.h"
#include "OrderScoring.h"
#include "RightOrder.h"
#include "LeftOrder.h"

#include <algorithm>
#include <iterator>

#include <boost/config.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/prim_minimum_spanning_tree.hpp>
#include <boost/property_map/property_map.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/topological_sort.hpp>
#include <boost/graph/breadth_first_search.hpp>

#include "edmonds/edmonds_optimum_branching.hpp"
#include "edmonds/edmonds_optimum_branching_impl.hpp"
// Define a directed graph type that associates a weight with each
// edge. We store the weights using internal properties as described
// in BGL.
typedef boost::property<boost::edge_weight_t, double> EdgeProperty;
typedef boost::adjacency_list<boost::listS,
                              boost::vecS,
                              boost::bidirectionalS, // OBS bidirectional to provide in_degree
                              boost::no_property,
                              EdgeProperty>
    DirectedGraph;
typedef boost::graph_traits<DirectedGraph>::vertex_descriptor Vertex;
typedef boost::graph_traits<DirectedGraph>::edge_descriptor Edge;

typedef boost::adjacency_list<boost::vecS,
                              boost::vecS,
                              boost::undirectedS,
                              boost::property<boost::vertex_distance_t, int>,
                              EdgeProperty>
    UndirectedGraph;
typedef boost::graph_traits<UndirectedGraph>::edge_descriptor UEdge;

vector<double> get_unrestricted_vec(size_t p, OrderScoring &scoring);

vector<double> all_restr(size_t p, OrderScoring &scoring);

LeftOrder extract_leftorder(RightOrder &ro, RightOrder &reference_order, OrderScoring &scoring);

UndirectedGraph prim(UndirectedGraph &g);

LeftOrder topo_left_order(DirectedGraph &g, RightOrder &ro, OrderScoring &scoring);

LeftOrder get_toporder_utree(UndirectedGraph &G, int n, RightOrder &ro, OrderScoring &scoring);

UndirectedGraph get_boost_ugraph(vector<vector<double>> &my_weights, RightOrder &ro);

DirectedGraph get_boost_dgraph(vector<vector<double>> &myweights, RightOrder &ro);

DirectedGraph edmonds(DirectedGraph &G);

vector<vector<double>> get_unrestr_mat(size_t p, OrderScoring &scoring);

vector<vector<double>> get_hard_restr_mat(size_t p, OrderScoring &scoring);

vector<int> get_suborder_prim(RightOrder &ro, vector<vector<double>> &M);

vector<RightOrder> prune_path(RightOrder &reference_order,
                              vector<RightOrder> &right_orders,
                              vector<vector<double>> &M,
                              vector<vector<double>> &H,
                              vector<double> &bottom_scores,
                              vector<double> &top_scores,
                              OrderScoring &scoring);

#endif