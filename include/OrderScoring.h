
#ifndef ORDERSCORING_H
#define ORDERSCORING_H

#include <cassert>
#include <chrono>
#include <iostream>
#include <cmath>
#include <Rcpp.h>

using namespace std;

class OrderScoring
{
private:
    vector<vector<int>> potential_plus1_parents;
    bool MAP;

public:
    vector<vector<int>> potential_parents;
    vector<Rcpp::IntegerVector> rowmaps_backwards;
    vector<Rcpp::IntegerVector> rowmaps_forward;
    vector<int> numparents;
    vector<vector<vector<double>>> scoretable;
    vector<vector<vector<double>>> scoresmatrices;
    vector<vector<vector<double>>> maxmatrix;
    vector<vector<vector<int>>> maxrow;
    vector<vector<vector<int>>> parenttable;

    /**
     * Re-calculation scores after swapping up node_a so that (a, b) -> (b, a).
     *
     * Order of node1 is lower than node2.
     *
     * Returns: (nodea_score, nodeb_score)
     */
    OrderScoring(vector<vector<int>> potential_parents, 
    vector<int> numparents, 
    vector<Rcpp::IntegerVector> rowmaps_backwards, 
    vector<Rcpp::IntegerVector> rowmaps_forward,
    vector<vector<int>> potential_plus1_parents, 
    vector<vector<vector<double>>> scoretable, 
    vector<vector<vector<double>>> scoresmatrices, 
    vector<vector<vector<double>>> maxmatrix, 
    vector<vector<vector<int>>> maxrow, 
    vector<vector<vector<int>>> parenttable, 
    bool MAP);

    tuple<double, double> swap_nodes(int nodea_index, int nodeb_index,
                                     vector<int> &ordering,
                                     const vector<double> &node_scores);

    vector<double> score(const vector<int> &ordering, const size_t &from_orderpos, const size_t &n_elements) const;

    double score_order(const vector<int> &ordering, const size_t &from_orderpos, const size_t &n_elements) const;

    /**
     * position is the index in the ordering of the node.
     */
    double score_pos(const vector<int> &ordering, const size_t &position) const;

    vector<int> get_opt_parents(const vector<int> &ordering, const size_t &position) const;

    double sum_log_probs(const vector<double> &log_probs) const;

    vector<int> get_plus1_indices(const int &position, const vector<int> &ordering) const;

    tuple<vector<int>, vector<int>> get_active_and_banned(const int &position, const vector<int> &ordering) const;

    //int get_f_bar_z(const int node, const vector<int> & parent_indices_banned_by_ordering) const;
    int get_f_bar_z(const int &, const vector<int> & ordering) const;

};


OrderScoring get_score(Rcpp::List ret);

#endif