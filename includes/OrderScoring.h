
#ifndef A2DD_H
#define A2DD_H

#include <cassert>
#include <chrono>
#include <iostream>
#include <cmath>

using namespace std;

class OrderScoring
{
private:
    vector<vector<int>> potential_parents;
    vector<Rcpp::IntegerVector> rowmaps_backwards;
    vector<vector<int>> potential_plus1_parents;

    bool MAP;

public:
    vector<int> numparents;
    vector<vector<vector<double>>> scoretable;
    vector<vector<vector<double>>> scoresmatrices;

    OrderScoring(
        vector<vector<int>> potential_parents,
        vector<int> numparents,
        vector<Rcpp::IntegerVector> rowmaps_backwards,
        vector<vector<int>> potential_plus1_parents,
        vector<vector<vector<double>>> scoretable,
        vector<vector<vector<double>>> scoresmatrices,
        bool MAP); 
        // : potential_parents(potential_parents),
        //             rowmaps_backwards(rowmaps_backwards),
        //             potential_plus1_parents(potential_plus1_parents),
        //             MAP(MAP),
        //             numparents(numparents),
        //             scoretable(scoretable),
        //             scoresmatrices(scoresmatrices);

    /**
     * Re-calculation scores after swapping up node_a so that (a, b) -> (b, a).
     *
     * Order of node1 is lower than node2.
     *
     * Returns: (nodea_score, nodeb_score)
     */
    tuple<double, double> swap_nodes(int nodea_index, int nodeb_index,
                                     vector<int> &ordering,
                                     const vector<double> &node_scores);
  

    vector<double> score(const vector<int> &ordering, const size_t &from_orderpos, const size_t &n_elements) const;

    /**
     * position is the index in the ordering of the node.
     */
    double score_pos(const vector<int> &ordering, const size_t &position) const;

    double sum_log_probs(const vector<double> &log_probs) const;

    vector<int> get_plus1_indices(const int &position, const vector<int> &ordering) const;

    int get_f_bar_z(const int &position, const vector<int> &ordering) const;
};

#endif

OrderScoring get_score(Rcpp::List ret);
