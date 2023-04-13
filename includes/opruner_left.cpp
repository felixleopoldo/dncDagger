#include <cassert>
#include <chrono>
#include <iostream>
#include "auxiliary.h"
#include "OrderScoring.h"
#include "LeftOrder.h"
#include "opruner_right.h"

tuple<vector<int>, double, size_t, size_t> opruner_left(OrderScoring &scoring)
{
}

vector<LeftOrder> prune_equal_sets(vector<LeftOrder> left_orders,
                                   bool right_type, double EPSILON)
{
    vector<vector<bool>> boolmat;
    vector<double> order_scores;

    // cout << "Creating boolmatrix" << endl;
    for (const LeftOrder &ro : left_orders)
    {
        vector<bool> boolvec = order_to_boolvec(ro.order, ro.n, right_type);
        boolmat.push_back(move(boolvec));
        order_scores.push_back(ro.order_score);
    }

    // cout << "Get indices of unique maximal scoring sets" << endl;
    vector<int> pruned_inds = unique_sets(boolmat, order_scores, EPSILON);
    vector<double> pruned_scores;
    vector<LeftOrder> kept_ros;

    for (const auto &ind : pruned_inds)
    {
        kept_ros.push_back(left_orders[ind]);
    }

    return (kept_ros);
}