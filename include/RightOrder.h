
#ifndef RIGHTORDER_H
#define RIGHTORDER_H

#include <cassert>
#include <chrono>
#include <iostream>
#include <cmath>
#include <vector>
/**
 *
 * Particle struct
 *
 */
class RightOrder
{
public:
    std::vector<int> order;
    double order_score;
    std::vector<double> node_scores;
    size_t n;
    std::vector<double> inserted_max_order_scores;
    std::vector<double> new_top_scores;
    std::vector<size_t> best_insert_pos;
    RightOrder(std::vector<int> &order,
               double order_score,
               std::vector<double> &node_scores,
               size_t n);

    int front() const;

    size_t front_ind() const;
};

#endif