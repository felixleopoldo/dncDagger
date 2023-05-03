
#ifndef RIGHTORDER_H
#define RIGHTORDER_H

#include <cassert>
#include <chrono>
#include <iostream>
#include <cmath>
#include <vector>
#include "LeftOrder.h"
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
    vector<int>::reverse_iterator begin();
    vector<int>::reverse_iterator end();
    vector<int>::iterator rbegin();
    vector<int>::iterator rend();
    vector<int>::reverse_iterator hidden_begin();
    vector<int>::reverse_iterator hidden_end();
};


RightOrder operator+( RightOrder &c1,  LeftOrder &c2);
ostream &operator<<(ostream &os, const RightOrder &ro);

#endif

