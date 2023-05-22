
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
    double upper_bound;
    double upper_bound_hidden;
    double hidden_top_score_sum;
   
    std::vector<double> node_scores;
    size_t n_nodes;
    std::vector<double> inserted_max_order_scores;
    std::vector<double> new_top_scores;
    std::vector<size_t> best_insert_pos;


    RightOrder(std::vector<int> &order,
               double order_score,
               std::vector<double> &node_scores,
               size_t n);

    int front() const;

    size_t front_ind() const;
    size_t size() const;
    size_t size_hidden() const;
    vector<int>::reverse_iterator begin();
    vector<int>::reverse_iterator end();
    vector<int>::iterator rbegin();
    vector<int>::iterator rend();
    vector<int>::iterator hidden_begin();
    vector<int>::iterator hidden_end();
};


RightOrder operator+( RightOrder &c1,  LeftOrder &c2);
RightOrder operator+(LeftOrder &lo, RightOrder &ro);
ostream &operator<<(ostream &os, const RightOrder &ro);

#endif

