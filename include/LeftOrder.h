
#ifndef LEFTORDER_H
#define LEFTORDER_H

#include <cassert>
#include <chrono>
#include <iostream>
#include <cmath>
#include<vector>

using namespace std;

class LeftOrder
{
public:
    vector<int> order;
    double order_score;
    vector<double> node_scores;
    size_t n;
    vector<double> inserted_max_order_scores;
    vector<double> new_back_scores;
    vector<size_t> best_insert_pos;
    LeftOrder(vector<int> &order,
              double order_score,
              vector<double> &node_scores,
              size_t n);

    int back() const;

    size_t back_ind() const;
    vector<int>::iterator begin();
    vector<int>::iterator end();
    vector<int>::iterator hidden_begin();
    vector<int>::iterator hidden_end();
};

#endif

ostream &operator<<(ostream &os, const LeftOrder &ro);
