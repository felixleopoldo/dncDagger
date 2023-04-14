#ifndef OPRUNER_LEFT_H
#define OPRUNER_LEFT_H

#include <iostream>
#include "auxiliary.h"
#include "OrderScoring.h"
#include "RightOrder.h"
#include "LeftOrder.h"

tuple<vector<int>, double, size_t, size_t> opruner_left(OrderScoring &scoring);

vector<LeftOrder> prune_equal_sets(vector<LeftOrder> left_orders, bool right_type, double);

std::vector<int> has_left_gaps(LeftOrder &lo, int new_node, std::vector<double> &bottom_scores, OrderScoring &scoring);

LeftOrder init_left_order(size_t node, OrderScoring &scoring);

bool optimal_back(const LeftOrder &lo, size_t new_node, OrderScoring &scoring);

std::tuple<std::vector<int>, double, size_t, size_t> sequential_opt_left(OrderScoring &scoring);

void sequential_opt_left_type(OrderScoring &scoring);

#endif