#ifndef OPRUNER_LEFT_H
#define OPRUNER_LEFT_H

#include <iostream>
#include "auxiliary.h"
#include "OrderScoring.h"
#include "RightOrder.h"
#include "LeftOrder.h"

tuple<vector<int>, double, size_t, size_t> opruner_left(OrderScoring &scoring);

vector<LeftOrder> prune_equal_sets(vector<LeftOrder> left_orders,
                                   bool right_type, double);

#endif