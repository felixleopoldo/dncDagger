#ifndef OPRUNER_LEFT_H
#define OPRUNER_LEFT_H

#include <iostream>
#include "auxiliary.h"
#include "OrderScoring.h"
#include "RightOrder.h"

void swap_nodes(const int lower, const int upper, RightOrder &ro, OrderScoring &scoring);

bool optimal_front(const RightOrder &ro, size_t new_node, OrderScoring &scoring);

bool optimal_front(const RightOrder &ro, OrderScoring &scoring);

bool independent_front(const RightOrder &ro, const vector<double> &top_scores, OrderScoring &scoring);

void make_visible(int from_index, int to_index, RightOrder &ro, OrderScoring &scoring);

double max_inserted_score(const RightOrder &ro, size_t node_ind, OrderScoring &scoring);

bool equal_and_unordered_top(const RightOrder &ro, int new_front, OrderScoring &scoring);

bool equal_and_unordered_top(RightOrder &ro, OrderScoring &scoring);

bool has_hidden_gap(const RightOrder &ro, const vector<double> &top_scores, OrderScoring &scoring);

bool unord_hidden_gap(const RightOrder &ro, const vector<double> &top_scores, OrderScoring &scoring);

RightOrder init_right_order(size_t node, OrderScoring &scoring);

void update_insertion_scores(RightOrder &ro, OrderScoring &scoring);

RightOrder add_node_in_front(const RightOrder &ro_prev, size_t index_of_el_to_insert, OrderScoring &scoring);

vector<bool> order_to_boolvec(const vector<int> &order, size_t n, bool right_type);

vector<RightOrder> prune_equal_sets(vector<RightOrder> right_orders, bool right_type);

vector<int> no_right_gaps(RightOrder &ro, vector<double> &top_scores, OrderScoring &scoring);

tuple<vector<int>, double, size_t, size_t> opruner_right(OrderScoring &scoring);

// [[Rcpp::plugins(cpp17)]]

// [[Rcpp::export]]

Rcpp::List r_opruner_right(Rcpp::List ret);

#endif