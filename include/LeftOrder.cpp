
#include "LeftOrder.h"

LeftOrder::LeftOrder(vector<int> &order,
                     double order_score,
                     vector<double> &node_scores,
                     size_t n) : order(order),
                                 order_score(order_score),
                                 node_scores(node_scores),
                                 n(n)
{
  inserted_max_order_scores = vector<double>(order.size());
  new_back_scores = vector<double>(order.size());
  best_insert_pos = vector<size_t>(order.size());
}

int LeftOrder::back() const
{
  return order[n - 1];
}

size_t LeftOrder::back_ind() const
{
  return n - 1;
}
