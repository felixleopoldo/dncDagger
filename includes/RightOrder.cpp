#include "RightOrder.h"
using namespace std;

RightOrder::RightOrder(vector<int> &order,
                       double order_score,
                       vector<double> &node_scores,
                       size_t n) : order(order),
                                   order_score(order_score),
                                   node_scores(node_scores),
                                   n(n)
{
  inserted_max_order_scores = vector<double>(order.size());
  new_top_scores = vector<double>(order.size());
  best_insert_pos = vector<size_t>(order.size());
}

int RightOrder::front() const
{
  return order[order.size() - n];
}

size_t RightOrder::front_ind() const
{
  return order.size() - n;
}

ostream &operator<<(ostream &os, const RightOrder &ro)
{
  size_t p = ro.order.size();

  if (ro.n != p)
  {
    os << "([...],";
  }
  else
  {
    os << "(";
  }

  for (size_t i = p - ro.n; i < p; i++)
  {
    if (i != p - 1)
    {
      os << ro.order[i] << ", ";
    }
    else
    {
      os << ro.order[i];
    }
  }

  os << ")";
  return os;
}