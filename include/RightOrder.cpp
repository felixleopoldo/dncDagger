#include "RightOrder.h"
#include "LeftOrder.h"
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

size_t RightOrder::size() const
{
  return n;
}

size_t RightOrder::size_hidden() const
{
  return order.size() - n;
}

vector<int>::reverse_iterator RightOrder::begin()
{
  return order.rbegin();
}

vector<int>::reverse_iterator RightOrder::end()
{
  return order.rbegin() + n;
}

vector<int>::iterator RightOrder::rbegin()
{
  return order.end() - n;
}
vector<int>::iterator RightOrder::rend()
{
  return order.end();
}

vector<int>::iterator RightOrder::hidden_begin()
{
  return order.begin() ;
}

vector<int>::iterator RightOrder::hidden_end()
{
  return order.begin() + order.size() - n;
}

RightOrder operator+(RightOrder &ro, LeftOrder &lo)
{
  vector<int> order;
  order.reserve(ro.n + lo.n);
  order.insert(order.end(), lo.begin(), lo.end());
  order.insert(order.end(), ro.rbegin(), ro.rend());

  double order_score = ro.order_score + lo.order_score;
  // This has to be fixe as well
  size_t n = ro.n + lo.n;

  vector<double> node_scores(n, -1);
  // loop trought all the nodes.
  for (auto i = lo.begin(); i != lo.end(); i++)
  {
    node_scores[*i] = lo.node_scores[*i];
  }
  for (auto i = ro.begin(); i != ro.end(); i++)
  {
    node_scores[*i] = ro.node_scores[*i];
  }

  RightOrder c(order, order_score, node_scores, n);
  return c;
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