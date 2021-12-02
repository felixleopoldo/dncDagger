#include <bits/stdc++.h>
#include <RInside.h>
#include <cassert>
#include <chrono>
#include <thread>
#include <iostream>
#include "thread_pool.hpp"
//#include "OrderScoring.cpp"
// #include <cstdio>
// #include <algorithm>

using namespace std::chrono;

double EPSILON = 0.0000001;

// // Function to find the nCr
// void printNcR(int n, int r)
// {

//     // p holds the value of n*(n-1)*(n-2)...,
//     // k holds the value of r*(r-1)...
//     long long p = 1, k = 1;

//     // C(n, r) == C(n, n-r),
//     // choosing the smaller value
//     if (n - r < r)
//         r = n - r;

//     if (r != 0)
//     {
//         while (r)
//         {
//             p *= n;
//             k *= r;

//             // gcd of p, k
//             long long m = __gcd(p, k);

//             // dividing by gcd, to simplify
//             // product division by their gcd
//             // saves from the overflow
//             p /= m;
//             k /= m;

//             n--;
//             r--;
//         }

//         // k should be simplified to 1
//         // as C(n, r) is a natural number
//         // (denominator should be 1 ) .
//     }

//     else
//         p = 1;

//     // if our approach is correct p = ans and k =1
//     std::cout << p << std::endl;
// }

bool approximatelyEqual(float a, float b, float epsilon)
{
    return fabs(a - b) <= ((fabs(a) < fabs(b) ? fabs(b) : fabs(a)) * epsilon);
}

bool essentiallyEqual(float a, float b, float epsilon)
{
    return fabs(a - b) <= ((fabs(a) > fabs(b) ? fabs(b) : fabs(a)) * epsilon);
}

bool definitelyGreaterThan(float a, float b, float epsilon)
{
    return (a - b) > ((fabs(a) < fabs(b) ? fabs(b) : fabs(a)) * epsilon);
}

bool definitelyLessThan(float a, float b, float epsilon)
{
    return (b - a) > ((fabs(a) < fabs(b) ? fabs(b) : fabs(a)) * epsilon);
}

int rand_int_in_range(const std::size_t &from, const std::size_t &to)
{
    return (from + (std::rand() % (to - from + 1)));
}

// CPP program To calculate The Value Of nCr
int fact(int n);

int nCr(int n, int r)
{
    return fact(n) / (fact(r) * fact(n - r));
}

// Returns factorial of n
int fact(int n)
{
    int res = 1;
    for (int i = 2; i <= n; i++)
        res = res * i;
    return res;
}

template <typename T>
void myswap(std::size_t i, std::size_t j, std::vector<T> &v)
{
    T tmp = v[i];
    v[i] = v[j];
    v[j] = tmp;
}

std::vector<std::size_t> parts(std::size_t m, std::size_t size)
{
    std::vector<std::size_t> res(size);

    for (std::size_t i = 0; i != m; ++i)
    {
        ++res[i % size];
    }
    return res;
}

std::size_t whichPart(std::size_t m, std::size_t size, std::size_t n)
{
    std::size_t index = 0;
    for (auto i : parts(m, size))
    {
        if (n < i)
        {
            return index;
        }
        ++index;
        n -= i;
    }
    throw std::runtime_error("invalid argument");
}

template <typename T>
void PrintVector(std::vector<T> &arr, int from = -1, int to = -1)
{
    if (from == -1)
    {
        copy(arr.begin(), arr.end(), std::ostream_iterator<T>(std::cout, ","));
    }
    else
    {
        copy(arr.begin() + from, arr.begin() + to, std::ostream_iterator<T>(std::cout, ","));
    }
    std::cout << std::endl;
}

template <typename T>
void PrintVector(const std::vector<T> &arr, int from = -1, int to = -1)
{
    const std::vector<T> vec(arr);
    PrintVector(vec, from, to);
}

/**
 *
 * Particle struct may be
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
               size_t n) : order(order),
                           order_score(order_score),
                           node_scores(node_scores),
                           n(n)
    {
        inserted_max_order_scores = std::vector<double>(order.size());
        new_top_scores = std::vector<double>(order.size());
        best_insert_pos = std::vector<size_t>(order.size());
    }

    int front() const
    {
        return order[order.size() - n];
    }

    size_t front_ind() const
    {
        return order.size() - n;
    }
};

std::ostream &operator<<(std::ostream &os, const RightOrder &ro)
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

class OrderScoring
{
private:
    std::vector<std::vector<int>> potential_parents;
    std::vector<Rcpp::IntegerVector> rowmaps_backwards;
    std::vector<std::vector<int>> potential_plus1_parents;

    bool MAP;

public:
    std::vector<int> numparents;
    std::vector<std::vector<std::vector<double>>> scoretable;
    std::vector<std::vector<std::vector<double>>> scoresmatrices;

    OrderScoring(
        std::vector<std::vector<int>> potential_parents,
        std::vector<int> numparents,
        std::vector<Rcpp::IntegerVector> rowmaps_backwards,
        std::vector<std::vector<int>> potential_plus1_parents,
        std::vector<std::vector<std::vector<double>>> scoretable,
        std::vector<std::vector<std::vector<double>>> scoresmatrices,
        bool MAP) : potential_parents(potential_parents),
                    rowmaps_backwards(rowmaps_backwards),
                    potential_plus1_parents(potential_plus1_parents),
                    MAP(MAP),
                    numparents(numparents),
                    scoretable(scoretable),
                    scoresmatrices(scoresmatrices)
    {
    }

    /**
     * Re-calculation scores after swapping up node_a so that (node_a, node_b) --> (node_b, node_a).
     *
     * Order of node1 is lower than node2.
     *
     * Returns: (nodea_score, nodeb_score)
     */
    std::tuple<double, double> swap_nodes(int nodea_index, int nodeb_index,
                                          std::vector<int> &ordering,
                                          const std::vector<double> &node_scores)
    {
        int node_a = ordering[nodea_index];
        int node_b = ordering[nodeb_index];
        double nodea_score = 0.0; // just to initialize
        double nodeb_score = 0.0; // just to initialize

        if (MAP == true)
        {
            myswap(nodea_index, nodeb_index, ordering);
            nodea_score = score_pos(ordering, nodeb_index);
            nodeb_score = score_pos(ordering, nodea_index);
            return (std::make_tuple(nodea_score, nodeb_score));
        }
        else
        {
            // Computing score for nodea, which is moved up
            // If b is a potential parent for a, we have to recompute the scores since b is now banned.
            if (std::find(potential_parents[node_a].begin(), potential_parents[node_a].end(), node_b) != potential_parents[node_a].end())
            {
                myswap(nodea_index, nodeb_index, ordering);
                nodea_score = score_pos(ordering, nodeb_index);
                myswap(nodea_index, nodeb_index, ordering);
            }
            else
            { // Since b is not a potential parent of a, f_bar_z is not altered.
                // Check if b was a plus1 parent for a
                std::vector<int>::iterator itr = std::find(potential_plus1_parents[node_a].begin(), potential_plus1_parents[node_a].end(), node_b);
                if (itr != potential_plus1_parents[node_a].end())
                {
                    // Subtract the bs plus1 score contibution.

                    myswap(nodea_index, nodeb_index, ordering);
                    int f_bar_z = get_f_bar_z(nodeb_index, ordering); // ok?
                    // find the correct index j and take it.
                    int plus1_parents_index = std::distance(potential_plus1_parents[node_a].begin(), itr);
                    double nodeb_as_plus1_score = scoresmatrices[node_a][f_bar_z][plus1_parents_index + 1];
                    myswap(nodea_index, nodeb_index, ordering);

                    // Sice there new additin is so small I get precision error, sice we get log(1-0.99999999999999) = -inf
                    double max = std::max(node_scores[node_a], nodeb_as_plus1_score); // I should use order_scores not node scores here. Or something. This is not right at least.
                    // std::cout << nodeb_as_plus1_score - node_scores[node_a] << " " << std::endl;
                    if (std::abs(nodeb_as_plus1_score - node_scores[node_a]) > 0.000001)
                    { // 0.000001 is arbitrary
                        // OK
                        // std::cout << "NO RECOMPUTE order score for node" << std::endl;
                        nodea_score = std::log(std::exp(node_scores[node_a] - max) - std::exp(nodeb_as_plus1_score - max)) + max; // gets inf... 0 at node_scores[node_a] but something at node_scores[node_b]
                    }
                    else
                    {
                        // round off error. Recompute.
                        // std::cout << "RECOMPUTE order score for node" << std::endl;
                        myswap(nodea_index, nodeb_index, ordering);
                        nodea_score = score_pos(ordering, nodeb_index);
                        myswap(nodea_index, nodeb_index, ordering);
                    }

                    // std::cout << "true score " << true_score << " calcuated score " << node_scores[node_a] << std::endl;
                    // assert(std::abs(nodea_score-true_score) < 0.001);
                }
            }

            // // Computing score for node_b, which is moved down
            if (std::find(potential_parents[node_b].begin(), potential_parents[node_b].end(), node_a) != potential_parents[node_b].end())
            {
                myswap(nodea_index, nodeb_index, ordering);
                nodeb_score = score_pos(ordering, nodea_index);
                myswap(nodea_index, nodeb_index, ordering);
            }
            else
            {
                std::vector<int>::iterator itr = std::find(potential_plus1_parents[node_b].begin(),
                                                           potential_plus1_parents[node_b].end(),
                                                           node_a);
                if (itr != potential_plus1_parents[node_b].cend())
                {
                    int plus1_parents_index = std::distance(potential_plus1_parents[node_b].begin(), itr); // since no parents is the first
                    myswap(nodea_index, nodeb_index, ordering);
                    int f_bar_z = get_f_bar_z(nodea_index, ordering);
                    double nodea_as_plus1_score = scoresmatrices[node_b][f_bar_z][plus1_parents_index + 1];
                    myswap(nodea_index, nodeb_index, ordering);
                    double max = std::max(node_scores[node_b], nodea_as_plus1_score);

                    nodeb_score = std::log(std::exp(node_scores[node_b] - max) + std::exp(nodea_as_plus1_score - max)) + max;
                }
            }

            // std::swap(ordering[nodea_index], ordering[nodeb_index]); // last swap
            myswap(nodea_index, nodeb_index, ordering);

            return (std::make_tuple(nodea_score, nodeb_score));
        }
    }

    std::vector<double> score(const std::vector<int> &ordering, const std::size_t &from_orderpos, const std::size_t &n_elements) const
    {
        std::size_t n = ordering.size();
        std::vector<double> orderscores = std::vector<double>(n, 0.0); // O(p)           // orderscores <- vector("double", n)
        std::vector<int> active_plus1_parents_indices;                 // active_plus1_parents_indices < -vector("list", n)
        int f_bar_z;

        for (std::size_t position = from_orderpos; position < from_orderpos + n_elements; ++position)
        {
            int node = ordering[position];
            if (position == n - 1)
            {
                // no parents allowed, i.e.only first row, only first list
                orderscores[node] = scoretable[node][0][0]; // orderscores[node] <- scoretable [[node]][[1]][1, 1]
                // active_plus1_parents_indices[node].assign({0});           // active_plus1_parents_indices [[node]] <- c(1)
                // f_bar_z[node] = std::pow(2, potential_parents[node].size()); // f_bar_z[node] <- c(2 ^ numparents[node])
            }
            else
            {
                f_bar_z = get_f_bar_z(position, ordering);
                active_plus1_parents_indices = get_plus1_indices(position, ordering);
                std::vector<double> plus1_parents_scores((active_plus1_parents_indices).size());
                for (std::size_t j = 0; j < plus1_parents_scores.size(); j++)
                {
                    plus1_parents_scores[j] = scoresmatrices[node][f_bar_z][(active_plus1_parents_indices)[j]]; // allowedscorelist is in numerical order
                }

                if (MAP == true)
                {
                    orderscores[node] = *std::max_element(plus1_parents_scores.begin(), plus1_parents_scores.end());
                }
                else
                {
                    orderscores[node] = sum_log_probs(plus1_parents_scores);
                }
            }
        }
        return (orderscores);
    }

    double score_pos(const std::vector<int> &ordering, const std::size_t &position) const
    {
        double orderscore(0.0); // O(p)           // orderscores <- vector("double", n)
        int f_bar_z;
        int node = ordering[position];

        if (position == ordering.size() - 1)
        {
            // no parents allowed, i.e.only first row, only first list
            orderscore = scoretable[node][0][0]; // orderscores[node] <- scoretable [[node]][[1]][1, 1]
        }
        else
        {
            f_bar_z = get_f_bar_z(position, ordering);
            std::vector<int> active_plus1_parents_indices = get_plus1_indices(position, ordering);
            std::vector<double> plus1_parents_scores((active_plus1_parents_indices).size());

            for (std::size_t j = 0; j < plus1_parents_scores.size(); j++)
            {
                plus1_parents_scores[j] = scoresmatrices[node][f_bar_z][active_plus1_parents_indices[j]]; // allowedscorelist is in numerical order
            }
            if (MAP == true)
            {
                orderscore = *std::max_element(plus1_parents_scores.begin(), plus1_parents_scores.end());
            }
            else
            {
                orderscore = sum_log_probs(plus1_parents_scores);
            }
        }

        return orderscore;
    }

    double sum_log_probs(const std::vector<double> &log_probs) const
    {
        double max_log_prob = *std::max_element(log_probs.begin(), log_probs.end());
        double score(0.0);
        for (auto &val : log_probs)
        {
            score += std::exp(val - max_log_prob);
        }
        return (max_log_prob + std::log(score));
    }

    std::vector<int> get_plus1_indices(const int &position, const std::vector<int> &ordering) const
    {

        // Find the possible parents after in the ordering (allowedscorelist)
        // potential_plus1_parents is a list of possible parent (sets?) (set with extra +1 parent) for node (node ocnsidering the ordering).
        // allowedscorelist is a filetered version of plius1listparents, removing the node before in the ordering
        const int node = ordering[position];
        std::vector<int> active_plus1_parents_indices;
        (active_plus1_parents_indices).push_back(0);                           // f(null)=0, no )=parents is always a possibility.?
        for (std::size_t j = 0; j < potential_plus1_parents[node].size(); j++) // Is j the plus1node? -No, but potential_plus1_parents[node][j] is.
        {
            if (std::find(ordering.begin() + position + 1, ordering.end(), potential_plus1_parents[node][j]) != ordering.end())
            {
                (active_plus1_parents_indices).push_back(j + 1); // +1 since the 0 is no +1 parent at all?
            }
        }
        return (active_plus1_parents_indices);
    }

    int get_f_bar_z(const int &position, const std::vector<int> &ordering) const
    {
        int f_bar_z;
        const int node = ordering[position];
        // Find the banned scores before in the ordering. (Double banned?)
        // potential_parents has the banned nodes for node.
        // bannednodes is a filtered version of potential_parents, removing the nodes after in the ordering.
        // it only purose is to compute f_bar_z[node], i.e. f(Z)
        // std::set<int> bannednodes(ordering.begin(), ordering.begin() + position + 1);
        std::vector<int> parent_indices_banned_by_ordering; // Index in the banned_parents[node] vector.
        for (std::size_t j = 0; j < potential_parents[node].size(); j++)
        {
            if (std::find(ordering.begin(), ordering.begin() + position + 1, potential_parents[node][j]) != ordering.begin() + position + 1)
            // if (bannednodes.find( potential_parents[node][j]) != bannednodes.end())
            {
                parent_indices_banned_by_ordering.push_back(j); // This has inly ints. It for computing f(Z)
            }
        }

        // Compute f(Z) (the labelling), where Z is the parents of node, accoring toe the paper.
        // I.e. f_bar_z[node] = f(Pa(node))
        // if (numparents[node] == 0 || active_potential_parents_indices.size() == 0)
        if (potential_parents[node].size() == 0 || parent_indices_banned_by_ordering.size() == 0)
        {
            f_bar_z = 0; // all parents allowed f(null) = 0
        }
        else
        {
            int indextmp = 0;
            for (auto &item : parent_indices_banned_by_ordering)
            {
                indextmp += std::pow(2, item); // compute f(Z). add 1 since nodes are labeled from 0.
            }
            f_bar_z = rowmaps_backwards[node][indextmp]; // I think indextmp is the acutl f_bar_z.
        }
        return (f_bar_z);
    }
};

void swap_nodes(const int lower, const int upper, RightOrder &ro, OrderScoring &scoring)
{
    int node1 = ro.order[lower];
    int node2 = ro.order[upper];
    const auto &[node1_scoretmp, node2_scoretmp] = scoring.swap_nodes(lower, upper, ro.order, ro.node_scores);
    ro.order_score = ro.order_score - (ro.node_scores[node1] + ro.node_scores[node2]) + (node1_scoretmp + node2_scoretmp);
    ro.node_scores[node1] = node1_scoretmp;
    ro.node_scores[node2] = node2_scoretmp;
}

/**
 *  Asserts index_to < index_from.?
 */
inline void move_element(std::vector<int> &v, const std::size_t &index_from, const std::size_t &index_to);
inline void move_element(std::vector<int> &v, const std::size_t &index_from, const std::size_t &index_to)
{
    int element_to_move = v[index_from];
    v.erase(v.begin() + index_from);
    v.insert(v.begin() + index_to, element_to_move);
    return;
}

bool optimal_front(const RightOrder &ro,
                   size_t new_node,
                   OrderScoring &scoring)
{

    // std::cout << ro.inserted_max_order_scores[new_node] << " > " << ro.new_top_scores[new_node] << "?" << std::endl;
    if (definitelyGreaterThan(ro.inserted_max_order_scores[new_node], ro.new_top_scores[new_node], EPSILON))
    {
        // std::cout << "YES, so " << new_node << " not optimal at top " << std::endl;
        return (false);
    }
    // std::cout << "NO, so " << new_node << "  optimal at top " << std::endl;
    return (true);
}

bool optimal_front(const RightOrder &ro,
                   OrderScoring &scoring)
{
    RightOrder ro_tmp(ro);
    int p = ro.order.size();

    for (int i = ro.front_ind(); i < p - 1; ++i)
    {
        swap_nodes(i, i + 1, ro_tmp, scoring);
        // std::cout << ro_tmp << std::endl;
        // std::cout << ro_tmp.order_score << " > " << ro.order_score << "?" << std::endl;
        if (definitelyGreaterThan(ro_tmp.order_score, ro.order_score, EPSILON)) // this implies worse performance..
        // if (ro_tmp.order_score > ro.order_score)
        {
            // std::cout << "YES, so " << ro.front() << " not optimal at top " << std::endl;
            return (false);
        }
    }

    // size_t x = ro.order[p - ro.n];
    // if (approximatelyEqual(top_scores[x], ro.node_scores[x], EPSILON))
    // {
    //     std::cout << x << " indep of rest in " << ro << std::endl;
    // }
    // std::cout << std::endl;
    //  Maybe just swap back instead of copying ro.
    return (true);
}

bool independent_front(const RightOrder &ro,
                       const std::vector<double> &top_scores,
                       OrderScoring &scoring)
{
    size_t x = ro.front();
    return (approximatelyEqual(top_scores[x], ro.node_scores[x], EPSILON));
}

/**
 * Moves an unsscored node to the right of the sub order and scores.
 * Eg if x is the node to score:
 * ([0,x,2,3],5,6,4) -> ([0,2,3],x,5,6,4)
 *
 * NOTE: It's not really visible since n does note change..
 */
void make_visible(int from_index,
                  int to_index,
                  RightOrder &ro,
                  OrderScoring &scoring)
{
    int p = ro.order.size();
    int node = ro.order[from_index];
    int n = ro.n;

    if (to_index < p - n)
    {
        // If the new node is to be put somewhere in the order.
        // [0,(1),2,3|5,6,4] -> [(1),0,2,3|5,6,4]
        move_element(ro.order, from_index, to_index);                 // put it in the back [0,(1),2,3,|5,6,4] -> [(1),0,2,3,|,5,6,4]
        ro.node_scores[node] = scoring.score_pos(ro.order, to_index); // O(p)? This should also be quite time consuming..
        ro.order_score += ro.node_scores[node];
    }

    if (to_index >= p - n)
    {
        // If the new node is to be put on the very left of the sub order
        // ([0,x,2,3],5,6,4) -s> ([0,2,3],x,5,6,4)
        move_element(ro.order, from_index, p - n - 1);                 // put it in the back [0,(1),2,3,|5,6,4] -> [0,2,3,|(1),5,6,4]
        ro.node_scores[node] = scoring.score_pos(ro.order, p - n - 1); // O(p)? This should also be quite time consuming..
        ro.order_score += ro.node_scores[node];

        // If the new node is put somewhere in the sub order we wave to shift it in from the left.
        // [0,(1),2,3,|5,6,4] -> [0,2,3,|(1),5,6,4] -> [0,2,3,|5,(1),6,4]
        for (int i = p - n - 1; i < to_index; ++i)
        {
            swap_nodes(i, i + 1, ro, scoring);
        }
    }
    ro.n++;
}
double max_inserted_score(const RightOrder &ro,
                          size_t node_ind,
                          OrderScoring &scoring)
{

    // std::cout << ro.inserted_max_order_scores[new_node] << " > " << ro.new_top_scores[new_node] << "?" << std::endl;

    RightOrder ro_tmp(ro);
    make_visible(node_ind, ro.front_ind() - 1, ro_tmp, scoring);
    double maxsc = ro_tmp.order_score;
    for (size_t i = ro_tmp.front_ind(); i < ro.order.size() - 1; i++)
    {
        swap_nodes(i, i + 1, ro_tmp, scoring);

        if (definitelyGreaterThan(ro_tmp.order_score, maxsc, EPSILON))
        {
            // std::cout << "YES, so " << new_node << " not optimal at top " << std::endl;
            maxsc = ro_tmp.order_score;
        }
    }
    // std::cout << "NO, so " << new_node << "  optimal at top " << std::endl;
    return (maxsc);
}

bool equal_and_unordered_top(RightOrder &ro, int new_front, OrderScoring &scoring)
{
    size_t p = ro.order.size();

    if (ro.best_insert_pos[new_front] == ro.front_ind())
    {
        // If new_node is best inserted as the second top ([..i..],a,b,c) before, and best pos for i is s([...],a,i,b,c)
        if (approximatelyEqual(ro.new_top_scores[new_front], ro.inserted_max_order_scores[new_front], EPSILON))
        {
            // Compare  ([...], i,a,b,c) to ([...], a,i,b,c)
            if (new_front > ro.front())
            {
                // If new_node can be added as new front and has higher value, prune.
                return (true);
            }
        }
    }
    return (false);
}

bool equal_and_unordered_top(RightOrder &ro, OrderScoring &scoring)
{
    int p = scoring.numparents.size();
    int n = ro.n;
    int node1 = ro.order[p - n];
    int node2 = ro.order[p - n + 1];
    double node1_score = ro.node_scores[node1];
    double node2_score = ro.node_scores[node2];

    const auto &[node1_score_swap, node2_score_swap] = scoring.swap_nodes(p - n, p - n + 1, ro.order, ro.node_scores); // This should update the orderscore directly maybe.
    myswap(p - n, p - n + 1, ro.order);                                                                                // swap back

    if (approximatelyEqual(node1_score + node2_score, node1_score_swap + node2_score_swap, EPSILON))
    {
        if (node1 > node2)
        {
            return (true);
        }
    }
    return (false);
}

bool has_gap_new(RightOrder &ro,
                 const std::vector<double> &top_scores,
                 OrderScoring &scoring)
{
    std::size_t p = ro.order.size();
    std::size_t n = ro.n;

    for (std::size_t node_index = 0; node_index < p - n; node_index++)
    {
        int inserted_node = ro.order[node_index];
        double score_at_very_top = ro.order_score + top_scores[inserted_node];
        if (ro.best_insert_pos[inserted_node] == ro.front_ind())
        {
            // Here best isert score is at the second top ([...],a,i,b,c)
            if (approximatelyEqual(ro.inserted_max_order_scores[inserted_node], score_at_very_top, EPSILON))
            {
                if (equal_and_unordered_top(ro, inserted_node, scoring))
                {
                    return (true);
                }
            }

            if (definitelyGreaterThan(ro.inserted_max_order_scores[inserted_node], score_at_very_top, EPSILON))
            {
                // Here s((i,[...],a,b,c)) = "best inserted score"
                return (true);
            }
        }
        else
        {
            if (definitelyGreaterThan(ro.inserted_max_order_scores[inserted_node], score_at_very_top, EPSILON))
            {
                // Here s((i,[...],a,b,c)) = "best inserted score"
                return (true);
            }
        }
    }

    return (false);
}

bool has_gap(RightOrder &ro,
             const std::vector<double> &top_scores,
             OrderScoring &scoring)
{
    std::size_t p = ro.order.size();
    std::size_t n = ro.n;

    // Each node (position) has a value corresponding to the maximal score when that node is inserted somewhere in the sub order
    // std::vector<double> inserted_max_order_scores(p);

    std::vector<int> indep_nodes;

    // PrintVector(ro.order);
    for (std::size_t node_index = 0; node_index < p - n; node_index++)
    {
        // Go though all of the un used nodes and inject them at second place (pos k)
        // of the sub order and at the very front of the order.
        int x = ro.order[node_index];
        RightOrder ro_bkp(ro);
        // Let say we have the sub order ([1,2,3,x],4,5), so 4 is the top.
        // Order score when x at very left (x,[1,2,3],4,5)
        double score_at_very_top = ro.order_score + top_scores[x];

        // Score order when x is inserted everywhere.
        std::vector<double> injected_order_scores;

        size_t top = p - n;
        size_t new_top = top - 1;

        make_visible(node_index, new_top, ro, scoring); // insert new node as new top.
        injected_order_scores.push_back(ro.order_score);

        // double max_rest_score = injected_order_scores[0];

        for (size_t i = 0; i < n; i++)
        {
            // compute order score when x is inserted at pos p-n-1,...,p-1 called 0,...,n.
            swap_nodes(new_top + i, new_top + i + 1, ro, scoring);
            injected_order_scores.push_back(ro.order_score); // TODO: The max of this should maybe be stored for optimel_front?

            if (definitelyLessThan(score_at_very_top, ro.order_score, EPSILON))
            {
                if (i == 0)
                {
                    // Maybe this shold be in the for loop above??
                    // scoring ([1,2,3],4,x,5)
                    if (approximatelyEqual(score_at_very_top, injected_order_scores[1], EPSILON))
                    {
                        // if s(x,[1,2,3],4,5) == s([1,2,3],4,x,5)
                        // then check if s([1,2,3],x,4,5) == s([1,2,3],4,x,5)
                        // if so, prune if x > 4.
                        if (approximatelyEqual(injected_order_scores[0], injected_order_scores[1], EPSILON))
                        {
                            if (ro.order[new_top] < ro.order[top])
                            {
                                return (true);
                            }
                        }
                    }
                }
                else
                {
                    return (true); // Prune
                }
            }
        }

        // auto max_score = std::max_element(injected_order_scores.begin, injected_order_scores.end());
        // inserted_max_order_scores[ro.order[node_index]] = *max_score;
        ro = ro_bkp;
    }

    // ro.inserted_max_order_scores = inserted_max_order_scores;

    return (false);
}
RightOrder init_right_order(size_t node, OrderScoring &scoring)
{
    size_t p = scoring.numparents.size();
    std::vector<int> order(p, 0);
    // It all vectors with all the nodes.
    // The makes it easier to keep track of which nodes are used.
    for (size_t i = 0; i < p; i++)
    {
        order[i] = i;
    }
    myswap(node, p - 1, order);
    // score nodes
    std::vector<double> node_scores = scoring.score(order, p - 1, 1);
    //  score order
    double order_score = node_scores[node];
    RightOrder ro(order, order_score, node_scores, 1);

    // Set insert node scores...
    size_t top_ind = p - 1;
    size_t new_top_ind = top_ind - 1;
    double order_score_bkp = ro.order_score;
    for (size_t i = 0; i < top_ind; i++)
    {
        size_t inserted_node = ro.order[i];
        make_visible(i, new_top_ind, ro, scoring); // ([..i..],a,b,c) -> ([...],i,a,b,c)
        ro.new_top_scores[inserted_node] = ro.order_score;

        swap_nodes(new_top_ind, top_ind, ro, scoring); // ([...],i,a,b,c) -> ([...],a,i,b,c)
        ro.inserted_max_order_scores[inserted_node] = ro.order_score;
        ro.best_insert_pos[inserted_node] = p - 1;

        move_element(ro.order, top_ind, i);  // Move back
        ro.order_score = order_score_bkp;    // Restore order score
        ro.node_scores[inserted_node] = 0.0; // restore node score
        ro.node_scores[node] = order_score;
        ro.n--;
    }

    return (ro);
}

RightOrder add_node_in_front(const RightOrder &ro_prev, size_t index_of_el_to_insert, OrderScoring &scoring)
{
    RightOrder ro(ro_prev);
    ro.n = ro_prev.n + 1;
    size_t n = ro.n;
    size_t p = ro.order.size();
    size_t node = ro_prev.order[index_of_el_to_insert];
    // Order score is the old score plus the new node score.
    move_element(ro.order, index_of_el_to_insert, p - n);   // put it in the front
    double node_score = scoring.score_pos(ro.order, p - n); // O(p)? This should also be quite time consuming..
    // double node_score = ro_prev.new_top_scores[node];
    ro.node_scores[node] = node_score;
    ro.order_score = ro.order_score + ro.node_scores[node];
    double order_score_bkp = ro.order_score;

    // Should also update the inserted node order scores.
    // 1. Remove node from the list.
    ro.inserted_max_order_scores[node] = 0;
    ro.best_insert_pos[node] = 0;

    if (n == p)
    {
        return (ro);
    }

    // 2. Update the other scores by adding score of node when also the inserted node is a child.
    //    The missing one is ([...],node,inserted_node,a,b,c), and it is calculated in step 3.
    //    This takes is O(p) time.
    // std::cout << "\nAdding " << node << std::endl;
    // std::cout << "(BEFORE ADDING) max insert score of " << 19 << " in " << ro_prev << ": " << ro_prev.inserted_max_order_scores[19] << std::endl;
    // std::cout << "score of " << node << ": " << node_score << std::endl;
    // std::cout << "best insert pos " << ro.best_insert_pos[19] << std::endl;
    for (size_t i = 0; i < p - n; i++) // <= ?
    {
        // Note that inserted_node must be put as child of node when calculating score of node.
        size_t inserted_node = ro.order[i];
        move_element(ro.order, i, ro.front_ind());
        ro.inserted_max_order_scores[inserted_node] += scoring.score_pos(ro.order, ro.front_ind() - 1);
        move_element(ro.order, ro.front_ind(), i);
    }

    // std::cout << "(updated) max insert score of 19 : " << ro.inserted_max_order_scores[19] << std::endl;
    //  3. Go through all other nodes and compare max inserted score to inserted score
    //     between node and a in ([...],node,a,b,c).
    size_t top_ind = p - n;
    assert(top_ind > 0);

    size_t new_top_ind = top_ind - 1;

    for (size_t i = 0; i <= new_top_ind; i++)
    {
        size_t inserted_node = ro.order[i];
        make_visible(i, new_top_ind, ro, scoring);                                     // ([..i..],a,b,c) -> ([...],i,a,b,c).
        ro.new_top_scores[inserted_node] = ro.order_score;                             // set the new top score i.e. s([...],i,a,b,c)
        swap_nodes(new_top_ind, top_ind, ro, scoring);                                 // ([...],i,a,b,c) -> ([...],a,i,b,c). This changes ro.order_score.
        double inserted_max_order_score = ro.inserted_max_order_scores[inserted_node]; // this should contain max score of the other positions.

        if (approximatelyEqual(ro.order_score, inserted_max_order_score, EPSILON) ||
            definitelyGreaterThan(ro.order_score, inserted_max_order_score, EPSILON)
            // ro.order_score >= inserted_max_order_score
        )
        {
            ro.inserted_max_order_scores[inserted_node] = ro.order_score;
            ro.best_insert_pos[inserted_node] = top_ind;
        }

        move_element(ro.order, top_ind, i);  // Move back
        ro.order_score = order_score_bkp;    // Restore order score
        ro.node_scores[inserted_node] = 0.0; // restore node score
        ro.node_scores[node] = node_score;
        ro.n--;
    }
    // std::cout << "max insert score of " << 19 << " in " << ro << ": " << ro.inserted_max_order_scores[19] << std::endl;
    // std::cout << "best insert pos " << ro.best_insert_pos[19] << std::endl;

    return (ro);
}

std::vector<bool> order_to_boolvec(const std::vector<int> &order, size_t n, bool right_type)
{
    std::vector<bool> boolvec(order.size(), false);
    std::size_t p = order.size();

    if (right_type)
    {
        for (std::size_t i = p - n; i < p; i++)
        {
            boolvec[order[i]] = true;
        }
    }
    else
    {
        for (std::size_t i = 0; i < n; i++)
        {
            boolvec[order[i]] = true;
        }
    }
    return (boolvec);
}

/**
 * Row matrix of bool vectors, representing sets or forderings, and selects the
 * indices with highest scores.
 *
 * It does is coulmb by columns and adding the indices corresonding to zeros
 * and ones respectively as arrays. in the next itetration, column n+1 these indices are used.
 *
 */
std::vector<int> unique_sets(const std::vector<std::vector<bool>> &mats,
                             const std::vector<double> &order_scores)
{
    std::size_t N = mats.size();
    std::size_t p = mats[0].size();
    std::vector<int> samesets;
    std::vector<std::vector<int>> q;
    // The starting indecis are all indices. The first element indicates which column of
    // mats these indices should be used. To start, this is column 0.
    std::vector<int> start(N + 1, 0);
    for (std::size_t i = 0; i < N; i++)
    {
        start[i + 1] = i;
    }
    q.push_back(std::move(start));

    while (q.size() > 0)
    {
        // take a set of indeces from the queue
        std::vector<int> inds = q.back();
        std::size_t n = inds[0]; // ge the column to look at
        q.pop_back();            // remove from the queue

        std::vector<int> ones(1, n + 1);  // get the indices correspoing to ones. And indicate that thes should be used at column 1.
        std::vector<int> zeros(1, n + 1); // same for zeros.
        std::vector<int>::iterator i;
        for (i = inds.begin() + 1; i != inds.end(); ++i)
        {
            // For each index, check if 0 or 1 and add to corresponding vector.
            if (mats[*i][n] == 0)
            {
                zeros.push_back(*i); // Isnt this alreadu n+1 long?????
            }
            else
            {
                ones.push_back(*i);
            }
        }
        if (n == p - 1)
        { // Add to the unique elements
            if (zeros.size() > 1)
            {
                // PrintVector(zeros);
                //  Get the max scoring index
                double maxscore = -INFINITY;
                int max_ind = -1;
                std::vector<int>::iterator i;
                for (i = zeros.begin() + 1; i != zeros.end(); ++i)
                {
                    // std::cout << "row " << *i << " score " << order_scores[*i] << std::endl;
                    if (order_scores[*i] > maxscore)
                    {
                        maxscore = order_scores[*i];
                        max_ind = *i;
                    }
                }
                samesets.push_back(max_ind);
            }

            if (ones.size() > 1)
            {
                // PrintVector(ones);
                //  Get the max scoring index
                double maxscore = -INFINITY;
                int max_ind = -1;
                std::vector<int>::iterator i;
                for (i = ones.begin() + 1; i != ones.end(); ++i)
                {
                    // std::cout << "row " << *i << " score " << order_scores[*i] << std::endl;
                    if (order_scores[*i] > maxscore)
                    {
                        maxscore = order_scores[*i];
                        max_ind = *i;
                    }
                }

                samesets.push_back(max_ind);
            }
        }
        else
        { // Add elements to the queue
            if (zeros.size() > 1)
                q.push_back(std::move(zeros));
            if (ones.size() > 1)
                q.push_back(std::move(ones));
        }
    }

    return (samesets);
}

std::vector<RightOrder> prune_equal_sets(std::vector<RightOrder> &right_orders,
                                         bool right_type)
{
    std::vector<std::vector<bool>> boolmat;
    std::vector<double> order_scores;

    // std::cout << "Creating boolmatrix" << std::endl;
    for (const RightOrder &ro : right_orders)
    {
        std::vector<bool> boolvec = order_to_boolvec(ro.order, ro.n, right_type);
        boolmat.push_back(std::move(boolvec));
        order_scores.push_back(ro.order_score);
    }

    // std::cout << "Get indices of unique maximal scoring sets" << std::endl;
    std::vector<int> pruned_inds = unique_sets(boolmat, order_scores);
    std::vector<double> pruned_scores;
    std::vector<RightOrder> kept_ros;

    for (const auto &ind : pruned_inds)
    {
        kept_ros.push_back(right_orders[ind]);
    }

    return (kept_ros);
}

std::vector<RightOrder> prune_indep_front(std::vector<RightOrder> &potential_orders,
                                          std::vector<double> &top_scores,
                                          OrderScoring &scoring)
{
    // get all vectors with inpenendent front.
    // keep the one with lowest index. - Highest??  / Felix
    int max_indep_node = 0;
    bool has_indep_node = false;
    size_t max_indep_ind = 0;
    size_t jj = 0;
    for (RightOrder &ro : potential_orders)
    {
        if (independent_front(ro, top_scores, scoring))
        {
            // std::cout << ro.front() << " indep of hidden nodes in " << ro << std::endl;
            if (ro.front() >= max_indep_node)
            {
                // std::cout << ro.front() << " >= " << max_indep_node << std::endl;
                has_indep_node = true;
                max_indep_node = ro.front();
                max_indep_ind = jj;
            }
        }
        jj++;
    }

    if (!has_indep_node)
    {
        return (potential_orders);
    }
    else
    {
        std::vector<RightOrder> ret = {potential_orders[max_indep_ind]};
        return (ret);
    }
}

std::tuple<std::vector<int>, double, size_t, size_t> sequential_opt(OrderScoring &scoring);
std::tuple<std::vector<int>, double, size_t, size_t> sequential_opt(OrderScoring &scoring)
{
    std::size_t p = scoring.numparents.size();
    std::cout << "Starting optimization" << std::endl;
    std::vector<RightOrder> right_orders;
    std::vector<RightOrder> right_orders_prev;

    /**
     * Compute S((x,[...])) for all x.
     */
    std::vector<double> top_scores(p);
    std::vector<int> order_tmp(p);

    for (size_t i = 0; i < p; i++)
    {
        order_tmp[i] = i;
    }

    for (size_t i = 0; i < p; i++)
    {
        move_element(order_tmp, i, 0);
        top_scores[i] = scoring.score_pos(order_tmp, 0);
        move_element(order_tmp, 0, i);
    }

    size_t orders1 = 0;
    size_t orders2 = 0;
    size_t orders3 = 0;
    size_t max_n_particles = 0;
    size_t tot_n_particles = 0;
    /**
     * Start build and prune
     */
    for (std::size_t n = 1; n <= p; n++)
    {

        if (n == 1)
        {
            std::vector<RightOrder> potential_orders;
            // std::cout << "Adding nodes to empty order" << std::endl;
            for (size_t node_index = 0; node_index <= p - n; node_index++)
            {

                RightOrder ro = init_right_order(node_index, scoring);
                if (!optimal_front(ro, scoring))
                {
                    continue;
                }
                ++orders1;

                if (has_gap_new(ro, top_scores, scoring))
                {
                    // assert(has_gap(ro, top_scores, scoring));
                    continue;
                }
                else
                {
                    // assert(!has_gap(ro, top_scores, scoring));
                }

                potential_orders.push_back(std::move(ro));
            }

            potential_orders = prune_indep_front(potential_orders, top_scores, scoring);
            right_orders = potential_orders;
            orders2 = right_orders.size();
            orders3 = orders2;
        }
        else
        {
            for (RightOrder &prev_order : right_orders_prev)
            {

                std::vector<RightOrder> potential_orders;
                for (size_t node_index = 0; node_index <= p - n; node_index++)
                {
                    int new_node = prev_order.order[node_index];

                    if (!optimal_front(prev_order, new_node, scoring))
                    {
                        continue;
                    }

                    if (!independent_front(prev_order, top_scores, scoring))
                    {
                        if (equal_and_unordered_top(prev_order, new_node, scoring))
                        {
                            continue;
                        }
                    }
                    RightOrder ro = add_node_in_front(prev_order, node_index, scoring);

                    potential_orders.push_back(std::move(ro));
                }
                // std::cout << "#potential orders" << potential_orders.size() << std::endl;
                potential_orders = prune_indep_front(potential_orders, top_scores, scoring);
                right_orders.insert(right_orders.end(), potential_orders.begin(), potential_orders.end());
            }

            // O(|right_orders_prev| * p) particles. space
            orders1 = right_orders.size();
            // std::cout << "Prune equal sets: #" << orders1 << std::endl;
            right_orders = prune_equal_sets(right_orders, true); // O(#particles * p) space and time
            orders2 = right_orders.size();

            // Remove if has gaps
            std::vector<RightOrder> right_orders_tmp;
            for (RightOrder &ro : right_orders)
            {
                if (has_gap_new(ro, top_scores, scoring))
                {
                    // assert(has_gap(ro, top_scores, scoring));
                    continue;
                }
                else
                {
                    // assert(!has_gap_new(ro, top_scores, scoring));
                }

                right_orders_tmp.push_back(std::move(ro));
            }
            right_orders = std::move(right_orders_tmp);
            orders3 = right_orders.size();
            // std::cout << "# orders: " << orders3 << std::endl;
        }

        // std::cout << "orders of size " << n << std::endl;
        // for (const RightOrder &opt_order : right_orders)
        // {
        //     if (n <= 3 || n >= 18)
        //     {
        //         std::cout << opt_order << ": " << opt_order.order_score << std::endl;
        //         // std::cout << "score: " << opt_order.order_score << std::endl;
        //     }
        // }

        // Print some statistics //
        /* Add to number of particles sum */
        tot_n_particles += orders3;
        /* Check that scores are correct */
        if (orders3 > max_n_particles)
        {
            max_n_particles = orders3;
        }
        
        auto max_ro = std::max_element(right_orders.begin(),
                                       right_orders.end(),
                                       [](const RightOrder &a, const RightOrder &b)
                                       { return a.order_score < b.order_score; });

        // std::cout << "max scoring sub order " << std::endl;
        // std::cout << *max_ro << std::endl;
        // std::cout << "score: " << max_ro->order_score << std::endl;

        // check correct score
        std::vector<double> sc = scoring.score(max_ro->order, p - n, n); // Take only the last n elements in the vector
        // PrintVector(sc);
        double max_score_check = std::accumulate(sc.begin(), sc.end(), 0.0);
        // std::cout << "correct max score: " << max_score_check << std::endl;
        assert(std::abs(max_ro->order_score - max_score_check) < EPSILON);

        std::cout << n << " & " << orders1 << " & " << orders2 << " & " << orders3 << " & " << max_ro->order_score << " \\\\" << std::endl;
        right_orders_prev = std::move(right_orders);
        // std::cout << "after move" << std::endl;
    }

    auto max_ro = std::max_element(right_orders_prev.begin(),
                                   right_orders_prev.end(),
                                   [](const RightOrder &a, const RightOrder &b)
                                   { return a.order_score < b.order_score; });

    std::cout << "MAX order " << *max_ro << std::endl;

    return (std::make_tuple(max_ro->order, max_ro->order_score, max_n_particles, tot_n_particles));
}

OrderScoring get_score(Rcpp::List ret)
{

    // Read MAP flag
    bool MAP = Rcpp::as<int>(ret["MAP"]);
    std::cout << "MAP:" << MAP << std::endl;

    // Read numparents
    std::vector<int> numparents = Rcpp::as<std::vector<int>>(ret["numparents"]);

    // Read scoretable
    std::vector<std::vector<std::vector<double>>> scoretable = Rcpp::as<std::vector<std::vector<std::vector<double>>>>(ret["scoretable"]);

    // Read parent table
    Rcpp::List parenttableR = Rcpp::as<Rcpp::List>(ret["parenttable"]);
    std::size_t p = parenttableR.size();
    std::vector<Rcpp::IntegerMatrix> parenttable;
    for (std::size_t i = 0; i < p; i++)
    {
        Rcpp::IntegerMatrix m = Rcpp::as<Rcpp::IntegerMatrix>(parenttableR[i]);
        parenttable.push_back(m);
    }

    // Read banned score
    Rcpp::List bannedscoreR = Rcpp::as<Rcpp::List>(ret["bannedscore"]);
    std::vector<std::vector<std::vector<double>>> bannedscore;
    for (std::size_t i = 0; i < p; i++)
    {
        Rcpp::NumericMatrix m = Rcpp::as<Rcpp::NumericMatrix>(bannedscoreR[i]);
        std::vector<std::vector<double>> mat(m.rows(), std::vector<double>(m.cols()));
        for (int j = 0; j < m.rows(); j++)
        {
            for (int k = 0; k < m.cols(); k++)
            {
                mat[j][k] = m(j, k);
            }
        }
        bannedscore.push_back(mat);
    }

    // Read rowmaps_backwards
    Rcpp::List rowmaps_backwardsR = Rcpp::as<Rcpp::List>(ret["rowmaps_backwards"]);
    std::vector<Rcpp::IntegerVector> rowmaps_backwards;
    for (std::size_t i = 0; i < p; i++)
    {
        Rcpp::IntegerVector m = Rcpp::as<Rcpp::IntegerVector>(rowmaps_backwardsR[i]);
        rowmaps_backwards.push_back(m);
    }

    // Read aliases
    Rcpp::List aliasesR = Rcpp::as<Rcpp::List>(ret["aliases"]);
    std::vector<std::vector<int>> aliases;
    for (std::size_t i = 0; i < p; i++)
    {
        std::vector<int> m = Rcpp::as<std::vector<int>>(aliasesR[i]);
        aliases.push_back(m);
    }

    // Read plus1listsparents
    Rcpp::List plus1listsparentsR = Rcpp::as<Rcpp::List>(ret["plus1listsparents"]);
    std::vector<std::vector<int>> plus1listsparents;
    for (std::size_t i = 0; i < p; i++)
    {
        Rcpp::IntegerVector m = Rcpp::as<Rcpp::IntegerVector>(plus1listsparentsR[i]);
        std::vector<int> tmp;
        for (auto e : m)
        {
            tmp.push_back(e);
        }
        plus1listsparents.push_back(tmp);
    }

    std::vector<int> scorepositions(p);
    for (std::size_t i = 0; i < p; ++i)
    {
        scorepositions[i] = i;
    }

    OrderScoring scoring(aliases,
                         numparents,
                         rowmaps_backwards,
                         plus1listsparents,
                         scoretable,
                         bannedscore,
                         MAP);

    return (scoring);
}

// [[Rcpp::plugins(cpp17)]]

// [[Rcpp::export]]
Rcpp::List r_sequential_opt(Rcpp::List ret)
{
    OrderScoring scoring = get_score(ret);

    const auto &[order, log_score, max_n_particles, tot_n_particles] = sequential_opt(scoring);

    Rcpp::List L = Rcpp::List::create(Rcpp::Named("order") = order,
                                      Rcpp::Named("log_score") = log_score,
                                      Rcpp::Named("max_n_particles") = max_n_particles,
                                      Rcpp::Named("tot_n_particles") = tot_n_particles);

    return (L);
}