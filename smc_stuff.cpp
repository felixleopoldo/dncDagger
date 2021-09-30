#include <bits/stdc++.h>
#include <RInside.h>
#include <cassert>
#include <chrono>
#include <thread>
#include <iostream>
#include "thread_pool.hpp"

using namespace std::chrono;
typedef std::tuple<std::vector<int>, int, int> cache_keytype;
typedef std::pair<std::vector<int>, std::set<int>> cache_keytype2;
typedef std::pair<std::vector<int>, int> cache_keytype3;

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
void PrintVector(std::vector<T> &arr)
{
    copy(arr.begin(), arr.end(), std::ostream_iterator<T>(std::cout, ","));
    std::cout << std::endl;
}

template <typename T>
void PrintVector(const std::vector<T> &arr)
{
    copy(arr.begin(), arr.end(),
         std::ostream_iterator<T>(std::cout, ","));
    std::cout << std::endl;
}

template <typename T>
void PrintVector(const std::valarray<T> &arr)
{
    std::for_each(std::begin(arr), std::end(arr), [](T c)
                  { std::cout << c << ','; });
    std::cout << std::endl;
}

class OrderScoring
{
private:
    std::vector<std::vector<int>> potential_parents;
    std::vector<Rcpp::IntegerVector> rowmaps_backwards;
    std::vector<std::vector<int>> potential_plus1_parents;
    std::vector<std::vector<std::vector<double>>> scoretable;
    std::vector<std::vector<std::vector<double>>> scoresmatrices;
    bool MAP;

public:
    std::vector<int> numparents;
    std::map<cache_keytype3, std::vector<double>> cache;
    int cache_hits;
    OrderScoring(
        std::vector<std::vector<int>> potential_parents,
        std::vector<int> numparents,
        std::vector<Rcpp::IntegerVector> rowmaps_backwards,
        std::vector<std::vector<int>> potential_plus1_parents,
        std::vector<std::vector<std::vector<double>>> scoretable,
        std::vector<std::vector<std::vector<double>>> scoresmatrices,
        bool MAP,
        std::map<cache_keytype3, std::vector<double>> cache) : potential_parents(potential_parents),
                                                               rowmaps_backwards(rowmaps_backwards),
                                                               potential_plus1_parents(potential_plus1_parents),
                                                               scoretable(scoretable),
                                                               scoresmatrices(scoresmatrices),
                                                               MAP(MAP),
                                                               numparents(numparents),
                                                               cache(cache),
                                                               cache_hits(0)
    {
    }

    /**
     * Re-calculation scores after swapping up node_a so that (node_a, node_b) --> (node_b, node_a).
     * 
     * Order of node1 is lower than node2.
     * 
     */
    std::tuple<double, double> swap_nodes(int nodea_index, int nodeb_index,
                                          std::vector<int> &ordering,
                                          const std::vector<double> &node_scores)
    {
        int node_a = ordering[nodea_index];
        int node_b = ordering[nodeb_index];
        double nodea_score = 0.0; // just to initialize
        double nodeb_score = 0.0; // just to initialize
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
                //std::cout << nodeb_as_plus1_score - node_scores[node_a] << " " << std::endl;
                if (std::abs(nodeb_as_plus1_score - node_scores[node_a]) > 0.000001)
                { // 0.000001 is arbitrary
                    // OK
                    //std::cout << "NO RECOMPUTE order score for node" << std::endl;
                    nodea_score = std::log(std::exp(node_scores[node_a] - max) - std::exp(nodeb_as_plus1_score - max)) + max; // gets inf... 0 at node_scores[node_a] but something at node_scores[node_b]
                }
                else
                {
                    // round off error. Recompute.
                    //std::cout << "RECOMPUTE order score for node" << std::endl;
                    myswap(nodea_index, nodeb_index, ordering);
                    nodea_score = score_pos(ordering, nodeb_index);
                    myswap(nodea_index, nodeb_index, ordering);
                }

                //std::cout << "true score " << true_score << " calcuated score " << node_scores[node_a] << std::endl;
                //assert(std::abs(nodea_score-true_score) < 0.001);
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
                //std::swap(ordering[nodea_index], ordering[nodeb_index]);
                myswap(nodea_index, nodeb_index, ordering);
                int f_bar_z = get_f_bar_z(nodea_index, ordering);
                double nodea_as_plus1_score = scoresmatrices[node_b][f_bar_z][plus1_parents_index + 1];
                //std::swap(ordering[nodea_index], ordering[nodeb_index]); // Swap back.
                myswap(nodea_index, nodeb_index, ordering);
                double max = std::max(node_scores[node_b], nodea_as_plus1_score);

                // Check that
                //if (std::abs(nodea_as_plus1_score - node_scores[node_b]) > 0.000001)
                //{ // -10 is arbitrary
                nodeb_score = std::log(std::exp(node_scores[node_b] - max) + std::exp(nodea_as_plus1_score - max)) + max;
                //}
                // else
                // {
                //     std::swap(ordering[nodea_index], ordering[nodeb_index]);
                //     nodeb_score = score_pos(ordering, nodea_index);
                //     std::swap(ordering[nodea_index], ordering[nodeb_index]);
                // }
            }
        }

        //std::swap(ordering[nodea_index], ordering[nodeb_index]); // last swap
        myswap(nodea_index, nodeb_index, ordering);

        return (std::make_tuple(nodea_score, nodeb_score));
    }

    /**
    * Score elements in scorenodes from scorepositions and n_elemens on.
    */
    // If in cache
    //PrintVector(ordering);
    //std::vector<int> nodeskey(ordering.begin() + from_orderpos, ordering.begin() + from_orderpos + n_elements);

    //PrintVector(nodeskey);

    //if(ordering.begin()+1+from_orderpos + n_elements == ordering.end() +1 ){
    //    std::cout << "end" << std::endl;
    //} else {

    //std::set<int> restkey(ordering.begin() + from_orderpos + n_elements, ordering.end());
    //std::set<vector> restkey(ordering.begin() + from_orderpos + n_elements, ordering.end());
    //std::vector<int> tmp(restkey.begin(), restkey.end());
    //PrintVector(tmp);

    //}

    ///std::cout << std::endl;
    //cache_keytype2 mykey = std::make_pair(nodeskey, restkey);

    // std::vector<int> nodeskey3(ordering.begin() + from_orderpos, ordering.end());
    // cache_keytype3 mykey = std::make_pair(nodeskey3, n_elements);
    // if (cache.count(mykey))
    // {
    //     cache_hits++;
    //     return (cache[mykey]);
    // }

    // cache should be the actual nodes in a list and the rst of the nodes in a set.
    //cache_keytype suborder = std::make_tuple(ordering, from_orderpos, n_elements);
    //if (cache.count(suborder))
    //{
    //std::cout << "cache hit" << std::endl;
    //    return (cache[suborder]);
    //}
    std::vector<double> *score(const std::vector<int> &ordering, const std::size_t &from_orderpos, const std::size_t &n_elements) const
    {
        std::size_t n = ordering.size();
        std::vector<double> *orderscores = new std::vector<double>(n, 0.0); // O(p)           // orderscores <- vector("double", n)
        std::vector<int> active_plus1_parents_indices;                      // active_plus1_parents_indices < -vector("list", n)
        int f_bar_z;

        for (std::size_t position = from_orderpos; position < from_orderpos + n_elements; ++position)
        {
            int node = ordering[position];
            if (position == n - 1)
            {
                // no parents allowed, i.e.only first row, only first list
                (*orderscores)[node] = scoretable[node][0][0]; // orderscores[node] <- scoretable [[node]][[1]][1, 1]
                //active_plus1_parents_indices[node].assign({0});           // active_plus1_parents_indices [[node]] <- c(1)
                //f_bar_z[node] = std::pow(2, potential_parents[node].size()); // f_bar_z[node] <- c(2 ^ numparents[node])
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
                    (*orderscores)[node] = *std::max_element(plus1_parents_scores.begin(), plus1_parents_scores.end());
                }
                else
                {
                    (*orderscores)[node] = sum_log_probs(plus1_parents_scores);
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
        //std::set<int> bannednodes(ordering.begin(), ordering.begin() + position + 1);
        std::vector<int> parent_indices_banned_by_ordering; // Index in the banned_parents[node] vector.
        for (std::size_t j = 0; j < potential_parents[node].size(); j++)
        {
            if (std::find(ordering.begin(), ordering.begin() + position + 1, potential_parents[node][j]) != ordering.begin() + position + 1)
            //if (bannednodes.find( potential_parents[node][j]) != bannednodes.end())
            {
                parent_indices_banned_by_ordering.push_back(j); // This has inly ints. It for computing f(Z)
            }
        }

        // Compute f(Z) (the labelling), where Z is the parents of node, accoring toe the paper.
        // I.e. f_bar_z[node] = f(Pa(node))
        //if (numparents[node] == 0 || active_potential_parents_indices.size() == 0)
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

std::pair<std::vector<int>, double> mh_swap_move(std::vector<int> &ordering,
                                                 std::vector<double> &node_scores,
                                                 const double &order_score,
                                                 OrderScoring &scoring)
{
    std::size_t p = ordering.size();
    int i = rand_int_in_range(1, p - 1);
    int nodea_index = i - 1;
    int nodeb_index = i;

    int node1 = ordering[i - 1];
    int node2 = ordering[i];

    double cur_node1_score = node_scores[i - 1];
    double cur_node2_score = node_scores[i];
    double cur_order_score = order_score;

    const auto &[node1_scoretmp, node2_scoretmp] = scoring.swap_nodes(i - 1, i, ordering, node_scores); // This should update the orderscore directly maybe.
    double prop_order_score = cur_order_score - (node_scores[node1] + node_scores[node2]) + (node1_scoretmp + node2_scoretmp);

    //std::cout << prop_order_score << std::endl;
    //std::cout << cur_order_score << std::endl;

    double ratio = std::exp(prop_order_score - cur_order_score);
    double alpha = std::min(1.0, ratio);
    double a = (float)std::rand() / RAND_MAX;

    //std::cout << "rand number: " << a << std::endl;
    //std::cout << "accept prob: " << alpha << std::endl;
    if (a < alpha)
    {
        //std::cout << "Accept" << std::endl;
        node_scores[node1] = node1_scoretmp;
        node_scores[node2] = node2_scoretmp;
        return (std::make_pair(ordering, prop_order_score));
    }
    else
    {
        //std::cout << "Reject" << std::endl;
        // swap back
        myswap(i - 1, i, ordering);
        return (std::make_pair(ordering, cur_order_score));
    }
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

inline std::vector<double> *dist_from_logprobs(const std::vector<double> &log_probs)
{
    std::size_t n = log_probs.size();
    double max = *std::max_element(log_probs.begin(), log_probs.end());
    std::vector<double> probs_rescaled(n);
    //std::vector<double> *probs_rescaled = new std::vector<double>(n);
    std::vector<double> *norm_prob = new std::vector<double>(n);
    double sum_exp = 0.0;
    for (std::size_t i = 0; i != n; ++i)
    {
        (probs_rescaled)[i] = std::exp(log_probs[i] - max);
        sum_exp += (probs_rescaled)[i];
    }

    for (std::size_t i = 0; i != n; ++i)
    {
        (*norm_prob)[i] = probs_rescaled[i] / sum_exp;
    }

    return (norm_prob);
    //return (probs_rescaled);
}

/**
 * Puts a random node at index 0 of the ordering. This means that its score is calulcated
 * based on all the other nodes as potential parents.
 */

double random_parent_init(OrderScoring &scoring,
                          std::vector<int> &output_order,
                          std::vector<double> &output_node_scores,
                          double &output_order_score,
                          std::default_random_engine &generator)
{
    int n = 0;
    int p = output_order.size();
    int index_of_el_to_insert = rand_int_in_range(n, p - 1); // Draw one of the remaining nodes

    if (index_of_el_to_insert != 0)
    {
        move_element(output_order, index_of_el_to_insert, n); // Move to node_index to index n=0
    }
    std::vector<double> *sc = scoring.score(output_order, n, 1); // OK, score only index 0
    output_node_scores = *sc;
    delete sc;
    output_order_score = output_node_scores[output_order[n]]; // std::accumulate(log_node_scores[i].begin(), log_node_scores[i].end(), 0.0)
    return (-std::log(p - n));
}

/**
 * Proposal which adds a new node to an exisitng ordering at a certain position 
 * with probability proportional to the new order. I.e. Score all orderings 
 * with the new node in any position then draw one of those orderings with 
 * corresonding relative probability.
 */

double neig_proportional_proposal(OrderScoring &scoring,
                                  const std::vector<int> &input_order,
                                  const std::vector<double> &nodes_scores,
                                  const double &input_order_score,
                                  std::vector<int> &output_order,
                                  std::vector<double> &output_node_scores,
                                  double &output_order_score,
                                  const int &first_n_elements,
                                  const int &index_of_el_to_insert,
                                  std::default_random_engine &generator);

double neig_proportional_proposal(OrderScoring &scoring,
                                  const std::vector<int> &input_order,
                                  const std::vector<double> &input_node_scores,
                                  const double &input_order_score,
                                  std::vector<int> &output_order,
                                  std::vector<double> &output_node_scores,
                                  double &output_order_score,
                                  const int &first_n_elements,

                                  std::default_random_engine &generator)
{
    int p = input_order.size();
    int index_of_el_to_insert = rand_int_in_range(first_n_elements, p - 1); // Draw one of the remaining nodes

    std::vector<double> order_scores(first_n_elements + 1);
    std::vector<double> node_scores(input_node_scores);
    std::vector<int> order(input_order);
    double order_score(input_order_score);
    int node = order[index_of_el_to_insert];

    // Calculate score for node at index first_n_elements. The other nodes have the
    // same score sine those are still befor in the ordering.
    move_element(order, index_of_el_to_insert, first_n_elements);   // put it in the back
    node_scores[node] = scoring.score_pos(order, first_n_elements); // O(p)? This should also be quite time consuming.
    order_scores[first_n_elements] = order_score + node_scores[node];

    for (int i = first_n_elements; i > 0; --i) // put it in the back instead
    {
        int node1 = order[i - 1];
        int node2 = order[i];
        const auto &[node1_scoretmp, node2_scoretmp] = scoring.swap_nodes(i - 1, i, order, node_scores); // This should update the orderscore directly maybe.
        order_scores[i - 1] = order_scores[i] - (node_scores[node1] + node_scores[node2]) + (node1_scoretmp + node2_scoretmp);
        node_scores[node1] = node1_scoretmp;
        node_scores[node2] = node2_scoretmp;
    }

    // Sample the new ordering from the order_scores distribution, where the index
    // specifies the position were node is inserted.
    std::vector<double> *neig_dist = dist_from_logprobs(order_scores);                  // Sample one order in the neighborhood
    std::discrete_distribution<int> distribution(neig_dist->begin(), neig_dist->end()); // Create a distribution. O(p+1)
    int insert_pos = distribution(generator);

    // Loop to the to get the correct node scorings. Ths is the fastest way I can think of right now.
    // Caching the node scores is O(p) for each iteration.
    output_order = input_order; // copy
    output_order_score = order_scores[insert_pos];
    output_node_scores = input_node_scores;                              // copy
    move_element(output_order, index_of_el_to_insert, first_n_elements); // put it in the back
    output_node_scores[node] = scoring.score_pos(output_order, first_n_elements);

    for (int i = first_n_elements; i > insert_pos; --i)
    {
        int node1 = output_order[i - 1];
        int node2 = output_order[i];
        const auto &[node1_scoretmp, node2_scoretmp] = scoring.swap_nodes(i - 1, i, output_order, output_node_scores);
        output_node_scores[node1] = node1_scoretmp;
        output_node_scores[node2] = node2_scoretmp;
    }

    double log_prob = std::log((*neig_dist)[insert_pos]);
    delete neig_dist;
    return log_prob - std::log(p - first_n_elements);
}

double max_at_end_proposal(OrderScoring &scoring,
                           const std::vector<int> &input_order,
                           const std::vector<double> &nodes_scores,
                           const double &input_order_score,
                           std::vector<int> &output_order,
                           std::vector<double> &output_node_scores,
                           double &output_order_score,
                           const int &first_n_elements,
                           const int &index_of_el_to_insert,
                           std::default_random_engine &generator);

double max_at_end_proposal(OrderScoring &scoring,
                           const std::vector<int> &input_order,
                           const std::vector<double> &input_node_scores,
                           const double &input_order_score,
                           std::vector<int> &output_order,
                           std::vector<double> &output_node_scores,
                           double &output_order_score,
                           const int &first_n_elements,
                           const int &index_of_el_to_insert,
                           std::default_random_engine &generator)
{
    std::vector<double> order_scores(first_n_elements + 1);
    std::vector<double> node_scores(input_node_scores);
    std::vector<int> order(input_order);
    int p = order.size();
    double order_score(input_order_score);
    int node = order[index_of_el_to_insert];

    // Calculate score for node at` index first_n_elements. The other nodes have the
    // same score sine those are still befor in the ordering.
    move_element(order, index_of_el_to_insert, p - 1);   // put it in the back
    node_scores[node] = scoring.score_pos(order, p - 1); // O(p)? This should also be quite time consuming..
    order_scores[first_n_elements] = order_score + node_scores[node];

    for (int i = p - 1; i > p - first_n_elements; --i) // put it in the back instead
    {
        int node1 = order[i - 1];
        int node2 = order[i];
        const auto &[node1_scoretmp, node2_scoretmp] = scoring.swap_nodes(i - 1, i, order, node_scores); // This should update the orderscore directly maybe.
        order_scores[i - 1] = order_scores[i] - (node_scores[node1] + node_scores[node2]) + (node1_scoretmp + node2_scoretmp);
        node_scores[node1] = node1_scoretmp;
        node_scores[node2] = node2_scoretmp;
    }

    // Take the ordering with the highest prob. Or only it is in the front?

    // Sample the new ordering from the order_scores distribution, where the index
    // specifies the position were node is inserted.
    std::vector<double> *neig_dist = dist_from_logprobs(order_scores);                  // Sample one order in the neighborhood
    std::discrete_distribution<int> distribution(neig_dist->begin(), neig_dist->end()); // Create a distribution. O(p+1)
    int insert_pos = distribution(generator);

    // Loop to the to get the correct node scorings. This is the fastest way I can think of right now.
    // Caching the node scores is O(p) for each iteration.
    output_order = input_order; // copy
    output_order_score = order_scores[insert_pos];
    output_node_scores = input_node_scores;                              // copy
    move_element(output_order, index_of_el_to_insert, first_n_elements); // put it in the back
    output_node_scores[node] = scoring.score_pos(output_order, first_n_elements);

    for (int i = first_n_elements; i > insert_pos; --i)
    {
        int node1 = output_order[i - 1];
        int node2 = output_order[i];
        const auto &[node1_scoretmp, node2_scoretmp] = scoring.swap_nodes(i - 1, i, output_order, output_node_scores);
        output_node_scores[node1] = node1_scoretmp;
        output_node_scores[node2] = node2_scoretmp;
    }

    double log_prob = std::log((*neig_dist)[insert_pos]);
    delete neig_dist;
    return log_prob;
}

std::tuple<std::vector<int>, double, std::vector<double>, int, bool> put_node_in_front(const std::vector<int> &input_order,
                                                                                       int n,
                                                                                       int index_of_el_to_insert,
                                                                                       const std::vector<double> &input_node_scores,
                                                                                       const double &input_order_score,
                                                                                       const int prev_order_number,
                                                                                       OrderScoring &scoring)
{
    //std::cout << "n=" << n << std::endl;
    // Try to p std::vector<double> order_scores(first_n_elements + 1);
    std::vector<double> node_scores(input_node_scores);
    std::vector<int> order(input_order);
    int p = order.size();
    //double order_score(input_order_score);
    int node = order[index_of_el_to_insert];

    int order_number = prev_order_number + node; // TODO: What is this? Seems strange..

    // Calculate score for node at index first_n_elements. The other nodes have the
    // same score sine those are still befor in the ordering.
    move_element(order, index_of_el_to_insert, p - n); // put it in the back
    std::vector<int> order_ret(order);

    double new_node_score = scoring.score_pos(order, p - n); // O(p)? This should also be quite time consuming..

    node_scores[node] = new_node_score;

    std::vector<double> order_scores(p);
    double order_score = input_order_score + node_scores[node];
    double max_score = order_score;
    bool optimal_at_front = true;

    //std::cout << p - n << ": " << order_score << ", ";
    for (int i = p - n; i < p - 1; ++i) // put it in the back instead?
    {
        int node1 = order[i];
        int node2 = order[i + 1];
        const auto &[node1_scoretmp, node2_scoretmp] = scoring.swap_nodes(i, i + 1, order, node_scores); // This should update the order score directly maybe.
        //std::cout << node2_scoretmp << ", " << node1_scoretmp << std::endl;
        order_score = order_score - (node_scores[node1] + node_scores[node2]) + (node1_scoretmp + node2_scoretmp);
        node_scores[node1] = node1_scoretmp;
        node_scores[node2] = node2_scoretmp;

        //std::cout << i + 1 << ": " << order_score << ", ";
        if (order_score > max_score)
        {
            //std::cout << "Better score at pos " << i + 1 << ", (" << order_score << ") instead of pos " << p - n << ". Break." << std::endl;
            optimal_at_front = false;
            break; // Just break the loop if some better osition if found.
        }
    }
    //std::cout << std::endl;
    //PrintVector(order);
    //PrintVector(order_ret);

    std::vector<double> node_scores_ret(input_node_scores);
    node_scores_ret[node] = new_node_score;
    //PrintVector(order_scores);
    // If not optimal at the front (p-n), nodes_scores is not correct
    return (std::make_tuple(order_ret, max_score, node_scores_ret, order_number, optimal_at_front));
}

void swap_nodes(const int lower, const int upper, std::vector<int> &order, std::vector<double> &node_scores, double &order_score, OrderScoring &scoring)
{
    int node1 = order[lower];
    int node2 = order[upper];
    const auto &[node1_scoretmp, node2_scoretmp] = scoring.swap_nodes(lower, upper, order, node_scores);
    //std::cout << node2_scoretmp << ", " << node1_scoretmp << std::endl;
    order_score = order_score - (node_scores[node1] + node_scores[node2]) + (node1_scoretmp + node2_scoretmp);
    node_scores[node1] = node1_scoretmp;
    node_scores[node2] = node2_scoretmp;
}

/**
 * Moves an unsscored node to the right of the sub order and scores. 
 * Eg if (1) is the node to score:
 * [0,(1),2,3,|5,6,4] -> [0,2,3,|(1),5,6,4]
 */
double insert_and_score_at(std::vector<int> &order,
                           int from_index,
                           int to_index,
                           int n,
                           std::vector<double> &node_scores,
                           double &order_score,
                           OrderScoring &scoring)
{
    int p = order.size();
    int node = order[from_index];

    if (to_index < p - n)
    {
        // If the new node is to be put somewhere in the order.
        // [0,(1),2,3|5,6,4] -> [(1),0,2,3|5,6,4]
        move_element(order, from_index, to_index);              // put it in the back [0,(1),2,3,|5,6,4] -> [(1),0,2,3,|,5,6,4]
        node_scores[node] = scoring.score_pos(order, to_index); // O(p)? This should also be quite time consuming..
        order_score += node_scores[node];
    }

    if (to_index >= p - n)
    {
        // If the new node is to be put on the very left of the sub order
        // [0,(1),2,3,|5,6,4] -> [0,2,3,|(1),5,6,4]
        move_element(order, from_index, p - n - 1);              // put it in the back [0,(1),2,3,|5,6,4] -> [0,2,3,|(1),5,6,4]
        node_scores[node] = scoring.score_pos(order, p - n - 1); // O(p)? This should also be quite time consuming..
        order_score += node_scores[node];

        // If the new node is put somewhere in the sub order we wave to shift it in from the left.
        // [0,(1),2,3,|5,6,4] -> [0,2,3,|(1),5,6,4] -> [0,2,3,|5,(1),6,4]
        for (int i = p - n - 1; i < to_index; ++i)
        {
            swap_nodes(i, i + 1, order, node_scores, order_score, scoring);
        }
    }

    return (order_score);
}

bool prune_right_type(std::vector<int> &order,
                      int n,
                      std::vector<double> &node_scores,
                      double &order_score,
                      OrderScoring &scoring)
{
    std::size_t p = order.size();
    for (std::size_t node_index = 0; node_index < p - n; node_index++)
    {
        // Go though all of the un used nodes and inject them at second place
        // of the sub order and at the very front of the order.

        double order_score_bkp = order_score;
        double score_at_top = insert_and_score_at(order, node_index, 0, n, node_scores, order_score, scoring);
        move_element(order, 0, node_index);   // move back the node.
        node_scores[order[node_index]] = 0.0; // reset the node score at pos 0
        order_score = order_score_bkp;        // reset order score

        double node_score_bkp = node_scores[order[p - n]];
        double score_at_pos_2 = insert_and_score_at(order, node_index, p - n, n, node_scores, order_score, scoring);
        move_element(order, p - n, node_index);     // move back the node.
        node_scores[order[node_index]] = 0.0;       // reset the node score at pos 0
        node_scores[order[p - n]] = node_score_bkp; // reset the node score at pos p-n

        order_score = order_score_bkp; // reset order score
        if (score_at_top < score_at_pos_2 + 0.0000001)
        {
            return (true); // Prune
        }
    }
    return (false);
}

void put_nodes_in_front(int n,
                        std::vector<std::vector<std::tuple<std::vector<int>, double, int, std::vector<double>>>> &opt_tuples,
                        OrderScoring &scoring)
{
    // Go through all orders with n-1 nodes. For each order, try to add the remaining nodes
    // at the front.
    std::size_t p = scoring.numparents.size();
    if (n == 1)
    {
        for (size_t node = 0; node < p; node++)
        {
            std::vector<int> order(p, 0);
            // It all vectors with all the nodes.
            // The makes it easier to keep track of which nodes are used.
            for (size_t i = 0; i < p; i++)
            {
                order[i] = i;
            }
            myswap(node, p - 1, order);
            // score nodes
            //PrintVector(order);
            std::vector<double> *sc = scoring.score(order, p - 1, 1); // OK, score only index p-1.
            std::vector<double> node_scores = *sc;
            //PrintVector(node_scores);
            delete sc;
            // score order
            double order_score = node_scores[node];
            int order_number = node;
            std::tuple<std::vector<int>, double, int, std::vector<double>>
                opt_tuple = std::make_tuple(order, order_score, order_number, node_scores);
            //std::cout << "score " << order_score << std::endl;

            //if (!prune_right_type(order, n, node_scores, order_score, scoring))
            {
                opt_tuples[n].push_back(opt_tuple);
            }
        }
        //std::cout << "order score " << std::get<double>(opt_tuples[1][0]) << std::endl;
    }
    else
    {
        //std::cout << "n=" << n << std::endl;
        //std::cout << "orders " << opt_tuples[n - 1].size() << std::endl;
        for (size_t i = 0; i < opt_tuples[n - 1].size(); i++)
        {
            //std::cout << "particle (order) " << i << " out of " << opt_tuples[n - 1].size() << std::endl;
            // take one of the orders.
            const auto &[prev_order, prev_score, prev_order_number, prev_node_scores] = opt_tuples[n - 1][i];
            //PrintVector(prev_order);
            for (size_t node_index = 0; node_index <= p - n; node_index++)
            {
                // Try to put m to the front

                auto [optimal_order, order_score, node_scores, order_number, is_optimal] = put_node_in_front(prev_order,
                                                                                                             n,
                                                                                                             node_index,
                                                                                                             prev_node_scores,
                                                                                                             prev_score,
                                                                                                             prev_order_number,
                                                                                                             scoring);
                // Pruning step I) Keep only if the new node (before at node_index) is optimal at the front.
                // Pruning step II) If it is better to move one of the nodes in the ordering to the front.
                // Pruning step III) If the top and next top nodes can be interchanged while the score is invariant,
                if (is_optimal == true)
                {
                    //if (!prune_right_type(optimal_order, n, node_scores, order_score, scoring))
                    {
                        // Check if it fits better as parent for all others (index 0, or very "left").
                        //PrintVector(optimal_order);
                        //std::cout << order_score << std::endl;
                        std::tuple<std::vector<int>, double, int, std::vector<double>> tuple = std::make_tuple(optimal_order, order_score, order_number, node_scores);
                        opt_tuples[n].push_back(tuple);
                    }
                }
            }
        }
    }
}

std::vector<bool> order_to_boolvec(const std::vector<int> &order, int n)
{
    std::vector<bool> boolvec(order.size(), false);
    std::size_t p = order.size();
    for (std::size_t i = p - n; i < p; i++)
    {
        boolvec[order[i]] = true;
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
std::vector<int> unique_sets(const std::vector<std::vector<bool>> &mats, const std::vector<double> &order_scores)
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
                zeros.push_back(*i);
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

                //PrintVector(zeros);
                // Get the max scoring index
                double maxscore = -INFINITY;
                int max_ind = -1;
                std::vector<int>::iterator i;
                for (i = zeros.begin() + 1; i != zeros.end(); ++i)
                {
                    //std::cout << "row " << *i << " score " << order_scores[*i] << std::endl;
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
                //PrintVector(ones);
                // Get the max scoring index
                double maxscore = -INFINITY;
                int max_ind = -1;
                std::vector<int>::iterator i;
                for (i = ones.begin() + 1; i != ones.end(); ++i)
                {
                    //std::cout << "row " << *i << " score " << order_scores[*i] << std::endl;
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

std::vector<std::tuple<std::vector<int>, double, int, std::vector<double>>> prune_equal_sets(std::vector<std::tuple<std::vector<int>, double, int, std::vector<double>>> &tuples, int n)
{
    std::vector<std::vector<bool>> boolmat;
    std::vector<double> order_scores;

    //std::cout << "Creating boolmatrix" << std::endl;
    for (const auto &t : tuples)
    {
        const auto &[optimal_order, order_score, order_number, node_scores] = t;
        std::vector<bool> boolvec = order_to_boolvec(optimal_order, n);
        boolmat.push_back(std::move(boolvec));
        order_scores.push_back(order_score);
    }

    //std::cout << "Get indices of unique maximal scoring sets" << std::endl;
    std::vector<int> pruned_inds = unique_sets(boolmat, order_scores);

    std::vector<double> pruned_scores;
    std::vector<std::tuple<std::vector<int>, double, int, std::vector<double>>> pruned_tuples;
    for (const auto &ind : pruned_inds)
    {
        const auto &[optimal_order, order_score, order_number, node_scores] = tuples[ind];
        pruned_tuples.push_back(tuples[ind]);
    }

    return (pruned_tuples);
}

std::vector<std::tuple<std::vector<int>, double, int, std::vector<double>>> prune_gaps(std::vector<std::tuple<std::vector<int>, double, int, std::vector<double>>> &tuples,
                                                                                       int n,
                                                                                       OrderScoring &scoring)
{
    std::vector<std::tuple<std::vector<int>, double, int, std::vector<double>>> pruned_tuples;
    int i = 0;
    for (auto &t : tuples)
    {
        auto &[optimal_order, order_score, order_number, node_scores] = t;

        if (!prune_right_type(optimal_order, n, node_scores, order_score, scoring))
        {
            pruned_tuples.push_back(t);
        }
    }
    return (pruned_tuples);
}

std::vector<std::tuple<std::vector<int>, double, int, std::vector<double>>> prune_subsample(std::vector<std::tuple<std::vector<int>, double, int, std::vector<double>>> &tuples)
{
    int seed = 1;
    std::srand(seed);
    std::random_device rd;
    std::mt19937 gen(rd());
    std::default_random_engine generator(seed);
    std::vector<double> pruned_scores;

    for (const auto &t : tuples)
    {
        const auto &[optimal_order, order_score, order_number, node_scores] = t;
        pruned_scores.push_back(order_score);
    }

    std::vector<double> *norm_w = dist_from_logprobs(pruned_scores);
    PrintVector(pruned_scores);
    std::cout << "Normalized vector of scores: " << std::endl;

    std::discrete_distribution<int> distribution(norm_w->begin(), norm_w->end());
    PrintVector(*norm_w);
    std::size_t NN = 100;
    std::size_t N = std::min(NN, tuples.size());
    std::vector<int> I(N);

    std::vector<std::tuple<std::vector<int>, double, int, std::vector<double>>> pruned_tuples;
    delete norm_w;
    // Resample particles
    for (std::size_t i = 0; i < N; i++)
    {
        I[i] = distribution(generator);
        std::cout << I[i] << ", ";
    }
    for (const auto &ind : I)
    {
        pruned_tuples.push_back(tuples[ind]);
    }

    std::cout << std::endl;
    std::cout << pruned_tuples.size() << std::endl;

    return (pruned_tuples);
}

void sequential_opt(OrderScoring &scoring);
void sequential_opt(OrderScoring &scoring)
{

    std::size_t p = scoring.numparents.size();
    std::cout << "Starting optimization" << std::endl;
    std::vector<std::vector<std::tuple<std::vector<int>, double, int, std::vector<double>>>> opt_tuples(p + 1);

    for (std::size_t n = 1; n <= p; n++)
    {
        std::cout << "\nn=" << n << std::endl;
        // This could be done in parallel.
        std::cout << "Add new nodes in the front " << std::endl;
        put_nodes_in_front(n, opt_tuples, scoring);
        // Prune opt_tuples. For all elements with the same set of nodes, keep only one with
        // the highest score.
        std::cout << "# orders: " << opt_tuples[n].size() << std::endl;
        double max_score = -INFINITY;
        std::vector<int> max_order;
        std::vector<double> max_nodescores;
        std::vector<std::vector<bool>> boolmat;
        std::vector<double> order_scores;

        std::cout << "Prune gaps " << std::endl;
        opt_tuples[n] = prune_gaps(opt_tuples[n], n, scoring);
        std::cout << "# orders: " << opt_tuples[n].size() << std::endl;

        std::cout << "Prune equal sets " << std::endl;
        opt_tuples[n] = prune_equal_sets(opt_tuples[n], n);
        std::cout << "# orders: " << opt_tuples[n].size() << std::endl;

        for (const auto &t : opt_tuples[n])
        {
            const auto &[optimal_order, order_score, order_number, node_scores] = t;
            if (order_score > max_score)
            {
                max_score = order_score;
                max_order = optimal_order;
                max_nodescores = node_scores;
            }
        }

        // Print some statistics //

        /* Check that scores are correct */
        std::vector<int> tmpv(max_order.end() - n, max_order.end());
        std::cout << "max scoring sub order " << std::endl;
        PrintVector(tmpv);
        // check correct score
        std::vector<double> *sc = scoring.score(max_order, p - n, n); // Take only the last n elements in the vector
        PrintVector(*sc);
        double max_score_check = std::accumulate(sc->begin(), sc->end(), 0.0);
        delete sc;
        assert(std::abs(max_score - max_score_check) < 0.00001);
        std::cout << "score: " << max_score << std::endl;
    }
}

std::tuple<std::vector<double>, std::vector<std::vector<int>>, std::vector<double>> smc(OrderScoring &scoring, std::size_t N, std::size_t p, std::default_random_engine &generator);
std::tuple<std::vector<double>, std::vector<std::vector<int>>, std::vector<double>> smc(OrderScoring &scoring, std::size_t N, std::size_t p, std::default_random_engine &generator)
{
    std::cout << "Starting SMC " << std::endl;

    std::vector<std::vector<double>> log_w(p, std::vector<double>(N, 0.0));
    std::vector<std::vector<std::vector<int>>> orders(p, std::vector<std::vector<int>>(N, std::vector<int>(p)));
    std::vector<int> I(N);
    int node_index;
    std::vector<int> order(p);
    std::vector<double> neig_scoring;
    std::vector<double> neig_dist;
    int insert_pos;
    double log_prop_prob;
    std::vector<std::vector<std::vector<double>>> log_node_scores(p, std::vector<std::vector<double>>(N, std::vector<double>(p, 0.0)));
    std::vector<std::vector<double>> log_order_scores(p, std::vector<double>(N, 0.0));
    std::vector<double> *norm_w; // = new(N, 0.0);
    std::discrete_distribution<int> distribution;

    for (size_t i = 0; i < N; i++)
    {
        for (size_t n = 0; n < p; n++)
        {
            for (size_t m = 0; m < p; m++)
            {
                orders[n][i][m] = m;
            }
        }
    }

    for (size_t n = 0; n < p; n++)
    {
        //std::cout << "\n\nn: " << n << " cache hits: " << scoring.cache_hits << std::endl;
        for (size_t i = 0; i < N; i++)
        {
            //std::cout << "\ni: " << i << std::endl;
            if (n == 0)
            {
                double log_initprop_prob = random_parent_init(scoring,
                                                              orders[n][i],
                                                              log_node_scores[n][i],
                                                              log_order_scores[n][i],
                                                              generator); // Propose from neighborhood.
                log_w[n][i] = log_order_scores[n][i] - log_initprop_prob;
            }
            else
            {
                double log_prop_prob = neig_proportional_proposal(scoring,
                                                                  orders[n - 1][I[i]],
                                                                  log_node_scores[n - 1][I[i]],
                                                                  log_order_scores[n - 1][I[i]],
                                                                  orders[n][i],
                                                                  log_node_scores[n][i],
                                                                  log_order_scores[n][i],
                                                                  n,
                                                                  generator); // Propose from neighborhood.
                log_w[n][i] = log_order_scores[n][i] - std::log(n + 1) - log_order_scores[n - 1][I[i]] - log_prop_prob;
            }
        }

        // rescale weights to
        norm_w = dist_from_logprobs(log_w[n]);
        std::discrete_distribution<int> distribution(norm_w->begin(), norm_w->end());
        delete norm_w;
        // Resample particles
        for (std::size_t i = 0; i < N; i++)
        {
            I[i] = distribution(generator);
        }
    }
    return (std::make_tuple(log_w[p - 1], orders[p - 1], log_order_scores[p - 1]));
}
std::tuple<std::vector<std::vector<int>>, std::vector<int>> backward_order_sampler(const std::vector<int> &order);
std::tuple<std::vector<std::vector<int>>, std::vector<int>> backward_order_sampler(const std::vector<int> &order)
{
    int p = order.size();
    std::vector<std::vector<int>> orders(p);
    std::vector<int> node_indices(p - 1);

    int node_index;
    orders[p - 1] = order;
    for (int n = p - 2; n >= 0; n--)
    {
        orders[n] = orders[n + 1];
        node_index = rand_int_in_range(0, n + 1);   // Draw random node
        move_element(orders[n], node_index, p - 1); // Move node at node_index to the back.
        node_indices[n] = node_index;
    }

    return std::make_tuple(orders, node_indices);
}

/**
 * 
 * Calculates the log proposal probabilities of the orders in asequence efficiently.
 */
std::vector<double> orders_log_prop_probs(const std::vector<std::vector<int>> &orders, std::vector<int> &new_node_inds_cond_orders, OrderScoring &scoring);
// std::vector<double> orders_log_prop_probs(const std::vector<std::vector<int>> &orders, std::vector<int> &new_node_inds_cond_orders, OrderScoring &scoring)
// {
//     std::size_t p = orders.size();
//     std::vector<double> order_log_prop_probs(orders.size());

//     for (std::size_t i = 0; i < p; i++)
//     {
//         int new_node_index = 1;
//         order_log_scores[i] = order_log_prop_prob(orders[i], i + 1, new_node_index, order_log_scores[i - 1], scoring);
//     }

//     return (order_log_scores);
// }

/**
 * 
 * Calculates the log proposal probability of an order.
 * The strategy is to place the new node in every position and calculate the order scores.
 * Since th node scores depends on the psostion also these as to be calculated.
 */
double order_log_prop_prob(const std::vector<int> &input_order, int prev_order_length, int new_node_index, double current_orderscore, OrderScoring &scoring);
double order_log_prop_prob(const std::vector<int> &input_order, int prev_order_length, int new_node_index, double current_orderscore, OrderScoring &scoring)
{

    std::vector<int> order(input_order);
    // compute these in a sequence n= 1,2,3,...
    int node = order[new_node_index];
    move_element(order, new_node_index, prev_order_length); // Move the element to the back at first. maybe one step left

    //for node at index first_n_elements. The other nodes have the
    // same score since those are still before in the ordering.
    std::vector<double> order_scores(prev_order_length + 1);
    std::vector<double> node_scores(input_order.size());

    node_scores[node] = scoring.score_pos(order, prev_order_length); // this should be the full scoring instead

    order_scores[prev_order_length] = current_orderscore + node_scores[node];

    for (int i = prev_order_length; i > 0; --i) // put it in the back instead
    {
        int node1 = order[i - 1];
        int node2 = order[i];
        const auto &[node1_scoretmp, node2_scoretmp] = scoring.swap_nodes(i - 1, i, order, node_scores); // This should update the orderscore directly maybe.
        order_scores[i - 1] = order_scores[i] - (node_scores[node1] + node_scores[node2]) + (node1_scoretmp + node2_scoretmp);
        node_scores[node1] = node1_scoretmp;
        node_scores[node2] = node2_scoretmp;
    }

    // Sample the new ordering from the order_scores distribution, where the index
    // specifies the position were node is inserted.
    //std::cout << "Scores to be dist" << std::endl;
    //PrintVector(order_scores);
    //std::cout << "Score for index "<< new_node_index <<": " << order_scores[new_node_index] << std::endl;
    std::vector<double> *neig_dist = dist_from_logprobs(order_scores); // TODO: dont need dist, just logs
    //PrintVector(*neig_dist);
    //std::cout << "prop prob" << (*neig_dist)[new_node_index] << std::endl;
    return std::log((*neig_dist)[new_node_index]); // NOTE: This can be -Inf sometimes.
}

void smc_cond_kernel(std::size_t n,
                     std::size_t i_from,
                     std::size_t i_to,
                     const std::vector<int> &I,
                     std::vector<std::vector<double>> &log_w,
                     std::vector<std::vector<std::vector<int>>> &orders,
                     std::vector<std::vector<std::vector<double>>> &log_node_scores,
                     std::vector<std::vector<double>> &log_order_scores,
                     const std::vector<int> &new_node_inds_cond_orders,
                     std::default_random_engine &generator,
                     OrderScoring &scoring)

{

    for (size_t i = i_from; i < i_to; i++)
    {
        int node_index;
        std::size_t p = scoring.numparents.size();
        if (i == 0) // Calculate probabilities for the fixed PGibbs particles.
        {
            std::vector<double> *sc = scoring.score(orders[n][0], 0, n + 1); // Score all the nodes in the sub order [0,..,n].
            log_node_scores[n][0] = *sc;
            delete sc;
            log_order_scores[n][0] = std::accumulate(log_node_scores[n][0].begin(), log_node_scores[n][0].end(), 0.0);
            if (n == 0)
            {
                log_w[n][0] = log_order_scores[n][0] - (-std::log(p));
            }
            else
            {
                double log_prop_prob = order_log_prop_prob(orders[n][0], n, new_node_inds_cond_orders[n - 1], log_order_scores[n - 1][0], scoring);
                if (!std::isfinite(log_prop_prob))
                {
                    log_prop_prob = std::numeric_limits<double>::min();
                }
                log_w[n][0] = log_order_scores[n][0] - std::log(n + 1) - log_order_scores[n - 1][0] - (log_prop_prob - std::log(p - n));
                // TODO: BUG?: Update log node scores!
            }
        }

        if (i > 0) // Update the rest of the particles.
        {
            if (n == 0)
            {
                double log_initprop_prob = random_parent_init(scoring,
                                                              orders[n][i],
                                                              log_node_scores[n][i],
                                                              log_order_scores[n][i],
                                                              generator); // Propose from neighborhood.
                log_w[n][i] = log_order_scores[n][i] - log_initprop_prob;
            }
            else
            {
                double log_prop_prob = neig_proportional_proposal(scoring,
                                                                  orders[n - 1][I[i]],
                                                                  log_node_scores[n - 1][I[i]],
                                                                  log_order_scores[n - 1][I[i]],
                                                                  orders[n][i],
                                                                  log_node_scores[n][i],
                                                                  log_order_scores[n][i],
                                                                  n,
                                                                  generator); // Propose from neighborhood.
                log_w[n][i] = log_order_scores[n][i] - std::log(n + 1) - log_order_scores[n - 1][I[i]] - log_prop_prob;
            }
        }
    }
}

void smc_cond_kernel2(std::size_t n,
                      std::size_t i,
                      const std::vector<int> &I,
                      std::vector<std::vector<double>> &log_w,
                      std::vector<std::vector<std::vector<int>>> &orders,
                      std::vector<std::vector<std::vector<double>>> &log_node_scores,
                      std::vector<std::vector<double>> &log_order_scores,
                      const std::vector<int> &new_node_inds_cond_orders,
                      std::default_random_engine &generator,
                      OrderScoring &scoring)
{

    smc_cond_kernel(n,
                    i,
                    i + 1,
                    I,
                    log_w,
                    orders,
                    log_node_scores,
                    log_order_scores,
                    new_node_inds_cond_orders,
                    generator,
                    scoring);
}

std::tuple<std::vector<double>, std::vector<std::vector<int>>, std::vector<double>> smc_cond(OrderScoring &scoring, std::size_t N,
                                                                                             const std::vector<std::vector<int>> &cond_orders,
                                                                                             const std::vector<int> &new_node_inds_cond_orders,
                                                                                             std::default_random_engine &generator,
                                                                                             thread_pool &pool);
std::tuple<std::vector<double>, std::vector<std::vector<int>>, std::vector<double>> smc_cond(OrderScoring &scoring, std::size_t N,
                                                                                             const std::vector<std::vector<int>> &cond_orders,
                                                                                             const std::vector<int> &new_node_inds_cond_orders,
                                                                                             std::default_random_engine &generator,
                                                                                             thread_pool &pool)
{
    std::size_t p = scoring.numparents.size();
    //std::cout << "Starting conditional SMC " << std::endl;
    std::vector<std::vector<double>> log_w(p, std::vector<double>(N, 0.0));
    std::vector<std::vector<std::vector<int>>> orders(p, std::vector<std::vector<int>>(N, std::vector<int>(p)));
    std::vector<int> I(N);
    int node_index;
    double log_prop_prob;
    std::vector<std::vector<std::vector<double>>> log_node_scores(p, std::vector<std::vector<double>>(N, std::vector<double>(p, 0.0)));
    std::vector<std::vector<double>> log_order_scores(p, std::vector<double>(N));
    std::vector<double> *norm_w;

    for (size_t i = 0; i < N; i++)
    {
        for (size_t n = 0; n < p; n++)
        {
            for (size_t m = 0; m < p; m++)
            {
                if (i == 0)
                {
                    orders[n][i][m] = cond_orders[n][m];
                }
                else
                {
                    orders[n][i][m] = m;
                }
            }
        }
    }

    for (size_t n = 0; n < p; n++)
    {

        // std::threads for a range seems to be faster than using thread_pool for each paticle or range.
        unsigned int maxthreads = N;
        std::size_t Num_Threads = std::min(std::thread::hardware_concurrency(), maxthreads);
        std::vector<std::thread> threads(Num_Threads);
        //std::cout << "num threads" << Num_Threads << std::endl;

        for (size_t i = 0; i < Num_Threads; i++)
        {
            int i_from = i * (N / Num_Threads);
            int i_to = (i + 1) * (N / Num_Threads);
            //std::cout << i_from << "-" << i_to << std::endl;

            if (i == Num_Threads - 1)
            {
                i_to = N;
            }

            threads[i] = std::thread(smc_cond_kernel,
                                     n,
                                     i_from,
                                     i_to,
                                     std::ref(I),
                                     std::ref(log_w),
                                     std::ref(orders),
                                     std::ref(log_node_scores),
                                     std::ref(log_order_scores),
                                     std::ref(new_node_inds_cond_orders),
                                     std::ref(generator),
                                     std::ref(scoring));
        }

        for (std::size_t i = 0; i < Num_Threads; i++)
        {
            threads[i].join();
        }

        //PrintVector(log_w[n]);
        // rescale weights to
        norm_w = dist_from_logprobs(log_w[n]);
        std::discrete_distribution<int> distribution(norm_w->begin(), norm_w->end());
        delete norm_w;
        // Resample particles
        for (std::size_t i = 0; i < N; i++)
        {
            I[i] = distribution(generator);
        }
    }

    return (std::make_tuple(log_w[p - 1], orders[p - 1], log_order_scores[p - 1]));
}

std::pair<std::vector<std::vector<int>>, std::vector<double>> pgibbs(std::size_t M, std::size_t N, OrderScoring &scoring, std::default_random_engine &generator, thread_pool &pool)
{
    std::cout << "Starting PGibbs" << std::endl;
    std::vector<std::vector<int>> orders(M);
    std::vector<double> log_scores(M);
    std::vector<double> *norm_w;
    std::size_t p = scoring.numparents.size();
    std::cout << "p " << p << std::endl;

    const auto &[smc_log_w, smc_orders, smc_log_scores] = smc(scoring, N, p, generator);

    // Sample order from log_w and put in orders.
    norm_w = dist_from_logprobs(smc_log_w);
    std::discrete_distribution<int> distribution(norm_w->begin(), norm_w->end());
    int sampled_index = distribution(generator);
    delete norm_w;
    orders[0] = smc_orders[sampled_index];
    log_scores[0] = smc_log_scores[sampled_index];
    //orders[0] = [6, 7, 0, 5, 1, 3, 4, 2]

    std::cout << "Order sampled from SMC " << std::endl;
    PrintVector(orders[0]);
    for (std::size_t j = 1; j < M; j++)
    {

        std::cout << "PGibbs sample  " << j << std::endl;
        // Backward sample orders.
        const auto &[orders_cond_traj, new_node_inds_cond_orders] = backward_order_sampler(orders[j - 1]);
        // std::cout << "Orders to condition on (generated by backward sampler)" << std::endl;
        // for(auto& o : orders_cond_traj){
        //     PrintVector(o);
        // }
        const auto &[smc_cond_log_w, smc_cond_orders, smc_cond_log_scores] = smc_cond(scoring, N, orders_cond_traj, new_node_inds_cond_orders, generator, pool);

        norm_w = dist_from_logprobs(smc_cond_log_w);

        std::discrete_distribution<int> distribution(norm_w->begin(), norm_w->end());
        sampled_index = distribution(generator);

        delete norm_w;
        // std::vector<int> pgibbsorder = smc_cond_orders[sampled_index];
        // float pgibbslogscore = smc_cond_log_scores[sampled_index];

        // // Add MH swaps about here with pgibbsorder
        // // Need the log node scores here
        // std::vector<double> *log_nodes_scores = scoring.score(pgibbsorder, 0, p);
        // //log_scores[0] = std::accumulate(log_nodes_scores->begin(), log_nodes_scores->end(), 0.0);
        // const auto &[mh_order, mh_log_score] = mh_swap_move(pgibbsorder, *log_nodes_scores, pgibbslogscore, scoring);

        // orders[j] = mh_order;
        // log_scores[j] = mh_log_score;

        orders[j] = smc_cond_orders[sampled_index];
        log_scores[j] = smc_cond_log_scores[sampled_index];
        std::cout << "PGibbs sample score " << log_scores[j] << std::endl;
    }

    return (std::make_pair(orders, log_scores));
}

std::pair<std::vector<int>, std::vector<double>> mh(std::size_t M, OrderScoring &scoring, std::default_random_engine &generator)
{
    std::size_t p = scoring.numparents.size();
    //std::vector<std::vector<int>> orders(M, std::vector<int>(p));
    //std::vector<std::vector<int>> orders(2, std::vector<int>(2));
    std::vector<double> log_scores(M);

    // Init vector

    std::vector<int> order = std::vector<int>(p);
    //orders[0] = std::vector<int>(p);
    for (std::size_t i = 0; i < p; i++)
    {
        //orders[0][i] = i;
        order[i] = i;
    }

    //std::vector<double> *log_nodes_scores = scoring.score(orders[0], 0, p);
    std::vector<double> *log_nodes_scores = scoring.score(order, 0, p);
    log_scores[0] = std::accumulate(log_nodes_scores->begin(), log_nodes_scores->end(), 0.0);
    double max_score = -INFINITY;
    std::vector<int> max_order;

    for (std::size_t j = 1; j < M; j++)
    {
        //const auto &[mh_order, mh_log_score] = mh_swap_move(orders[j - 1], *log_nodes_scores, log_scores[j - 1], scoring);
        const auto &[mh_order, mh_log_score] = mh_swap_move(order, *log_nodes_scores, log_scores[j - 1], scoring);
        //orders[j] = mh_order;
        order = mh_order;
        log_scores[j] = mh_log_score;
        if (mh_log_score > max_score)
        {
            max_score = mh_log_score;
            max_order = order;
        }
    }

    delete log_nodes_scores;

    return (std::make_pair(max_order, log_scores));
}

OrderScoring get_score(Rcpp::List ret)
{

    // Read MAP flag
    bool MAP = Rcpp::as<bool>(ret["MAP"]);
    std::cout << MAP << std::endl;

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
    //std::vector<Rcpp::NumericMatrix> bannedscore;
    std::vector<std::vector<std::vector<double>>> bannedscore;
    for (std::size_t i = 0; i < p; i++)
    {
        Rcpp::NumericMatrix m = Rcpp::as<Rcpp::NumericMatrix>(bannedscoreR[i]);
        //std::vector<std::vector<double>> mat;
        std::vector<std::vector<double>> mat(m.rows(), std::vector<double>(m.cols()));
        for (int j = 0; j < m.rows(); j++)
        {
            for (int k = 0; k < m.cols(); k++)
            {
                mat[j][k] = m(j, k);
                //NumericVector v = m( _ , j );
                //std::vector<double> row;
                //for(auto val: v){
                //    row.push_back(val);
                //}
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

    std::map<cache_keytype3, std::vector<double>> cache;
    OrderScoring scoring(aliases,
                         numparents,
                         rowmaps_backwards,
                         plus1listsparents,
                         scoretable,
                         bannedscore,
                         MAP, cache);

    return (scoring);
}

// [[Rcpp::plugins(cpp17)]]

// [[Rcpp::export]]
Rcpp::List r_pgibbs(Rcpp::List ret, int N, int M, int seed)
{
    OrderScoring scoring = get_score(ret);

    std::srand(seed);
    std::random_device rd;
    std::mt19937 gen(rd());
    std::default_random_engine generator(seed);

    thread_pool pool;
    //pool.push_task(task, arg1, arg2);
    //pool.wait_for_tasks();

    const auto &[pgibbs_orders, pgibbs_log_scores] = pgibbs(M, N, scoring, generator, pool);

    Rcpp::List L = Rcpp::List::create(Rcpp::Named("pgibbs_orders") = pgibbs_orders,
                                      Rcpp::Named("pgibbs_log_scores") = pgibbs_log_scores);

    return (L);
}

// [[Rcpp::plugins(cpp17)]]

// [[Rcpp::export]]
Rcpp::List r_mh(Rcpp::List ret, int M, int seed)
{
    OrderScoring scoring = get_score(ret);

    std::srand(seed);
    std::random_device rd;
    std::mt19937 gen(rd());
    std::default_random_engine generator(seed);

    const auto &[mh_orders, mh_log_scores] = mh(M, scoring, generator);

    Rcpp::List L = Rcpp::List::create(Rcpp::Named("mh_orders") = mh_orders,
                                      Rcpp::Named("mh_log_scores") = mh_log_scores);

    return (L);
}

// std::vector<double> orderscorePlus1(std::size_t n,
//                                     std::vector<int> &scorenodes,
//                                     std::vector<int> &scorepositions,
//                                     std::vector<std::vector<int>> &aliases,
//                                     std::vector<int> &numparents,
//                                     std::vector<Rcpp::IntegerVector> &rowmaps_backwards,
//                                     std::vector<Rcpp::IntegerVector> &plus1listsparents,
//                                     std::vector<std::vector<std::vector<double>>> &scoretable,
//                                     std::vector<Rcpp::NumericMatrix> &scoresmatrices,
//                                     std::vector<int> &permy);

// std::vector<double> orderscorePlus1(std::size_t n,
//                                     std::vector<int> &scorenodes,
//                                     std::vector<int> &scorepositions,
//                                     std::vector<std::vector<int>> &aliases,
//                                     std::vector<int> &numparents,
//                                     std::vector<Rcpp::IntegerVector> &rowmaps_backwards,
//                                     std::vector<Rcpp::IntegerVector> &plus1listsparents,
//                                     std::vector<std::vector<std::vector<double>>> &scoretable,
//                                     std::vector<Rcpp::NumericMatrix> &scoresmatrices,
//                                     std::vector<int> &permy)
// {

//     std::vector<double> orderscores(n);                 // orderscores <- vector("double", n)
//     std::vector<std::vector<int>> allowedscorelists(n); // allowedscorelists < -vector("list", n)
//     std::vector<int> therows(n);                        // therows <- vector("integer", n) # what is this? / Felix
//     size_t k = 0;                                       // k <- 1

//     for (int i : scorenodes)
//     {
//         std::cout << "node " << i << std::endl;
//         size_t position = scorepositions[k];
//         std::cout << "position " << position << std::endl;
//         if (position == n - 1)
//         {
//             // no parents allowed, i.e.only first row, only first list
//             orderscores[i] = scoretable[i][0][0];    // orderscores[i] <- scoretable [[i]][[1]][1, 1]
//             allowedscorelists[i].assign({0});        // allowedscorelists [[i]] <- c(1)
//             therows[i] = std::pow(2, numparents[i]); // therows[i] <- c(2 ^ numparents[i])
//         }
//         else
//         {
//             std::vector<int> bannednodes(permy.begin(), permy.begin() + position + 1); // bannednodes <- permy[1:position]
//             std::vector<int> allowednodes(permy.begin() + position + 1, permy.end());  // allowednodes < -permy [(position + 1):n]
//             std::vector<int> bannedpool;                                               // bannedpool <- which(aliases [[i]] % in % bannednodes)
//             for (int j = 0; j < aliases[i].size(); j++)
//             {
//                 if (std::find(bannednodes.begin(), bannednodes.end(), aliases[i][j]) != bannednodes.end())
//                 {
//                     bannedpool.push_back(j);
//                 }
//             }

//             if (numparents[i] == 0 || bannedpool.size() == 0)
//             {
//                 therows[i] = 0; // all parents allowed
//             }
//             else
//             {
//                 int indextmp = 0;
//                 for (auto &item : bannedpool)
//                 {
//                     indextmp += std::pow(2, item + 1); // add 1 since nodes are labeled from 0
//                 }
//                 indextmp = indextmp / 2;
//                 therows[i] = rowmaps_backwards[i][indextmp]; // this might be list I guess. rowmaps_backwards[i][std::sum(2 ^ bannedpool) / 2 + 1]
//             }

//             allowedscorelists[i].push_back(0); // allowedscorelists [[i]] <- c(1, which(plus1lists$parents [[i]] % in % allowednodes) + 1)
//             // Optimize when optimizing order.
//             for (int j = 0; j < plus1listsparents[i].size(); j++)
//             {
//                 if (std::find(allowednodes.begin(), allowednodes.end(), plus1listsparents[i][j]) != allowednodes.end())
//                 {
//                     allowedscorelists[i].push_back(j + 1);
//                 }
//             }

//             std::vector<double> scoresvec;
//             for (int allowedscore : allowedscorelists[i])
//             {
//                 scoresvec.push_back(scoresmatrices[i](therows[i], allowedscore));
//             }

//             double maxallowed = *std::max_element(scoresvec.begin(), scoresvec.end()); // maxallowed <- max(scoresvec)

//             for (auto &val : scoresvec)
//             {
//                 orderscores[i] += std::exp(val - maxallowed);
//             }
//             orderscores[i] = maxallowed + std::log(orderscores[i]);
//         }
//         k++;
//     }
//     //scores < -list()
//     //scores$therow < -therows
//     //scores$allowedlists < -allowedscorelists
//     //scores$totscores < -orderscores
//     return (orderscores);
// }

// std::valarray<int> numbering(int n, std::vector<int> inds, const std::valarray<bool> &mats, std::valarray<int> &vec)
//     std::valarray<int> numbering(const std::valarray<bool> &mats, int n, std::valarray<int> &vec)
// {
//     int N = vec.size();
//     // Find indices for zeros
//     std::valarray<bool> one_inds(N);
//     std::valarray<bool> zero_inds(N);

//     for (size_t i = 0; i < N; i++)
//     {
//         one_inds[i] = mats[std::slice(n * N, N, 1)] && true;
//         zero_inds[i] = mats[std::slice(n * N, N, 1)] && false;
//     }

//     if (n == vec.size() - 1)
//     {
//         vec[one_inds] += zero_inds.sum();
//         //vec[zero_inds] += 0;
//     }
//     else
//     {
//         // The indices of the 1's.
//         std::valarray<int> one_numbers = numbering(mats[std::slice(n * N, N, 1)][one_inds], n + 1, vec[one_inds]);
//         vec[one_inds] += zero_inds.size() + one_numbers;
//         // Add the number of zeros to the one inds.
//         std::valarray<int> zero_numbers = numbering(mats[std::slice(n * N, N, 1)][zero_inds], n + 1, vec[zero_inds]);
//         vec[zero_inds] += zero_numbers;
//     }
//     return (vec);
// }

// std::tuple<std::vector<int>, double, double> proposal(const OrderScoring &scoring,
//                                                       std::vector<int> input_order,
//                                                       std::vector<double> node_scores,
//                                                       double order_score,
//                                                       int first_n_elements, int index_of_el_to_insert, std::default_random_engine generator)
// {

//     std::vector<double> order_scores = score_sub_order_neigh(scoring, input_order, node_scores, order_score, int_)
//         std::vector<double>
//             neig_dist = dist_from_logprobs(order_scores);                             // Sample one order in the neighborhood
//     std::discrete_distribution<int> distribution(neig_dist.begin(), neig_dist.end()); // Create a distribution. O(p+1)
//     int insert_pos = distribution(generator);                                         // O(N)

//     // move mode to index insert_pos
//     std::vector<int> ret_order(input_order);

//     move_element(ret_order, index_of_el_to_insert, insert_pos); // Move the element to the selectd position
//     // TODO: Need to return log_node_scores
//     return (std::make_tuple(ret_order, std::log(neig_dist[insert_pos]), order_scores[insert_pos]));
// }