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
    copy(arr.begin(), arr.end(), std::ostream_iterator<T>(std::cout, ", "));
    std::cout << std::endl;
}

template <typename T>
void PrintVector(const std::vector<T> &arr)
{
    copy(arr.begin(), arr.end(),
         std::ostream_iterator<T>(std::cout, ", "));
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
        std::map<cache_keytype3, std::vector<double>> cache) : potential_parents(potential_parents),
                                                               numparents(numparents),
                                                               rowmaps_backwards(rowmaps_backwards),
                                                               potential_plus1_parents(potential_plus1_parents),
                                                               scoretable(scoretable),
                                                               scoresmatrices(scoresmatrices),
                                                               cache(cache),
                                                               cache_hits(0)
    {
    }

    /**
     * Re-calculation scores after swapping up node so that (node_a, node_b) --> (node_b, node_a).
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
        double nodea_score;
        double nodeb_score;
        // Computing score for nodea, which is moved up
        // If b is a potential parent for a, we have to recompute the scores since b i now banned.
        if (std::find(potential_parents[node_a].begin(), potential_parents[node_a].end(), node_b) != potential_parents[node_a].end())
        {
            std::swap(ordering[nodea_index], ordering[nodeb_index]);
            nodea_score = score_pos(ordering, nodeb_index);
            std::swap(ordering[nodea_index], ordering[nodeb_index]); // Swap back since mor computations has do be done
        }
        else
        { // Since b is not a potential parent of a, f_bar_z is not altered.
            // Check if b was a plus1 parent for a
            std::vector<int>::iterator itr = std::find(potential_plus1_parents[node_a].begin(), potential_plus1_parents[node_a].end(), node_b);
            if (itr != potential_plus1_parents[node_a].end())
            {
                // Subtract the bs plus1 score contibution.
                std::swap(ordering[nodea_index], ordering[nodeb_index]);
                int f_bar_z = get_f_bar_z(nodeb_index, ordering); // ok?
                // find the correct index j and take it.
                int plus1_parents_index = std::distance(potential_plus1_parents[node_a].begin(), itr);
                //double nodeb_as_plus1_score = scoresmatrices[node_a][f_bar_z][(active_plus1_parents_indices)[plus1_parents_index+1]];
                double nodeb_as_plus1_score = scoresmatrices[node_a][f_bar_z][plus1_parents_index + 1];
                std::swap(ordering[nodea_index], ordering[nodeb_index]); // Swap back

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
                    std::swap(ordering[nodea_index], ordering[nodeb_index]);
                    nodea_score = score_pos(ordering, nodeb_index);
                    std::swap(ordering[nodea_index], ordering[nodeb_index]); // Swap back since more computations has do be done
                }

                //std::cout << "true score " << true_score << " calcuated score " << node_scores[node_a] << std::endl;
                //assert(std::abs(nodea_score-true_score) < 0.001);
            }
        }

        // // Computing score for node_b, which is moved down
        if (std::find(potential_parents[node_b].begin(), potential_parents[node_b].end(), node_a) != potential_parents[node_b].end())
        {
            std::swap(ordering[nodea_index], ordering[nodeb_index]);
            nodeb_score = score_pos(ordering, nodea_index);
            std::swap(ordering[nodea_index], ordering[nodeb_index]);
        }
        else
        {
            std::vector<int>::iterator itr = std::find(potential_plus1_parents[node_b].begin(),
                                                       potential_plus1_parents[node_b].end(),
                                                       node_a);
            if (itr != potential_plus1_parents[node_b].cend())
            {
                int plus1_parents_index = std::distance(potential_plus1_parents[node_b].begin(), itr); // since no parents is the first
                std::swap(ordering[nodea_index], ordering[nodeb_index]);

                int f_bar_z = get_f_bar_z(nodea_index, ordering);
                double nodea_as_plus1_score = scoresmatrices[node_b][f_bar_z][plus1_parents_index + 1];
                std::swap(ordering[nodea_index], ordering[nodeb_index]); // Swap back.
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

        std::swap(ordering[nodea_index], ordering[nodeb_index]); // last swap

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
    std::vector<double> *score(const std::vector<int> &ordering, const int &from_orderpos, const int &n_elements) const
    {
        std::size_t n = ordering.size();
        std::vector<double> *orderscores = new std::vector<double>(n, 0.0); // O(p)           // orderscores <- vector("double", n)
        std::vector<int> active_plus1_parents_indices;                      // active_plus1_parents_indices < -vector("list", n)
        int f_bar_z;

        for (int position = from_orderpos; position < from_orderpos + n_elements; ++position)
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

                for (int j = 0; j < plus1_parents_scores.size(); j++)
                {
                    plus1_parents_scores[j] = scoresmatrices[node][f_bar_z][(active_plus1_parents_indices)[j]]; // allowedscorelist is in numerical order
                }
                (*orderscores)[node] = sum_log_probs(plus1_parents_scores);
                //delete active_plus1_parents_indices;
            }
        }
        return (orderscores);
    }

    double score_pos(const std::vector<int> &ordering, const int &position) const
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

            for (int j = 0; j < plus1_parents_scores.size(); j++)
            {
                plus1_parents_scores[j] = scoresmatrices[node][f_bar_z][active_plus1_parents_indices[j]]; // allowedscorelist is in numerical order
            }
            orderscore = sum_log_probs(plus1_parents_scores);
            //delete active_plus1_parents_indices;
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
        (active_plus1_parents_indices).push_back(0);                   // f(null)=0, no )=parents is always a possibility.?
        for (int j = 0; j < potential_plus1_parents[node].size(); j++) // Is j the plus1node? -No, but potential_plus1_parents[node][j] is.
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
        for (int j = 0; j < potential_parents[node].size(); j++)
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

/**
 *  Asserts index_to < index_from.
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
    int n = log_probs.size();
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
 * Score neigborhood of ordering.
 * when the element at index index_of_el_to_insert is inserted at positions 0,..,first_n_elements.
 */

double score_sub_order_neigh(OrderScoring &scoring,
                             const std::vector<int> &input_order,
                             const std::vector<double> &nodes_scores,
                             const double &input_order_score,
                             std::vector<int> &output_order,
                             std::vector<double> &output_node_scores,
                             double &output_order_score,
                             const int &first_n_elements,
                             const int &index_of_el_to_insert,
                             std::default_random_engine &generator);

double score_sub_order_neigh(OrderScoring &scoring,
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
    double order_score(input_order_score);
    int node = order[index_of_el_to_insert];

    // Calculate score for node at index first_n_elements. The other nodes have the
    // same score sine those are still befor in the ordering.
    move_element(order, index_of_el_to_insert, first_n_elements); // put it in the back
    node_scores[node] = scoring.score_pos(order, first_n_elements);
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
    return log_prob;
}

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

int rand_int_in_range(const std::size_t &from, const std::size_t &to)
{
    return (from + (std::rand() % (to - from + 1)));
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
                // This is the initialisation
                node_index = rand_int_in_range(n, p - 1); // Draw random node
                if (node_index != 0)
                {
                    move_element(orders[n][i], node_index, n); // Move to node_index to index n=0
                }
                std::vector<double> *sc = scoring.score(orders[n][i], n, 1); // OK, score only index 0
                log_node_scores[n][i] = *sc;
                delete sc;
                log_order_scores[n][i] = log_node_scores[n][i][orders[n][i][n]]; // std::accumulate(log_node_scores[i].begin(), log_node_scores[i].end(), 0)
                log_w[n][i] = log_order_scores[n][i] + std::log(p);              // First weight
            }
            else
            {
                node_index = rand_int_in_range(n, p - 1); // Draw one of the remaining nodes
                log_prop_prob = score_sub_order_neigh(scoring,
                                                      orders[n - 1][I[i]],
                                                      log_node_scores[n - 1][I[i]],
                                                      log_order_scores[n - 1][I[i]],
                                                      orders[n][i],
                                                      log_node_scores[n][i],
                                                      log_order_scores[n][i],
                                                      n,
                                                      node_index,
                                                      generator); // Propose from neighborhood.
                log_w[n][i] = log_order_scores[n][i] - std::log(n + 1) - log_order_scores[n - 1][I[i]] - (log_prop_prob - std::log(p - n));
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
    // if (prev_order_length == 0)
    // {
    //     return (-std::log(input_order.size()));
    // }

    std::vector<int> order(input_order);
    // compute these in a sequence n= 1,2,3,...
    int node = order[new_node_index];
    move_element(order, new_node_index, prev_order_length); // Move the element to the back at first. maybe one step left

    //for node at index first_n_elements. The other nodes have the
    // same score sine those are still before in the ordering.
    std::vector<double> order_scores(prev_order_length + 1);
    //std::vector<double> node_scores(first_n_elements);
    std::vector<double> node_scores(input_order.size());

    node_scores[node] = scoring.score_pos(order, prev_order_length); // this should be the full scoring instead
    //order_scores[first_n_elements] = order_score + node_scores[node];
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

// std::tuple<std::vector<double>, std::vector<std::vector<int>>, std::vector<double>> smc_cond(OrderScoring &scoring, std::size_t N,
//                                                                                              const std::vector<std::vector<int>> &cond_orders,
//                                                                                              const std::vector<int> &new_node_inds_cond_orders,
//                                                                                              std::default_random_engine &generator);
// std::tuple<std::vector<double>, std::vector<std::vector<int>>, std::vector<double>> smc_cond(OrderScoring &scoring, std::size_t N,
//                                                                                              const std::vector<std::vector<int>> &cond_orders,
//                                                                                              const std::vector<int> &new_node_inds_cond_orders,
//                                                                                              std::default_random_engine &generator)
// {
//     std::size_t p = scoring.numparents.size();
//     std::cout << "Starting conditional SMC " << std::endl;
//     std::vector<std::vector<double>> log_w(p, std::vector<double>(N, 0.0));
//     std::vector<std::vector<std::vector<int>>> orders(p, std::vector<std::vector<int>>(N, std::vector<int>(p)));
//     std::vector<int> I(N);
//     int node_index;
//     double log_prop_prob;
//     std::vector<std::vector<std::vector<double>>> log_node_scores(p, std::vector<std::vector<double>>(N, std::vector<double>(p, 0.0)));
//     std::vector<std::vector<double>> log_order_scores(p, std::vector<double>(N));
//     std::vector<double> *norm_w;

//     for (size_t i = 0; i < N; i++)
//     {
//         for (size_t n = 0; n < p; n++)
//         {
//             for (size_t m = 0; m < p; m++)
//             {
//                 if (i == 0)
//                 {
//                     orders[n][i][m] = cond_orders[n][m];
//                 }
//                 else
//                 {
//                     orders[n][i][m] = m;
//                 }
//             }
//         }
//     }

//     for (size_t n = 0; n < p; n++)
//     {
//         //std::cout << "\n\nn: " << n << " cache hits: " << scoring.cache_hits << " pgibbs order (first " << n + 1 << " elements) " << std::endl;
//         //PrintVector(orders[n][0]);
//         // Compute the conditional orders for particle Gibbs.
//         std::vector<double> *sc = scoring.score(orders[n][0], 0, n + 1); // Score all the nodes in the sub order [0,..,n].
//         log_node_scores[n][0] = *sc;

//         //std::cout << "Node scores " << std::endl;
//         //PrintVector(log_node_scores[n][0]);
//         delete sc;
//         log_order_scores[n][0] = std::accumulate(log_node_scores[n][0].begin(), log_node_scores[n][0].end(), 0.0); // TODO: this may be too slow.

//         //std::cout << "fixed order score "  << log_order_scores[n][0] << std::endl;
//         //std::cout << "Order score " << log_order_scores[n][0] << std::endl;

//         if (n == 0)
//         {
//             log_w[n][0] = log_order_scores[n][0] - (-std::log(p));
//         }
//         else
//         {
//             // std::cout << "New node inds" << std::endl;
//             // PrintVector(new_node_inds_cond_orders);
//             double log_prop_prob, log_prop_prob2 = order_log_prop_prob(orders[n][0], n, new_node_inds_cond_orders[n - 1], log_order_scores[n - 1][0], scoring);
//             //std::cout << "Log prop prob " << log_prop_prob << std::endl;
//             assert(std::exp(log_prop_prob) >= 0 && std::exp(log_prop_prob) <= 1);

//             if (!std::isfinite(log_prop_prob))
//             {
//                 log_prop_prob = std::numeric_limits<double>::min();
//                 //log_w[n][0] =
//             }
//             log_w[n][0] = log_order_scores[n][0] - std::log(n + 1) - log_order_scores[n - 1][0] - (log_prop_prob - std::log(p - n));
//             //std::cout << "w " << log_w[n][0] << std::endl;
//             //assert(std::isfinite(log_prop_prob2));
//             // TODO: BUG: Update log node scores!
//         }

//         // Update the rest of the particles.
//         for (size_t i = 1; i < N; i++)
//         {
//             //std::cout << "\ni: " << i << std::endl;
//             if (n == 0)
//             {
//                 // This is the initialisation
//                 node_index = rand_int_in_range(n, p - 1); // Draw random node

//                 if (node_index != 0)
//                 {
//                     move_element(orders[n][i], node_index, n); // Move to node_index to index n=0
//                 }
//                 std::vector<double> *sc = scoring.score(orders[n][i], n, 1); // OK, score only index 0
//                 log_node_scores[n][i] = *sc;
//                 delete sc;
//                 log_order_scores[n][i] = log_node_scores[n][i][orders[n][i][n]]; // std::accumulate(log_node_scores[i].begin(), log_node_scores[i].end(), 0.0)
//                 log_w[n][i] = log_order_scores[n][i] - (-std::log(p));           // First weight
//             }
//             else
//             {
//                 node_index = rand_int_in_range(n, p - 1); // Draw one of the remaining nodes
//                 log_prop_prob = score_sub_order_neigh(scoring,
//                                                       orders[n - 1][I[i]],
//                                                       log_node_scores[n - 1][I[i]],
//                                                       log_order_scores[n - 1][I[i]],
//                                                       orders[n][i],
//                                                       log_node_scores[n][i],
//                                                       log_order_scores[n][i],
//                                                       n,
//                                                       node_index,
//                                                       generator); // Propose from neighborhood.

//                 //std::cout << "order score " << log_order_scores[n][i] << std::endl;
//                 log_w[n][i] = log_order_scores[n][i] - std::log(n + 1) - log_order_scores[n - 1][I[i]] - (log_prop_prob - std::log(p - n));
//             }
//         }
//         //PrintVector(log_w[n]);
//         // rescale weights to
//         norm_w = dist_from_logprobs(log_w[n]);
//         std::discrete_distribution<int> distribution(norm_w->begin(), norm_w->end());
//         delete norm_w;
//         // Resample particles
//         for (std::size_t i = 0; i < N; i++)
//         {
//             I[i] = distribution(generator);
//         }
//     }

//     return (std::make_tuple(log_w[p - 1], orders[p - 1], log_order_scores[p - 1]));
// }

void smc_cond_kernel(std::size_t n,
                     std::size_t i,
                     //std::size_t i_from,
                     //std::size_t i_to,
                     const std::vector<int> &I,
                     std::vector<std::vector<double>> &log_w,
                     std::vector<std::vector<std::vector<int>>> &orders,
                     std::vector<std::vector<std::vector<double>>> &log_node_scores,
                     std::vector<std::vector<double>> &log_order_scores,
                     const std::vector<int> &new_node_inds_cond_orders,
                     std::default_random_engine &generator,
                     OrderScoring &scoring)
{

    //for (size_t i = i_from; i < i_to; i++)
    {
        int node_index;
        std::size_t p = scoring.numparents.size();
        if (i == 0)
        {
            std::vector<double> *sc = scoring.score(orders[n][0], 0, n + 1); // Score all the nodes in the sub order [0,..,n].
            log_node_scores[n][0] = *sc;

            //std::cout << "Node scores " << std::endl;
            //PrintVector(log_node_scores[n][0]);
            delete sc;
            log_order_scores[n][0] = std::accumulate(log_node_scores[n][0].begin(), log_node_scores[n][0].end(), 0.0); // TODO: this may be too slow.

            //std::cout << "fixed order score "  << log_order_scores[n][0] << std::endl;
            //std::cout << "Order score " << log_order_scores[n][0] << std::endl;

            if (n == 0)
            {
                log_w[n][0] = log_order_scores[n][0] - (-std::log(p));
            }
            else
            {
                // std::cout << "New node inds" << std::endl;
                // PrintVector(new_node_inds_cond_orders);
                double log_prop_prob, log_prop_prob2 = order_log_prop_prob(orders[n][0], n, new_node_inds_cond_orders[n - 1], log_order_scores[n - 1][0], scoring);
                //std::cout << "Log prop prob " << log_prop_prob << std::endl;
                //assert(std::exp(log_prop_prob) >= 0 && std::exp(log_prop_prob) <= 1);

                if (!std::isfinite(log_prop_prob))
                {
                    log_prop_prob = std::numeric_limits<double>::min();
                    //log_w[n][0] =
                }
                log_w[n][0] = log_order_scores[n][0] - std::log(n + 1) - log_order_scores[n - 1][0] - (log_prop_prob - std::log(p - n));
                //std::cout << "w " << log_w[n][0] << std::endl;
                //assert(std::isfinite(log_prop_prob2));
                // TODO: BUG: Update log node scores!
            }
        }

        if (i > 0)
        // Update the rest of the particles.
        //for (size_t i = 1; i < N; i++)
        {
            //std::cout << "\ni: " << i << std::endl;
            if (n == 0)
            {
                // This is the initialisation
                node_index = rand_int_in_range(n, p - 1); // Draw random node

                if (node_index != 0)
                {
                    move_element(orders[n][i], node_index, n); // Move to node_index to index n=0
                }
                std::vector<double> *sc = scoring.score(orders[n][i], n, 1); // OK, score only index 0
                log_node_scores[n][i] = *sc;
                delete sc;
                log_order_scores[n][i] = log_node_scores[n][i][orders[n][i][n]]; // std::accumulate(log_node_scores[i].begin(), log_node_scores[i].end(), 0.0)
                log_w[n][i] = log_order_scores[n][i] - (-std::log(p));           // First weight
            }
            else
            {
                node_index = rand_int_in_range(n, p - 1); // Draw one of the remaining nodes
                double log_prop_prob = score_sub_order_neigh(scoring,
                                                             orders[n - 1][I[i]],
                                                             log_node_scores[n - 1][I[i]],
                                                             log_order_scores[n - 1][I[i]],
                                                             orders[n][i],
                                                             log_node_scores[n][i],
                                                             log_order_scores[n][i],
                                                             n,
                                                             node_index,
                                                             generator); // Propose from neighborhood.

                //std::cout << "order score " << log_order_scores[n][i] << std::endl;
                log_w[n][i] = log_order_scores[n][i] - std::log(n + 1) - log_order_scores[n - 1][I[i]] - (log_prop_prob - std::log(p - n));
            }
        }
    }
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

        //std::cout << "\n\nn: " << n << " cache hits: " << scoring.cache_hits << " pgibbs order (first " << n + 1 << " elements) " << std::endl;
        //PrintVector(orders[n][0]);
        // Compute the conditional orders for particle Gibbs.
        //std::vector<std::thread> threads(N);

        for (size_t i = 0; i < N; i++)
        {
            pool.push_task(smc_cond_kernel,
                           n,
                           i,
                           std::ref(I),
                           std::ref(log_w),
                           std::ref(orders),
                           std::ref(log_node_scores),
                           std::ref(log_order_scores),
                           std::ref(new_node_inds_cond_orders),
                           std::ref(generator),
                           std::ref(scoring));
        }
        pool.wait_for_tasks();

        // int Num_Threads = std::thread::hardware_concurrency();
        // std::vector<std::thread> threads(Num_Threads);
        // //std::cout << "num threads" << Num_Threads << std::endl;
        // for (size_t i = 0; i < Num_Threads; i++)
        // {
        //     int i_from = i * (N / Num_Threads);
        //     int i_to = (i + 1) * (N / Num_Threads);
        //     //std::cout << i_from << "-" << i_to << std::endl;

        //     if (i == Num_Threads - 1)
        //     {
        //         i_to = N;
        //     }

        //     pool.push_task(smc_cond_kernel,
        //                     n,
        //                     i_from,
        //                     i_to,
        //                     std::ref(I),
        //                     std::ref(log_w),
        //                     std::ref(orders),
        //                     std::ref(log_node_scores),
        //                     std::ref(log_order_scores),
        //                     std::ref(new_node_inds_cond_orders),
        //                     std::ref(generator),
        //                     std::ref(scoring));
        // threads[i] = std::thread(smc_cond_kernel,
        //                          n,
        //                          i_from,
        //                          i_to,
        //                          std::ref(I),
        //                          std::ref(log_w),
        //                          std::ref(orders),
        //                          std::ref(log_node_scores),
        //                          std::ref(log_order_scores),
        //                          std::ref(new_node_inds_cond_orders),
        //                          std::ref(generator),
        //                          std::ref(scoring));
        //}
        //pool.wait_for_tasks();
        // for (size_t i = 0; i < Num_Threads; i++)
        // {
        //     threads[i].join();
        // }

        // smc_cond_kernel(n,
        //                 0,
        //                 N,
        //                 I,
        //                 log_w,
        //                 orders,
        //                 log_node_scores,
        //                 log_order_scores,
        //                 new_node_inds_cond_orders,
        //                 generator,
        //                 scoring);

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
        //std::cout << "Order sampled from conditinal SMC " << std::endl;

        //std::cout << "smc_cond_log_scores " << std::endl;
        //PrintVector(smc_cond_log_scores);
        //std::cout << "smc_cond_log_w " << std::endl;
        //PrintVector(smc_cond_log_w);
        // Sample order from log_w and put in orders.
        norm_w = dist_from_logprobs(smc_cond_log_w);

        std::discrete_distribution<int> distribution(norm_w->begin(), norm_w->end());
        sampled_index = distribution(generator);
        //std::cout << "norm_w " << std::endl;
        //PrintVector(*norm_w);

        delete norm_w;
        orders[j] = smc_cond_orders[sampled_index];
        log_scores[j] = smc_cond_log_scores[sampled_index];
        std::cout << "PGibbs sample score " << log_scores[j] << std::endl;
    }

    return (std::make_pair(orders, log_scores));
}

OrderScoring get_score(Rcpp::List ret)
{

    // Read order
    std::vector<int> order = Rcpp::as<std::vector<int>>(ret["order"]);

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
        for (std::size_t j = 0; j < m.rows(); j++)
        {
            for (std::size_t k = 0; k < m.cols(); k++)
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
    for (int i = 0; i < p; ++i)
    {
        scorepositions[i] = i;
    }

    std::map<cache_keytype3, std::vector<double>> cache;
    OrderScoring scoring(aliases,
                         numparents,
                         rowmaps_backwards,
                         plus1listsparents,
                         scoretable,
                         bannedscore, cache);

    return (scoring);
}

// [[Rcpp::plugins(cpp17)]]

// [[Rcpp::export]]
Rcpp::List r_pgibbs(Rcpp::List ret, int N, int M)
{
    OrderScoring scoring = get_score(ret);

    int seed = 2;
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
