
#include <bits/stdc++.h>
#include <RInside.h>
#include <cassert>
#include <chrono>
#include <thread>
#include <iostream>
#include "thread_pool.hpp"
// #include <cstdio>
// #include <algorithm>

using namespace std::chrono;

double EPSILON = 0.00001;

inline std::vector<double> *dist_from_logprobs(const std::vector<double> &log_probs)
{
    std::size_t n = log_probs.size();
    double max = *std::max_element(log_probs.begin(), log_probs.end());
    std::vector<double> probs_rescaled(n);
    // std::vector<double> *probs_rescaled = new std::vector<double>(n);
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
}

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

    // std::cout << prop_order_score << std::endl;
    // std::cout << cur_order_score << std::endl;

    double ratio = std::exp(prop_order_score - cur_order_score);
    double alpha = std::min(1.0, ratio);
    double a = (float)std::rand() / RAND_MAX;

    // std::cout << "rand number: " << a << std::endl;
    // std::cout << "accept prob: " << alpha << std::endl;
    if (a < alpha)
    {
        // std::cout << "Accept" << std::endl;
        node_scores[node1] = node1_scoretmp;
        node_scores[node2] = node2_scoretmp;
        return (std::make_pair(ordering, prop_order_score));
    }
    else
    {
        // std::cout << "Reject" << std::endl;
        //  swap back
        myswap(i - 1, i, ordering);
        return (std::make_pair(ordering, cur_order_score));
    }
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
    std::vector<double> sc = scoring.score(output_order, n, 1); // OK, score only index 0
    // output_node_scores = *sc;
    // delete sc;
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
        // std::cout << "\n\nn: " << n << " cache hits: " << scoring.cache_hits << std::endl;
        for (size_t i = 0; i < N; i++)
        {
            // std::cout << "\ni: " << i << std::endl;
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
            std::vector<double> sc = scoring.score(orders[n][0], 0, n + 1); // Score all the nodes in the sub order [0,..,n].
            log_node_scores[n][0] = sc;
            // delete sc;
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

    // for node at index first_n_elements. The other nodes have the
    //  same score since those are still before in the ordering.
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
    // std::cout << "Scores to be dist" << std::endl;
    // PrintVector(order_scores);
    // std::cout << "Score for index "<< new_node_index <<": " << order_scores[new_node_index] << std::endl;
    std::vector<double> *neig_dist = dist_from_logprobs(order_scores); // TODO: dont need dist, just logs
    // PrintVector(*neig_dist);
    // std::cout << "prop prob" << (*neig_dist)[new_node_index] << std::endl;
    return std::log((*neig_dist)[new_node_index]); // NOTE: This can be -Inf sometimes.
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
    // std::cout << "Starting conditional SMC " << std::endl;
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
        // std::cout << "num threads" << Num_Threads << std::endl;

        for (size_t i = 0; i < Num_Threads; i++)
        {
            int i_from = i * (N / Num_Threads);
            int i_to = (i + 1) * (N / Num_Threads);
            // std::cout << i_from << "-" << i_to << std::endl;

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

        // PrintVector(log_w[n]);
        //  rescale weights to
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
    // orders[0] = [6, 7, 0, 5, 1, 3, 4, 2]

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
    // std::vector<std::vector<int>> orders(M, std::vector<int>(p));
    // std::vector<std::vector<int>> orders(2, std::vector<int>(2));
    std::vector<double> log_scores(M);

    // Init vector

    std::vector<int> order = std::vector<int>(p);
    // orders[0] = std::vector<int>(p);
    for (std::size_t i = 0; i < p; i++)
    {
        // orders[0][i] = i;
        order[i] = i;
    }

    // std::vector<double> *log_nodes_scores = scoring.score(orders[0], 0, p);
    std::vector<double> log_nodes_scores = scoring.score(order, 0, p);
    log_scores[0] = std::accumulate(log_nodes_scores.begin(), log_nodes_scores.end(), 0.0);
    double max_score = -INFINITY;
    std::vector<int> max_order;

    for (std::size_t j = 1; j < M; j++)
    {
        // const auto &[mh_order, mh_log_score] = mh_swap_move(orders[j - 1], *log_nodes_scores, log_scores[j - 1], scoring);
        const auto &[mh_order, mh_log_score] = mh_swap_move(order, log_nodes_scores, log_scores[j - 1], scoring);
        // orders[j] = mh_order;
        order = mh_order;
        log_scores[j] = mh_log_score;
        if (mh_log_score > max_score)
        {
            max_score = mh_log_score;
            max_order = order;
        }
    }

    // delete log_nodes_scores;

    return (std::make_pair(max_order, log_scores));
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
    // pool.push_task(task, arg1, arg2);
    // pool.wait_for_tasks();

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
