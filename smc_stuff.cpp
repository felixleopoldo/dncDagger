#include <bits/stdc++.h>
#include <RInside.h>
#include <cassert>

typedef std::tuple<std::vector<int>, int, int> cache_keytype;

class OrderScoring
{
private:
    std::vector<int> scorepositions;
    std::vector<std::vector<int>> aliases;
    std::vector<int> numparents;
    std::vector<Rcpp::IntegerVector> rowmaps_backwards;
    std::vector<Rcpp::IntegerVector> plus1listsparents;
    std::vector<std::vector<std::vector<double>>> scoretable;
    std::vector<Rcpp::NumericMatrix> scoresmatrices;
    std::vector<int> permy;
    std::map<cache_keytype, std::vector<double>> cache;

public:
    OrderScoring(
        std::vector<std::vector<int>> aliases,
        std::vector<int> numparents,
        std::vector<Rcpp::IntegerVector> rowmaps_backwards,
        std::vector<Rcpp::IntegerVector> plus1listsparents,
        std::vector<std::vector<std::vector<double>>> scoretable,
        std::vector<Rcpp::NumericMatrix> scoresmatrices,
        std::map<cache_keytype, std::vector<double>> cache) : aliases(aliases), numparents(numparents),
                                                               rowmaps_backwards(rowmaps_backwards),
                                                               plus1listsparents(plus1listsparents), 
                                                               scoretable(scoretable), scoresmatrices(scoresmatrices),
                                                               cache(cache)

    {
    }

    /**
         * Score elements in scorenodes from scorepositions and n_elemens on.
         */
    std::vector<double> score(const std::vector<int> &scorenodes, const int &from_orderpos, const int &n_elements) 
    {
        // If in cache
        cache_keytype suborder = std::make_tuple(scorenodes, from_orderpos, n_elements);
        if (cache.count(suborder))
        {
            //std::cout << "cache hit" << std::endl;
            return (cache[suborder]);
        }
        std::size_t n = scorenodes.size();
        std::vector<double> orderscores(n, 0.0);            // orderscores <- vector("double", n)
        std::vector<std::vector<int>> allowedscorelists(n); // allowedscorelists < -vector("list", n)
        std::vector<int> therows(n);                        // therows <- vector("integer", n) # what is this? / Felix
        size_t k = 0;                                       // k <- 1

        for (int i = from_orderpos; i < from_orderpos + n_elements; ++i)
        {
            int node = scorenodes[i];
            size_t position = i;

            if (position == n - 1)
            {
                // no parents allowed, i.e.only first row, only first list
                // orderscores[node] <- scoretable [[node]][[1]][1, 1]
                orderscores[node] = scoretable[node][0][0];
                // allowedscorelists [[node]] <- c(1)
                allowedscorelists[node].assign({0});
                // therows[node] <- c(2 ^ numparents[node])
                therows[node] = std::pow(2, numparents[node]);
            }
            else
            {
                //std::vector<int> bannednodes(scorenodes.begin(), scorenodes.begin() + position + 1); // This variable is not necessary
                //std::vector<int> allowednodes(scorenodes.begin() + position + 1, scorenodes.end()); // This variable is not necessary
                std::vector<int> bannedpool;
                for (int j = 0; j < aliases[node].size(); j++)
                {
                    //if (std::find(bannednodes.begin(), bannednodes.end(), aliases[node][j]) != bannednodes.end())
                    if (std::find(scorenodes.begin(), scorenodes.begin() + position + 1, aliases[node][j]) != scorenodes.begin() + position + 1)
                    {
                        bannedpool.push_back(j);
                    }
                }

                if (numparents[node] == 0 || bannedpool.size() == 0)
                {
                    therows[node] = 0; // all parents allowed
                }
                else
                {
                    int indextmp = 0;
                    for (auto &item : bannedpool)
                    {
                        indextmp += std::pow(2, item + 1); // add 1 since nodes are labeled from 0
                    }
                    indextmp = indextmp / 2;
                    // this might be list I guess. rowmaps_backwards[node][std::sum(2 ^ bannedpool) / 2 + 1]
                    therows[node] = rowmaps_backwards[node][indextmp];
                }

                // allowedscorelists [[node]] <- c(1, which(plus1lists$parents [[node]] % in % allowednodes) + 1)
                allowedscorelists[node].push_back(0);
                // Optimize when optimizing order.
                for (int j = 0; j < plus1listsparents[node].size(); j++)
                {
                    //if (std::find(allowednodes.begin(), allowednodes.end(), plus1listsparents[node][j]) != allowednodes.end())
                    if (std::find(scorenodes.begin() + position + 1, scorenodes.end(), plus1listsparents[node][j]) != scorenodes.end())
                    {
                        allowedscorelists[node].push_back(j + 1);
                    }
                }

                //std::vector<double> scoresvec;
                std::vector<double> scoresvec(allowedscorelists[node].size());
                for (int j = 0; j < scoresvec.size(); j++)
                //for (int allowedscore : allowedscorelists[node])
                {
                    //scoresvec.push_back(scoresmatrices[node](therows[node], allowedscore));
                    scoresvec[j] = scoresmatrices[node](therows[node], allowedscorelists[node][j]);
                }

                double maxallowed = *std::max_element(scoresvec.begin(), scoresvec.end()); // maxallowed <- max(scoresvec)

                for (auto &val : scoresvec)
                {
                    orderscores[node] += std::exp(val - maxallowed);
                }
                orderscores[node] = maxallowed + std::log(orderscores[node]);
            }
            k++;
        }
        cache[suborder] = orderscores;
        return (orderscores);
    }
};

std::vector<double> orderscorePlus1(std::size_t n,
                                    std::vector<int> &scorenodes,
                                    std::vector<int> &scorepositions,
                                    std::vector<std::vector<int>> &aliases,
                                    std::vector<int> &numparents,
                                    std::vector<Rcpp::IntegerVector> &rowmaps_backwards,
                                    std::vector<Rcpp::IntegerVector> &plus1listsparents,
                                    std::vector<std::vector<std::vector<double>>> &scoretable,
                                    std::vector<Rcpp::NumericMatrix> &scoresmatrices,
                                    std::vector<int> &permy);

std::vector<double> orderscorePlus1(std::size_t n,
                                    std::vector<int> &scorenodes,
                                    std::vector<int> &scorepositions,
                                    std::vector<std::vector<int>> &aliases,
                                    std::vector<int> &numparents,
                                    std::vector<Rcpp::IntegerVector> &rowmaps_backwards,
                                    std::vector<Rcpp::IntegerVector> &plus1listsparents,
                                    std::vector<std::vector<std::vector<double>>> &scoretable,
                                    std::vector<Rcpp::NumericMatrix> &scoresmatrices,
                                    std::vector<int> &permy)
{

    std::vector<double> orderscores(n);                 // orderscores <- vector("double", n)
    std::vector<std::vector<int>> allowedscorelists(n); // allowedscorelists < -vector("list", n)
    std::vector<int> therows(n);                        // therows <- vector("integer", n) # what is this? / Felix
    size_t k = 0;                                       // k <- 1

    for (int i : scorenodes)
    {
        std::cout << "node " << i << std::endl;
        size_t position = scorepositions[k];
        std::cout << "position " << position << std::endl;
        if (position == n - 1)
        {
            // no parents allowed, i.e.only first row, only first list
            orderscores[i] = scoretable[i][0][0];    // orderscores[i] <- scoretable [[i]][[1]][1, 1]
            allowedscorelists[i].assign({0});        // allowedscorelists [[i]] <- c(1)
            therows[i] = std::pow(2, numparents[i]); // therows[i] <- c(2 ^ numparents[i])
        }
        else
        {
            std::vector<int> bannednodes(permy.begin(), permy.begin() + position + 1); // bannednodes <- permy[1:position]
            std::vector<int> allowednodes(permy.begin() + position + 1, permy.end());  // allowednodes < -permy [(position + 1):n]
            std::vector<int> bannedpool;                                               // bannedpool <- which(aliases [[i]] % in % bannednodes)
            for (int j = 0; j < aliases[i].size(); j++)
            {
                if (std::find(bannednodes.begin(), bannednodes.end(), aliases[i][j]) != bannednodes.end())
                {
                    bannedpool.push_back(j);
                }
            }

            if (numparents[i] == 0 || bannedpool.size() == 0)
            {
                therows[i] = 0; // all parents allowed
            }
            else
            {
                int indextmp = 0;
                for (auto &item : bannedpool)
                {
                    indextmp += std::pow(2, item + 1); // add 1 since nodes are labeled from 0
                }
                indextmp = indextmp / 2;
                therows[i] = rowmaps_backwards[i][indextmp]; // this might be list I guess. rowmaps_backwards[i][std::sum(2 ^ bannedpool) / 2 + 1]
            }

            allowedscorelists[i].push_back(0); // allowedscorelists [[i]] <- c(1, which(plus1lists$parents [[i]] % in % allowednodes) + 1)
            // Optimize when optimizing order.
            for (int j = 0; j < plus1listsparents[i].size(); j++)
            {
                if (std::find(allowednodes.begin(), allowednodes.end(), plus1listsparents[i][j]) != allowednodes.end())
                {
                    allowedscorelists[i].push_back(j + 1);
                }
            }

            std::vector<double> scoresvec;
            for (int allowedscore : allowedscorelists[i])
            {
                scoresvec.push_back(scoresmatrices[i](therows[i], allowedscore));
            }

            double maxallowed = *std::max_element(scoresvec.begin(), scoresvec.end()); // maxallowed <- max(scoresvec)

            for (auto &val : scoresvec)
            {
                orderscores[i] += std::exp(val - maxallowed);
            }
            orderscores[i] = maxallowed + std::log(orderscores[i]);
        }
        k++;
    }
    //scores < -list()
    //scores$therow < -therows
    //scores$allowedlists < -allowedscorelists
    //scores$totscores < -orderscores
    return (orderscores);
}

template <typename T>
void PrintVector(std::vector<T> &arr)
{
    copy(arr.begin(), arr.end(),
         std::ostream_iterator<T>(std::cout, ", "));
    std::cout << std::endl;
}

/**
 *  Assets index_to < index_from.
 */
void move_element(std::vector<int> &v, std::size_t index_from, std::size_t index_to);

void move_element(std::vector<int> &v, std::size_t index_from, std::size_t index_to)
{
    int element_to_move = v[index_from];
    v.erase(v.begin() + index_from);
    v.insert(v.begin() + index_to, element_to_move);
    return;
}

std::vector<double> dist_from_logprobs(std::vector<double> &log_probs)
{
    int n = log_probs.size();
    double max = *std::max_element(log_probs.begin(), log_probs.end());
    std::vector<double> log_probs_rescaled(n, 0.0);

    double sum_exp = 0.0;
    for (std::size_t i = 0; i != n; ++i)
    {
        log_probs_rescaled[i] = log_probs[i] - max;
        sum_exp += std::exp(log_probs_rescaled[i]);
    }

    std::vector<double> norm_prob(n);
    for (std::size_t i = 0; i != n; ++i)
    {
        norm_prob[i] = std::exp(log_probs_rescaled[i]) / sum_exp;
    }

    return (norm_prob);
}
/**
 * Score neigborhood of ordering.
 * when the element at index index_of_el_to_insert is inserted at positions 0,..,first_n_elements.
 */
std::tuple<std::vector<int>, double, std::vector<double>, double> score_sub_order_neigh( OrderScoring &scoring,
                                                                                        const std::vector<int> &input_order,
                                                                                        const std::vector<double> &nodes_scores,
                                                                                        double input_order_score,
                                                                                        int first_n_elements,
                                                                                        int index_of_el_to_insert, std::default_random_engine generator);

std::tuple<std::vector<int>, double, std::vector<double>, double> score_sub_order_neigh( OrderScoring &scoring,
                                                                                        const std::vector<int> &input_order,
                                                                                        const std::vector<double> &input_node_scores,
                                                                                        double input_order_score,
                                                                                        int first_n_elements,
                                                                                        int index_of_el_to_insert, std::default_random_engine generator)
{
    std::vector<double> order_scores(first_n_elements + 1);
    std::vector<double> new_node_scores;
    std::vector<double> node_scores_tmp;
    std::vector<double> node_scores(input_node_scores);
    std::vector<int> order(input_order);
    double order_score(input_order_score);

    //cache_keytype neig_cache_key = std::make_tuple(input_order, input_order_score, index_of_el_to_insert);
    //typedef std::tuple<std::discrete_distribution<int>, std::vector<std::vector<double> > > dist_nodescores;
    //std::map<cache_keytype, dist_nodescores > neig_cache;
    // if(neig_cache.count(neig_cache_key)){
    //     int pos = neig_cache[neig_cache_key].distribution(generator);
    //     node_scores = neig_cache[neig_cache_key][pos].nodescores[]
    //     return (std::make_tuple(order, std::log(neig_dist[insert_pos]), node_scores, order_scores[insert_pos]));
    // }

    int node = order[index_of_el_to_insert];
    move_element(order, index_of_el_to_insert, first_n_elements); // put it in the back
    // Calculate score for node at index first_n_elements. The other nodes have the
    // same score sine those are still befor in the ordering.
    new_node_scores = scoring.score(order, first_n_elements, 1); // new element is on index first_n_elements
    node_scores[node] = new_node_scores[node];                   // update the score for node
    order_scores[first_n_elements] = order_score + node_scores[node];

    for (int i = first_n_elements; i > 0; --i) // put it in the back instead
    {
        int node1 = order[i - 1];
        int node2 = order[i];
        std::swap(order[i - 1], order[i]);
        new_node_scores = scoring.score(order, i - 1, 2); // Only need to rescore the swapped nodes
        order_scores[i - 1] = order_scores[i] - node_scores[node1] - node_scores[node2] + new_node_scores[node1] + new_node_scores[node2];
        node_scores[node1] = new_node_scores[node1];
        node_scores[node2] = new_node_scores[node2];

        //node_scores_tmp = scoring.score(order, 0, first_n_elements + 1);
        //double order_score_check = std::accumulate(node_scores_tmp.begin(), node_scores_tmp.end(), 0.0);
        //assert(std::abs(order_score_check - order_scores[i - 1]) < 0.00001 );
    }
    // order now has node as index 0

    // Sample the new ordering from the order_scores distribution, where the index
    // specifies the position were node is inserted.

    std::vector<double> neig_dist = dist_from_logprobs(order_scores);                 // Sample one order in the neighborhood
    std::discrete_distribution<int> distribution(neig_dist.begin(), neig_dist.end()); // Create a distribution. O(p+1)
    int insert_pos = distribution(generator);

    // Loop to the to get the correct node scorings. Ths is the fastest way I can think of right now.
    // Caching the node scores is O(p) for each iteration.

    order = input_order;
    node_scores = input_node_scores;

    move_element(order, index_of_el_to_insert, first_n_elements); // put it in the back
    new_node_scores = scoring.score(order, first_n_elements, 1);  // new element is on index first_n_elements
    node_scores[node] = new_node_scores[node];                    // update the score for node

    for (int i = first_n_elements; i > insert_pos; --i)
    {
        int node1 = order[i - 1];
        int node2 = order[i];
        std::swap(order[i - 1], order[i]);
        new_node_scores = scoring.score(order, i - 1, 2); // Only need to rescore the swapped nodes
        node_scores[node1] = new_node_scores[node1];
        node_scores[node2] = new_node_scores[node2];
    }

    return (std::make_tuple(order, std::log(neig_dist[insert_pos]), node_scores, order_scores[insert_pos]));
}

// std::tuple<std::vector<int>, double, double> proposal(const OrderScoring &scoring,
//                                           std::vector<int> input_order,
//                                           std::vector<double> node_scores,
//                                           double order_score,
//                                           int first_n_elements, int index_of_el_to_insert, std::default_random_engine generator)
// {

//     std::vector<double> order_scores = score_sub_order_neigh(scoring, input_order, node_scores, order_score, int_)
//     std::vector<double> neig_dist = dist_from_logprobs(order_scores);                                     // Sample one order in the neighborhood
//     std::discrete_distribution<int> distribution(neig_dist.begin(), neig_dist.end()); // Create a distribution. O(p+1)
//     int insert_pos = distribution(generator);                                             // O(N)

//     // move mode to index insert_pos
//     std::vector<int> ret_order(input_order);

//     move_element(ret_order, index_of_el_to_insert, insert_pos); // Move the element to the selectd position
//     // TODO: Need to return log_node_scores
//     return (std::make_tuple(ret_order, std::log(neig_dist[insert_pos]), order_scores[insert_pos]));
// }

int rand_int_in_range(std::size_t from, std::size_t to)
{
    return (from + (std::rand() % (to - from + 1)));
}

void smc( OrderScoring &scoring, std::size_t N, std::size_t p);
void smc( OrderScoring &scoring, std::size_t N, std::size_t p)
{
    std::random_device rd;
    std::mt19937 gen(rd());
    std::default_random_engine generator(1);
    std::vector<std::vector<double>> log_w(p, std::vector<double>(N, 0.0));
    std::vector<std::vector<int>> new_orders(N);
    std::vector<std::vector<int>> old_orders(N);

    for (size_t i = 0; i < N; i++)
    {
        new_orders[i] = std::vector<int>(p);
        for (size_t n = 0; n < p; n++)
        {
            new_orders[i][n] = n;
        }
    }

    std::cout << "Starting SMC " << std::endl;

    std::vector<int> I(N);
    int node_index;
    std::vector<int> order(p);
    std::vector<double> neig_scoring;
    std::vector<double> neig_dist;
    int insert_pos;

    std::vector<std::vector<double>> log_node_scores(N, std::vector<double>(p, 0.0));
    std::vector<std::vector<double>> log_node_scores_prev(N, std::vector<double>(p, 0.0));
    std::vector<double> log_order_scores(N, 0.0);
    std::vector<double> log_order_scores_prev(N, 0.0);

    std::vector<double> log_node_score(p);

    double order_score;
    std::vector<double> norm_w(N, 0.0);
    std::discrete_distribution<int> distribution;

    for (size_t n = 0; n < p; n++)
    {
        //std::cout << "\n\nn: " << n << std::endl;
        for (size_t i = 0; i < N; i++)
        {
            //std::cout << "\ni: " << i << std::endl;
            if (n == 0)
            {
                // This is the initialisation
                node_index = rand_int_in_range(n, p - 1); // Draw random node
                if (node_index != 0)
                {
                    move_element(new_orders[i], node_index, n); // Move to node_index to index n=0
                }
                //PrintVector(new_orders[i]);
                log_node_scores[i] = scoring.score(new_orders[i], n, 1);    // OK, score oncly index 0
                log_order_scores[i] = log_node_scores[i][new_orders[i][n]]; //std::accumulate(log_node_scores[i].begin(), log_node_scores[i].end(), 0)
                log_w[n][i] = log_order_scores[i] - std::log(p);            // First weight
            }
            else
            {
                node_index = rand_int_in_range(n, p - 1); // Draw one of the remaining nodes
                auto [prop_order, log_prop_prob, log_node_score, log_order_score] = score_sub_order_neigh(scoring,
                                                                                                          old_orders[I[i]],            // O(p)
                                                                                                          log_node_scores_prev[I[i]],  // O(p)
                                                                                                          log_order_scores_prev[I[i]], // O(1)
                                                                                                          n,
                                                                                                          node_index, generator); // Scores the neighbohood
                log_node_scores[i] = log_node_score;
                log_order_scores[i] = log_order_score;
                log_w[n][i] = log_order_scores[i] - log_order_scores_prev[I[i]] - (log_prop_prob - std::log(p - n));
                new_orders[i] = prop_order;
            }
        }

        // rescale weights to
        norm_w = dist_from_logprobs(log_w[n]);
        std::discrete_distribution<int> distribution(norm_w.begin(), norm_w.end());

        // Resample particles
        for (std::size_t i = 0; i < N; i++)
        {
            I[i] = distribution(generator);
        }

        old_orders = new_orders; // Copy in some way
        log_node_scores_prev = log_node_scores;
        log_order_scores_prev = log_order_scores;
    }

    int maxElementIndex = std::max_element(log_order_scores.begin(), log_order_scores.end()) - log_order_scores.begin();

    std::map<std::vector<int>, double> orders_probs;
    std::set<std::vector<int>> distinct_orders;

    for (auto o : new_orders)
    {
        distinct_orders.insert(o);
    }

    for (int i = 0; i < N; i++)
    {
        if (orders_probs.count(new_orders[i]))
        {
            orders_probs[new_orders[i]] += norm_w[i];
        }
        else
        {
            //std::cout << orders_probs[new_orders[i]] << std::endl;
            orders_probs[new_orders[i]] = norm_w[i];
        }
    }

    for (auto o : distinct_orders)
    {
        PrintVector(o);
        std::cout << orders_probs[o] << std::endl;
    }

    std::cout << "number of distinct orders: " << distinct_orders.size() << std::endl;
    std::cout << "index: " << maxElementIndex << std::endl;
    PrintVector(new_orders[maxElementIndex]);
    std::cout << "prob: " << orders_probs[new_orders[maxElementIndex]] << std::endl;
    std::cout << "score: " << log_order_scores[maxElementIndex] << std::endl;
}
