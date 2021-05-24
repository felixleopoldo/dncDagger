#include <bits/stdc++.h>
#include <RInside.h>
#include <cassert>
#include <chrono>
using namespace std::chrono;
typedef std::tuple<std::vector<int>, int, int> cache_keytype;
typedef std::pair<std::vector<int>, std::set<int>> cache_keytype2;
typedef std::pair<std::vector<int>, int> cache_keytype3;

template <typename T>
void PrintVector(std::vector<T> &arr)
{
    copy(arr.begin(), arr.end(),
         std::ostream_iterator<T>(std::cout, ", "));
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
    std::vector<int> scorepositions;
    std::vector<std::vector<int>> potential_parents;
    std::vector<int> numparents;
    std::vector<Rcpp::IntegerVector> rowmaps_backwards;
    std::vector<std::vector<int>> potential_plus1_parents;
    std::vector<std::vector<std::vector<double>>> scoretable;
    std::vector<std::vector<std::vector<double>>> scoresmatrices;
    std::vector<int> permy;

public:
    std::map<cache_keytype3, std::vector<double>> cache;
    int cache_hits;
    OrderScoring(
        std::vector<std::vector<int>> potential_parents,
        std::vector<int> numparents,
        std::vector<Rcpp::IntegerVector> rowmaps_backwards,
        std::vector<std::vector<int>> potential_plus1_parents,
        std::vector<std::vector<std::vector<double>>> scoretable,
        std::vector<std::vector<std::vector<double>>> scoresmatrices,
        std::map<cache_keytype3, std::vector<double>> cache) : potential_parents(potential_parents), numparents(numparents),
                                                               rowmaps_backwards(rowmaps_backwards),
                                                               potential_plus1_parents(potential_plus1_parents),
                                                               scoretable(scoretable), scoresmatrices(scoresmatrices),
                                                               cache(cache), cache_hits(0)
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
                //double max = std::max(node_scores[node_a], nodeb_as_plus1_score);                                                 // I should use order_scores not node scores here. Or something. This is not right at least.
                //std::cout << nodeb_as_plus1_score - node_scores[node_a] << " " << std::endl;
                if(std::abs(nodeb_as_plus1_score - node_scores[node_a]) > 0) { // -10 is arbitrary
                    // OK
                    //std::cout << "NO RECOMPUTE order score for node" << std::endl;
                    nodea_score = std::log(std::exp(node_scores[node_a]) - std::exp(nodeb_as_plus1_score)); // gets inf... 0 at node_scores[node_a] but something at node_scores[node_b]
                } else {
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
                std::swap(ordering[nodea_index], ordering[nodeb_index]);                                // Swap back.
                double max = std::max(node_scores[node_b], nodea_as_plus1_score);
                
                if(std::abs(nodea_as_plus1_score - node_scores[node_b]) > 0) { // -10 is arbitrary
                    nodeb_score = std::log(std::exp(node_scores[node_b] - max) + std::exp(nodea_as_plus1_score - max)) + max;
                } else {
                    std::swap(ordering[nodea_index], ordering[nodeb_index]);
                    nodeb_score = score_pos(ordering, nodea_index);
                    std::swap(ordering[nodea_index], ordering[nodeb_index]);

                }
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
 *  Assets index_to < index_from.
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
    std::vector<double> *norm_prob = new std::vector<double>(n);
    double sum_exp = 0.0;
    for (std::size_t i = 0; i != n; ++i)
    {
        probs_rescaled[i] = std::exp(log_probs[i] - max);
        sum_exp += probs_rescaled[i];
    }

    for (std::size_t i = 0; i != n; ++i)
    {
        (*norm_prob)[i] = probs_rescaled[i] / sum_exp;
    }

    return (norm_prob);
}
/**
 * Score neigborhood of ordering.
 * when the element at index index_of_el_to_insert is inserted at positions 0,..,first_n_elements.
 */
//std::tuple<std::vector<int>, double, std::vector<double>, double> score_sub_order_neigh(const OrderScoring &scoring,
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

//std::tuple<std::vector<int>, double, std::vector<double>, double> score_sub_order_neigh(const OrderScoring &scoring,
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
    double node1_score;
    double node2_score;
    int node = order[index_of_el_to_insert];

    // Calculate score for node at index first_n_elements. The other nodes have the
    // same score sine those are still befor in the ordering.
    move_element(order, index_of_el_to_insert, first_n_elements); // put it in the back
    node_scores[node] = scoring.score_pos(order, first_n_elements);
    order_scores[first_n_elements] = order_score + node_scores[node];

    for (int i = first_n_elements; i > 0; --i) // put it in the back instead
    {
        //PrintVector(order);
        int node1 = order[i - 1];
        int node2 = order[i];

        //std::vector<int> ordertmp(order);
        //std::vector<double> node_scorestmp(node_scores);
        //const auto& [node1_scoretmp, node2_scoretmp] = scoring.swap_nodes(i - 1, i, ordertmp, node_scorestmp); // This should update the orderscore directly maybe.
        const auto& [node1_scoretmp, node2_scoretmp] = scoring.swap_nodes(i - 1, i, order, node_scores); // This should update the orderscore directly maybe.
        node1_score = node1_scoretmp;
        node2_score = node2_scoretmp;


        // std::swap(order[i - 1], order[i]);
        // node1_score = scoring.score_pos(order, i);     
        // node2_score = scoring.score_pos(order, i - 1); 


        double node1_score_diff = node1_score - node_scores[node1];
        double node2_score_diff = node2_score - node_scores[node2];
        order_scores[i - 1] = order_scores[i] + node1_score_diff + node2_score_diff; //- (node_scores[node1] + node_scores[node2]) + (node1_score + node2_score);        
        node_scores[node1] = node1_score;
        node_scores[node2] = node2_score;

        // std::cout << "node1 (" << node1 << ") new " << node1_scoretmp << " true: " << node_scores[node1] << " node2 (" << node2 << ") new " << node2_scoretmp << " true " << node_scores[node2] << std::endl;
        // assert(std::abs(node1_scoretmp - node_scores[node1]) < 0.001);
        // assert(std::abs(node2_scoretmp - node_scores[node2]) < 0.001);
        // assert(std::abs(node2_score_tmp - node2_score)<0.001);
    }

    // order now has node as index 0
    // Sample the new ordering from the order_scores distribution, where the index
    // specifies the position were node is inserted.
    std::vector<double> *neig_dist = dist_from_logprobs(order_scores);                  // Sample one order in the neighborhood
    std::discrete_distribution<int> distribution(neig_dist->begin(), neig_dist->end()); // Create a distribution. O(p+1)
    int insert_pos = distribution(generator);

    // Loop to the to get the correct node scorings. Ths is the fastest way I can think of right now.
    // Caching the node scores is O(p) for each iteration.
    output_order = input_order; // copy
    output_order_score = order_scores[insert_pos];
    output_node_scores = input_node_scores;                              //copy
    move_element(output_order, index_of_el_to_insert, first_n_elements); // put it in the back
    output_node_scores[node] = scoring.score_pos(output_order, first_n_elements);

    for (int i = first_n_elements; i > insert_pos; --i)
    {
        int node1 = output_order[i - 1];
        int node2 = output_order[i];
        std::swap(output_order[i - 1], output_order[i]);
        node1_score = scoring.score_pos(output_order, i);
        node2_score = scoring.score_pos(output_order, i - 1);
        output_node_scores[node1] = node1_score;
        output_node_scores[node2] = node2_score;
    }

    double log_prob = std::log((*neig_dist)[insert_pos]);

    delete neig_dist;
    return log_prob;
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

int rand_int_in_range(const std::size_t &from, const std::size_t &to)
{
    return (from + (std::rand() % (to - from + 1)));
}

void smc(OrderScoring &scoring, std::size_t N, std::size_t p);
void smc(OrderScoring &scoring, std::size_t N, std::size_t p)
{
    std::cout << "Starting SMC " << std::endl;
    auto start = high_resolution_clock::now();
    std::random_device rd;
    std::mt19937 gen(rd());
    std::default_random_engine generator(1);
    std::vector<std::vector<double>> log_w(p, std::vector<double>(N, 0.0));
    std::vector<std::vector<std::vector<int>>> orders_full(N, std::vector<std::vector<int>>(p, std::vector<int>(p)));
    std::vector<int> I(N);
    int node_index;
    std::vector<int> order(p);
    std::vector<double> neig_scoring;
    std::vector<double> neig_dist;
    int insert_pos;
    double log_prop_prob;
    std::vector<std::vector<std::vector<double>>> log_node_scores_full(N, std::vector<std::vector<double>>(p, std::vector<double>(p, 0.0)));
    std::vector<std::vector<double>> log_order_scores_full(p, std::vector<double>(N, 0.0));
    std::vector<double> *norm_w; // = new(N, 0.0);
    std::discrete_distribution<int> distribution;

    for (size_t i = 0; i < N; i++)
    {
        for (size_t n = 0; n < p; n++)
        {
            for (size_t m = 0; m < p; m++)
            {
                orders_full[i][n][m] = m;
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
                    move_element(orders_full[i][n], node_index, n); // Move to node_index to index n=0
                }
                std::vector<double> *sc = scoring.score(orders_full[i][n], n, 1); // OK, score only index 0
                log_node_scores_full[i][n] = *sc;
                delete sc;
                log_order_scores_full[n][i] = log_node_scores_full[i][n][orders_full[i][n][n]]; // std::accumulate(log_node_scores[i].begin(), log_node_scores[i].end(), 0)
                log_w[n][i] = log_order_scores_full[n][i] + std::log(p);                        // First weight
            }
            else
            {
                node_index = rand_int_in_range(n, p - 1); // Draw one of the remaining nodes
                log_prop_prob = score_sub_order_neigh(scoring,
                                                      orders_full[I[i]][n - 1],
                                                      log_node_scores_full[I[i]][n - 1],
                                                      log_order_scores_full[n - 1][I[i]],
                                                      orders_full[i][n],
                                                      log_node_scores_full[i][n],
                                                      log_order_scores_full[n][i],
                                                      n,
                                                      node_index,
                                                      generator); // Scores the neighbohood

                log_w[n][i] = log_order_scores_full[n][i] - std::log(n + 1) - log_order_scores_full[n - 1][I[i]] - (log_prop_prob - std::log(p - n));
            }
        }

        // rescale weights to
        norm_w = dist_from_logprobs(log_w[n]);
        std::discrete_distribution<int> distribution(norm_w->begin(), norm_w->end());

        // Resample particles
        for (std::size_t i = 0; i < N; i++)
        {
            I[i] = distribution(generator);
        }
    }
    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<milliseconds>(stop - start);

    //int maxElementIndex = std::max_element(log_order_scores.begin(), log_order_scores.end()) - log_order_scores.begin();
    int maxElementIndex = std::max_element(log_order_scores_full[p - 1].begin(), log_order_scores_full[p - 1].end()) - log_order_scores_full[p - 1].begin();

    //PrintVector(log_order_scores_full[p-1]);
    std::map<std::vector<int>, double> orders_probs;
    std::set<std::vector<int>> distinct_orders;

    for (int i = 0; i < N; i++)
    {
        distinct_orders.insert(orders_full[i][p - 1]);
        //distinct_orders.insert(o);
    }

    for (int i = 0; i < N; i++)
    {
        if (orders_probs.count(orders_full[i][p - 1]))
        {
            orders_probs[orders_full[i][p - 1]] += (*norm_w)[i];
        }
        else
        {
            orders_probs[orders_full[i][p - 1]] = (*norm_w)[i];
        }
    }

    for (auto o : distinct_orders)
    {
        PrintVector(o);
        std::vector<double> *scr = scoring.score(o, 0, p);
        double sc = std::accumulate(scr->begin(), scr->end(), 0.0);
        std::cout << orders_probs[o] << " score " << sc << std::endl;
        delete scr;
    }

    std::cout << "number of distinct orders: " << distinct_orders.size() << std::endl;
    std::cout << "index: " << maxElementIndex << std::endl;
    //PrintVector(new_orders[maxElementIndex]);
    PrintVector(orders_full[maxElementIndex][p - 1]);

    //std::cout << "prob: " << orders_probs[new_orders[maxElementIndex]] << std::endl;
    std::cout << "prob: " << orders_probs[orders_full[maxElementIndex][p - 1]] << std::endl;
    //std::cout << "score: " << log_order_scores[maxElementIndex] << std::endl;
    std::cout << "score: " << log_order_scores_full[p - 1][maxElementIndex] << std::endl;
    //PrintVector(log_order_scores_full[p - 1]);
    // To get the value of duration use the count()
    // member function on the duration object
    std::cout << duration.count() << " ms." << std::endl;
    delete norm_w;
}

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