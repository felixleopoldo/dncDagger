#include <bits/stdc++.h>
#include <RInside.h>

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
    int n;

public:
    OrderScoring(
        std::vector<std::vector<int>> aliases,
        std::vector<int> numparents,
        std::vector<Rcpp::IntegerVector> rowmaps_backwards,
        std::vector<Rcpp::IntegerVector> plus1listsparents,
        std::vector<std::vector<std::vector<double>>> scoretable,
        std::vector<Rcpp::NumericMatrix> scoresmatrices,
        std::vector<int> permy) : aliases(aliases), numparents(numparents), rowmaps_backwards(rowmaps_backwards),
                                  plus1listsparents(plus1listsparents), scoretable(scoretable), scoresmatrices(scoresmatrices),
                                  permy(permy)
    {
        n = permy.size();
    }

    /**
         * Score elements in scorenodes from scorepositions and n_elemens on.
         */
    std::vector<double> score(std::vector<int> &scorenodes, int from_orderpos, int n_elements)
    {
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
                orderscores[node] = scoretable[node][0][0];    // orderscores[node] <- scoretable [[node]][[1]][1, 1]
                allowedscorelists[node].assign({0});           // allowedscorelists [[node]] <- c(1)
                therows[node] = std::pow(2, numparents[node]); // therows[node] <- c(2 ^ numparents[node])
            }
            else
            {
                std::vector<int> bannednodes(permy.begin(), permy.begin() + position + 1); // bannednodes <- permy[1:position]
                std::vector<int> allowednodes(permy.begin() + position + 1, permy.end());  // allowednodes < -permy [(position + 1):n]
                std::vector<int> bannedpool;                                               // bannedpool <- which(aliases [[node]] % in % bannednodes)
                for (int j = 0; j < aliases[node].size(); j++)
                {
                    if (std::find(bannednodes.begin(), bannednodes.end(), aliases[node][j]) != bannednodes.end())
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
                    therows[node] = rowmaps_backwards[node][indextmp]; // this might be list I guess. rowmaps_backwards[node][std::sum(2 ^ bannedpool) / 2 + 1]
                }

                allowedscorelists[node].push_back(0); // allowedscorelists [[node]] <- c(1, which(plus1lists$parents [[node]] % in % allowednodes) + 1)
                // Optimize when optimizing order.
                for (int j = 0; j < plus1listsparents[node].size(); j++)
                {
                    if (std::find(allowednodes.begin(), allowednodes.end(), plus1listsparents[node][j]) != allowednodes.end())
                    {
                        allowedscorelists[node].push_back(j + 1);
                    }
                }

                std::vector<double> scoresvec;
                for (int allowedscore : allowedscorelists[node])
                {
                    scoresvec.push_back(scoresmatrices[node](therows[node], allowedscore));
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
        //scores < -list()
        //scores$therow < -therows
        //scores$allowedlists < -allowedscorelists
        //scores$totscores < -orderscores
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

/**
 * Score neigborhood of ordering.
 * when the element at index index_of_el_to_insert is inserted at positions 0,..,first_n_elements.
 */
std::vector<double> score_sub_order_neigh(OrderScoring &scoring, std::vector<int> order, int first_n_elements, int index_of_el_to_insert);

std::vector<double> score_sub_order_neigh(OrderScoring &scoring, std::vector<int> order, int first_n_elements, int index_of_el_to_insert)
{
    std::vector<int> scorepositions(order.size());
    for (int i = 0; i < order.size(); ++i)
    {
        scorepositions[i] = i;
    }

    std::vector<double> scores(first_n_elements + 1);
    std::vector<double> node_scores;

    move_element(order, index_of_el_to_insert, 0);

    node_scores = scoring.score(order, 0, first_n_elements + 1);

    //PrintVector(order);
    //PrintVector(node_scores);
    // TODO: optimize to only sum the needed part.
    scores[0] = std::accumulate(node_scores.begin(), node_scores.end(), 0);

    for (int i = 0; i < first_n_elements ; ++i)
    {
        std::swap(order[i], order[i + 1]);
        node_scores = scoring.score(order, 0, first_n_elements + 1);
        //PrintVector(order);
        //PrintVector(node_scores);
        // TODO: optimize to only sum the needed part.
        scores[i+1] = std::accumulate(node_scores.begin(), node_scores.end(), 0);
    }

    return (scores);
}

std::vector<double> dist_from_logprobs(std::vector<double> &log_probs)
{
    int n = log_probs.size();
    double max = *std::max_element(log_probs.begin(), log_probs.end());
    std::vector<double> log_probs_rescaled(n, 0.0);
    
    double sum_exp = 0.0;
    for (std::size_t i = 0; i !=n; ++i) {
        log_probs_rescaled[i] = log_probs[i] - max;
        sum_exp += std::exp(log_probs_rescaled[i]);
    }

    std::vector<double> norm_prob(n);
    for (std::size_t i = 0; i !=n; ++i) {
        norm_prob[i] = std::exp(log_probs_rescaled[i]) / sum_exp;
    }

    return(norm_prob);
}

int rand_int_in_range(std::size_t from, std::size_t to) {
    return(from + ( std::rand() % ( to - from + 1 ) ));
}


void smc(OrderScoring scoring, std::size_t N, std::size_t p);
void smc(OrderScoring scoring, std::size_t N, std::size_t p)
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

    std::cout << " hohooho "<< std::endl;

    std::vector<int> I(N);
    int node_index;
    std::vector<int> order(p);
    std::vector<double> neig_scoring;
    std::vector<double> neig_dist;
    int insert_pos;
    std::vector<double> log_node_scores;
    double order_score;

    for (size_t n = 0; n < p; n++)
    {

        std::cout << "\n\nn: " << n << std::endl;
        if (n > 0)
        {
            // rescale weights to
            std::vector<double> norm_w = dist_from_logprobs(log_w[n]);
            std::discrete_distribution<int> distribution(norm_w.begin(), norm_w.end());

            // Resample particles
            for (std::size_t i = 0; i < N; i++)
            {
                I[i] = distribution(generator);
            }
        }

        for (size_t i = 0; i < N; i++)
        {
            std::cout << "\ni: " << i << std::endl;
            
            if (n == 0)
            {
                // This is the initialisation
                node_index = rand_int_in_range(n, p - 1);
                if (node_index != 0)
                {
                    move_element(new_orders[i], node_index, n);
                }
                //PrintVector(new_orders[i]);

                log_w[n][i] = 0.0;
            }
            else
            {
                
                node_index = rand_int_in_range(n, p - 1);

                std::vector<int> order(old_orders[I[i]]); // Same every time even though I[i] differ.
                neig_scoring = score_sub_order_neigh(scoring, order, n - 1, node_index);
                std::vector<double> neig_dist = dist_from_logprobs(neig_scoring);
                std::discrete_distribution<int> distribution(neig_dist.begin(), neig_dist.end());

                insert_pos = distribution(generator);

                // PrintVector(order);
                // std::cout << "stemming from particle: " << I[i] << std::endl;
                // PrintVector(neig_dist);
                // std::cout << "move node " << order[node_index] << " to pos " << insert_pos << std::endl;

                move_element(order, node_index, insert_pos);
                //PrintVector(order);
                log_node_scores = scoring.score(order, 0, p); // true unnormalised node posteriors
                order_score = std::accumulate(log_node_scores.begin(), log_node_scores.end(), 0) - std::log(p - n);
                log_w[n][i] = order_score - neig_dist[insert_pos];
                new_orders[i] = order;
            }
        }
        old_orders = new_orders; // Copy in some way
    }

    std::vector<double> last_dist = dist_from_logprobs(log_w[p-1]);
    for (size_t i = 0; i < N; i++) {
        PrintVector(new_orders[i]);
        std::cout << "weight: " << last_dist[i] << std::endl;
    }
    
}
