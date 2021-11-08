#include <bits/stdc++.h>
#include <RInside.h>
#include <cassert>
#include <chrono>
#include <thread>
#include <iostream>

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
    // std::map<cache_keytype3, std::vector<double>> cache;
    // int cache_hits;
    OrderScoring(
        std::vector<std::vector<int>> potential_parents,
        std::vector<int> numparents,
        std::vector<Rcpp::IntegerVector> rowmaps_backwards,
        std::vector<std::vector<int>> potential_plus1_parents,
        std::vector<std::vector<std::vector<double>>> scoretable,
        std::vector<std::vector<std::vector<double>>> scoresmatrices,
        bool MAP
        // std::map<cache_keytype3, std::vector<double>> cache
        ) : potential_parents(potential_parents),
            rowmaps_backwards(rowmaps_backwards),
            potential_plus1_parents(potential_plus1_parents),
            scoretable(scoretable),
            scoresmatrices(scoresmatrices),
            MAP(MAP),
            numparents(numparents)
    // cache(cache),
    // cache_hits(0)
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