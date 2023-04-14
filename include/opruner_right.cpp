
#include <cassert>
#include <chrono>
#include <iostream>
#include <Rcpp.h>
#include "auxiliary.h"
#include "OrderScoring.h"
#include "RightOrder.h"

double EPSILON = 0.0000001;

using namespace std;

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
 * Compares the max score when new_node is inserted somewhere to when its put in at the top.
 */
bool optimal_front(const RightOrder &ro,
                   size_t new_node,
                   OrderScoring &scoring)
{

    // cout << ro.inserted_max_order_scores[new_node] << " > " << ro.new_top_scorses[new_node] << "?" << endl;
    if (definitelyGreaterThan(ro.inserted_max_order_scores[new_node], ro.new_top_scores[new_node], EPSILON)) // a > b
    {
        // cout << "YES, so " << new_node << " not optimal at top " << endl;
        return (false);
    }
    // cout << "NO, so " << new_node << "  optimal at top " << endl;
    return (true);
}

/**
 * Checks if the front is optimal where it is or if it fits better somewhere else.
 */
bool optimal_front(const RightOrder &ro,
                   OrderScoring &scoring)
{
    RightOrder ro_tmp(ro);
    int p = ro.order.size();

    for (int i = ro.front_ind(); i < p - 1; ++i)
    {
        swap_nodes(i, i + 1, ro_tmp, scoring);
        if (definitelyGreaterThan(ro_tmp.order_score, ro.order_score, EPSILON)) // this implies worse performance..
        {
            return (false);
        }
    }
    return (true);
}

/**
 *
 */
bool independent_front(const RightOrder &ro,
                       const vector<double> &top_scores,
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
 * NOTE: It's not really visible since n does not change..
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
/// @brief The maximum score when the node at index node_ind is inserted among the visible ones.
/// @param ro RightOrder
/// @param node_ind
/// @param scoring
/// @return
double max_inserted_score(const RightOrder &ro,
                          size_t node_ind,
                          OrderScoring &scoring)
{
    RightOrder ro_tmp(ro);
    make_visible(node_ind, ro.front_ind() - 1, ro_tmp, scoring);
    double maxsc = ro_tmp.order_score;
    for (size_t i = ro_tmp.front_ind(); i < ro.order.size() - 1; i++)
    {
        swap_nodes(i, i + 1, ro_tmp, scoring);

        if (definitelyGreaterThan(ro_tmp.order_score, maxsc, EPSILON))
        {
            maxsc = ro_tmp.order_score;
        }
    }
    return (maxsc);
}

bool equal_and_unordered_top(const RightOrder &ro, int new_front, OrderScoring &scoring)
{
    size_t p = ro.order.size();

    if (ro.best_insert_pos[new_front] == ro.front_ind()) // So we can only compare when its best at right after the top.. sort of..
    {
        // If new_node is best inserted as the second top ([..i..],a,b,c) before, and best pos for i is s([...],a,i,b,c)
        if (approximatelyEqual(ro.new_top_scores[new_front], ro.inserted_max_order_scores[new_front], EPSILON))
        {
            // Compare a > i in ([...], i,a,b,c) and ([...], a,i,b,c)
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

bool has_hidden_gap(const RightOrder &ro,
                    const vector<double> &top_scores,
                    OrderScoring &scoring)
{
    size_t p = ro.order.size();
    size_t n = ro.n;

    // Checking insertion of each of the hidden nodes.
    for (size_t node_index = 0; node_index < p - n; node_index++)
    {
        int inserted_node = ro.order[node_index];
        double score_at_very_top = ro.order_score + top_scores[inserted_node];
        if (definitelyGreaterThan(ro.inserted_max_order_scores[inserted_node], score_at_very_top, EPSILON))
        {
            // Here s(([...],a,b,i,c)) = "best inserted score"
            // cout << inserted_node <<" better at knd place "<< endl;
            return (true);
        }
    }

    return (false);
}

bool unord_hidden_gap(const RightOrder &ro,
                      const vector<double> &top_scores,
                      OrderScoring &scoring)
{
    size_t p = ro.order.size();
    size_t n = ro.n;

    // Checking insertion of each of the remaining nodes
    for (size_t node_index = 0; node_index < p - n; node_index++)
    {
        int inserted_node = ro.order[node_index];
        // cout <<  "Insert "<< inserted_node<< endl;
        //  Total score when the new node is at the top
        double score_at_very_top = ro.order_score + top_scores[inserted_node];
        // Special case if the best position is the second top (i.e ro's top), since then we have to check if the
        // top score is invariant under swapping the top 2 nodes and if so, prune if they are not ordered.
        if (ro.best_insert_pos[inserted_node] == ro.front_ind())
        {
            // Here, best insert score is at the second top ([...],a,i,b,c)
            // Need to compare S(i,[...],a,b,c) == S([...],a,i,b,c)
            if (approximatelyEqual(ro.inserted_max_order_scores[inserted_node], score_at_very_top, EPSILON))
            {
                // We also have to check if the top is ordered.
                // If not, we prune, i.e. return True.
                // Evaluate S([...],i,a,b,c) == S([...],a,i,b,c) and i > a
                if (inserted_node > ro.front())
                {
                    // cout <<  " Equal unordered top "<< endl;
                    // cout << ro.inserted_max_order_scores[inserted_node]<< " = "<< score_at_very_top<< endl;
                    return (true);
                }
            }
        }
    }

    return (false);
}

RightOrder init_right_order(size_t node, OrderScoring &scoring)
{
    size_t p = scoring.numparents.size();
    vector<int> order(p, 0);
    // It all vectors with all the nodes.
    // The makes it easier to keep track of which nodes are used.
    for (size_t i = 0; i < p; i++)
    {
        order[i] = i;
    }
    myswap(node, p - 1, order);
    // score nodes
    vector<double> node_scores = scoring.score(order, p - 1, 1);
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

/**
 * Note that the ro should have the insertions scores etc from its parent when entering this function.
 */
void update_insertion_scores(RightOrder &ro, OrderScoring &scoring)
{

    size_t n = ro.n;
    size_t p = ro.order.size();

    double order_score_bkp = ro.order_score; // I added this
    size_t node = ro.front();
    double node_score = ro.node_scores[node];

    // If the last one
    if (n == p)
    {
        return;
    }

    for (size_t i = 0; i < p - n; i++)
    {
        // Note that inserted_node must be put as child of node when calculating score of node. ([...],node,inserted_node,a,b,c) ?
        size_t inserted_node = ro.order[i];
        move_element(ro.order, i, ro.front_ind());                                                      // ([..,inserted_node,..],node,a,b,c) -> ([...],node,inserted_node,a,b,c)
        ro.inserted_max_order_scores[inserted_node] += scoring.score_pos(ro.order, ro.front_ind() - 1); // O(p). This can be O(1) if inserted_node is not a plus1 parent of node.
        move_element(ro.order, ro.front_ind(), i);
    }

    //  3. Go through all other nodes and compare max inserted score to inserted score
    //     between node and a in ([...],node,a,b,c).
    size_t top_ind = p - n; // order.front_ind()
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
}

RightOrder add_node_in_front(const RightOrder &ro_prev, size_t index_of_el_to_insert, OrderScoring &scoring)
{
    RightOrder ro(ro_prev); // O(p)
    ro.n = ro_prev.n + 1;
    size_t n = ro.n;
    size_t p = ro.order.size();
    size_t node = ro_prev.order[index_of_el_to_insert];

    move_element(ro.order, index_of_el_to_insert, p - n);     // put node in the front: ([...],node,a,b,c)
    double orderscore_new = ro_prev.new_top_scores[node];     // O(1) this is for the whole ro with node added in the front/top.
    double node_score = orderscore_new - ro_prev.order_score; // O(1)

    ro.node_scores[node] = node_score;
    ro.order_score = orderscore_new;
    double order_score_bkp = ro.order_score;

    // Should also update the inserted node order scores.
    // 1. Remove node from the list, since this is now part of the order.
    ro.inserted_max_order_scores[node] = 0;
    ro.best_insert_pos[node] = 0;
    return (ro);
}

/// @brief
/// @param order
/// @param n
/// @param right_type
/// @return
vector<bool> order_to_boolvec(const vector<int> &order, size_t n, bool right_type)
{
    vector<bool> boolvec(order.size(), false);
    size_t p = order.size();

    if (right_type)
    {
        for (size_t i = p - n; i < p; i++)
        {
            boolvec[order[i]] = true;
        }
    }
    else // left order
    {
        for (size_t i = 0; i < n; i++)
        {
            boolvec[order[i]] = true;
        }
    }
    return (boolvec);
}

vector<RightOrder> prune_equal_sets(vector<RightOrder> right_orders,
                                    bool right_type)
{
    vector<vector<bool>> boolmat;
    vector<double> order_scores;

    // cout << "Creating boolmatrix" << endl;
    for (const RightOrder &ro : right_orders)
    {
        vector<bool> boolvec = order_to_boolvec(ro.order, ro.n, right_type);
        boolmat.push_back(move(boolvec));
        order_scores.push_back(ro.order_score);
    }

    // cout << "Get indices of unique maximal scoring sets" << endl;
    vector<int> pruned_inds = unique_sets(boolmat, order_scores, EPSILON);
    vector<double> pruned_scores;
    vector<RightOrder> kept_ros;

    for (const auto &ind : pruned_inds)
    {
        kept_ros.push_back(right_orders[ind]);
    }

    return (kept_ros);
}

vector<int> no_right_gaps(RightOrder &ro,
                          vector<double> &top_scores,
                          OrderScoring &scoring)

{
    // get all vectors with indendent front.
    int max_indep_node = 0;
    bool has_indep_node = false;
    size_t max_indep_ind = 0;
    size_t p = ro.order.size();
    size_t n = ro.n + 1;

    vector<int> indices_to_consider;

    for (size_t node_index = 0; node_index <= p - n; node_index++) // O(p)
    {
        int new_node = ro.order[node_index];

        if (approximatelyEqual(top_scores[new_node] + ro.order_score, ro.new_top_scores[new_node], EPSILON))
        {
            // new node is independent of hidden nodes.
            has_indep_node = true;
            if (new_node > max_indep_node)
            {
                max_indep_node = new_node;
                max_indep_ind = node_index;
            }
        }
    }

    // Add the nodes with higher number that the maximal that is indep of hidden.
    if (has_indep_node)
    {
        for (size_t node_index = 0; node_index <= p - n; node_index++)
        {

            int new_node = ro.order[node_index];
            if (new_node >= max_indep_node)
            {
                indices_to_consider.push_back(node_index);
            }
        }
    }
    else // No node is indep of hidden, so add all.
    {
        for (size_t node_index = 0; node_index <= p - n; node_index++)
        {
            indices_to_consider.push_back(node_index);
        }
    }

    // cout << "has indep of hodden node " << has_indep_node << endl;

    return (indices_to_consider);
}

tuple<vector<int>, double, size_t, size_t> opruner_right(OrderScoring &scoring)
{
    size_t p = scoring.numparents.size();
    // cout << "Starting optimization" << endl;
    vector<RightOrder> right_orders;
    vector<RightOrder> right_orders_prev;

    /**
     * Compute S((x,[...])) for all x.
     */
    vector<double> top_scores(p);
    vector<int> order_tmp(p);

    for (size_t i = 0; i < p; i++)
    {
        order_tmp[i] = i;
    }
    // top_scores has scores for individual nodes.
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
    for (size_t n = 1; n <= p; n++)
    {

        if (n == 1)
        {
            // vector<RightOrder> potential_orders;
            //  cout << "Adding nodes to empty order" << endl;
            for (size_t node_index = 0; node_index <= p - n; node_index++)
            {

                RightOrder ro = init_right_order(node_index, scoring);
                if (!optimal_front(ro, scoring))
                {
                    continue;
                }
                ++orders1;

                if (has_hidden_gap(ro, top_scores, scoring))
                {
                    continue;
                }

                if (unord_hidden_gap(ro, top_scores, scoring))
                {
                    // assert(has_gap(ro, top_scores, scoring));
                    continue;
                }

                right_orders.push_back(move(ro));
                // potential_orders.push_back(move(ro));
            }

            // potential_orders = prune_indep_front(potential_orders, top_scores, scoring);
            // right_orders = potential_orders;
            orders2 = right_orders.size();
            orders3 = orders2;
        }
        else
        {
            // Loop over all the sub orders from the previous round. O( (p CR n-1) )
            for (RightOrder &prev_order : right_orders_prev)
            {
                // If some new node is independent of the hidden, we prune all with lower numbers using no_right_gaps
                // than the maximal one of these.
                vector<int> inds_to_consider = no_right_gaps(prev_order, top_scores, scoring); // O(p)

                // PrintVector(inds_to_consider);
                for (auto node_index : inds_to_consider)
                {
                    int new_node = prev_order.order[node_index];

                    // Check if the new node is optimal at front, is S(([hidden], n_node, visible)) S(([hidden], visible and n)).
                    // If not, don't consider this node for addition.
                    // This is O(1)
                    if (!optimal_front(prev_order, new_node, scoring)) // O(1)
                    {
                        continue;
                    }

                    // the top can be swapped, i.e. S(([hidden],n,m,rest-m)) = S(([hidden],m,n,rest-m))
                    // This is O(1)
                    // OBS! Here we already knoiw that new_node is optimal in the front or equly good somewhere else.
                    // so just have to check if the best insertion position is the second top and that the score at that position is equal.
                    // If so we check the ordering of the top two nodes.
                    if (equal_and_unordered_top(prev_order, new_node, scoring)) // O(1)
                    {
                        continue;
                    }

                    // If the new node survived the checks above, we add it.
                    // Update the scores for inserting other nodes after equal set pruning,
                    // to avoid O(p).
                    RightOrder ro = add_node_in_front(prev_order, node_index, scoring); // O(p)

                    right_orders.push_back(move(ro));
                }
            }

            // O(|right_orders_prev| * p) particles. space
            orders1 = right_orders.size();
            // cout << "# orders after add in front and prune indep front: " << orders1 << endl;

            // For orders with the same nodes, keep only one having the maximal score
            // O(#particles * p) = O( (p CR n) * n * p) space and time.
            right_orders = prune_equal_sets(right_orders, true);

            orders2 = right_orders.size();

            // updated the insertion scores for all the right orders
            for (RightOrder &ro : right_orders)
            {
                update_insertion_scores(ro, scoring);
            }
            // cout << "# orders after prune equal sets: " << orders2 << endl;
            //  Remove any order that has gaps.

            vector<RightOrder> right_orders_tmp;
            // After equal set pruning there are O((p CR n)) particles left so
            // this is O((p CR n))

            for (RightOrder &ro : right_orders)
            {
                // cout << "Check for gaps in: " << ro << endl;
                //  This is O(p-n) = O(p)
                if (has_hidden_gap(ro, top_scores, scoring))
                {
                    // assert(has_gap(ro, top_scores, scoring));
                    continue;
                }

                if (unord_hidden_gap(ro, top_scores, scoring))
                {
                    continue;
                }

                right_orders_tmp.push_back(move(ro));
            }
            right_orders = move(right_orders_tmp);

            orders3 = right_orders.size();
            // cout << "# orders after has gap prune: " << orders3 << endl;
        }

        // cout << "orders of size " << n << endl;

        // Print some statistics //
        /* Add to number of particles sum */
        tot_n_particles += orders3;
        /* Check that scores are correct */
        if (orders3 > max_n_particles)
        {
            max_n_particles = orders3;
        }

        // This is O((p CR n))

        // cout << "# of orders: " << right_orders.size() << endl;
        auto max_ro = max_element(right_orders.begin(),
                                  right_orders.end(),
                                  [](const RightOrder &a, const RightOrder &b)
                                  { return a.order_score < b.order_score; });

        // cout << "max scoring sub order " << endl;
        // cout << *max_ro << endl;
        // cout << "score: " << max_ro->order_score << endl;

        // check correct score
        vector<double> sc = scoring.score(max_ro->order, p - n, n); // Take only the last n elements in the vector
        // PrintVector(sc);
        double max_score_check = accumulate(sc.begin(), sc.end(), 0.0);
        // cout << "correct max score: " << max_score_check << endl;
        assert(abs(max_ro->order_score - max_score_check) < EPSILON);

        // cout << n << " & " << orders1 << " & " << orders2 << " & " << orders3 << " & " << max_ro->order_score << " \\\\" << endl;
        right_orders_prev = move(right_orders);
        // cout << "after move" << endl;
    }

    auto max_ro = max_element(right_orders_prev.begin(),
                              right_orders_prev.end(),
                              [](const RightOrder &a, const RightOrder &b)
                              { return a.order_score < b.order_score; });

    // cout << "MAX order " << *max_ro << endl;

    return (make_tuple(max_ro->order, max_ro->order_score, max_n_particles, tot_n_particles));
}

// [[Rcpp::plugins(cpp17)]]

// [[Rcpp::export]]

Rcpp::List r_opruner_right(Rcpp::List ret)
{
    OrderScoring scoring = get_score(ret);

    const auto &[order, log_score, max_n_particles, tot_n_particles] = opruner_right(scoring);

    Rcpp::List L = Rcpp::List::create(Rcpp::Named("order") = order,
                                      Rcpp::Named("log_score") = log_score,
                                      Rcpp::Named("max_n_particles") = max_n_particles,
                                      Rcpp::Named("tot_n_particles") = tot_n_particles);

    return (L);
}