#include <cassert>
#include <chrono>
#include <iostream>
#include "auxiliary.h"
#include "OrderScoring.h"
#include "RightOrder.h"
#include "LeftOrder.h"
#include "path_pruning.h"
#include <boost/graph/graphviz.hpp>
#include <queue>

LeftOrder extract_leftorder(RightOrder &ro, RightOrder &reference_order, OrderScoring &scoring)
{
  size_t p = ro.order.size();

  // set the internal order to be the same as the reference order for the hidden nodes.
  vector<int> left_order(p, -1);
  int i = 0;
  for (auto &node : reference_order.order)
  {
    // if node is hidden, add it to the left order.
    if (find(ro.hidden_begin(), ro.hidden_end(), node) != ro.hidden_end())
    {
      left_order[i] = node; // reference_order.order[node];
      i++;
    }
  }

  // score the left order.
  vector<double> left_node_scores = scoring.score(left_order, 0, p - ro.n); // O(p^2) ?
  double left_order_score = accumulate(left_node_scores.begin(), left_node_scores.end(), 0.0);

  LeftOrder left_order_obj(left_order, left_order_score, left_node_scores, p - ro.n); // O(p)
  return (left_order_obj);
}

UndirectedGraph prim(UndirectedGraph &g)
{

  using namespace boost;
  vector<graph_traits<UndirectedGraph>::vertex_descriptor> MST(num_vertices(g));          // PredecessorMap
  property_map<UndirectedGraph, edge_weight_t>::type weightmap = get(edge_weight_t(), g); // weights

  prim_minimum_spanning_tree(g, &MST[0]);

  // create a new graph with only the edges from prim
  UndirectedGraph MST_graph(num_vertices(g));

  for (size_t i = 0; i != MST.size(); ++i)
  {
    if (MST[i] != i)
    {
      UEdge ed = edge(MST[i], i, g).first;
      // Take from the original
      double weight = get(weightmap, ed);           //
      add_edge(i, MST[i], 1.0 / weight, MST_graph); // these weights are the originals.
    }
  }
  return MST_graph;
}

LeftOrder topo_left_order(DirectedGraph &g, RightOrder &ro, OrderScoring &scoring)
{
  // Create left order from topological order and score it.
  typedef vector<Vertex> container;
  container c;
  topological_sort(g, std::back_inserter(c));
  vector<int> top_order;
  // for (container::reverse_iterator ii = c.rbegin(); ii != c.rend(); ++ii) {
  for (container::iterator ii = c.begin(); ii != c.end(); ++ii)
  {
    top_order.push_back(ro.order[*ii]); // TODO: verify that this is correct.
  }
  top_order.insert(top_order.end(), ro.begin(), ro.end());                         // add the ro nodes to the right
  vector<double> topo_node_scores = scoring.score(top_order, 0, ro.size_hidden()); // O(p^2)?
  double topo_order_score = accumulate(topo_node_scores.begin(), topo_node_scores.end(), 0.0);
  LeftOrder topo_left_order = LeftOrder(top_order, topo_order_score, topo_node_scores, ro.size_hidden());

  return topo_left_order;
}

vector<vector<int>> get_toporders(UndirectedGraph &G)
{
  vector<vector<int>> toporders;
  return toporders;
}

vector<int> get_topo_order(DirectedGraph &G)
{
  vector<int> visited;
  queue<int> to_visit;
  vector<int> children;
  // first add all the source nodes to the queue

  // These firs are choildredn o god sort of
  BOOST_FOREACH (Vertex v, vertices(G))
  {
    if (in_degree(v, G) == 0)
    {
      children.push_back(v);
    }
  }
  // then order them
  sort(children.begin(), children.end());
  // add them to the queue
  for (auto &child : children)
  {
    to_visit.push(child);
  }

  // Now visit each od them on order
  while (to_visit.size() > 0)
  {
    // pop from the queue
    int v = to_visit.front();
    to_visit.pop();
    visited.push_back(v);

    // get the children of v and sort them
    children.clear();
    BOOST_FOREACH (Vertex child, adjacent_vertices(v, G))
    {
      children.push_back(child);
    }
    // sort them
    sort(children.begin(), children.end());

    // add them to the queue
    for (auto &child : children)
    {
      to_visit.push(child);
    }
  }
  reverse(visited.begin(), visited.end());
  return visited;
}

DirectedGraph edmonds(DirectedGraph &G)
{
  using namespace boost;
  // This is how we can get a property map that gives the weights of
  // the edges.
  property_map<DirectedGraph, edge_weight_t>::type weights =
      get(edge_weight_t(), G);

  // This is how we can get a property map mapping the vertices to
  // integer indices.
  property_map<DirectedGraph, vertex_index_t>::type vertex_indices =
      get(vertex_index_t(), G);

  // Print the graph (or rather the edges of the graph).
  // cout << "This is the graph:\n";
  // BOOST_FOREACH (Edge e, edges(G))
  // {
  //     cout << "(" << source(e, G) << ", "
  //          << target(e, G) << ")\t"
  //          << get(weights, e) << "\n";
  // }

  // Find the maximum branching.
  vector<Edge> edmonds_edges;
  edmonds_optimum_branching<true, false, false>(G,
                                                vertex_indices,
                                                weights,
                                                static_cast<Vertex *>(0),
                                                static_cast<Vertex *>(0),
                                                back_inserter(edmonds_edges));

  // typedef adjacency_list<vecS, vecS, directedS> CGraph;
  // typedef graph_traits<CGraph>::vertex_descriptor CVertex;
  // typedef graph_traits<CGraph>::edge_descriptor CEdge;

  DirectedGraph cgraph(num_vertices(G));
  for (auto e : edmonds_edges)
  {
    add_edge(source(e, G), target(e, G), cgraph);
  }

  return cgraph;
  //  return edmonds_edges;

  // Print the edges of the maximum branching
  // cout << "This is the maximum branching (for negative weights)\n";
  // BOOST_FOREACH (Edge e, branching)
  // {
  //     cout << "(" << source(e, G) << ", "
  //          << target(e, G) << ")\t"
  //          << get(weights, e) << "\n";
  // }

  // Compute the weight of the maximum branching.

  // tuple<double, vector<Edge>> result = make_tuple(-weight, branching);
}

vector<int> get_suborder_prim(RightOrder &ro, vector<vector<double>> &M)
{
  // Graph hard_rest_graph = get_boost_ugraph(H);

  vector<int> ret(ro.order.size() - ro.n);
  return (ret);
}

vector<RightOrder> prune_path(RightOrder &reference_order,
                              vector<RightOrder> &right_orders,
                              vector<vector<double>> &M,
                              vector<vector<double>> &H,
                              vector<double> &bottom_scores,
                              vector<double> &top_scores,
                              OrderScoring &scoring)
{

  using namespace boost;
  vector<RightOrder> kept_right_orders;

  for (auto &ro : right_orders)
  {
    // The hidden nodes
    std::vector<std::string> name(ro.size_hidden());
    for (size_t i = 0; i < ro.size_hidden(); i++)
    {
      name[i] = to_string(ro.order[i]);
    }

    cout << endl;
    cout << "Considering right order " << ro << endl;

    /************* Edmonds ********/
    // Edmond finds MAX spanning tree, so we set negative weights.
    DirectedGraph loose_rest_graph = get_boost_dgraph(M, ro);
    // negative weights so we get the MAX spanning tree
    DirectedGraph edmond_tree = edmonds(loose_rest_graph);

    // ********** Prim MST **********
    // We set inverse weights, since prime finds the minimum spanning tree
    UndirectedGraph hard_rest_graph = get_boost_ugraph(H, ro);
    UndirectedGraph prim_tree = prim(hard_rest_graph);

    property_map<UndirectedGraph, edge_weight_t>::type prim_weights = get(edge_weight_t(), prim_tree);
    cout << "prim dot tree" << endl;

    // write_graphviz(cout, prim_tree, make_label_writer(&name[0]));

    cout << "prim tree" << endl;
    BOOST_FOREACH (UEdge e, boost::edges(prim_tree))
    {
      int from = boost::source(e, prim_tree);
      int to = boost::target(e, prim_tree);
      cout << ro.order[from] << "-" << ro.order[to] << ": weight " << get(prim_weights, e) << endl;
    }

    // cout << "edmond dot tree" << endl;
    // write_graphviz(cout, edmond_tree, make_label_writer(&name[0]));
    cout << "edmond tree" << endl;
    property_map<DirectedGraph, edge_weight_t>::type edmond_weights = get(edge_weight_t(), edmond_tree); // there are no weights for this one..

    BOOST_FOREACH (Edge e, boost::edges(edmond_tree))
    {
      int from = boost::source(e, edmond_tree);
      int to = boost::target(e, edmond_tree);
      cout << ro.order[to] << "<-" << ro.order[from] << ": weight " << get(edmond_weights, e) << endl;
    }
    /************ Edmond topological ************/
    LeftOrder edmond_topo = topo_left_order(edmond_tree, ro, scoring);
    vector<int> edmond_topo2 = get_topo_order(edmond_tree);
    // print the topological order
    cout << "edmond topo order: "  << endl;
    for (auto &i : edmond_topo2)
    {
      cout << ro.order[i] << " ";
    }
    cout << endl;
    /********** Hidden nodes put in the same order as in the reference_order **********/
    LeftOrder left_order = extract_leftorder(ro, reference_order, scoring);

    // ********* Compute upper bound for the current right order (ro) *****/
    // get weight of min spanning tree
    property_map<DirectedGraph, edge_weight_t>::type weights = get(edge_weight_t(), edmond_tree);
    double min_span_tree_weight = 0;
    BOOST_FOREACH (Edge e, edges(edmond_tree))
    {
      min_span_tree_weight -= get(weights, e); // OBS!!! Negative, again, as when defining the weights.
    }
    // get the unrestricted weight sum
    double hidden_unrestr_sum = 0.0;
    for (size_t i = 0; i < ro.size_hidden(); i++)
    {
      hidden_unrestr_sum += top_scores[ro.order[i]];
    }
    double edmond_upper_bound_hidden = min_span_tree_weight + hidden_unrestr_sum;
    double ro_upper_bound = edmond_upper_bound_hidden + ro.order_score;

    /************ Compute Prim lower bound *************/
    double prim_lower_bound = 0.0;
    BOOST_FOREACH (UEdge e, boost::edges(prim_tree))
    {
      prim_lower_bound += get(prim_weights, e);
    }
    // Add the diagonal element scores
    double hard_restree_weight = 0.0;
    for (size_t i = 0; i < ro.size_hidden(); i++)
    {
      hard_restree_weight += bottom_scores[ro.order[i]];
    }
    prim_lower_bound += hard_restree_weight;
    double ro_lower_bound = prim_lower_bound + ro.order_score;

    // Get linear extensions of the Prim tree
    // vector<LeftOrder> prim_left_orders = get_all_left_orders(prim_tree, ro, scoring);

    /************ Pruning *************/
    cout << "left order bounds [prim, edmond]" << endl;
    cout << "[" << prim_lower_bound << ", " << edmond_upper_bound_hidden << "]" << endl;

    cout << "edmond_topo (hidden) " << edmond_topo << endl;

    cout << "extracted left order (of hidden nodes) " << left_order.order_score << endl;
    cout << "topo_order_score (of hidden nodes) " << edmond_topo.order_score << endl;

    cout << "full order bounds [prim, edmond]" << endl;
    cout << "[" << ro_lower_bound << ", " << ro_upper_bound << "]" << endl;

    cout << "edmond_topo (full order) " << edmond_topo.order_score + ro.order_score << endl;

    cout << "reference_order.order_score " << reference_order.order_score << endl;

    if (ro_upper_bound <= reference_order.order_score)
    {
      cout << "******************************** Edmond Pruning " << ro << endl;
      continue;
    }
    else
    {
      kept_right_orders.push_back(ro);
    }

    if (approximatelyEqual(prim_lower_bound, edmond_upper_bound_hidden, 0.0000001))
    {
      cout << "******************************** Upper bound = lower bound " << endl;
      cout << "******************************** Prim Pruning " << ro << endl;
      continue;
    }

    /************** Updating the current best (v*) *************/

    if (edmond_topo.order_score + ro.order_score > reference_order.order_score)
    {
      cout << "Replacing reference order by Edmonds topologocal order" << endl;
      reference_order = ro + edmond_topo;
    }
    if (left_order.order_score + ro.order_score > reference_order.order_score)
    {
      cout << "Replacing reference order by Edmonds topologocal order" << endl;
      reference_order = ro + left_order;
    }
  }
  // print number of kept right orders and the nomber of right orders
  cout << "right_orders.size() " << right_orders.size() << endl;
  cout << "kept_right_orders.size() " << kept_right_orders.size() << endl;
  return (kept_right_orders);
}

UndirectedGraph get_boost_ugraph(vector<vector<double>> &my_weights, RightOrder &ro)
{
  int N = ro.size_hidden();
  cout << "N " << N << endl;
  UndirectedGraph G(N);

  for (int i = 0; i < N - 1; i++)
  {
    for (int j = i + 1; j < N; j++)
    {
      // OBS!!! inverse weights as Prim finds minimal spanning tree,
      // but we want to maximum spanning tree.
      add_edge(i,
               j,
               1.0 / my_weights[ro.order[i]][ro.order[j]],
               // my_weights[ro.order[i]][ro.order[j]],
               G);
    }
  }

  return (G);
}

DirectedGraph get_boost_dgraph(vector<vector<double>> &myweights, RightOrder &ro)
{
  int N = ro.size_hidden();
  DirectedGraph G(N);

  for (int i = 0; i < N; i++)
  {
    for (int j = 0; j < N; j++)
    {
      // OBS!!! Negative weights as edmond finds maximal matching!
      add_edge(i,
               j,
               -myweights[ro.order[i]][ro.order[j]],
               G);
    }
  }

  return (G);
}

vector<vector<double>> get_unrestr_mat(size_t p, OrderScoring &scoring)
{
  // vector of integers from 0 to p-1
  vector<int> order(p);
  iota(order.begin(), order.end(), 0);
  size_t k = p;

  // create a matrix of scores for all suborders of size k.
  vector<vector<double>> M(k, vector<double>(k));
  // create symmetric matrix of scores for all suborders of size k.
  vector<double> unrestr(k, 0);
  // working with the hidden nodes here
  for (size_t i = 0; i < k; i++)
  {
    move_element(order, i, 0); // put i in the front (i,...)

    unrestr[i] = scoring.score_pos(order, 0);
    for (size_t j = 0; j < k; j++)
    {
      if (i == j)
      {
        M[i][j] = unrestr[i];
        continue;
      }
      // put j in the front so that i is at the second place (j,i,...)
      size_t j_ind = j;
      if (i > j)
      {
        j_ind = j + 1; // +1 if j>i since the indices gets shifted
      }
      move_element(order, j_ind, 0);         //
      M[i][j] = scoring.score_pos(order, 1); // score of i in (j,i,...)
      // move back j.
      move_element(order, 0, j_ind);
    }
    // move back i.
    move_element(order, 0, i);
  }

  // subtract each row by its diagonal elemets.
  for (size_t i = 0; i < k; i++)
  {
    for (size_t j = 0; j < k; j++)
    {
      M[i][j] -= unrestr[i];
      assert(M[i][j] <= 0);
    }
  }

  // for each element and its transpose, keep the maximal value.
  double EPSILON = 0.0000001;
  for (size_t i = 0; i < k; i++)
  {
    for (size_t j = i + 1; j < k; j++)
    {
      double m = max(M[i][j], M[j][i]);
      if (definitelyLessThan(M[i][j], m, EPSILON))
      { // maybe approximate here? EPSILON
        M[i][j] = 0;
      }

      if (definitelyLessThan(M[j][i], m, EPSILON))
      { // maybe approximate here? EPSILON
        M[j][i] = 0;
      }
    }
  }

  return (M); // should keep track of the order of the nodes.
}

vector<vector<double>> get_hard_restr_mat(size_t p, OrderScoring &scoring)
{
  vector<int> order(p);
  iota(order.begin(), order.end(), 0);
  size_t end_ind = p - 1;

  // create a matrix of scores
  vector<vector<double>> M(p, vector<double>(p));
  // vector of unrestriced scores
  vector<double> allrestr(p, 0);
  // working with the hidden nodes here
  for (size_t i = 0; i < p; i++)
  {
    move_element(order, i, end_ind); // put i in the back (...,i)

    allrestr[i] = scoring.score_pos(order, end_ind);
    for (size_t j = 0; j < p; j++)
    {
      if (i == j)
      {
        M[i][j] = allrestr[i]; // score i in (...,i)
        continue;
      }
      // put j in the back so that i is at the second last place (...,i,j).
      size_t j_ind = j;
      if (i < j)
      {
        j_ind = j - 1;
      }
      move_element(order, j_ind, end_ind);
      M[i][j] = scoring.score_pos(order, end_ind - 1); // score i in (...,i,j)
      // move back j.
      move_element(order, end_ind, j_ind);
    }
    // move back i.
    move_element(order, end_ind, i);
  }

  // subtract each row by its diagonal elemets.
  // The diagonal will be zero of course.
  for (size_t i = 0; i < p; i++)
  {
    for (size_t j = 0; j < p; j++)
    {
      M[i][j] -= allrestr[i];
      //            assert(M[i][j] >= 0);
    }
  }

  return (M); // should keep track of the order of the nodes.
}

vector<double> get_unrestricted_vec(size_t p, OrderScoring &scoring)
{
  /**
   * Compute S((x,[...])) for nodes 0,...,p-1.
   */
  vector<int> order_tmp(p);
  iota(order_tmp.begin(), order_tmp.end(), 0);
  vector<double> top_scores(p);
  // top_scores has scores for individual nodes.
  for (size_t i = 0; i < p; i++)
  {
    move_element(order_tmp, i, 0);
    top_scores[i] = scoring.score_pos(order_tmp, 0);
    move_element(order_tmp, 0, i);
  }
  return (top_scores);
}

/**
 * Compute S(([...],x,)) for nodes 0,...,p-1.
 */
vector<double> all_restr(size_t p, OrderScoring &scoring)
{
  vector<int> order_tmp(p);
  iota(order_tmp.begin(), order_tmp.end(), 0);
  vector<double> all_restr(p);
  // bottom scores for individual nodes.
  for (size_t i = 0; i < p; i++)
  {
    move_element(order_tmp, i, p - 1);
    all_restr[i] = scoring.score_pos(order_tmp, p - 1);
    move_element(order_tmp, p - 1, i);
  }
  return (all_restr);
}