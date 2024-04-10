#include <RInside.h>
#include <Rcpp.h>
#include <chrono>
#include "include/opruner_right.h"
#include "include/OrderScoring.h"
#include "include/path_pruning.h"

using namespace std::chrono;

// for string delimiter
std::vector<std::string> split(std::string s, std::string delimiter)
{
    size_t pos_start = 0, pos_end, delim_len = delimiter.length();
    std::string token;
    std::vector<std::string> res;

    while ((pos_end = s.find(delimiter, pos_start)) != std::string::npos)
    {
        token = s.substr(pos_start, pos_end - pos_start);
        pos_start = pos_end + delim_len;
        res.push_back(token);
    }

    res.push_back(s.substr(pos_start));
    return res;
}

/**
 * Get DAG associated with order.
*/
vector<vector<int>> get_dag(const vector<int> &order, const OrderScoring &scoring)
{
    vector<vector<int>> dag(scoring.numparents.size(), vector<int>(scoring.numparents.size(), 0));
    for (size_t i = 0; i < order.size(); i++)
    {
        int node = order[i];
        
        // What should we do now ?:)
        scoring.maxmatrix[node];
        scoring.maxrow[node];

    }
    return dag;
}

int main(int argc, char **argv)
{
    std::string datafilename;
    std::string scoretype;
    std::string am = "1";
    std::string chi = "0.5";
    std::string aw = "NULL";
    std::string edgepf = "2";
    bool leftorder = false;

    if (argc < 4)
    {
        std::cout << "usage: ./run_opruner --filename datafilename --scoretype [bde|bge] [--am am|--chi chi] [--aw aw|--edgepf edgepf]" << std::endl;
        std::cout << "example: ./run_opruner --filename data/p20n300gaussdata.csv --scoretype bge " << std::endl;
        return (0);
    }

    datafilename = *argv;

    for (int i = 1; i < argc; i++)
    {
        if (i + 1 != argc)
        {
            if (strcmp(argv[i], "--filename") == 0)
            {
                datafilename = argv[i + 1];
                i++;
            }

            if (strcmp(argv[i], "--scoretype") == 0)
            {
                scoretype = argv[i + 1];
                i++;
            }

            if (strcmp(argv[i], "--aw") == 0)
            {
                aw = argv[i + 1];
                i++;
            }

            if (strcmp(argv[i], "--am") == 0)
            {
                am = argv[i + 1];
                i++;
            }

            if (strcmp(argv[i], "--chi") == 0)
            {
                chi = argv[i + 1];
                i++;
            }

            if (strcmp(argv[i], "--edgepf") == 0)
            {
                edgepf = argv[i + 1];
                i++;
            }
        }
    }

    std::cout << datafilename << std::endl;
    std::string r_code;
    std::string plus1it = "1";
    if (scoretype == "bge")
    {
        r_code = "source('R/scoring.R'); set.seed(1); ret <- get_scores('" + datafilename + "', scoretype='" + scoretype + "', bgepar=list(am=" + am + ", aw=" + aw + "), plus1it=" + plus1it + "); ret";
    }
    else if (scoretype == "bde")
    {
        r_code = "source('R/scoring.R'); set.seed(1); ret <- get_scores('" + datafilename + "', scoretype='" + scoretype + "', bdepar=list(chi=" + chi + ", edgepf=" + edgepf + "), plus1it=" + plus1it + "); ret";
    }

    RInside R(argc, argv);
    Rcpp::List ret = R.parseEval(r_code);

    OrderScoring scoring = get_score(ret);
    size_t p = scoring.numparents.size();
    vector<double> top_scores = get_unrestricted_vec(p, scoring);

    vector<RightOrder> initial_right_orders = {};
    
    //RightOrder initial_ro = init_right_order(top_scores, scoring);
    // initial_ro = add_node_in_front(initial_ro, 3, top_scores, scoring);
    // update_insertion_scores(initial_ro, scoring);
    // initial_ro = add_node_in_front(initial_ro, 1, top_scores, scoring);
    // update_insertion_scores(initial_ro, scoring);
    // initial_ro = add_node_in_front(initial_ro, 3, top_scores, scoring);
    // update_insertion_scores(initial_ro, scoring);
    // initial_ro = add_node_in_front(initial_ro, 1, top_scores, scoring);
    // update_insertion_scores(initial_ro, scoring);

    // For each of the hidden nodes, we want to set the best
    // insertion position to be the last position.
    // Manipulate the best insertion position to be the last position:
    // Set the best insertion position to be the last position
    // Set inserted_max_order_scores to be the same as new_top_scores
    // for (size_t i = 0; i < initial_ro.size_hidden(); i++)
    // {
    //     size_t hidden_node = initial_ro.order[i];
    //     // initial_ro.best_insert_pos[hidden_node] = initial_ro.front_ind()-1;
    //     // Make the best inserted score worse than adding it to the front
    //     initial_ro.inserted_max_order_scores[hidden_node] = initial_ro.new_top_scores[hidden_node] - 100000;
    //     cout << "hidden node: " << hidden_node << endl;
    //     cout << initial_ro.best_insert_pos[hidden_node] << endl;
    //     cout << initial_ro.new_top_scores[hidden_node] << endl;
    //     cout << initial_ro.inserted_max_order_scores[hidden_node] << endl;
    // }

    // initial_right_orders.push_back(initial_ro);

    auto start = high_resolution_clock::now();
    const auto &[order, log_score, node_scores, max_n_particles, tot_n_particles] = opruner_right(scoring, initial_right_orders);
    cout << log_score << endl;
    // Get the DAG associated with order
    vector<vector<int>> dag = get_dag(order, scoring);


    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<milliseconds>(stop - start);
}
