#include <RInside.h>
#include <Rcpp.h>
#include <chrono>
#include "include/opruner_right.h"
#include "include/OrderScoring.h"
#include "include/dnc.h"
//#include "include/path_pruning.h"

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
    std::string plus1it = "2";
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

    auto start = high_resolution_clock::now();
 
    // // convert these into vector<vector<double>>
    vector<vector<bool>> H_max_adj(p, vector<bool>(p));
    vector<vector<bool>> H_min_adj(p, vector<bool>(p));
  
    // print H_min_adj and H_max_adj
    tie(H_min_adj, H_max_adj) = get_diff_matrices(scoring);

    vector<int> dnc_order = dnc(scoring, H_min_adj, H_max_adj);
    vector<vector<bool>> dnc_mat = order_to_dag(dnc_order, scoring);

    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<milliseconds>(stop - start);
    cout << "Time taken by dnc function: " << duration.count() << " milliseconds" << endl;
    
    print_matrix(dnc_mat);

    vector<RightOrder> initial_right_orders = {}; 

     start = high_resolution_clock::now();
    const auto &[order, log_score, node_scores, max_n_particles, tot_n_particles] = opruner_right(scoring, initial_right_orders);
    cout << log_score << endl;

    // print the order
    cout << "Order: ";
    for (size_t i = 0; i < order.size(); i++)
    {
        cout << order[i] << " ";
    }
    cout << endl;

    vector<vector<bool>> dag = order_to_dag(order, scoring);
    // print as matrix
    cout << "DAG: " << endl;
    for (size_t i = 0; i < dag.size(); i++)
    {
        for (size_t j = 0; j < dag[i].size(); j++)
        {
            cout << dag[i][j] << " ";
        }
        cout << endl;
    }    

     stop = high_resolution_clock::now();
     duration = duration_cast<milliseconds>(stop - start);
}
