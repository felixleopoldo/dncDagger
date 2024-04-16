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

    // read the R matrix ret$H_max from into a c++ matrix    
    Rcpp::NumericMatrix H_max_r = ret["H_max_adj"];
    Rcpp::NumericMatrix H_min_r = ret["H_min_adj"];

    // convert these into vector<vector<double>>
    vector<vector<bool>> H_max_adj(p, vector<bool>(p));
    vector<vector<bool>> H_min_adj(p, vector<bool>(p));
    // read the values from th Rcpp matrices
    for (size_t i = 0; i < p; i++)
    {
        for (size_t j = 0; j < p; j++)
        {
            if (i == j)
            {
                H_max_adj[i][i] = 0;
                continue;
            }
            H_max_adj[i][j] = H_max_r(i, j) == true;
            H_min_adj[i][j] = H_min_r(i, j) == true;

        }
    }
    // // print the matrices
    // cout << "H_max: " << endl;
    // for (size_t i = 0; i < p; i++)
    // {
    //     for (size_t j = 0; j < p; j++)
    //     {
    //         cout << H_max_adj[i][j] << " ";
    //     }
    //     cout << endl;
    // }

    // cout << "H_min: " << endl;
    // for (size_t i = 0; i < p; i++)
    // {
    //     for (size_t j = 0; j < p; j++)
    //     {
    //         cout << H_min_adj[i][j] << " ";
    //     }
    //     cout << endl;
    // }


    vector<vector<bool>> optmat = dnc(scoring, H_min_adj, H_max_adj);

    vector<RightOrder> initial_right_orders = {}; 

    auto start = high_resolution_clock::now();
    const auto &[order, log_score, node_scores, max_n_particles, tot_n_particles] = opruner_right(scoring, initial_right_orders);
    cout << log_score << endl;

    // print the order
    cout << "Order: ";
    for (size_t i = 0; i < order.size(); i++)
    {
        cout << order[i] << " ";
    }
    cout << endl;

    vector<vector<int>> dag = order_to_dag(order, scoring);
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

    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<milliseconds>(stop - start);
}

/*
    c++ version of the following code.

get_diff_matrices <- function(rowmaps, scoretable, aliases, var_labels){

    nvars <- length(rowmaps)

    H_min = matrix(, nrow = nvars, ncol = nvars)
    H_max = matrix(, nrow = nvars, ncol = nvars)
    colnames(H_min) <- var_labels
    colnames(H_max) <- var_labels
    rownames(H_min) <- var_labels
    rownames(H_max) <- var_labels

    for (i in seq(nvars)) {
        var <- rowmaps[[i]]
        # print("#############")
        # print("var:")
        # print(var_labels[[i]])
        # print("aliases:")
        # print(aliases[[i]])
        # print("forward hashs:")
        # print(var$forward)

        # total number of possible parents         
        n_pos_parents <- length(aliases[[i]])# sqrt(length(var$forward))
        
        # For the possible parents, i.e. not plus1 parents
        # We exclude each parent in turn and see how the score changes
        for (parent_ind in seq(n_pos_parents)){
            # if no possible parents, skip
            if(n_pos_parents == 0) next # since seq is weird
            parent <- aliases[[i]][[parent_ind]] # parent to exclude
            for (hash in var$forward){                
                
                # Check if the parent is in the hash
                check <- (hash-1) %% 2^(parent_ind)
                if (check < 2^(parent_ind-1)){
                    
                    # Compute the hash with the parent excluded
                    hash_with_parent <- hash + 2^(parent_ind-1)                                        
                    # Using the no plus1 score table, i.e. index 1
                    score_diff <- scoretable[[i]][[1]][var$backwards[[hash_with_parent]]] - scoretable[[i]][[1]][var$backwards[[hash]]]                    
                    # Update the H matrices
                    if (is.na(H_max[i, parent])){                        
                        H_max[i,parent] <- score_diff
                    } else {
                        H_max[i, parent] <- max(H_max[i, parent], score_diff)
                    }

                    if (is.na(H_min[i, parent])){                        
                        H_min[i,parent] <- score_diff
                    } else {
                        H_min[i, parent] <- min(H_min[i, parent], score_diff)
                    }
                 }
            }
        }

        ## For the plus1 parents        
        # plus1parents are those parents that are not in the aliases
        plus1parents <- c()
        plus1parent_inds <- c()

        j <- 1
        for (label in var_labels){
            if (j == i) {  # skip itself
                j <- j + 1
                next
            }
            if (!(label %in% labels(aliases[[i]]))){
                # label not in aliases so it is a plus1 parent
                plus1parents <- c(plus1parents, label)
                plus1parent_inds <- c(plus1parent_inds, j)
            }
            j <- j + 1
        }

        for(j in seq(1, length(plus1parents))){ 
            score_diffs <- scoretable[[i]][[j+1]] - scoretable[[i]][[1]] # the first one is the no plus1 score table. Subtracting all at once.
            H_max[i, plus1parent_inds[j]] <- max(score_diffs) 
            H_min[i, plus1parent_inds[j]] <- min(score_diffs)
        }
    }
    # print("H_max:")
    # print((H_max > 0) * 1)
    # print("H_min:")
    # print((H_min > 0) * 1)
    return(list(H_min = H_min, H_max = H_max))
}

*/

// pair<vector<vector<double>>, vector<vector<double>>> get_diff_matrices(vector<RowMap> rowmaps, vector<vector<vector<double>>> scoretable, vector<vector<int>> aliases, vector<string> var_labels)
// {
//     size_t nvars = rowmaps.size();
//     vector<vector<double>> H_min(nvars, vector<double>(nvars, 0));
//     vector<vector<double>> H_max(nvars, vector<double>(nvars, 0));

//     for (size_t i = 0; i < nvars; i++)
//     {
//         RowMap var = rowmaps[i];
//         size_t n_pos_parents = aliases[i].size();

//         for (size_t parent_ind = 0; parent_ind < n_pos_parents; parent_ind++)
//         {
//             if (n_pos_parents == 0)
//                 continue;
//             int parent = aliases[i][parent_ind];
//             for (size_t hash : var.forward)
//             {
//                 int check = (hash - 1) % (1 << parent_ind);
//                 if (check < (1 << (parent_ind - 1)))
//                 {
//                     int hash_with_parent = hash + (1 << (parent_ind - 1));
//                     double score_diff = scoretable[i][0][var.backwards[hash_with_parent]] - scoretable[i][0][var.backwards[hash]];
//                     if (isnan(H_max[i][parent]))
//                     {
//                         H_max[i][parent] = score_diff;
//                     }
//                     else
//                     {
//                         H_max[i][parent] = max(H_max[i][parent], score_diff);
//                     }

//                     if (isnan(H_min[i][parent]))
//                     {
//                         H_min[i][parent] = score_diff;
//                     }
//                     else
//                     {
//                         H_min[i][parent] = min(H_min[i][parent], score_diff);
//                     }
//                 }
//             }
//         }

//         vector<int> plus1parents;
//         vector<int> plus1parent_inds;

//         for (size_t j = 0; j < var_labels.size(); j++)
//         {
//             if (j == i)
//             {
//                 continue;
//             }
//             if (find(aliases[i].begin(), aliases[i].end(), j) == aliases[i].end())
//             {
//                 // add label to plus1parents
//                 //plus1parents.push_back(j);
//                 plus1parents.push_back(var_labels[j]);
//                 plus1parent_inds.push_back(j);

//             }
//         }

//         for (size_t j = 0; j < plus1parents.size(); j++)
//         {
//             vector<double> score_diffs;
//             for (size_t k = 0; k < scoretable[i].size(); k++)
//             {
//                 score_diffs.push_back(scoretable[i][j + 1][k] - scoretable[i][0][k]);
//             }
//             H_max[i][plus1parent_inds[j]] = *max_element(score_diffs.begin(), score_diffs.end());
//             H_min[i][plus1parent_inds[j]] = *min_element(score_diffs.begin(), score_diffs.end());
//         }
//     }
//     return {H_min, H_max};

// }







