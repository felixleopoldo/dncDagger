#include <RInside.h>
#include <Rcpp.h>
#include <chrono>
#include "include/opruner_right.h"
#include "include/OrderScoring.h"
#include "include/dnc.h"
#include <fstream>
// #include "include/path_pruning.h"

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

// function that writes the matrix to a csv file using std::ofstream
// Function to write a boolean matrix to a CSV file
void writeBoolMatrixToCSV(const vector<vector<bool>>& matrix, const string& filename, const vector<string>& labels) {
    std::ofstream file(filename);

    if (!file.is_open()) {
        cerr << "Could not open the file!" << endl;
        return;
    }
    // write the labels
    for (size_t i = 0; i < labels.size(); ++i) {
        file << labels[i];
        if (i != labels.size() - 1) {
            file << ",";
        }
    }
    file << "\n";

    for (const auto& row : matrix) {
        for (size_t col = 0; col < row.size(); ++col) {
            file << (row[col] ? "1" : "0");
            if (col != row.size() - 1) {
                file << ",";
            }
        }
        file << "\n";
    }

    file.close();
}


int main(int argc, char **argv)
{
    std::string datafilename;
    std::string scoretype;
    std::string am = "1";
    std::string chi = "0.5";
    std::string aw = "NULL";
    std::string edgepf = "2";
    std::string output_csv = "output.csv";
    bool leftorder = false;

    if (argc < 4)
    {
        std::cout << "usage: ./run_opruner --filename datafilename --scoretype [bde|bge] [--am am|--chi chi] [--aw aw|--edgepf edgepf]" << std::endl;
        std::cout << "example: ./run_opruner --filename data/asia.csv --scoretype bge " << std::endl;
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
             
            if (strcmp(argv[i], "--output_csv") == 0)
            {
                output_csv = argv[i + 1];
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

    auto start = high_resolution_clock::now();

    const auto &[dnc_order, dnc_score, dnc_max_particles, dnc_tot_n_particles] = dnc(scoring);
    vector<vector<bool>> dnc_mat = order_to_dag(dnc_order, scoring);
    
    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<milliseconds>(stop - start);
    cout << "Time taken by dnc function: " << duration.count() << " milliseconds" << endl;

    print_matrix(dnc_mat);
    // write to csv
    // get labels from the datafilename csv
    std::ifstream file(datafilename);
    std::string line;
    std::getline(file, line);
    std::vector<std::string> labels = split(line, ",");
    


    writeBoolMatrixToCSV(dnc_mat, output_csv, labels);


    // cout << "Dnc score: " << dnc_score << endl;
    // cout << "Dnc max particles: " << dnc_max_particles << endl;
    // cout << "Dnc total particles: " << dnc_tot_n_particles << endl;

    // vector<RightOrder> initial_right_orders = {};

    // start = high_resolution_clock::now();
    // const auto &[order, log_score, node_scores, max_n_particles, tot_n_particles] = opruner_right(scoring, initial_right_orders);
    // cout << log_score << endl;

    // // print the order
    // cout << "Order: ";
    // for (size_t i = 0; i < order.size(); i++)
    // {
    //     cout << order[i] << " ";
    // }
    // cout << endl;
    // cout << "Score: " << log_score << endl;
    // cout << "Max particles: " << max_n_particles << endl;
    // cout << "Total particles: " << tot_n_particles << endl;

    // vector<vector<bool>> dag = order_to_dag(order, scoring);
    // // print as matrix    
    // stop = high_resolution_clock::now();
    // cout << "DAG: " << endl;
    // print_matrix(dag);
    // duration = duration_cast<milliseconds>(stop - start);
    // cout << "Time taken by opruner_right function: " << duration.count() << " milliseconds" << endl;
}
