#include <RInside.h>
#include <Rcpp.h>
#include <chrono>
#include "include/opruner_right.h"
#include "include/OrderScoring.h"

using namespace std::chrono;

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

    RightOrder ro = init_right_order(0, scoring);



    auto start = high_resolution_clock::now();
    const auto &[order, log_score, max_n_particles, tot_n_particles] = opruner_right(scoring);
    cout << log_score << endl;
    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<milliseconds>(stop - start);

}
