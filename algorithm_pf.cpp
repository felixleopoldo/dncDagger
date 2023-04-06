//============================================================================
// MCKL/example/algorithm/src/algorithm_pf.cpp
//----------------------------------------------------------------------------
// MCKL: Monte Carlo Kernel Library
//----------------------------------------------------------------------------
// Copyright (c) 2013-2018, Yan Zhou
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//   Redistributions of source code must retain the above copyright notice,
//   this list of conditions and the following disclaimer.
//
//   Redistributions in binary form must reproduce the above copyright notice,
//   this list of conditions and the following disclaimer in the documentation
//   and/or other materials provided with the distribution.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS AS IS
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
// SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
// INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
// CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
// ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.
// ============================================================================

//#include "algorithm_pf.hpp"
#include <RInside.h>
#include <Rcpp.h>
#include "seq_opt.cpp"
#include "thread_pool.hpp"
#include <experimental/any>
#include <chrono>

//#include "OrderScoring.cpp"
using namespace std::chrono;

int main(int argc, char **argv)
{
    //--argc;        
    //++argv; // For the program name

    //std::size_t N = 3;
    std::string datafilename;
    std::string scoretype;
    //double am = 1;
    std::string am = "1";
    //double chi = 0.5;
    std::string chi = "0.5";
    //double aw = NULL;
    std::string aw = "NULL";
    //double edgepf = 2;
    std::string edgepf = "2";

    if(argc < 4){
        std::cout << "usage: ./orderpruner --filename datafilename --scoretype [bde|bge] [--am am|--chi chi] [--aw aw|--edgepf aw]" << std::endl;
        return(0);
    }

    datafilename = *argv;
    // ++argv;

    // scoretype = *argv;
    // ++argv;
    

    // for (size_t i = 0; i < argc+1; i++){

    //     std::cout << "hej" << argv[i];
    // }

    for (int i = 1; i < argc; i++)
    {  
        if (i + 1 != argc)
        {
            if (strcmp(argv[i], "--filename") == 0) // This is your parameter name
            {                 
                datafilename = argv[i + 1];    // The next value in the array is your value
                i++;    // Move to the next flag
                std::cout << datafilename << std::endl;
            }

            if (strcmp(argv[i], "--scoretype") == 0) // This is your parameter name
            {                 
                scoretype = argv[i + 1];    // The next value in the array is your value
                i++;    // Move to the next flag
                std::cout << scoretype << std::endl;
            }

            if (strcmp(argv[i], "--aw") == 0) // This is your parameter name
            {                 
                //char* aw = argv[i + 1];    // The next value in the array is your value
                //aw = static_cast<double>(std::atof(argv[i+1]));
                aw = argv[i+1];
                i++;    // Move to the next flag
                std::cout << aw << std::endl;
            }

            if (strcmp(argv[i], "--am") == 0) // This is your parameter name
            {                 
                //char* aw = argv[i + 1];    // The next value in the array is your value
                //am = static_cast<double>(std::atof(argv[i+1]));
                am = argv[i+1];
                i++;    // Move to the next flag
                std::cout << am << std::endl;
            }

            if (strcmp(argv[i], "--chi") == 0) // This is your parameter name
            {                 
                //char* aw = argv[i + 1];    // The next value in the array is your value
                //chi = static_cast<double>(std::atof(argv[i+1]));
                chi = argv[i+1];
                i++;    // Move to the next flag
                std::cout << chi << std::endl;
            }

            if (strcmp(argv[i], "--edgepf") == 0) // This is your parameter name
            {                 
                //char* aw = argv[i + 1];    // The next value in the array is your value
                //edgepf = static_cast<double>(std::atof(argv[i+1]));
                edgepf = argv[i+1];
                i++;    // Move to the next flag
                std::cout << edgepf << std::endl;
            }
        }
    }

    std::cout << datafilename << std::endl;
    std::string r_code;
    if(scoretype == "bge") {
        r_code = "source('helper_functions.R'); ret <- get_scores('" + datafilename + "', scoretype='"+scoretype+"', bgepar=list(am="+am+", aw="+aw+")); ret";
    } else if(scoretype == "bde")
    {
        r_code = "source('helper_functions.R'); ret <- get_scores('" + datafilename + "', scoretype='"+scoretype+"', bdepar=list(chi="+chi+", edgepf="+edgepf+")); ret";
    } 
    RInside R(argc, argv);
    
    Rcpp::List ret = R.parseEval(r_code);
    OrderScoring scoring = get_score(ret);

    for (auto &mat : scoring.scoresmatrices)
    {
        for (size_t i = 0; i < mat.size(); i++)
        {
            for (size_t j = 0; j < mat[i].size(); j++)
            {
                //mat[i][j] = 0.0;
            }
        }
    }

    for (auto &mat : scoring.scoretable)
    {
        for (size_t i = 0; i < mat.size(); i++)
        {
            for (size_t j = 0; j < mat[i].size(); j++)
            {
                //mat[i][j] = 0.0;
            }
        }
    }

    // int seed = 1;
    // std::srand(seed);
    // std::random_device rd;
    // std::mt19937 gen(rd());
    // std::default_random_engine generator(seed);

    // std::vector<int> myorder = {1,3,6,7,5,2,0,4}; // myasiandata.csv
    ////std::vector<int> myorder = {7,6,2,1,5,0,4,3}; // gobnilp asia
    //std::vector<int> myorder = {1,3,6,7,5,2,0,4}; // gobnilp asia
    // std::vector<double> sc = scoring.score(myorder,0,8);
    // double score_check = std::accumulate(sc.begin(), sc.end(), 0.0);
    // PrintVector(myorder);
    // PrintVector(sc, myorder);
    // PrintVector(sc, {0,1,2,3,4,5,6,7});
    // //PrintVector(sc);
    // std::cout << score_check << std::endl;

    auto start = high_resolution_clock::now();

    const auto &[order, log_score, max_n_particles, tot_n_particles] = sequential_opt(scoring);
    std::cout << log_score << std::endl;
    // std::cout << definitelyGreaterThan(0.0003, 0.0002, 0.001) << std::endl;
    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<milliseconds>(stop - start);
    //std::cout << "Execution time: " << duration.count() << std::endl;
    // int M = 100000000;
    // const auto &[max_order, mh_log_scores] = mh(M, scoring, generator);
    // int mh_max_score_ind = std::max_element(mh_log_scores.begin(), mh_log_scores.end()) - mh_log_scores.begin();
    // std::cout << "MH only swap move " << M << " iterations" << std::endl;
    // std::cout << "max order" << std::endl;
    // PrintVector(max_order);
    // std::cout << "score: " << mh_log_scores[mh_max_score_ind] << std::endl;

    // thread_pool pool;
    // //pool.push_task(task, arg1, arg2);
    // //pool.wait_for_tasks();

    // const auto &[pgibbs_orders, pgibbs_log_scores] = pgibbs(M, N, scoring, generator, pool);

    // std::cout << "PGibbs log scores " << std::endl;
    // PrintVector(pgibbs_log_scores);

    // for (auto &o : pgibbs_orders)
    // {
    //     PrintVector(o);
    // }
    // int pgibbs_max_score_ind = std::max_element(pgibbs_log_scores.begin(), pgibbs_log_scores.end()) - pgibbs_log_scores.begin();
    // std::cout << "PGibbs max log score " << pgibbs_log_scores[pgibbs_max_score_ind] << std::endl;
}

// auto start = high_resolution_clock::now();
// const auto &[smc_log_w, smc_orders, smc_log_scores] = smc(scoring, N, order.size(), generator);
// auto stop = high_resolution_clock::now();
// auto duration = duration_cast<milliseconds>(stop - start);

// std::vector<double> * norm_w = dist_from_logprobs(smc_log_w);
// std::discrete_distribution<int> distribution(norm_w->begin(), norm_w->end());

// int maxElementIndex = std::max_element(smc_log_scores.begin(), smc_log_scores.end()) - smc_log_scores.begin();

// //PrintVector(log_order_scores[p-1]);
// std::map<std::vector<int>, double> orders_probs;
// std::set<std::vector<int>> distinct_orders;

// for (int i = 0; i < N; i++)
// {
//     distinct_orders.insert(smc_orders[i]);
// }

// for (int i = 0; i < N; i++)
// {
//     if (orders_probs.count(smc_orders[i]))
//     {
//         orders_probs[smc_orders[i]] += (*norm_w)[i];
//     }
//     else
//     {
//         orders_probs[smc_orders[i]] = (*norm_w)[i];
//     }
// }

// for (auto o : distinct_orders)
// {
//     PrintVector(o);
//     std::vector<double> *scr = scoring.score(o, 0, p);
//     double sc = std::accumulate(scr->begin(), scr->end(), 0.0);
//     //std::cout << orders_probs[o] << " score " << sc << std::endl;
//     delete scr;
// }

// std::cout << "number of distinct orders: " << distinct_orders.size() << std::endl;
// std::cout << "index: " << maxElementIndex << std::endl;
// //PrintVector(new_orders[maxElementIndex]);
// PrintVector(smc_orders[maxElementIndex]);

// //std::cout << "prob: " << orders_probs[new_orders[maxElementIndex]] << std::endl;
// std::cout << "prob: " << orders_probs[smc_orders[maxElementIndex]] << std::endl;
// //std::cout << "score: " << log_order_scores[maxElementIndex] << std::endl;
// std::cout << "score: " << smc_log_scores[maxElementIndex] << std::endl;
// //PrintVector(smc_log_scores);
// // To get the value of duration use the count()
// // member function on the duration object
// std::cout << duration.count() << " ms." << std::endl;
// delete norm_w;
