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
#include "smc_stuff.cpp"

int main(int argc, char **argv)
{
    RInside R(argc, argv);
    // R["txt"] = "Hello, world!\n"; // assign a char* (string) to 'txt'

    // R.parseEvalQ("print(txt)");

    --argc;
    ++argv;

    std::size_t N = 3;
    if (argc > 0)
    {
        N = static_cast<std::size_t>(std::atoi(*argv));
        --argc;
        ++argv;
    }

    //std::string r_code = "source(\"felixtestar.R\"); ret";
    std::string r_code = "ret <- readRDS('jackdata.csv.rds'); ret";
    //std::string r_code = "ret <- readRDS('jackdata.csv.rds'); print(ret$bannedscore); print('aliases'); print(ret$aliases); print('rowmaps_backwards'); print(ret$rowmaps_backwards); ret";
    //std::string r_code = "ret <- readRDS('myasiandata.csv.rds'); print(ret$bannedscore); print('aliases'); print(ret$aliases); print('rowmaps_backwards'); print(ret$rowmaps_backwards); print('potential plus1 parents'); print(ret$plus1listsparents); ret";

    //std::string r_code = "source(\"readtables.R\"); ret";

    Rcpp::List ret = R.parseEval(r_code);

    // Read order
    std::vector<int> order = Rcpp::as<std::vector<int>>(ret["order"]);

    // Read numparents
    std::vector<int> numparents = Rcpp::as<std::vector<int>>(ret["numparents"]);

    // Read scoretable
    std::vector<std::vector<std::vector<double>>> scoretable = Rcpp::as<std::vector<std::vector<std::vector<double>>>>(ret["scoretable"]);

    // Read parent table
    Rcpp::List parenttableR = Rcpp::as<Rcpp::List>(ret["parenttable"]);
    std::size_t p = parenttableR.size();
    std::vector<Rcpp::IntegerMatrix> parenttable;
    for (std::size_t i = 0; i < p; i++)
    {
        Rcpp::IntegerMatrix m = Rcpp::as<Rcpp::IntegerMatrix>(parenttableR[i]);
        parenttable.push_back(m);
    }

    // Read banned score
    Rcpp::List bannedscoreR = Rcpp::as<Rcpp::List>(ret["bannedscore"]);
    //std::vector<Rcpp::NumericMatrix> bannedscore;
    std::vector<std::vector<std::vector<double>>> bannedscore;
    for (std::size_t i = 0; i < p; i++)
    {
        Rcpp::NumericMatrix m = Rcpp::as<Rcpp::NumericMatrix>(bannedscoreR[i]);
        //std::vector<std::vector<double>> mat;
        std::vector<std::vector<double>> mat(m.rows(), std::vector<double>(m.cols()));
        for (std::size_t j = 0; j < m.rows(); j++){
            for (std::size_t k = 0; k < m.cols(); k++){
                mat[j][k] = m(j,k);
                //NumericVector v = m( _ , j );
                //std::vector<double> row;
                //for(auto val: v){
                //    row.push_back(val);
                //}
            }
        }
        bannedscore.push_back(mat);
    }

    // Read rowmaps_backwards
    Rcpp::List rowmaps_backwardsR = Rcpp::as<Rcpp::List>(ret["rowmaps_backwards"]);
    std::vector<Rcpp::IntegerVector> rowmaps_backwards;
    for (std::size_t i = 0; i < p; i++)
    {
        Rcpp::IntegerVector m = Rcpp::as<Rcpp::IntegerVector>(rowmaps_backwardsR[i]);
        rowmaps_backwards.push_back(m);
    }

    // Read aliases
    Rcpp::List aliasesR = Rcpp::as<Rcpp::List>(ret["aliases"]);
    std::vector<std::vector<int>> aliases;
    for (std::size_t i = 0; i < p; i++)
    {
        std::vector<int> m = Rcpp::as<std::vector<int>>(aliasesR[i]);
        aliases.push_back(m);
    }

    // Read plus1listsparents
    Rcpp::List plus1listsparentsR = Rcpp::as<Rcpp::List>(ret["plus1listsparents"]);
    std::vector<std::vector<int>> plus1listsparents;
    for (std::size_t i = 0; i < p; i++)
    {
        Rcpp::IntegerVector m = Rcpp::as<Rcpp::IntegerVector>(plus1listsparentsR[i]);
        std::vector<int> tmp;
        for (auto e : m){
            tmp.push_back(e);
        }
        plus1listsparents.push_back(tmp);
    }

    std::vector<int> scorepositions(p);
    for (int i = 0; i < p; ++i)
    {
        scorepositions[i] = i;
    }

    std::map<cache_keytype3, std::vector<double>> cache;
    OrderScoring scoring(aliases,
                         numparents,
                         rowmaps_backwards,
                         plus1listsparents,
                         scoretable,
                         bannedscore, cache);

    int seed = 1;
    std::srand(seed);
   
    std::random_device rd;
    std::mt19937 gen(rd());
    std::default_random_engine generator(seed);




    // std::vector<int> v = {0,1,2,3,4,5,6};
    // const auto&[orders, new_nodes] = backward_order_sampler(v);

    // for(auto o : orders){
    //     PrintVector(o);
    // }



    auto start = high_resolution_clock::now();
    const auto &[smc_log_w, smc_orders, smc_log_scores] = smc(scoring, N, order.size(), generator);
    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<milliseconds>(stop - start);

    std::vector<double> * norm_w = dist_from_logprobs(smc_log_w);
    std::discrete_distribution<int> distribution(norm_w->begin(), norm_w->end());

    int maxElementIndex = std::max_element(smc_log_scores.begin(), smc_log_scores.end()) - smc_log_scores.begin();

    //PrintVector(log_order_scores[p-1]);
    std::map<std::vector<int>, double> orders_probs;
    std::set<std::vector<int>> distinct_orders;

    for (int i = 0; i < N; i++)
    {
        distinct_orders.insert(smc_orders[i]);
    }

    for (int i = 0; i < N; i++)
    {
        if (orders_probs.count(smc_orders[i]))
        {
            orders_probs[smc_orders[i]] += (*norm_w)[i];
        }
        else
        {
            orders_probs[smc_orders[i]] = (*norm_w)[i];
        }
    }

    for (auto o : distinct_orders)
    {
        PrintVector(o);
        std::vector<double> *scr = scoring.score(o, 0, p);
        double sc = std::accumulate(scr->begin(), scr->end(), 0.0);
        //std::cout << orders_probs[o] << " score " << sc << std::endl;
        delete scr;
    }

    std::cout << "number of distinct orders: " << distinct_orders.size() << std::endl;
    std::cout << "index: " << maxElementIndex << std::endl;
    //PrintVector(new_orders[maxElementIndex]);
    PrintVector(smc_orders[maxElementIndex]);

    //std::cout << "prob: " << orders_probs[new_orders[maxElementIndex]] << std::endl;
    std::cout << "prob: " << orders_probs[smc_orders[maxElementIndex]] << std::endl;
    //std::cout << "score: " << log_order_scores[maxElementIndex] << std::endl;
    std::cout << "score: " << smc_log_scores[maxElementIndex] << std::endl;
    //PrintVector(smc_log_scores);
    // To get the value of duration use the count()
    // member function on the duration object
    std::cout << duration.count() << " ms." << std::endl;
    delete norm_w;

    int M = 1000;

    const auto &[pgibbs_orders, pgibbs_log_scores] = pgibbs(M, N, scoring, generator);

    std::cout << "PGibbs log scores " << std::endl;
    PrintVector(pgibbs_log_scores);

    for (auto & o: pgibbs_orders){
        PrintVector(o);
    }
    int pgibbs_max_score_ind = std::max_element(pgibbs_log_scores.begin(), pgibbs_log_scores.end()) - pgibbs_log_scores.begin();
    std::cout << "PGibbs max log score " << pgibbs_log_scores[pgibbs_max_score_ind] << std::endl;
    exit(0);
}
