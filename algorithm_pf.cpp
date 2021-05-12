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
//============================================================================

//#include "algorithm_pf.hpp"
#include <RInside.h>
#include "smc_stuff.cpp"
#include <chrono>
using namespace std::chrono;

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

    std::string r_code = "source(\"felixtestar.R\"); ret";

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
    std::vector<Rcpp::NumericMatrix> bannedscore;

    for (std::size_t i = 0; i < p; i++)
    {
        Rcpp::NumericMatrix m = Rcpp::as<Rcpp::NumericMatrix>(bannedscoreR[i]);
        bannedscore.push_back(m);
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
    std::vector<Rcpp::IntegerVector> plus1listsparents;
    for (std::size_t i = 0; i < p; i++)
    {
        Rcpp::IntegerVector m = Rcpp::as<Rcpp::IntegerVector>(plus1listsparentsR[i]);
        plus1listsparents.push_back(m);
    }

    std::vector<int> scorepositions(p);
    for (int i = 0; i < p; ++i)
    {
        scorepositions[i] = i;
    }

    // int nsmall = 4;

    // std::vector<int> sub_order(order.begin() + nsmall, order.end());
    // std::vector<int> sub_scorepositions(scorepositions.begin() + nsmall, scorepositions.end());

    // std::vector<double>
    //     orderscores = orderscorePlus1(p,
    //                                   sub_order,
    //                                   sub_scorepositions,
    //                                   aliases,
    //                                   numparents,
    //                                   rowmaps_backwards,
    //                                   plus1listsparents,
    //                                   scoretable,
    //                                   bannedscore,
    //                                   order);

    // PrintVector(order);
    // PrintVector(orderscores);

    // std::cout << "New scoring " << std::endl;
    std::map<cache_keytype, std::vector<double>> cache;
    OrderScoring scoring(aliases,
                         numparents,
                         rowmaps_backwards,
                         plus1listsparents,
                         scoretable,
                         bannedscore, cache);

    // std::vector<double> orderscores2 = scoring.score(order, 4, 4);
    // PrintVector(orderscores2);

    // int cur_order_length = 4;
    // int tot_order_length = order.size();
    // int to = order.size() - 1;

    int seed = 1;
    std::srand(seed);

    // int node_index = rand_int_in_range(cur_order_length, to);
    // std::cout << "node index " << node_index << ", node " << order[node_index] << std::endl;

    // double order_score = std::accumulate(orderscores2.begin(), orderscores2.end(), 0);
    // std::vector<double> neig_scoring = score_sub_order_neigh(scoring, order, orderscores2,
    //                                                          order_score, cur_order_length, node_index);
    // PrintVector(neig_scoring);

    // std::vector<double> neig_dist = dist_from_logprobs(neig_scoring);
    // PrintVector(neig_dist);

    std::default_random_engine generator(seed);
    // std::discrete_distribution<int> distribution(neig_dist.begin(), neig_dist.end());

    // int insert_pos = distribution(generator);

    // std::cout << "insert pos " << insert_pos << std::endl;

    // PrintVector(order);
    // move_element(order, node_index, insert_pos);

    // PrintVector(order);
    auto start = high_resolution_clock::now();
    smc(scoring, N, order.size());
    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<milliseconds>(stop - start);

    // To get the value of duration use the count()
    // member function on the duration object
    std::cout << duration.count() << " ms."<< std::endl;
    //algorithm_pf(N);

    exit(0);
}
