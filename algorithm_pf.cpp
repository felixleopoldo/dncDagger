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

std::vector<double> orderscorePlus1(std::size_t n,
                                    std::vector<int> &scorenodes,
                                    std::vector<int> &scorepositions,
                                    std::vector<std::vector<int>> aliases,
                                    std::vector<int> numparents,
                                    std::vector<Rcpp::IntegerVector> rowmaps_backwards,
                                    std::vector<Rcpp::IntegerVector> plus1listsparents,
                                    std::vector<std::vector<std::vector<double>>> scoretable,
                                    std::vector<Rcpp::NumericMatrix> scoresmatrices,
                                    std::vector<int> &permy)
{

    std::vector<double> orderscores(n); // orderscores <- vector("double", n)
    std::vector<std::vector<int>> allowedscorelists(n); // allowedscorelists < -vector("list", n)
    std::vector<int> therows(n);                        // therows <- vector("integer", n) # what is this? / Felix
    size_t k = 0;                                       // k <- 1

    for (int i : scorenodes)
    {
        size_t position = scorepositions[k];
        if (position == n - 1)
        {
            // no parents allowed, i.e.only first row, only first list
            orderscores[i] = scoretable[i][0][0];    // orderscores[i] <- scoretable [[i]][[1]][1, 1]
            allowedscorelists[i].assign({0});        // allowedscorelists [[i]] <- c(1)
            therows[i] = std::pow(2, numparents[i]); // therows[i] <- c(2 ^ numparents[i])
        }
        else
        {
            std::vector<int> bannednodes(permy.begin(), permy.begin() + position + 1); // bannednodes <- permy[1:position]
            std::vector<int> allowednodes(permy.begin() + position + 1, permy.end());  // allowednodes < -permy [(position + 1):n]
            std::vector<int> bannedpool;                                               // bannedpool <- which(aliases [[i]] % in % bannednodes)
            for (int j = 0; j < aliases[i].size(); j++)
            {
                if (std::find(bannednodes.begin(), bannednodes.end(), aliases[i][j]) != bannednodes.end())
                {
                    bannedpool.push_back(j);
                }
            }

            if (numparents[i] == 0 || bannedpool.size() == 0)
            {
                therows[i] = 0; // all parents allowed
            }
            else
            {
                int indextmp = 0;
                for (auto &item : bannedpool)
                {
                    indextmp += std::pow(2, item + 1); // add 1 since nodes are labeled from 0
                }
                indextmp = indextmp / 2;
                therows[i] = rowmaps_backwards[i][indextmp]; // this might be list I guess. rowmaps_backwards[i][std::sum(2 ^ bannedpool) / 2 + 1]
            }

            allowedscorelists[i].push_back(0); // allowedscorelists [[i]] <- c(1, which(plus1lists$parents [[i]] % in % allowednodes) + 1)
            for (int j = 0; j < plus1listsparents[i].size(); j++)
            {
                if (std::find(allowednodes.begin(), allowednodes.end(), plus1listsparents[i][j]) != allowednodes.end())
                {
                    allowedscorelists[i].push_back(j + 1);
                }
            }

            std::vector<double> scoresvec;
            for (int allowedscore : allowedscorelists[i])
            {
                scoresvec.push_back(scoresmatrices[i](therows[i], allowedscore));
            }

            double maxallowed = *std::max_element(scoresvec.begin(), scoresvec.end()); // maxallowed <- max(scoresvec)

            for (auto &val : scoresvec)
            {
                orderscores[i] += std::exp(val - maxallowed);
            }
            orderscores[i] = maxallowed + std::log(orderscores[i]);
        }
        k++;
    }
    //scores < -list()
    //scores$therow < -therows
    //scores$allowedlists < -allowedscorelists
    //scores$totscores < -orderscores
    return (orderscores);
}

template <typename T>
void PrintVector(std::vector<T> &arr)
{
    copy(arr.begin(), arr.end(),
         std::ostream_iterator<T>(std::cout, "; "));
    std::cout << std::endl;
}

/* Return 1 if arr2[] is a subset of arr1[] */
bool isSubset(int arr1[], int arr2[], int m, int n)
{
    int i = 0, j = 0;

    if (m < n)
        return 0;

    // Sort both the arrays
    std::sort(arr1, arr1 + m);
    std::sort(arr2, arr2 + n);

    // Iterate till they donot exceed their sizes
    while (i < n && j < m)
    {
        // If the element is smaller than
        // Move aheaad in the first array
        if (arr1[j] < arr2[i])
            j++;
        // If both are equal, then move
        // both of them forward
        else if (arr1[j] == arr2[i])
        {
            j++;
            i++;
        }

        // If we donot have a element smaller
        // or equal to the second array then break
        else if (arr1[j] > arr2[i])
            return 0;
    }

    return (i < n) ? false : true;
}

int main(int argc, char **argv)
{
    RInside R(argc, argv);
    R["txt"] = "Hello, world!\n"; // assign a char* (string) to 'txt'

    R.parseEvalQ("print(txt)");

    --argc;
    ++argv;

    std::size_t N = 1000;
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

    std::vector<double> orderscores = orderscorePlus1(p,
                                                      order,
                                                      scorepositions,
                                                      aliases,
                                                      numparents,
                                                      rowmaps_backwards,
                                                      plus1listsparents,
                                                      scoretable,
                                                      bannedscore,
                                                      order);

    PrintVector(order);
    PrintVector(orderscores);

    //smc(10)
    //algorithm_pf(N);

    exit(0);
}

// void smc(std::size_t N){

//     for (size_t j = 0; j < p; j++)
//         {
//             if (j>0)
//             {
//                 norm_w =
//             }

//             std::vector<int> I = // resample weights
//             for(size_t i=0; i<N; i++)
//             {
//                 if (j==0)
//                 {
//                     // This is the initialisation
//                     int node = j; //draw random node
//                     l = std::list<int>;
//                     l.push_front(node)
//                     indperms[i, j] = l;
//                     log_w = 0.0;
//                 }
//                 else{
//                     indperms[i, j] = sdt:list<int>(indperms[I[i], j]);
//                     log_w = 0.0;
//                 }
//             }
//         }

// }
