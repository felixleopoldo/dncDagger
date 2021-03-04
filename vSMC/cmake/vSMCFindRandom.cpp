//============================================================================
// vSMC/cmake/vSMCFindRandom.cpp
//----------------------------------------------------------------------------
//                         vSMC: Scalable Monte Carlo
//----------------------------------------------------------------------------
// Copyright (c) 2013,2014, Yan Zhou
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

#include <cassert>
#include <vsmc/cxx11/random.hpp>
#include <vsmc/rng/rng_set.hpp>

#if VSMC_HAS_RANDOM123
#include <Random123/philox.h>
#include <Random123/conventional/Engine.hpp>
#endif

int main ()
{
    const int N = 1000;

    vsmc::cxx11::uniform_real_distribution<> ru01(0, 1);
    double u;
    {vsmc::cxx11::mt19937       eng; u = ru01(eng); assert(u >=0 && u <= 1);}
    {vsmc::cxx11::mt19937_64    eng; u = ru01(eng); assert(u >=0 && u <= 1);}
    {vsmc::cxx11::minstd_rand0  eng; u = ru01(eng); assert(u >=0 && u <= 1);}
    {vsmc::cxx11::minstd_rand   eng; u = ru01(eng); assert(u >=0 && u <= 1);}
    {vsmc::cxx11::ranlux24_base eng; u = ru01(eng); assert(u >=0 && u <= 1);}
    {vsmc::cxx11::ranlux48_base eng; u = ru01(eng); assert(u >=0 && u <= 1);}
    {vsmc::cxx11::ranlux24      eng; u = ru01(eng); assert(u >=0 && u <= 1);}
    {vsmc::cxx11::ranlux48      eng; u = ru01(eng); assert(u >=0 && u <= 1);}
    {vsmc::cxx11::knuth_b       eng; u = ru01(eng); assert(u >=0 && u <= 1);}

#if VSMC_HAS_RANDOM123
    vsmc::RngSet<r123::Engine<r123::Philox2x32>, vsmc::Vector> eng(N);
#else
    vsmc::traits::RngSetTypeTrait<vsmc::NullType>::type eng(N);
#endif

    vsmc::cxx11::uniform_real_distribution<> runif(0, 1);
    for (int i = 0; i != N; ++i) {
        double u = runif(eng[i]);
        assert(u >= 0 && u <= 1);
    }

    vsmc::cxx11::uniform_int_distribution<std::size_t> ruint(1, 100);
    for (int i = 0; i != N; ++i) {
        std::size_t u = ruint(eng[i]);
        assert(u >= 1 && u <= 100);
    }

    vsmc::cxx11::bernoulli_distribution rbern(0.5);
    for (int i = 0; i != N; ++i){
        int b = rbern(eng[i]);
        assert(b == 0 || b == 1);
    }

    vsmc::cxx11::binomial_distribution<> rbinom(20, 0.7);
    for (int i = 0; i != N; ++i) {
        int b = rbinom(eng[i]);
        assert(b >= 0 && b <= 20);
    }

    vsmc::cxx11::gamma_distribution<> rgamma(1, 1);
    for (int i = 0; i != N; ++i) {
        double g = rgamma(eng[i]);
        assert(g >= 0);
    }

    vsmc::cxx11::lognormal_distribution<> rlnom(1, 1);
    for (int i = 0; i != N; ++i) {
        double l = rlnom(eng[i]);
        assert(l >= 0);
    }

    vsmc::cxx11::normal_distribution<> rnorm(0, 1);
    for (int i = 0; i != N; ++i) {
        double n = rnorm(eng[i]);
    }

    return 0;
}
