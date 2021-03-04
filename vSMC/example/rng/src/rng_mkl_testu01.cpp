//============================================================================
// vSMC/example/rng/src/rng_mkl_testu01.cpp
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

#include "rng_testu01.hpp"
#include <vsmc/rng/mkl.hpp>

VSMC_RNG_TESTU01_FUNCTION(MKL_MCG59)
VSMC_RNG_TESTU01_FUNCTION(MKL_MT19937)
VSMC_RNG_TESTU01_FUNCTION(MKL_MT19937_64)
VSMC_RNG_TESTU01_FUNCTION(MKL_MT2203)
VSMC_RNG_TESTU01_FUNCTION(MKL_MT2203_64)
VSMC_RNG_TESTU01_FUNCTION(MKL_SFMT19937)
VSMC_RNG_TESTU01_FUNCTION(MKL_SFMT19937_64)
#if VSMC_HAS_RDRAND
VSMC_RNG_TESTU01_FUNCTION(MKL_NONDETERM)
VSMC_RNG_TESTU01_FUNCTION(MKL_NONDETERM_64)
#endif

int main (int argc, char **argv)
{
    VSMC_RNG_TESTU01_OPTION_PRE;

    VSMC_RNG_TESTU01_OPTION(MKL_MCG59);
    VSMC_RNG_TESTU01_OPTION(MKL_MT19937);
    VSMC_RNG_TESTU01_OPTION(MKL_MT19937_64);
    VSMC_RNG_TESTU01_OPTION(MKL_MT2203);
    VSMC_RNG_TESTU01_OPTION(MKL_MT2203_64);
    VSMC_RNG_TESTU01_OPTION(MKL_SFMT19937);
    VSMC_RNG_TESTU01_OPTION(MKL_SFMT19937_64);
#if VSMC_HAS_RDRAND
    VSMC_RNG_TESTU01_OPTION(MKL_NONDETERM);
    VSMC_RNG_TESTU01_OPTION(MKL_NONDETERM_64);
#endif

    VSMC_RNG_TESTU01_OPTION_POST;

    VSMC_RNG_TESTU01(MKL_MCG59);
    VSMC_RNG_TESTU01(MKL_MT19937);
    VSMC_RNG_TESTU01(MKL_MT19937_64);
    VSMC_RNG_TESTU01(MKL_MT2203);
    VSMC_RNG_TESTU01(MKL_MT2203_64);
    VSMC_RNG_TESTU01(MKL_SFMT19937);
    VSMC_RNG_TESTU01(MKL_SFMT19937_64);
#if VSMC_HAS_RDRAND
    VSMC_RNG_TESTU01(MKL_NONDETERM);
    VSMC_RNG_TESTU01(MKL_NONDETERM_64);
#endif

    return 0;
}
