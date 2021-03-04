//============================================================================
// vSMC/include/vsmc/rng/rng.hpp
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

#ifndef VSMC_RNG_RNG_HPP
#define VSMC_RNG_RNG_HPP

#include <vsmc/internal/config.hpp>

#include <vsmc/rng/rng_set.hpp>
#include <vsmc/rng/seed.hpp>

#include <vsmc/rng/generator_wrapper.hpp>

#include <vsmc/rng/xor_combine_engine.hpp>
#include <vsmc/rng/xorshift.hpp>

#include <vsmc/rng/philox.hpp>
#include <vsmc/rng/threefry.hpp>

#if VSMC_HAS_AES_NI
#include <vsmc/rng/aes.hpp>
#include <vsmc/rng/aes_ni.hpp>
#include <vsmc/rng/ars.hpp>
#endif

#if VSMC_HAS_GSL
#include <vsmc/rng/gsl.hpp>
#endif

#if VSMC_HAS_MKL
#include <vsmc/rng/mkl.hpp>
#endif

#if VSMC_HAS_RDRAND
#include <vsmc/rng/rdrand.hpp>
#endif

#include <vsmc/rng/stable_distribution.hpp>
#include <vsmc/rng/u01.hpp>
#include <vsmc/rng/uniform_real_distribution.hpp>

#if VSMC_HAS_SSE2
#include <vsmc/rng/m128i.hpp>
#endif

#endif // VSMC_RNG_RNG_HPP
