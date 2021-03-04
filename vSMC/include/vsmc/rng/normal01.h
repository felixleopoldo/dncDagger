//============================================================================
// vSMC/include/vsmc/rng/normal01.h
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

#ifndef VSMC_RNG_NORMAL01_H
#define VSMC_RNG_NORMAL01_H

/// \page normal01 Normal distribution
///
/// - Header: `<vsmc/rng/normal01.h>`
/// - Distribution: Standard Normal
/// \f$ f(x) = \frac{1}{\sqrt{2\pi}} e^{-x^2/2} \f$
/// - Method: Box-Muller
///
/// ### Types and functions
///
/// ~~~{.c}
/// normal01_<N>x<W>_<F>;
/// void normal01_<N>x<W>_<F>_init (normal01_<N>x<W>_<F> *, cburng<N>x<W>_rng_t *, <FT> shape);
/// <FT> normal01_<N>x<W>_<F>_init (normal01_<N>x<W>_<F> *, cburng<N>x<W>_rng_t *);
/// ~~~
///
/// ### Macros
///
/// ~~~{.c}
/// NORMAL01_<N>x<W>
/// NORMAL01_<N>x<W>_INIT
/// NORMAL01_<N>x<W>_RAND
/// ~~~
///
/// ### Examples
///
/// ~~~{.c}
/// #define VSMC_HAS_OPENCL_DOUBLE 1
/// #include <vsmc/rng/normal01.h>
///
/// cburng4x32_rng_t rng;
/// cburng4x32_init(&rng);
///
/// double shape = 2;
/// normal01_4x32_53 rnorm_53;
/// normal01_4x32_53_init(&rnorm_53, &rng, shape);
/// double r_53 = normal01_4x32_53_rand(&rnorm_53, &rng);
///
/// NORMAL01_4x32 rnorm;
/// NORMAL01_4x32_INIT(&rnorm, &rng, shape);
/// double r = NORMAL01_4x32_RAND(&rnorm, &rng);
/// ~~~

#include <vsmc/rng/internal/common.h>
#include <vsmc/rng/urng.h>
#include <vsmc/rng/u01.h>

#if VSMC_HAS_OPENCL_DOUBLE

#define NORMAL01_2x32      normal01_2x32_53
#define NORMAL01_2x32_INIT normal01_2x32_53_init
#define NORMAL01_2x32_RAND normal01_2x32_53_rand

#define NORMAL01_2x64      normal01_2x64_53
#define NORMAL01_2x64_INIT normal01_2x64_53_init
#define NORMAL01_2x64_RAND normal01_2x64_53_rand

#define NORMAL01_4x32      normal01_4x32_53
#define NORMAL01_4x32_INIT normal01_4x32_53_init
#define NORMAL01_4x32_RAND normal01_4x32_53_rand

#define NORMAL01_4x64      normal01_4x64_53
#define NORMAL01_4x64_INIT normal01_4x64_53_init
#define NORMAL01_4x64_RAND normal01_4x64_53_rand

#else // VSMC_HAS_OPENCL_DOUBLE

#define NORMAL01_2x32      normal01_2x32_24
#define NORMAL01_2x32_INIT normal01_2x32_24_init
#define NORMAL01_2x32_RAND normal01_2x32_24_rand

#define NORMAL01_4x32      normal01_4x32_24
#define NORMAL01_4x32_INIT normal01_4x32_24_init
#define NORMAL01_4x32_RAND normal01_4x32_24_rand

#endif // VSMC_HAS_OPENCL_DOUBLE

#define VSMC_DEFINE_RNG_NORMAL01(N, W, F, FT) \
    typedef struct {                                                         \
        FT u1;                                                               \
        FT u2;                                                               \
        unsigned char saved;                                                 \
    } normal01_##N##x##W##_##F;

#define VSMC_DEFINE_RNG_NORMAL01_INIT(N, W, F, FT) \
    VSMC_STATIC_INLINE void normal01_##N##x##W##_##F##_init (                \
            normal01_##N##x##W##_##F *rnorm, cburng##N##x##W##_rng_t *rng)   \
    {                                                                        \
        rnorm->u1 = u01_open_closed_##W##_##F(cburng##N##x##W##_rand(rng));  \
        rnorm->u2 = u01_open_closed_##W##_##F(cburng##N##x##W##_rand(rng));  \
        rnorm->saved = 1;                                                    \
    }

#define VSMC_DEFINE_RNG_NORMAL01_RAND(N, W, F, FT) \
    VSMC_STATIC_INLINE FT normal01_##N##x##W##_##F##_rand (                  \
            normal01_##N##x##W##_##F *rnorm, cburng##N##x##W##_rng_t *rng)   \
    {                                                                        \
        const FT c_2pi = 6.2831853071795865;                                 \
        if (rnorm->saved) {                                                  \
            rnorm->saved = 0;                                                \
            return sqrt(-2 * log(rnorm->u1)) * cos(c_2pi * rnorm->u2);       \
        } else {                                                             \
            normal01_##N##x##W##_##F##_init(rnorm, rng);                     \
            return sqrt(-2 * log(rnorm->u1)) * sin(c_2pi * rnorm->u2);       \
        }                                                                    \
    }

/// \ingroup CLRNG
VSMC_DEFINE_RNG_NORMAL01(2, 32, 24, float)
/// \ingroup CLRNG
VSMC_DEFINE_RNG_NORMAL01(4, 32, 24, float)

/// \ingroup CLRNG
VSMC_DEFINE_RNG_NORMAL01_INIT(2, 32, 24, float)
/// \ingroup CLRNG
VSMC_DEFINE_RNG_NORMAL01_INIT(4, 32, 24, float)

/// \ingroup CLRNG
VSMC_DEFINE_RNG_NORMAL01_RAND(2, 32, 24, float)
/// \ingroup CLRNG
VSMC_DEFINE_RNG_NORMAL01_RAND(4, 32, 24, float)

#if VSMC_HAS_OPENCL_DOUBLE

/// \ingroup CLRNG
VSMC_DEFINE_RNG_NORMAL01(2, 32, 53, double)
/// \ingroup CLRNG
VSMC_DEFINE_RNG_NORMAL01(4, 32, 53, double)

/// \ingroup CLRNG
VSMC_DEFINE_RNG_NORMAL01_INIT(2, 32, 53, double)
/// \ingroup CLRNG
VSMC_DEFINE_RNG_NORMAL01_INIT(4, 32, 53, double)

/// \ingroup CLRNG
VSMC_DEFINE_RNG_NORMAL01_RAND(2, 32, 53, double)
/// \ingroup CLRNG
VSMC_DEFINE_RNG_NORMAL01_RAND(4, 32, 53, double)

/// \ingroup CLRNG
VSMC_DEFINE_RNG_NORMAL01(2, 64, 53, double)
/// \ingroup CLRNG
VSMC_DEFINE_RNG_NORMAL01(4, 64, 53, double)

/// \ingroup CLRNG
VSMC_DEFINE_RNG_NORMAL01_INIT(2, 64, 53, double)
/// \ingroup CLRNG
VSMC_DEFINE_RNG_NORMAL01_INIT(4, 64, 53, double)

/// \ingroup CLRNG
VSMC_DEFINE_RNG_NORMAL01_RAND(2, 64, 53, double)
/// \ingroup CLRNG
VSMC_DEFINE_RNG_NORMAL01_RAND(4, 64, 53, double)

#endif // VSMC_HAS_OPENCL_DOUBLE

#endif // VSMC_RNG_NORMAL01_H
