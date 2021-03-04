//============================================================================
// vSMC/include/vsmc/resample/residual_systematic.hpp
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

#ifndef VSMC_RESAMPLE_RESIDUAL_SYSTEMATIC_HPP
#define VSMC_RESAMPLE_RESIDUAL_SYSTEMATIC_HPP

#include <vsmc/resample/common.hpp>

namespace vsmc {

namespace internal {

typedef cxx11::integral_constant<ResampleScheme, ResidualSystematic>
    ResampleResidualSystematic;

} // namespace vsmc::internal

/// \brief Residual systematic resampling
/// \ingroup Resample
template <>
class Resample<internal::ResampleResidualSystematic>
{
    public :

    template <typename IntType, typename RngType>
    void operator() (std::size_t M, std::size_t N, RngType &rng,
            const double *weight, IntType *replication)
    {
        using std::modf;

        residual_.resize(M);
        integral_.resize(M);
        double *const rptr = &residual_[0];
        double *const iptr = &integral_[0];
        for (std::size_t i = 0; i != M; ++i)
            rptr[i] = modf(N * weight[i], iptr + i);
        double rsum = 0;
        for (std::size_t i = 0; i != M; ++i)
            rsum += rptr[i];
        for (std::size_t i = 0; i != M; ++i)
            rptr[i] /= rsum;

        IntType R = 0;
        for (std::size_t i = 0; i != M; ++i)
            R += static_cast<IntType>(iptr[i]);
        std::size_t NN = N - static_cast<std::size_t>(R);
        u01_.resize(NN);
        double *const uptr = &u01_[0];
        cxx11::uniform_real_distribution<double> runif(0, 1);
        const double delta = 1.0 / N;
        const double u = runif(rng) * delta;
        for (std::size_t i = 0; i != NN; ++i)
            uptr[i] = u + i * delta;
        inversion_(M, NN, rptr, uptr, replication, true);

        for (std::size_t i = 0; i != N; ++i)
            replication[i] += static_cast<IntType>(iptr[i]);
    }

    private :

    internal::Inversion inversion_;
    std::vector<double, AlignedAllocator<double> > residual_;
    std::vector<double, AlignedAllocator<double> > integral_;
    std::vector<double, AlignedAllocator<double> > u01_;
}; // Residual systematic resampling

} // namespace vsmc

#endif // VSMC_RESAMPLE_RESIDUAL_SYSTEMATIC_HPP
