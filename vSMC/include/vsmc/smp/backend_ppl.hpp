//============================================================================
// vSMC/include/vsmc/smp/backend_ppl.hpp
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

#ifndef VSMC_SMP_BACKEND_PPL_HPP
#define VSMC_SMP_BACKEND_PPL_HPP

#include <vsmc/smp/backend_base.hpp>
#include <ppl.h>

namespace vsmc {

VSMC_DEFINE_SMP_FORWARD(PPL)

/// \brief Particle::value_type subtype using Parallel Pattern Library
/// \ingroup PPL
template <typename BaseState>
class StatePPL : public BaseState
{
    public :

    typedef typename traits::SizeTypeTrait<BaseState>::type size_type;

    explicit StatePPL (size_type N) : BaseState(N) {}

    template <typename IntType>
    void copy (size_type N, const IntType *copy_from)
    {
        VSMC_RUNTIME_ASSERT_SMP_BACKEND_BASE_COPY_SIZE_MISMATCH(PPL);

        ::concurrency::parallel_for(static_cast<size_type>(0), N,
                copy_work_<IntType>(this, copy_from));
    }

    private :

    template <typename IntType>
    struct copy_work_
    {
        copy_work_ (StatePPL<BaseState> *state, const IntType *copy_from) :
            state_(state), copy_from_(copy_from) {}

        void operator() (size_type to) const
        {state_->copy_particle(copy_from_[to], to);}

        private :

        StatePPL<BaseState> *const state_;
        const IntType *const copy_from_;
    }; // class copy_work_
}; // class StatePPL

/// \brief Sampler<T>::init_type subtype using Parallel Pattern Library
/// \ingroup PPL
template <typename T, typename Derived>
class InitializePPL : public InitializeBase<T, Derived>
{
    public :

    std::size_t operator() (Particle<T> &particle, void *param)
    {
        typedef typename Particle<T>::size_type size_type;
        const size_type N = static_cast<size_type>(particle.size());
        this->initialize_param(particle, param);
        this->pre_processor(particle);
        ::concurrency::combinable<std::size_t> accept(accept_init_);
        ::concurrency::parallel_for(static_cast<size_type>(0), N,
                work_(this, &particle, &accept));
        this->post_processor(particle);

        return accept.combine(accept_accu_);
    }

    protected :

    VSMC_DEFINE_SMP_IMPL_COPY(PPL, Initialize)

    private :

    struct work_
    {
        typedef typename Particle<T>::size_type size_type;

        work_ (InitializePPL<T, Derived> *init,
                Particle<T> *particle,
                ::concurrency::combinable<std::size_t> *accept) :
            init_(init), particle_(particle), accept_(accept) {}

        void operator() (size_type i) const
        {
            accept_->local() += init_->initialize_state(
                    SingleParticle<T>(i, particle_));
        }

        private :

        InitializePPL<T, Derived> *const init_;
        Particle<T> *const particle_;
        ::concurrency::combinable<std::size_t> *const accept_;
    }; // class work_

    static std::size_t accept_init_ () {return 0;}
    static std::size_t accept_accu_ (std::size_t a, std::size_t b)
    {return a + b;}
}; // class InitializePPL

/// \brief Sampler<T>::move_type subtype using Parallel Pattern Library
/// \ingroup PPL
template <typename T, typename Derived>
class MovePPL : public MoveBase<T, Derived>
{
    public :

    std::size_t operator() (std::size_t iter, Particle<T> &particle)
    {
        typedef typename Particle<T>::size_type size_type;
        const size_type N = static_cast<size_type>(particle.size());
        this->pre_processor(iter, particle);
        ::concurrency::combinable<std::size_t> accept(accept_init_);
        ::concurrency::parallel_for(static_cast<size_type>(0), N,
                work_(this, iter, &particle, &accept));
        this->post_processor(iter, particle);

        return accept.combine(accept_accu_);
    }

    protected :

    VSMC_DEFINE_SMP_IMPL_COPY(PPL, Move)

    private :

    struct work_
    {
        typedef typename Particle<T>::size_type size_type;

        work_ (MovePPL<T, Derived> *move, std::size_t iter,
                Particle<T> *particle,
                ::concurrency::combinable<std::size_t> *accept):
            move_(move), particle_(particle), accept_(accept), iter_(iter) {}

        void operator() (size_type i) const
        {
            accept_->local() += move_->move_state(iter_,
                    SingleParticle<T>(i, particle_));
        }

        private :

        MovePPL<T, Derived> *const move_;
        Particle<T> *const particle_;
        ::concurrency::combinable<std::size_t> *const accept_;
        const std::size_t iter_;
    }; // class work_

    static std::size_t accept_init_ () {return 0;}
    static std::size_t accept_accu_ (std::size_t a, std::size_t b)
    {return a + b;}
}; // class MovePPL

/// \brief Monitor<T>::eval_type subtype using Parallel Pattern Library
/// \ingroup PPL
template <typename T, typename Derived>
class MonitorEvalPPL : public MonitorEvalBase<T, Derived>
{
    public :

    void operator() (std::size_t iter, std::size_t dim,
            const Particle<T> &particle, double *res)
    {
        typedef typename Particle<T>::size_type size_type;
        const size_type N = static_cast<size_type>(particle.size());
        this->pre_processor(iter, particle);
        ::concurrency::parallel_for(static_cast<size_type>(0), N,
                work_(this, iter, dim, &particle, res));
        this->post_processor(iter, particle);
    }

    protected :

    VSMC_DEFINE_SMP_IMPL_COPY(PPL, MonitorEval)

    private :

    struct work_
    {
        typedef typename Particle<T>::size_type size_type;

        work_ (MonitorEvalPPL<T, Derived> *monitor,
                std::size_t iter, std::size_t dim,
                const Particle<T> *particle, double *res) :
            monitor_(monitor), particle_(particle), res_(res),
            iter_(iter), dim_(dim) {}

        void operator() (size_type i) const
        {
            monitor_->monitor_state(iter_, dim_,
                    ConstSingleParticle<T>(i, particle_), res_ + i * dim_);
        }

        private :

        MonitorEvalPPL<T, Derived> *const monitor_;
        const Particle<T> *const particle_;
        double *const res_;
        const std::size_t iter_;
        const std::size_t dim_;
    }; // class work_
}; // class MonitorEvalPPL

/// \brief Path<T>::eval_type subtype using Parallel Pattern Library
/// \ingroup PPL
template <typename T, typename Derived>
class PathEvalPPL : public PathEvalBase<T, Derived>
{
    public :

    double operator() (std::size_t iter, const Particle<T> &particle,
            double *res)
    {
        typedef typename Particle<T>::size_type size_type;
        const size_type N = static_cast<size_type>(particle.size());
        this->pre_processor(iter, particle);
        ::concurrency::parallel_for(static_cast<size_type>(0), N,
                work_(this, iter, &particle, res));
        this->post_processor(iter, particle);

        return this->path_grid(iter, particle);
    }

    protected :

    VSMC_DEFINE_SMP_IMPL_COPY(PPL, PathEval)

    private :

    struct work_
    {
        typedef typename Particle<T>::size_type size_type;

        work_ (PathEvalPPL<T, Derived> *path, std::size_t iter,
                const Particle<T> *particle, double *res) :
            path_(path), particle_(particle), res_(res), iter_(iter) {}

        void operator() (size_type i) const
        {
            res_[i] = path_->path_state(iter_,
                    ConstSingleParticle<T>(i, particle_));
        }

        private :

        PathEvalPPL<T, Derived> *const path_;
        const Particle<T> *const particle_;
        double *const res_;
        const std::size_t iter_;
    }; // class work_
}; // PathEvalPPL

} // namespace vsmc

#endif // VSMC_SMP_BACKEND_PPL_HPP
