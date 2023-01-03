//============================================================================
// MCKL/example/algorithm/include/algorithm_pf.hpp
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

#ifndef MCKL_EXAMPLE_ALGORITHM_PF_HPP
#define MCKL_EXAMPLE_ALGORITHM_PF_HPP

#include <mckl/algorithm/smc.hpp>
#include <mckl/random/normal_distribution.hpp>
#include <mckl/smp.hpp>
#include <mckl/utility/stop_watch.hpp>
#include <list>

// Felix:
// This says that the state is a StateMatrix.
// and that the state has 4 double values stored with the ColMajor layout.

// Since I dont know yet how to use lists (orders) of varying size,
// I use the StateMatrix with ints of maximal size. 
// The iter can be uses to see which elements to look at.

// This is an order 
//using AlgorithmPFBase = mckl::StateMatrix<mckl::ColMajor, int, 10>; 

using AlgorithmPFBase = mckl::StateMatrix<mckl::ColMajor, int, 4>;

class AlgorithmPF : public AlgorithmPFBase
{
  public:
    AlgorithmPF(std::size_t N) : AlgorithmPFBase(N)
    {
        double x = 0;
        double y = 0;
        std::ifstream data("algorithm_pf.data");
        while (data >> x >> y) {
            x_.push_back(x);
            y_.push_back(y);
        }
        data.close();
    }

    std::size_t n() const { return 10; }
    double x(std::size_t iter) { return x_[iter]; }
    double y(std::size_t iter) { return y_[iter]; }

  private:
    mckl::Vector<double> x_;
    mckl::Vector<double> y_;
}; // class AlgorithmPF

template <typename Backend>
class AlgorithmPFSelection : public mckl::SMCSamplerEvalSMP<AlgorithmPF, AlgorithmPFSelection<Backend>, Backend>
{
  public:
    void operator()(std::size_t iter, mckl::Particle<AlgorithmPF> &particle)
    {
        // run is run(std::size_t iter, Particle<T> &particle, std::size_t, Args &&...)
        //std::cout << "RUNNING " << std::endl;
        this->run(iter, particle, 1000); // why 1000?
    }

    void eval_each(std::size_t iter, mckl::ParticleIndex<AlgorithmPF> idx)
    {
        std::cout << "AlgorithmPFSelection: eval_eachccccccc " << std::endl;
        // eval_dispatch(iter, idx, &Derived::eval_each);
    }

    void eval_first(std::size_t, mckl::Particle<AlgorithmPF> &particle)
    {
        w_.resize(particle.size());
    }

    void eval_last(std::size_t, mckl::Particle<AlgorithmPF> &particle)
    {
        particle.weight().add_log(w_.data()); // Here
    }

    // eval each?

    // How does MCKL know what is the weight? - Right above in eval last, it seems.. /Felix
    void eval_range(
        std::size_t iter, const mckl::ParticleRange<AlgorithmPF> &range)
    {
        // what uis a range? - An iterator (collection) of particles. Or an interator for the particle object.
        // What is a state? - Mixt be the x,y AlgorithmPF class.. object.. mckl::StateMatrix<mckl::ColMajor, double, 4>;
        // What is a col_data?
        // What is a particle?
        auto rng = range.begin().rng();

        const std::size_t N = range.size();

        int p = range.begin().state().dim();
        // How can we run data() on double Vector? / Felix
        double *const w = w_.data() + range.ibegin(); //index begin

        std::cout << iter << std::endl << std::endl;
        
        //Print the current state
        for (std::size_t j = 0; j != N; ++j) {
            for (std::size_t i = 0; i != p; ++i) {
                std::cout << *range.particle().at(j).state().col_data(i) << std::endl;
            }
        }
        if (iter == 0) {
            for (std::size_t j = 0; j != N; ++j) {
                for (std::size_t i = 0; i != p; ++i) {
                    *range.particle().at(j).state().col_data(i) = i+100;          
                }         
            }
        } else {
            for (std::size_t j = 0; j != N; ++j) {
                for (std::size_t i = 0; i != p; ++i) {
                    // std::cout << *range.particle().at(j).state().col_data(i) << std::endl;
                    // *range.particle().at(j).state().col_data(i) = i+100                
                }         
            }

        }
        
        range.begin().rng() = rng;
    }

  private:
    mckl::Vector<double> w_; // What are these? - Noise in SSM / Felix
}; // AlgorithmPFSelection

template <typename Backend>
class AlgorithmPFPos : public mckl::SMCEstimatorEvalSMP<AlgorithmPF, AlgorithmPFPos<Backend>, Backend>
{
  public:
    void eval_each(std::size_t, std::size_t, mckl::ParticleIndex<AlgorithmPF> idx, double *r)
    {
        std::cout << "AlgorithmPFPos eval_each: where is this called.." << std::endl;
        r[0] = idx(0); // What is this? /Felix
        r[1] = idx(1);
        std::cout << "r[0]" << r[0]<< " r[1]" << r[1] << std::endl;
    }
}; // class AlgorithmPFPos

template <typename Backend>
inline std::string algorithm_pf_name();

template <>
inline std::string algorithm_pf_name<mckl::BackendSEQ>()
{
    return "SEQ";
}

template <>
inline std::string algorithm_pf_name<mckl::BackendSTD>()
{
    return "STD";
}

template <>
inline std::string algorithm_pf_name<mckl::BackendOMP>()
{
    return "OMP";
}

template <>
inline std::string algorithm_pf_name<mckl::BackendTBB>()
{
    return "TBB";
}

template <typename Backend>
inline void algorithm_pf(std::size_t N)
{
    // AlgorithmPF is defined above and derives from 
    // AlgorithmPFBase = mckl::StateMatrix<mckl::ColMajor, double, 4>
    // SMCSampler derives from Sampler<SMCSampler<T, U>>
    // see smc.hpp? should take 2 args...
    mckl::SMCSampler<AlgorithmPF> sampler(N); // greates sampler object in the stack.
    // AlgorithmPFSelection is defined above and derives from 
    // mckl::SMCSamplerEvalSMP<AlgorithmPF, AlgorithmPFSelection<Backend>, Backend>
    // defined in backend_omp.hpp
    sampler.selection(AlgorithmPFSelection<Backend>()); //see above
    sampler.resample(mckl::Stratified);
    sampler.resample_threshold(0.5);
    sampler.selection_estimator(mckl::SMCEstimator<AlgorithmPF>(2, AlgorithmPFPos<Backend>()));

    std::cout << "resample estimator "<< std::endl;
    sampler.resample_estimator(mckl::SMCEstimator<AlgorithmPF>(2, AlgorithmPFPos<Backend>()));
    sampler.mutation_estimator(mckl::SMCEstimator<AlgorithmPF>(2, AlgorithmPFPos<Backend>()));

    std::cout << "iterate "<< std::endl;
    sampler.iterate(sampler.particle().state().n());
    sampler.clear();


    // std::cout << "iterate with watch "<< std::endl;
    // mckl::StopWatch watch;
    // watch.start();
    // sampler.iterate(sampler.particle().state().n());
    // watch.stop();

    const std::string name = algorithm_pf_name<Backend>();
    // std::cout << name << ": " << watch.seconds() << 's' << std::endl;

    std::ofstream save("algorithm_pf_" + name + ".save");
    save << sampler << std::endl;
    save.close();
}

inline void algorithm_pf(std::size_t N)
{
    std::cout << std::fixed << std::setprecision(3);

    algorithm_pf<mckl::BackendSEQ>(N);
    //algorithm_pf<mckl::BackendSTD>(N);
    //algorithm_pf<mckl::BackendOMP>(N);
#if MCKL_HAS_TBB
    algorithm_pf<mckl::BackendTBB>(N);
#endif
}

#endif // MCKL_EXAMPLE_ALGORITHM_PF_HPP
