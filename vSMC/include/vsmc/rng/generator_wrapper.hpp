//============================================================================
// vSMC/include/vsmc/rng/generator_wrapper.hpp
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

#ifndef VSMC_RNG_GENERATOR_WRAPPER_HPP
#define VSMC_RNG_GENERATOR_WRAPPER_HPP

#include <vsmc/rng/internal/common.hpp>

#define VSMC_STATIC_ASSERT_RNG_GENERATOR_WRAPPER_RESULT_TYPE(ResultType) \
    VSMC_STATIC_ASSERT((                                                     \
                cxx11::is_same<ResultType, uint16_t>::value ||               \
                cxx11::is_same<ResultType, uint32_t>::value ||               \
                cxx11::is_same<ResultType, uint64_t>::value),                \
            USE_GeneratorWrapper_WITH_RESULT_TYPE_OTHER_THAN_uint16_t_OR_uint32_t_OR_uint64_t)

#define VSMC_STATIC_ASSERT_RNG_GENERATOR_WRAPPER \
    VSMC_STATIC_ASSERT_RNG_GENERATOR_WRAPPER_RESULT_TYPE(ResultType);

namespace vsmc {

namespace traits {

/// \brief Default traits of GeneratorWrapper
/// \ingroup Traits
template <typename ResultType, typename = void>
struct GeneratorWrapperMinMaxTrait
{
    static VSMC_CONSTEXPR const ResultType _Min = 0;
    static VSMC_CONSTEXPR const ResultType _Max = static_cast<ResultType>(
            ~(static_cast<ResultType>(0)));

    static VSMC_CONSTEXPR ResultType min VSMC_MNE () {return _Min;}
    static VSMC_CONSTEXPR ResultType max VSMC_MNE () {return _Max;}
}; // struct GeneratorWrapperMinMaxTrait

} // namespace traits

/// \brief A thin wrapper over any RNG Generator for use with C++11 API
/// \ingroup RNGAdapter
///
/// \details
/// Generator need to have a member function named `generate` which returns
/// uniform integers on the range of `ResultType`, which has to be one of
/// `uint16_t`, `uint32_t` and `uint64_t`. Note that, the generator does not
/// necessarily returns one of these types. For example, on some platforms,
/// some RNG returns `unsigned long long`, which is 64-bits, while `uint64_t`
/// is also 64-bits but a typedef of `unsigned long`. It is ok as long as the
/// `static_cast` from the returned value to `uint64_t` does not change, add,
/// or lose any bits.
///
/// For most RNG generators, use only need to write a thin wrapper to use this
/// class. Also note that the `operator==` will always return false as this
/// wrapper assume that the internal states of the RNG cannot be determined by
/// this library (otherwise will have to require more interfaces of Generator).
/// And all member functions, except `operator()`, does nothing. To seed and
/// change the engine, use the Generator type object.
template <typename ResultType, class Generator, typename Traits =
    traits::GeneratorWrapperMinMaxTrait<ResultType, Generator> >
class GeneratorWrapper : public Traits
{
    public :

    typedef ResultType result_type;

    explicit GeneratorWrapper (result_type = 0)
    {VSMC_STATIC_ASSERT_RNG_GENERATOR_WRAPPER;}

    template <typename SeedSeq>
    explicit GeneratorWrapper (SeedSeq &,
            typename cxx11::enable_if<internal::is_seed_seq<SeedSeq,
            result_type, GeneratorWrapper<ResultType, Generator, Traits>
            >::value>::type * = VSMC_NULLPTR)
    {VSMC_STATIC_ASSERT_RNG_GENERATOR_WRAPPER;}

    void seed (result_type) {}

    template <typename SeedSeq>
    void seed (SeedSeq &,
            typename cxx11::enable_if<internal::is_seed_seq<SeedSeq,
            result_type, GeneratorWrapper<ResultType, Generator, Traits>
            >::value>::type * = VSMC_NULLPTR) {}

    result_type operator() ()
    {return static_cast<result_type>(generator_.generate());}

    void discard (std::size_t nskip)
    {
        for (std::size_t i = 0; i != nskip; ++i)
            operator()();
    }

    Generator &generator () {return generator_;}

    const Generator &generator () const {return generator_;}

    friend inline bool operator== (
            const GeneratorWrapper<ResultType, Generator, Traits> &,
            const GeneratorWrapper<ResultType, Generator, Traits> &)
    {return false;}

    friend inline bool operator!= (
            const GeneratorWrapper<ResultType, Generator, Traits> &,
            const GeneratorWrapper<ResultType, Generator, Traits> &)
    {return true;}

    template <typename CharT, typename CharTraits>
    friend inline std::basic_ostream<CharT, CharTraits> &operator<< (
            std::basic_ostream<CharT, CharTraits> &os,
            const GeneratorWrapper<ResultType, Generator, Traits> &)
    {return os;}

    template <typename CharT, typename CharTraits>
    friend inline std::basic_istream<CharT, CharTraits> &operator>> (
            std::basic_istream<CharT, CharTraits> &is,
            GeneratorWrapper<ResultType, Generator, Traits> &)
    {return is;}

    private :

    Generator generator_;
}; // clss GeneratorWrapper

} // namespace vsmc

#endif // VSMC_RNG_GENERATOR_WRAPPER_HPP
