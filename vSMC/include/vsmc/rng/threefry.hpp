//============================================================================
// vSMC/include/vsmc/rng/threefry.hpp
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

#ifndef VSMC_RNG_THREEFRY_HPP
#define VSMC_RNG_THREEFRY_HPP

#include <vsmc/rng/internal/common.hpp>

#define VSMC_STATIC_ASSERT_RNG_THREEFRY_RESULT_TYPE(ResultType) \
    VSMC_STATIC_ASSERT(                                                      \
            (cxx11::is_same<ResultType, uint32_t>::value ||                  \
             cxx11::is_same<ResultType, uint64_t>::value),                   \
            USE_ThreefryEngine_WITH_INTEGER_TYPE_OTHER_THAN_uint32_t_OR_uint64_t)

#define VSMC_STATIC_ASSERT_RNG_THREEFRY_SIZE(K) \
    VSMC_STATIC_ASSERT((K == 2 || K == 4),                                   \
            USE_ThreefryEngine_WITH_SIZE_OTHER_THAN_2_OR_4)

#define VSMC_STATIC_ASSERT_RNG_THREEFRY \
        VSMC_STATIC_ASSERT_RNG_THREEFRY_RESULT_TYPE(ResultType);             \
        VSMC_STATIC_ASSERT_RNG_THREEFRY_SIZE(K);

#define VSMC_DEFINE_RNG_THREEFRY_ROTATE_CONSTANT(T, K, N, I, val) \
    template <> struct ThreefryRotateConstantValue < T, K, N, I > :          \
        public cxx11::integral_constant< unsigned, val > {};

/// \brief ThreefryEngine default rounds
/// \ingroup Config
#ifndef VSMC_RNG_THREEFRY_ROUNDS
#define VSMC_RNG_THREEFRY_ROUNDS 20
#endif

namespace vsmc {

namespace internal {

template <typename> struct ThreefryKSConstantValue;

template <> struct ThreefryKSConstantValue<uint32_t> :
    public cxx11::integral_constant<uint32_t,
           UINT32_C(0x1BD11BDA)> {};

template <> struct ThreefryKSConstantValue<uint64_t> :
    public cxx11::integral_constant<uint64_t,
           UINT64_C(0x1BD11BDAA9FC1A22)> {};

template <typename, std::size_t, std::size_t, std::size_t>
struct ThreefryRotateConstantValue;

VSMC_DEFINE_RNG_THREEFRY_ROTATE_CONSTANT(uint32_t, 2, 0, 0, 13)
VSMC_DEFINE_RNG_THREEFRY_ROTATE_CONSTANT(uint32_t, 2, 1, 0, 15)
VSMC_DEFINE_RNG_THREEFRY_ROTATE_CONSTANT(uint32_t, 2, 2, 0, 26)
VSMC_DEFINE_RNG_THREEFRY_ROTATE_CONSTANT(uint32_t, 2, 3, 0,  6)
VSMC_DEFINE_RNG_THREEFRY_ROTATE_CONSTANT(uint32_t, 2, 4, 0, 17)
VSMC_DEFINE_RNG_THREEFRY_ROTATE_CONSTANT(uint32_t, 2, 5, 0, 29)
VSMC_DEFINE_RNG_THREEFRY_ROTATE_CONSTANT(uint32_t, 2, 6, 0, 16)
VSMC_DEFINE_RNG_THREEFRY_ROTATE_CONSTANT(uint32_t, 2, 7, 0, 24)

VSMC_DEFINE_RNG_THREEFRY_ROTATE_CONSTANT(uint32_t, 4, 0, 0, 10)
VSMC_DEFINE_RNG_THREEFRY_ROTATE_CONSTANT(uint32_t, 4, 1, 0, 11)
VSMC_DEFINE_RNG_THREEFRY_ROTATE_CONSTANT(uint32_t, 4, 2, 0, 13)
VSMC_DEFINE_RNG_THREEFRY_ROTATE_CONSTANT(uint32_t, 4, 3, 0, 23)
VSMC_DEFINE_RNG_THREEFRY_ROTATE_CONSTANT(uint32_t, 4, 4, 0,  6)
VSMC_DEFINE_RNG_THREEFRY_ROTATE_CONSTANT(uint32_t, 4, 5, 0, 17)
VSMC_DEFINE_RNG_THREEFRY_ROTATE_CONSTANT(uint32_t, 4, 6, 0, 25)
VSMC_DEFINE_RNG_THREEFRY_ROTATE_CONSTANT(uint32_t, 4, 7, 0, 18)

VSMC_DEFINE_RNG_THREEFRY_ROTATE_CONSTANT(uint32_t, 4, 0, 1, 26)
VSMC_DEFINE_RNG_THREEFRY_ROTATE_CONSTANT(uint32_t, 4, 1, 1, 21)
VSMC_DEFINE_RNG_THREEFRY_ROTATE_CONSTANT(uint32_t, 4, 2, 1, 27)
VSMC_DEFINE_RNG_THREEFRY_ROTATE_CONSTANT(uint32_t, 4, 3, 1,  5)
VSMC_DEFINE_RNG_THREEFRY_ROTATE_CONSTANT(uint32_t, 4, 4, 1, 20)
VSMC_DEFINE_RNG_THREEFRY_ROTATE_CONSTANT(uint32_t, 4, 5, 1, 11)
VSMC_DEFINE_RNG_THREEFRY_ROTATE_CONSTANT(uint32_t, 4, 6, 1, 10)
VSMC_DEFINE_RNG_THREEFRY_ROTATE_CONSTANT(uint32_t, 4, 7, 1, 20)

VSMC_DEFINE_RNG_THREEFRY_ROTATE_CONSTANT(uint64_t, 2, 0, 0, 16)
VSMC_DEFINE_RNG_THREEFRY_ROTATE_CONSTANT(uint64_t, 2, 1, 0, 42)
VSMC_DEFINE_RNG_THREEFRY_ROTATE_CONSTANT(uint64_t, 2, 2, 0, 12)
VSMC_DEFINE_RNG_THREEFRY_ROTATE_CONSTANT(uint64_t, 2, 3, 0, 31)
VSMC_DEFINE_RNG_THREEFRY_ROTATE_CONSTANT(uint64_t, 2, 4, 0, 16)
VSMC_DEFINE_RNG_THREEFRY_ROTATE_CONSTANT(uint64_t, 2, 5, 0, 32)
VSMC_DEFINE_RNG_THREEFRY_ROTATE_CONSTANT(uint64_t, 2, 6, 0, 24)
VSMC_DEFINE_RNG_THREEFRY_ROTATE_CONSTANT(uint64_t, 2, 7, 0, 21)

VSMC_DEFINE_RNG_THREEFRY_ROTATE_CONSTANT(uint64_t, 4, 0, 0, 14)
VSMC_DEFINE_RNG_THREEFRY_ROTATE_CONSTANT(uint64_t, 4, 1, 0, 52)
VSMC_DEFINE_RNG_THREEFRY_ROTATE_CONSTANT(uint64_t, 4, 2, 0, 23)
VSMC_DEFINE_RNG_THREEFRY_ROTATE_CONSTANT(uint64_t, 4, 3, 0,  5)
VSMC_DEFINE_RNG_THREEFRY_ROTATE_CONSTANT(uint64_t, 4, 4, 0, 25)
VSMC_DEFINE_RNG_THREEFRY_ROTATE_CONSTANT(uint64_t, 4, 5, 0, 46)
VSMC_DEFINE_RNG_THREEFRY_ROTATE_CONSTANT(uint64_t, 4, 6, 0, 58)
VSMC_DEFINE_RNG_THREEFRY_ROTATE_CONSTANT(uint64_t, 4, 7, 0, 32)

VSMC_DEFINE_RNG_THREEFRY_ROTATE_CONSTANT(uint64_t, 4, 0, 1, 16)
VSMC_DEFINE_RNG_THREEFRY_ROTATE_CONSTANT(uint64_t, 4, 1, 1, 57)
VSMC_DEFINE_RNG_THREEFRY_ROTATE_CONSTANT(uint64_t, 4, 2, 1, 40)
VSMC_DEFINE_RNG_THREEFRY_ROTATE_CONSTANT(uint64_t, 4, 3, 1, 37)
VSMC_DEFINE_RNG_THREEFRY_ROTATE_CONSTANT(uint64_t, 4, 4, 1, 33)
VSMC_DEFINE_RNG_THREEFRY_ROTATE_CONSTANT(uint64_t, 4, 5, 1, 12)
VSMC_DEFINE_RNG_THREEFRY_ROTATE_CONSTANT(uint64_t, 4, 6, 1, 22)
VSMC_DEFINE_RNG_THREEFRY_ROTATE_CONSTANT(uint64_t, 4, 7, 1, 32)

template <typename ResultType, unsigned N> struct ThreefryRotateImpl;

template <unsigned N>
struct ThreefryRotateImpl<uint32_t, N>
{static uint32_t eval (uint32_t x) {return x << N | x >> (32 - N);}};

template <unsigned N>
struct ThreefryRotateImpl<uint64_t, N>
{static uint64_t eval (uint64_t x) {return x << N | x >> (64 - N);}};

template <typename ResultType, std::size_t K, std::size_t N, bool = (N > 0)>
struct ThreefryRotate {static void eval (Array<ResultType, K> &) {}};

template <typename ResultType, std::size_t N>
struct ThreefryRotate<ResultType, 2, N, true>
{
    static void eval (Array<ResultType, 2> &state)
    {
        state[Position<0>()] += state[Position<1>()];
        state[Position<1>()] =
            ThreefryRotateImpl<ResultType, ThreefryRotateConstantValue<
            ResultType, 2, r_, 0>::value>::eval(state[Position<1>()]);
        state[Position<1>()] ^= state[Position<0>()];
    }

    private :

    static VSMC_CONSTEXPR const unsigned r_ = (N - 1) % 8;
}; // struct ThreefryRotate

template <typename ResultType, std::size_t N>
struct ThreefryRotate<ResultType, 4, N, true>
{
    static void eval (Array<ResultType, 4> &state)
    {
        state[Position<0>()] += state[Position<i0_>()];
        state[Position<i0_>()] =
            ThreefryRotateImpl<ResultType, ThreefryRotateConstantValue<
            ResultType, 4, r_, 0>::value>::eval(state[Position<i0_>()]);
        state[Position<i0_>()] ^= state[Position<0>()];

        state[Position<2>()] += state[Position<i2_>()];
        state[Position<i2_>()] =
            ThreefryRotateImpl<ResultType, ThreefryRotateConstantValue<
            ResultType, 4, r_, 1>::value>::eval(state[Position<i2_>()]);
        state[Position<i2_>()] ^= state[Position<2>()];
    }

    private :

    static VSMC_CONSTEXPR const std::size_t i0_ = N % 2 ? 1 : 3;
    static VSMC_CONSTEXPR const std::size_t i2_ = N % 2 ? 3 : 1;
    static VSMC_CONSTEXPR const unsigned r_ = (N - 1) % 8;
}; // struct ThreefryRotate

template <typename ResultType, std::size_t K, std::size_t N,
         bool = (N % 4 == 0)>
struct ThreefryInsertKey
{
    static void eval (Array<ResultType, K> &,
            const Array<ResultType, K + 1> &) {}
}; // struct ThreefryInsertKey

template <typename ResultType, std::size_t N>
struct ThreefryInsertKey<ResultType, 2, N, true>
{
    static void eval (Array<ResultType, 2> &state,
            const Array<ResultType, 3> &par)
    {
        state[Position<0>()] += par[Position<i0_>()];
        state[Position<1>()] += par[Position<i1_>()];
        state[Position<1>()] += inc_;
    }

    private :

    static VSMC_CONSTEXPR const std::size_t inc_ = N / 4;
    static VSMC_CONSTEXPR const std::size_t i0_ = (inc_ + 0) % 3;
    static VSMC_CONSTEXPR const std::size_t i1_ = (inc_ + 1) % 3;
}; // struct ThreefryInsertKey

template <typename ResultType, std::size_t N>
struct ThreefryInsertKey<ResultType, 4, N, true>
{
    static void eval (Array<ResultType, 4> &state,
            const Array<ResultType, 5> &par)
    {
        state[Position<0>()] += par[Position<i0_>()];
        state[Position<1>()] += par[Position<i1_>()];
        state[Position<2>()] += par[Position<i2_>()];
        state[Position<3>()] += par[Position<i3_>()];
        state[Position<3>()] += inc_;
    }

    private :

    static VSMC_CONSTEXPR const std::size_t inc_ = N / 4;
    static VSMC_CONSTEXPR const std::size_t i0_ = (inc_ + 0) % 5;
    static VSMC_CONSTEXPR const std::size_t i1_ = (inc_ + 1) % 5;
    static VSMC_CONSTEXPR const std::size_t i2_ = (inc_ + 2) % 5;
    static VSMC_CONSTEXPR const std::size_t i3_ = (inc_ + 3) % 5;
}; // struct ThreefryInsertKey

} // namespace vsmc::internal

/// \brief Threefry RNG engine reimplemented
/// \ingroup R123RNG
///
/// \details
/// This is a reimplementation of the algorithm Threefry as described in
/// [Parallel Random Numbers: As Easy as 1, 2, 3][r123paper] and implemented in
/// [Random123][r123lib].
///
/// [r123paper]:http://sc11.supercomputing.org/schedule/event_detail.php?evid=pap274
/// [r123lib]: https://www.deshawresearch.com/resources_random123.html
///
/// Depending on the compilers, processors and RNG configurations, it might be
/// slightly faster or slower than the original implementation. At most
/// two-folds performace difference (both faster and slower) were observed.
///
/// This implementation is slightly more flexible in the sense that it does not
/// limit the number of rounds. However, larger number of rounds can have
/// undesired effects. To say the least, currently all loops are unrolled,
/// which can slow down significantly when the number of rounds is large.
///
/// Compared to `r123:Engine<r123::Threefry4x32>` etc., when using the default
/// constructor or the one with a single seed, the output shall be exactly the
/// same for the first \f$2^n\f$ iterations, where \f$n\f$ is the number of
/// bits (32 or 64).  Further iterations may produce different results, as vSMC
/// increment the counter slightly differently, but it still cover the same
/// range and has the same period as the original.
template <typename ResultType, std::size_t K,
         std::size_t Rounds = VSMC_RNG_THREEFRY_ROUNDS>
class ThreefryEngine
{
    public :

    typedef ResultType result_type;
    typedef Array<ResultType, K> buffer_type;
    typedef Array<ResultType, K> ctr_type;
    typedef Array<ResultType, K> key_type;

    private :

    typedef Counter<ctr_type> counter;

    public :

    explicit ThreefryEngine (result_type s = 0) : index_(K)
    {
        VSMC_STATIC_ASSERT_RNG_THREEFRY;
        seed(s);
    }

    template <typename SeedSeq>
    explicit ThreefryEngine (SeedSeq &seq,
            typename cxx11::enable_if<internal::is_seed_seq<SeedSeq,
            result_type, key_type, ThreefryEngine<ResultType, K, Rounds>
            >::value>::type * = VSMC_NULLPTR) : index_(K)
    {
        VSMC_STATIC_ASSERT_RNG_THREEFRY;
        seed(seq);
    }

    ThreefryEngine (const key_type &k) : index_(K)
    {
        VSMC_STATIC_ASSERT_RNG_THREEFRY;
        seed(k);
    }

    void seed (result_type s)
    {
        counter::reset(ctr_);
        key_type k;
        k.fill(0);
        k.front() = s;
        init_par(k);
        index_ = K;
    }

    template <typename SeedSeq>
    void seed (SeedSeq &seq,
            typename cxx11::enable_if<internal::is_seed_seq<SeedSeq,
            result_type, key_type, ThreefryEngine<ResultType, K, Rounds>
            >::value>::type * = VSMC_NULLPTR)
    {
        counter::reset(ctr_);
        key_type k;
        seq.generate(k.begin(), k.end());
        init_par(k);
        index_ = K;
    }

    void seed (const key_type &k)
    {
        counter::reset(ctr_);
        init_par(k);
        index_ = K;
    }

    ctr_type ctr () const {return ctr_;}

    key_type key () const
    {
        key_type k;
        for (std::size_t i = 0; i != K; ++i)
            k[i] = par_[i];

        return k;
    }

    void ctr (const ctr_type &c)
    {
        counter::set(ctr_, c);
        index_ = K;
    }

    void key (const key_type &k)
    {
        init_par(k);
        index_ = K;
    }

    result_type operator() ()
    {
        if (index_ == K) {
            counter::increment(ctr_);
            generate_buffer(ctr_, buffer_);
            index_ = 0;
        }

        return buffer_[index_++];
    }

    /// \brief Generate a buffer of random bits given a counter using the
    /// current key
    buffer_type operator() (const ctr_type &c) const
    {
        buffer_type buf;
        generate_buffer(c, buf);

        return buf;
    }

    /// \brief Generate random bits in a pre-allocated buffer given a counter
    /// using the current key
    void operator() (const ctr_type &c, buffer_type &buf) const
    {generate_buffer(c, buf);}

    void discard (result_type nskip)
    {
        std::size_t n = static_cast<std::size_t>(nskip);
        if (index_ + n <= K) {
            index_ += n;
            return;
        }

        n -= K - index_;
        if (n <= K) {
            index_ = K;
            operator()();
            index_ = n;
            return;
        }

        counter::increment(ctr_, static_cast<result_type>(n / K));
        index_ = K;
        operator()();
        index_ = n % K;
    }

    static VSMC_CONSTEXPR const result_type _Min = 0;
    static VSMC_CONSTEXPR const result_type _Max = static_cast<result_type>(
            ~(static_cast<result_type>(0)));

    static VSMC_CONSTEXPR result_type min VSMC_MNE () {return _Min;}
    static VSMC_CONSTEXPR result_type max VSMC_MNE () {return _Max;}

    friend inline bool operator== (
            const ThreefryEngine<ResultType, K, Rounds> &eng1,
            const ThreefryEngine<ResultType, K, Rounds> &eng2)
    {
        return
            eng1.index_ == eng2.index_ &&
            eng1.ctr_ == eng2.ctr_ &&
            eng1.par_ == eng2.par_;
    }

    friend inline bool operator!= (
            const ThreefryEngine<ResultType, K, Rounds> &eng1,
            const ThreefryEngine<ResultType, K, Rounds> &eng2)
    {return !(eng1 == eng2);}

    template <typename CharT, typename Traits>
    friend inline std::basic_ostream<CharT, Traits> &operator<< (
            std::basic_ostream<CharT, Traits> &os,
            const ThreefryEngine<ResultType, K, Rounds> &eng)
    {
        if (!os.good())
            return os;

        os << eng.ctr_ << ' ';
        os << eng.par_ << ' ';
        os << eng.buffer_ << ' ';
        os << eng.index_;

        return os;
    }

    template <typename CharT, typename Traits>
    friend inline std::basic_istream<CharT, Traits> &operator>> (
            std::basic_istream<CharT, Traits> &is,
            ThreefryEngine<ResultType, K, Rounds> &eng)
    {
        if (!is.good())
            return is;

        ThreefryEngine<ResultType, K, Rounds> eng_tmp;
        is >> std::ws >> eng_tmp.ctr_;
        is >> std::ws >> eng_tmp.par_;
        is >> std::ws >> eng_tmp.buffer_;
        is >> std::ws >> eng_tmp.index_;

        if (is.good()) {
#if VSMC_HAS_CXX11_RVALUE_REFERENCES
            eng = cxx11::move(eng_tmp);
#else
            eng = eng_tmp;
#endif
        }

        return is;
    }

    private :

    ctr_type ctr_;
    Array<ResultType, K + 1> par_;
    buffer_type buffer_;
    std::size_t index_;

    void generate_buffer (const ctr_type &c, buffer_type &buf) const
    {
        buf = c;
        generate_buffer<0>(buf, cxx11::true_type());
    }

    template <std::size_t>
    void generate_buffer (buffer_type &, cxx11::false_type) const {}

    template <std::size_t N>
    void generate_buffer (buffer_type &buf, cxx11::true_type) const
    {
        internal::ThreefryRotate<ResultType, K, N>::eval(buf);
        internal::ThreefryInsertKey<ResultType, K, N>::eval(buf, par_);
        generate_buffer<N + 1>(buf,
                cxx11::integral_constant<bool, N < Rounds>());
    }

    void init_par (const key_type &key)
    {
        par_.back() = internal::ThreefryKSConstantValue<ResultType>::value;
        par_xor<0>(key, cxx11::integral_constant<bool, 0 < K>());
    }

    template <std::size_t>
    void par_xor (const key_type &, cxx11::false_type) {}

    template <std::size_t N>
    void par_xor (const key_type &key, cxx11::true_type)
    {
        par_[Position<N>()] = key[Position<N>()];
        par_.back() ^= key[Position<N>()];
        par_xor<N + 1>(key, cxx11::integral_constant<bool, N + 1 < K>());
    }
}; // class ThreefryEngine

/// \brief Threefry2x32 RNG engine reimplemented
/// \ingroup R123RNG
typedef ThreefryEngine<uint32_t, 2> Threefry2x32;

/// \brief Threefry4x32 RNG engine reimplemented
/// \ingroup R123RNG
typedef ThreefryEngine<uint32_t, 4> Threefry4x32;

/// \brief Threefry2x64 RNG engine reimplemented
/// \ingroup R123RNG
typedef ThreefryEngine<uint64_t, 2> Threefry2x64;

/// \brief Threefry4x64 RNG engine reimplemented
/// \ingroup R123RNG
typedef ThreefryEngine<uint64_t, 4> Threefry4x64;

/// \brief The default 32-bits Threefry engine
/// \ingroup R123RNG
typedef Threefry4x32 Threefry;

/// \brief The default 64-bits Threefry engine
/// \ingroup R123RNG
typedef Threefry4x64 Threefry_64;

} // namespace vsmc

#endif // VSMC_RNG_THREEFRY_HPP
