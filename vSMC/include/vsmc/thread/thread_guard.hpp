//============================================================================
// vSMC/include/vsmc/thread/thread_guard.hpp
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

#ifndef VSMC_THREAD_THREAD_GUARD_HPP
#define VSMC_THREAD_THREAD_GUARD_HPP

#include <vsmc/internal/common.hpp>

namespace vsmc {

/// \brief Strictly scope-based thread ownership wrapper
/// \ingroup Thread
template <typename ThreadType>
class ThreadGuard
{
    public :

    typedef ThreadType thread_type;

    ThreadGuard () VSMC_NOEXCEPT {}

    ThreadGuard (const ThreadGuard &) = delete;

    ThreadGuard &operator= (const ThreadGuard &) = delete;

    ThreadGuard (thread_type &&thr) VSMC_NOEXCEPT :
        thread_(cxx11::move(thr)) {}

    ThreadGuard (ThreadGuard &&other) VSMC_NOEXCEPT :
        thread_(cxx11::move(other.thread_)) {}

    ThreadGuard &operator= (ThreadGuard &&other) VSMC_NOEXCEPT
    {thread_ = cxx11::move(other.thread_); return *this;}

    ~ThreadGuard () VSMC_NOEXCEPT {if (thread_.joinable()) thread_.join();}

    private :

    thread_type thread_;
}; // class ThreadGuard

} // namespace vsmc

#endif // VSMC_THREAD_THREAD_GUARD_HPP
