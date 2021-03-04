# Replicate Results in the JSS paper

# Introduction

This document provide detailed information on how to replicate the results
shown in the JSS paper. In the following, we will assume that working directory
start with the top level of the vSMC source tree. To begin, one need CMake
2.8.3 or later.

~~~sh
mkdir build
cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
make paper gmm
~~~

# Particle filter example

To reproduce Figure 1 in the paper, one need `R` and the `ggplot2` package.

~~~sh
cd build/example/paper
./paper_pf
Rscript paper_pf.R
~~~

# Gaussian mixture model example

To run the example and see the normalizing constant estimator,

~~~sh
cd build/example/paper
./paper_gmm_seq
~~~

To reproduce Figure 2 and 3 in the paper,

~~~sh
cd build/example/gmm
Rscript gmm_smc_benchmark.R
~~~

Below are some details.

The program used to produce the benchmark results are not exactly the same as
the one described in the paper, which is available as `paper/paper_gmm`. The
only difference is that it supports some command options which help the
automated benchmark process. Using `paper_gmm` one can also produce the same
results, but it will be quite tiresome. On the other hand, describing the
benchmark programs will make the paper unnecessarily long. So I opt to provide
one version that shows the implementation of the algorithm (in paper) and
another with additional functionalities for providing benchmarking. `paper_gmm`
is build when running 'make paper', and again, like `paper_pf`, it has minimal
dependency and thus has the best chance to work out of box.

To exactly reproduce the results, one need the following,

* Intel C++ compiler 14 or later, in C++11 mode, used with a recent version GCC
  in a recent Linux distribution.
* `libdispatch` and its dependencies on Linux
* OpenCL implementations

I have tried my best to make the automated benchmark script as smooth as
possible, but one still might run into some platform specific problems since it
will try to test all possible parallelizations (in reality, one or two are
perhaps enough for the user, for example, I myself only use TBB in work. But
different people hardly use the same one).

However, without following the above dependencies, one can get some benchmark
results, only with the backends not detected by CMake missing from the figures.

Depending on the computer and which backends are available, the benchmark can
take from half an hour to a couple of hours for the CPU part. For the GPU part,
it will take longer or shorter depending on how powerful the OpenCL hardware
is.
