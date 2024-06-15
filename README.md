<div id="top"></div>
<!--
*** Thanks for checking out the Best-README-Template. If you have a suggestion
*** that would make this better, please fork the repo and create a pull request
*** or simply open an issue with the tag "enhancement".
*** Don't forget to give the project a star!
*** Thanks again! Now go create something AMAZING! :D
-->


[![License](https://img.shields.io/badge/License-Apache_2.0-blue.svg)](https://opensource.org/licenses/Apache-2.0)


<!-- ABOUT THE PROJECT -->
## About The Project

The D&C order pruner is an exact structure learning algorithm for Bayesian networks, operating on the space of topological orders and using a divide and conquer technique.

<p align="right">(<a href="#top">back to top</a>)</p>


<!-- GETTING STARTED -->
# Getting Started

 Clone the repository and make

```sh
git clone https://github.com/felixleopoldo/orderpruner.git
```
Easiest way to run the algorithms is to build a Docker image from the Docker file and run the algorithm in a container.

## Docker

```bash
    $ docker build . -t dnc
```
Run the container and mount the current directory, with the repo
```bash
    $ docker run -w /mnt -v $(pwd):/mnt -it dnc
```
Inside the container run 
```bash
make
```
Now follow the instructions below

## Native

### Requirements

- Boost 1.82.0 C++ library

#### R packages
- RInside
- Rcpp
- Jmisc
- argparser
- igraph
- BiDAG 2.0.0

#### For plotting:
- ggplot2
- dplyr
- latex2exp
- patchwork

### Installation


After cloning th repo, type make 

```sh
make
```

<p align="right">(<a href="#top">back to top</a>)</p>


<!-- USAGE EXAMPLES -->
## Usage

To run the algorithm with the data file *data/p20n300gaussdata.csv* type
```sh
./run_opruner --filename data/asiadata.csv --scoretype bge --am 0.1 --aw NULL --output_csv dag_adjmat.csv
```
The estimate dag is then as an adjacency matrix in the CSV file *dag_adjmat.csv*.
If M_ij=1 in the adjacency matrix, then there is an edge i->j in the DAG.


## Benchmarks 
To generate benchmarks (with the setting specified in *R/run_opruner.R*) type

```sh
Rscript R/benchmark_opruner.R --output_dir results --filename res.csv --seeds_from 1 --seeds_to 5
```

This produces joined results in the file *res.csv* which can be analysed by typing
```R
source('R/plotting.R')
```


<p align="right">(<a href="#top">back to top</a>)</p>



<!-- LICENSE -->
## License

Distributed under the Apache 2.0 License. See `LICENSE.txt` for more information.


<p align="right">(<a href="#top">back to top</a>)</p>




<!-- MARKDOWN LINKS & IMAGES -->
<!-- https://www.markdownguide.org/basic-syntax/#reference-style-links -->
[contributors-shield]: https://img.shields.io/github/contributors/othneildrew/Best-README-Template.svg?style=for-the-badge
[contributors-url]: https://github.com/othneildrew/Best-README-Template/graphs/contributors
[forks-shield]: https://img.shields.io/github/forks/othneildrew/Best-README-Template.svg?style=for-the-badge
[forks-url]: https://github.com/othneildrew/Best-README-Template/network/members
[stars-shield]: https://img.shields.io/github/stars/othneildrew/Best-README-Template.svg?style=for-the-badge
[stars-url]: https://github.com/othneildrew/Best-README-Template/stargazers
[issues-shield]: https://img.shields.io/github/issues/othneildrew/Best-README-Template.svg?style=for-the-badge
[issues-url]: https://github.com/othneildrew/Best-README-Template/issues
[license-shield]: https://img.shields.io/github/license/othneildrew/Best-README-Template.svg?style=for-the-badge
[license-url]: https://github.com/othneildrew/Best-README-Template/blob/master/LICENSE.txt
[linkedin-shield]: https://img.shields.io/badge/-LinkedIn-black.svg?style=for-the-badge&logo=linkedin&colorB=555
