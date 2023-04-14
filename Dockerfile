FROM r-base:4.0.5

RUN apt-get -y update
RUN apt-get -y install cmake
RUN apt-get install -y nasm 
RUN apt-get install libblas-dev liblapack-dev -y
RUN apt-get install gdb -y

RUN R -e "install.packages(\"RInside\", repos=\"https://cran.rstudio.com\")" --no-save
RUN R -e "if (!requireNamespace(\"BiocManager\", quietly = TRUE)) {install.packages(\"BiocManager\", repo=\"http://cran.rstudio.com/\")}" --no-save
RUN R -e "BiocManager::install()"  --no-save
RUN R -e "BiocManager::install(c(\"gRbase\", \"RBGL\", \"Rgraphviz\", \"gRain\"))" --no-save
RUN R -e "install.packages(\"pcalg\", repos=\"https://cran.rstudio.com\")" --no-save
RUN R -e "packageurl <- \"https://cran.r-project.org/src/contrib/Archive/BiDAG/BiDAG_2.0.0.tar.gz\" ; install.packages(packageurl, repos=NULL, type=\"source\")" --no-save
RUN R -e "install.packages(c(\"pcalg\", \"ggplot2\", \"Jmisc\", \"argparser\"))" --no-save
RUN R -e "install.packages(c(\"testit\"))" --no-save
RUN R -e "install.packages(c(\"dplyr\"))" --no-save
WORKDIR /order_pruning

RUN make