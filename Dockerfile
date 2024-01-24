FROM bpimages/bidag:2.1.3

RUN R -e "install.packages(\"RInside\", repos=\"https://cran.rstudio.com\")"
RUN R -e "install.packages(\"Jmisc\")" 
RUN R -e "install.packages(\"argparser\")"
RUN R -e "install.packages(\"igraph\")" 

RUN wget https://boostorg.jfrog.io/artifactory/main/release/1.82.0/source/boost_1_82_0.tar.gz
RUN tar xvf boost_1_82_0.tar.gz
WORKDIR boost_1_82_0
RUN ./bootstrap.sh --prefix=/usr/
RUN ./b2 install

RUN apt-get update
RUN apt install git -y
RUN git clone https://github.com/atofigh/edmonds-alg.git

WORKDIR /orderpruner
COPY . .
# RUN make



#RUN R -e "install.packages(\"ggplot2\")" --no-save
#RUN R -e "install.packages(c(\"testit\"))" --no-save
#RUN R -e "install.packages(c(\"dplyr\"))" --no-save

#RUN apt-get -y update
#RUN apt-get install gdb -y