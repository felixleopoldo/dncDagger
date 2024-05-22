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

# Installing Apptainer, could be useful if Docker is not available.
RUN wget https://github.com/apptainer/apptainer/releases/download/v1.3.1/apptainer_1.3.1_amd64.deb
RUN apt install -y ./apptainer_1.3.1_amd64.deb

# install some R pakcages fro the plotting
RUN R -e "install.packages(\"ggplot2\")" 
RUN R -e "install.packages(\"dplyr\")" 
RUN R -e "install.packages(\"latex2exp\")" 
RUN R -e "install.packages(\"patchwork\")" 
# RUN R -e "install.packages(\"quantreg\")" # not compatible with R 4.2.3 used here

# set default servers for apptainer
RUN apptainer remote add --no-login SylabsCloud cloud.sycloud.io

# Order pruner:
#WORKDIR /orderpruner
#COPY . .
  
# can run make on the singularity/apptainer container

  