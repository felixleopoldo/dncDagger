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


RUN wget https://github.com/apptainer/apptainer/releases/download/v1.3.1/apptainer_1.3.1_amd64.deb
RUN apt install -y ./apptainer_1.3.1_amd64.deb
# Setuid install
#RUN echo setuid installation
#RUN wget https://github.com/apptainer/apptainer/releases/download/v1.3.1/apptainer-suid_1.3.1_amd64.deb
#RUN dpkg -i ./apptainer-suid_1.3.1_amd64.deb

# Fix namespace creation. Needed? 
#RUN sh -c 'echo kernel.unprivileged_userns_clone=1 >/etc/sysctl.d/90-unprivileged_userns.conf'
#RUN sysctl -p /etc/sysctl.d /etc/sysctl.d/90-unprivileged_userns.conf

# Order pruner:
WORKDIR /orderpruner
COPY . .
# run make on the singularity container
