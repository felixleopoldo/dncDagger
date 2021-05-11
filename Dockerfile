FROM debian:stretch 

RUN apt-get -y update && apt-get install cmake -y
RUN apt-get install -y nasm
RUN apt-get install libblas-dev liblapack-dev -y
RUN apt-get install gdb -y

#RUN apt-get install libboost-all-dev -y
RUN apt install r-base
RUN "R install RInside"

COPY ./MCKL/include/mckl/ /usr/local/include/mckl
COPY ./vSMC/include/vsmc/ /usr/local/include/vsmcß

ENV CXXFLAGS -std=c++14
#WORKDIR /usr/src/myapp/MCKL/build

WORKDIR /usr/src/myapp 

#RUN apt-get install build-essential

LABEL Name=ordersmc Version=1.0
