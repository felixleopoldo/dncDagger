FROM bpimages/bidag:2.1.3

RUN R -e "install.packages(\"RInside\", repos=\"https://cran.rstudio.com\")" --no-save
RUN R -e "install.packages(\"Jmisc\")" --no-save
RUN R -e "install.packages(\"argparser\")" --no-save

WORKDIR /orderpruner
COPY . .
RUN make

#RUN R -e "install.packages(\"ggplot2\")" --no-save
#RUN R -e "install.packages(c(\"testit\"))" --no-save
#RUN R -e "install.packages(c(\"dplyr\"))" --no-save

#RUN apt-get -y update
#RUN apt-get install gdb -y