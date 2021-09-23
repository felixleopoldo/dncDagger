library("Rcpp")
library("Jmisc")
sourceAll(path = "R")
filename <- 'data/myasiandata.csv.rds'

ret <- readRDS(filename); 

Sys.setenv("PKG_CXXFLAGS"="-Wall -pipe -Wno-unused -pedantic -Wall -L /usr/lib/R/lib -l R  -pthread -lpcre -llzma -lbz2 -lz -lrt -ldl -lm -L /usr/local/lib/R/site-library/RInside/lib/ -l RInside -Wl,-rpath,/usr/local/lib/R/site-library/RInside/lib  -llapack -lblas  -I /usr/local/lib/R/site-library/RInside/include/ -I /usr/local/lib/R/site-library/Rcpp/include/ -I /usr/share/R/include/ -I thread-pool-master/ -std=c++17 -O3")
sourceCpp("smc_stuff.cpp", verbose=TRUE)

a <- r_mh(ret, 10000000, 1)

saveRDS(a, "mh_output.rds")
