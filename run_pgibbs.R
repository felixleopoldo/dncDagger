library("Rcpp")
library("Jmisc")
sourceAll(path = "R")
ret <- readRDS('jackdata.csv.rds'); 

Sys.setenv("PKG_CXXFLAGS"="-Wall -pipe -Wno-unused -pedantic -Wall -L /usr/lib/R/lib -l R  -pthread -lpcre -llzma -lbz2 -lz -lrt -ldl -lm -L /usr/local/lib/R/site-library/RInside/lib/ -l RInside -Wl,-rpath,/usr/local/lib/R/site-library/RInside/lib  -llapack -lblas  -I /usr/local/lib/R/site-library/RInside/include/ -I /usr/local/lib/R/site-library/Rcpp/include/ -I /usr/share/R/include/ -I thread-pool-master/ -std=c++17 -O3")
sourceCpp("smc_stuff.cpp", verbose=TRUE)

a <- r_pgibbs(ret, 10, 10000)

saveRDS(a, "pgibbs_output.rds")

#scoretable
#vars

#omcmcres <- orderMCMC(myscore, scoreout = TRUE, chainout = TRUE, plus1 = TRUE, MAP = FALSE,
#                      startorder = startorder, scoretable = scoretable,
#                      startspace = startspace, iterations = 1000000)
#omcmcres$maxorder
#omcmcres$score
#plot(omcmcres$trace, type="l")

#max(omcmcres$traceadd$orderscores)
#plot(omcmcres$traceadd$orderscores[600:1000])
#print(as.integer(omcmcres$maxorder) - 1)
#print(omcmcres$score)

# omcmc_orderindex <- c()
# for (v in omcmcres$maxorder) {
#   omcmc_orderindex <- c(omcmc_orderindex, which(vars == v)[[1]] - 1)
# }
# omcmc_orderindex
# smcorder <- c(32, 48, 37, 50, 8, 7, 27, 55, 52, 28, 66, 53, 25, 33, 11, 56, 12, 69, 31, 13, 49, 62, 6, 38, 30, 21, 46, 34, 24, 58, 2, 67, 20, 65, 63, 60, 40, 57, 36, 61, 64, 5, 18, 44, 35, 14, 9, 1, 4, 0, 39, 29, 17, 59, 3, 22, 54, 43, 45, 68, 10, 26, 47, 42, 16, 19, 51, 15, 41, 23)+1

# vars[smcorder]


