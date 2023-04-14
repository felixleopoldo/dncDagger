CFLAGS=-Wall -pipe -Wno-unused -pedantic -Wall -L /usr/lib/R/lib -l R -L /usr/local/lib/R/site-library/RInside/lib/ -l RInside -Wl,-rpath,/usr/local/lib/R/site-library/RInside/lib -I /usr/local/lib/R/site-library/RInside/include/ -I /usr/local/lib/R/site-library/Rcpp/include/ -I /usr/share/R/include/ -std=c++17
INCLUDE=include/auxiliary.cpp include/RightOrder.cpp include/OrderScoring.cpp include/opruner_right.cpp include/opruner_left.cpp

all: ; c++ -o run_opruner run_opruner.cpp $(INCLUDE) -O3 $(CFLAGS)
debug: ; c++ -o run_opruner run_opruner.cpp $(INCLUDE) -g $(CFLAGS)
gprof: ; c++  -o run_opruner run_opruner.cpp $(INCLUDE) -pg -no-pie -fno-builtin $(CFLAGS)
