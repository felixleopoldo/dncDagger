R_LIBS=/usr/local/lib/R/site-library

CFLAGS=-Wall -pipe -Wno-unused -pedantic -Wall -L /usr/lib/R/lib -l R -L $(R_LIBS)/RInside/lib/ -l RInside -Wl,-rpath,$(R_LIBS)/RInside/lib -I $(R_LIBS)/RInside/include/  -I $(R_LIBS)/Rcpp/include/ -I /usr/share/R/include/ -std=c++17 -fconcepts
INCLUDE=include/auxiliary.cpp include/LeftOrder.cpp include/dnc.cpp include/RightOrder.cpp include/OrderScoring.cpp include/opruner_right.cpp include/path_pruning.cpp

allnoopt: ; c++ -o dncpruner dncpruner.cpp $(INCLUDE) $(CFLAGS)
all: ; c++ -o dncpruner dncpruner.cpp $(INCLUDE) -03 $(CFLAGS)
all_arm: ; c++ -o dncpruner dncpruner.cpp $(INCLUDE) $(CFLAGS_ARM)
sandbox: ; c++ -o sandbox sandbox.cpp $(INCLUDE) -O3 $(CFLAGS)
debug: ; c++ -o dncpruner dncpruner.cpp $(INCLUDE) -g $(CFLAGS)
gprof: ; c++  -o dncpruner dncpruner.cpp $(INCLUDE) -pg -no-pie -fno-builtin $(CFLAGS)
clean: ; rm -f dncpruner sandbox
