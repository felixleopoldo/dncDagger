R_LIBS=/usr/local/lib/R/site-library
R_LIBS_ARM=/Library/Frameworks/R.framework/Versions/4.3-arm64/Resources/library

CFLAGS_ARM=-Wall -pipe -Wno-unused -pedantic -Wall -L /Library/Frameworks/R.framework/Versions/4.3-arm64/Resources/lib -l R  -I /Library/Frameworks/R.framework/Versions/4.3-arm64/Resources/include/ \
			  -L $(R_LIBS_ARM)/RInside/lib/ -l RInside -Wl,-rpath,$(R_LIBS_ARM)/RInside/lib -I $(R_LIBS_ARM)/RInside/include/ -I $(R_LIBS_ARM)/Rcpp/include/ -std=c++17 

CFLAGS=-Wall -pipe -Wno-unused -pedantic -Wall -L /usr/lib/R/lib -l R -L $(R_LIBS)/RInside/lib/ -l RInside -Wl,-rpath,$(R_LIBS)/RInside/lib -I $(R_LIBS)/RInside/include/  -I $(R_LIBS)/Rcpp/include/ -I /usr/share/R/include/ -std=c++17 -fconcepts
INCLUDE=include/auxiliary.cpp include/LeftOrder.cpp include/dnc.cpp include/RightOrder.cpp include/OrderScoring.cpp include/opruner_right.cpp include/path_pruning.cpp
#INCLUDE=include/*.cpp

allopt: ; c++ -o run_opruner run_opruner.cpp $(INCLUDE) -O3 $(CFLAGS)
all: ; c++ -o run_opruner run_opruner.cpp $(INCLUDE) $(CFLAGS)
all_arm: ; c++ -o run_opruner run_opruner.cpp $(INCLUDE) $(CFLAGS_ARM)
sandbox: ; c++ -o sandbox sandbox.cpp $(INCLUDE) -O3 $(CFLAGS)
debug: ; c++ -o run_opruner run_opruner.cpp $(INCLUDE) -g $(CFLAGS)
gprof: ; c++  -o run_opruner run_opruner.cpp $(INCLUDE) -pg -no-pie -fno-builtin $(CFLAGS)
clean: ; rm -f run_opruner sandbox

