#CFLAGS=-I/usr/share/R/include/ -I/usr/local/lib/R/site-library/Rcpp/include/ -I/usr/local/lib/R/site-library/VineCopula/include -dynamiclib -Wl,-headerpad_max_install_names -shared -L/usr/lib/R/lib -lR 
#CFLAGS2=-I/usr/share/R/include/ -I/usr/local/lib/R/site-library/Rcpp/include/ -I/usr/local/lib/R/site-library/VineCopula/include 
CFLAGS=-Wall -pipe -Wno-unused -pedantic -Wall -L /usr/lib/R/lib -l R  -pthread -lpcre -llzma -lbz2 -lz -lrt -ldl -lm -L /usr/local/lib/R/site-library/RInside/lib/ -l RInside -Wl,-rpath,/usr/local/lib/R/site-library/RInside/lib  -llapack -lblas  -I /usr/local/lib/R/site-library/RInside/include/ -I /usr/local/lib/R/site-library/Rcpp/include/ -I /usr/share/R/include/ -I thread-pool-master/ -std=c++17

#LDFLAGS=-DNDEBUG -fpic  -g -O2 -fstack-protector-strong -Wformat -Werror=format-security -Wdate-time -D_FORTIFY_SOURCE=2 -g


all: ; c++ -DMCKL_USE_ASM_LIB=1 -o run_opruner run_opruner.cpp includes/auxiliary.cpp includes/RightOrder.cpp includes/OrderScoring.cpp -O3 $(CFLAGS)
debug: ; c++ -DMCKL_USE_ASM_LIB=1 -o run_opruner run_opruner.cpp -g $(CFLAGS)
gprof: ; c++ -DMCKL_USE_ASM_LIB=1 -o run_opruner run_opruner.cpp -pg -no-pie -fno-builtin $(CFLAGS)
sandbox: ; c++ -o main main.cpp  -g -Wall
docker_build: ;  docker build -t onceltuca/orderpruning .
docker_push: ;  docker push onceltuca/orderpruning
singu_build: ; singularity pull --dir ~/singularity_images/ docker://onceltuca/orderpruning:1.2.1
sung_push: ; singularity push -U ~/singularity_images/orderpruning_1.2.1.sif library://felixleopoldo/bn/orderpruning:1.2.1