CXXFLAGS = -g -std=c++11 -I. Eigen #After update method is implemented and debugged, use these flags for optimization: CXXFLAGS = -O3 -std=c++11 -I. Eigen
PROGS = metropolis compute_observables
all: $(PROGS) 
metropolis: su3.o
compute_observables: su3.o
.PHONY: clean
clean:
	${RM} su3.o $(PROGS)


