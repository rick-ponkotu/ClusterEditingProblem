CC = /usr/bin/g++
CFLAGS = -I${HOME}/local/include -std=c++11 -g3 -Wall
#CPLEX = -I/Applications/CPLEX_Studio1210/cplex/include/ -I/Applications/CPLEX_Studio1210/concert/include/ -DIL_STD
#CPLEXLIB = -L/Applications/CPLEX_Studio1210/cplex/lib/x86-64_osx/static_pic -L/Applications/CPLEX_Studio1210/concert/lib/x86-64_osx/static_pic

CPLEX = -I/Applications/CPLEX_Studio201/cplex/include/ -I/Applications/CPLEX_Studio201/concert/include/ -DIL_STD
CPLEXLIB = -L/Applications/CPLEX_Studio201/cplex/lib/x86-64_osx/static_pic -L/Applications/CPLEX_Studio201/concert/lib/x86-64_osx/static_pic

# COMPILE
all: main.o graph.o solver.o reduction.o lp.o exactBranching.o heuristicBranching.o randomPivot.o lpUsedBranching.o
	$(CC) $(CPLEXLIB) -o clusteredit main.o graph.o solver.o reduction.o lp.o exactBranching.o heuristicBranching.o  randomPivot.o lpUsedBranching.o -lm -lilocplex -lcplex -lconcert

main.o: main.cpp main.h graph.h solver.h exactBranching.h
	$(CC) $(CFLAGS) -c main.cpp

graph.o: graph.cpp main.h graph.h
	$(CC) $(CFLAGS) -c graph.cpp

solver.o: solver.cpp solver.h main.h graph.h
	$(CC) $(CFLAGS) $(CPLEX) -c solver.cpp

reduction.o: reduction.cpp reduction.h main.h graph.h
	$(CC) $(CFLAGS) -c reduction.cpp

lp.o: lp.cpp lp.h main.h graph.h
	$(CC) $(CFLAGS) $(CPLEX) -c lp.cpp

exactBranching.o: exactBranching.cpp exactBranching.h main.h graph.h
	$(CC) $(CFLAGS) $(CPLEX) -c exactBranching.cpp

heuristicBranching.o: heuristicBranching.cpp heuristicBranching.h main.h graph.h
	$(CC) $(CFLAGS) -c heuristicBranching.cpp

randomPivot.o: randomPivot.cpp randomPivot.o main.h graph.h
	$(CC) $(CFLAGS) -c randomPivot.cpp

lpUsedBranching.o: lpUsedBranching.cpp lpUsedBranching.h main.h graph.h randomPivot.h
	$(CC) $(CFLAGS) -c lpUsedBranching.cpp

clean:
	$(RM) -f  *.o

test: all
	./clusteredit

