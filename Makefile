CC = g++
FLAGS = 

CPPS = node.cpp nodesupp.cpp mathfunc.cpp main.cpp

EXENAME = fmm2d

$(EXENAME):
		$(CC) -o $(EXENAME) $(CPPS) $(FLAGS)
