CC = g++
FLAGS = 

CPPS = node.cpp mathfunc.cpp main.cpp

EXENAME = fmm2d

$(EXENAME):
		$(CC) -o $(EXENAME) $(CPPS) $(FLAGS)
