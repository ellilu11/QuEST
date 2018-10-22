CC = g++
FLAGS = -fopenmp

CPPS = node.cpp nodefunc.cpp mathfunc.cpp out.cpp main.cpp

EXENAME = fmm2d

$(EXENAME):
		$(CC) -o $(EXENAME) $(CPPS) $(FLAGS)
