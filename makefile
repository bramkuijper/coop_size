EXE=coop_size.eexe
CPP=individual.cpp simulation.cpp main.cpp 
HPP=individual.hpp simulation.hpp parameters.hpp 
CXX=g++
CXXFLAGS=-Wall -O3 -ggdb -std=c++17

$(EXE) : $(CPP) $(HPP)
	$(CXX) $(CXXFLAGS) -o $(EXE) $(CPP)

clean :
	rm -rf $(EXE)
