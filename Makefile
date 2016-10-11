CXXFLAGS = -std=c++11 -Wall -Werror -O3 -g -Wno-unused-function
CXX = g++

ECFP = ecfp.o molecule.o stringhash.o 
SCORE = score.o molecule.o stringhash.o

default: ecfp score

ecfp: $(ECFP)
	$(CXX) $(CXXFLAGS) -o $@ $^

score: $(SCORE)
	$(CXX) $(CXXFLAGS) -o $@ $^

ecfp.o: stringhash.cc molecule.cc

score.o: stringhash.cc molecule.cc

clean:
	rm -f ecfp score *.o *~ output/*
