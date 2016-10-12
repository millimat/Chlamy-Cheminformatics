CXXFLAGS = -std=c++11 -Wall -Werror -O3 -g
CXX = g++

ECFP = ecfp.o molecule.o stringhash.o 
ECFP-TEST = tests/ecfp-test.o molecule.o stringhash.o
STRINGHASH-TEST = tests/stringhash-test.o stringhash.o

default: ecfp tests

tests: tests/ecfp-test tests/stringhash-test

ecfp: $(ECFP)
	$(CXX) $(CXXFLAGS) -o $@ $^

tests/ecfp-test: $(ECFP-TEST)
	$(CXX) $(CXXFLAGS) -o $@ $^

tests/stringhash-test: $(STRINGHASH-TEST)
	$(CXX) $(CXXFLAGS) -o $@ $^

ecfp.o: stringhash.cc molecule.cc
tests/ecfp-test.o: stringhash.cc molecule.cc
tests/stringhash-test.o: stringhash.cc

clean:
	rm -f ecfp *.o *~ tests/ecfp-test tests/stringhash-test tests/*.o tests/*~ output/*
