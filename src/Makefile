CXX ?= g++
CXXFLAGS ?= -O3 -DNDEBUG -std=c++20 -Wall -Wextra -Wpedantic
OPENMP ?= -fopenmp
SOURCES = main.cpp compute.cpp timestep.cpp read.cpp print.cpp checkpoint.cpp dipole.cpp
OBJECTS = $(SOURCES:.cpp=.o)
.PHONY: all clean test benchmark
all: point_vortex
point_vortex: $(OBJECTS)
	$(CXX) $(CXXFLAGS) $(OPENMP) -o $@ $(OBJECTS)
%.o: %.cpp
	$(CXX) $(CXXFLAGS) $(OPENMP) -MMD -MP -c $< -o $@
tests: tests.cpp compute.cpp timestep.cpp read.cpp checkpoint.cpp dipole.cpp
	$(CXX) $(CXXFLAGS) $(OPENMP) -o $@ $^
test: tests
	./tests
benchmark: benchmark.cpp compute.cpp
	$(CXX) $(CXXFLAGS) $(OPENMP) -o $@ $^
clean:
	$(RM) $(OBJECTS) $(OBJECTS:.o=.d) point_vortex tests benchmark
-include $(OBJECTS:.o=.d)
