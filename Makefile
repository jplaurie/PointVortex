CXX ?= g++
CPPFLAGS ?= -Isrc -Iinitial_conditions
CXXFLAGS ?= -O3 -DNDEBUG -std=c++20 -Wall -Wextra -Wpedantic
OPENMP ?= -fopenmp
MPICXX ?= mpicxx
NVCC ?= nvcc
CUDA_ARCH ?= sm_120
BUILD_DIR ?= build/make
OBJ_DIR := $(BUILD_DIR)/obj

LIB_SOURCES = src/compute.cpp src/timestep.cpp src/read.cpp src/checkpoint.cpp src/dipole.cpp \
	src/backend_common.cpp initial_conditions/initial_condition.cpp
APP_SOURCES = src/main.cpp src/print.cpp
LIB_OBJECTS = $(patsubst %.cpp,$(OBJ_DIR)/%.o,$(LIB_SOURCES))
APP_OBJECTS = $(patsubst %.cpp,$(OBJ_DIR)/%.o,$(APP_SOURCES))
CPU_OBJECT = $(OBJ_DIR)/src/backend_cpu.o
CUDA_OBJECT = $(OBJ_DIR)/src/backend_cuda.o

.PHONY: all cpu mpi cuda initial test benchmark clean
all: cpu initial
cpu: $(BUILD_DIR)/point_vortex $(BUILD_DIR)/point_vortex_cpu
mpi: $(BUILD_DIR)/point_vortex_mpi
cuda: $(BUILD_DIR)/point_vortex_cuda
initial: $(BUILD_DIR)/point_vortex_initial
benchmark: $(BUILD_DIR)/point_vortex_benchmark

$(BUILD_DIR) $(OBJ_DIR):
	mkdir -p $@

$(OBJ_DIR)/%.o: %.cpp | $(OBJ_DIR)
	mkdir -p $(dir $@)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(OPENMP) -MMD -MP -c $< -o $@

$(CUDA_OBJECT): src/backend_cuda.cu | $(OBJ_DIR)
	mkdir -p $(dir $@)
	$(NVCC) $(CPPFLAGS) -O3 -std=c++20 -arch=$(CUDA_ARCH) -Xcompiler=$(OPENMP) -c $< -o $@

$(BUILD_DIR)/point_vortex $(BUILD_DIR)/point_vortex_cpu: $(APP_OBJECTS) $(CPU_OBJECT) $(LIB_OBJECTS) | $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) $(OPENMP) -o $@ $^

$(BUILD_DIR)/point_vortex_mpi: $(APP_SOURCES) src/backend_mpi.cpp $(LIB_SOURCES) | $(BUILD_DIR)
	$(MPICXX) $(CPPFLAGS) $(CXXFLAGS) $(OPENMP) -o $@ $^

$(BUILD_DIR)/point_vortex_cuda: $(APP_OBJECTS) $(CUDA_OBJECT) $(LIB_OBJECTS) | $(BUILD_DIR)
	$(NVCC) -arch=$(CUDA_ARCH) -Xcompiler=$(OPENMP) -o $@ $^

$(BUILD_DIR)/point_vortex_initial: initial_conditions/generate_initial.cpp \
		initial_conditions/initial_condition.cpp src/read.cpp | $(BUILD_DIR)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -o $@ $^

$(BUILD_DIR)/point_vortex_tests: tests/tests.cpp $(LIB_SOURCES) | $(BUILD_DIR)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(OPENMP) -o $@ $^

test: $(BUILD_DIR)/point_vortex_tests
	$(BUILD_DIR)/point_vortex_tests

$(BUILD_DIR)/point_vortex_benchmark: src/benchmark.cpp src/compute.cpp | $(BUILD_DIR)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(OPENMP) -o $@ $^

clean:
	$(RM) -r $(BUILD_DIR)

-include $(LIB_OBJECTS:.o=.d) $(APP_OBJECTS:.o=.d) $(CPU_OBJECT:.o=.d)
