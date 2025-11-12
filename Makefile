# ==== Compiler and Flags ====
CXX       := mpic++
CXXFLAGS  := -std=c++17 -O2 -Wall -I./linalglib -I./solverlib

# ==== Directories ====
LIN_DIR       := linalglib
GIBBS_DIR	  := gibbslib
TEST_DIR      := testing
BUILD_DIR     := build
PROGRAMS_DIR  := programs
SOD_DIR       := explicit_sod
CFDSOLVER_DIR := cfdsolver
SOLVER_DIR 	  := solverlib

# ==== Files ====
LIN_SRC    := $(LIN_DIR)/linalg.cpp
LIN_OBJ    := $(BUILD_DIR)/linalg.o

GIBBS_SRC	:= $(GIBBS_DIR)/gibbs.cpp
GIBBS_OBJ	:= $(BUILD_DIR)/gibbs.o

SOLVER_SRC := $(SOLVER_DIR)/solver.cpp
SOLVER_OBJ := $(BUILD_DIR)/solver.o 

CFDSOLVER_SRC := $(CFDSOLVER_DIR)/main.cpp
CFDSOLVER_OBJ := $(BUILD_DIR)/cfdsolver.o
CFDSOLVER_TARGET := $(PROGRAMS_DIR)/cfdsolver

# Parallel Sod
SOD_SRC    := $(SOD_DIR)/parallel.cpp
SOD_OBJ    := $(BUILD_DIR)/sod.o
SOD_TARGET := $(PROGRAMS_DIR)/psod

# Serial Sod
SSOD_SRC    := $(SOD_DIR)/serial.cpp
SSOD_OBJ    := $(BUILD_DIR)/ssod.o
SSOD_TARGET := $(PROGRAMS_DIR)/ssod

TEST_SRC   := $(TEST_DIR)/main.cpp
TEST_OBJ   := $(BUILD_DIR)/test.o
TEST_TARGET := $(PROGRAMS_DIR)/test

# ==== Default target ====
all: $(SOD_TARGET) $(SSOD_TARGET) $(TEST_TARGET) $(CFDSOLVER_TARGET)

# ==== Ensure build and programs dirs exist ====
$(BUILD_DIR):
	mkdir -p $(BUILD_DIR)

$(PROGRAMS_DIR):
	mkdir -p $(PROGRAMS_DIR)

# ==== Object file compilation ====
$(LIN_OBJ): $(LIN_SRC) | $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(GIBBS_OBJ): $(GIBBS_SRC) | $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(SOD_OBJ): $(SOD_SRC) | $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(SSOD_OBJ): $(SSOD_SRC) | $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@	

$(TEST_OBJ): $(TEST_SRC) | $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(SOLVER_OBJ): $(SOLVER_SRC) | $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(CFDSOLVER_OBJ): $(CFDSOLVER_SRC) | $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

# ==== Link targets ====
$(CFDSOLVER_TARGET): $(CFDSOLVER_OBJ) $(LIN_OBJ) $(GIBBS_OBJ) $(SOLVER_OBJ) | $(PROGRAMS_DIR)
	$(CXX) $(CXXFLAGS) $^ -o $@

$(SOD_TARGET): $(SOD_OBJ) $(LIN_OBJ) | $(PROGRAMS_DIR)
	$(CXX) $(CXXFLAGS) $^ -o $@

$(SSOD_TARGET): $(SSOD_OBJ) $(LIN_OBJ) | $(PROGRAMS_DIR)
	$(CXX) $(CXXFLAGS) $^ -o $@

$(TEST_TARGET): $(TEST_OBJ) $(LIN_OBJ) $(GIBBS_OBJ) $(SOLVER_OBJ) | $(PROGRAMS_DIR)
	$(CXX) $(CXXFLAGS) $^ -o $@

# ==== Clean ====
clean:
	rm -rf $(BUILD_DIR) $(PROGRAMS_DIR)

.PHONY: all clean
