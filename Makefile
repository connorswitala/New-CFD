# ==== Compiler and Flags ====
CXX       := mpic++
CXXFLAGS  := -std=c++17 -O2 -Wall -I./linalglib

# ==== Directories ====
LIN_DIR       := linalglib
TEST_DIR      := testing
BUILD_DIR     := build
PROGRAMS_DIR  := programs
SOD_DIR       := explicit_sod

# ==== Files ====
LIN_SRC    := $(LIN_DIR)/linalg.cpp
LIN_OBJ    := $(BUILD_DIR)/linalg.o

SOD_SRC    := $(SOD_DIR)/main.cpp
SOD_OBJ    := $(BUILD_DIR)/sod.o
SOD_TARGET := $(PROGRAMS_DIR)/sod

TEST_SRC   := $(TEST_DIR)/main.cpp
TEST_OBJ   := $(BUILD_DIR)/test.o
TEST_TARGET := $(PROGRAMS_DIR)/test

# ==== Default target ====
all: $(SOD_TARGET) $(TEST_TARGET)

# ==== Ensure build and programs dirs exist ====
$(BUILD_DIR):
	mkdir -p $(BUILD_DIR)

$(PROGRAMS_DIR):
	mkdir -p $(PROGRAMS_DIR)

# ==== Object file compilation ====
$(LIN_OBJ): $(LIN_SRC) | $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(SOD_OBJ): $(SOD_SRC) | $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(TEST_OBJ): $(TEST_SRC) | $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

# ==== Link targets ====
$(SOD_TARGET): $(SOD_OBJ) $(LIN_OBJ) | $(PROGRAMS_DIR)
	$(CXX) $(CXXFLAGS) $^ -o $@

$(TEST_TARGET): $(TEST_OBJ) $(LIN_OBJ) | $(PROGRAMS_DIR)
	$(CXX) $(CXXFLAGS) $^ -o $@

# ==== Clean ====
clean:
	rm -rf $(BUILD_DIR) $(PROGRAMS_DIR)

.PHONY: all clean
