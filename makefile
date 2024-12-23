# Compiler and flags
CXX = g++
CXXFLAGS = -Wall -Wextra -std=c++20 -I. -I1D -I2D -I3D

# Executable names
EXEC_1D = 1D_program
EXEC_2D = 2D_program
EXEC_3D = 3D_program
EXEC_3D_2 = new_3D_program
EXEC_test = testvalues

# Source file directories
DIR_1D = 1D
DIR_2D = 2D
DIR_3D = 3D

# Source and object files
SRCS_1D = $(wildcard $(DIR_1D)/*.cpp)
OBJS_1D = $(SRCS_1D:.cpp=.o)
DEPS_1D = $(OBJS_1D:.o=.d)

SRCS_2D = $(wildcard $(DIR_2D)/*.cpp)
OBJS_2D = $(SRCS_2D:.cpp=.o)
DEPS_2D = $(OBJS_2D:.o=.d)

SRCS_3D = $(wildcard $(DIR_3D)/*.cpp)
OBJS_3D = $(SRCS_3D:.cpp=.o)
DEPS_3D = $(OBJS_3D:.o=.d)

SRCS_3D_2 = 3dCombined.cpp
OBJS_3D_2 = $(SRCS_3D_2:.cpp=.o)
DEPS_3D_2 = $(OBJS_3D_2:.o=.d)

SRCS_test = testvalues.cpp 3D/SolveRiemann.cpp
OBJS_test = $(SRCS_test:.cpp=.o)
DEPS_test = $(OBJS_test:.o=.d)

# Default target
all: $(EXEC_1D) $(EXEC_2D) $(EXEC_3D) $(EXEC_3D_2) $(EXEC_test)

# Rules to build individual executables
$(EXEC_1D): $(OBJS_1D)
	$(CXX) $(CXXFLAGS) -o $@ $^

$(EXEC_2D): $(OBJS_2D)
	$(CXX) $(CXXFLAGS) -o $@ $^

$(EXEC_3D): $(OBJS_3D)
	$(CXX) $(CXXFLAGS) -o $@ $^

$(EXEC_3D_2): $(OBJS_3D_2)
	$(CXX) $(CXXFLAGS) -o $@ $^

$(EXEC_test): $(OBJS_test)
	$(CXX) $(CXXFLAGS) -o $@ $^

# Rule to build object files and generate dependencies
%.o: %.cpp
	$(CXX) $(CXXFLAGS) -MMD -MP -c $< -o $@

# Include dependency files
-include $(DEPS_1D) $(DEPS_2D) $(DEPS_3D) $(DEPS_3D_2) $(DEPS_test)

# Targets to build only one executable
1D: $(EXEC_1D)

2D: $(EXEC_2D)

3D: $(EXEC_3D)

new3D: $(EXEC_3D_2)

test_values: $(EXEC_test)

# Clean up build files
clean:
	rm -f $(OBJS_1D) $(OBJS_2D) $(OBJS_3D) $(DEPS_1D) $(DEPS_2D) $(DEPS_3D) $(EXEC_1D) $(EXEC_2D) $(EXEC_3D) $(OBJS_3D_2) $(DEPS_3D_2) $(EXEC_3D_2) $(OBJS_test) $(DEPS_test) $(EXEC_test)

# Phony targets
.PHONY: all clean 1D 2D 3D new3D test_values
