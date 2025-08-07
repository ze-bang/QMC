# Alternative build system using Make
# Use this if CMake is not available

CXX = g++
CXXFLAGS = -std=c++17 -O3 -march=native -Wall -Wextra -fopenmp
INCLUDES = -Iinclude
SRCDIR = src
OBJDIR = obj
SOURCES = $(wildcard $(SRCDIR)/*.cpp)
OBJECTS = $(SOURCES:$(SRCDIR)/%.cpp=$(OBJDIR)/%.o)
TARGET = sse_qmc

.PHONY: all clean debug

all: $(TARGET)

$(TARGET): $(OBJECTS)
	$(CXX) $(OBJECTS) -o $@ $(CXXFLAGS)

$(OBJDIR)/%.o: $(SRCDIR)/%.cpp | $(OBJDIR)
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

$(OBJDIR):
	mkdir -p $(OBJDIR)

debug: CXXFLAGS = -std=c++17 -g -O0 -Wall -Wextra -fopenmp
debug: $(TARGET)

clean:
	rm -rf $(OBJDIR) $(TARGET)

# Run examples
run-square-heisenberg:
	./$(TARGET) --lattice square --model heisenberg --Lx 8 --Ly 8 --J 1.0 --beta 10.0 --therm 10000 --meas 50000

run-triangular-heisenberg:
	./$(TARGET) --lattice triangular --model heisenberg --Lx 6 --Ly 6 --J 1.0 --beta 10.0 --therm 10000 --meas 50000

run-chain-ising:
	./$(TARGET) --lattice chain --model ising --Lx 32 --J 1.0 --beta 5.0 --therm 10000 --meas 50000

run-xxz-model:
	./$(TARGET) --lattice square --model xxz --Lx 8 --Ly 8 --J 1.0 --Delta 0.5 --beta 10.0 --therm 10000 --meas 50000
