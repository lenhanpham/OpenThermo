# Makefile for OpenThermo
# Auto-detect compiler version

# Directory structure
SRC_DIR = src
BUILD_DIR = build

# Compiler settings
# Auto-detect compiler. Prefers Intel compilers (icpx, icpc, icc) over GCC (g++) and Clang.
COMPILER_LIST := icpx icpc icc clang++ g++
CXX := $(firstword $(foreach c,$(COMPILER_LIST),$(if $(shell command -v $(c)),$(c))))

# Fallback to a default if no compiler is found in PATH and print a warning.
ifeq ($(CXX),)
    $(warning "No supported compiler (icpx, icpc, icc, g++) found in PATH. Defaulting to g++.")
    CXX = g++
endif
$(info Using compiler: $(CXX))

# Base flags
CXXFLAGS = -std=c++17 -Wall -Wextra -O2
LDFLAGS =

# Add Intel-specific flags if using Intel compiler
ifneq ($(filter icpx icpc icc,$(CXX)),)
    CXXFLAGS += -fp-model=precise
    # Intel compilers may require TBB for some features, but OpenThermo doesn't use threading
    # LDFLAGS += -ltbb  # Uncomment if needed
    $(info Intel compiler detected. TBB linking commented out as not required for OpenThermo.)
endif

DEBUGFLAGS = -g -DDEBUG_BUILD -fsanitize=address -fno-omit-frame-pointer

# Platform detection
UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Linux)
    LDFLAGS += -lrt -lstdc++fs
endif
ifeq ($(UNAME_S),Darwin)
    # macOS specific flags if needed
endif

# Source files
SOURCES = $(SRC_DIR)/main.cpp \
          $(SRC_DIR)/atommass.cpp \
          $(SRC_DIR)/calc.cpp \
          $(SRC_DIR)/loadfile.cpp \
          $(SRC_DIR)/symmetry.cpp \
          $(SRC_DIR)/util.cpp \
          $(SRC_DIR)/help_utils.cpp

HEADERS = $(SRC_DIR)/atommass.h \
           $(SRC_DIR)/calc.h \
           $(SRC_DIR)/chemsys.h \
           $(SRC_DIR)/loadfile.h \
           $(SRC_DIR)/symmetry.h \
           $(SRC_DIR)/util.h \
           $(SRC_DIR)/help_utils.h

OBJECTS = $(SOURCES:$(SRC_DIR)/%.cpp=$(BUILD_DIR)/%.o)
TARGET = $(BUILD_DIR)/OpenThermo

# Ensure build directory exists
$(shell mkdir -p $(BUILD_DIR))

# Default target
all: $(TARGET)

# Build the main executable
$(TARGET): $(OBJECTS)
	$(CXX) $(OBJECTS) -o $@ $(LDFLAGS)

# Compile source files to object files
$(BUILD_DIR)/%.o: $(SRC_DIR)/%.cpp $(HEADERS)
	@mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Debug build with additional safety checks
debug: CXXFLAGS += $(DEBUGFLAGS)
debug: LDFLAGS += -fsanitize=address
debug: clean $(TARGET)

# Release build with optimizations
release: CXXFLAGS += -O3 -DNDEBUG -march=native
release: clean $(TARGET)

# Clean build artifacts
clean:
	rm -rf $(BUILD_DIR) $(TARGET)

# Help target
help:
	@echo "OpenThermo C++ Makefile - Available targets:"
	@echo ""
	@echo "  all          - Build the program (default)"
	@echo "  debug        - Build with debug symbols and AddressSanitizer"
	@echo "  release      - Build optimized release version"
	@echo "  clean        - Remove build artifacts"
	@echo "  help         - Show this help message"
	@echo ""
	@echo "Compiler Support:"
	@echo "  Auto-detects: icpx, icpc, icc, g++ (in order of preference)"
	@echo "  Force GCC:    make CXX=g++"
	@echo "  Intel notes:  TBB linking is commented out as OpenThermo doesn't require it"
	@echo ""
	@echo "Usage examples:"
	@echo "  make                    # Build with default settings"
	@echo "  make CXX=g++            # Force GCC compilation"
	@echo "  make debug              # Build debug version"
	@echo "  make release            # Build optimized version"
	@echo "  make clean              # Clean build"

# Declare phony targets
.PHONY: all debug release clean help