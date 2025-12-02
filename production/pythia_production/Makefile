# Makefile for HNL Production Simulation
#
# Requires PYTHIA8 to be installed. Set PYTHIA8 environment variable
# or modify the path below.

# Pythia8 location (modify if needed)
PYTHIA8 ?= /usr/local/pythia8
PYTHIA8_INCLUDE = $(PYTHIA8)/include
PYTHIA8_LIB = $(PYTHIA8)/lib

# Compiler settings
CXX = g++
CXXFLAGS = -O2 -std=c++17 -Wall -I$(PYTHIA8_INCLUDE)
LDFLAGS = -L$(PYTHIA8_LIB) -Wl,-rpath,$(PYTHIA8_LIB) -lpythia8 -ldl

# Targets
TARGET = main_hnl_production
SOURCES = main_hnl_production.cc

.PHONY: all clean help test

all: $(TARGET)

$(TARGET): $(SOURCES)
	$(CXX) $(CXXFLAGS) -o $@ $< $(LDFLAGS)

clean:
	rm -f $(TARGET) *.csv

help:
	@echo "HNL Production Simulation Makefile"
	@echo ""
	@echo "Usage:"
	@echo "  make          - Build the simulation"
	@echo "  make clean    - Remove built files"
	@echo "  make test     - Run a quick test"
	@echo ""
	@echo "Environment variables:"
	@echo "  PYTHIA8       - Path to Pythia8 installation"
	@echo "                  (default: /usr/local/pythia8)"
	@echo ""
	@echo "Example:"
	@echo "  PYTHIA8=/opt/pythia8312 make"

test: $(TARGET)
	@echo "Running quick test (1000 events)..."
	./$(TARGET) 2.0 muon 1000
	@echo ""
	@echo "Test complete. Check output CSV file."
