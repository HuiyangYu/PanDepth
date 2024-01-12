BIN_DIR := $(shell pwd)/bin
LIB_DIR := $(shell pwd)/lib
SRC_DIR := $(shell pwd)/src
INCLUDE_DIR := $(shell pwd)/include

.PHONY: all clean

all: $(BIN_DIR)/pandepth

$(BIN_DIR)/pandepth: $(SRC_DIR)/PanDepth.cpp | $(BIN_DIR)
	g++ -L$(LIB_DIR) -Wl,-rpath=$(LIB_DIR) --std=c++11 -g -O3 $< -lhts -ldeflate -lz -pthread -I$(INCLUDE_DIR) -o $@

$(BIN_DIR):
	mkdir -p $(BIN_DIR)

clean:
	rm -rf $(BIN_DIR)


