# Name of the compiled executable
EXEC = engulfment-integrator
BUILD_DIR = ./build
# Create build directory if it doesn't exist
$(shell mkdir -p $(BUILD_DIR))
#$(info $(shell mkdir -p $(BUILD_DIR)))

# Compiler
CC = gcc

SRC_DIR = ./src
SRC = $(wildcard $(SRC_DIR)/*.c)
CFLAGS = -O3 -Wall

OBJ = $(patsubst $(SRC_DIR)/%.c,$(BUILD_DIR)/%.o,$(SRC))
LIBS = -lm -lgsl -lgslcblas -lhdf5

$(BUILD_DIR)/$(EXEC): $(OBJ)
	$(CC) $(OBJ) $(LIBS) -o $(BUILD_DIR)/$(EXEC)

$(OBJ): $(BUILD_DIR)/%.o : $(SRC_DIR)/%.c
	$(CC) $(CFLAGS) $(INCL) -c $< -o $@ 

.PHONY : clean

clean:
	rm -rf $(BUILD_DIR)
	rm -rf docs/sphinx/build
	rm -rf docs/doxygen/build
