CC = g++
MPICC = mpicc

debugFlag= -g
LIB = -lm 
C_FLAGS= -O3 $(LIB) ${debugFlag}
# C_FLAGS= -O3 -march=znver1 -mavx2 -fopenmp $(LIB) ${debugFlag}

INCLUDEPATH = feGRASS
SRC_DIR = feGRASS
BUILD_DIR = build/bin

SRCS = $(wildcard $(SRC_DIR)/*.cpp)
OBJ = $(SRCS:.cpp=.o)
FILENAME = $(SRCS:.cpp=)

.DEFAULT_GOAL := all
all : ${SRCS}
	echo "compiling $(SRC_DIR) ${FILENAME}"
	$(CC) $^ $(C_FLAGS) -o $(BUILD_DIR)/main

mpi: ${SRCS}
	echo "compiling $(SRC_DIR) ${FILENAME}"
	$(MPICC) $^ $(C_FLAGS) -o $(BUILD_DIR)/main

debugMpi: ${SRCS}
	echo "compiling $(SRC_DIR) ${FILENAME}"
	$(MPICC) -DDEBUG -DTIME $^ $(C_FLAGS) -o $(BUILD_DIR)/main

timeMpi: ${SRCS}
	echo "compiling $(SRC_DIR) ${FILENAME}"
	$(MPICC) -DTIME $^ $(C_FLAGS) -o $(BUILD_DIR)/main




checkdirs: $(BUILD_DIR)
$(BUILD_DIR):
	@mkdir -p $@

.PHONY: clean
clean:
	rm -rf $(BUILD_DIR)/main *.o