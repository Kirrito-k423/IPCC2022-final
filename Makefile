CC = g++
MPICC = mpicc

C_FLAGS= -O3 -fopenmp 
LIB = -lgomp
# C_FLAGS= -fopenmp $(LIB) ${debugFlag}
# C_FLAGS= -O3 -march=znver1 -mavx2 -fopenmp $(LIB) ${debugFlag}

INCLUDEPATH = feGRASS
SRC_DIR = feGRASS
BUILD_DIR = build/bin

SRCS = $(wildcard $(SRC_DIR)/*.cpp)
OBJ = $(SRCS:.cpp=.o)
DEBUG_OBJ = $(SRCS:.cpp=_debug.o)
TIME_OBJ = $(SRCS:.cpp=_time.o)
FILENAME = $(SRCS:.cpp=)

.DEFAULT_GOAL := all
all : ${OBJ}
	echo "compiling $(SRC_DIR) ${FILENAME}"
	$(CC) $^ $(LIB) -o $(BUILD_DIR)/main

debugPrint: ${DEBUG_OBJ}
	echo "compiling $(SRC_DIR) ${FILENAME}"
	$(CC) $^ $(LIB) -o $(BUILD_DIR)/main

timePrint: ${TIME_OBJ}
	echo "compiling $(SRC_DIR) ${FILENAME}"
	$(CC) $^ $(LIB) -o $(BUILD_DIR)/main

mpi: ${SRCS}
	echo "compiling $(SRC_DIR) ${FILENAME}"
	$(MPICC) $^ $(C_FLAGS) -o $(BUILD_DIR)/main

debugMpi: ${SRCS}
	echo "compiling $(SRC_DIR) ${FILENAME}"
	$(MPICC) -DDEBUG -DTIME $^ $(C_FLAGS) -o $(BUILD_DIR)/main

timeMpi: ${SRCS}
	echo "compiling $(SRC_DIR) ${FILENAME}"
	$(MPICC) -DTIME $^ $(C_FLAGS) -o $(BUILD_DIR)/main

%.o: %.cpp
	$(CC) $(C_FLAGS) -c $< -o $@ 
%_debug.o: %.cpp
	$(CC) -DDEBUG -DTIME $(C_FLAGS) -c $< -o $@ 
%_time.o: %.cpp
	$(CC) -DTIME $(C_FLAGS) -c $< -o $@ 


checkdirs: $(BUILD_DIR)
$(BUILD_DIR):
	@mkdir -p $@

.PHONY: clean
clean:
	rm -rf $(BUILD_DIR)/main ${SRC_DIR}/*.o