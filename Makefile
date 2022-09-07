CC = g++
MPICC = mpicc

MST_CFLAGS = -mcx16 -std=c++17 -DNDEBUG -DPARLAY_OPENMP -I"parallelKruskal"
MST_LFLAGS = -DPARLAY_OPENMP -ldl

INCLUDE = -IfeGRASS
C_FLAGS= -O3 -fopenmp $(INCLUDE) $(MST_CFLAGS)
L_FLAGS = -fopenmp $(MST_LFLAGS)

SRC_DIR = feGRASS
BUILD_DIR = build/bin
MST_SRC = parallelKruskal

SRCS = $(wildcard $(SRC_DIR)/*.cpp)
OBJ = $(SRCS:.cpp=.o) $(MST_SRC)/MST.o
DEBUG_OBJ = $(SRCS:.cpp=_debug.o) $(MST_SRC)/MST.o
TIME_OBJ = $(SRCS:.cpp=_time.o) $(MST_SRC)/MST.o
FILENAME = $(SRCS:.cpp=)

.DEFAULT_GOAL := main
main: ${OBJ}
	$(CC) $^ $(L_FLAGS) -o $(BUILD_DIR)/main

$(MST_SRC)/MST.o:
	make -C parallelKruskal

debugPrint: ${DEBUG_OBJ}
	$(CC) $^ $(L_FLAGS) -o $(BUILD_DIR)/main

timePrint: ${TIME_OBJ}
	$(CC) $^ $(L_FLAGS) -o $(BUILD_DIR)/main

mpi: ${SRCS}
	$(MPICC) $^ $(C_FLAGS) -o $(BUILD_DIR)/main

debugMpi: ${SRCS}
	$(MPICC) -DDEBUG -DTIME $^ $(C_FLAGS) -o $(BUILD_DIR)/main

timeMpi: ${SRCS}
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
	make clean -C parallelKruskal