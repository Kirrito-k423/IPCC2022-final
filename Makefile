CC = g++
MPICC = mpicc

C_FLAGS= -O3 -fopenmp 
LIB = -lgomp
# C_FLAGS= -fopenmp $(LIB) $(debugFlag)
# C_FLAGS= -O3 -march=znver1 -mavx2 -fopenmp $(LIB) $(debugFlag)

INCLUDEPATH = feGRASS
SRC_DIR = feGRASS
BUILD_DIR = build/bin

SRCS = $(wildcard $(SRC_DIR)/*.cpp)
HEADERS = $(wildcard $(SRC_DIR)/*.h)
OBJ = $(SRCS:.cpp=.o)
DEBUG_OBJ = $(SRCS:.cpp=_debug.o)
TIME_OBJ = $(SRCS:.cpp=_time.o)
FILENAME = $(SRCS:.cpp=)

$(info    SRCS is: $(SRCS))
$(info    HEADERS is: $(HEADERS))

.DEFAULT_GOAL := main
main : $(OBJ)
	$(CC) $(OBJ) $(LIB) -o $(BUILD_DIR)/main

debugPrint: $(DEBUG_OBJ)
	$(CC) $(DEBUG_OBJ) $(LIB) -o $(BUILD_DIR)/main

timePrint: $(TIME_OBJ)
	$(CC) $(TIME_OBJ) $(LIB) -o $(BUILD_DIR)/main

mpi: $(SRCS)
	$(MPICC) $^ $(C_FLAGS) -o $(BUILD_DIR)/main

debugMpi: $(SRCS)
	$(MPICC) -DDEBUG -DTIME $^ $(C_FLAGS) -o $(BUILD_DIR)/main

timeMpi: $(SRCS)
	$(MPICC) -DTIME $^ $(C_FLAGS) -o $(BUILD_DIR)/main

%.o: %.cpp $(HEADERS)
	$(CC) $(C_FLAGS) -c $< -o $@ 
%_debug.o: %.cpp $(HEADERS)
	$(CC) -DDEBUG -DTIME $(C_FLAGS) -c $< -o $@ 
%_time.o: %.cpp $(HEADERS)
	$(CC) -DTIME $(C_FLAGS) -c $< -o $@ 


checkdirs: $(BUILD_DIR)
$(BUILD_DIR):
	@mkdir -p $@

.PHONY: clean
clean:
	rm -rf $(BUILD_DIR)/main $(SRC_DIR)/*.o