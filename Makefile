CC = g++
MPICC = mpicc

C_FLAGS= -O3 -fopenmp -g
LIB = -lgomp
# C_FLAGS= -fopenmp $(LIB) $(debugFlag)
# C_FLAGS= -O3 -march=znver1 -mavx2 -fopenmp $(LIB) $(debugFlag)

INCLUDEPATH = feGRASS
SRC_DIR = feGRASS
BUILD_DIR = build/bin
OBJ_DIR = build/obj

SRCS = $(wildcard $(SRC_DIR)/*.cpp)
TMP = $(patsubst %.cpp,${OBJ_DIR}/%.cpp,$(notdir ${SRCS}))
HEADERS = $(wildcard $(SRC_DIR)/*.h)
OBJ = $(TMP:.cpp=.o)
DEBUG_OBJ = $(TMP:.cpp=_debug.o)
TIME_OBJ = $(TMP:.cpp=_time.o)
FILENAME = $(TMP:.cpp=)

$(info    SRCS is: $(SRCS))
$(info    OBJ is: $(OBJ))
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

${OBJ_DIR}/%.o: ${SRC_DIR}/%.cpp $(HEADERS)
	$(CC) $(C_FLAGS) -c $< -o $@ 
${OBJ_DIR}/%_debug.o: ${SRC_DIR}/%.cpp $(HEADERS)
	$(CC) -DDEBUG -DTIME $(C_FLAGS) -c $< -o $@ 
${OBJ_DIR}/%_time.o: ${SRC_DIR}/%.cpp $(HEADERS)
	$(CC) -DTIME $(C_FLAGS) -c $< -o $@ 


checkdirs: $(BUILD_DIR)
$(BUILD_DIR):
	@mkdir -p $(BUILD_DIR)
	@mkdir -p $(OBJ_DIR)


.PHONY: clean
clean:
	rm -rf $(BUILD_DIR)/main $(OBJ_DIR)/*.o