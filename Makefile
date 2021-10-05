SRC_DIR := src
HEADER_DIR := include
OBJ_DIR := obj
BIN_DIR := bin

SOURCES := $(filter-out $(SRC_DIR)/Tests.cpp,$(wildcard $(SRC_DIR)/*.cpp))

OBJECTS_DBG := $(patsubst $(SRC_DIR)/%.cpp,$(OBJ_DIR)/%_dbg.o,$(SOURCES))
OBJECTS_REL := $(patsubst $(SRC_DIR)/%.cpp,$(OBJ_DIR)/%_rel.o,$(SOURCES))
OBJECTS_TEST := $(patsubst $(SRC_DIR)/%.cpp,$(OBJ_DIR)/%_test.o,$(SOURCES))

OBJECTS_TEST := $(subst Main,Tests,$(OBJECTS_TEST))

DEPENDS_DBG := $(OBJECTS_DBG:.o=.d)
DEPENDS_REL := $(OBJECTS_REL:.o=.d)
DEPENDS_TEST := $(OBJECTS_TEST:.o=.d)

TARGET_DBG := $(BIN_DIR)/main_dbg
TARGET_REL := $(BIN_DIR)/main
TARGET_TEST := $(BIN_DIR)/tests
TARGETS := $(TARGET_DBG) $(TARGET_REL) $(TARGET_TEST)

CXX_EXPERIMENTAL_FILESYSTEM := $(shell expr `g++ -dumpversion | cut -f1 -d.` \<= 7)

CXX := g++
CPPFLAGS := -I$(HEADER_DIR) `pkg-config --cflags yaml-cpp` -Wall -Werror -Wextra -Wundef -Wfloat-equal -Wpointer-arith -Wcast-align -Wwrite-strings -Wcast-qual -Wswitch-default -Wswitch-enum -Wunreachable-code -Wimplicit-fallthrough -Wstringop-overflow=4 -Wconversion -MMD
CXXFLAGS := -std=c++17

ifeq "$(CXX_EXPERIMENTAL_FILESYSTEM)" "1"
	LDFLAGS := `pkg-config --libs yaml-cpp` -lstdc++fs
else
	LDFLAGS := `pkg-config --libs yaml-cpp`
endif

CPPFLAGS_DBG := $(CPPFLAGS) -g -Og
CPPFLAGS_REL := $(CPPFLAGS) -flto -Ofast -march=native
CPPFLAGS_TEST := $(CPPFLAGS) -fsanitize=address -fsanitize=undefined -g -Og

CXXFLAGS_DBG := $(CXXFLAGS) -D_GLIBCXX_ASSERTIONS -DTRACE_MATRIX_OPTIMIZATION
CXXFLAGS_REL := $(CXXFLAGS) -D_GLIBCXX_ASSERTIONS
CXXFLAGS_TEST := $(CXXFLAGS) -D_GLIBCXX_DEBUG

LDFLAGS_DBG := $(LDFLAGS) -g -Og
LDFLAGS_REL := $(LDFLAGS) -flto -Ofast -mtune=native
LDFLAGS_TEST := $(LDFLAGS) -fsanitize=address -fsanitize=undefined -g -Og

.PHONY: rel all dbg test clean veryclean

rel: $(TARGET_REL)
all: $(TARGETS)
dbg: $(TARGET_DBG)
test: $(TARGET_TEST)
test_internal: $(TARGET_TEST)
	$(TARGET_TEST)

clean:
	-@rm -f $(OBJECTS_DBG) $(OBJECTS_REL) $(OBJECTS_TEST) $(TARGETS)
veryclean: clean
	-@rm -fr $(BIN_DIR) $(OBJ_DIR)

$(OBJ_DIR)/%_dbg.o: $(SRC_DIR)/%.cpp
	@mkdir -p $(OBJ_DIR)
	$(CXX) -c $(CPPFLAGS_DBG) $(CXXFLAGS_DBG) -o $@ $<

$(OBJ_DIR)/%_rel.o: $(SRC_DIR)/%.cpp $(PERF_FILE)
	@mkdir -p $(OBJ_DIR)
	$(CXX) -c $(CPPFLAGS_REL) $(CXXFLAGS_REL) -o $@ $<

$(OBJ_DIR)/%_test.o: $(SRC_DIR)/%.cpp
	@mkdir -p $(OBJ_DIR)
	$(CXX) -c $(CPPFLAGS_TEST) $(CXXFLAGS_TEST) -o $@ $<

$(TARGET_DBG): $(OBJECTS_DBG)
	@mkdir -p $(BIN_DIR)
	$(CXX) $^ $(LDFLAGS_DBG) -o $@

$(TARGET_REL): $(OBJECTS_REL)
	@mkdir -p $(BIN_DIR)
	$(CXX) $^ $(LDFLAGS_REL) -o $@

$(TARGET_TEST): $(OBJECTS_TEST)
	@mkdir -p $(BIN_DIR)
	$(CXX) $^ $(LDFLAGS_TEST) -o $@


-include $(DEPENDS_DBG)
-include $(DEPENDS_REL)
-include $(DEPENDS_TEST)
