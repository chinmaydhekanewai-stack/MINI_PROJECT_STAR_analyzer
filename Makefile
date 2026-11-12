
CC       ?= gcc
CXX      ?= g++

CFLAGS   ?= -std=c99 -Wall -Wextra -Wpedantic -O2
CXXFLAGS ?= -std=c++14 -Wall -Wextra -O2
CPPFLAGS += -I. -MMD -MP
LDLIBS   += -lm
GTEST_INC  ?=
GTEST_LIBS ?= -lgtest -lgtest_main -pthread

BUILD := build
BIN   := bin

LIB_SRC  := data_access.c parameters.c sky_position.c zones.c
APP_SRC  := main.c
TEST_SRC := tests/test_stellarcalc.cpp

LIB_OBJ  := $(LIB_SRC:%.c=$(BUILD)/%.o)
APP_OBJ  := $(APP_SRC:%.c=$(BUILD)/%.o)
TEST_OBJ := $(TEST_SRC:tests/%.cpp=$(BUILD)/%.o)

.PHONY: all objs tests test clean

all: $(BIN)/stellarcalc

# compile only: every .o, nothing linked
objs: $(LIB_OBJ) $(APP_OBJ)

$(BIN)/stellarcalc: $(LIB_OBJ) $(APP_OBJ) | $(BIN)
	$(CC) $(CFLAGS) $^ -o $@ $(LDLIBS)

tests: $(BIN)/run_tests

$(BIN)/run_tests: $(LIB_OBJ) $(TEST_OBJ) | $(BIN)
	$(CXX) $(CXXFLAGS) $^ -o $@ $(GTEST_LIBS) $(LDLIBS)

test: tests
	./$(BIN)/run_tests

$(BUILD)/%.o: %.c | $(BUILD)
	$(CC) $(CPPFLAGS) $(CFLAGS) -c $< -o $@

$(BUILD)/%.o: tests/%.cpp | $(BUILD)
	$(CXX) $(CPPFLAGS) $(GTEST_INC) $(CXXFLAGS) -c $< -o $@

$(BUILD) $(BIN):
	mkdir -p $@

clean:
	rm -rf $(BUILD) $(BIN)

-include $(wildcard $(BUILD)/*.d)
