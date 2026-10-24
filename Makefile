MAKEFLAGS += --warn-undefined-variables

CC	:= gcc
SRC_DIR	:= src
TEST_DIR	:= tests
BUILD	:= build

SRCS	:= $(wildcard $(SRC_DIR)/*.c)
TESTS	:= $(wildcard $(TEST_DIR)/*.c)
TEST_NAMES	:= $(basename $(notdir $(TESTS)))

COMMON	:= -fopenmp -Wall -Wextra -I$(SRC_DIR) -MMD -MP
DBG_FLAGS	:= -g -O0 -fsanitize=address,undefined
OPT_FLAGS	:= -O3 -march=native
LDLIBS	:= -lm

DBG	:= $(BUILD)/dbg
OPT	:= $(BUILD)/opt

DBG_SRC_OBJ	:= $(SRCS:$(SRC_DIR)/%.c=$(DBG)/src/%.o)
OPT_SRC_OBJ	:= $(SRCS:$(SRC_DIR)/%.c=$(OPT)/src/%.o)
DBG_TEST_OBJ	:= $(TESTS:$(TEST_DIR)/%.c=$(DBG)/tests/%.o)
OPT_TEST_OBJ	:= $(TESTS:$(TEST_DIR)/%.c=$(OPT)/tests/%.o)
DBG_BINS	:= $(TEST_NAMES:%=$(DBG)/bin/%)
OPT_BINS	:= $(TEST_NAMES:%=$(OPT)/bin/%)

.DEFAULT_GOAL	:= all
all: $(DBG_BINS)
opt: $(OPT_BINS)

# Debug build

$(DBG)/src/%.o: $(SRC_DIR)/%.c
	@mkdir -p $(dir $@)
	$(CC) $(COMMON) $(DBG_FLAGS) -c $< -o $@

$(DBG)/tests/%.o: $(TEST_DIR)/%.c
	@mkdir -p $(dir $@)
	$(CC) $(COMMON) $(DBG_FLAGS) -c $< -o $@

$(DBG)/bin/%: $(DBG)/tests/%.o $(DBG_SRC_OBJ)
	@mkdir -p $(dir $@)
	$(CC) $(DBG_FLAGS) -fopenmp -o $@ $^ $(LDLIBS)

# Optimized build

$(OPT)/src/%.o: $(SRC_DIR)/%.c
	@mkdir -p $(dir $@)
	$(CC) $(COMMON) $(OPT_FLAGS) -c $< -o $@

$(OPT)/tests/%.o: $(TEST_DIR)/%.c
	@mkdir -p $(dir $@)
	$(CC) $(COMMON) $(OPT_FLAGS) -c $< -o $@

$(OPT)/bin/%: $(OPT)/tests/%.o $(OPT_SRC_OBJ)
	@mkdir -p $(dir $@)
	$(CC) $(OPT_FLAGS) -fopenmp -o $@ $^ $(LDLIBS)

run-%: $(DBG)/bin/%
	./$<

bench-%: $(OPT)/bin/%
	./$<

clean:
	rm -rf $(BUILD)

.SECONDARY:
.PHONY: test bench run clean

-include $(DBG_SRC_OBJ:.o=.d) $(DBG_TEST_OBJ:.o=.d) $(OPT_SRC_OBJ:.o=.d) $(OPT_TEST_OBJ:.o=.d)
