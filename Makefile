#
# Makefile for DRSA C code.
# @author Guilherme Oliveira Chagas (guilherme.o.chagas[a]gmail.com)
# @date 27/03/2018

CXX = gcc
CXXFLAGS = -g -Wall -Wextra -O3
LDFLAGS = -g
RM = rm -rf

LDLIBS = -lm

SRC_DIR = src
OBJ_DIR = obj

SRCS = $(SRC_DIR)/drsa_bmp.c $(SRC_DIR)/mmio.c
OBJS = $(addprefix $(OBJ_DIR)/,$(notdir $(subst .c,.o, $(SRCS) ) ) )
OUT_RELEASE = exec_drsa

all: out_release

out_release: $(OBJS)
	$(CXX) $(LDFLAGS) $(OBJS) -o $(OUT_RELEASE) $(LDLIBS)

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.c
	$(CXX) $(CXXFLAGS) -c $^ -o $@

clean:
	$(RM) $(OBJS) $(OUT_RELEASE)

.PHONY: clean