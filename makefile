# Problem to compile
PROB = torus

SYSTEM_LIBDIR = /lib64

CC=gcc
# Example CFLAGS for going fast with GCC
CFLAGS = -std=gnu99 -O3 -march=native -mtune=native -flto -fopenmp -funroll-loops
MATH_LIB = -lm

# Name of the executable
EXE = harm

MAKEFILE_PATH := $(dir $(abspath $(firstword $(MAKEFILE_LIST))))

LINK = $(CC)
LDFLAGS = $(CFLAGS)

CORE_DIR := $(MAKEFILE_PATH)/core/
PROB_DIR := $(MAKEFILE_PATH)/prob/$(PROB)/
VPATH = $(CORE_DIR):$(PROB_DIR)
ARC_DIR := $(CURDIR)/build_archive/

SRC := $(wildcard $(CORE_DIR)/*c) $(wildcard $(PROB_DIR)/*.c)
HEAD := $(wildcard $(CORE_DIR)/*h) $(wildcard $(PROB_DIR)/*.h)

HEAD_ARC := $(addprefix $(ARC_DIR), $(notdir $(HEAD)))
OBJ := $(addprefix $(ARC_DIR), $(notdir $(SRC:%.c=%.o)))

INC = -I$(ARC_DIR)
LIBDIR =
LIB = $(MATH_LIB)
ifneq ($(strip $(SYSTEM_LIBDIR)),)
	LIBDIR += -L$(SYSTEM_LIBDIR)
endif

default: build

build: $(EXE)
	@echo -e "Completed build of prob: $(PROB)"
	@echo -e "CFLAGS: $(CFLAGS)"

clean:
	@echo "Cleaning build files..."
	@rm -f $(EXE) $(OBJ)

distclean: clean
	@echo "Cleaning config files..."
	@rm -rf build_archive

$(EXE): $(ARC_DIR)$(EXE)
	@cp $(ARC_DIR)$(EXE) .

$(ARC_DIR)$(EXE): $(OBJ)
	@$(LINK) $(LDFLAGS) $(OBJ) $(LIBDIR) $(LIB) -o $(ARC_DIR)$(EXE)
	@rm $(OBJ)
	

$(OBJ): $(ARC_DIR)%.o: $(ARC_DIR)%.c $(HEAD_ARC)
	@echo -e "\tCompiling $(notdir $<)"
	@$(CC) $(CFLAGS) $(INC) -c $< -o $@

$(ARC_DIR)%: % | $(ARC_DIR)
	@cp $< $(ARC_DIR)

#$(ARC_DIR)%.c: $(CORE_DIR)%.c $(PROB_DIR)%.c $(HEAD) | $(ARC_DIR)
#	@cp $^ $(ARC_DIR)

$(ARC_DIR):
	@mkdir $(ARC_DIR)
