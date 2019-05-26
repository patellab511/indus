### indus Makefile ###

#############
### Input ###
#############

XDR_DIR=$(HOME)/programs/xdrfile/1.1.4

# Compilers
CC=mpicc
CXX=mpic++

# TODO detect MPI automatically
is_mpi_enabled=1

# Output binary's name
PROJECT=indus

# Where to install 
INSTALL_DIR=$(HOME)/programs/indus


### !!! Shouldn't have to change anything below this line !!! ###


# Testing options
# - How to invoke mpirun for testing purposes
mpirun_cmd="mpirun --use-hwthread-cpus -np 4"
# - Whether to echo the differences between the tested files
echo_failed_diffs=1


#############
### Flags ###
#############

# GROMACS XDR Library
XDR_LIB     := $(XDR_DIR)/lib
XDR_INCLUDE := $(XDR_DIR)/include/xdrfile

# Compiler flags
# - Misc: -Wno-comment -Wno-sign-compare -DPLUMED_MODE
CXXFLAGS += -g -std=c++11 -DCPLUSPLUS -I$(XDR_INCLUDE) -Wall -DINDUS_STANDALONE_MODE
# - Optimizations
CXXFLAGS += -ffast-math -march=native -O3

# Check for whether MPI is enabled
ifeq ($(is_mpi_enabled),1)
CXXFLAGS += -DMPI_ENABLED
endif

# Linking
LDFLAGS += -g -L$(XDR_LIB)
LIBS    += -lm -lxdrfile

#version_flag = -mmacosx-version-min=10.6
#CXXFLAGS += $(version_flag)
#LDFLAGS  += $(version_flag)
#endif

#############################
### Directories and Files ###
#############################

START_DIR := $(PWD)

SRC_DIR         := src
BUILD_DIR       := build
BUILD_BIN_DIR   := $(BUILD_DIR)/bin
INSTALL_BIN_DIR := $(INSTALL_DIR)/bin

# Get source and target objects
SRC_FILES = $(shell find $(SRC_DIR) -name '*.cpp')
SRC_DIRS  = $(shell find $(SRC_DIR) -type d | sed 's/$(SRC_DIR)/./g' )
OBJECTS   = $(patsubst $(SRC_DIR)/%.cpp,$(BUILD_DIR)/%.o,$(SRC_FILES))


####################
### Main Program ###
####################

all : buildrepo $(PROJECT)

# Make directories 
buildrepo :
	@{ \
	for dir in $(BUILD_DIR) $(BUILD_BIN_DIR) $(INSTALL_DIR)/bin ; do \
		if [ ! -d $$dir ]; then \
			echo "Making directory $$dir ..." ;\
			mkdir -p $$dir ;\
		fi ;\
	done ;\
	}

# Compile
$(BUILD_DIR)/%.o : $(SRC_DIR)/%.cpp
	$(CXX) -o $@ $(CXXFLAGS) -c $<

# Link project
$(PROJECT) : $(OBJECTS)
	$(CXX) -o $(BUILD_BIN_DIR)/$@ $(LDFLAGS) $(OBJECTS) $(LIBS)

.PHONY : install
install :
	@{ \
	if [ ! -d $(INSTALL_BIN_DIR) ]; then \
		echo "Creating $(INSTALL_BIN_DIR) ..." ;\
		mkdir -p $(INSTALL_BIN_DIR) ;\
	fi ;\
	echo "Installing at $(INSTALL_DIR) ..." ;\
	cp $(BUILD_BIN_DIR)/* $(INSTALL_BIN_DIR) ;\
	echo "Done" ;\
	}


###############
### Cleanup ###
###############

# Clean up build directory
.PHONY : clean	
clean :
	rm -rf $(BUILD_DIR)/*

# Clean up install directory
.PHONY : clean_install
clean_install :
	rm -rf $(INSTALL_DIR)/*


####################
### Test drivers ###
####################

# Get a list of objects excluding main.o
OBJECTS_EXCEPT_MAIN = $(filter-out $(BUILD_DIR)/main.o,$(OBJECTS))


########################
### Regression Tests ###
########################

.PHONY: test
test:
	@{ \
	cd test ;\
	./run_tests.sh $(START_DIR)/$(BUILD_BIN_DIR)/$(PROJECT) $(echo_failed_diffs) $(is_mpi_enabled) $(mpirun_cmd) ;\
	}

