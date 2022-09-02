# LIBRARY SETTINGS - SET AS NECESSARY
# 
# For recent versions of BamTools (>= 2.5):
#
BOOST_PATH ?= /usr/include/boost
BAMTOOLS_PATH ?= /home/dilthey/bamtools/install
BOOST_INCLUDE = $(BOOST_PATH)/include
BOOST_LIB = $(BOOST_PATH)/lib
BAMTOOLS_INCLUDE = $(BAMTOOLS_PATH)/include/bamtools
BAMTOOLS_SRC = $(BAMTOOLS_PATH)/../src
BAMTOOLS_LIB = $(BAMTOOLS_PATH)/lib
BAMTOOLS_LIB2 = $(BAMTOOLS_PATH)/lib64

INCS = -I$(BOOST_INCLUDE) -I$(BAMTOOLS_INCLUDE) -I$(BAMTOOLS_SRC)
LIBS = -L$(BOOST_LIB) -L$(BAMTOOLS_LIB) -L$(BAMTOOLS_LIB2) -lboost_random -lboost_filesystem -lboost_system  -lbamtools -lz -lboost_serialization

# use the following for older versions of BamTools (e.g. commit https://github.com/pezmaster31/bamtools/commit/2d7685d2aeedd11c46ad3bd67886d9ed65c30f3e)
# ... in which case you probably also have to comment in the '#include "utils/bamtools_utilities.h"' include statement in mapper/processBAM.cpp.
# 
# BOOST_PATH ?= /data/projects/phillippy/software/boost_1_60_0/
# BAMTOOLS_PATH ?= /data/projects/phillippy/software/bamtools
# BOOST_INCLUDE = $(BOOST_PATH)/include
# BOOST_LIB = $(BOOST_PATH)/lib
# BAMTOOLS_INCLUDE = $(BAMTOOLS_PATH)/include
# BAMTOOLS_SRC = $(BAMTOOLS_PATH)/src
# BAMTOOLS_LIB = $(BAMTOOLS_PATH)/lib

# INCS = -I$(BOOST_INCLUDE) -I$(BAMTOOLS_INCLUDE) -I$(BAMTOOLS_SRC)
# LIBS = -L$(BOOST_LIB) -L$(BAMTOOLS_LIB) -lboost_random -lboost_filesystem -lboost_system  -lbamtools -lbamtools-utils -lz -lboost_serialization


MKDIR_P = mkdir -p

.PHONY: directories
	
# END LIBRARY SETTINGS

#
# object and binary dirs  
#

DIR_OBJ = ../obj
DIR_BIN = ../bin

#CXX    = g++ 
COPTS  = -ggdb -O2 -fopenmp -std=gnu++0x -fstack-protector-all
CFLAGS = 
COMPILE = $(CXX) $(INCS) $(CFLAGS) $(COPTS)
VPATH = Graph:simulator:mapper:mapper/reads:mapper/aligner:mapper/bwa:mapper/bowtie2:Graph/graphSimulator:hla:linearALTs:fullLengthHMM
        
OBJS = \
        $(DIR_OBJ)/Edge.o \
        $(DIR_OBJ)/Graph.o \
        $(DIR_OBJ)/GraphAndEdgeIndex.o \
        $(DIR_OBJ)/Node.o \
        $(DIR_OBJ)/HaplotypePanel.o \
        $(DIR_OBJ)/simpleGraphSimulator.o \
        $(DIR_OBJ)/LocusCodeAllocation.o \
        $(DIR_OBJ)/simulator.o \
        $(DIR_OBJ)/readSimulator.o \
        $(DIR_OBJ)/trueReadLevels.o \
        $(DIR_OBJ)/processBAM.o \
        $(DIR_OBJ)/protoSeeds.o \
        $(DIR_OBJ)/oneRead.o \
        $(DIR_OBJ)/oneReadPair.o \
        $(DIR_OBJ)/oneReadPairwithSeedChains.o \
        $(DIR_OBJ)/alignerBase.o \
        $(DIR_OBJ)/alignmentContext.o \
        $(DIR_OBJ)/extensionAligner.o \
        $(DIR_OBJ)/BWAmapper.o \
        $(DIR_OBJ)/Bowtie2mapper.o \
        $(DIR_OBJ)/statistics.o \
        $(DIR_OBJ)/VirtualNWUnique.o \
        $(DIR_OBJ)/PRGContigBAMAlignment.o \
        $(DIR_OBJ)/verboseSeedChain.o \
        $(DIR_OBJ)/HLATyper.o \
        $(DIR_OBJ)/pathFinder.o \
        $(DIR_OBJ)/linearALTs.o \
        $(DIR_OBJ)/oneExonPosition.o \
        $(DIR_OBJ)/seedChain.o \
        $(DIR_OBJ)/fullLengthHMM.o \
        $(DIR_OBJ)/Utilities.o \
        
#
# list executable file names
#
EXECS = HLA-LA sam2alignment

OUT_DIR = ../obj ../bin

directories: ${OUT_DIR}


#
# compile and link
#
default:
	@echo
	@echo " to build:"
	@echo "    make all"
	@echo
	@echo " to clean:"
	@echo "    make clean"
	@echo "    make realclean"
	@echo

all: directories $(EXECS)

HLA-LA: $(OBJS)
	$(foreach EX, HLA-LA, $(COMPILE) $(EX).cpp -c -o $(DIR_OBJ)/$(EX).o;)
	$(foreach EX, HLA-LA, $(COMPILE) $(OBJS) $(DIR_OBJ)/$(EX).o -o $(DIR_BIN)/$(EX) $(LIBS);)

sam2alignment:
	$(foreach EX, sam2alignment, $(COMPILE) $(EX).cpp -c -o $(DIR_OBJ)/$(EX).o;)
	$(foreach EX, sam2alignment, $(COMPILE) $(DIR_OBJ)/$(EX).o -o $(DIR_BIN)/$(EX);)

$(DIR_OBJ)/%.o: %.cpp %.h
	$(COMPILE) $< -c -o $@


#
# odds and ends
#
clean:
	/bin/rm $(DIR_OBJ)/*

realclean: clean
	/bin/rm $(DIR_BIN)/*

${OUT_DIR}:
	${MKDIR_P} ${OUT_DIR}

