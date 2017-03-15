# LIBRARY SETTINGS - SET AS NECESSARY
# 
# The library settings typically require some tinkering - for reasons beyond me, sometimes one has to include
# the shared object files (.so), and sometimes the .a files (particularly for bamtools).
# Also, curiously, sometimes bamtools requires the explicit inclusion of libz (either as 
# file or just via -lz)
# The following values work for me (see below for an alternative):
#
BOOST_INCLUDE = /home/dilthey/boost/boost_1_63_0/install/include
BOOST_LIB = /home/dilthey/boost/boost_1_63_0/install/lib
BAMTOOLS_INCLUDE = /home/dilthey/bamtools/bamtools/include
BAMTOOLS_SRC = /home/dilthey/bamtools/bamtools/src
BAMTOOLS_LIB = /home/dilthey/bamtools/bamtools/lib

INCS = -I$(BOOST_INCLUDE) -I$(BAMTOOLS_INCLUDE) -I$(BAMTOOLS_SRC)
LIBS = -L$(BOOST_LIB) -L$(BAMTOOLS_LIB) -lboost_random -lboost_filesystem -lboost_system  -lbamtools -lbamtools-utils -lz -lboost_serialization


# an alternative line (courtesy Peter Humburg, not working for me but for him) is
# LIBS = /home/dilthey/PnP/libs/boost_1_59_0/lib/lib/libboost_random.so /home/dilthey/PnP/libs/boost_1_59_0/lib/lib/libboost_filesystem.so /home/dilthey/PnP/libs/boost_1_59_0/lib/lib/libboost_system.so /home/dilthey/bamtools/bamtools/lib/libbamtools.so /home/dilthey/bamtools/bamtools/lib/libbamtools-utils.a -lz

MKDIR_P = mkdir -p

.PHONY: directories
	
# END LIBRARY SETTINGS

#
# object and binary dirs  
#

DIR_OBJ = ../obj
DIR_BIN = ../bin

CXX    = g++
COPTS  = -ggdb -O2 -fopenmp -std=gnu++0x -fstack-protector-all
CFLAGS = 
COMPILE = $(CXX) $(INCS) $(CFLAGS) $(COPTS)
VPATH = Graph:simulator:mapper:mapper/reads:mapper/aligner:mapper/bwa:mapper/bowtie2:Graph/graphSimulator:hla:linearALTs
        
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
        $(DIR_OBJ)/Utilities.o \
        
#
# list executable file names
#
EXECS = HLA-PRG-LA

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

$(EXECS): $(OBJS)
	$(foreach EX, $(EXECS), $(COMPILE) $(EX).cpp -c -o $(DIR_OBJ)/$(EX).o;)
	$(foreach EX, $(EXECS), $(COMPILE) $(OBJS) $(DIR_OBJ)/$(EX).o -o $(DIR_BIN)/$(EX) $(LIBS);)

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

