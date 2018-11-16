# Makefile for datasim library

LIBDIR = ../lib/

ifeq  ($(NODETYPE),powerpc)
        MFLAG = -m64
	LIB = $(LIBDIR)libmdrppc.a
else
	LIB = $(LIBDIR)libmdr.a
endif

CC = g++
# LIB= $(LIBDIR)libmdr.a
CFLAGS= -O2 -Wall -Wno-deprecated
CFLAGS= -g -Wall -Wno-deprecated
OBJS = ComboGenerator.o ConfigFileReader.o Config.o DatafileReader.o \
	Dataset.o FlatIndex.o FormattedOutput.o MapFileReader.o \
	MatrixOutput.o MDR.o MDRExcept.o Model.o OutputWriter.o \
	PartitionOutput.o RulesOutput.o TabDelimOutput.o TrainList.o \
	LogOutput.o Analysis.o Checkpoint.o StatCreator.o LRInteract.o \
	LogistReg.o FileManager.o BinDataFileReader.o Stringmanip.o \
	BioComboGenerator.o KnuthComboGenerator.o BioTextFileAccess.o \
	BioBinFileAccess.o BioFileAccess.o AllModelOutput.o ModelFileReader.o

INTERFACES = ComboGenerator.h ConfigFileReader.h Config.h DatafileReader.h \
	Dataset.h FlatIndex.h FormattedOutput.h MapFileReader.h \
	MatrixOutput.h MDR.h MDRExcept.h Model.h OutputWriter.h \
	PartitionOutput.h RulesOutput.h TabDelimOutput.h TrainList.h \
	LogOutput.h  rbtree.h Analysis.h Checkpoint.h Stat.h \
	BalancedAccStat.h ErrorStat.h OddsRatioStat.h fitness.h \
	StatCreator.h LRInteract.h LogistReg.h FileReader.h \
	FileManager.h BinDataFileReader.h Stringmanip.h BioComboGenerator.h \
	KnuthComboGenerator.h BioFileAccess.h BioTextFileAccess.h \
	BioBinFileAccess.h AllModelOutput.h ModelFileReader.h

# GNU compiler flags (works with g++ and egcs)
CC          = g++
PRELINK     =
LD          = g++ -w
AR          = ar rv
INSTALL     = install -c
RM          = rm
INC_DIRS    = -I.

%.o: %.cpp $(INTERFACES)
	$(CC) $(CFLAGS) $(MFLAG) $(INC_DIRS) -c $< -o $@
$(LIB): $(OBJS)
	$(PRELINK)
	$(AR) $(LIB) $?
	ranlib $(LIB)
	@echo $(LIB) is now up-to-date
.cpp.a:;
clean:
	rm -f *.o *~ core 
