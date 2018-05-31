# -----------------------------------------------------------------------------
# Make pMDR
# -----------------------------------------------------------------------------

CC = mpicxx
EXECUTABLE = pMDR
LIBS = -lmdr

# compiler
# CC = g++

# compiler flags
MDRDIR = ../mdrlib
INCLUDEDIR = -I. -I$(MDRDIR)
CCFLAGS = -Wall -O2 $(INCLUDEDIR)

LIBDIR = -L../lib

PROG_NAME = pMDR

PMDROBJS = pMDR.o ParallelPerm.o MasterAnalysis.o SlaveAnalysis.o pBioComboGenerator.o \
	pKnuthComboGenerator.o 

PMDRHEADERS = ParallelPerm.h MasterAnalysis.h SlaveAnalysis.h pBioComboGenerator.h \
	pKnuthComboGenerator.h pComboTransfer.h

all: $(PROG_NAME)

$(PROG_NAME): $(PMDROBJS) $(PMDRHEADERS)
	$(CC) -o ../../bin/$(EXECUTABLE) $(PMDROBJS) $(LIBDIR) $(LIBS)

%.o: %.cpp $(PMDRHEADERS)
	$(CC) -c $(CCFLAGS) $(MFLAG) $< -o $@

tarball:
	tar cvfz sMDR_source_`date '+%Y%m%d'`.tar.gz \
	*.cpp *.h* Makefile*

clean:
	rm -f *.o *~ core
