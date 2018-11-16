# -----------------------------------------------------------------------------
# Make sMDR
# -----------------------------------------------------------------------------

# compiler
CC = g++

ifeq  ($(NODETYPE),powerpc)
        MFLAG = -m64
        EXECUTABLE = sMDR.ppc64
        LIBS = -lmdrppc
else
        EXECUTABLE = sMDR
        LIBS = -lmdr
endif

# compiler flags
MDRDIR = ../mdrlib
INCLUDEDIR = -I. -I$(MDRDIR)
CCFLAGS = -Wall -O2 $(INCLUDEDIR)
CCFLAGS = -Wall -g $(INCLUDEDIR)

#LIBDIR = -L/usr/lib/mysql -L/usr/local/lib -L../lib
LIBDIR = -L../lib
# LIBS = -lmdr

PROG_NAME = sMDR

SMDROBJS = sMDR.o PermutationTester.o sAnalysis.o

SMDRHEADERS = PermutationTester.h sAnalysis.h

all: $(PROG_NAME)

$(PROG_NAME): $(SMDROBJS) $(SMDRHEADERS)
	$(CC) -o ../../bin/$(EXECUTABLE) $(SMDROBJS) $(LIBDIR) $(LIBS)

%.o: %.cpp $(SMDRHEADERS)
	$(CC) -c $(CCFLAGS) $(MFLAG) $< -o $@

tarball:
	tar cvfz sMDR_source_`date '+%Y%m%d'`.tar.gz \
	*.cpp *.h* Makefile*

clean:
	rm -f *.o *~ core
