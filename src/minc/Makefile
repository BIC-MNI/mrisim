##############################################################################
# MRISIM MINC Makefile
#
# To use:  set up directories and make compiler selection
#
# May 16, 1996
# R. Kwan 
#
##############################################################################

#
# Directories
#

# Where are MINC header files?
MINC_INCLUDE      = /usr/local/include/
# Where are mrisim library headers?
MRISIM_SIGNAL_DIR = ../signal
# Where is the unit test directory?
TD                = ./Tests
# Where is this makefile?
CURRENT_DIR       = ../mrisrc/minc
# Where is the example program directory?
PROG_DIR          = ./Progs

#
# SGI CC Compiler
#
SGI_DEBUG_FLAGS   = -g -DDEBUG -fullwarn -I/usr/include/ -I$(MINC_INCLUDE) -I../
SGI_RELEASE_FLAGS = -O -fullwarn -I/usr/include/ -I$(MINC_INCLUDE) -I../
SGI_LIBS          = -lminc -lnetcdf -lsun -lc_s -lm
SGI_DEBUG_CXX     = CC $(SGI_DEBUG_FLAGS)
SGI_DEBUG_CC      = cc $(SGI_DEBUG_FLAGS)
SGI_RELEASE_CXX   = CC $(SGI_RELEASE_FLAGS)
SGI_RELEASE_CC    = cc $(SGI_RELEASE_FLAGS)

#
# GNU C++ Compiler
#
GNU_DEBUG_FLAGS   = -gstabs -DDEBUG -I$(MINC_INCLUDE) -I../
GNU_RELEASE_FLAGS = -O -I$(MINC_INCLUDE) -I../
GNU_LIBS          = -lg++ -lminc -lnetcdf -lsun -lc_s -lm
GNU_DEBUG_CXX     = g++ $(GNU_DEBUG_FLAGS)
GNU_DEBUG_CC      = gcc $(GNU_DEBUG_FLAGS)
GNU_RELEASE_CXX   = g++ $(GNU_RELEASE_FLAGS)
GNU_RELEASE_CC    = gcc $(GNU_RELEASE_FLAGS)

##############################################################################
# --- Compiler selection ---
# Set CXX to one of (SGI_DEBUG_CXX, SGI_RELEASE_CXX, 
#                    GNU_DEBUG_CXX, GNU_RELEASE_CXX)
# Set CC  to one of (SGI_DEBUG_CC, SGI_RELEASE_CC,
#                    GNU_DEBUG_CC, GNU_RELEASE_CC)
##############################################################################

CXX  = $(SGI_RELEASE_CXX)
CC   = $(SGI_RELEASE_CC)
LIBS = $(SGI_LIBS)
GET  = co

##############################################################################
# Object models
##############################################################################

MINC_OBJS     = mincicv.o mincfile.o imincfile.o omincfile.o iomincfile.o \
                time_stamp.o 
MVOL_OBJS     = mrimatrix.o fourn.o mrivolume.o mristring.o mriimage.o \
                mrilabel.o chirp.o
TESTS         = mincinfo minccopy mincstat testmat testchirp
OBJS          = $(MVOL_OBJS) $(MINC_OBJS)

all:   $(TESTS) 
lib:   libmriminc.a

#
# Objects
#

mriminc.h:
	$(GET) mriminc.h

mincicv.h:
	$(GET) mincicv.h
mincicv.cxx:
	$(GET) mincicv.cxx
mincicv.o:	mincicv.h mincicv.cxx
	$(CXX) -c mincicv.cxx -o mincicv.o

mincfile.h:
	$(GET) mincfile.h
mincfile.cxx:
	$(GET) mincfile.cxx
mincfile.o:	mincfile.h mincfile.cxx mincicv.o
	$(CXX) -c mincfile.cxx -o mincfile.o

imincfile.h:
	$(GET) imincfile.h
imincfile.cxx:
	$(GET) imincfile.cxx
imincfile.o:	imincfile.h imincfile.cxx mincfile.o $(MVOL_OBJS)
	$(CXX) -c imincfile.cxx -o imincfile.o

omincfile.h:
	$(GET) omincfile.h
omincfile.cxx:
	$(GET) omincfile.cxx
omincfile.o:	omincfile.h omincfile.cxx mincfile.o
	$(CXX) -c omincfile.cxx -o omincfile.o

iomincfile.h:
	$(GET) iomincfile.h
iomincfile.cxx:
	$(GET) iomincfile.cxx
iomincfile.o:	iomincfile.h iomincfile.cxx imincfile.o omincfile.o
	$(CXX) -c iomincfile.cxx -o iomincfile.o

mrimatrix.h:
	$(GET) mrimatrix.h
mrimatrix.cxx:
	$(GET) mrimatrix.cxx
mrimatrix.o:	mrimatrix.h mrimatrix.cxx fourn.o
	$(CXX) -c mrimatrix.cxx -o mrimatrix.o

mrivolume.h:
	$(GET) mrivolume.h
mrivolume.cxx:
	$(GET) mrivolume.cxx
mrivolume.o:	mrivolume.h mrivolume.cxx
	$(CXX) -c mrivolume.cxx -o mrivolume.o

mristring.h:
	$(GET) mristring.h
mristring.cxx:
	$(GET) mristring.cxx
mristring.o:	mristring.h mristring.cxx
	$(CXX) -c mristring.cxx -o mristring.o

mriimage.h:
	$(GET) mriimage.h
mriimage.cxx:
	$(GET) mriimage.cxx
mriimage.o:	mriimage.h mriimage.cxx mrimatrix.o
	$(CXX) -c mriimage.cxx -o mriimage.o

mrilabel.h:
	$(GET) mrilabel.h
mrilabel.cxx:
	$(GET) mrilabel.cxx
mrilabel.o:	mrilabel.h mrilabel.cxx mrimatrix.o
	$(CXX) -c mrilabel.cxx -o mrilabel.o

fourn.h:
	$(GET) fourn.h
fourn.c: 
	$(GET) fourn.c 
fourn.o:	fourn.h fourn.c  
	$(CC) -DTRACE -c fourn.c -o fourn.o

chirp.h:
	$(GET) chirp.h
chirp.cxx:
	$(GET) chirp.cxx
chirp.o:	chirp.h chirp.cxx fourn.o
	$(CXX) -c chirp.cxx -o chirp.o

time_stamp.h:
	$(GET) time_stamp.h
time_stamp.c:
	$(GET) time_stamp.c
time_stamp.o:	time_stamp.h time_stamp.c
	$(CC) -c time_stamp.c -o time_stamp.o

parseargv.o:	../mrisim/ParseArgv.h ../mrisim/ParseArgv.c
	$(CC) -c ../mrisim/ParseArgv.c -o parseargv.o

#
# Example application programs
#

mincinfo:	$(PROG_DIR)/mincinfo.cxx $(MINC_OBJS) $(MVOL_OBJS)
	$(CXX) $(PROG_DIR)/mincinfo.cxx $(MINC_OBJS) $(MVOL_OBJS) \
               -o mincinfo $(LIBS)

minccopy:	$(PROG_DIR)/minccopy.cxx $(MINC_OBJS) $(MVOL_OBJS)
	$(CXX) $(PROG_DIR)/minccopy.cxx $(MINC_OBJS) $(MVOL_OBJS) \
               -o minccopy $(LIBS)

mincdiff:	$(PROG_DIR)/mincdiff.cxx $(MINC_OBJS) $(MVOL_OBJS)
	$(CXX) $(PROG_DIR)/mincdiff.cxx $(MINC_OBJS) $(MVOL_OBJS) \
               -o mincdiff $(LIBS)

maketest:	$(PROG_DIR)/maketest.cxx $(MINC_OBJS) $(MVOL_OBJS)
	$(CXX) $(PROG_DIR)/maketest.cxx $(MINC_OBJS) $(MVOL_OBJS) \
               -o maketest $(LIBS)

mincstat:	$(PROG_DIR)/mincstat.cxx parseargv.o $(MINC_OBJS) $(MVOL_OBJS)
	$(CXX) $(PROG_DIR)/mincstat.cxx $(MINC_OBJS) $(MVOL_OBJS)  \
               parseargv.o -o mincstat $(LIBS)

minctest:	$(TD)/minctest.cxx mrimatrix.o $(MINC_OBJS) $(MVOL_OBJS)
	$(CXX) $(TD)/minctest.cxx $(MINC_OBJS) $(MVOL_OBJS)  -o minctest $(LIBS)

testmat:	$(TD)/testmat.cxx $(MVOL_OBJS)
	$(CXX) $(TD)/testmat.cxx $(MVOL_OBJS) -o testmat $(LIBS)

testvol:	$(TD)/testvol.cxx $(MINC_OBJS) $(MVOL_OBJS)
	$(CXX) $(TD)/testvol.cxx $(MINC_OBJS) $(MVOL_OBJS) -o testvol $(LIBS)

testchirp:	$(TD)/testchirp.cxx mrimatrix.o chirp.o fourn.o
	$(CXX) $(TD)/testchirp.cxx mrimatrix.o chirp.o fourn.o \
               -o testchirp $(LIBS)

testfczt:	$(TD)/testchirp.cxx mrimatrix.o chirp.o fourn.o
	$(CXX) -DFOURN $(TD)/testchirp.cxx mrimatrix.o chirp.o fourn.o \
               -o testfczt $(LIBS)

tracechirp:	$(TD)/tracechirp.cxx fourn.o
	$(CXX) $(TD)/tracechirp.cxx fourn.o -o tracechirp $(LIBS)

loadphan:	$(TD)/loadphan.cxx $(MINC_OBJS) $(MVOL_OBJS)
	$(CXX) $(TD)/loadphan.cxx $(MINC_OBJS) $(MVOL_OBJS) -o loadphan $(LIBS)

minc:		$(TD)/minc.c
	$(CC) $(TD)/minc.c -o minc $(LIBS)

libmriminc.a:	$(OBJS)
	rm -f libmriminc.a
	ar ruv libmriminc.a $(OBJS)
	ranlib libmriminc.a
clean:
	rm -f *.o *~ *.a
