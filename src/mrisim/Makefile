##############################################################################
# MRISIM Makefile
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
MRISIM_MINC_DIR   = ../minc
MRISIM_SIGNAL_DIR = ../signal
# Where is the unit test directory?
TD                = ./Tests
# Where is this makefile?
CURRENT_DIR       = ../mrisim

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

MRISIM_MINC_LIB   = $(MRISIM_MINC_DIR)/libmriminc.a
MRISIM_SIGNAL_LIB = $(MRISIM_SIGNAL_DIR)/libsignal.a
MRLIBS            = $(MRISIM_MINC_LIB) $(MRISIM_SIGNAL_LIB)

SCANNER  = mriscanner.o scanner_output.o
RF_COIL  = rf_coil.o intrinsic_coil.o image_snr_coil.o percent_coil.o
PHAN     = phantom.o tissue_phantom.o
RF_PHAN  = rf_tissue_phantom.o 
DISCRETE = $(PHAN) discrete_label_phantom.o discrete_phantom.o
FUZZY    = $(DISCRETE) fuzzy_label_phantom.o fuzzy_phantom.o
DISCR_RF = $(DISCRETE) $(RF_PHAN) discrete_rf_phantom.o
FUZZY_RF = $(FUZZY) $(RF_PHAN) discrete_rf_phantom.o fuzzy_rf_phantom.o
SUPPORT  = paramfile.o ParseArgv.o mrisimargs.o

OBJS     = $(FUZZY_RF) $(RF_COIL) $(SCANNER) $(SUPPORT)

TD       = ./Tests

UNIT_TESTS =  d.test drf.test f.test frf.test 
INT_TESTS  =  dpv.test select.test linselect.test

##############################################################################
# --- Dependencies ---
##############################################################################

all:   $(OBJS)

scanner:   $(SCANNER)
coil:      $(RF_COIL)
phantom:   $(DISCRETE) $(DISCR_RF)

unit_test: $(UNIT_TESTS)

clean:
	rm -f *.o *~ *.a

clean_unit:
	rm -f $(TD)/*.test


#
# --- Libraries --- 
#

$(MRISIM_MINC_LIB):	
	cd $(MRISIM_MINC_DIR); make lib; cd $(CURRENT_DIR)

$(MRISIM_SIGNAL_LIB):
	cd $(MRISIM_SIGNAL_DIR); make lib; cd $(CURRENT_DIR)
#
# --- Objects --- 
#

mrisim.h:
	$(GET) mrisim.h

$(MRISIM_MINC_DIR)/mriminc.h:
	cd $(MRISIM_MINC_DIR); make mriminc.h; cd $(CURRENT_DIR)

# --- SCANNER --- 

mriscanner.h:
	$(GET) mriscanner.h
mriscanner.cxx:
	$(GET) mriscanner.cxx
mriscanner.o:	mriscanner.h mriscanner.cxx rf_coil.h phantom.h
	$(CXX) -c mriscanner.cxx -o mriscanner.o

scanner_output.h:
	$(GET) scanner_output.h
scanner_output.cxx:
	$(GET) scanner_output.cxx
scanner_output.o: scanner_output.h scanner_output.cxx mriscanner.o 
	$(CXX) -c scanner_output.cxx -o scanner_output.o

# --- RF_COIL ---

rf_coil.h:
	$(GET) rf_coil.h
rf_coil.cxx:
	$(GET) rf_coil.cxx
rf_coil.o:	rf_coil.h rf_coil.cxx
	$(CXX) -c rf_coil.cxx -o rf_coil.o

intrinsic_coil.h:
	$(GET) intrinsic_coil.h
intrinsic_coil.cxx:
	$(GET) intrinsic_coil.cxx
intrinsic_coil.o:       intrinsic_coil.h intrinsic_coil.cxx rf_coil.o
	$(CXX) -c intrinsic_coil.cxx -o intrinsic_coil.o

image_snr_coil.h:
	$(GET) image_snr_coil.h
image_snr_coil.cxx:
	$(GET) image_snr_coil.cxx
image_snr_coil.o:       image_snr_coil.h image_snr_coil.cxx rf_coil.o
	$(CXX) -c image_snr_coil.cxx -o image_snr_coil.o

percent_coil.h:
	$(GET) percent_coil.h
percent_coil.cxx:       
	$(GET) percent_coil.cxx
percent_coil.o: percent_coil.h percent_coil.cxx rf_coil.o
	$(CXX) -c percent_coil.cxx -o percent_coil.o

# --- PHANTOM ---

phantom.h:
	$(GET) phantom.h
phantom.cxx:
	$(GET) phantom.cxx
phantom.o:	phantom.h phantom.cxx
	$(CXX) -c phantom.cxx -o phantom.o

tissue_phantom.h:
	$(GET) tissue_phantom.h
tissue_phantom.cxx:
	$(GET) tissue_phantom.cxx
tissue_phantom.o:	tissue_phantom.h tissue_phantom.cxx
	$(CXX) -c tissue_phantom.cxx -o tissue_phantom.o

rf_tissue_phantom.h:
	$(GET) rf_tissue_phantom.h
rf_tissue_phantom.cxx:
	$(GET) rf_tissue_phantom.cxx
rf_tissue_phantom.o:	rf_tissue_phantom.h rf_tissue_phantom.cxx rf_coil.o
	$(CXX) -c rf_tissue_phantom.cxx -o rf_tissue_phantom.o

discrete_label_phantom.h:
	$(GET) discrete_label_phantom.h
discrete_label_phantom.cxx:
	$(GET) discrete_label_phantom.cxx
discrete_label_phantom.o:	discrete_label_phantom.h \
                        discrete_label_phantom.cxx \
                        phantom.o tissue_phantom.o
	$(CXX) -c discrete_label_phantom.cxx -o discrete_label_phantom.o

discrete_phantom.h:
	$(GET) discrete_phantom.h
discrete_phantom.cxx:
	$(GET) discrete_phantom.cxx
discrete_phantom.o:	discrete_phantom.h discrete_phantom.cxx \
                        phantom.o tissue_phantom.o discrete_label_phantom.o
	$(CXX) -c discrete_phantom.cxx -o discrete_phantom.o

fuzzy_label_phantom.h:
	$(GET) fuzzy_label_phantom.h
fuzzy_label_phantom.cxx:
	$(GET) fuzzy_label_phantom.cxx
fuzzy_label_phantom.o: fuzzy_label_phantom.h fuzzy_label_phantom.cxx \
                       phantom.o tissue_phantom.o
	$(CXX) -c fuzzy_label_phantom.cxx -o fuzzy_label_phantom.o

fuzzy_phantom.h:
	$(GET) fuzzy_phantom.h
fuzzy_phantom.cxx:
	$(GET) fuzzy_phantom.cxx
fuzzy_phantom.o: fuzzy_phantom.h fuzzy_phantom.cxx \
                 phantom.o tissue_phantom.o fuzzy_label_phantom.o
	$(CXX) -c fuzzy_phantom.cxx -o fuzzy_phantom.o

discrete_rf_phantom.h:
	$(GET) discrete_rf_phantom.h
discrete_rf_phantom.cxx:
	$(GET) discrete_rf_phantom.cxx
discrete_rf_phantom.o:	discrete_rf_phantom.h discrete_rf_phantom.cxx \
                        discrete_label_phantom.o rf_tissue_phantom.o
	$(CXX) -c discrete_rf_phantom.cxx -o discrete_rf_phantom.o

fuzzy_rf_phantom.h:
	$(GET) fuzzy_rf_phantom.h
fuzzy_rf_phantom.cxx:
	$(GET) fuzzy_rf_phantom.cxx
fuzzy_rf_phantom.o: fuzzy_rf_phantom.h fuzzy_rf_phantom.cxx \
                    fuzzy_label_phantom.o rf_tissue_phantom.o
	$(CXX) -c fuzzy_rf_phantom.cxx -o fuzzy_rf_phantom.o

# --- SUPPORT ---

paramfile.h:
	$(GET) paramfile.h
paramfile.cxx:
	$(GET) paramfile.cxx
paramfile.o:    paramfile.h paramfile.cxx
	$(CXX) -c paramfile.cxx -o paramfile.o

ParseArgv.h:
	$(GET) ParseArgv.h
ParseArgv.c:
	$(GET) ParseArgv.c 
ParseArgv.o:    ParseArgv.h ParseArgv.c
	$(CC) -c ParseArgv.c -o ParseArgv.o

mrisimargs.h:
	$(GET) mrisimargs.h
mrisimargs.cxx:
	$(GET) mrisimargs.cxx
mrisimargs.o:   mrisimargs.h mrisimargs.cxx ParseArgv.o
	$(CXX) -c mrisimargs.cxx -o mrisimargs.o

# --- MAIN EXECUTABLE ---

mrisim_main.h:  $(MRISIM_MINC_DIR)/mriminc.h
	$(GET) mrisim_main.h
mrisim_main.cxx:
	$(GET) mrisim_main.cxx
mrisim:         $(MRLIBS) mrisim_main.h mrisim_main.cxx $(OBJS)
	$(CXX) mrisim_main.cxx $(OBJS) -o mrisim $(MRLIBS) $(LIBS)

#
# --- UNIT TESTS ---
#

d.test:	$(TD)/d.cxx $(DISCRETE)
	$(CXX) $(TD)/d.cxx $(DISCRETE) $(MRLIBS) $(LIBS) -o $(TD)/d.test  

drf.test:	$(TD)/drf.cxx $(DISCR_RF) 
	$(CXX) $(TD)/drf.cxx $(RF_COIL) $(DISCR_RF) $(MRLIBS) $(LIBS) \
               -o $(TD)/drf.test  

f.test:	$(TD)/f.cxx $(FUZZY)
	$(CXX) $(TD)/f.cxx $(FUZZY) $(MRLIBS) $(LIBS) -o $(TD)/f.test  

frf.test:	$(TD)/frf.cxx $(FUZZY_RF)
	$(CXX) $(TD)/frf.cxx $(RF_COIL) $(FUZZY_RF) $(MRLIBS) $(LIBS) \
               -o $(TD)/frf.test  

dpv.test:	$(TD)/dpv.cxx $(DISCRETE)
	$(CXX) $(TD)/dpv.cxx $(DISCRETE) $(MRLIBS) $(LIBS) -o $(TD)/dpv.test  

select.test:	$(TD)/select.cxx $(DISCRETE)
	$(CXX) $(TD)/select.cxx $(DISCRETE) $(MRLIBS) $(LIBS) -o $(TD)/select.test  

linselect.test:	$(TD)/linselect.cxx $(DISCRETE)
	$(CXX) $(TD)/linselect.cxx $(DISCRETE) $(MRLIBS) $(LIBS) -o $(TD)/linselect.test  
