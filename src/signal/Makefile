##############################################################################
# MRISIM signal Makefile
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
# Where are the mrisim header files?
MRISIM_DIR        = ../mrisim
# Where is the unit test directory?
TD                = ./Tests
# Where is this makefile?
CURRENT_DIR       = ../signal

#
# SGI CC Compiler
#
SGI_DEBUG_FLAGS   = -g -DDEBUG -fullwarn -I/usr/include/ -I../
SGI_RELEASE_FLAGS = -O -fullwarn -I/usr/include/ -I../
SGI_LIBS          = -lm
SGI_DEBUG_CXX     = CC $(SGI_DEBUG_FLAGS)
SGI_DEBUG_CC      = cc $(SGI_DEBUG_FLAGS)
SGI_RELEASE_CXX   = CC $(SGI_RELEASE_FLAGS)
SGI_RELEASE_CC    = cc $(SGI_RELEASE_FLAGS)

#
# GNU C++ Compiler
#
GNU_DEBUG_FLAGS   = -gstabs -DDEBUG -I../
GNU_RELEASE_FLAGS = -O -I../
GNU_LIBS          = -lg++ -lm
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

QUICK    = quickseq.o quick_model.o \
           se.o ir.o ffe.o spoiled_flash.o fisp.o flash.o ce_fast.o

CUSTOM   = customseq.o vector_model.o isochromat_model.o \
           fast_iso_model.o rf_pulse.o repeat.o spoiler.o

SUPPORT  = $(MRISIM_MINC_DIR)/mristring.o tissue.o vector.o \
           spin_model.o \
           pulseseq.o event.o sample.o

OBJS     = $(SUPPORT) $(QUICK) $(CUSTOM)

TESTS    = $(TD)/spinecho $(TD)/invrecovery $(TD)/ffetest

##############################################################################
# --- Dependencies ---
##############################################################################

all:    $(OBJS)
tests:  $(TESTS)
lib:    libsignal.a

#
# Class dependencies
#

# Basic utility classes

$(MRISIM_MINC_DIR)/mristring.o:
	cd $(MRISIM_MINC_DIR); make mristring.o; cd $(CURRENT_DIR)

$(MRISIM_MINC_DIR)/fourn.o:
	cd $(MRISIM_MINC_DIR); make fourn.o; cd $(CURRENT_DIR)

$(MRISIM_DIR)/mrisim.h:
	cd $(MRISIM_DIR); make mrisim.h; cd $(CURRENT_DIR)

vector.h:
	$(GET) vector.h
vector.cxx:
	$(GET) vector.cxx
vector.o:     vector.cxx vector.h;
	$(CXX) -c vector.cxx -o vector.o

tissue.h:
	$(GET) tissue.h
tissue.cxx:
	$(GET) tissue.cxx
tissue.o:     tissue.cxx tissue.h $(MRISIM_DIR)/mrisim.h ../minc/mristring.o;
	$(CXX) -c tissue.cxx -o tissue.o

# Old magnetization spin models
# Obsolete
#magvec.o:     magvec.cxx magvec.h vector.h tissue.h;
#$(CXX) -c magvec.cxx -o magvec.o
#
#magsys.o:     magsys.cxx magvec.h vector.h tissue.h;
#$(CXX) -c magsys.cxx -o magsys.o

# Magnetization spin models

spin_model.h:
	$(GET) spin_model.h
spin_model.cxx:
	$(GET) spin_model.cxx
spin_model.o:   spin_model.cxx spin_model.h;
	$(CXX) -c spin_model.cxx -o spin_model.o

quick_model.h:
	$(GET) quick_model.h
quick_model.cxx:
	$(GET) quick_model.cxx
quick_model.o:  quick_model.cxx quick_model.h spin_model.o tissue.o;
	$(CXX) -c quick_model.cxx -o quick_model.o

vector_model.h:
	$(GET) vector_model.h
vector_model.cxx:
	$(GET) vector_model.cxx
vector_model.o:	vector_model.cxx vector_model.h spin_model.o tissue.o;
	$(CXX) -c vector_model.cxx -o vector_model.o

isochromat_model.h:
	$(GET) isochromat_model.h
isochromat_model.cxx:
	$(GET) isochromat_model.cxx
isochromat_model.o:	isochromat_model.cxx isochromat_model.h spin_model.o tissue.o;
	$(CXX) -c isochromat_model.cxx -o isochromat_model.o

fast_iso_model.h:
	$(GET) fast_iso_model.h
fast_iso_model.cxx:
	$(GET) fast_iso_model.cxx
fast_iso_model.o:	fast_iso_model.cxx fast_iso_model.h spin_model.o tissue.o;
	$(CXX) -c fast_iso_model.cxx -o fast_iso_model.o

# Programmable pulse sequence events

event.h:
	$(GET) event.h
event.cxx:
	$(GET) event.cxx
event.o:      event.h vector.o spin_model.o;
	$(CXX) -c event.cxx -o event.o

rf_pulse.h:
	$(GET) rf_pulse.h
rf_pulse.cxx:
	$(GET) rf_pulse.cxx
rf_pulse.o:   rf_pulse.cxx rf_pulse.h event.o;
	$(CXX) -c rf_pulse.cxx -o rf_pulse.o

spoiler.h:
	$(GET) spoiler.h
spoiler.cxx:
	$(GET) spoiler.cxx
spoiler.o:    spoiler.cxx spoiler.h event.o;
	$(CXX) -c spoiler.cxx -o spoiler.o

sample.h:	
	$(GET) sample.h
sample.cxx:
	$(GET) sample.cxx
sample.o:     sample.cxx sample.h event.o;
	$(CXX) -c sample.cxx -o sample.o 

repeat.h:
	$(GET) repeat.h
repeat.cxx:
	$(GET) repeat.cxx
repeat.o:     repeat.cxx repeat.h event.o;
	$(CXX) -c repeat.cxx -o repeat.o 

# Pulse Sequence classes

pulseseq.h:
	$(GET) pulseseq.h
pulseseq.cxx:
	$(GET) pulseseq.cxx 
pulseseq.o:   pulseseq.cxx pulseseq.h event.o;
	$(CXX) -c pulseseq.cxx -o pulseseq.o

customseq.h:
	$(GET) customseq.h
customseq.cxx:
	$(GET) customseq.cxx
customseq.o:  customseq.cxx customseq.h pulseseq.o;
	$(CXX) -c customseq.cxx -o customseq.o

quickseq.h:
	$(GET) quickseq.h
quickseq.cxx:
	$(GET) quickseq.cxx
quickseq.o:   quickseq.cxx quickseq.h pulseseq.o quick_model.h sample.h;
	$(CXX) -c quickseq.cxx -o quickseq.o

se.h:
	$(GET) se.h
se.cxx:	
	$(GET) se.cxx
se.o:         se.cxx se.h quickseq.o;
	$(CXX) -c se.cxx -o se.o

ir.h:
	$(GET) ir.h
ir.cxx: 
	$(GET) ir.cxx
ir.o:         ir.cxx ir.h quickseq.o;
	$(CXX) -c ir.cxx -o ir.o

ffe.h:
	$(GET) ffe.h
ffe.cxx:
	$(GET) ffe.cxx
ffe.o:         ffe.cxx ffe.h quickseq.o;
	$(CXX) -c ffe.cxx -o ffe.o

spoiled_flash.h:
	$(GET) spoiled_flash.h
spoiled_flash.cxx:
	$(GET) spoiled_flash.cxx
spoiled_flash.o:         spoiled_flash.cxx spoiled_flash.h ffe.o;
	$(CXX) -c spoiled_flash.cxx -o spoiled_flash.o

fisp.h:
	$(GET) fisp.h
fisp.cxx:
	$(GET) fisp.cxx
fisp.o:         fisp.cxx fisp.h ffe.o;
	$(CXX) -c fisp.cxx -o fisp.o

flash.h:	
	$(GET) flash.h
flash.cxx:
	$(GET) flash.cxx
flash.o:         flash.cxx flash.h ffe.o;
	$(CXX) -c flash.cxx -o flash.o

ce_fast.h:	
	$(GET) ce_fast.h
ce_fast.cxx:
	$(GET) ce_fast.cxx
ce_fast.o:         ce_fast.cxx ce_fast.h ffe.o;
	$(CXX) -c ce_fast.cxx -o ce_fast.o

signal.h:
	$(GET) signal.h

# Signal calculation library

libsignal.a: $(OBJS) signal.h
	rm -f libsignal.a
	ar ruv libsignal.a $(OBJS)
	ranlib libsignal.a

#
# TESTS
#

$(TD)/quick:	$(TD)/quick.cxx $(QUICK) $(SUPPORT)
	$(CXX) $(TD)/quick.cxx $(QUICK) $(SUPPORT) $(LIBS) -o $(TD)/quick

$(TD)/spinecho:   $(TD)/spinecho.cxx $(OBJS) $(MRISIM_MINC_DIR)/fourn.o
	$(CXX) $(TD)/spinecho.cxx $(OBJS) $(MRISIM_MINC_DIR)/fourn.o \
               $(LIBS) -o $(TD)/spinecho

$(TD)/invrecovery:   $(TD)/invrecovery.cxx $(MRISIM_MINC_DIR)/mristring.o \
                     $(OBJS)
	$(CXX) $(TD)/invrecovery.cxx $(MRISIM_MINC_DIR)/mristring.o \
               $(OBJS) $(LIBS) -o $(TD)/invrecovery

$(TD)/ffetest:    $(TD)/ffetest.cxx $(QUICK) $(SUPPORT)
	$(CXX) $(TD)/ffetest.cxx $(QUICK) $(SUPPORT) $(LIBS) -o $(TD)/ffetest

$(TD)/isotest: $(TD)/isotest.cxx $(OBJS) $(MRISIM_MINC_DIR)/fourn.o
	$(CXX) $(TD)/isotest.cxx $(OBJS) $(MRISIM_MINC_DIR)/fourn.o \
               $(LIBS) -o $(TD)/isotest

$(TD)/fftsim: $(TD)/fftsim.cxx $(OBJS) $(MRISIM_MINC_DIR)/fourn.o
	$(CXX) $(TD)/fftsim.cxx $(OBJS) $(MRISIM_MINC_DIR)/fourn.o \
               $(LIBS) -o $(TD)/fftsim

$(TD)/fastiso: $(TD)/fastiso.cxx $(OBJS) $(MRISIM_MINC_DIR)/fourn.o
	$(CXX) $(TD)/fastiso.cxx $(OBJS) $(MRISIM_MINC_DIR)/fourn.o \
               $(LIBS) -o $(TD)/fastiso

#
#  CLEAN
#

clean: 
	rm -f *.o *~ *.a

rcsclean:
	rcsclean
clean_tests:
	rm -f $(TESTS); rm -f $(TD)/*~
