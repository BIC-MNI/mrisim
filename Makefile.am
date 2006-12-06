AUTOMAKE_OPTIONS = check-news subdir-objects
ACLOCAL_AMFLAGS = -I m4

AM_CPPFLAGS =
AM_CFLAGS = 
AM_CXXFLAGS = -Wno-deprecated

SUBDIRS = .

# Add these to the include searches.
#
INCLUDES = -I$(top_srcdir)/src

EXTRA_DIST = $(m4_files) \
	$(coil_files) \
	$(doc_files) \
	$(field_files) \
	$(sequence_files) \
	epm-header.in

m4_files = \
	m4/mni_REQUIRE_LIB.m4 \
	m4/mni_REQUIRE_MNILIBS.m4 \
	m4/smr_WITH_BUILD_PATH.m4

coil_files = \
	coils/coil_pd_pn2_1mm.rf \
	coils/coil.rf \
	coils/coil_t1_intrinsic.rf \
	coils/coil_t1_pn2_1mm.rf \
	coils/coil_t2_pn2_1mm.rf

doc_files = \
	doc/mrisim.man \
	doc/mrisim.ps \
	doc/switch.txt

field_files = \
	fields/Philips/makerf.zip \
	fields/Philips/qbc_id.mnc.Z \
	fields/Philips/qbc_q450.mnc.Z \
	fields/Philips/qbc_q86.mnc.Z \
	fields/Philips/qhc.mnc.Z

sequence_files = \
	sequences/pd_ai.seq \
	sequences/pd_icbm.seq \
	sequences/t1_ai.seq \
	sequences/t1_icbm_avg.seq \
	sequences/t1_icbm.seq \
	sequences/t2_ai.seq \
	sequences/t2_icbm_avg.seq \
	sequences/t2_icbm.seq

noinst_HEADERS = \
	src/minc/chirp.h \
	src/minc/fourn.h \
	src/minc/imincfile.h \
	src/minc/iomincfile.h \
	src/minc/mincfile.h \
	src/minc/mincicv.h \
	src/minc/mriimage.h \
	src/minc/mrilabel.h \
	src/minc/mrimatrix.h \
	src/minc/mriminc.h \
	src/minc/mristring.h \
	src/minc/mrivolume.h \
	src/minc/omincfile.h \
	src/minc/time_stamp.h \
	src/mrisim/discrete_label_phantom.h \
	src/mrisim/discrete_phantom.h \
	src/mrisim/discrete_rf_phantom.h \
	src/mrisim/fuzzy_label_phantom.h \
	src/mrisim/fuzzy_phantom.h \
	src/mrisim/fuzzy_rf_phantom.h \
	src/mrisim/image_snr_coil.h \
	src/mrisim/intrinsic_coil.h \
	src/mrisim/mriscanner.h \
	src/mrisim/mrisimargs.h \
	src/mrisim/mrisim.h \
	src/mrisim/mrisim_main.h \
	src/mrisim/paramfile.h \
	src/mrisim/ParseArgv.h \
	src/mrisim/percent_coil.h \
	src/mrisim/phantom.h \
	src/mrisim/rf_coil.h \
	src/mrisim/rf_tissue_phantom.h \
	src/mrisim/scanner_output.h \
	src/mrisim/tissue_phantom.h \
	src/signal/ce_fast.h \
	src/signal/customseq.h \
	src/signal/event.h \
	src/signal/fast_iso_model.h \
	src/signal/ffe.h \
	src/signal/fisp.h \
	src/signal/flash.h \
	src/signal/ir.h \
	src/signal/isochromat_model.h \
	src/signal/pulseseq.h \
	src/signal/quick_model.h \
	src/signal/quickseq.h \
	src/signal/repeat.h \
	src/signal/rf_pulse.h \
	src/signal/sample.h \
	src/signal/se.h \
	src/signal/signal.h \
	src/signal/spin_model.h \
	src/signal/spoiled_flash.h \
	src/signal/spoiler.h \
	src/signal/tissue.h \
	src/signal/vector.h \
	src/signal/vector_model.h



# Programs we build and install.
#
bin_PROGRAMS = \
	mrisim

mrisim_SOURCES = \
	src/minc/chirp.cxx \
	src/minc/fourn.c \
	src/minc/imincfile.cxx \
	src/minc/iomincfile.cxx \
	src/minc/mincfile.cxx \
	src/minc/mincicv.cxx \
	src/minc/mriimage.cxx \
	src/minc/mrilabel.cxx \
	src/minc/mrimatrix.cxx \
	src/minc/mristring.cxx \
	src/minc/mrivolume.cxx \
	src/minc/omincfile.cxx \
	src/minc/time_stamp.c \
	src/mrisim/discrete_label_phantom.cxx \
	src/mrisim/discrete_phantom.cxx \
	src/mrisim/discrete_rf_phantom.cxx \
	src/mrisim/fuzzy_label_phantom.cxx \
	src/mrisim/fuzzy_phantom.cxx \
	src/mrisim/fuzzy_rf_phantom.cxx \
	src/mrisim/image_snr_coil.cxx \
	src/mrisim/intrinsic_coil.cxx \
	src/mrisim/mriscanner.cxx \
	src/mrisim/mrisimargs.cxx \
	src/mrisim/mrisim_main.cxx \
	src/mrisim/paramfile.cxx \
	src/mrisim/ParseArgv.c \
	src/mrisim/percent_coil.cxx \
	src/mrisim/phantom.cxx \
	src/mrisim/rf_coil.cxx \
	src/mrisim/rf_tissue_phantom.cxx \
	src/mrisim/scanner_output.cxx \
	src/mrisim/tissue_phantom.cxx \
	src/signal/ce_fast.cxx \
	src/signal/customseq.cxx \
	src/signal/event.cxx \
	src/signal/fast_iso_model.cxx \
	src/signal/ffe.cxx \
	src/signal/fisp.cxx \
	src/signal/flash.cxx \
	src/signal/ir.cxx \
	src/signal/isochromat_model.cxx \
	src/signal/pulseseq.cxx \
	src/signal/quick_model.cxx \
	src/signal/quickseq.cxx \
	src/signal/repeat.cxx \
	src/signal/rf_pulse.cxx \
	src/signal/sample.cxx \
	src/signal/se.cxx \
	src/signal/spin_model.cxx \
	src/signal/spoiled_flash.cxx \
	src/signal/spoiler.cxx \
	src/signal/tissue.cxx \
	src/signal/vector.cxx \
	src/signal/vector_model.cxx
