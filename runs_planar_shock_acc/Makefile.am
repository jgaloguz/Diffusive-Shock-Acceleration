## Process this file with automake to produce Makefile.in

AM_CXXFLAGS = $(MPI_CFLAGS)
AM_LDFLAGS = $(MPI_LIBS)

bin_PROGRAMS = dsa_forward_acc dsa_likelihood_expectation

SPBL_COMMON_DIR = ../common
SPBL_SOURCE_DIR = ../src

dsa_forward_acc_SOURCES = dsa_forward_acc.cc dsa_common.hh \
   $(SPBL_SOURCE_DIR)/background_smooth_shock.cc \
   $(SPBL_SOURCE_DIR)/background_smooth_shock.hh \
   $(SPBL_SOURCE_DIR)/background_shock.cc \
   $(SPBL_SOURCE_DIR)/background_shock.hh \
   $(SPBL_SOURCE_DIR)/background_base.cc \
   $(SPBL_SOURCE_DIR)/background_base.hh \
   $(SPBL_SOURCE_DIR)/diffusion_other.cc \
   $(SPBL_SOURCE_DIR)/diffusion_other.hh \
   $(SPBL_SOURCE_DIR)/diffusion_base.cc \
   $(SPBL_SOURCE_DIR)/diffusion_base.hh \
   $(SPBL_COMMON_DIR)/matrix.cc \
   $(SPBL_COMMON_DIR)/matrix.hh \
   $(SPBL_COMMON_DIR)/vectors.cc \
   $(SPBL_COMMON_DIR)/vectors.hh \
   $(SPBL_COMMON_DIR)/multi_index.hh \
   $(SPBL_COMMON_DIR)/simple_array.hh \
   $(SPBL_COMMON_DIR)/arithmetic.hh \
   $(SPBL_COMMON_DIR)/params.cc \
   $(SPBL_COMMON_DIR)/params.hh \
   $(SPBL_COMMON_DIR)/physics.cc \
   $(SPBL_COMMON_DIR)/physics.hh \
   $(SPBL_COMMON_DIR)/data_container.hh \
   $(SPBL_COMMON_DIR)/data_container.cc \
   $(SPBL_COMMON_DIR)/spatial_data.hh \
   $(SPBL_COMMON_DIR)/definitions.hh \
   $(SPBL_COMMON_DIR)/random.hh \
   $(SPBL_COMMON_DIR)/print_warn.hh

dsa_forward_acc_LDADD = $(MPI_LIBS) $(GSL_LIBS)

dsa_likelihood_expectation_SOURCES = dsa_likelihood_expectation.cc dsa_common.hh \
   $(SPBL_SOURCE_DIR)/background_smooth_shock.cc \
   $(SPBL_SOURCE_DIR)/background_smooth_shock.hh \
   $(SPBL_SOURCE_DIR)/background_shock.cc \
   $(SPBL_SOURCE_DIR)/background_shock.hh \
   $(SPBL_SOURCE_DIR)/background_base.cc \
   $(SPBL_SOURCE_DIR)/background_base.hh \
   $(SPBL_SOURCE_DIR)/diffusion_other.cc \
   $(SPBL_SOURCE_DIR)/diffusion_other.hh \
   $(SPBL_SOURCE_DIR)/diffusion_base.cc \
   $(SPBL_SOURCE_DIR)/diffusion_base.hh \
   $(SPBL_COMMON_DIR)/matrix.cc \
   $(SPBL_COMMON_DIR)/matrix.hh \
   $(SPBL_COMMON_DIR)/vectors.cc \
   $(SPBL_COMMON_DIR)/vectors.hh \
   $(SPBL_COMMON_DIR)/multi_index.hh \
   $(SPBL_COMMON_DIR)/simple_array.hh \
   $(SPBL_COMMON_DIR)/arithmetic.hh \
   $(SPBL_COMMON_DIR)/params.cc \
   $(SPBL_COMMON_DIR)/params.hh \
   $(SPBL_COMMON_DIR)/physics.cc \
   $(SPBL_COMMON_DIR)/physics.hh \
   $(SPBL_COMMON_DIR)/data_container.hh \
   $(SPBL_COMMON_DIR)/data_container.cc \
   $(SPBL_COMMON_DIR)/spatial_data.hh \
   $(SPBL_COMMON_DIR)/definitions.hh \
   $(SPBL_COMMON_DIR)/random.hh \
   $(SPBL_COMMON_DIR)/print_warn.hh

dsa_likelihood_expectation_LDADD = $(MPI_LIBS) $(GSL_LIBS)