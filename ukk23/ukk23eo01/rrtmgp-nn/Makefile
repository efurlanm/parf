#!/bin/make
# RTE_DIR = rte
RRTMGP_DIR = rrtmgp

#
# Compiler variables FC, FCFLAGS can be set in the environment or in Makefile.conf
#
# -include Makefile.conf
#
# Choose kernels depending on environment variable
#

ifeq ($(USE_OPENACC),1)
  # RTE_KERNEL_DIR    = rte/kernels-openacc
  RRTMGP_KERNEL_DIR = rrtmgp/kernels-openacc
  FCFLAGS += -DUSE_OPENACC
endif

# RTE_KERNEL_DIR += rte/kernels
RRTMGP_KERNEL_DIR += rrtmgp/kernels
NEURAL_DIR = neural
VPATH = $(RRTMGP_DIR):$(RRTMGP_KERNEL_DIR):$(NEURAL_DIR)

LIBRRTMGP = ../lib/librrtmgp.a ../lib/libneural.a
all: $(LIBRRTMGP)

COMPILE = $(FC) $(FCFLAGS) $(FCINCLUDE) -c
%.o: %.F90
	$(COMPILE) $<

# include $(RTE_DIR)/Make.depends
include $(RRTMGP_DIR)/Make.depends
include $(NEURAL_DIR)/Make.depends

# ../lib/librte.a: $(RTE_SRC)
# 	ar r ../lib/librte.a $(RTE_SRC)

../lib/librrtmgp.a: $(RRTMGP_SRC)
	ar r ../lib/librrtmgp.a $(RRTMGP_SRC)

../lib/libneural.a: $(NEURAL_SRC)
	ar r ../lib/libneural.a $(NEURAL_SRC)

clean:
	rm -f *.optrpt *.mod *.o $(LIBRRTMGP)

