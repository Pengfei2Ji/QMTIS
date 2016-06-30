#
# Top Makefile for linear algebra libraries
#

include config.mk

INSTALL = ./install-sh

all_targets all:
	cd blas && $(MAKE)
	cd lapack && $(MAKE)

install:
	$(INSTALL) -d -m 755 $(prefix)/lib
	$(INSTALL) -m 644 blas/libblas.a $(prefix)/lib
	$(INSTALL) -m 644 lapack/liblapack.a $(prefix)/lib

clean:
	cd blas && $(MAKE) clean
	cd lapack && $(MAKE) clean

