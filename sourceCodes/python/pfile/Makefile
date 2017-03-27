# Makefile for generating Python PFile library
#
# Copyright 2011 <fill in later>

# There are two kinds of active 64-bit machines in the cluster: as41-amd64
# and as50-amd64. The former are Nikola 4 boxes, and only have Python 2.5
# binaries; the latter are Nikola 5 boxes, and only have Python 2.6 binaries.
# To get a library that works on both types of machines, run this makefile on
# a as41-amd64 machine.
#
# To find the as41-amd64 machines in the cluster, run the following command:
# > condor_status -const 'NikolaNarch == "as41-amd64"' | less

CXX = g++
SWIG = $(HOME)/swig/bin/swig
PYVER = $(shell ./pyver.py)
NIKVER = $(shell nikver -n)
PYINCLUDE = /usr/nikola/pkgs/python/include/python$(PYVER)

objects = error.o general.o rand.o pfile.o vbyteswapping.o libpfile_wrap.o

_libpfile.so: $(objects)
	@echo "INFO: Building on Nikola $(NIKVER)"
ifeq "$(NIKVER)" "as41-amd64"
	@echo "INFO: Library will work on Nikola 4 and Nikola 5 machines."
else
	@echo "INFO: Library will *not* work on Nikola 4 machines."
endif
	$(CXX) -shared -Wl $(objects) -o _libpfile.so

.cc.o:
	$(CXX) -DHAVE_FSEEKO -c -fPIC $<

libpfile_wrap.o: libpfile_wrap.cxx
	$(CXX) -DHAVE_FSEEKO -c -fPIC -I $(PYINCLUDE) libpfile_wrap.cxx

libpfile_wrap.cxx: libpfile.i
	$(SWIG) -c++ -python libpfile.i

clean:
	-rm -f $(objects) libpfile.pyc _libpfile.so *~

distclean: clean
	-rm -f libpfile_wrap.cxx libpfile.py
