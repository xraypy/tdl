###############################################################################
# Run make in the each of the src directories
#
# Note this requires the environment variable TDLDIR be defined
#   - TDLDIR = /path/to/tdl
#
# Also assume that header files for gsl are located under /usr/local/
#
# Execute this file as:
# >>make build    # to build
# >>make rebuild  # to rebuild (clean then build)
# >>make clean    # to clean
# >>make          # default is to build
###############################################################################
TDLLIBDIR  = $(TDLDIR)/lib
TDLSRCDIR  = $(TDLLIBDIR)/src

export COMPILER   = gcc
export COMPFLAGS  = -Wall -O
export TDLINCLUDE = -I $(TDLSRCDIR) -I /usr/local/ -I /usr/include

AR      = ar -rv
CCS     = $(COMPILER) -shared
GSLLIB  = /usr/lib/libgsl.a /usr/lib/libgslcblas.a

LIBOBJ  = utils.o numfcns.o lmdif.o fit.o peak.o xrr.o
WRAPOBJ = wrap_xrr.o wrap_hello.o
OBJ     = $(LIBOBJ) $(WRAPOBJ)
LIBTDL  = libtdl.a
###############################################################################

# default built all dll's
build: _xrr.so _hello.so

# clean then build
rebuild: clean build

# generate _xrr.dll from wrap_xrr.o
_xrr.so: $(LIBTDL)
	$(CCS) -o $@  wrap_xrr.o $(LIBTDL)  $(GSLLIB)
	cp $@  $(TDLLIBDIR)/.

# generate _hello.dll from wrap_hello.o
_hello.so: $(LIBTDL)
	$(CCS) -o $@ wrap_hello.o $(LIBTDL) 
	cp $@  $(TDLLIBDIR)/.

# generate LIBTDL from OBJ
$(LIBTDL): $(LIBOBJ)
	$(AR) $(LIBTDL) $(OBJ)

$(OBJ):
	$(MAKE) -C  $(TDLSRCDIR)/fit 
	$(MAKE) -C  $(TDLSRCDIR)/peak
	$(MAKE) -C  $(TDLSRCDIR)/utils 
	$(MAKE) -C  $(TDLSRCDIR)/wrap
	$(MAKE) -C  $(TDLSRCDIR)/xrr 
	cd $(MAKEDIR)

#clean
clean:
	$(MAKE) -C  $(TDLSRCDIR)/fit    clean
	$(MAKE) -C  $(TDLSRCDIR)/peak   clean
	$(MAKE) -C  $(TDLSRCDIR)/utils  clean
	$(MAKE) -C  $(TDLSRCDIR)/wrap   clean
	$(MAKE) -C  $(TDLSRCDIR)/xrr    clean
	cd $(MAKEDIR)
	rm -f *.o; rm -f *.obj; 
	rm -f *.so; rm -f *.dll
	rm -f *.a; rm -f *.lib

