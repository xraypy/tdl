##############################################################
# NOTE the following environment variable must be set:
#   - COMPILER: Which compiler to be used
#   - COMPFLAGS: Compiler flags
#   - TDLINCLUDE: Include statement for compiler
#
# Execute this file as:
# >>make build    # to build the obj code
# >>make rebuild  # to rebuild the obj code
# >>make clean    # to delete obj code
##############################################################
CC       = $(COMPILER)
CFLAGS   = $(COMPFLAGS)
INCLUDE  = $(TDLINCLUDE)

SRCFILES = utils.c numfcns.c
##############################################################

build: $(SRCFILES)
	$(CC) $(CFLAGS) $(INCLUDE) -c $(SRCFILES)
	cp *.o ../../build/.

clean:
	rm -f *.o

rebuild: clean build
