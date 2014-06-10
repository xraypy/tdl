#
# Makefile for XRR C library
#

include Mk.config
include Mk.install

INSTALLNAME = xrr
STATICLIB   = $(DLLPRE)$(INSTALLNAME).a
SHAREDLIB   = $(DLLPRE)$(INSTALLNAME).$(DLLEXT)

.f.o:
	$(F77) -c $*.f

.PHONY: clean install default

default: $(SHAREDLIB)

LINK_FFLAGS  =
LINK         = gcc -O2 $(ARCH_TARGS) -fPIC -shared
XRRHEADERS = xrr.h wrap.h numfcns.h utils.h complex.h mconf.h protos.h

OBJ  = xrr.o  numfcns.o utils.o clog.o cmplx.o

SHARED = $(CCOPTS) -shared
LINKXRR = -L$(libdir) -L. -lxrr

xrr.o: xrr.h utils.h numfcns.h complex.h
clog.o: complex.h mconf.h

$(SHAREDLIB): $(OBJ)
	@rm -f $@
	$(LINK) -o $(SHAREDLIB) $(OBJ) $(LINK_FFLAGS) $(LINK_CFLAGS)

$(STATICLIB): $(OBJ)
	$(AR) cru $(STATICLIB) $(OBJ)
	$(RANLIB) $(STATICLIB)

clean :
	rm -f $(OBJ)  $(STATICLIB) $(SHAREDLIB)

install: $(SHAREDLIB)
	mkdir -p $(INSTALLDIR)
	cp $(SHAREDLIB) $(INSTALLDIR)/.
