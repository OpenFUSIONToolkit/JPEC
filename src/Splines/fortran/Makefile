#include ../install/DEFAULTS.inc

# === CONFIGURABLE ===
FC = gfortran  # You can override this: `make FC=ifort`
FFLAGS ?= -fPIC
LDFLAGS ?= -shared
LIBS ?= -framework Accelerate
OUTLIB = libspline.dylib
RECURSFLAG ?= -frecursive

F90 = $(FC) $(FFLAGS) $(RECURSFLAG)

.f.o:
	$(F90) -c $< -o $@


OBJS = \
	defs.o \
	spline.o \
	cspline.o \
	bicube.o \
	spline_c_api.o

# targets

all : $(OUTLIB) 

$(OUTLIB): $(OBJS)
	$(FC) $(LDFLAGS) -o $@ $^ $(LIBS)
	mv $(OUTLIB) ../../../deps/

clean:
	rm -f *.o *.mod *.original *.dylib

# dependencies

spline.o: defs.o
cspline.o: defs.o spline.o
bicube.o: defs.o spline.o
spline_c_api.o: defs.o spline.o cspline.o bicube.o
