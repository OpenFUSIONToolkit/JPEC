#include ../install/DEFAULTS.inc

# === CONFIGURABLE ===
FC = gfortran  # You can override this: `make FC=ifort`
FFLAGS ?= -fPIC -g -O0 -fcheck=all -fbacktrace
LDFLAGS ?= -shared
LIBS ?= -llapack -lblas
OUTLIB = libspline.dll
RECURSFLAG ?= -frecursive

F90 = $(FC) $(FFLAGS) $(RECURSFLAG)

.f.o:
	$(F90) -c $< -o $@


OBJS = \
	defs.o \
	spline.o \
	cspline.o \
	bicube.o \
	fft.o \
	fspline.o \
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
fft.o: defs.o
fspline.o: defs.o spline.o cspline.o fft.o
spline_c_api.o: defs.o spline.o cspline.o bicube.o
