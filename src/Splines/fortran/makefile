FC ?= gfortran  # You can override this: `make FC=ifort`
FFLAGS ?= -fPIC
LDFLAGS ?= -shared
LIBS ?= -framework Accelerate
RECURSFLAG ?= -frecursive
LIBSUFFIX ?= .dylib

OUTLIB = libspline$(LIBSUFFIX)

F90 = $(FC) $(FFLAGS) $(RECURSFLAG)

.f.o:
	$(F90) -c $< -o $@

# 3.  Object list
OBJS = \
	defs.o \
	spline.o \
	cspline.o \
	bicube.o \
	fft.o \
	fspline.o \
	spline_c_api.o

# 4.  Targets

all : $(OUTLIB) 

$(OUTLIB): $(OBJS)
	$(FC) $(LDFLAGS) -o $@ $^ $(LIBS)
	mv $(OUTLIB) ../../../deps/

clean:
	rm -f *.o *.mod *.original *.dylib
	
# 5. dependencies
spline.o: defs.o
cspline.o: defs.o spline.o
bicube.o: defs.o spline.o
fft.o: defs.o
fspline.o: defs.o spline.o cspline.o fft.o
spline_c_api.o: defs.o spline.o cspline.o bicube.o
