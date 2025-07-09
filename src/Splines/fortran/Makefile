# === CONFIGURABLE ===
FC      = gfortran                 
FFLAGS  ?= -fPIC                    
RECURS  ?= -frecursive
LDFLAGS ?= -shared     
DEPSDIR = ../../../deps

# 1.  Simple OS/CPU detection (works on MSYS/MinGW, Linux, macOS)

OS  ?= $(OS)
CPU ?= $(PROCESSOR_ARCHITECTURE)

ifeq ($(strip $(OS)),)             
  OS  := $(shell uname -s 2>/dev/null)
endif
ifeq ($(strip $(CPU)),)
  CPU := $(shell uname -p 2>/dev/null)
endif

ifneq ($(filter MINGW%,$(OS)),)
  OS := Windows_NT
endif



ifeq ($(OS),Windows_NT)                        # -------- Windows ----------
  SOEXT   := dll
  FFLAGS  += -g -O0 -fcheck=all -fbacktrace
  LIBS    ?= -llapack -lblas                  
else ifeq ($(OS),Darwin)                       # -------- macOS ------------
  SOEXT   := dylib
  LIBS    ?= -framework Accelerate
else                                           # -------- Linux -----------
  SOEXT   := so
  LIBS    ?= -llapack -lblas
endif


# 2.  Derived variables & pattern rules
OUTLIB  := libspline.$(SOEXT)
F90     := $(FC) $(FFLAGS) $(RECURS)

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
all: $(OUTLIB)

$(OUTLIB): $(OBJS)
	$(FC) $(LDFLAGS) -o $@ $^ $(LIBS)
	@mkdir -p $(DEPSDIR)
	mv -f $(OUTLIB)  $(DEPSDIR)

clean:
	rm -f *.o *.mod *.original *.dylib *.dll *.so

# 5. dependencies
spline.o: defs.o
cspline.o: defs.o spline.o
bicube.o: defs.o spline.o
fft.o: defs.o
fspline.o: defs.o spline.o cspline.o fft.o
spline_c_api.o: defs.o spline.o cspline.o bicube.o
