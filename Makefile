# Start oF the makeFile DeFining variables
lib_f := $(wildcard Library/*.f)
#lib_f90 := $(wildcard Library/*.f90)
lib_f90 := Library/timing.o \
           Library/general_utility.o\
	   Library/precin.o \
	   Library/Gauss_Quadrature.o\
	   Library/Lagrange.o\
	   Library/Lagrange_weights.o\
	   Library/hermite.o\
	   
lib_objects := $(lib_f90:.f90=.o) $(lib_f:.f=.o) 

objects := \
	$(lib_objects) \
	fconfig.o \
	parameters.o \
	grid.o \
	potential.o \
	banded_matrices.o\
	pulse.o \
	propagator.o\
	integral_method.o\
	examples.o\
  main.o 


F90comp:=gfortran
NOLINK = -c
#LDFLAGS = -L/c:/lib/ -lblas -llapack

FCFLAGS  = -g -fbacktrace

######################################
# SET COMPILER FLAGS DEPENDING ON WHICH ONE IS BEING USED
######################################
ifeq  ($(F90comp),gfortran)
LDFLAGS = -L/usr/lib/ -lblas -llapack
else
ifeq  ($(F90comp),ifort)
LDFLAGS = -mkl#sequential
else
echo "please provide the lapack link of libraries"
endif
endif

## MakeFile
run: $(objects)
	$(F90comp) $(OMPFLAGS) $(FCFLAGS) -o run_volterra.exe $(objects) $(LDFLAGS)
%: %.o
	$(F90comp) $(FCFLAGS) -o $@ $^ $(LDFLAGS)


%.o: %.f90
	$(F90comp) -o $@ $(FCFLAGS) $(OMPFLAGS) -c $<
%.o: %.F90
	$(F90comp) -o $@ $(FCFLAGS) $(OMPFLAGS) -c $<
%.o: %.f
	$(F90comp) -o $@ $(FCFLAGS) $(OMPFLAGS) -c $<


# Cleaning everything
clean:
	rm -fR *.x
	rm -f $(objects)
	rm -f $(objects:.o=.mod)