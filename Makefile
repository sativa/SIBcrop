# Makfile for SiB3crop
#
# You will need GNU's version of Make, which is the default on Linux, OS X
# but not necessarily other Unixes.
#
# If you want to change the compiler or optimization without editing the
# makefile, set the COMPILER or OPT environment variables on the command
# line like so:
#
# user@host:~$ COMPILER=ifort OPT=opt make clean all
#
# Values for COMPILER: gcc*, pgi, ifort
# Values for OPT: debug*, opt
#
# * = default
#
# To name the SiB binary with the platform, define the FULL_SUFFIX environment
# like this before running make:
#
# user@host:~$ export FULL_SUFFIX=1                        # bourne shells
#
# user@host:~$ setenv FULL_SUFFIX 1                        # c shell
#
# You can tell the Makefile where to find NetCDF and LAPACK/BLAS libraries
# by setting the NETCDF_ROOT and LAPACK_ROOT environment variables as well:
#
# user@host:~$ export NETCDF_ROOT=/usr/local
#
# the includes will be searched for in $(NETCDF_ROOT)/include and
# $(LAPACK_ROOT)/include, and the libraries will be searched for in
# $(NETCDF_ROOT)/lib and $(LAPACK_ROOT)/lib.

OS = $(shell uname -s)
PROC = $(shell uname -p)

# Uncomment one of these to select your compiler. Use gcc for gfortran
# and pgi for the portland group compilers.
ifndef COMPILER
  COMPILER = gcc
  #COMPILER = pgi
  #COMPILER = ifort
endif

# Uncomment one of these to select optimize vs. debug compile modes.
ifndef OPT
  #OPT = opt
  OPT = debug
endif

# Non-System Headers and Libraries
# -------------------------------------------------------------------------
# BLAS (Basic Linear Algebra System) and LAPACK (Linear Algebra Package)

# Directory containing netcdf's include/ and lib/ directories.
ifndef NETCDF_ROOT
  NETCDF_ROOT = /usr/local/netcdf3-$(COMPILER)
  #NETCDF_ROOT = /usr/local/
endif

ifndef LAPACK_ROOT
	LAPACK_ROOT = /usr/local/atlas
endif

ifeq ($(OS),Darwin)
  ifdef FULL_SUFFIX
    SUFFIX	= -mac$(COMPILER)-$(OPT)
  endif
  INCLUDES 	= -I$(NETCDF_ROOT)/include
  LALIB    	= -framework vecLib
  LIBS  	= $(LALIB) -L$(NETCDF_ROOT)/lib -lnetcdf
endif

ifeq ($(OS),Linux)
  ifdef FULL_SUFFIX
    SUFFIX	= -lx$(COMPILER)-$(OPT)
  endif
  LALIB 	= -L$(LAPACK_ROOT)/lib -llapack -lblas -lcblas -latlas
  INCLUDES 	= -I$(NETCDF_ROOT)/include
  LIBS     	= $(LALIB) -L$(NETCDF_ROOT)/lib -lnetcdf
endif

ifeq ($(COMPILER),gcc)
  F90 = gfortran
  ifeq ($(OPT),opt)
    F90FLAGS = -O2
    F90FLAGS_NOOPT =
  else ifeq ($(OPT),debug)
    F90FLAGS = -g -fbounds-check -ffpe-trap=invalid,zero,overflow
    F90FLAGS_NOOPT = $(F90FLAGS)
  endif
  ifeq ($(PROC),powerpc)
    F90FLAGS += -fconvert=little-endian
  endif
  F90FLAGS += -DPGF -fimplicit-none -Wsurprising $(INCLUDES)
  F90FLAGS_NOOPT += -DPGF -fimplicit-none -Wsurprising $(INCLUDES)
  LFLAGS = $(LIBS)
endif

ifeq ($(COMPILER),pgi)
  F90 = pgf90
  ifeq ($(OPT),opt)
    F90FLAGS = -fast -Mnoframe
    F90FLAGS_NOOPT =
  else ifeq ($(OPT),debug)
    F90FLAGS = -g -Mbounds -Ktrap=fp
    F90FLAGS_NOOPT = $(F90FLAGS)
  endif
  F90FLAGS += -DPGF=1 -Minfo=loop,inline -Minform=inform $(INCLUDES)
  F90FLAGS_NOOPT += -DPGF=1 -Minfo=loop,inline -Minform=inform $(INCLUDES)
  LFLAGS   = -v -Minform=inform $(LIBS)
endif

ifeq ($(COMPILER),ifort)
  F90 = ifort
  ifeq ($(OPT),opt)
    F90FLAGS = -O2
    F90FLAGS_NOOPT = -O0
  else ifeq ($(OPT),debug)
    F90FLAGS = -g -check bounds -fpe0
    F90FLAGS_NOOPT = $(F90FLAGS)
  endif
  ifeq ($(PROC),powerpc)
    F90FLAGS += -convert little_endian
  endif
  F90FLAGS += -DPGF -implicitnone -warn all $(INCLUDES)
  F90FLAGS_NOOPT += -DPGF -implicitnone -warn all $(INCLUDES)
  LFLAGS = $(LIBS)
endif

# Objects (DO NOT EDIT - Unless You Add or Remove Files)
# -------------------------------------------------------------------------
# Variable objects

VAR_OBJS  = kinds.o \
	    physical_parameters.o \
	    eau_params.o  \
        sib_const_module.o \
        sib_io_module.o \
	    sib_bc_module.o \
        sibtype.o \
	    timetype.o

# Pre-requisite objects (Required for the Scientific objects)
PRE_OBJS  = 

# Scientific objects
SCI_OBJS  = addinc.o \
	    balan.o \
	    dtess_eau.o \
	    begtem.o \
	    cfrax.o \
	    clm_combo.o \
	    combine_snow.o \
	    compact_snow.o \
	    cycalc.o \
	    delef.o \
	    delhf.o \
	    dellwf.o \
	    hydro_canopy.o \
	    hydro_snow.o \
	    hydro_soil.o \
	    qsat_eau.o \
	    ess_eau.o \
	    netrad.o \
	    phosib.o \
	    rada2.o \
	    rbrd.o \
	    respsib.o \
	    rnload.o \
	    soilwater.o \
	    sortin.o \
	    subdivide_snow.o \
	    tridiag_solver.o \
	    vmfcalz.o \
	    vmfcalzo.o \
	    vntlat.o \
	    sibslv.o \
	    update.o \
	    resp_control.o \
	    sib.o # sib.o NEEDS TO BE LAST

# Netcdf Objects
NCDF_OBJS = nc_util.o

# SiB Drive objects
DRV_OBJS  = init_grid.o \
	    init_var.o\
	    zenith.o \
	    phen_mapper.o \
	    crop_accum.o \
	    leaf_weight.o \
	    init_sibdrv.o \
	    init_crop.o \
        time_init.o \
        time_manager.o \
        sibdrv_read_single.o \
	    sibdrv_read_ecmwf.o \
        sibdrv_read_ncep2.o \
        sibdrv_read_geos4.o \
	sibdrv_read_narr.o \
	    read_ti.o \
	    set_ti.o \
	    mapper.o \
	    phen_mapper.o \
	    calculate_td.o \
	    read_ndvi.o \
	    previous_bc.o \
	    update_bc.o \
        rtape_sib.o \
        diagnostic_output.o \
        pbpwrite.o \
        qpwrite.o \
        output_control.o \
        sibdrv_interp.o \
        bc_interp.o \
	    SiBDRV.o

# Compilation (DO NOT EDIT BELOW THIS LINE)
# -------------------------------------------------------------------------
%.o: %.f90
	$(F90) $(F90FLAGS) -c $< -o $@

%.o: %.F90
	$(F90) $(F90FLAGS) -c $< -o $@

all: SiBD3crop$(SUFFIX) sibmerge$(SUFFIX)

remake: clean all

SiBD3crop$(SUFFIX): $(VAR_OBJS) $(PRE_OBJS) $(SCI_OBJS) $(NCDF_OBJS) $(DRV_OBJS)
	$(F90) $^ -o $@ $(LFLAGS)
	@echo -e "\n"
	@date

sibmerge$(SUFFIX): sibmerge.o $(NCDF_OBJS)
	$(F90) sibmerge.o $(NCDF_OBJS) -o $@ $(LFLAGS)

crop_accum.o: crop_accum.F90
	$(F90) $(F90FLAGS_NOOPT) -c $< -o $@

leaf_weight.o: leaf_weight.F90
	$(F90) $(F90FLAGS_NOOPT) -c $< -o $@

clean:
	rm -f SiBD3crop$(SUFFIX) sibmerge$(SUFFIX) *.o *.mod *.stb *~

distclean: clean
	rm -rf SiBD3crop-* sibmerge-* #*#

help:
	@echo ""
	@echo "Supported Rules:"
	@echo "1) [nothing], SiBD3 and sibmerge: compile sib3 model, " \
	"recompiling dependencies as necessary.  Then, compile sibmerge program"
	@echo "2) SiBD3: compile sib3 model only, recompiling dependencies " \
	"as necessary."
	@echo "3) sibmerge: compile sibmerge program"
	@echo "4) <object name>.o: compile specific object. matches the file" \
	"name."
	@echo "5) remake: delete all objects, and re-compile the sib3 model " \
	"and sibmerge."
	@echo "6) clean: remove all compiler-generated files and executable."
	@echo "7) help: print this message."
	@echo ""

# Dependencies

$(VAR_OBJS):

$(PRE_OBJS): $(VAR_OBJS)

$(SCI_OBJS): $(VAR_OBJS) $(PRE_OBJS)

$(NCDF_OBJS): $(VAR_OBJS)

$(DRV_OBJS): $(VAR_OBJS) $(NCDF_OBJS)
