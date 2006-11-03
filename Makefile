# Makfile for sib3
# All options are set at the top of this file

# Non-System Headers and Libraries 
# -------------------------------------------------------------------------
# BLAS (Basic Linear Algebra System) and LAPACK (Linear Algebra Package)
#LALIB     = /usr/local/pgi/linux86/lib
# BLAS-(64-BIT)
LALIB     = /opt/pgi-6.1-4/linux86-64/6.1/lib

# Netcdf
NETCDFINC = /usr/local/pgilibs/include
NETCDFLIB = /usr/local/pgilibs/lib


# Includes and Links
# -------------------------------------------------------------------------
INCLUDES  = -I$(NETCDFINC)
LINKS     = -L$(NETCDFLIB) -lnetcdf -lnetcdff -L$(LALIB) -llapack -lblas -lpgmp -lpthread
#LINKS     = -L$(NETCDFLIB) -lnetcdf -L$(LALIB) -llapack -lblas -lg2c -lpthread

#PROCESSOR INFO: determine your processor with the command
#  more /proc/cpuinfo


# Fortran90 Compiler
# -------------------------------------------------------------------------
# Portland Group
F90       = pgf90

#flags for running on nodes 1-9 of the cluster
F90FLAGS  =  -Mnoframe -DPGF=1 -tp p5 -Ktrap=fp -Minfo=loop,inline -Minform=inform $(INCLUDES)
#flags for running on nodes 20-29 of the cluster
#F90FLAGS  =  -Mnoframe -DPGF=1 -tp p7 -Ktrap=fp -Minfo=loop,inline -Minform=inform $(INCLUDES)

# WithDebugging and Profiling
#insert the code for your local processor instead of athlonxp
#F90FLAGS  = -g -DPGF=1 -Mnoframe -tp athlonxp -Ktrap=fp -Minfo=loop,inline -Minform,inform $(INCLUDES)

#Opteron 
#F90FLAGS  = -DPGF=1 -Mnoframe -tp k8-64 -Ktrap=fp -Minfo=loop,inline -Minform=inform $(INCLUDES)
#F90FLAGS  =    -DPGF=1 -Mnoframe -tp k8-64 -Ktrap=fp -Minfo=loop,inline -Minform=inform $(INCLUDES)

# Intel 8.0 compiler
#F90       = ifort
#F90FLAGS  = -g -DPGF=1 -fast -fpe0 $(INCLUDES)

# F90Linker Arguments
# -------------------------------------------------------------------------
# Portland Group
LNKFLAGS  = -v -fast -Minform=inform $(LINKS)
# With Debugging and Profiling
#LNKFLAGS  = -g -v -Minform,inform $(LINKS)
# Intel 8.0 compiler
#LNKFLAGS  = -g -fast $(LINKS)


# Objects (DO NOT EDIT - Unless You Add or Remove Files)
# -------------------------------------------------------------------------
# Variable objectsr

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
NCDF_OBJS = handle_err.o 

# SiB Drive objects
DRV_OBJS  = init_grid.o \
	    init_var.o\
	    zenith.o \
	    crop_accum.o \
	    init_sibdrv.o \
            time_init.o \
            time_manager.o \
            sibdrv_read_single.o \
	    sibdrv_read_ecmwf.o \
            sibdrv_read_ncep2.o \
            sibdrv_read_geos4.o \
	    read_ti.o \
	    mapper.o \
	    calculate_td.o \
	    read_ndvi.o \
	    previous_bc.o \
            rtape_sib.o \
            diagnostic_output.o \
            pbpwrite.o \
            qpwrite.o \
            output_control.o \
            sibdrv_interp.o \
            bc_interp.o \
	    SiBDRV.o

#sibmerge objects
SIBMERGE_OBJS = sibmerge.o

# Compilation (DO NOT EDIT BELOW THIS LINE)
# -------------------------------------------------------------------------
all:	SiBcrop sibmerge

SiBcrop: $(VAR_OBJS) $(PRE_OBJS) $(SCI_OBJS) $(NCDF_OBJS) $(DRV_OBJS)
	$(F90) $^  $(LNKFLAGS) -o  $@
	@echo -e "\n"
	@date

sibmerge: $(SIBMERGE_OBJS)
	$(F90) $^  $(LNKFLAGS) -o  $@ 
	@echo -e "\n"
	@date

$(VAR_OBJS):

$(PRE_OBJS): $(VAR_OBJS)

$(SCI_OBJS): $(VAR_OBJS) $(PRE_OBJS)

$(NCDF_OBJS): $(VAR_OBJS)

$(DRV_OBJS): $(VAR_OBJS) $(NCDF_OBJS)

%.o: %.f90
	$(F90) -c $(F90FLAGS)  $< -o $@
	@echo ""
%.o: %.F90
	$(F90) -c $(F90FLAGS) $< -o $@
	@echo ""
remake: clean all

clean:
	rm -f SiBcrop sibmerge *.o *.mod *.stb

help:
	@echo ""
	@echo "Supported Rules:"
	@echo "1) [nothing], SiBcrop and sibmerge: compile sib3 model, " \
	"recompiling dependencies as necessary.  Then, compile sibmerge program"
	@echo "2) SiBcrop: compile sib3 model only, recompiling dependencies " \
	"as necessary."
	@echo "3) sibmerge: compile sibmerge program"
	@echo "4) <object name>.o: compile specific object. matches the file" \
	"name."
	@echo "5) remake: delete all objects, and re-compile the sib3 model " \
	"and sibmerge."
	@echo "6) clean: remove all compiler-generated files and executable."
	@echo "7) help: print this message."
	@echo ""
