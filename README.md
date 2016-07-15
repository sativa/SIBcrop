# SIBcrop

SiB is is the Simple Biosphere model begun by Piers Sellers. It models soil and plant level biochemical processes affecting the atmospheric carbon cycle. It is now actively maintained by Scott Denning's Biocycle Research Group at CSU. There are a number of branches of development, but we are only considering SiBcrop at this time.

Bibliography

Sellers, P.J., D.A. Randall, G.J. Collatz, J.A. Berry, C.B. Field, D.A. Dazlich, C. Zhang, G.D. Collelo, and L. L. Bounoua, 1996a: A Revised Land Surface Parameterization (SiB2) for Atmospheric GCMs. Part I: Model Formulation. J. Climate, 9, 676-705.

Sellers, P.J., S.O. Los, C.J. Tucker, C.O. Justice, D.A. Dazlich, G.J. Collatz, and D.A. Randall, 1996b: A Revised Land Surface Parameteriztion (SiB2) for Atmospheric GCMs. Part II: The Generation of Global Fields of Terrestrial Biophysical Parameters from Satellite Data. J. Climate, 9,706-737.

Baker, I.T. and Denning, A.S. 2008. SiB3 modeled global 1-degree hourly biosphere-atmosphere carbon flux, 1998-2006. Data set. Oak Ridge National Laboratory Distributed Active Archive Center, Oak Ridge, Tennessee, U.S.A., http://daac.ornl.gov/.

SiBcrop

This is SiB with crop phenology added by Erandi Lokupitiya.

Bibliography

Lokupitiya, E., Denning, A.S., Paustian, K., Baker, I.T., Schaefer, K. and co-authors. 2009. Incorporation of crop phenology in Simple Biosphere Model (SiBcrop) to improve land-atmosphere carbon exchanges from croplands. Biogeosciences 6, 969-986.

Supported Platforms

Compilers:
GCC Gfortran,
Portland Group (pgf90), and
Intel Fortran (ifort)
Operating Systems:
Linux,
OS X
We work with SiB on Linux and OS X with Gfortran, and Portland Group's pgf90. Sometimes we test with Intel's ifort. Other Fortran 90 compilers and POSIX operating systems may work, but you're on your own. If some simple changes will make SiB run on your platform of choice, consider sending them to us for inclusion in our source code repository.

Prerequisites

You will need the libraries:
NetCDF version 3,
LAPACK
BLAS
All versions of SiB depend on NetCDF version 3, LAPACK and BLAS libraries. LAPACK and BLAS may come with your compiler, and on OS X the routines are in the vecLib framework. Otherwise, you will need to compile them yourself. We recommend building the ATLAS implementation of BLAS with netlib.org's LAPACK distribution.


Compile

We provide an example Makefile (Makefile.example) to compile SiB. It requires GNU make, which is the default on the supported platforms. It attempts to detect your operating system and configure itself accordingly, but you may still want to make some changes. Copy Makefile.example to Makefile and make your site-specific modifications to the latter. You will usually want to make these changes:
select the default compiler by uncommenting one of the COMPILER = lines,
select the optimization level by uncommenting one of the OPT = lines,
change NETCDF_ROOT to the place you installed NetCDF, and
change LAPACK_ROOT to your install of LAPACK and BLAS (not needed if you use -framework veclib).
You can also change these values without editing the Makefile by setting them as environment variables before issuing the make command. Example:

    $ export NETCDF_ROOT=/usr/local
    $ COMPILER=ifort OPT=opt make clean all
See the comments in the Makefile for details.

Now run make. If everything is configured correctly this will build the SiB3crop and sibmerge programs; if not, make sure everything is setup correctly in the Makefile.

Download Datasets

Get one of these datasets if you do not have your own.

Bondville—Bondville, Illinois, USA; corn and soy crops over two years.
Ponca—Ponca, Oklahoma, USA; winter wheat.
Run

Get a dataset (one we have provided or one you have created). There should be a namel_sibdrv file which tells SiB where all these input files and output directory are. You may need to edit these to change the paths. We recommend you keep your sib code in the SiB3crop/ directory that Subversion creates by default, and your data directories in the same directory SiB3crop/ is in, i.e.:

    bishop@chysis:~/src$ ls my_sib_project/
    SiB3crop/  data.bondville/  data.ponca.tar.gz  data.ponca/
Once you got your directory structure in order, run SiBD3crop from the data directory with the namel_sibdrv:
    bishop@chysis:~/src/my_sib_project/data.bondville$ ../SiB3crop/SiBD3crop
     INIT_GRID:
        reading sib inlist
        reading sib i/olst
        reading subgrid values
        reading sib pbplst
        reading sib_control_lst
        SiB time step (s) =                   600
        SiB out written (months) =                     1
        SiB restart written (months) =                     1
        nsib=            1
        drvr_type= single  
        reading parameter file: param/sib_param_TI.nc
     INIT_SIBDRV:
        diagnostics initialized
        initialize time variables
        reading time-invariant boundary conditions
           reading param/sib_param_TI.nc
        obtaining previous month time-variant boundary conditions
           reading file: param/sib_param_1998
        reading in initial conditions: param/sib_restart.nc
           opening ic fileparam/sib_restart.nc
           reading in slabs...
           load data into the structure
        setting soil properties 
        reading in initial time-step driver data
        reading in respFactor
        initializing crops
        initializing solar declination
     January                       1                 1999
     January                       2                 1999
     January                       3                 1999
     January                       4                 1999
     January                       5                 1999
     January                       6                 1999
     January                       7                 1999
     .
     .
     .
    
As shown above, you should see initialization information, then as the model runs it prints out each day it completes. When the run is done, the output will be in a series of NetCDF files. Look in the output/ directory in our datasets. Read the data with your favorite NetCDF tools.
