# Makefile for G.A.S. MPI version

N_PROC = 4

LIBDIR         = ./lib/
MAINDIR        = ./main/
PARAM_FILE_DIR = ../param_files/

F90C    = mpif90
DEBUG   = -g -Wall -Warray-bounds -Wunderflow -fbounds-check -fbacktrace -frepack-arrays -finit-local-zero -ffree-line-length-none
OPTIM   = -O2 -frepack-arrays -finit-local-zero -ffree-line-length-none 
FFLAGS  = -cpp $(OPTIM)
LIB     = libGalICS.a -lcfitsio

MAIN    = GAS.f90
PROG    = GAS_MPI

#Compilation Options
##############################################################
# -DRECOVERY_MODE              : Generate FITS main output files from tmp files
# -DPRINTALL                   : Print run information on screen
# -DNOT_REMOVE_TMP_FILES       : Keep tmp files in the ouput directory
# -DNO_OUTPUTS                 : Do not generate main FITS output files
#                                WARNING, this option has to be use WITH: NOT_REMOVE_TMP_FILES but is imcompatible with RECOVERY_MODE
# -DPRINT_WARNING              : Print extension of warning messages
# -DLUMINOUS_PROCESSES         : Use MPI physical/luminous combination to compute stellar spectra and magnitudes
# -DONLY_LUM                   : CReate only luminous process (use with -DRECOVERY_MODE to write output files)
# -DGAL_SPECTRA                : Save young, old stellar population and galaxy spectra
# -DAGN_SPECTRUM               : Take into account AGN SED (Fritz+2006) in the complete galaxy SED
# -DFOLLOW                     : Run only on followed halos
# -DHORIZON                    : Use Horizon dm-simulation
# -DBOLSHOI                    : Use Bolshoi dm-simulation (may be in a future)
# -DUNLINKED_TREES             : Use unlinked tree, 
# 				 In these merger-trees, substructure are always defined by the link with their main halos is broken
#                                WARNING, this option is not compatible with the SUB_TRANSFER option !! 
# -TWOMYRS                     : Use 2Myr stellar timestep libraries
# -DCLEAN                      : apply cleaning algorithm on merger tree
# -DTREE_EVOLVE                : Compute baryonic evolution of haloes (build galaxies). 
#                                If not present only the dark matter properties will be saved
# -DPHOTOIONISATION            : Use (Gnedin+00 and Kravtsov+04) or Okamoto+08 model photoionization model 
#   -GNEDIN_2000               : If Gnedin and Kravtsov
#   -OKAMOTO_2008              : If Okamoto 
# -DREACCRETION                : Post accretion of photoionised/ejected gas
# -DCOLD_STREAMS               : Split baryonic accretion in cold and hot part
# -DSUB_QUENCHING              : No accretion onto substructure
# -DSUB_TRANSFER               : Apply gas transfer between sub and main structures (no available with -DUNLINKED_TREES)
# -DSUB_STRIPPING              : Take into account gas stripping on the diffuse gas contained in the disc (To be used with -DSUB_TRANSFER)
# -DNO_EVAP                    : force to keep all the hot gas into the atmosphere --> no evaporation
# -DPOLYTROPIC                 : use polytropic description of the hot atmosphere
# -DSN_FEEDBACK_PROP_TO_SFR    : use a SN ejecta rate proporional to the star formation rate
# -DTHERMAL_INSTABILITIES      : Take into account thermal instabilities in the cooling process (reduce the effective cooling rate)
##############################################################

# low level option
#OPTIONS = -DHORIZON -DUNLINKED_TREES -DCLEAN
# global model
OPTIONS = -DHORIZON -DUNLINKED_TREES -DCLEAN -DTREE_EVOLVE -DPHOTOIONISATION -DOKAMOTO_2008 -DCOLD_STREAMS -DSUB_QUENCHING -DREACCRETION -DTHERMAL_INSTABILITIES
# 
#OPTIONS = -DHORIZON -DUNLINKED_TREES -DCLEAN 

all: module prog clean

clean:
	rm -f *.o
	rm -f *.mod
	rm -f *.out
	rm -f *.a
	rm -f *~
	rm -f $(LIBDIR)*~
	rm -f $(MAINDIR)*~

module: 
	$(F90C) $(FFLAGS) $(OPTIONS) -c $(LIBDIR)global_variables.f90
	$(F90C) $(FFLAGS) $(OPTIONS) -c $(LIBDIR)PrDi.f90
	$(F90C) $(FFLAGS) $(OPTIONS) -c $(LIBDIR)IO.f90
	$(F90C) $(FFLAGS) $(OPTIONS) -c $(LIBDIR)utilities.f90
	$(F90C) $(FFLAGS) $(OPTIONS) -c $(LIBDIR)gas.f90
	$(F90C) $(FFLAGS) $(OPTIONS) -c $(LIBDIR)cooling.f90
	$(F90C) $(FFLAGS) $(OPTIONS) -c $(LIBDIR)stellar_population_library.f90
	$(F90C) $(FFLAGS) $(OPTIONS) -c $(LIBDIR)filters.f90
	$(F90C) $(FFLAGS) $(OPTIONS) -c $(LIBDIR)dust.f90
	$(F90C) $(FFLAGS) $(OPTIONS) -c $(LIBDIR)dm.f90
	$(F90C) $(FFLAGS) $(OPTIONS) -c $(LIBDIR)stars.f90
	$(F90C) $(FFLAGS) $(OPTIONS) -c $(LIBDIR)agn.f90
	$(F90C) $(FFLAGS) $(OPTIONS) -c $(LIBDIR)bulge.f90	
	$(F90C) $(FFLAGS) $(OPTIONS) -c $(LIBDIR)disc.f90
	$(F90C) $(FFLAGS) $(OPTIONS) -c $(LIBDIR)galaxy.f90
	$(F90C) $(FFLAGS) $(OPTIONS) -c $(LIBDIR)baryon_halo.f90
	$(F90C) $(FFLAGS) $(OPTIONS) -c $(LIBDIR)halo.f90
	$(F90C) $(FFLAGS) $(OPTIONS) -c $(LIBDIR)tree.f90
	$(F90C) $(FFLAGS) $(OPTIONS) -c $(LIBDIR)IO_next.f90

	ar -r libGalICS.a *.o 

prog:	
	${F90C} ${FFLAGS} ${OPTIONS} $(MAINDIR)${MAIN} -o $(PROG) $(LIB)
run:
	rm -f *.log
#	qsub -q mpi ./run_onto_mib
	mpiexec -np $(N_PROC) ./$(PROG) $(PARAM_FILE_DIR)GAS_param.in
