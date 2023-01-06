# Makefile for G.A.S. MPI version

N_PROC = 4

MODULESDIR   = ./src/modules/
UTILSDIR = ./src/utils/
MAINDIR  = ./src/main/
TESTSDIR = ./src/tests/
BINDIR   = ./bin/
LOGDIR   = ./logs/
VALIDDIR = ./validation/
PARAMFILE_DIR = ./param_files/
PARAMFILE = GAS_param.in

FORT    = gfortran
DEBUG   = -g -Wall -Warray-bounds -Wunderflow -fbounds-check -fbacktrace -frepack-arrays -finit-local-zero -ffree-line-length-none
OPTIM   = -O2 -frepack-arrays -finit-local-zero -ffree-line-length-none 
FFLAGS  = -cpp $(DEBUG)
LIB     = libGAS.a -llapack

MAIN_FILE = GAS.f08
TEST_FILE = tests.f08
MAIN_PROG = GAS
TEST_PROG = GAS_TESTS

# Compilation Options
##############################################################
##############################################################

OPTIONS = 

all: utils modules test_modules lib tests main clean

gas_tests: utils modules test_modules lib tests clean

gas: utils modules lib main clean

clean:
	rm -f *~
	rm -f $(UTILSDIR)*~
	rm -f $(LIBDIR)*~
	rm -f $(MAINDIR)*~
	rm -f $(LOGDIR)*.log
	rm -f $(VALIDDIR)*.log
	rm -f $(VALIDDIR)*.dat

utils:
	$(FORT) $(FFLAGS) $(OPTIONS) -c $(UTILSDIR)parameters.f08
	$(FORT) $(FFLAGS) $(OPTIONS) -c $(UTILSDIR)PrDi.f08
	$(FORT) $(FFLAGS) $(OPTIONS) -c $(UTILSDIR)log.f08
	$(FORT) $(FFLAGS) $(OPTIONS) -c $(UTILSDIR)config.f08
	$(FORT) $(FFLAGS) $(OPTIONS) -c $(UTILSDIR)model.f08

modules: 
	$(FORT) $(FFLAGS) $(OPTIONS) -c $(MODULESDIR)gas.f08
	$(FORT) $(FFLAGS) $(OPTIONS) -c $(MODULESDIR)status.f08
	$(FORT) $(FFLAGS) $(OPTIONS) -c $(MODULESDIR)scale.f08
	$(FORT) $(FFLAGS) $(OPTIONS) -c $(MODULESDIR)gsh.f08
	$(FORT) $(FFLAGS) $(OPTIONS) -c $(MODULESDIR)ssp.f08
	$(FORT) $(FFLAGS) $(OPTIONS) -c $(MODULESDIR)sp.f08

test_modules:
	$(FORT) $(FFLAGS) $(OPTIONS) -c $(TESTSDIR)gas_tests.f08
	$(FORT) $(FFLAGS) $(OPTIONS) -c $(TESTSDIR)scale_tests.f08
	$(FORT) $(FFLAGS) $(OPTIONS) -c $(TESTSDIR)gsh_tests.f08
	$(FORT) $(FFLAGS) $(OPTIONS) -c $(TESTSDIR)sp_tests.f08

lib:
	mv *.o ${BINDIR}
	ar -r ${BINDIR}libGAS.a ${BINDIR}*.o
	rm -f ${BINDIR}*.o

tests:
	${FORT} ${FFLAGS} ${OPTIONS} $(TESTSDIR)${TEST_FILE} -o ${BINDIR}$(TEST_PROG) ${BINDIR}$(LIB)
	
main:	
	${FORT} ${FFLAGS} ${OPTIONS} $(MAINDIR)${MAIN_FILE} -o ${BINDIR}$(MAIN_PROG) ${BINDIR}$(LIB)

run:
	rm -f *.log
#	qsub -q mpi ./run_onto_mib
	mpiexec -np $(N_PROC) ./$(MAIN_PROG) $(PARAMFILE_DIR)$(PARAMFILE)
