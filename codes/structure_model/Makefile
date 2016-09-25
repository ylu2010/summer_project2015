#=======================================================
# makefile for galaxy structure model code
#OPT	+= -DMASS
#OPT	+= -DVEL
#=======================================================
CC	= gcc
CXX = g++
AR	= ar r

#
EXEC	= gstructure
PRESAMP_EXEC	= sample_gstructure
TESTMCMC_EXEC	= test_mcmc

GAL_SOURCE	= evolve_galaxy.c init.c model_halo_accretion_history.c model_starformation_molecule.c cosmo.c \
			model_hotgas_accretion.c model_coldgas_accretion.c model_reionization.c model_disc_size.c model_adiabatic_contraction.c \
			cooling.c misc.c setup_parameters.c \
			run_galaxy.c 
			
GAL_OBJS	= evolve_galaxy.o init.o model_halo_accretion_history.o model_starformation_molecule.o cosmo.o \
			model_hotgas_accretion.o model_coldgas_accretion.o model_reionization.o model_disc_size.o model_adiabatic_contraction.o \
			cooling.o misc.o setup_parameters.o \
			run_galaxy.o 
			
GAL_INCL	= cosmo.h variables.h  proto.h

IHS_SOURCE	= ihs.cpp
IHS_INCL	= ihs.hpp
IHS_OBJS	= ihs.o

PRESAMP_SRC	= main_run_pre_sample.c
PRESAMP_OBJS	= main_run_pre_sample.o

RUN_SOURCE	= main_run_single.c			
RUN_OBJS	= main_run_single.o 

MCMC_SRC	= interface_likelihood.c
MCMC_OBJS	= interface_likelihood.o
			
#
OPTIMIZE	= -O3 -Wall
LIBS	= -lm -lgsl -lgslcblas 
CFLAGS =   -g $(OPT) $(OPTIMIZE)

#
$(EXEC):	$(GAL_OBJS) $(GAL_INCL) $(IHS_OBJS) $(IHS_INCL) $(RUN_OBJS) $(PRESAMP_OBJS) $(MCMC_OBJS) Makefile
	$(CXX) $(RUN_OBJS) $(GAL_OBJS) $(IHS_OBJS) $(LIBS)  -o  $(EXEC)
	$(CXX) $(PRESAMP_OBJS) $(GAL_OBJS) $(IHS_OBJS) $(LIBS)  -o  $(PRESAMP_EXEC)
	$(CXX) $(MCMC_OBJS) $(GAL_OBJS) $(IHS_OBJS) $(LIBS)  -o  $(TESTMCMC_EXEC)

$(GAL_OBJS):	$(GAL_SOURCE) $(GAL_INCL)
	$(CC) -c $(OPTIMIZE) $(GAL_SOURCE)
	
$(IHS_OBJS):	$(IHS_SOURCE) $(IHS_INCL)	
	$(CXX) -c $(IHS_SOURCE)
	
$(RUN_OBJS): $(RUN_SOURCE)
	$(CC) -c $(RUN_SOURCE)

$(PRESAMP_OBJS): $(PRESAMP_SRC)
	$(CC) -c $(PRESAMP_SRC)
	
$(MCMC_OBJS): $(MCMC_SRC)
	$(CC) -c $(MCMC_SRC)

#========================================================
all:	$(GAL_OBJS) $(GAL_INCL) $(IHS_OBJS) $(IHS_INCL) $(RUN_OBJS) $(PRESAMP_OBJS) $(MCMC_OBJS) Makefile
	$(CXX) $(OPTIMIZE) $(RUN_OBJS) $(GAL_OBJS) $(IHS_OBJS) $(LIBS)  -o  $(EXEC)
	$(CXX) $(PRESAMP_OBJS) $(GAL_OBJS) $(IHS_OBJS) $(LIBS)  -o  $(PRESAMP_EXEC)
	$(CXX) $(MCMC_OBJS) $(GAL_OBJS) $(IHS_OBJS) $(LIBS)  -o  $(TESTMCMC_EXEC)
		
clean:
	rm -f $(GAL_OBJS) $(IHS_OBJS) $(RUN_OBJS) $(PRESAMP_OBJS) $(MCMC_OBJS)

tidy:
	rm -f $(GAL_OBJS) $(IHS_OBJS) $(RUN_OBJS) $(PRESAMP_OBJS) $(MCMC_OBJS) $(EXEC) $(PRESAMP_EXEC) $(TESTMCMC_EXEC)
	
#libihs.a: ihs.o
#	$(AR) $@ $^ 

#libgal.a: $(GAL_OBJS)
#	$(AR) $@ $^ 

