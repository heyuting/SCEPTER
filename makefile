# Start of the makefile
# Defining variables

FC            = gfortran
# FC            = ifort
FFLAGS = -O2 -Wall -cpp $(CPFLAGS)

CPFLAGS       = 
CPFLAGS       += -Dno_intr_findloc # need to use in cluster
# CPFLAGS       += -Dshow_PSDiter # showing iteration process during PSD calculation
# CPFLAGS       += -Dparallel_ON # testing parallelization
# CPFLAGS       += -Dnpar_in=1 # number of threads for parallelization 
# CPFLAGS       += -Dnpar_in=1 # number of threads for parallelization 
# CPFLAGS       += -Dparpsd_chk # checking parallelization results
# CPFLAGS       += -Dksld_chk # checking rate consts for sld species
CPFLAGS       += -DolddustPSD # using old PSD for dust (not user input but prescribed one)
# CPFLAGS       += -Ddisp_cnst=6.59754e-2 # forcing constant and common dispersion for aqueous species
# CPFLAGS       += -Dnrec_prof_in=120 # number of profile records
# CPFLAGS       += -Derrmtx_printout # 
CPFLAGS       += -Dmod_basalt_cmp # using basalt composition defined in <basalt_define.h>
# CPFLAGS       += -Ddef_flx_save_alltime # flux reported each integration (costs lots of bites)
# CPFLAGS       += -Dfull_flux_report # output all cumulative flux
# CPFLAGS       += -Ddisp_lim # limiting the display of results
# CPFLAGS       += -Ddiss_only # not allowing precipitation of minerals
# CPFLAGS       += -Dlim_minsld # limiting mineral lowest conc. 
# CPFLAGS       += -Dporoiter # do iteration for porosity  
# CPFLAGS       += -Dcalcw_full # fully coupled w calcuation  
# CPFLAGS       += -Ddispiter # showing PSD flux in each iteration   
# CPFLAGS       += -DdispPSDiter # showing PSD flux in each iteration   
# CPFLAGS       += -Dcalcporo_full # fully coupled poro calcuation  
# CPFLAGS       += -Diwtypein=0 # uplift type 0--cnst w, 1-- cnst poro*w, 2-- cnst (1-poro)*w, 3--- w-flexible, if not defined 0 is taken

ifeq ($(FC),gfortran)
  # CFLAGS        = -fcheck=all -g -O3  
  # CFLAGS        = -Wall -O3 -g -fcheck=all -ffpe-trap=invalid,zero,overflow -fbacktrace
  # CFLAGS        = -Wall -O3 -g -fcheck=all -fbacktrace
  CFLAGS        = -fimplicit-none  -Wall  -Wline-truncation  -Wcharacter-truncation  -Wsurprising  \
	  -Waliasing  -Wimplicit-interface  -Wunused-parameter  -fwhole-file  -fcheck=all  -std=gnu  -pedantic  -fbacktrace -O3

endif

ifeq ($(FC),ifort)
  CFLAGS        = -O3 -heap-arrays -g -traceback -check bounds -fp-stack-check -gen-interfaces -warn interfaces -check arg_temp_created 
endif 

# LDFLAGS       = -L/usr/local/lib
LDFLAGS       = 

LIBS          = -lopenblas

ifneq (,$(findstring -Dmod_basalt_cmp,$(CPFLAGS)))
  # Found -Dmod_basalt_cmp
  INC          = -I/home/yhs5/project/SCEPTER/data 
else
  # Not found
  INC          = 
endif

ifneq (,$(findstring -Dparallel_ON,$(CPFLAGS)))
  # Found -Dparallel_ON
  CFLAGS        += -fopenmp
else
  # Not found
endif

SRC           = \
                scepter_constants.f90 \
                scepter_variables.f90 \
                scepter_IO.f90 \
                scepter_findloc.f90 \
                scepter_input.f90 \
                scepter_thermodynamics.f90 \
                scepter_concentration.f90 \
                scepter_equilibrium.f90 \
                scepter_kinetics.f90 \
                scepter_psd.f90 \
                scepter_eq_charge.f90 \
                $(wildcard scepter*.f90)
PROGRAM       = scepter
OBJS          = $(SRC:.f90=.o)

all:            $(PROGRAM)

# Linking step
$(PROGRAM):     $(OBJS)
	$(FC) $(OBJS) -o $(PROGRAM) -cpp $(CPFLAGS) $(CFLAGS) $(LIBS) $(LDFLAGS) $(INC)

# Compilation rule
%.o: %.f90
	$(FC) $< -c -cpp $(CPFLAGS) $(CFLAGS) $(LIBS) $(LDFLAGS) $(INC)

clean:;         rm -f *.o  *.mod *~ $(PROGRAM)
blank:;         truncate -s 0 *.out
cleanall:;         rm -f *.o *.out *~ $(PROGRAM)

