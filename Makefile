#OPT           =-02
OPT           =
 
CXX           =
ObjSuf        = o
SrcSuf        = cxx
ExeSuf        =
DllSuf        = so
OutPutOpt     = -o

ROOTCFLAGS   := $(shell root-config --cflags)
ROOTLIBS     := $(shell root-config --libs)
ROOTGLIBS    := $(shell root-config --glibs)

CXX           = g++
CXXFLAGS      = -Wall -fPIC
LD            = g++
LDFLAGS       =
SOFLAGS       = -shared

CXXFLAGS     += $(ROOTCFLAGS)

# Directories to look for fortran compiler
BINDIRS	      = . /bin /usr/bin /usr/local/bin

#Fortran
PRIMARYF77    = gfortran
ALTF77	      =  g77

F77           = $(PRIMARYF77)
F77LIBSO      = $(shell gfortran -print-file-name=libgfortran.so)
ALTF77LIBSO   = $(shell g77 -print-file-name=libg2c.so)

F77FLAGS      = -fPIC
F77OPT        = $(OPT)
CXXOUT        = -o # keep whitespace after "-o"

#Set fortran compiler 
EXIST = $(strip $(foreach dir,$(BINDIRS),$(wildcard $(dir)/$(PRIMARYF77) ) ) )
ifeq ($(EXIST),)
F77 = $(ALTF77)
F77LIBSO = $(ALTF77LIBSO)
EXIST = $(strip $(foreach dir,$(BINDIRS),$(wildcard $(dir)/$(ALTF77) ) ))
endif
#Alt 
ifeq ($(EXIST),)
$(error Fortran compiler not found. Either $(PRIMARYF77) or $(ALTF77) should be present)
endif

LIBS          = $(ROOTLIBS) $(SYSLIBS) 
GLIBS         = $(ROOTGLIBS) $(SYSLIBS)

#--------------HYDJET-----------------

HYDJETO        = Particle.o HadronDecayer.o InitialState.o InitialStateHydjet.o\
               UKUtility.o StrangeDensity.o StrangePotential.o EquationSolver.o\
	       HankelFunction.o GrandCanonical.o RunHadronSource.o RandArrayFunction.o\
	       DecayChannel.o ParticlePDG.o DatabasePDG.o\
	       PYQUEN/pythia-6.4.24.o PYQUEN/pyquen1_5.o PYQUEN/progs_fortran.o 
              
HYDJET	       = HYDJET

#--------------HYDJET_HISTO-----------------

HYDJETO_H      = Particle.o HadronDecayer.o InitialState.o InitialStateHydjet.o\
               UKUtility.o StrangeDensity.o StrangePotential.o EquationSolver.o\
	       HankelFunction.o GrandCanonical.o RunHadronSourceHISTO.o RandArrayFunction.o\
	       DecayChannel.o ParticlePDG.o DatabasePDG.o\
	       PYQUEN/pythia-6.4.24.o PYQUEN/pyquen1_5.o PYQUEN/progs_fortran.o 

HYDJET_H      = HYDJET_HISTO

#--------------------------------------------------------------------------------------------
all:	tree
tree:	$(HYDJET)

$(HYDJET):      $(HYDJETO) 
		$(LD)  -O2 $(LDFLAGS) $^ $(LIBS) $(F77LIBSO) $(OutPutOpt) $@ 
		@echo "$@ done"

histo:	$(HYDJET_H)

$(HYDJET_H):    $(HYDJETO_H)
		$(LD)  -O2 $(LDFLAGS) $^  $(LIBS) $(F77LIBSO) $(OutPutOpt) $@
		@echo "$@ done"		
		
clean:
		@rm -f $(HYDJETO) $(HYDJETO_H) $(HYDJET) $(HYDJET_H)

%.o : %.cxx
	$(CXX) -O2 $(CXXFLAGS) -c $<

%.o: %.f
	$(F77) $(F77OPT) $(F77FLAGS) $(CXXOUT)$@ -c $<
