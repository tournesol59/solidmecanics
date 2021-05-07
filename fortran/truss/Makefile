####################################################################
#   Makefile, erzeugt die Library "libTrussPkkum.a"                 #
####################################################################

#There are many ways to write makefiles, but a good, simple model to start from for exploring pattern rules would be a combination of your first and second examples:

#all: launch.exe

#launch.exe: source.o 2ndfile.o
#    g++ -o $@ $^

#%.o: %.cpp
#    g++ -c $< -std=gnu++11

####################################################################
# moegliche Dateiendungen werden bekanntgegeben
#
.SUFFIXES: .d .mod .o .a .f .f90 


####################################################################
# Macros fuer den Bibliotheksaufbau 
#
LIB = libTrussPkkum.a
AR = ar
ARFLAGS = -crU

#Compiler Wahl fuer Linux:
#FCOMP="/opt/intel/compiler70/ia32/bin/ifc"
#FCOMP="/opt/intel_fc_80/bin/ifort"
FCOMP=gfortran

#Compiler Optionen fuer Linux:
#FC_FLAGS="-O2 -autodouble"
FC_FLAGS=-O2 -g -cpp
LDFLAGS=

####################################################################
# Maschinenabhaengige Einstellungen
# (Compiler und Optionen)
#
#SGI_FLAGS = \
        F90="f90" \
        F90_CFLAGS="-C -O2 -I Modules" \
	MOD_FILES="if test -e *.mod ; then mv -f *.mod Modules; fi " \
##	ALLOCATEFIELDS="Modules/allocateFields_mod.mod" \
##	TYPES="Modules/TYPES.mod" \

ifeq ($(FCOMP),"notgfortran")
LINUX_FLAGS = \
        F90=$(FCOMP) \
        F90_CFLAGS=$(FC_FLAGS) \
	MOD_FILES=""\
	TYPESBALKEN="TYPESBALKEN.mod" \
	INPUTTRUSS="INPUTTRUSS_mod.mod" \
        CREATETRUSS="CREATETRUSS_mod.mod" \
	PRINTABSTRACTSTRUCT="PRINTABSTRACTSTRUCT_mod.mod" \
	RIGIDTRUSSCONTACTS="RIGIDTRUSSCONTACTS_mod.mod" \
        CREATEBEAMTRIANGLES="CREATEBEAMTRIANGLES_mod.mod" \
        FORCESLIMITS="FORCESLIMITS_mod" \
        SETTIMOBEAMELEMENT="SETTIMOBEAMELEMENT"
else
LINUX_FLAGS = \
        F90=$(FCOMP) \
        F90_CFLAGS=$(FC_FLAGS) \
	MOD_FILES="" \
	TYPESBALKEN="typesbalken.mod" \
	INPUTTRUSS="inputtruss_mod.mod" \
        CREATETRUSS="createtruss_mod.mod" \
	PRINTABSTRACTSTRUCT="printabstractstruct_mod.mod" \
	RIGIDTRUSSCONTACTS="rigidtrusscontacts_mod.mod" \
        CREATEBEAMTRIANGLES="createbeamtriangles_mod.mod" \
        FORCESLIMITS="forceslimits_mod.mod" \
        SETTIMOBEAMELEMENT="settimobeamelement"
endif
####################################################################
# Targets fuer alle Maschinen, die zugelassen sind
#
linux:
#	touch modules.tar
#	tar xf modules.tar
	make $(LINUX_FLAGS) praktikum
	tar cf modules.tar *.mod *.o
	rm -f *.mod *.o

sgi:
	if ! test -d Modules ; then mkdir Modules ; fi
	make $(SGI_FLAGS) praktikum

######################################################
# Targets, die letztendlich das Programm Herstellen
#
library: $(LIB)(typesbalken.o) \
         $(LIB)(inputtruss.o) \
         $(LIB)(createtruss.o) \
	 $(LIB)(printabstractstruct.o) \
         $(LIB)(rigidtrusscontacts.o) \
         $(LIB)(createbeamtriangles.o) \
         $(LIB)(forceslimits.o) \
         $(LIB)(settimobeamelement.o)

OBJS: typesbalken.o \
      inputtruss.o \
      createtruss.o \
      printabstractstruct.o \
      rigidtrusscontacts.o \
      createbeamtriangles.o \
      forceslimits.o \
      settimobeamelement.o \
      mainTruss.o

# General rule for building prog from prog.o; $^ (GNU extension) is
# used in order to list additional object files on which the
# executable depends

praktikum: $(OBJS) mainTruss.o
#	$(FCOMP) $(FC_FLAGS) -o $@ $^ $(LDFLAGS)
# works only with a proper installation of lapack (not my debian)
	$(FCOMP) $(FC_CFLAGS) mainTruss.o -o Praktikum -L. -lTrussPkkum -g
#	$(FCOMP) $(FC_CFLAGS) main.o -o Praktikum -L. -lTrussPkkum -llapack -g
	rm -f main.o 
	@echo "============================================="
	@echo "   Program Praktikum successfully compiled   "
	@echo "============================================="

####################################################################
# Regeln zum Uebersetzen 
#
# frederic: try make NameOfTheModule.o each Module after each (yes the .o is important!)
#%.o : %.f90
#	$(F90) -c $(F90_CFLAGS) $<mv

%.o: %.f90
	$(FCOMP) $(FC_FLAGS) -o $@ -c $< 

#%.a : %.o
#	/usr/bin/ar rcv $@ $<
#	$(MOD_FILES)
##	rm -f $<

#Regeln fuer Module
$(TYPESBALKEN): typesbalken.o
$(INPUTTRUSS): inputtruss.o
$(CREATETRUSS): createtruss.o
$(PRINTABSTRACTSTRUCT): printabstractstruct.o
${RIGIDTRUSSCONTACTS}: rigidtrusscontacts.o
${CREATEBEAMTRIANGLES}: createbeamtriangles.o
${FORCESLIMITS}: forceslimits.o
${SETTIMOBEAMELEMENT}: settimobeamelement.o


######################################################
# Targets zum Aufraeumen
#
clean:
	rm -f *.o *.mod *.d modules.tar
	@echo "make clean completed..."

cleanall:
	rm -f Praktikum *.o *.mod *.d modules.tar 
	rm -f *.out *.dat *.plt *.dx *~ work.pc work.pcl ifc*
	rm -f $(LIB)
	rm -fr Modules
	@echo "make cleanall completed..."


######################################################
# Abhaengigkeiten
#
$(LIB)(printabstractstruct.o): $(TYPESBALKEN)
$(LIB)(rigidtrussstruct.o):$(PRINTABSTRACTSTRUCT) $(TYPESBALKEN) $(INPUTTRUSS)
$(LIB)(createtruss.o): $(TYPESBALKEN)
$(LIB)(createbeamtriangles.o): $(TYPESBALKEN)
$(LIB)(forceslimits.o): $(TYPESBALKEN)
$(LIB)(settimobeamelement.o): $(TYPESBALKEN)
##main.o: $(RIGIDTRUSSCONTACTS) $(SETTIMOBEAMELEMENT)
main.o: $(TYPESBALKEN) $(PRINTABTRACTSTRUCT) $(RIGIDTRUSSCONTACTS) $(INPUTTRUSS) $(CREATETRUSS) \
                $(CREATEBEAMTRIANGLES) $(FORCESLIMITS) $(SETTIMOBEAMELEMENT)
	@echo "make cleanall completed..."

