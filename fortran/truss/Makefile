####################################################################
#   Makefile, erzeugt die Library "libTrussPkkum.a"                 #
####################################################################


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
FCOMP="gfortran"

#Compiler Optionen fuer Linux:
#FC_FLAGS="-O2 -autodouble"
FC_FLAGS="-O2"

####################################################################
# Maschinenabhaengige Einstellungen
# (Compiler und Optionen)
#
SGI_FLAGS = \
        F90="f90" \
        F90_CFLAGS="-C -O2 -I Modules" \
	MOD_FILES="if test -e *.mod ; then mv -f *.mod Modules; fi " \
#	ALLOCATEFIELDS="Modules/allocateFields_mod.mod" \
#	TYPES="Modules/TYPES.mod" \

ifeq ($(FCOMP),"notgfortran")
LINUX_FLAGS = \
        F90=$(FCOMP) \
        F90_CFLAGS=$(FC_FLAGS) \
	MOD_FILES=""\
	TYPESBALKEN="TYPESBALKEN.mod" \
	INPUTTRUSS="INPUTTRUSS_mod.mod" \
        CREATETRUSS="CREATETRUSS_mod.mod" \
	PRINTABSTRACTSTRUCT="PRINTABSTRACTSTRUCT_mod.mod" \
	RIGIDTRUSSCONTACTS="RIGIDTRUSSCONTACTS_mod.mod" 
else
LINUX_FLAGS = \
        F90=$(FCOMP) \
        F90_CFLAGS=$(FC_FLAGS) \
	MOD_FILES="" \
	TYPESBALKEN="typesbalken.mod" \
	INPUTTRUSS="inputtruss_mod.mod" \
        CREATETRUSS="createtruss_mod.mod" \
	PRINTABSTRACTSTRUCT="printabstractstruct_mod.mod" \
	RIGIDTRUSSCONTACTS="rigidtrusscontacts_mod.mod"
endif

####################################################################
# Regeln zum Uebersetzen 
#
%.o : %.f90
	$(F90) -c $(F90_CFLAGS) $<mv

%.a : %.o
	/usr/bin/ar rcv $@ $<
	$(MOD_FILES)
#	rm -f $<

#Regeln fuer Module
$(TYPESBALKEN): typesbalken.o
$(INPUTTRUSS): inputtruss.o
$(CREATETRUSS): createtruss.o
$(PRINTABSTRACTSTRUCT): printabstractstruct.o
${RIGIDTRUSSCONTACTS}: rigidtrusscontacts.o

####################################################################
# Targets fuer alle Maschinen, die zugelassen sind
#
linux:
	touch modules.tar
	tar xf modules.tar
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
	 $(LIB)(rigidtrusscontacts.o)

praktikum:library main.o
# works only with a proper installation of lapack (not my debian)
	$(FCOMP) $(FC_CFLAGS) main.o -o Praktikum -L. -lTrussPkkum -g
#	$(FCOMP) $(FC_CFLAGS) main.o -o Praktikum -L. -lTrussPkkum -llapack -g
	rm -f main.o 
	@echo "============================================="
	@echo "   Program Praktikum successfully compiled   "
	@echo "============================================="



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
main.o: $(TYPESBALKEN) $(PRINTABTRACTSTRUCT) $(RIGIDTRUSSCONTACTS) $(INPUTTRUSS) $(CREATETRUSS)

	@echo "make cleanall completed..."





