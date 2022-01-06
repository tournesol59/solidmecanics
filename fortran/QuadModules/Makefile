####################################################################
#   Makefile, erzeugt die Library "libPraktikum.a"                 #
####################################################################


####################################################################
# moegliche Dateiendungen werden bekanntgegeben
#
.SUFFIXES: .d .mod .o .a .f .f90 


####################################################################
# Macros fuer den Bibliotheksaufbau 
#
LIB = libPraktikum.a
AR = ar
ARFLAGS = -cr

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
	TYPES="Modules/TYPES.mod" \
        PRINTABSTRACTSTRUCT="Modules/printabstractstruct_mod" \
	INPUT_MOD="Modules/Input_mod.mod" \
	ALLOCATEFIELDS="Modules/allocateFields_mod.mod" \
	CREATEMESH="Modules/createmesh_mod.mod" \
	INTERPOL2D="Modules/Interpol2d_mod.mod" \
        PATCHTEST2D5EL="Modules/Patchtest2d5el_mod.mod

ifeq ($(FCOMP),"notgfortran")
LINUX_FLAGS = \
        F90=$(FCOMP) \
        F90_CFLAGS=$(FC_FLAGS) \
	MOD_FILES=""\
	TYPES="TYPES.mod" \
        PRINTABSTRACTSTRUCT="PRINTABSTRACTSTRUCT_mod" \
	INPUT_MOD="INPUT_mod.mod" \
	ALLOCATEFIELDS="ALLOCATEFIELDS_mod.mod" \
	CREATEMESH="CREATEMESH_mod.mod" \
	INTERPOL2D="INTERPOL2D_mod.mod" \
        PATCHTEST2D5EL="PATCHTEST2D5EL_mod.mod
else
LINUX_FLAGS = \
        F90=$(FCOMP) \
        F90_CFLAGS=$(FC_FLAGS) \
	MOD_FILES="" \
	TYPES="types.mod" \
        PRINTABSTRACTSTRUCT="printabstractstruct_mod" \
	INPUT_MOD="Modules/Input_mod.mod" \
	ALLOCATEFIELDS="allocateFields_mod.mod" \
	CREATEMESH="createmesh_mod.mod" \
	INTERPOL2D="Interpol2d_mod.mod" \
        PATCHTEST2D5EL="Patchtest2d5el_mod.mod
endif

####################################################################
# Regeln zum Uebersetzen 
#
%.o : %.f90
	$(F90) -c $(F90_CFLAGS) $<

%.a : %.o
	/usr/bin/ar rcv $@ $<
	$(MOD_FILES)
#	rm -f $<

#Regeln fuer Module
$(TYPES): typesdef.o
$(PRINTABSTRACTSTRUCT): printabstractstruct.o
${INPUT_MOD}: input_mod.o
$(ALLOCATEFIELDS): allocatefields.o
${CREATEMESH}: createmesh.o
{INTERPOL2D}: interpol2d.o
{PATCHTEST2D5EL}: Patchtest2d5el.o

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
library: $(LIB)(typesdef.o) \
	 $(LIB)(input_mod.o) \
	 $(LIB)(allocatefields.o) \
	 $(LIB)(createmesh.o) \
	 $(LIB)(interpol2d.o) \
         $(LIB)(patchtest2d5el.o) \

praktikum:library main.o
# works only with a proper installation of lapack (not my debian)
	$(FCOMP) $(FC_CFLAGS) main.o -o Praktikum -L. -lPraktikum -llapack
#	$(FCOMP) $(FC_CFLAGS) main.o -o Praktikum -L. -lPraktikum
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
$(LIB)(printabstractstruct.o): $(TYPES)
$(LIB)(Input_mod.o): $(TYPES)
$(LIB)(createMesh.o): $(TYPES)
$(LIB)(allocateFields.o): $(TYPES)
$(LIB)(interpol2D.o): $(TYPES)
$(LIB)(createMesh.o): $(TYPES) 
$(LIB)(patchtest2d5el.o): $(TYPES)
main.o: $(TYPES) $(ALLOCATEFIELDS) $(CREATEMESH) $(INTERPOL2D) \
        $(PATCHTEST2D5EL)

	@echo "make cleanall completed..."





