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
	ALLOCATEFIELDS="Modules/allocateFields_mod.mod" \
	TYPES="Modules/TYPES.mod" \

ifeq ($(FCOMP),"notgfortran")
LINUX_FLAGS = \
        F90=$(FCOMP) \
        F90_CFLAGS=$(FC_FLAGS) \
	MOD_FILES=""\
	TYPES="TYPES.mod" \
	INPUT_MOD="Input_mod.mod" \
	ALLOCATEFIELDS="ALLOCATEFIELDS_mod.mod" \
	CREATEMESH="CREATEMESH_mod.mod" \
	INTERPOL2D="INTERPOL2D_mod.mod" \
        LAGRANGE2D="LAGRANGE2D_mod.mod"
else
LINUX_FLAGS = \
        F90=$(FCOMP) \
        F90_CFLAGS=$(FC_FLAGS) \
	MOD_FILES="" \
	TYPES="types.mod" \
	INPUT_MOD="Input_mod.mod" \
	ALLOCATEFIELDS="allocatefields_mod.mod" \
	CREATEMESH="createmesh_mod.mod" \
	INTERPOL2D="interpol2d_mod.mod" \
        LAGRANGE2D="lagrange2d_mod.mod"
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
$(TYPES): TypesDef.o
${INPUT_MOD}: Input_mod.o
$(ALLOCATEFIELDS): allocateFields.o
${CREATEMESH}: createMesh.o
{INTERPOL2D}: Interpol2D.o
{LAGRANGE2D}: Lagrange2D.o

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
library: $(LIB)(TypesDef.o) \
	 $(LIB)(Input_mod.o) \
	 $(LIB)(allocateFields.o) \
	 $(LIB)(createMesh.o) \
	 $(LIB)(Interpol2D.o) \
	 $(LIB)(Lagrange2D.o)


praktikum:library main.o
	$(F90) $(F90_CFLAGS) main.o -o Praktikum -L. -lPraktikum -llapack
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
$(LIB)(Input_mod.o): $(TYPES)
$(LIB)(createMesh.o): $(TYPES)
$(LIB)(allocateFields.o): $(TYPES)
$(LIB)(Interpol2D.o): $(TYPES)
$(LIB)(createMesh.o): $(TYPES)
$(LIB)(Lagrange2D.o): $(TYPES) $(INTERPOL2D) 
main.o: $(TYPES) $(ALLOCATEFIELDS) $(CREATEMESH) $(INTERPOL2D) \
         $(LAGRANGE2D)






