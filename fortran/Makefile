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
	TYPESROUND="TYPESROUND.mod" \
	ALLOCATEFIELDSROUND="ALLOCATEFIELDSROUND_mod.mod" \
	CREATEMESHROUND="CREATEMESHROUND_mod.mod" \
        LAGRANGEGRAD2DROUND="LAGRANGEGRAD2DROUND_mod.mod"
else
LINUX_FLAGS = \
        F90=$(FCOMP) \
        F90_CFLAGS=$(FC_FLAGS) \
	MOD_FILES="" \
	TYPESROUND="typesround.mod" \
	ALLOCATEFIELDSROUND="allocatefieldsround_mod.mod" \
	CREATEMESHROUND="createmeshround_mod.mod" \
        LAGRANGEGRAD2DROUND="Lagrangegrad2dround_mod.mod"
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
$(TYPESROUND): typesrounddef.o
$(ALLOCATEFIELDSROUND): allocatefieldsround.o
${CREATEMESHROUND}: createmeshround.o
{LAGRANGEGRAD2DROUND}: Lagrangegrad2dround.o

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
library: $(LIB)(typesrounddef.o) \
	 $(LIB)(allocatefieldsround.o) \
	 $(LIB)(createmeshround.o) \
	 $(LIB)(Lagrangegrad2dround.o) 

praktikum:library main.o
# works only with a proper installation of lapack (not my debian)
	$(FCOMP) $(FC_CFLAGS) main.o -o Praktikum -L. -lPraktikum
#	$(FCOMP) $(FC_CFLAGS) main.o -o Praktikum -L. -lPraktikum -llapack
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
$(LIB)(createmeshround.o): $(TYPESROUND)
$(LIB)(allocatefieldsround.o): $(TYPESROUND)
$(LIB)(Lagrangegrad2dround.o): $(TYPESROUND)
main.o: $(TYPESROUND) $(ALLOCATEFIELDSROUND) $(CREATEROUNDMESH) \
        $(LAGRANGEGRAD2DROUND) 

	@echo "make cleanall completed..."





