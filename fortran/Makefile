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
ARFLAGS = -crU

#Compiler Wahl fuer Linux:
#FCOMP="/opt/intel/compiler70/ia32/bin/ifc"
#FCOMP="/opt/intel_fc_80/bin/ifort"
FCOMP="gfortran"

#Compiler Optionen fuer Linux:
#FC_FLAGS="-O2 -autodouble"
FC_CFLAGS="-O2"

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
        INPUT="INPUT_mod" \
	TYPESROUND="TYPESROUND.mod" \
	ALLOCATEFIELDSROUND="ALLOCATEFIELDSROUND_mod.mod" \
	CREATEMESHROUND="CREATEMESHROUND_mod.mod" \
        LAGRANGEGRAD2DROUND="LAGRANGEGRAD2DROUND_mod.mod" \
        HERMITE2DDKT="HERMITE2DDKT_mod.mod" \
        SPACEGEOMETRY="SPACEGEOMETRY_mod.mod" \
        VERIFYMESHROUND="VERIFYMESHROUND_mod.mod" \
        TESTCONNRIGIDMAT="TESTCONNRIGIDMAT_mod.mod" \
        LAGRANGEROUND2D="LAGRANGEROUND2D_mod.mod"
else
LINUX_FLAGS = \
        F90=$(FCOMP) \
        F90_CFLAGS=$(FC_FLAGS) \
	MOD_FILES="" \
        INPUT="input_mod" \
	TYPESROUND="typesround.mod" \
	ALLOCATEFIELDSROUND="allocatefieldsround_mod.mod" \
	CREATEMESHROUND="createmeshround_mod.mod" \
        LAGRANGEGRAD2DROUND="Lagrangegrad2dround_mod.mod" \
        HERMITE2DDKT="hermite2dDKT_mod.mod" \
        SPACEGEOMETRY="spacegeometry_mod.mod" \
        VERIFYMESHROUND="verifymeshround_mod.mod" \
        TESTCONNRIGIDMAT="testconnrigidmat_mod.mod" \
        LAGRANGEROUND2D="lagrangeround2d_mod.mod"
endif

####################################################################
# Regeln zum Uebersetzen 
#
%.o : %.f90
	$(FCOMP) $(FC_CFLAGS) -g -cpp -o $@ -c $<

%.a : %.o
	/usr/bin/ar rcv $@ $<
	$(MOD_FILES)
#	rm -f $<

#Regeln fuer Module
$(TYPESROUND): typesrounddef.o
$(INPUT): input_mod.o
$(ALLOCATEFIELDSROUND): allocatefieldsround.o
${CREATEMESHROUND}: createmeshround.o
${LAGRANGEGRAD2DROUND}: Lagrangegrad2dround.o
${HERMITE2DDKT}: hermite2dDKT.o
${SPACEGEOMETRY}: spacegeometry.o
${VERIFYMESHROUND}: verifymeshround.o 
${TESTCONNRIGIDMAT}: testconnrigidmat.o
${LAGRANGEROUND2D}: lagrangeround2d.o
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
         $(LIB)(input_mod.o) \
	 $(LIB)(allocatefieldsround.o) \
	 $(LIB)(createmeshround.o) \
	 $(LIB)(Lagrangegrad2dround.o) \
	 $(LIB)(hermite2dDKT.o) \
         $(LIB)(verifymeshround.o) \
         $(LIB)(spacegeometry.o) \
         $(LIB)(testconnrigidmat.o) \
         $(LIB)(lagrangeround2d.o)

main:    mainRound.f90
	$(FCOMP) $(FC_CFLAGS) -g -cpp -c $< -o main.o
praktikum:library main.o
# works only with a proper installation of lapack (not my debian)
	$(FCOMP) $(FC_CFLAGS) main.o -o Praktikum -L. -lPraktikum -g -lgfortran
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
$(LIB)(input_mod.o): $(INPUT)
$(LIB)(createmeshround.o): $(TYPESROUND)
$(LIB)(allocatefieldsround.o): $(TYPESROUND)
$(LIB)(Lagrangegrad2dround.o): $(TYPESROUND)
$(LIB)(hermite2dDKT.o): $(TYPESROUND)
$(LIB)(spacegeometry.o): $(TYPESROUND)
$(LIB)(verifymeshround.o): $(SPACEGEOMETRY) $(TYPESROUND)
$(LIB)(testconnrigidmat.o): $(TYPESROUND)
$(LIB)(lagrangeround2d.o): $(LAGRANGEROUND2D)

main.o: $(TYPESROUND) $(INPUT) $(ALLOCATEFIELDSROUND) $(CREATEROUNDMESH) \
        $(HERMITE2DDKT) $(LAGRANGEGRAD2DROUND) $(TESTCONNRIGIDMAT) \
        $(LAGRANGEROUND2D) $(VERIFYMESHROUND)

	@echo "make cleanall completed..."





