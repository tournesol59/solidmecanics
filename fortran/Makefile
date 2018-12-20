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
        ACCESSORIES="ACCESSORIES_mod.mod" \
	CREATEMESH="CREATEMESH_mod.mod" \
	INTERPOL2D="INTERPOL2D_mod.mod" \
        USEDTRIANGLE2D="USEDTRIANGLE2D_mod.mod" \
        FLEXIONGRAD2D="FLEXIONGRAD2D_mod.mod"
else
LINUX_FLAGS = \
        F90=$(FCOMP) \
        F90_CFLAGS=$(FC_FLAGS) \
	MOD_FILES="" \
	TYPES="types.mod" \
	INPUT_MOD="Input_mod.mod" \
	ALLOCATEFIELDS="allocatefields_mod.mod" \
        ACCESSORIES="accessories_mod.mod" \
	CREATEMESH="createmesh_mod.mod" \
	INTERPOL2D="interpol2d_mod.mod" \
        USEDTRIANGLE2D="usedtriangled2D_mod.mod" \
        FLEXIONGRAD2D="flexiongrad2d_mod.mod"
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
${INPUT_MOD}: input_mod.o
$(ALLOCATEFIELDS): allocatefields.o
$(ACCESSORIES): accessories.o
${CREATEMESH}: createmesh.o
{INTERPOL2D}: interpol2d.o
{USEDTRIANGLE2D}: usedtriangled2d.o
{FLEXIONGRAD2D}: flexiongrad2d.o 

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
         $(LIB)(accessories.o) \
	 $(LIB)(createmesh.o) \
	 $(LIB)(interpol2d.o) \
	 $(LIB)(usedtriangle2d.o) \
	 $(LIB)(flexiongrad2d.o) \

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
$(LIB)(Input_mod.o): $(TYPES)
$(LIB)(createMesh.o): $(TYPES)
$(LIB)(allocateFields.o): $(TYPES)
$(LIB)(accessories.o): $(TYPES)
$(LIB)(interpol2D.o): $(TYPES)
$(LIB)(createMesh.o): $(TYPES)
$(LIB)(usedtriangle2d.o): $(TYPES)
$(LIB)(flexiongrad2d.o): $(TYPES)
main.o: $(TYPES) $(ALLOCATEFIELDS) $(ACCESSORIES) $(CREATEMESH) $(INTERPOL2D) \
         $(USEDTRIANGLE2D) $(FLEXIONQUAD2D)

	@echo "make cleanall completed..."





