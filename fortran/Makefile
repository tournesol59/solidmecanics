####################################################################
#   Makefile, erzeugt das Program "solidmain(.exe)"                #
####################################################################

####################################################################
# moegliche Dateiendungen werden bekanntgegeben
#
.SUFFIXES: .d .mod .o .a .f .f90 

COMP=gfortran
LIB = libPraktikum.a
AR = ar
ARFLAGS = -cr

LINUX_FLAGS =\
        F90=$(COMP) \
        F90_CFLAGS=$(FLAGS) \
        FLAGS= -fdefault-real-8 -O0 -g -fbounds-check -Wall \
  -ffpe-trap=invalid,zero,overflow,underflow \
  -fbacktrace -ftrapv -fimplicit-none \
	MOD_FILES=""\
	TYPES="types.mod" \
	ALLOCATEFIELDS="allocatefields_mod.mod" \
	INPUT_MOD="input_mod.mod" \
	TRIDIAG="tridiag_mod.mod" \
        CREATEMESH="createmesh_mod.mod" \
        INTERPOL2D="interpol2d_mod.mod" \
        CALCGRADIENTSURF="calcgradientsurf_mod.mod" 

#FINCLUDE=../include/

LAPACKELIB=/usr/lib/lapack/lib/liblapacke.a

INCLAPACKE=/usr/lib/lapack/LAPACKE/include
TMGLIB=/usr/lib/lapack/lib/libtmglib.a
LAPACKLIB=/usr/lib/lapack/lib/liblapack.a
BLASLIB=/usr/lib/lapack/lib/librefblas.a

LIBS=-lm

####################################################################
# Regeln zum Uebersetzen 
#
%.o : %.f90
	$(F90) -c $(FLAGS) $<

%.a : %.o
	/usr/bin/ar rcv $@ $<
	$(MOD_FILES)
#	rm -f $<

#Regeln fuer Module
$(TYPES): TypesDef.o
$(INPUT_MOD): Input_mod.o
$(ALLOCATEFIELDS): allocateFields.o
$(CREATEMESH): createMesh.o
$(INTERPOL2D): interpol2d.o 
$(CALCGRADIENTSURF): calcgradientsurf.o
$(SOLIDMAIN): solidmain.o
#$(QUELLEN): Quellen.o
#$(AUSGABE): ausgabe.o
#$(TRIDIAG): tridiag.o

####################################################################
# Targets fuer alle Maschinen, die zugelassen sind
#
linux:
	touch modules.tar
	tar xf modules.tar
	make $(LINUX_FLAGS) solidpraktikum
	tar cf modules.tar *.mod *.o
	rm -f *.mod *.o

######################################################
# Targets, die letztendlich das Programm Herstellen
#
library: $(LIB)(TypesDef.o) \
	 $(LIB)(input_mod.o) \
	 $(LIB)(allocateFields.o) \
	 $(LIB)(createMesh.o) \
	 $(LIB)(interpol2d.o) \
	 $(LIB)(calcgradientsurf.o) \
	 $(LIB)(solidmain.o) 

solidmain: library solidmain.o
	$(F90) $(LINUX_FLAGS) $(TYPES) $(INPUT_MOD) $(ALLOCATEFIELDS) $(CREATEMESH) $(INTERPOL2D) $(CALCGRADIENTSURF) $(SOLIDMAIN) -o solidmain -L. -lsolidpraktikum -ltecio
	rm -f solidmain.o 
	@echo "============================================="
	@echo "   Program solidmain successfully compiled   "
	@echo "============================================="

clean:
	rm -rf *.o *.mod solidmain


######################################################
# Abhaengigkeiten
#
$(Input_mod.o): $(TYPES)
$(allocateFields.o): $(TYPES)
$(CREATEMESH): $(TYPES) $(ALLOCATEFIELDS)
$(INTERPOL2D): $(TYPES)
$(CALCGRADIENTSURF): $(TYPES) $(INTERPOL2D)
$(SOLIDMAIN):: $(TYPES) $(ALLOCATEFIELDS) $(INTERPOL2D)


