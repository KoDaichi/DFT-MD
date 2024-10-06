include ../make.inc
OBJ = commons.o phase_commons.o tspace_defines.o \
    container_commons.o container_psicoef.o tspace_commons.o container_lattice.o \
    degeneracy.o transd_vector.o character_table.o \
    reduce_band.o compatibility.o classification.o

LIB = -I. $(OBJ) ../tspace/tsp98.o
BIN=..
all: install
install: $(OBJ) main.o main
	mv main $(BIN)/band_symm
