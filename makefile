FC = ifort
FFLAGS = -g #-Wall -ffree-line-length-none -fcheck=all
FOPT = -O3
FLIBS = -llapack -lblas # -lccplib
EXEC = welding

SRC = main.f90 modelo_completo.f90
OBJ = main.o modelo_completo.o


all: $(OBJ)
	$(FC) $^ -o $(EXEC) $(FFLAGS) $(FLIBS)
$(OBJ): $(SRC)
	$(FC) $(FFLAGS) $(FOPT) $^ -c
