# Paths
DSRC = src
DOBJ = build
DEXE = app
DMOD = mod
EXEN = fgrav
all: $(DEXE)/$(EXEN)
# Flags
LIBS = 
FLAGS = -O3 -I$(DOBJ) -I$(DMOD) -fcheck=all -fbacktrace -g -ffree-line-length-none -fimplicit-none
CC = gfortran $(FLAGS) -J$(DMOD) $(LIBS) -c
CCL = gfortran -o

# Objects
OBJECTS = $(DOBJ)/constants.o $(DOBJ)/body_module.o $(DOBJ)/sim.o $(DOBJ)/out.o
MAIN_OBJ = $(DOBJ)/main.o


$(DOBJ)/main.o: $(DSRC)/main.f90 $(DOBJ)/constants.o $(DOBJ)/sim.o $(DOBJ)/body_module.o $(DOBJ)/out.o
$(DOBJ)/sim.o: $(DSRC)/sim.f90 $(DOBJ)/constants.o $(DOBJ)/body_module.o
$(DOBJ)/body_module.o: $(DSRC)/body_module.f90 $(DOBJ)/constants.o
$(DOBJ)/constants.o: $(DSRC)/constants.f90
$(DOBJ)/out.o: $(DSRC)/out.f90 $(DOBJ)/constants.o $(DOBJ)/body_module.o
# Default target


$(DOBJ)/%.o: $(DSRC)/%.f90 | $(DOBJ) $(DMOD)
	$(CC) $< -o $@


# Ensure required directories exist
$(DOBJ) $(DEXE) $(DMOD):
	mkdir -p $@

# Targets
$(DEXE)/$(EXEN): $(MAIN_OBJ) $(OBJECTS) | $(DEXE)
	$(CCL) $@ $(MAIN_OBJ) $(OBJECTS) $(LIBS)

main: $(DEXE)/$(EXEN)

run: $(DEXE)/$(EXEN)
	$(DEXE)/$(EXEN) $(ARG1) $(ARG2)



clean:
	rm -rf $(DOBJ)/*.o $(DEXE)/* $(DMOD)/*.mod *.dat 

.PHONY: all clean run main
