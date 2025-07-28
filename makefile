# Paths
DSRC = src
DOBJ = build
DEXE = app
DMOD = mod
EXEN = fgrav
EXEFMM = fmm
all: $(DEXE)/$(EXEN)

VPATH = $(DSRC):$(DTEST):$(DSRC)/fmm
# Flags
LIBS = 
FLAGS = -O3 -I$(DOBJ) -I$(DMOD) -fcheck=all -fbacktrace -g -ffree-line-length-none -fimplicit-none
CC = gfortran $(FLAGS) -J$(DMOD) $(LIBS) -c
CCL = gfortran -o

# Objects
OBJECTS = $(DOBJ)/constants.o $(DOBJ)/body_module.o $(DOBJ)/sim.o $(DOBJ)/out.o $(DOBJ)/fmm.o $(DOBJ)/fmm_math.o

MAIN_OBJ = $(DOBJ)/main.o

FMM_OBJ = $(DOBJ)/main_fmm.o

$(DOBJ)/main.o: $(DSRC)/main.f90 $(DOBJ)/constants.o $(DOBJ)/sim.o $(DOBJ)/body_module.o $(DOBJ)/out.o
$(DOBJ)/sim.o: $(DSRC)/sim.f90 $(DOBJ)/constants.o $(DOBJ)/body_module.o
$(DOBJ)/body_module.o: $(DSRC)/body_module.f90 $(DOBJ)/constants.o
$(DOBJ)/constants.o: $(DSRC)/constants.f90
$(DOBJ)/out.o: $(DSRC)/out.f90 $(DOBJ)/constants.o $(DOBJ)/body_module.o

$(DOBJ)/main_fmm.o: $(DSRC)/fmm/main_fmm.f90 $(DOBJ)/constants.o $(DOBJ)/fmm.o
$(DOBJ)/fmm.o: $(DSRC)/fmm/fmm.f90 $(DOBJ)/constants.o
$(DOBJ)/fmm_math.o: $(DSRC)/fmm/fmm_math.f90 $(DOBJ)/constants.o $(DOBJ)/fmm.o
# Default target


# Build rules
$(DOBJ)/%.o: %.f90 | $(DOBJ) $(DMOD)
	$(CC) $< -o $@


# Ensure required directories exist
$(DOBJ) $(DEXE) $(DMOD):
	mkdir -p $@

# Targets
$(DEXE)/$(EXEN): $(MAIN_OBJ) $(OBJECTS) | $(DEXE)
	$(CCL) $@ $(MAIN_OBJ) $(OBJECTS) $(LIBS)

# Targets
$(DEXE)/$(EXEFMM): $(FMM_OBJ) $(OBJECTS) | $(DEXE)
	$(CCL) $@ $(FMM_OBJ) $(OBJECTS) $(LIBS)

main: $(DEXE)/$(EXEN)

run: $(DEXE)/$(EXEN)
	$(DEXE)/$(EXEN)

fmm: $(DEXE)/$(EXEFMM)
	$(DEXE)/$(EXEFMM)


clean:
	rm -rf $(DOBJ)/*.o $(DEXE)/* $(DMOD)/*.mod *.dat 

.PHONY: all clean run main
