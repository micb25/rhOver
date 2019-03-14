
# additional libraries
LIBS=-llapack

# symbols
SYMBOLS=-DDEBUG

# target name
TARGET = rhover

# global stuff
SRC  = global_c.f90 util_rot.f90 util_ylm.f90 util_str.f90 cgto.f90 int_overlap.f90
# libboys
SRC += libboys/libboys_data.f90 libboys/libboys.f90 
# grid
SRC += grid_lebedev.f90 grid.f90 radialwf.f90 sievers.f90 coulomb.f90 util_cf.f90 rhover.f90

OBJ = $(SRC:.f90=.o)

default:
	$(MAKE) -j1 global_c.o cgto.o grid_lebedev.o grid.o coulomb.o
	$(MAKE) $(OBJ)
	$(MAKE) -j1 $(TARGET)
	
%.o %.mod: %.f90
	$(FC) $(FLAGS) $(SYMBOLS) $(LIBS) -c -o $@ $<
	
$(TARGET): $(OBJ)
	$(FC) $(FLAGS) $(SYMBOLS) $(LIBS) -o $@ $^

clean:
	rm -f $(TARGET) *.mod *.o

