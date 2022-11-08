# Name:     makefile
# Function: Provides the make utility with instructions for
#           compiling the aorsa3d code into an executable module.
#
# 11/2021 - Revised for compilation in MN4.


# Names and directories
EXEC = xaorsa3d.mn
SRC_DIR = ./src
INCLUDE_DIR = ./src
OBJ_DIR = ./obj

OBJFILES = \
 $(OBJ_DIR)/kron3_mod.o \
 $(OBJ_DIR)/aorsa3din_mod.o \
 $(OBJ_DIR)/swim_global_data_mod.o \
 $(OBJ_DIR)/profile_mod.o \
 $(OBJ_DIR)/vmec_setup.o \
 $(OBJ_DIR)/qlsum.o \
 $(OBJ_DIR)/ql_myra.o \
 $(OBJ_DIR)/mets2aorsa.o \
 $(OBJ_DIR)/cauchy_ppart.o \
 $(OBJ_DIR)/vlog.o \
 $(OBJ_DIR)/aorsa3dMain.o \
 $(OBJ_DIR)/aorsaSubs.o \
 $(OBJ_DIR)/sigma.o \
 $(OBJ_DIR)/zfunction.o \
 $(OBJ_DIR)/ztable.o \
 $(OBJ_DIR)/current.o \
 $(OBJ_DIR)/fourier.o \
 $(OBJ_DIR)/assert.o \
 $(OBJ_DIR)/setupblacs.o \
 $(OBJ_DIR)/bessel.o \
 $(OBJ_DIR)/check.o

INC_FILES = $(SRC_DIR)/kron3_mod.f

# Compilers and compiler flags
COMMON_FLAGS = -O2 -save -r8 -i8 -I $(INCLUDE_DIR)

FLAGS77 =
FLAGS90 =
FLAGSEX = -mkl=cluster

COMPILE77 = mpiifort -c $(COMMON_FLAGS) $(FLAGS77)
COMPILE90 = mpiifort -c $(COMMON_FLAGS) $(FLAGS90)
COMPILEEX = mpiifort    $(COMMON_FLAGS) $(FLAGSEX)

# External libraries
PGPLOT_LIB = ./pgplot
INCS = -I$(NETCDF_INC)
LIBS = -L$(NETCDF_LIB) -lnetcdf -lnetcdff -L$(PGPLOT_LIB) -lpgplot


###########################################################################
# Compile the program
$(EXEC):   $(OBJFILES)
	$(COMPILEEX) -o ./$(EXEC)  $(OBJFILES) $(LIBS) $(FLAGSEX)


# Compiling objects
$(OBJ_DIR)/mets2aorsa.o:     $(SRC_DIR)/mets2aorsa.f90 $(INC_FILES)
	$(COMPILE90) -o $(OBJ_DIR)/mets2aorsa.o $(SRC_DIR)/mets2aorsa.f90

$(OBJ_DIR)/qlsum.o:  $(SRC_DIR)/qlsum.f90 $(INC_FILES)
	$(COMPILE90) -o $(OBJ_DIR)/qlsum.o $(SRC_DIR)/qlsum.f90

$(OBJ_DIR)/cauchy_ppart.o:   $(SRC_DIR)/cauchy_ppart.f90 $(INC_FILES)
	$(COMPILE90) -o $(OBJ_DIR)/cauchy_ppart.o $(SRC_DIR)/cauchy_ppart.f90

$(OBJ_DIR)/vlog.o:           $(SRC_DIR)/vlog.f90 $(INC_FILES)
	$(COMPILE77) -o $(OBJ_DIR)/vlog.o $(SRC_DIR)/vlog.f90

$(OBJ_DIR)/kron3_mod.o:      $(SRC_DIR)/kron3_mod.f $(INC_FILES)
	$(COMPILE77) -o $(OBJ_DIR)/kron3_mod.o $(SRC_DIR)/kron3_mod.f

$(OBJ_DIR)/aorsa3din_mod.o:  $(SRC_DIR)/aorsa3din_mod.f $(INC_FILES)
	$(COMPILE77) -o $(OBJ_DIR)/aorsa3din_mod.o $(SRC_DIR)/aorsa3din_mod.f

$(OBJ_DIR)/swim_global_data_mod.o:  $(SRC_DIR)/swim_global_data_mod.f $(INC_FILES)
	$(COMPILE90) -o $(OBJ_DIR)/swim_global_data_mod.o $(SRC_DIR)/swim_global_data_mod.f

$(OBJ_DIR)/profile_mod.o:    $(SRC_DIR)/profile_mod.f $(INC_FILES)
	$(COMPILE90) -o $(OBJ_DIR)/profile_mod.o $(SRC_DIR)/profile_mod.f

$(OBJ_DIR)/vmec_setup.o:     $(SRC_DIR)/vmec_setup.f90 $(INC_FILES)
	$(COMPILE90)  -o $(OBJ_DIR)/vmec_setup.o $(SRC_DIR)/vmec_setup.f90 $(INCS)

$(OBJ_DIR)/aorsa3dMain.o:    $(SRC_DIR)/aorsa3dMain.f $(INC_FILES)
	$(COMPILE77) -o $(OBJ_DIR)/aorsa3dMain.o $(SRC_DIR)/aorsa3dMain.f

$(OBJ_DIR)/aorsaSubs.o:      $(SRC_DIR)/aorsaSubs.f $(INC_FILES)
	$(COMPILE77) -o $(OBJ_DIR)/aorsaSubs.o $(SRC_DIR)/aorsaSubs.f

$(OBJ_DIR)/sigma.o:          $(SRC_DIR)/sigma.f $(INC_FILES)
	$(COMPILE77) -o $(OBJ_DIR)/sigma.o $(SRC_DIR)/sigma.f

$(OBJ_DIR)/zfunction.o:      $(SRC_DIR)/zfunction.f $(INC_FILES)
	$(COMPILE77) -o $(OBJ_DIR)/zfunction.o $(SRC_DIR)/zfunction.f

$(OBJ_DIR)/ztable.o:         $(SRC_DIR)/ztable.f90 $(INC_FILES)
	$(COMPILE90) -o $(OBJ_DIR)/ztable.o $(SRC_DIR)/ztable.f90

$(OBJ_DIR)/bessel.o:         $(SRC_DIR)/bessel.f $(INC_FILES)
	$(COMPILE77) -o $(OBJ_DIR)/bessel.o $(SRC_DIR)/bessel.f

$(OBJ_DIR)/current.o:        $(SRC_DIR)/current.f $(INC_FILES)
	$(COMPILE77) -o $(OBJ_DIR)/current.o $(SRC_DIR)/current.f

$(OBJ_DIR)/ql_myra.o:        $(SRC_DIR)/ql_myra.f $(INC_FILES)
	$(COMPILE77) -o $(OBJ_DIR)/ql_myra.o $(SRC_DIR)/ql_myra.f

$(OBJ_DIR)/fourier.o:        $(SRC_DIR)/fourier.f $(INC_FILES)
	$(COMPILE77) -o $(OBJ_DIR)/fourier.o $(SRC_DIR)/fourier.f

$(OBJ_DIR)/assert.o:         $(SRC_DIR)/assert.f90 $(INC_FILES)
	$(COMPILE77) -o $(OBJ_DIR)/assert.o $(SRC_DIR)/assert.f90

$(OBJ_DIR)/descinit.o:       $(SRC_DIR)/descinit.f90 $(INC_FILES)
	$(COMPILE77) -o $(OBJ_DIR)/descinit.o $(SRC_DIR)/descinit.f90

$(OBJ_DIR)/infog2l.o:        $(SRC_DIR)/infog2l.f90 $(INC_FILES)
	$(COMPILE77) -o $(OBJ_DIR)/infog2l.o $(SRC_DIR)/infog2l.f90

$(OBJ_DIR)/pxerbla.o:        $(SRC_DIR)/pxerbla.f90 $(INC_FILES)
	$(COMPILE77) -o $(OBJ_DIR)/pxerbla.o $(SRC_DIR)/pxerbla.f90

$(OBJ_DIR)/setupblacs.o:     $(SRC_DIR)/setupblacs.f $(INC_FILES)
	$(COMPILE77) -o $(OBJ_DIR)/setupblacs.o $(SRC_DIR)/setupblacs.f

$(OBJ_DIR)/check.o:          $(SRC_DIR)/check.f $(INC_FILES)
	$(COMPILE77) -o $(OBJ_DIR)/check.o $(SRC_DIR)/check.f

$(OBJ_DIR)/cauchy_mod.o:     $(SRC_DIR)/cauchy_mod.f $(INC_FILES)
	$(COMPILE77) -o $(OBJ_DIR)/cauchy_mod.o $(SRC_DIR)/cauchy_mod.f


all: $(EXEC)

clean:
	rm *.mod $(EXEC)
	rm $(OBJ_DIR)/*.o
