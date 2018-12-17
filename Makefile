############################################################################
# Makefile for building the binary file for: SparseFI, J-SparseFI or J-SparseFI-HM
# Author: Claas Grohnfeldt
############################################################################

# Directories and files
DIR_BASE = .
DIR_OBJ = $(DIR_BASE)/obj
DIR_SRC = $(DIR_BASE)/src
DIR_ARC = $(DIR_BASE)/archive
DIR_BIN = $(DIR_BASE)/bin

SRC = $(DIR_SRC)/dataIO.cpp \
      $(DIR_SRC)/filter.cpp \
      $(DIR_SRC)/JS.cpp \
      $(DIR_SRC)/auxFcts.cpp \
      $(DIR_SRC)/eval_alg.cpp \
      $(DIR_SRC)/mpi_counter.cpp \
      $(DIR_SRC)/userSettings.cpp \
      $(DIR_SRC)/paths.cpp \
      $(DIR_SRC)/JSparseFIHM_alg.cpp \
      $(DIR_SRC)/nnls.cpp

SRC_MAIN = $(DIR_SRC)/main.cpp
ARC      = $(DIR_ARC)/JSparseFIHM
EXE      = $(DIR_BIN)/JSparseFIHM

HDR = $(wildcard ${DIR_SRC}/*.h)
OBJ = $(addprefix $(DIR_OBJ)/,$(notdir $(SRC:.cpp=.o)))

### Libraries
PATH_LIBRARY_GDAL  = ./lib/gdal/lib
PATH_INCLUDE_EIGEN = ./lib/eigen
PATH_INCLUDE_GDAL  = ./lib/gdal/include
# inc
FLAGS_INC = -I $(PATH_INCLUDE_EIGEN) -I $(PATH_INCLUDE_GDAL) -I./ 
# lib
LIB_GDAL = -lgdal
LIB_DIR  = -L $(PATH_LIBRARY_GDAL)
FLAGS_LD=$(LIB_DIR) $(LIB_GDAL)

### Compiler preferences
CXX = mpic++
#   = mpiCC   <- on SuperMUC
#   = mpic++  <- more common on PCs
FLAGS_COMP = -O3 -w

########################################################################
#                                                                      #
#                  DON'T MODIFY BEYOND THIS LINE                       #
#                                                                      #
########################################################################
.PHONY: clean all
 
all: clean $(EXE)  

$(OBJ): $(DIR_OBJ)/%.o: $(DIR_SRC)/%.cpp
	@$(CXX) $(FLAGS_COMP) $(FLAGS_INC) -c $< -o $@

$(EXE): $(OBJ)
	$(CXX) $(SRC_MAIN) $(FLAGS_COMP) $(FLAGS_INC) $(OBJ) $(FLAGS_LD) -o $(EXE)

clean:
	rm -f $(OBJ) $(EXE) *~ Depends

archive_snapshot:
	tar czvf $(ARC)$(shell date +_%Y_%m_%d).tar.gz $(EXE) $(SRC_MAIN) $(SRC) $(HDR) Makefile 

.SUFFIXES : .o .cpp

Depends depend:
	$(CXX) -E -MM $(FLAGS_INC) $(SRC_MAIN) $(SRC) > Depends