############################################################################
# Makefile for building the binary file for: SparseFI, J-SparseFI or J-SparseFI-HM
# Author: Claas Grohnfeldt
############################################################################

### local settings
# path to eigen library (you might need to adapt the version number 'f562a193118d' in the path)
EIGEN_HEADER_PATH = ./lib/eigen/eigen-eigen-f562a193118d/
# path to GDAL library
GDAL_INCLUDE_PATH = ./lib/gdal/inc
GDAL_LIBRARY_PATH = ./lib/gdal/lib

### internal path variables
COREDIR = .
OBJDIR = $(COREDIR)/obj
SRCDIR = $(COREDIR)/src
SUBDIR = $(COREDIR)/archive
BINDIR = $(COREDIR)/bin

### files
SRC = $(SRCDIR)/dataIO.cpp \
	$(SRCDIR)/filter.cpp \
	$(SRCDIR)/JS.cpp \
	$(SRCDIR)/auxFcts.cpp \
	$(SRCDIR)/eval_alg.cpp \
	$(SRCDIR)/mpi_counter.cpp \
	$(SRCDIR)/userSettings.cpp \
	$(SRCDIR)/paths.cpp \
	$(SRCDIR)/JSparseFI_alg.cpp \
	$(SRCDIR)/nnls.cpp

SRCMAIN = $(SRCDIR)/JSparseFI.cpp
EXE = $(BINDIR)/JSparseFIHM_supermuc

SUB = $(SUBDIR)/JSparseFIHM

HDR = $(wildcard ${SRCDIR}/*.h)
OBJ = $(addprefix $(OBJDIR)/,$(notdir $(SRC:.cpp=.o)))

### includes
INCFLAGS = -I $(EIGEN_HEADER_PATH) -I $(GDAL_INCLUDE_PATH) -I./ 

### libraries
LIB_GDAL = -lgdal
LIB_DIR  = -L $(GDAL_LIBRARY_PATH)

### compiler
CXX      = mpic++
#        = mpiCC   <- on SuperMUC
#        = mpic++  <- alternative
CFLAGS   = -O3 -w 
LDFLAGS=$(LIB_DIR) $(LIB_GDAL)

### exec
RUN      = mpiexec
RUNFLAGS = -n 2

########################################
#                                      #
#          program arguments           #
#                                      #
########################################

# 1: job name (string roughly describing this job)
LOADL_JOB_NAME=JSparseFIHM_demo_665211108350_Pavia_University

# 2: job ID (6 digit integer - on the SuperMUC this number will be automatically set. On a local machine it is set manually to an arbitrary string)
LOADL_PID=xxxxxx

# 3: dataset ID (string) 
datasetID=665211108350
# examples (see src/paths.cpp):
#__________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________
# data set ID           | orig. sensor    | LR sensor       | HR sensor       | fusion type | filter kernel | scene            | size ID       | fDS         | SNR          | redundancy digit # directory/link name
# --------------------- |------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# HS-MS (J-SparseFI-HM) 
#   665211108350    <=  | 6 (ROSIS)       - 6 (ROSIS)       - 5 (Quickbird)   - 2 (HS-MS)   - 1 (gauss)     - 1 (Pavia Univ.)  - 1 (560x320)   - 08 (fDS=8)  - 35 (SNR=35db)- 0                # 665211108350_ROSIS_Pavia_Univeristy
#   11119211105350  <=  | 11 (GRSSDFC2013)- 11 (GRSSDFC2013)- 9 (Sentinel2)   - 2 (HS-MS)   - 1 (gauss)     - 1 (Housten Univ.)- 1 (320x540)   - 05 (fDS=5)  - 35 (SNR=inf) - 0                # 11119211105350_2013IEEEGRSSDFC_Sentinel2_Univ
#   774212106350    <=  | 7 (Headwall)    - 7 (Headwall)    - 4 (WorldView-2) - 2 (HS-MS)   - 1 (gauss)     - 2 (Chikusei n.u.)- 1 (540x420)   - 06 (fDS=6)  - 35 (SNR=35db)- 0                # 774212106350_Headwall_Chikusei_nonUrban
#   885211404350    <=  | 8 (HYDICE)      - 8 (HYDICE)      - 5 (Quickbird)   - 2 (HS-MS)   - 1 (gauss)     - 1 (Wash.DC Mall) - 4 (420x300)   - 04 (fDS=4)  - 35 (SNR=35db)- 0                # 885211404350_HYDICE_WashDC_Mall
#   335213304350    <=  | 3 (Aviris)      - 3 (Aviris)      - 5 (Quickbird)   - 2 (HS-MS)   - 1 (gauss)     - 3 (Moffett Field)- 3 (360x360)   - 04 (fDS=4)  - 35 (SNR=35db)- 0                # 335213304350_Aviris_Moffett_Field
#   3315211304350   <=  | 3 (Aviris)      - 3 (Aviris)     - 15 (WorldView-3) - 2 (HS-MS)   - 1 (gauss)     - 1 (Ind.Pines)    - 3 (360x360)   - 04 (fDS=4)  - 35 (SNR=35db)- 0                # 3315211304350_Aviris_IndianPines_WV3_VNIR_SWIR
#   3315212405350   <=  | 3 (Aviris)      - 3 (Aviris)     - 15 (WorldView-3) - 2 (HS-MS)   - 1 (gauss)     - 2 (Cuprite Sc03) - 4 (420x360)   - 05 (fDS=5)  - 35 (SNR=35db)- 0                # 3315212405350_Aviris_Cuprite_sc03_WV3_VNIR_SWIR
#   109211103990    <=  |   (HyMap)       - 10 (EnMAP)      - 9 (Sentinel2)   - 2 (HS-MS)   - 1 (gauss)     - 1 (Rodalquilar)  - 1 (261x867)   - 03 (fDS=3)  - 99 (SNR=inf) - 0                # 109211103990_EnMAP_Sentinel2
#                       
# MS-PAN (J-SparseFI)   
# 155111203350 HySpex 3600x1200 QB
# 2444101104000_Real_WV2_scene
# 2444102104000_Real_WV2_HongKong_from_Naoto
#   114211504350    <=  | 1 (HySpex)      - 1 (HySpex)      - 4 (WorldView-2) - 2 (HS-MS)   - 1 (gauss)     - 1 (MUC_Oly)      - 5 (540x540)   - 04 (fDS=4)  - 35 (SNR=35db)- 0                #
#   114211510350    <=  | 1 (HySpex)      - 1 (HySpex)      - 4 (WorldView-2) - 2 (HS-MS)   - 1 (gauss)     - 1 (MUC_Oly)      - 5 (540x540)   - 10 (fDS=10) - 35 (SNR=35db)- 0                #
#   114211515350    <=  | 1 (HySpex)      - 1 (HySpex)      - 4 (WorldView-2) - 2 (HS-MS)   - 1 (gauss)     - 1 (MUC_Oly)      - 5 (540x540)   - 15 (fDS=15) - 35 (SNR=35db)- 0                #
# --------------------- |------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# 4: fusion method (string)
alg=JSparseFIHM
# options: 
# - JSparseFIHM
# - GroupedJSparseFI
# - SparseFI

#########################################################
# local-non-local processing module
#########################################################

# 5: regularization parameter trading the weighting of the l_2,1 norm term in the joint sparsity optimization problem (float)
lambda=1e0

#*************************************
# dictionary generation
#*************************************

# 6: use simulated high resolution image X for dictionary learning (bool)
useSimulatedImXforDictLearn=1

# 7: ImX simulation mode (int) 
ImX_sim_mode=2
# option:
# 0: correlation based
# 1: unconstrained least-squares based 
# 2: non-negative least-squares based

# 8: Select coupled LR and HR Pan dictionaries according to:
dictselect=8
# option: 
# 0: Dictionary contains ONLY the current patch -> (NDP=1) & Alpha is calculated by least squares
# 1: Nearest Neighbor
# 2: PanLR norm
# 3: SRF approximate PanLR  norm
# 4: PanHR norm (POSITIVE PanHR correlation)   (depricated)
# 5: PanLR-PanHR joint ranking                 (depricated)
# 6: ABSOLUTE PanHR correlation
# 7: PanHR uncorrelation, including current patch as first atom
# 8: Random, including current patch as first atom
# 9: PanHR self uncorrelated basis approximation, including current patch as first atom

# 9: number of dictionary atoms/patches (int)
NDP=900

# 10: patchsize, measured in low resolution pixels (int)
psz=3

# 11: patch overlap, measured in low resolution pixels (int)
ovrlp=0
# options:
# any int between 0 and psz-1

#*************************************
# coefficient estimation
#*************************************

# 12~14: regularization parameters (double)
mu_X=1e-5 # 316
mu_Y=1e-5 # 0.316
mu_Z=1

#*************************************
# correlation-based Hyperspectral Grouping (CorHySpeG)
#*************************************

# 15: maximum size of spectral group above groups will be double-checked and perhaps split into subgoups (int) 
Nc_max=25

# 16: minimum cross-correlation within spectral goups (double)
CC_min=0.96

# 17: size of window around patch: Must be odd if patchsize is odd and even if patchsize is even, in order to have both centers matched; Used for correlation calculations (int)
winSize=6

#######################################################
# Global processing module                             
# (for full-image optimization)                        
#######################################################

# 18: regularization parameter trading the relative weighting of the high resolution input image I_X (double)
mu_X_prime=1e0

# 19: regularization parameter trading the relative weighting of the low resolution input image I_Y (double)
mu_Y_prime=1e0

# 20: maximum number of iterations in the GS step to solve the least squares problem on the final image level (int) 
maxiter_globalOpt=1500 

# 21: error tolerance (double) 
tol_r_globalOpt=1e-12 

# 22: do full image optimization on the full image space (not on a subspace); 1: in the full image optimization step, operate on a subspace (optimize for A, where  ImZ = E*A and E is known) (only effenctive if LQ_post_opt_im is set to 1 above)  (bool)
fullImOptOnSubspace=1

# 23: subspace transformation type
subspace_transform_type=SVD
# options:
#                      =PCA
#                      =SVD
#                      =VCA
#                      =none

# 24: subspace dimension
subspace_dim=10

# 25: calc. coeff. in full image opt. eq. via SNR calc. of ImX and ImY (bool)
SNR_normalization=1

######################################################### 
# initialization and interplay between local-non-local  
# and global processing module 
######################################################### 

# 26: use estimated SRFs instead of apriori given ones (bool)
use_estimated_SRFs=1

# 27: type of initial high resolution image ImZ_init; flag (int)
ImZ_init_type=2
# options:
# 0: # mu_Z=0 in 1st iter (no initial image)
# 1: # upsampled and bilinearly interpolated low resolution image ImY
# 2: # reconstruction result of another algorithm (e.g. CNMF or HySure)
# 3: # reference image (for tests only)

# 28: use/activate global processing module (bool) 
use_global_proc_module=1

# 29: skip local-non-local processing module, i.e. jump straight to global proceccing module. (bool)
use_ONLY_global_proc_module=0

# 30: number of coupled ImZ calculations, i.e. local-non-local / glocal processing iterations 
iterMain=1

########################################################
# Output settings                                      #
########################################################

# 31: evaluate reconstructed image(s) (bool)
eval=1
# options:
# 1: evaluation will be done immetiately after the reconstruction 
# 0: evaluation will have to be done in a separate post-processing step)

# 32: evaluate initial image (bool)
evaluate_ImZ_init=1

# 33: write fused image in file (1: create file and write resulting image in file; 0: to not write image in file (useful for analyses only)) (bool)
writeImageFile=1

# 34: write all intermediate image fusion resulta (after every iteration) (1: create file and write resulting image in file; 0: to not write image in file (useful for analyses only)) (bool)
writeImageFileAfterEveryIter=1

# 35: save output in double format (64bit) instead of uint16 (bool)
saveAsDouble=1


########################################################################
#                                                                      #
#                  DON'T MODIFY BEYOND THIS LINE                       #
#                                                                      #
########################################################################
.PHONY: clean all
 
all: clean $(EXE)  

$(OBJ): $(OBJDIR)/%.o: $(SRCDIR)/%.cpp
	@$(CXX) $(CFLAGS) $(INCFLAGS) -c $< -o $@

$(EXE): $(OBJ)
	$(CXX) $(SRCMAIN) $(CFLAGS) $(INCFLAGS) $(OBJ) $(LDFLAGS) -o $(EXE)

run_JSparseFI:
		$(RUN) $(RUNFLAGS) -x LD_LIBRARY_PATH=$(GDAL_LIBRARY_PATH) $(EXE) $(LOADL_JOB_NAME) $(LOADL_PID) \
		$(datasetID) $(alg) $(lambda) $(useSimulatedImXforDictLearn) $(ImX_sim_mode) \
		$(dictselect) $(NDP) $(psz) $(ovrlp) $(mu_X) $(mu_Y) $(mu_Z) $(Nc_max) $(CC_min) \
		$(winSize) $(mu_X_prime) $(mu_Y_prime) $(maxiter_globalOpt) $(tol_r_globalOpt) \
		$(fullImOptOnSubspace) $(subspace_transform_type) $(subspace_dim) $(SNR_normalization) \
		$(use_estimated_SRFs) $(ImZ_init_type) $(use_global_proc_module) $(use_ONLY_global_proc_module) \
		$(iterMain) $(eval) $(evaluate_ImZ_init) $(writeImageFile) $(writeImageFileAfterEveryIter) $(saveAsDouble) 
		
		# $(datasetID) $(lambda) $(NDP) $(twoStep) $(Nc) $(No) $(psz) $(ovrlp) $(alg) $(eval) \
		# $(tol_SRF) $(store_patches_tmp_on_drive) $(parWrNumProc) $(saveAlphas) $(pFirstAlpha) $(pLastAlpha) \
		# $(saveDicts) $(pFirstDict) $(pLastDict) $(chBundleFirst) $(chBundleLast) \
		# $(uLFirst) $(uLLast) $(vLFirst) $(vLLast) \
		# $(numProcPerPatch) $(workStealingTurns) $(dictselect) \
		# $(writeImageFile) $(delete_tmp_patch_folders) $(imageConstructionOnly) $(contUnfinishedRec) \
		# $(PathToIncompletePatchSetCSV) $(dir_tmp_patches_additional_num) $(matrixNorm) $(addMeanPixelwise) \
		# $(maxiter_out) $(tol) $(LQ_post_opt) $(lambdaX) $(lambdaY) $(maxiter_CGLS) $(tol_r_CGLS) $(fix_Alpha) $(fix_delta_m) \
		# $(LQ_post_opt_im) $(lambdaX_im) $(lambdaY_im) $(maxiter_CGLS_im) $(tol_r_CGLS_im) \
		# $(useNewMethodForCalculatingZ) $(useSimulatedImXforDictLearn) $(lambdaX_ABC) $(lambdaY_ABC) $(lambdaZ_ABC) $(lambdaZ_ABC_in_1st_iter) $(ImZ_init_type) \
		# $(doFullImOptWithoutPatRec) $(iterMain) $(Nc_max) $(CC_min) $(winSize) \
		# $(evaluate_ImZ_init) $(set_neg_to_0) $(use_estimated_SRFs) $(writeImageFileAfterEveryIter) $(fullImOptOnSubspace) \
		# $(platformID) $(subspace_transform_type) $(subspace_dim) $(ImX_sim_mode) $(SNR_normalization) $(balance_ImX_term_coef) $(saveAsDouble) $(use_LRnorm_for_dic_normalization) $(load_DictHR_and_DictLR)\
		# $(dir_tmp_patches_additional_1) $(dir_tmp_patches_additional_2) $(dir_tmp_patches_additional_3)

clean:
	rm -f $(OBJ) $(EXE) $(SUB).tar.gz *~ Depends

submit:
	tar czvf $(SUB)$(shell date +_%Y_%m_%d).tar.gz $(EXE) $(SRCMAIN) $(SRC) $(HDR) Makefile 

.SUFFIXES : .o .cpp

Depends depend:
	$(CXX) -E -MM $(INCFLAGS) $(SRCMAIN) $(SRC) > Depends
