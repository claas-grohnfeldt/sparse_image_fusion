############################################################################
# Makefile for building the binary file for: SparseFI, J-SparseFI or J-SparseFI-HM
# Author: Claas Grohnfeldt
############################################################################

### local settings
include make.inc

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
EXE = $(BINDIR)/JSparseFIHM

SUB = $(SUBDIR)/JSparseFIHM

HDR = $(wildcard ${SRCDIR}/*.h)
OBJ = $(addprefix $(OBJDIR)/,$(notdir $(SRC:.cpp=.o)))

### includes
INCFLAGS = -I $(EIGEN_HEADER_PATH) -I $(GDAL_INCLUDE_PATH) -I./ 

### libraries
LIB_GDAL = -lgdal
LIB_DIR  = -L $(GDAL_LIBRARY_PATH)

### compiler
#CXX      = mpic++
CXX      = mpiCC
CFLAGS   = -O3 -w 
LDFLAGS=$(LIB_DIR) $(LIB_GDAL)

### exec
RUN      = mpiexec
RUNFLAGS = -n 2

########################################
###        program arguments         ###
########################################
#######################  1: job name (string roughly describing this job)
THIS_JOB_NAME=demo_665211108350_iter0_T_subspaceDim
#######################  2: job ID (6 digit integer - on the SuperMUC this number will be automatically set. On a local machine it is set manually to an arbitrary string)
THIS_JOB_ID=xxxxxx
#######################  3: dataset ID (int) 
datasetID=DEFINEDBELOW
LIST_datasetID = 665211108350
# examples:
# 
# HS-MS:
# 109211103990_EnMAP_Sentinel2
# 11119211105350_2013IEEEGRSSDFC_Sentinel2_Univ
# 3315211304350_Aviris_IndianPines_WV3_VNIR_SWIR
# 3315212405350_Aviris_Cuprite_sc03_WV3_VNIR_SWIR
# 335213304350_Aviris_Moffett_Field
# 665211108350_ROSIS_Pavia_Univeristy
# 774212106350_Headwall_Chikusei_nonUrban
# 885211404350_HYDICE_WashDC_Mall
#
# MS-PAN
# 155111203350 HySpex 3600x1200 QB 
# 2444101104000_Real_WV2_scene
# 2444102104000_Real_WV2_HongKong_from_Naoto

#_______________________________________________________________________________________________________________________________________________________________________________________________________________
#                     | platform     | orig. sensor    | LR sensor       | HR sensor       | fusion type | filter kernel | scene         | size ID       | fDS         | SNR          | redundant digit (>0 for non-regular datasets #
# Test data           |--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#   2114211004350     | 2 (CG-PC)    - 1 (HySpex)      - 1 (HySpex)      - 4 (WorldView-2) - 2 (HS-MS)   - 1 (gauss)     - 1 (MUC_Oly)   - 0 (60x80)     - 04 (fDS=4)  - 35 (SNR=35)  - 0                                            #	
#   2114211108000  <= | 2 (CG-PC)    - 1 (HySpex)      - 1 (HySpex)      - 4 (WorldView-2) - 2 (HS-MS)   - 1 (gauss)     - 1 (MUC_Oly)   - 1 (240x120)   - 08 (fDS=8)  - 00 (SNR=inf) - 0                                            #
#   2115211105000  <= | 2 (CG-PC)    - 1 (HySpex)      - 1 (HySpex)      - 5 (Quickbird)   - 2 (HS-MS)   - 1 (gauss)     - 1 (MUC_Oly)   - 1 (240x120)   - 05 (fDS=5)  - 00 (SNR=inf) - 0                                            #
#   2144111104000  <= | 2 (CG-PC)    - 1 (HySpex)      - 4 (WorldView-2) - 4 (WorldView-2) - 1 (MS-Pan)  - 1 (gauss)     - 1 (MUC_Oly)   - 1 (240x120)   - 04 (fDS=4)  - 00 (SNR=inf) - 0                                            # 
# JSpFI TGRS paper
#
#   2144111210000  <= | 2 (CG-PC)    - 1 (HySpex)      - 4 (WorldView-2) - 4 (WorldView-2) - 1 (MS-Pan)  - 1 (gauss)     - 1 (MUC_Oly)   - 2 (3600x1200) - 10 (fDS=10) - 00 (SNR=inf) - 0                                            #
#
#   2144111210100   <= | 2 (CG-PC)    - 1 (HySpex)      - 4 (WorldView-2) - 4 (WorldView-2) - 1 (MS-Pan)  - 1 (gauss)     - 1 (MUC_Oly)   - 2 (3600x1200) - 10 (fDS=10) - 10 (SNR=10)  - 0                                            #
#   2144111210150   <= | 2 (CG-PC)    - 1 (HySpex)      - 4 (WorldView-2) - 4 (WorldView-2) - 1 (MS-Pan)  - 1 (gauss)     - 1 (MUC_Oly)   - 2 (3600x1200) - 10 (fDS=10) - 15 (SNR=15)  - 0                                            #
#   2144111210200   <= | 2 (CG-PC)    - 1 (HySpex)      - 4 (WorldView-2) - 4 (WorldView-2) - 1 (MS-Pan)  - 1 (gauss)     - 1 (MUC_Oly)   - 2 (3600x1200) - 10 (fDS=10) - 20 (SNR=20)  - 0                                            #
#   2144111210300   <= | 2 (CG-PC)    - 1 (HySpex)      - 4 (WorldView-2) - 4 (WorldView-2) - 1 (MS-Pan)  - 1 (gauss)     - 1 (MUC_Oly)   - 2 (3600x1200) - 10 (fDS=10) - 30 (SNR=30)  - 0                                            #
#   2144111210400   <= | 2 (CG-PC)    - 1 (HySpex)      - 4 (WorldView-2) - 4 (WorldView-2) - 1 (MS-Pan)  - 1 (gauss)     - 1 (MUC_Oly)   - 2 (3600x1200) - 10 (fDS=10) - 40 (SNR=40)  - 0                                            #
#
#   2144111210002   <= | 2 (CG-PC)    - 1 (HySpex)      - 4 (WorldView-2) - 4 (WorldView-2) - 1 (MS-Pan)  - 1 (gauss)     - 1 (MUC_Oly)   - 2 (3600x1200) - 10 (fDS=10) - 00 (SNR=inf) - 2  (with Pan-band replaced by rec. band 2)   #
#   2144111210006   <= | 2 (CG-PC)    - 1 (HySpex)      - 4 (WorldView-2) - 4 (WorldView-2) - 1 (MS-Pan)  - 1 (gauss)     - 1 (MUC_Oly)   - 2 (3600x1200) - 10 (fDS=10) - 00 (SNR=inf) - 6  (with Pan-band replaced by rec. band 6)   #
#-----------------------
#   2144111204000   <= | 2 (CG-PC)    - 1 (HySpex)      - 4 (WorldView-2) - 4 (WorldView-2) - 1 (MS-Pan)  - 1 (gauss)     - 1 (MUC_Oly)   - 2 (3600x1200) - 04 (fDS=4)  - 00 (SNR=inf) - 0                                            #
#
#   2144111204002   <= | 2 (CG-PC)    - 1 (HySpex)      - 4 (WorldView-2) - 4 (WorldView-2) - 1 (MS-Pan)  - 1 (gauss)     - 1 (MUC_Oly)   - 2 (3600x1200) - 04 (fDS=4)  - 00 (SNR=inf) - 2  (with Pan-band replaced by rec. band 2)   #
#   2144111204006   <= | 2 (CG-PC)    - 1 (HySpex)      - 4 (WorldView-2) - 4 (WorldView-2) - 1 (MS-Pan)  - 1 (gauss)     - 1 (MUC_Oly)   - 2 (3600x1200) - 04 (fDS=4)  - 00 (SNR=inf) - 6  (with Pan-band replaced by rec. band 6)   #
#-----------------------
#   2444111104000   <= | 2 (CG-PC)    - 4 (WorldView-2) - 4 (WorldView-2) - 4 (WorldView-2) - 1 (MS-Pan)  - 1 (gauss)     - 1 (MUC)       - 1 (960x1320)  - 04 (fDS=4)  - 00 (SNR=inf) - 0                                            #
#
#   2444111104002   <= | 2 (CG-PC)    - 4 (WorldView-2) - 4 (WorldView-2) - 4 (WorldView-2) - 1 (MS-Pan)  - 1 (gauss)     - 1 (MUC)       - 1 (960x1320)  - 04 (fDS=4)  - 00 (SNR=inf) - 2 (with Pan-band replaced by rec. band 2)    #
#   2444111104006   <= | 2 (CG-PC)    - 4 (WorldView-2) - 4 (WorldView-2) - 4 (WorldView-2) - 1 (MS-Pan)  - 1 (gauss)     - 1 (MUC)       - 1 (960x1320)  - 04 (fDS=4)  - 00 (SNR=inf) - 6 (with Pan-band replaced by rec. band 6)    #
#-----------------------
# Standard data (151105)
#   2114211504350   <= | 2 (CG-PC)    - 1 (HySpex)      - 1 (HySpex)      - 4 (WorldView-2) - 2 (HS-MS)   - 1 (gauss)     - 1 (MUC_Oly)      - 5 (540x540)   - 04 (fDS=4)  - 35 (SNR=35db)- 0                                            #
#   2114211510350   <= | 2 (CG-PC)    - 1 (HySpex)      - 1 (HySpex)      - 4 (WorldView-2) - 2 (HS-MS)   - 1 (gauss)     - 1 (MUC_Oly)      - 5 (540x540)   - 10 (fDS=10) - 35 (SNR=35db)- 0                                            #
#   2114211515350   <= | 2 (CG-PC)    - 1 (HySpex)      - 1 (HySpex)      - 4 (WorldView-2) - 2 (HS-MS)   - 1 (gauss)     - 1 (MUC_Oly)      - 5 (540x540)   - 15 (fDS=15) - 35 (SNR=35db)- 0                                            #
#   2334211304350   <= | 2 (CG-PC)    - 3 (Aviris)      - 3 (Aviris)      - 4 (WorldView-2) - 2 (HS-MS)   - 1 (gauss)     - 1 (Ind.Pines)    - 3 (360x360)   - 04 (fDS=4)  - 35 (SNR=35db)- 0                                            #
#   2334212404350   <= | 2 (CG-PC)    - 3 (Aviris)      - 3 (Aviris)      - 4 (WorldView-2) - 2 (HS-MS)   - 1 (gauss)     - 2 (Cuprite Sc03) - 4 (420x360)   - 04 (fDS=4)  - 35 (SNR=35db)- 0                                            #
#   2334212405350   <= | 2 (CG-PC)    - 3 (Aviris)      - 3 (Aviris)      - 4 (WorldView-2) - 2 (HS-MS)   - 1 (gauss)     - 2 (Cuprite Sc03) - 4 (420x360)   - 05 (fDS=5)  - 35 (SNR=35db)- 0                                            #
#   2335213304350   <= | 2 (CG-PC)    - 3 (Aviris)      - 3 (Aviris)      - 5 (Quickbird)   - 2 (HS-MS)   - 1 (gauss)     - 3 (Moffett Field)- 3 (360x360)   - 04 (fDS=4)  - 35 (SNR=35db)- 0                                            #
#   2665211108350   <= | 2 (CG-PC)    - 6 (ROSIS)       - 6 (ROSIS)       - 5 (Quickbird)   - 2 (HS-MS)   - 1 (gauss)     - 1 (Pavia Univ.)  - 1 (560x320)   - 08 (fDS=8)  - 35 (SNR=35db)- 0                                            #
#   2665211108990   <= | 2 (CG-PC)    - 6 (ROSIS)       - 6 (ROSIS)       - 5 (Quickbird)   - 2 (HS-MS)   - 1 (gauss)     - 1 (Pavia Univ.)  - 1 (560x320)   - 08 (fDS=8)  - 99 (SNR=inf) - 0                                            #
#   2774211106350   <= | 2 (CG-PC)    - 7 (Headwall)    - 7 (Headwall)    - 4 (WorldView-2) - 2 (HS-MS)   - 1 (gauss)     - 1 (Chikusei)     - 1 (540x420)   - 06 (fDS=6)  - 35 (SNR=35db)- 0                                            #
#   2774212106350   <= | 2 (CG-PC)    - 7 (Headwall)    - 7 (Headwall)    - 4 (WorldView-2) - 2 (HS-MS)   - 1 (gauss)     - 2 (Chikusei n.u.)- 1 (540x420)   - 06 (fDS=6)  - 35 (SNR=35db)- 0                                            #
#   2885211404350   <= | 2 (CG-PC)    - 8 (HYDICE)      - 8 (HYDICE)      - 5 (Quickbird)   - 2 (HS-MS)   - 1 (gauss)     - 1 (Wash.DC Mall) - 4 (420x300)   - 04 (fDS=4)  - 35 (SNR=35db)- 0                                            #
#   22109211103990  <= | 2 (CG-PC)    - 2 (HyMap)       - 10 (EnMAP)      - 9 (Sentinel2)   - 2 (HS-MS)   - 1 (gauss)     - 1 (Rodalquilar)  - 1 (261x867)   - 03 (fDS=3)  - 99 (SNR=inf) - 0                                            #
#   211119211105350 <= | 2 (CG-PC)    - 11 (GRSSDFC2013)- 11 (GRSSDFC2013)- 9 (Sentinel2)   - 2 (HS-MS)   - 1 (gauss)     - 1 (Housten Univ.)- 1 (320x540)   - 05 (fDS=5)  - 35 (SNR=inf) - 0                                            #

#######################  4: regularization parameter trading the weighting of the l_2,1 norm term in the joint sparsity optimization problem (float)
lambda=1000 # 1
#######################  5: number of dictionary atoms/patches (int)
NDP=300 # 900
#######################  6: two_step_estimation [default=0] (bool)
twoStep=0
#######################  7: number of jointly processed channels, i.e. number of channels per bundle / per optimization problem (int)
Nc=1
#######################  8: number of overlapping channels, i.e. number of adjacent channels in which adjacend channel bundles overlap (int)
No=0
#######################  9: patchsize, measured in low resolution pixels (int)
psz=3
#######################  10: overlap, measured in low resolution pixels (int)
ovrlp=0
#######################  11: fusion method (either of SparseFI, JSparseFI, GroupedJSparseFI or JJSparseFIHM) (string)
alg=JJSparseFIHM
#######################  12: evaluate.. (1: evaluation will be done immetiately after the reconstruction during the same job; 0: evaluation will have to be done in a separate after processing step) (bool)
eval=1
#######################  13: relative tolerace [double between 0 and 0.5] for the decision matrix / spectral response functions (float)
tol_SRF=0.7
#######################  14: (bool)
store_patches_tmp_on_drive=0
#######################  15: maximum number of processors to be used for writing final image in parallel. [only affective if parWr==1] (int)
parWrNumProc=1
#######################  16: save coefficient vectors (alphas) (bool)
saveAlphas=0
#######################  17: first patch coefficient vector to be saved #### [applicable ONLY IF (saveAlpha==TRUE)] (int)  
pFirstAlpha=0
#######################  18: last patch coefficient vector to be saved  #### [applicable ONLY IF (saveAlpha==TRUE)] (int)
pLastAlpha=9999999
#######################  19: save local dictionary coordinates (uP,vP) (bool)
saveDicts=0
#######################  20: first patch which of you want the local dictionary coordinates to be saved #### [applicable ONLY IF (saveDicts==TRUE)] (int)
pFirstDict=0
#######################  21: last patch of which you want the local dictionary coordinates to be saved  #### [applicable ONLY IF (saveDicts==TRUE)] (int)
pLastDict=9999999
#######################  22: first Y band to be processed (int)
chBundleFirst=0
#######################  23: last Y band to be processed (int)
chBundleLast=999
#######################  24: low resolution vertival pixel coordinate from where the reconstruction should start (default: 0) (int)
uLFirst=0
#######################  25: low resolution vertival pixel coordinate from where the reconstruction should end (default: 9999999) (int)
uLLast=99999
#######################  26: low resolution horizontal pixel coordinate from where the reconstruction should start (default: 0) (int)
vLFirst=0 # 
#######################  27: low resolution horizontal pixel coordinate from where the reconstruction should end (default: 9999999) (int)
vLLast=99999
#######################  28: number of processes working on one patch. 1 <= numProcPerPatch <= numProbPerPatch, where numProbPerPatch is the number of reconstruction problems per patch (only greater than 1 for hyperspectral image enhancement) (default: 9999999) (int)
numProcPerPatch=1
#######################  29: parallelization strategy (<0 no work stealing, =0 work stealing only for the spare patches, = 1,2,3,4, ...  work stealing for the spare patches plus number of turns (multiples of concurrently processed patches) before end to start work stealing)) (default: 0) (int)
workStealingTurns=-1
#######################  30: Select coupled LR and HR Pan dictionaries according to:
# Dictionary contains ONLY the current patch -> (NDP=1) & Alpha is calculated by least squares    (0)
# Nearest Neighbor                                                                                (1)
# PanLR norm                                                                                      (2)
# SRF approximate PanLR  norm                                                                     (3)
# PanHR norm (POSITIVE PanHR correlation)   (depricated)                                          (4)
# PanLR-PanHR joint ranking                 (depricated)                                          (5)
# ABSOLUTE PanHR correlation                                                                      (6)
# PanHR uncorrelation, including current patch as first atom                                      (7)
# Random, including current patch as first atom                                                   (8)
# PanHR self uncorrelated basis approximation, including current patch as first atom              (9)
dictselect=8
#######################  31: write fused image in file (1: create file and write resulting image in file; 0: to not write image in file (useful for analyses only)) (bool)
writeImageFile=1
#######################  32: delete tmp patch folders after the fusion process? (bool, default=1)
delete_tmp_patch_folders=1
#######################  33: Construct the image exclusively from previously processed TMP patches. No patch processing will be done, which makes most of the parameters set above inaffective (bool)
                             #### TMP patches must be located in the directories given as last arguments in this list
imageConstructionOnly=0
#######################  34: continue unfinished reconstruction (bool)
contUnfinishedRec=0
#######################  35: path to .csv file that contains the list of linear patch IDs (iP numbers) which remain to be processed in oder to finish incomplete reconstruction
                             #### [applicable ONLY IF (contUnfinishedRec==TRUE)] (string)
PathToIncompletePatchSetCSV=path_to_csv_file #./tmp.csv
#######################  36: number of additional directories that contain previously processed tmp patches needed to reconstruct full image (int)
                             #### !! Must be the number of strings that are specified as last program arguments
dir_tmp_patches_additional_num=0
############ TMP >>>>>>
#######################  37: matrix (dictionary) normalization norm: 1: spectral norm, 0: Frobenious norm (bool)  
matrixNorm=1
#######################  38: instead of adding the mean of the LR patches to the (zero-mean) patch reconstruction result, the LR pixel values are added to the corrensponding set of HR pixels  (bool)  
addMeanPixelwise=0
############ <<<<<< TMP 
#######################################################
# J-P-FISTA settings for joint sparse reconstruction  #
#######################################################
#######################  39: maximum number of iterations to run the algorithm (int)    
maxiter_out=200000
#######################  40: tolerance (double)
tol=1e-12
########################################################
# for Least squares post processing on the patch level #
########################################################
#######################  41: flag used to decide whether or not the least square post-minimization is activated (bool) 
LQ_post_opt=0
#######################  42: regularization parameter trading the relative weighting of the high resolution input patch xHR (double) 
lambdaX=DEFINEDBELOW
LIST_lambdaX = 1e3
#######################  43: regularization parameter trading the relative weighting of the low resolution input patch yLR (double) 
lambdaY=DEFINEDBELOW
LIST_lambdaY = 0
#######################  44: maximum number of iterations in the GS step to solve the least squares problem (int) 
maxiter_CGLS=1000
#######################  45: error tolerance (double) 
tol_r_CGLS=1e-12
#######################  46: decides whether or not the coefficients get updates via least squares (bool) 
fix_Alpha=DEFINEDBELOW
LIST_fix_Alpha = 1
#######################  47: decides whether or not the mean values of Z are set to the same mean values as Y (i.e. either delta_m remains the initial zero vector or it gets updated via least squares) (bool) 
fix_delta_m=DEFINEDBELOW
LIST_fix_delta_m = 1
########################################################
# for Least squares post processing on the image level #
########################################################
#######################  48: flag used to decide whether or not the least square post-minimization (of the final image) is activated (bool)
LQ_post_opt_im=1
#######################  49: regularization parameter trading the relative weighting of the high resolution input image I_X (double)
lambdaX_im=DEFINEDBELOW
LIST_lambdaX_im= 1e-5
#######################  50: regularization parameter trading the relative weighting of the low resolution input image I_Y (double)
lambdaY_im=DEFINEDBELOW
LIST_lambdaY_im=1e-5
#######################  51: maximum number of iterations in the GS step to solve the least squares problem on the final image level (int)
maxiter_CGLS_im=1500
#######################  52: error tolerance (double)
tol_r_CGLS_im=1e-8 # 1e-12
#########################################################
# for coefficient estimation
#########################################################
########################  53: use new method for calculating Z based on step-wise least squares optimization [added in October 2015] (bool)
useNewMethodForCalculatingZ=1
########################  54: use simulated high resolution image X for dictionary learning [added in October 2015] (bool)
useSimulatedImXforDictLearn=1
########################  55~57: regularization parameters for new coefficient estimation [added in October 2015] (double)
lambdaX_ABC=1e-5 # 316
lambdaY_ABC=1e-5 # 0.316
lambdaZ_ABC=1
####################### 58 set lambdaZ_ABC to this number only in the first iteration. A low value can be helpful e.g. if the initial image ImZ_init is not very good/trustworthy (double)
lambdaZ_ABC_in_1st_iter=1 #Set_to_1
####################### 59  type of initial high resolution image ImZ_init; flag (int) 
#ImZ_init_type=0 # lambdaZ_ABZ=0 in 1st iter (no initial image)
#ImZ_init_type=1 # upsampled and bilinearly interpolated low resolution image ImY
ImZ_init_type=2 # reconstruction result of another algorithm (e.g. HySure, Bayesian Sparse, CNMF or MAPSMM. Depends on dataset)
#ImZ_init_type=3  # reference image (for tests only)
####################### 60 jump to the full image optimization (of the initial image) without doing the patch-wise imge reconstruction. (bool)
doFullImOptWithoutPatRec=0
####################### 61 number of coupled ImZ calculations/iterations (int)
iterMain=2
####################### 62 maximum size of spectral group above groups will be double-checked and perhaps split into subgoups (int) 
Nc_max=25
####################### 63 minimum cross-correlation within spectral goups (double)
CC_min=0.96
####################### 64 size of window around patch: Must be odd if patchsize is odd and even if patchsize is even, in order to have both centers matched; Used for correlation calculations (int)
winSize=6
####################### 65 evaluate initial image (bool)
evaluate_ImZ_init=1
####################### 66 set negative values to zero (int)
# set_neg_to_0=0    # set negative values to zero only at the very end, before writing the final image
# set_neg_to_0=1    # set negative values to zero only after patch reconstruction
# set_neg_to_0=2    # set negative values to zero only after full image optimization
set_neg_to_0=3    # set negative values to zero both after patch reconstruction and after full image optimization
####################### 67 use estimated SRFs instead of apriori given ones (bool)
use_estimated_SRFs_LIST=1
####################### 68: write all intermediate image fusion resulta (after every iteration) (1: create file and write resulting image in file; 0: to not write image in file (useful for analyses only)) (bool)
writeImageFileAfterEveryIter=1
####################### 69: 0: do full image optimization on the full image space (not on a subspace); 1: in the full image optimization step, operate on a subspace (optimize for A, where  ImZ = E*A and E is known) (only effenctive if LQ_post_opt_im is set to 1 above)  (bool)
fullImOptOnSubspace_LIST=1
####################### 70: e.g. =1 for SuperMUC and =2 for CG local CP (int)
platformID=2
####################### 71: subspace transformation type
subspace_transform_type_LIST=SVD
#                           =PCA
#                           =SVD
#                           =VCA
#                           =none
####################### 6: subspace dimension
subspace_dim_LIST=10
####################### 73: ImX simulation mode: 0: correlation based; 1: unconstrained LS based; 2: NNLS based
ImX_sim_mode=2
####################### 74: calc. coeff. in full image opt. eq. via SNR calc. of ImX and ImY (bool)
SNR_normalization_LIST=1
####################### 75: set the coeff. of |R*Z-X| in full image opt. eq. to NChY/NChX [only relevant if SNR_normalization==1 ] (bool)
balance_ImX_term_coef_LIST=0
####################### 76: save output in double format (64bit) instead of uint16 (bool)
saveAsDouble=1
####################### 77: use LR (low resolution) patch norm for normalization of corresponding LR and HR patch in coupled dictionaries. If set to 0 the HR nor is used by default. (bool)
use_LRnorm_for_dic_normalization=1
####################### 78: load and use a-priori calculated dictionaries. (bool)
load_DictHR_and_DictLR=0

#######################  LAST ARGUMENT(S): [optional for reparation purposes] additional directories that contain possibly inclomplete sets of previously processed tmp patches possibly needed to reconstruct full image (string(s))
                             #### !!! This list of strings must be at least as long as specified by the argument 'dir_tmp_patches_additional_num'  !!!
                             ####     These MUST be the very last arguments !!!                             
####dir_tmp_patches_additional_1=path_of_first_directory_to_be_searched_if_patch_not_found
dir_tmp_patches_additional_1=/users/groh_cl/supermuc/gpfs/work/ga39yoz2/tmp/patches/150710_150814_846510
dir_tmp_patches_additional_2=path_of_second_directory_to_be_searched_if_patch_not_found
dir_tmp_patches_additional_3=path_of_third_directory_to_be_searched_if_patch_not_found

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
		$(foreach datasetID, $(LIST_datasetID), \
		$(foreach fix_delta_m, $(LIST_fix_delta_m), \
		$(foreach lambdaX_im, $(LIST_lambdaX_im), \
		$(foreach lambdaY_im, $(LIST_lambdaY_im), \
		$(foreach use_estimated_SRFs, $(use_estimated_SRFs_LIST), \
		$(foreach fullImOptOnSubspace, $(fullImOptOnSubspace_LIST), \
		$(foreach subspace_transform_type, $(subspace_transform_type_LIST), \
		$(foreach subspace_dim, $(subspace_dim_LIST), \
		$(foreach balance_ImX_term_coef, $(balance_ImX_term_coef_LIST), \
		$(foreach SNR_normalization, $(SNR_normalization_LIST), \
		$(RUN) $(RUNFLAGS) -x LD_LIBRARY_PATH=$(GDAL_LIBRARY_PATH) $(EXE) $(THIS_JOB_NAME) $(THIS_JOB_ID) \
		$(datasetID) $(lambda) $(NDP) $(twoStep) $(Nc) $(No) $(psz) $(ovrlp) $(alg) $(eval) \
		$(tol_SRF) $(store_patches_tmp_on_drive) $(parWrNumProc) $(saveAlphas) $(pFirstAlpha) $(pLastAlpha) \
		$(saveDicts) $(pFirstDict) $(pLastDict) $(chBundleFirst) $(chBundleLast) \
		$(uLFirst) $(uLLast) $(vLFirst) $(vLLast) \
		$(numProcPerPatch) $(workStealingTurns) $(dictselect) \
		$(writeImageFile) $(delete_tmp_patch_folders) $(imageConstructionOnly) $(contUnfinishedRec) \
		$(PathToIncompletePatchSetCSV) $(dir_tmp_patches_additional_num) $(matrixNorm) $(addMeanPixelwise) \
		$(maxiter_out) $(tol) $(LQ_post_opt) $(lambdaX) $(lambdaY) $(maxiter_CGLS) $(tol_r_CGLS) $(fix_Alpha) $(fix_delta_m) \
		$(LQ_post_opt_im) $(lambdaX_im) $(lambdaY_im) $(maxiter_CGLS_im) $(tol_r_CGLS_im) \
		$(useNewMethodForCalculatingZ) $(useSimulatedImXforDictLearn) $(lambdaX_ABC) $(lambdaY_ABC) $(lambdaZ_ABC) $(lambdaZ_ABC_in_1st_iter) $(ImZ_init_type) \
		$(doFullImOptWithoutPatRec) $(iterMain) $(Nc_max) $(CC_min) $(winSize) \
		$(evaluate_ImZ_init) $(set_neg_to_0) $(use_estimated_SRFs) $(writeImageFileAfterEveryIter) $(fullImOptOnSubspace) \
		$(platformID) $(subspace_transform_type) $(subspace_dim) $(ImX_sim_mode) $(SNR_normalization) $(balance_ImX_term_coef) $(saveAsDouble) $(use_LRnorm_for_dic_normalization) $(load_DictHR_and_DictLR)\
		$(dir_tmp_patches_additional_1) $(dir_tmp_patches_additional_2) $(dir_tmp_patches_additional_3);))))))))))

clean:
	rm -f $(OBJ) $(EXE) $(SUB).tar.gz *~ Depends

submit:
	tar czvf $(SUB)$(shell date +_%Y_%m_%d).tar.gz $(EXE) $(SRCMAIN) $(SRC) $(HDR) Makefile 

.SUFFIXES : .o .cpp

Depends depend:
	$(CXX) -E -MM $(INCFLAGS) $(SRCMAIN) $(SRC) > Depends
