############################################################################
# Makefile for building the binary file for: SparseFI, J-SparseFI or J-SparseFI-HM
# Author: Claas Grohnfeldt
############################################################################

### exec
RUN="mpiexec"
RUNFLAGS="-n 6"

EXE="./bin/JSparseFIHM"

LIBRARY_PATH_GDAL="./lib/gdal/lib"

########################################
#                                      #
#          program arguments           #
#                                      #
########################################

# 1: job name (string roughly describing this job)
JOB_NAME="SparseFIHM_demo_665211108350_Pavia_University"
#       =LOADL_JOB_NAME

# OBSOLETE
# 2: job ID (6 digit integer - on the SuperMUC this number will be automatically set. On a local machine it is set manually to an arbitrary string)
LOADL_PID=xxxxxx

# OBSOLETE
# 3: dataset ID (string) 
datasetID=665211108350
# examples (see src/paths.cpp):

# 4: fusion method (string)
alg=JSparseFIHM
# - JSparseFIHM
# - GroupedJSparseFI
# - SparseFI
# options: 

#########################################################
# local-non-local processing module
#########################################################

# 5: regularization parameter trading the weighting of the l_2,1 norm term in the joint sparsity optimization problem (float)
lambda=1e1

#*************************************
# dictionary generation
#*************************************

# 6: use simulated high resolution image X for dictionary learning (bool)
useSimulatedImXforDictLearn=1

# 7: ImX simulation mode (int) 
# options:
# 0: correlation based
# 1: unconstrained least-squares based 
# 2: non-negative least-squares based
ImX_sim_mode=2

# 8: Select coupled LR and HR Pan dictionaries according to:
# options: 
# 0: Dictionary contains ONLY the current patch -> (N_a=1) & Alpha is calculated by least squares
# 1: Nearest Neighbor
# 2: PanLR norm
# 3: SRF approximate PanLR  norm
# 4: PanHR norm (POSITIVE PanHR correlation)   (depricated)
# 5: PanLR-PanHR joint ranking                 (depricated)
# 6: ABSOLUTE PanHR correlation
# 7: PanHR uncorrelation, including current patch as first atom
# 8: Random, including current patch as first atom
# 9: PanHR self uncorrelated basis approximation, including current patch as first atom
dictselect=8

# 9: number of dictionary atoms/patches (int)
N_a=900

# 10: patchsize, measured in low resolution pixels (int)
psz=4

# 11: patch overlap, measured in low resolution pixels (int)
ovrlp=2
# options:
# any int between 0 and psz-1

#*************************************
# coefficient estimation
#*************************************

# 12~13: regularization parameters (double)
mu_X=1e0
mu_Y=1e0

#*************************************
# correlation-based Hyperspectral Grouping (CorHySpeG)
#*************************************

# 14: maximum size of spectral group above groups will be double-checked and perhaps split into subgoups (int) 
N_c=25

# 15: minimum cross-correlation within spectral goups (double)
theta=0.96

# 16: size of window around patch: Must be odd if patchsize is odd and even if patchsize is even, in order to have both centers matched; Used for correlation calculations (int)
winSize=6

#######################################################
# Global processing module                             
# (for full-image optimization)                        
#######################################################

# 17: regularization parameter trading the relative weighting of the high resolution input image I_X (double)
mu_X_prime=1e0

# 18: regularization parameter trading the relative weighting of the low resolution input image I_Y (double)
mu_Y_prime=1e0

# 19: maximum number of iterations in the GS step to solve the least squares problem on the final image level (int) 
maxiter_globalOpt=1500 

# 20: error tolerance (double) 
tol_r_globalOpt=1e-12 

# 21: do full image optimization on the full image space (not on a subspace); 1: in the full image optimization step, operate on a subspace (optimize for A, where  ImZ = E*A and E is known) (only effenctive if LQ_post_opt_im is set to 1 above)  (bool)
fullImOptOnSubspace=1

# 22: subspace transformation type
# options:
# - PCA
# - SVD
# - VCA
# - none
subspace_transform_type=SVD

# 23: subspace dimension
subspace_dim=10

# 24: include SNR normalization to compensate for colored (band-dependend) noise (bool)
SNR_normalization=1

# 25: scale the coefficient of |R*Z-X| term by a factor of NChY/NChX (experimenatal, depricated) (bool)
balance_ImX_term_coef=1

######################################################### 
# initialization and interplay between local-non-local  
# and global processing module 
######################################################### 

# 26: use estimated SRFs instead of apriori given ones (bool)
use_estimated_SRFs=1

# 27: type of initial high resolution image ImZ_init; flag (int)
# options:
# 0: # mu_Z=0 in 1st iter (no initial image)
# 1: # upsampled and bilinearly interpolated low resolution image ImY
# 2: # reconstruction result of another algorithm (e.g. CNMF or HySure)
# 3: # reference image (for tests only)
ImZ_init_type=2

# 28: use/activate global processing module (bool) 
use_global_proc_module=1

# 29: skip local-non-local processing module, i.e. jump straight to global proceccing module. (bool)
use_ONLY_global_proc_module=0

# 30: number of coupled ImZ calculations, i.e. local-non-local / glocal processing iterations 
iterMain=2

########################################################
# Output settings                                      #
########################################################

# 31: evaluate reconstructed image(s) (bool)
# options:
# 1: evaluation will be done immetiately after the reconstruction 
# 0: evaluation will have to be done in a separate post-processing step)
eval=1

# 32: evaluate initial image (bool)
evaluate_ImZ_init=1

# 33: write fused image in file (1: create file and write resulting image in file; 0: to not write image in file (useful for analyses only)) (bool)
writeImageFile=1

# 34: write all intermediate image fusion results/images (after every iteration) (bool) 
# options:
# 0: do not write image in file (useful for analyses only)
# 1: create file and write resulting image in file
writeImageFileAfterEveryIter=1

# 35: save output in double format (64bit) instead of uint16 (bool)
saveAsDouble=1

# DIR_DATA="./data/HS_MS/demo_small_dataset"
DIR_DATA="./data/HS_MS/ROSIS_Pavia_Univeristy"
# 36:
# fname_ImX="${DIR_DATA}/InputData/multispectral_highRes/ROSIS_Pavia_University_100x80_Sentinel2_10m_bands_MSHR_SNR35.dat"
fname_ImX="${DIR_DATA}/InputData/multispectral_highRes/ROSIS_Pavia_University_synthesized_QuickBird_MSHR_SNR35.dat"
# 37:
# fname_ImY="${DIR_DATA}/InputData/hyperspectral_lowRes/ROSIS_Pavia_University_100x80_HSLR_fDS4_SNR35.dat"
fname_ImY="${DIR_DATA}/InputData/hyperspectral_lowRes/ROSIS_Pavia_University_synthesized_HSLR_fDS8_SNR35.dat"
# 38:
# fname_ImZ_init="${DIR_DATA}/FusionResults/bilinearly_interpolated/ROSIS_Pavia_University_100x80_HSLR_fDS4_US_via_bilinear_SNR35.dat"
fname_ImZ_init="${DIR_DATA}/FusionResults/CNMF/ROSIS_Pavia_Univeristy_FusionResult_CNMF.dat"
# 39:
# fname_ImZ_ref="${DIR_DATA}/ReferenceData/ROSIS_Pavia_Univeristy_hyperspectral_highRes_reference_subarea.dat"
fname_ImZ_ref="${DIR_DATA}/ReferenceData/ROSIS_Pavia_Univeristy_hyperspectral_highRes_reference.dat"
# 40:
fname_SRF="${DIR_DATA}/InputData/SRFs/SRFs_of_QuickBird_sampled_to_centers_of_SRFs_of_ROSIS_for_scene_Pavia_Univeristy.csv"
# 41:
dir_out="${DIR_DATA}/FusionResults/JSparseFIHM"
# 42:
dir_tmp="./tmp"

# ------------------------------------------------------------------------

PROG_ARGS="$JOB_NAME $LOADL_PID $datasetID $alg $lambda $useSimulatedImXforDictLearn $ImX_sim_mode $dictselect $N_a $psz $ovrlp $mu_X $mu_Y $N_c $theta $winSize $mu_X_prime $mu_Y_prime $maxiter_globalOpt $tol_r_globalOpt $fullImOptOnSubspace $subspace_transform_type $subspace_dim $SNR_normalization $balance_ImX_term_coef $use_estimated_SRFs $ImZ_init_type $use_global_proc_module $use_ONLY_global_proc_module $iterMain $eval $evaluate_ImZ_init $writeImageFile $writeImageFileAfterEveryIter $saveAsDouble $fname_ImX $fname_ImY $fname_ImZ_init $fname_ImZ_ref $fname_SRF $dir_out $dir_tmp"

# ------------------------------------------------------------------------

$RUN $RUNFLAGS -x LD_LIBRARY_PATH=$LIBRARY_PATH_GDAL $EXE $PROG_ARGS