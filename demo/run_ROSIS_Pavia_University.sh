#!/bin/bash

######################################################################
# Shell script for running the program "sparse_image_fusion" for 
# multi-sensor image super-resolution, which incoporates an HPC 
# implementation of the algorithm J-SparseFI-HM, which is published in
# "C. Grohnfeldt, 'Multi-sensor Data Fusion for Multi- and Hyper-
# spectral Resolution Enhancement Based on Sparse Representations', 
# Ph.D. Dissertation, Technical University of Munich, 2017; 
# doi:10.14459/2017md1366768"
#
# Author: Claas Grohnfeldt
######################################################################

RUN='mpiexec'
RUNFLAGS="-n $1"

LIBRARY_PATH_GDAL='./lib/gdal/lib'

EXE='./bin/sparse_image_fusion'

#---------------------------------------------------
#
#                program arguments
#
#---------------------------------------------------

#----------------
# 1: Arbitrary string roughly describing this job (string)
#----------------
job_name='JSparseFIHM_ROSIS_Pavia_University'

#----------------------------------
# paths
#----------------------------------

DIR_DATA_BASE='./data/HS_MS/ROSIS_Pavia_University'

#----------------
# 2: Filename of image with higher spatial resolution (string) 
#----------------
filename_ImX="${DIR_DATA_BASE}/InputData/multispectral_highRes/ROSIS_\
Pavia_University_synthesized_QuickBird_MSHR_SNR35.dat"

#----------------
# 3: Filename of image with lower spatial resolution higher spectral 
#    resolution (string)
#----------------
filename_ImY="${DIR_DATA_BASE}/InputData/hyperspectral_lowRes/ROSIS_P\
avia_University_synthesized_HSLR_fDS8_SNR35.dat"

#----------------
# 4: Filename of initial image (string)
#    (The final fusion product is the result of an alternating opti-
#    mization process, which requires initialization. A simple initia-
#    lization can be optained by bilinearly interpolating ImY - the
#    lower-resolution input image)
#----------------
filename_ImZ_init="${DIR_DATA_BASE}/FusionResults/CNMF/ROSIS_Pavia_Un\
iveristy_FusionResult_CNMF.dat"

#----------------
# 5: Boolean flag to specify whether or not a high-resolution 
#    reference ("ground truth") image (ImZ_ref) is available. (bool)
#----------------
reference_image_available='1'

#----------------
# 6: Filename of image reference ("ground truth") image. (bool)
#----------------
filename_ImZ_ref="${DIR_DATA_BASE}/ReferenceData/ROSIS_Pavia_Univeris\
ty_hyperspectral_highRes_reference.dat"
#----------------
# 7: Boolean flag to specify whether the program should use SRFs that 
#    are (1) estimated directly from the data, or (0) known from the 
#    sensors' specifications a priori. In the latter case, those known 
#    SRFs must be loaded from an existing .csv file (see next program 
#    argument below). (bool)
#----------------
estimate_SRFs_from_data='1'

#----------------
# 8: File name of .csv file containing a priory known SRFs. Only used
#    (required) if the flag 'estimage_SRFs_from_data' is set to '0'.
#    (bool)
#----------------
filename_SRF="${DIR_DATA_BASE}/InputData/SRFs/SRFs_of_QuickBird_sampl\
ed_to_centers_of_SRFs_of_ROSIS_for_scene_Pavia_Univeristy.csv"

#----------------
# 9: Arbitrary directory into which the program's output including 
#    image fusion results will be stored.
#----------------
dirname_output="${DIR_DATA_BASE}/FusionResults/JSparseFIHM"

#----------------------------------
# local-non-local processing module
#----------------------------------

#----------------
# 10: patch size, measured in low resolution pixels. Images patches 
#     are set to be square and contain patch_size^2 pixels at the low-
#     resolution scale. (int)
#----------------
patch_size='3'

#----------------
# 11: Patch overlap, measured in low resolution pixels.
#     patch_overlap can be set to any number between 0 and 
#     patch_size-1. Usually, the higher the patch overlap, the better 
#     the fusion results. (int)
#----------------
patch_overlap='2'

#----------------
# 12: Number of dictionary atoms/patches (int)
#----------------
N_a='900'

#----------------
# 13: Select coupled LR and HR Pan dictionaries according to:
# options: '0': Dictionary contains ONLY the current patch 
#               Hence, (N_a=1) & Alpha is calculated by least squares
#          '1': Nearest Neighbors
#          '2': PanLR norm
#          '3': SRF approximate PanLR  norm
#          '4': PanHR norm (POSITIVE PanHR correlation)   (depricated)
#          '5': PanLR-PanHR joint ranking                 (depricated)
#          '6': ABSOLUTE PanHR correlation
#          '7': PanHR uncorrelation, including current patch as first 
#               atom
#          '8': Random, including current patch as first atom
#          '9': PanHR self uncorrelated basis approximation, including
#               current patch as first atom
#----------------
dictionary_selection_method='8'

#----------------
# 14: Regularization parameter weighting the l_2,1 norm term in the 
#     joint sparsity optimization problem. (float)
#----------------
lambda='1e0'

#----------------
# 15: Minimum cross-correlation value within spectral groups. (double)
#----------------
theta='0.96'

#----------------
# 16: Spectral group size threshold above which spectral channel
#     groups will be double-checked and possibly split into channel
#     sub-groups. (int)
#----------------
N_c='25'

#----------------
# 17-18: regularization parameters for weighing the impact of ImX and
#        ImY, respectively, to the product (double)
#----------------
mu_X='3.16'
mu_Y='0.316'

#----------------
# 19: ImX simulation mode flag (int) 
# options: '0': correlation based
#          '1': unconstrained least-squares based 
#          '2': non-negative least-squares based
#----------------
ImX_sim_mode='1'

#----------------
# 20: Size of window around patch: For symmetry reasons, it should be
#     be set to an odd number if patch_size is odd and to en even if 
#     patch_size is even. Should be set to a number larger than or
#     equal to patch_size. (int)
#----------------
winSize='9'

#----------------------------------
# Global-non-local processing module
# (full image optimization)
#----------------------------------

#----------------
# 21-22: Regularization parameters trading the weighting of the high-
#        and low-resolution input images, I_X and I_Y, respectively.
#        (double)
#----------------
mu_X_prime='1.0'
mu_Y_prime='1.0'

#----------------
# 23: Subspace dimension
#----------------
subspace_dim='10'

#----------------------------------
# Concerning both local-non-local and global processing modules
#----------------------------------

#----------------
# 24: Processing module selection flag (int)
# options: '0': run only local-non-local processing module
#          '1': run only global processing module
#          '2': run both processing module in an alternating manner
#               (recommended)
#----------------
processing_module_flag='2'

#----------------
# 25: Maximum number of times each of the two processing module should
#     be run. More specifically, maximum number of iteration in the 
#     alternating calculation of ImZ using the local-non-local and 
#     global processing mdules. (int) 
#----------------
iterMain='1'

#----------------------------------
# Output settings
#----------------------------------

#----------------
# 26: Boolean flag to specify whether the quality of the fusion result
#     (image ImZ) should be assessed relative to the supposedly 
#     provided high-resolution reference ("ground truth") image, i.e.
#     ImZ_ref. (bool)
#----------------
perform_quality_assessment_of_fusion_result='1'

#----------------
# 27: Boolean flag to specify whether or not the fusion result (super-
#     resolved image) should bewritten to a file. (bool)
#----------------
writeImageFile='1'

#----------------
# 28: Output image format flag.
# options: '0': 16-bit unsigned integer (UInt16) 
#          '1': 64-bit float (Float64)
#----------------
output_image_format_flag='1'

#----------------
# DEPRICATED
# 29: scale the coefficient of |R*Z-X| term by a factor of NChY/NChX 
#     (experimenatal) (bool)
balance_ImX_term_coef='0'
# ---------------

#---------------------------------------------------------------------
# User settings end here. Don't modify beyond this line.
#---------------------------------------------------------------------

PROG_ARGS="$job_name $filename_ImX $filename_ImY $filename_ImZ_init \
$reference_image_available $filename_ImZ_ref \
$estimate_SRFs_from_data $filename_SRF $dirname_output $patch_size \
$patch_overlap $N_a $dictionary_selection_method $lambda $theta $N_c \
$mu_X $mu_Y $ImX_sim_mode $winSize $mu_X_prime $mu_Y_prime \
$subspace_dim $processing_module_flag $iterMain \
$perform_quality_assessment_of_fusion_result $writeImageFile \
$output_image_format_flag $balance_ImX_term_coef"

#----------------------------------
# Run program
#----------------------------------
export LD_LIBRARY_PATH=$LIBRARY_PATH_GDAL
$RUN $RUNFLAGS $EXE $PROG_ARGS