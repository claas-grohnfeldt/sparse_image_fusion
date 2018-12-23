/*
 * userSettings.cpp
 *
 *  Created on: Jan 20, 2014
 *      Author: Grohnfeldt, Claas
 */

#include "userSettings.h"


void getUserSettings(SpEOPaths *paths, SpEODataIOSetting *dSetting, SpEOFusionSetting *fSetting, SpEOOutputSetting *oSetting, SpEOSolverSetting *sSetting, SpEOParallelSetting *pSetting, int argc, char **argv){
    
	/*=================================================================*
	 *                          User Settings                          *
	 *=================================================================*/

	int my_processes, my_rank;
	MPI_Comm_size(MPI_COMM_WORLD, &my_processes);
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

	// MPI / parallelization settings
	// total number of processes
	pSetting->numProcTot = my_processes;

	//----------------
	// 1: Arbitrary string roughly describing this job (string)
	//----------------
	dSetting->jobName = argv[1]; // job_name

	//----------------------------------
	// paths
	//----------------------------------

	//----------------
	// 2: Filename of image with higher spatial resolution (string) 
	//----------------
	paths->fname_ImX = argv[2]; // filename_ImX

	//----------------
	// 3: Filename of image with lower spatial resolution higher spectral 
	//    resolution (string)
	//----------------
	paths->fname_ImY = argv[3]; // filename_ImY

	//----------------
	// 4: Filename of initial image (string)
	//    (The final fusion product is the result of an alternating opti-
	//    mization process, which requires initialization. A simple initia-
	//    lization can be optained by bilinearly interpolating ImY - the
	//    lower-resolution input image)
	//----------------
	paths->fname_ImZ_init = argv[4]; // filename_ImZ_init

	//----------------
	// 5: Boolean flag to specify whether or not a high-resolution 
	//    reference ("ground truth") image (ImZ_ref) is available. (bool)
	//----------------
	fSetting->ImZ_ref_avlbl = (bool)atoi(argv[5]); // reference_image_available

	//----------------
	// 6: Filename of image reference ("ground truth") image. (bool)
	//----------------
	paths->fname_ImZ_ref = argv[6]; // filename_ImZ_ref

	//----------------
	// 7: Boolean flag to specify whether the program should use SRFs that 
	//    are (1) estimated directly from the data, or (0) known from the 
	//    sensors' specifications a priori. In the latter case, those known 
	//    SRFs must be loaded from an existing .csv file (see next program 
	//    argument below). (bool)
	//----------------
	fSetting->use_estimated_SRFs = (bool)atoi(argv[7]); // estimate_SRFs_from_data

	//----------------
	// 8: File name of .csv file containing a priory known SRFs. Only used
	//    (required) if the flag 'estimage_SRFs_from_data' is set to '0'.
	//    (bool)
	//----------------
	paths->fname_SRF = argv[8]; // filename_SRF

	//----------------
	// 9: Arbitrary directory into which the program's output including 
	//    image fusion results will be stored.
	//----------------
	paths->dir_out = argv[9]; // dirname_output

	paths->fname_ImZ_out = "JSparesFIHM_fusion_result";

	//----------------------------------
	// local-non-local processing module
	//----------------------------------

	//----------------
	// 10: patch size, measured in low resolution pixels. Images patches 
	//     are set to be square and contain patch_size^2 pixels at the low-
	//     resolution scale. (int)
	//----------------
	fSetting->patchsize = atoi(argv[10]); // patch_size

	//----------------
	// 11: Patch overlap, measured in low resolution pixels.
	//     patch_overlap can be set to any number between 0 and 
	//     patch_size-1. Usually, the higher the patch overlap, the better 
	//     the fusion results. (int)
	//----------------
	fSetting->overlap = atoi(argv[11]); // patch_overlap

	//----------------
	// 12: Number of dictionary atoms/patches (int)
	//----------------
	fSetting->NDP = atoi(argv[12]); // N_a

	//----------------
	// 13: Select coupled LR and HR Pan dictionaries according to:
	// options: '0': Dictionary contains ONLY the current patch 
	//               Hence, (N_a=1) & Alpha is calculated by least squares
	//          '1': Nearest Neighbors
	//          '2': PanLR norm
	//          '3': SRF approximate PanLR  norm
	//          '4': PanHR norm (POSITIVE PanHR correlation)   (depricated)
	//          '5': PanLR-PanHR joint ranking                 (depricated)
	//          '6': ABSOLUTE PanHR correlation
	//          '7': PanHR uncorrelation, including current patch as first 
	//               atom
	//          '8': Random, including current patch as first atom
	//          '9': PanHR self uncorrelated basis approximation, including
	//               current patch as first atom
	//----------------
	fSetting->dictselect = atoi(argv[13]); // dictselect

	//----------------
	// 14: Regularization parameter weighting the l_2,1 norm term in the 
	//     joint sparsity optimization problem. (float)
	//----------------
	fSetting->lambda = atof(argv[14]); // lambda

	//----------------
	// 15: minimum cross-correlation within spectral goups (double)
	//----------------
	fSetting->theta = (double)atof(argv[15]); // theta

	//----------------
	// 16: maximum size of spectral group above groups will be double-checked and perhaps split into subgoups (int) 
	//----------------
	fSetting->Nc_max = atoi(argv[16]); // N_c

	//----------------
	// 17-18: regularization parameters for weighing the impact of ImX and
	//        ImY, respectively, to the product (double)
	//----------------
	fSetting->lambdaX_ABC = (double)atof(argv[17]); // mu_X
	fSetting->lambdaY_ABC = (double)atof(argv[18]); // mu_Y

	//----------------
	// 19: ImX simulation mode flag (int) 
	// options: '0': correlation based
	//          '1': unconstrained least-squares based 
	//          '2': non-negative least-squares based
	//----------------
	fSetting->ImX_sim_mode = atoi(argv[19]); // ImX_sim_mode

	//----------------
	// 20: Size of window around patch: For symmetry reasons, it should be
	//     be set to an odd number if patch_size is odd and to en even if 
	//     patch_size is even. Should be set to a number larger than or
	//     equal to patch_size. (int)
	//----------------
	fSetting->winSize = atoi(argv[20]); // winSize

	//----------------------------------
	// Global-non-local processing module
	// (full image optimization)
	//----------------------------------

	//----------------
	// 21-22: Regularization parameters trading the weighting of the high-
	//        and low-resolution input images, I_X and I_Y, respectively.
	//        (double)
	//----------------
	fSetting->lambdaX_im = (double)atof(argv[21]); //mu_X_prime
	fSetting->lambdaY_im = (double)atof(argv[22]); // mu_Y_prime

	//----------------
	// 23: Subspace dimension
	//----------------
	fSetting->subspace_dim = atoi(argv[23]); // subspace_dim

	//----------------------------------
	// Concerning both local-non-local and global processing modules
	//----------------------------------

	//----------------
	// 24: Processing module selection flag (int)
	// options: '0': run only local-non-local processing module
	//          '1': run only global processing module
	//          '2': run both processing module in an alternating manner
	//               (recommended)
	//----------------
	int processing_module_flag = atoi(argv[24]);
	switch(processing_module_flag){
		case 0: {
			fSetting->doFullImOptWithoutPatRec=false;
			fSetting->LQ_post_opt_im = false;
			break;
		}case 1: {
			fSetting->doFullImOptWithoutPatRec=true;
			fSetting->LQ_post_opt_im = true;
			break;
		}case 2: {
			fSetting->doFullImOptWithoutPatRec=false;
			fSetting->LQ_post_opt_im = true;
			break;
		}
	}
	//----------------
	// 25: Maximum number of times each of the two processing module should
	//     be run. More specifically, maximum number of iteration in the 
	//     alternating calculation of ImZ using the local-non-local and 
	//     global processing mdules. (int) 
	//----------------
	fSetting->iterMain = atoi(argv[25]); // iterMain

	//----------------------------------
	// Output settings
	//----------------------------------

	//----------------
	// 26: Boolean flag to specify whether the quality of the fusion result
	//     (image ImZ) should be assessed relative to the supposedly 
	//     provided high-resolution reference ("ground truth") image, i.e.
	//     ImZ_ref. (bool)
	//----------------
	fSetting->evaluate = (bool)atoi(argv[26]); // perform_quality_assessment_of_fusion_result

	fSetting->evaluate_ImZ_init = fSetting->evaluate;

	//----------------
	// 27: write fused image in file (1: create file and write resulting image in file; 0: to not write image in file (useful for analyses only)) (bool)
	//----------------
	oSetting->writeImageFile = (bool)atoi(argv[27]); // writeImageFile

	//----------------
	// write all intermediate image fusion results/images (after every iteration) (bool) 
	// options:
	// 0: do not write image in file (useful for analyses only)
	// 1: create file and write resulting image in file
	// CAN BE PUT BACK TO PROGRAM ARGUMENTS IF NEEDED 
	//----------------
	oSetting->writeImageFileAfterEveryIter = false; // writeImageFileAfterEveryIter

	//----------------
	// 28: Output image format flag.
	// options: '0': 16-bit unsigned integer (UInt16) 
	//          '1': 64-bit float (Float64)
	//----------------
	dSetting->saveAsDouble = (bool)atoi(argv[28]);// output_image_format_flag

	//----------------
	// DEPRICATED
	// 29: scale the coefficient of |R*Z-X| term by a factor of NChY/NChX 
	//     (experimenatal) (bool)
	fSetting->balance_ImX_term_coef = (bool)atoi(argv[29]); // balance_ImX_term_coef
	// ---------------




	//########################################################
    //**********************************
    // default parameters
    // advanced: change with caution!
    //**********************************
	//########################################################

	// use simulated high resolution image X for dictionary learning (bool)
	//fSetting->useSimulatedImXforDictLearn = true; //(bool)atoi(argv[6]);
	// number of processes to work on one sparse reconstruction problem (solver parallelization)
	pSetting->numProcGrp          = 1;
	// Choose if input images should be normalized before calculation
	fSetting->nrmlIm              = false; // to be set to false
	// Choose if the columns of the dictionaries should be normalized for calc.
	fSetting->nrmlDicts           = true;
	// Choose if the dictionary's columns should have zero means.
	fSetting->substrMean          = true;
	// Specify whether or not original high resolution MS image is available.
	// If so, the fusion results will be evaluated.
	oSetting->prec                = 10;
	fSetting->Nc                  = 1; // relevant for older version of J-SparseFI-HM
	fSetting->No                  = 0; // relevant for older version of J-SparseFI-HM
	fSetting->tol_SRF             = 0.5; // relevant for older version of J-SparseFI-HM
	fSetting->two_step_estimation = false;
	pSetting->store_patches_tmp_on_drive = false;
	pSetting->parWrNumProc       = 1; 
	// output settings
	oSetting->saveAlphas         = false;
	oSetting->pFirstAlpha        = 0;
  	oSetting->pLastAlpha         = 999999999;
	oSetting->saveDicts          = false;
	oSetting->pFirstDict         = 0;
  	oSetting->pLastDict          = 999999999;
  	dSetting->chBundleFirst      = 0;
  	dSetting->chBundleLast       = 999999999;
	dSetting->uLFirst            = 0;
	dSetting->uLLast             = 999999999;
	dSetting->vLFirst            = 0;
	dSetting->vLLast             = 999999999;

	pSetting->numProcPerPatch    = 1;
	// patch parallelization strategy: work stealing (true) or fixed work scheduling (false)
	pSetting->workStealingTurns  = -1;

	// calc. coeff. in full image opt. eq. via SNR calc. of ImX and ImY (bool).
	// include SNR normalization to compensate for colored (band-dependend) noise (bool)
	fSetting->SNR_normalization = true;
//  matrix (dictionary) normalization norm: 0: spectral norm, 1: Frobenious norm (bool)
    fSetting->matrixNorm        = true;
    fSetting->addMeanPixelwise  = false;
	//#######################################################
	//# J-P-FISTA settings for joint sparse reconstruction  #
    //#######################################################
	// Choose a solver for the sparse reconstruction (optimization) problem (currently JPFISTA only)
    sSetting->solver             = JPFISTA;
    // maximum number of iterations to run the algorithm (int)
	sSetting->maxiter_out        = 200000;
	// tolerance
	sSetting->tol                = 0.000000000001;
	// decides whether or not the mean values of Z are set to the same mean values as Y (i.e. either delta_m remains the initial zero vector or it gets updated via least squares) (bool)
	sSetting->fix_delta_m  = false;

    //set lambdaZ_ABC to this number only in the first iteration. A low value can be helpful e.g. if the initial image ImZ_init is not very good/trustworthy (double)
	fSetting->lambdaZ_ABC_in_1st_iter = 1.0;

    // use LR (low resolution) patch norm for normalization of corresponding LR and HR patch in coupled dictionaries. If set to 0 the HR nor is used by default. (bool)
	fSetting->use_LRnorm_for_dic_normalization = true; //
	
	fSetting->lambdaZ_ABC = 1.0;

	// set negative values to zero (int)
	fSetting->set_neg_to_0 = 3;
    //			               = 0  => set negative values to zero only at the very end, before writing the final image
    //	                       = 1  => set negative values to zero only after patch reconstruction
    //	                       = 2  => set negative values to zero only after full image optimization
    //	                       = 3  => set negative values to zero both after patch reconstruction and after full image optimization

	// maximum number of iterations in the GS step to solve the least squares problem on the final image level (int)
	sSetting->maxiter_CGLS_im = 1500; //atoi(argv[19]);
	// error tolerance (double)
	sSetting->tol_r_CGLS_im   = (double)atof("1e-12");

	// write all intermediate image fusion results (after every iteration) (1: create file and write resulting image in file; 0: to not write image in file (useful for analyses only)) (bool)
	fSetting->fullImOptOnSubspace = true; // (bool)atoi(argv[21]);
    // subspace transformation type
    // options:
    // - PCA
    // - SVD
    // - VCA
    // - none
    fSetting->subspace_transform_type = "SVD"; //argv[22];

	//##############################################//
	//##         Check & Correct Input            ##//
	//##############################################//

	if(fSetting->dictselect == 0){
		fSetting->NDP = 1;
		fSetting->fMethod = LeastSquares;
	}
	if(oSetting->writeImageFile && pSetting->store_patches_tmp_on_drive && pSetting->parWrNumProc>my_processes){
		if(my_rank == 0){
			cout << endl << "WARNING: number of processors used for parallel writing (parWrNumProc) was set to a number larger than the total number of processors! It got corrected from " << pSetting->parWrNumProc << " to " << my_processes << "!"  << endl;
		}
		pSetting->parWrNumProc=my_processes;
	}
	if(oSetting->writeImageFile && pSetting->store_patches_tmp_on_drive && pSetting->parWrNumProc<1){
		int tmp = min(my_processes,100);
		if(my_rank == 0){
			cout << endl << "WARNING: number of processors used for parallel writing (parWrNumProc) has to be set to a positive number! It got corrected from " << pSetting->parWrNumProc << " to "<< tmp << "!"  << endl;
		}
		pSetting->parWrNumProc = tmp;
	}
}

