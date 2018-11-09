/*
 * userSettings.cpp
 *
 *  Created on: Jan 20, 2014
 *      Author: Grohnfeldt, Claas
 */

#include "userSettings.h"


void getUserSettings(SpEODataIOSetting *dSetting, SpEOFusionSetting *fSetting, SpEOOutputSetting *oSetting, SpEOSolverSetting *sSetting, SpEOParallelSetting *pSetting, int argc, char **argv){

    
	/*=================================================================*
	 *                          User Settings                          *
	 *=================================================================*/

	int my_processes, my_rank;
	MPI_Comm_size(MPI_COMM_WORLD, &my_processes);
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

	// MPI / parallelization settings
	// total number of processes
	pSetting->numProcTot = my_processes;
	dSetting->jobName             = argv[1];
	dSetting->jobID               = argv[2];

	// Choose a sparse image fusion methods
	//fSetting->fMethod             = JSparseFI;
	//                            = SparseFI
	//                         or = JSparseFI
	//                         or = JSparseFIHM
	if(strcmp(argv[4],"JSparseFIHM") == 0){
		fSetting->fMethod             = JSparseFIHM;
	    // use new method for calculating Z based on step-wise least squares optimization (bool)
	    fSetting->useNewMethodForCalculatingZ = true;
	//}else if(strcmp(argv[4],"JSparseFI") == 0){  <- old implementation of J-SparesFI
	//	fSetting->fMethod             = JSparseFI;
	}else if(strcmp(argv[4],"GroupedJSparseFI") == 0){
		fSetting->fMethod             = GroupedJSparseFI; // <- new implementation of J-SparseFI for WorldView-2
	    // use new method for calculating Z based on step-wise least squares optimization (bool)
	    fSetting->useNewMethodForCalculatingZ = false;
	}else if(strcmp(argv[4],"SparseFI") == 0){
		fSetting->fMethod             = SparseFI;
	    // use new method for calculating Z based on step-wise least squares optimization (bool)
	    fSetting->useNewMethodForCalculatingZ = false;
	}else{
		cerr << "ERROR: UNKNOWN Fusion Method";
		exit(2);
	}

	//########################################################
	// local-non-local processing module
	//########################################################
	// Set Lagrangian multiplier
	fSetting->lambda              = atof(argv[5]);

	//*************************************
	// dictionary settings
	//*************************************
	// use simulated high resolution image X for dictionary learning (bool)
	fSetting->useSimulatedImXforDictLearn = (bool)atoi(argv[6]);
    // ImX simulation mode: 0: correlation based; 1: unconstrained LS based; 2: NNLS based
    fSetting->ImX_sim_mode        = atoi(argv[7]);
    // dictionary selection method
	fSetting->dictselect 		  = atoi(argv[8]);
	// Set number of nearest patches to be selected in order to form the reduced
	// dictionaries. Set npp=9999999 in order to use the full dictionary (all patches).
	fSetting->NDP                 = atoi(argv[9]);
	// Set size of (quadratic) patch, i.e. number of pixels in one direction.
	fSetting->patchsize           = atoi(argv[10]);
	// Set overlap, measured in pixels at the spatial scale of the LR MS image.
	fSetting->overlap             = atoi(argv[11]);

	//*************************************
	// coefficient estimation
	//*************************************
	// regularization parameters for new coefficient estimation (double)
	fSetting->lambdaX_ABC = (double)atof(argv[12]);
	fSetting->lambdaY_ABC = (double)atof(argv[13]);
	fSetting->lambdaZ_ABC = (double)atof(argv[14]);

	//*************************************
	// correlation-based Hyperspectral Grouping (CorHySpeG)
	//*************************************
	// maximum size of spectral group above groups will be double-checked and perhaps split into subgoups (int)
	fSetting->Nc_max = atoi(argv[15]);
	// minimum cross-correlation within spectral goups (double)
	fSetting->CC_min = (double)atof(argv[16]);
	// size of window around patch: Must have the same sign as patchsize in order to have both centers matched; Used for correlation calculations (int)
	fSetting->winSize = atoi(argv[17]);


	//########################################################
    //# Global processing module                             #
	//# (for full-image optimization)                        #
	//########################################################
	// >>>// flag used to decide whether or not the least square post-minimization (of the final image) is activated (bool)
	// >>>fSetting->LQ_post_opt_im  = (bool)atoi(argv[19]);
	// regularization parameter trading the relative weighting of the high resolution input image I_X (double)
	fSetting->lambdaX_im      = (double)atof(argv[18]);
	// regularization parameter trading the relative weighting of the low resolution input image I_Y (double)
	fSetting->lambdaY_im      = (double)atof(argv[19]);
	// maximum number of iterations in the GS step to solve the least squares problem on the final image level (int)
	sSetting->maxiter_CGLS_im = atoi(argv[20]);
	// error tolerance (double)
	sSetting->tol_r_CGLS_im   = (double)atof(argv[21]);

	// write all intermediate image fusion resulta (after every iteration) (1: create file and write resulting image in file; 0: to not write image in file (useful for analyses only)) (bool)
	fSetting->fullImOptOnSubspace = (bool)atoi(argv[22]);
        // subspace transformation type
        fSetting->subspace_transform_type = argv[23];
        // subspace dimension
        fSetting->subspace_dim = atoi(argv[24]);
	// calc. coeff. in full image opt. eq. via SNR calc. of ImX and ImY (bool)
	fSetting->SNR_normalization = (bool)atoi(argv[25]);


    //########################################################
    // initialization and interplay between local-non-local 
    // and global processing module
    //########################################################
    // use estimated SRFs instead of apriori given ones (bool)
	fSetting->use_estimated_SRFs = (bool)atoi(argv[26]);
	// type of initial high resolution image ImZ_init; flag (int)
	fSetting->ImZ_init_type=atoi(argv[27]);
	       // ImZ_init_type=0: lambdaZ_ABZ=0 in 1st iter (no initial image)
	       // ImZ_init_type=1: upsampled and bilinearly interpolated low resolution image ImY
	       // ImZ_init_type=2: reconstruction result of another algorithm (e.g. Bayesian Sparse or CNMF. Depends on dataset)
	
	// flag used to decide whether or not the least square post-minimization (of the final image) is activated (bool)
	fSetting->LQ_post_opt_im  = (bool)atoi(argv[28]);
    // jump to the full image optimization (of the initial image) without doing the patch-wise imge reconstruction.
	fSetting->doFullImOptWithoutPatRec=(bool)atoi(argv[29]);
	// number of coupled ImZ calculations iterations
	fSetting->iterMain = atoi(argv[30]);

	//########################################################
    //# Output settings                                      #
	//########################################################

	fSetting->evaluate           = (bool)atoi(argv[31]);
    //65 evaluate initial image (bool)
	fSetting->evaluate_ImZ_init = (bool)atoi(argv[32]);

	oSetting->writeImageFile      = (bool)atoi(argv[33]);

	// write all intermediate image fusion resulta (after every iteration) (1: create file and write resulting image in file; 0: to not write image in file (useful for analyses only)) (bool)
	oSetting->writeImageFileAfterEveryIter = (bool)atoi(argv[34]);
	// save output in double format (64bit) instead of uint16 (bool)
	dSetting->saveAsDouble = (bool)atoi(argv[35]);	


	//########################################################
    //**********************************
    // default parameters
    // advanced: change with caution!
    //**********************************
	//########################################################

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
	fSetting->ImZ_ref_avlbl       = true;
	oSetting->prec                = 10;

	//fSetting->pFirst             = 0; // atoi(argv[10]);
	//fSetting->pLast              = 999999999; // atoi(argv[11]);

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
	pSetting->workStealingTurns  = 1;
	
    dSetting->delete_tmp_patch_folders       = true;
	dSetting->imageConstructionOnly          = false;
	dSetting->contUnfinishedRec              = false;

//	Attention: set in getPaths.cpp!!!
//	if(dSetting->contUnfinishedRec){
//		paths->PathToIncompletePatchSetCSV = argv[35];
//	}

	dSetting->dir_tmp_patches_additional_num = 0; //atoi(argv[36]);

//  matrix (dictionary) normalization norm: 0: spectral norm, 1: Frobenious norm (bool)
    fSetting->matrixNorm        = true;
    fSetting->addMeanPixelwise  = false;
	//#######################################################
	//# J-P-FISTA settings for joint sparse reconstruction  #
    //#######################################################
	// Choose a solver for the sparse reconstruction (optimization) problem (currently JPFISTA only)
    sSetting->solver             = JPFISTA;
    // 39: maximum number of iterations to run the algorithm (int)
	sSetting->maxiter_out        = 200000;
	// 40: tolerance
	sSetting->tol                = 0.000000000001;
    //########################################################
	//# for Least squares post processing on the patch level #
	//########################################################
	// 41: flag used to decide whether or not the least square post-minimization is activated (bool)
	fSetting->LQ_post_opt  = false;
	// 42: regularization parameter trading the relative weighting of the high resolution input patch xHR (double)
	fSetting->lambdaX      = 1.0;
	// 43: regularization parameter trading the relative weighting of the low resolution input patch yLR (double)
	fSetting->lambdaY      = 1.0; 
	// 44: maximum number of iterations in the GS step to solve the least squares problem (int)
	sSetting->maxiter_CGLS = 1000;
	// 45: error tolerance (double)
	sSetting->tol_r_CGLS   = 0.000000000001;
	// 46: decides whether or not the coefficients get updates via least squares (bool)
	sSetting->fix_Alpha    = 1;
	// 47: decides whether or not the mean values of Z are set to the same mean values as Y (i.e. either delta_m remains the initial zero vector or it gets updated via least squares) (bool)
	sSetting->fix_delta_m  = true;

    //58 set lambdaZ_ABC to this number only in the first iteration. A low value can be helpful e.g. if the initial image ImZ_init is not very good/trustworthy (double)
	fSetting->lambdaZ_ABC_in_1st_iter = 1.0;

	//70: e.g. =0 for generic; =1 for SuperMUC-pr45ne and =2 for CG local (int)
	dSetting->platformID = 0;

	//75: set the coeff. of |R*Z-X| in full image opt. eq. to NChY/NChX [only relevant if SNR_normalization==1 ] (bool)
	fSetting->balance_ImX_term_coef = false;

    //77: use LR (low resolution) patch norm for normalization of corresponding LR and HR patch in coupled dictionaries. If set to 0 the HR nor is used by default. (bool)
	fSetting->use_LRnorm_for_dic_normalization = true; //
	//78: load and use a-priori calculated dictionaries. (bool)
	fSetting->load_DictHR_and_DictLR = false;
	

	// set negative values to zero (int)
	fSetting->set_neg_to_0 = 3;
    //			               = 0  => set negative values to zero only at the very end, before writing the final image
    //	                       = 1  => set negative values to zero only after patch reconstruction
    //	                       = 2  => set negative values to zero only after full image optimization
    //	                       = 3  => set negative values to zero both after patch reconstruction and after full image optimization


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
	if(fSetting->LQ_post_opt && pSetting->numProcPerPatch>1){
		if(my_rank==0){
			cerr << endl << "ERROR: The least-squares based post processing step (program argument LQ_post_opt) is set active. This step is only supported for ONE PROCESSOR PER PATCH! However, the program argument numProcPerPatch is set to a number greater than one!" << endl << endl;
		}
		MPI_Barrier(MPI_COMM_WORLD);
		exit(2);
	}
	if(pSetting->store_patches_tmp_on_drive && dSetting->contUnfinishedRec ){
		if(my_rank==0){
			cerr << endl << "ERROR: The two flags 'store_patches_tmp_on_drive' and 'contUnfinishedRec' are not compatible!" << endl << endl;
		}
		MPI_Barrier(MPI_COMM_WORLD);
		exit(2);
	}
//}




}

