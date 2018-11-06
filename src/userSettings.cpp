/*
 * userSettings.cpp
 *
 *  Created on: Jan 20, 2014
 *      Author: Grohnfeldt, Claas
 */

#include "userSettings.h"

// for nnls tests >>>>>
//#include <Eigen/Eigen>
//#include "nnls.h"
//using namespace Eigen;
//<<<<<<<<<<<<<<<<<<<<<


void getUserSettings(SpEODataIOSetting *dSetting, SpEOFusionSetting *fSetting, SpEOOutputSetting *oSetting, SpEOSolverSetting *sSetting, SpEOParallelSetting *pSetting, int argc, char **argv){

    
	/*=================================================================*
	 *                          User Settings                          *
	 *=================================================================*/

	//MPI_Status recv_status;
	int my_processes, my_rank;
	MPI_Comm_size(MPI_COMM_WORLD, &my_processes);
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    //**********************************
    // default parameters
    //**********************************
	// MPI / parallelization settings
	// total number of processes
	pSetting->numProcTot = my_processes;
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
	dSetting->jobName             = argv[1];
	dSetting->jobID               = argv[2];
	// Set Lagrangian multiplier
	fSetting->lambda              = atof(argv[4]);

	// Set number of nearest patches to be selected in order to form the reduced
	// dictionaries. Set npp=9999999 in order to use the full dictionary (all patches).
	fSetting->NDP                 = atoi(argv[5]);
	fSetting->two_step_estimation = (bool)atoi(argv[6]); // false
	fSetting->Nc                  = atoi(argv[7]); // 5;	// [More relevant for J-SparseFI-HM]
	fSetting->No                  = atoi(argv[8]); // 4;	// [More relevant for J-SparseFI-HM]
	// Set size of (quadratic) patch, i.e. number of pixels in one direction.
	fSetting->patchsize           = atoi(argv[9]); // 5;
	// Set overlap, measured in pixels at the spatial scale of the LR MS image.
	fSetting->overlap             = atoi(argv[10]); // 4;

	// Choose a sparse image fusion methods
	//fSetting->fMethod             = JSparseFI;
	//                            = SparseFI
	//                         or = JSparseFI
	//                         or = JSparseFIHM
	if(strcmp(argv[11],"GroupedJSparseFI") == 0){
		fSetting->fMethod             = GroupedJSparseFI;
	}else if(strcmp(argv[11],"JJSparseFIHM") == 0){
		fSetting->fMethod             = JSparseFIHM;
	}else if(strcmp(argv[11],"JSparseFI") == 0){
		fSetting->fMethod             = JSparseFI;
	}else if(strcmp(argv[11],"SparseFI") == 0){
		fSetting->fMethod             = SparseFI;
	}else{
//	if(~strcmp(argv[11],"SparseFI")){
//		fSetting->fMethod             = SparseFI;
//	}else if(~strcmp(argv[11],"JSparseFI")){
//		fSetting->fMethod             = JSparseFI;
//	}else if(~strcmp(argv[11],"JJSparseFIHM")){
//		fSetting->fMethod             = JSparseFIHM;
//	}else if(~strcmp(argv[11],"GroupedJSparseFI")){
//		fSetting->fMethod             = GroupedJSparseFI;
//	}else{
		cerr << "ERROR: UNKNOWN Fusion Method";
		exit(2);
	}

//	fSetting->pFirst             = atoi(argv[10]);
//	fSetting->pLast              = atoi(argv[11]);
	fSetting->evaluate           = (bool)atoi(argv[12]);
	//	fSetting->tol_SRF             = 0.15;	// [More relevant for J-SparseFI-HM]
	fSetting->tol_SRF            = atof(argv[13]);

	pSetting->store_patches_tmp_on_drive = (bool)atoi(argv[14]);
	pSetting->parWrNumProc       = atoi(argv[15]);
	// output settings
	oSetting->saveAlphas         = (bool)atoi(argv[16]);
	oSetting->pFirstAlpha        = atoi(argv[17]);
  	oSetting->pLastAlpha         = atoi(argv[18]);

	oSetting->saveDicts          = (bool)atoi(argv[19]);
	oSetting->pFirstDict         = atoi(argv[20]);
  	oSetting->pLastDict          = atoi(argv[21]);
  	dSetting->chBundleFirst      = atoi(argv[22]);
  	dSetting->chBundleLast       = atoi(argv[23]);

	dSetting->uLFirst            = atoi(argv[24]);
	dSetting->uLLast             = atoi(argv[25]);
	dSetting->vLFirst            = atoi(argv[26]);
	dSetting->vLLast             = atoi(argv[27]);

	pSetting->numProcPerPatch    = atoi(argv[28]);
	// patch parallelization strategy: work stealing (true) or fixed work scheduling (false)
	pSetting->workStealingTurns  = atoi(argv[29]);

	fSetting->dictselect 		  = atoi(argv[30]);

	oSetting->writeImageFile      = (bool)atoi(argv[31]);

	dSetting->delete_tmp_patch_folders       = (bool)atoi(argv[32]);

	dSetting->imageConstructionOnly          = (bool)atoi(argv[33]);

	dSetting->contUnfinishedRec              = (bool)atoi(argv[34]);

//	Attention: set in getPaths.cpp!!!
//	if(dSetting->contUnfinishedRec){
//		paths->PathToIncompletePatchSetCSV = argv[35];
//	}

	dSetting->dir_tmp_patches_additional_num = atoi(argv[36]);

//  TMP >>>>>>
//  37: matrix (dictionary) normalization norm: 0: spectral norm, 1: Frobenious norm (bool)
    fSetting->matrixNorm        = (bool)atoi(argv[37]);
    fSetting->addMeanPixelwise  = (bool)atoi(argv[38]);
    //  <<<<<< TMP
	//#######################################################
	//# J-P-FISTA settings for joint sparse reconstruction  #
    //#######################################################
	// Choose a solver for the sparse reconstruction (optimization) problem (currently JPFISTA only)
    sSetting->solver             = JPFISTA;
    //#######################  39: maximum number of iterations to run the algorithm (int)
	sSetting->maxiter_out        = atoi(argv[39]); // 400000;//200000;
	//#######################  40: tolerance
	sSetting->tol                = (double)atof(argv[40]); // 1e-14; //1e-12; // tolerance (double)
    //########################################################
	//# for Least squares post processing on the patch level #
	//########################################################
	//#######################  41: flag used to decide whether or not the least square post-minimization is activated (bool)
	fSetting->LQ_post_opt  = (bool)atoi(argv[41]);
	//#######################  42: regularization parameter trading the relative weighting of the high resolution input patch xHR (double)
	fSetting->lambdaX      = (double)atof(argv[42]);
	//#######################  43: regularization parameter trading the relative weighting of the low resolution input patch yLR (double)
	fSetting->lambdaY      = (double)atof(argv[43]);
	//#######################  44: maximum number of iterations in the GS step to solve the least squares problem (int)
	sSetting->maxiter_CGLS = atoi(argv[44]);
	//#######################  45: error tolerance (double)
	sSetting->tol_r_CGLS   = (double)atof(argv[45]);
	//#######################  46: decides whether or not the coefficients get updates via least squares (bool)
	sSetting->fix_Alpha    = atoi(argv[46]);
	//#######################  47: decides whether or not the mean values of Z are set to the same mean values as Y (i.e. either delta_m remains the initial zero vector or it gets updated via least squares) (bool)
	sSetting->fix_delta_m  = (bool)atoi(argv[47]);
	//########################################################
	//# for Least squares post processing on the image level #
	//########################################################
	//#######################  48: flag used to decide whether or not the least square post-minimization (of the final image) is activated (bool)
	fSetting->LQ_post_opt_im  = (bool)atoi(argv[48]);
	//#######################  49: regularization parameter trading the relative weighting of the high resolution input image I_X (double)
	fSetting->lambdaX_im      = (double)atof(argv[49]);
	//#######################  50: regularization parameter trading the relative weighting of the low resolution input image I_Y (double)
	fSetting->lambdaY_im      = (double)atof(argv[50]);
	//#######################  51: maximum number of iterations in the GS step to solve the least squares problem on the final image level (int)
	sSetting->maxiter_CGLS_im = atoi(argv[51]);
	//#######################  52: error tolerance (double)
	sSetting->tol_r_CGLS_im   = (double)atof(argv[52]);
	//########################################################
	// for coefficient estimation
	//########################################################
	//#######################  53: use new method for calculating Z based on step-wise least squares optimization [added in October 2015] (bool)
	fSetting->useNewMethodForCalculatingZ = (bool)atoi(argv[53]);
	//#######################  54: use simulated high resolution image X for dictionary learning [added in October 2015] (bool)
	fSetting->useSimulatedImXforDictLearn = (bool)atoi(argv[54]);
	//#######################  55~57: regularization parameters for new coefficient estimation  [added in October 2015] (double)
	fSetting->lambdaX_ABC = (double)atof(argv[55]); // fSetting->lambdaX_im;//(double)atof(argv[55]);
	fSetting->lambdaY_ABC = (double)atof(argv[56]); // fSetting->lambdaY_im;//(double)atof(argv[56]);
	fSetting->lambdaZ_ABC = (double)atof(argv[57]); // 1.0;                 //(double)atof(argv[57]);
    //####################### 58 set lambdaZ_ABC to this number only in the first iteration. A low value can be helpful e.g. if the initial image ImZ_init is not very good/trustworthy (double)
	fSetting->lambdaZ_ABC_in_1st_iter = (double)atof(argv[58]); // 1.0;//
	//####################### 59  type of initial high resolution image ImZ_init; flag (int)
	fSetting->ImZ_init_type=atoi(argv[59]);
	       // ImZ_init_type=0: lambdaZ_ABZ=0 in 1st iter (no initial image)
	       // ImZ_init_type=1: upsampled and bilinearly interpolated low resolution image ImY
	       // ImZ_init_type=2: reconstruction result of another algorithm (e.g. Bayesian Sparse or CNMF. Depends on dataset)
	//####################### 60 jump to the full iimage optimization (of the initial image) without doing the patch-wise imge reconstruction.
	fSetting->doFullImOptWithoutPatRec=(bool)atoi(argv[60]);
	//####################### 61 number of coupled ImZ calculations iterations
	fSetting->iterMain = atoi(argv[61]);
	//####################### 62 maximum size of spectral group above groups will be double-checked and perhaps split into subgoups (int)
	fSetting->Nc_max = atoi(argv[62]);
	//####################### 63 minimum cross-correlation within spectral goups (double)
	fSetting->CC_min = (double)atof(argv[63]);
	//####################### 64 size of window around patch: Must have the same sign as patchsize in order to have both centers matched; Used for correlation calculations (int)
	fSetting->winSize = atoi(argv[64]);

    //####################### 65 evaluate initial image (bool)
	fSetting->evaluate_ImZ_init = (bool)atoi(argv[65]);
	//####################### 66 set negative values to zero (int)
	fSetting->set_neg_to_0 = atoi(argv[66]);
//			               = 0  => set negative values to zero only at the very end, before writing the final image
//	                       = 1  => set negative values to zero only after patch reconstruction
//	                       = 2  => set negative values to zero only after full image optimization
//	                       = 3  => set negative values to zero both after patch reconstruction and after full image optimization
    //####################### 67 use estimated SRFs instead of apriori given ones (bool)
	fSetting->use_estimated_SRFs = (bool)atoi(argv[67]);

	//####################### 68: write all intermediate image fusion resulta (after every iteration) (1: create file and write resulting image in file; 0: to not write image in file (useful for analyses only)) (bool)
	oSetting->writeImageFileAfterEveryIter = (bool)atoi(argv[68]);
	//####################### 69: write all intermediate image fusion resulta (after every iteration) (1: create file and write resulting image in file; 0: to not write image in file (useful for analyses only)) (bool)
	fSetting->fullImOptOnSubspace = (bool)atoi(argv[69]);
	//####################### 70: e.g. =1 for SuperMUC and =2 for CG local CP (int)
	dSetting->platformID = atoi(argv[70]);
        //####################### 71: subspace transformation type
        fSetting->subspace_transform_type = argv[71];
        //####################### 72: subspace dimension
        fSetting->subspace_dim = atoi(argv[72]);
        //####################### 73: ImX simulation mode: 0: correlation based; 1: unconstrained LS based; 2: NNLS based
        fSetting->ImX_sim_mode = atoi(argv[73]);
	//####################### 74: calc. coeff. in full image opt. eq. via SNR calc. of ImX and ImY (bool)
	fSetting->SNR_normalization = (bool)atoi(argv[74]);
	//####################### 75: set the coeff. of |R*Z-X| in full image opt. eq. to NChY/NChX [only relevant if SNR_normalization==1 ] (bool)
	fSetting->balance_ImX_term_coef = (bool)atoi(argv[75]);
	//####################### 76: save output in double format (64bit) instead of uint16 (bool)
	dSetting->saveAsDouble = (bool)atoi(argv[76]);	
        //####################### 77: use LR (low resolution) patch norm for normalization of corresponding LR and HR patch in coupled dictionaries. If set to 0 the HR nor is used by default. (bool)
	fSetting->use_LRnorm_for_dic_normalization = (bool)atoi(argv[77]);
	//####################### 78: load and use a-priori calculated dictionaries. (bool)
	fSetting->load_DictHR_and_DictLR = (bool)atoi(argv[78]);
	

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

