/*#############################################################################
##                                                                           ##
##    File:        JSparseFI.cpp                                             ##
##                                                                           ##
##    Author:      Claas Grohnfeldt                                          ##
##                                                                           ##
##    Contributing co-authors:                                               ##
##                 - Steffen Peter                                           ##
##                 - Xiaoxiang Zhu                                           ##
##                                                                           ##
##    Contact:     Claas Grohnfeldt                                          ##
##                 E-mail: claas.grohnfeldt@gmail.com                        ##
##                                                                           ##
##    Version:     1.2.2 (11/2018)                                           ##
##                                                                           ##
##    Last modification: November 9, 2018                                    ##
##                                                                           ##
##    Copyright:   Claas Grohnfeldt                                          ##
##                                                                           ##
##    Description: Main file for the image fusion methods SparseFI,          ##
##                 J-SparseFI, and J-SparseFI-HM.                            ##
##                                                                           ##
##    Input:       ---                                                       ##
##                                                                           ##
##    Output:      ---                                                       ##
##                                                                           ##
#############################################################################*/

#include "JSparseFI.h"
#include "nnls.h"


using namespace std;
using namespace Eigen;

int main(int argc, char* argv[]) {

	//==================================================================//
	//                        MPI initialization                        //
	//==================================================================//
#ifdef _OPENMP
      int iprovided;
      MPI_Init_thread(&argc,&argv,MPI_THREAD_FUNNELED,&iprovided);
#else
      MPI_Init(&argc,&argv);
#endif
      //MPI_Status recv_status;
      int my_rank; int my_processes;
      MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
      MPI_Comm_size(MPI_COMM_WORLD, &my_processes);

  	//==================================================================//
  	//                         Get user settings                        //
  	//==================================================================//
      SpEODataIOSetting dSetting;
      SpEOFusionSetting fSetting;
      SpEOOutputSetting oSetting;
      SpEOSolverSetting sSetting;
      SpEOParallelSetting pSetting;

      getUserSettings(&dSetting, &fSetting, &oSetting, &sSetting, &pSetting, argc, argv);
      
    //==================================================================//
	//                         Get data set paths                       //
	//==================================================================//
	  SpEOPaths paths;
	  getPaths(&paths,&dSetting,&pSetting,argc, argv);

    //==================================================================//
  	//          Initialize report and broadcast output file name        //
  	//==================================================================//
  	  SpEOReport report;
  	  if(my_rank == 0){
  		  report.initialize(&paths, &dSetting, &fSetting, &oSetting, &sSetting, &pSetting, argc, argv);
  	  }

  	  // broadcast the output directory name to all processes
  	  int mylen = 0;
  	  if(my_rank==0){
  		  mylen = paths.dir_out.length();
  	  }
  	  MPI_Bcast(&mylen, 1, MPI_INT, 0,MPI_COMM_WORLD);
  	  char buf2[mylen+1];
  	  if(my_rank==0){
  		sprintf (buf2, "%s", paths.dir_out.c_str());
  	  }
  	  MPI_Bcast(&buf2, sizeof(buf2), MPI_CHAR, 0,MPI_COMM_WORLD);
  	  string mystringstr(buf2);
  	  paths.dir_out = mystringstr;


	//==================================================================//
	//                          Load Images                             //
	//==================================================================//
  	  SpEODataset *ImX, *ImY, *ImZ, *ImZ_init;
	  ImX        = new SpEODataset(HR, imFlag_X);
	  ImY        = new SpEODataset(LR, imFlag_Y);
	  ImZ        = new SpEODataset(HR, imFlag_Z);
	  ImZ_init   = new SpEODataset(HR, imFlag_Z_init);
	  ImX->dataRead(     paths.fname_ImX,      &report);
	  ImY->dataRead(     paths.fname_ImY,      &report);

	  switch(fSetting.ImZ_init_type){
		  case 0:{
			  fSetting.lambdaZ_ABC_in_1st_iter=0;
			  ImZ_init->dataRead(paths.fname_ImZ_ref, &report);
			  break;
		  }case 1:{
			  ImZ_init->dataRead(paths.fname_ImZ_init_ImY_US, &report);
			  break;
		  }case 2:{
			  ImZ_init->dataRead(paths.fname_ImZ_init_rec, &report);
			  break;
		  }case 3:{
			  ImZ_init->dataRead(paths.fname_ImZ_ref, &report);
			  break;
		  }
	  }

	//=======================================================================//
	//  possibly correct ImY and ImX to selected spectral and spatial subset //
	//=======================================================================//
	  SpEOGlobalParams glPrms;
	  glPrms.NChX_orig = ImX->get_NCh();
	  glPrms.NChY_orig = ImY->get_NCh();
	  cutRelevantInput(ImX, ImY, &dSetting);
	  if(fSetting.subspace_dim > ImY->get_NCh()){
		fSetting.subspace_dim = ImY->get_NCh();
	  }

	//==================================================================//
	//                    Calculate global parameters                   //
	//                      and check input settings                    //
	//==================================================================//
	  //SpEOGlobalParams glPrms;
	  calcGlobalParams(&glPrms, &paths, &dSetting, &fSetting, ImY, ImX);
	  setMetaInfo(ImZ, ImY, ImX, &dSetting, &glPrms);

	  if(!pSetting.store_patches_tmp_on_drive){
		  if( (fSetting.fMethod == GroupedJSparseFI) &&  (    (dSetting.uLLast - dSetting.uLFirst +1 < ImY->get_sizeU())
				  	  	  	  	  	  	  	  	  	  	   || (dSetting.vLLast - dSetting.vLFirst +1 < ImY->get_sizeV()) ) ){
			  if(my_rank==0){
				cerr << endl << "ERROR: The algorithm 'GroupedJSparseFI' is not compatible to calculating only a spatial sub-image, because, for the dictionary generation, the intermediate pan-images (which are reconstruction results) need to be of the same size as the original Pan image!" << endl << endl;
			  }
			  MPI_Barrier(MPI_COMM_WORLD);
			  exit(2);
		  }
	  }

	//=======================================================================================//
	//  Read subspace transformation matrix file (e.g. endmember matrix in case of unmixing) //
	//=======================================================================================//
	  // declare, initialize and possibly load endmember matrix
	  
	  SpEOMatrixD SubspaceTransformMat(ImY->get_NCh(),fSetting.subspace_dim);
	  if(fSetting.subspace_transform_type == "none"){
	     SubspaceTransformMat =  SpEOMatrixD::Identity(ImY->get_NCh(),ImY->get_NCh());
	  }else if(fSetting.subspace_transform_type == "VCA"){
	     if(fSetting.fullImOptOnSubspace){
	             std::string fname_SubspaceTransformMat = paths.fname_SubspaceTransformMat;
	             if (my_rank==0) cout << "read subspace transformation matrix from file: " << endl << fname_SubspaceTransformMat << endl;
	     	  int stat_CSV_read = read_CSV(&SubspaceTransformMat, fname_SubspaceTransformMat.c_str(),',',0);
	             if(stat_CSV_read==-1){
                            if(my_rank==0){
                   	        cout << "["<< my_rank << "] ERROR: The file '" << fname_SubspaceTransformMat.c_str() << "' does not exist" << endl;
	                    }
                            MPI_Barrier(MPI_COMM_WORLD);
                            if(my_rank==0){
                   	        cerr << "["<< my_rank << "] ERROR: The file '" << fname_SubspaceTransformMat.c_str() << "' does not exist" << endl;
	                    }
                   	 exit(2);
	             }	  

		     JacobiSVD<SpEOMatrixD> svd(SubspaceTransformMat, ComputeThinU);
		     SpEOMatrixD SVD_mat_U = svd.matrixU();
		     fSetting.subspace_dim = SubspaceTransformMat.cols();
		     SubspaceTransformMat = SVD_mat_U.leftCols(fSetting.subspace_dim);
	     }
	  }else if(fSetting.subspace_transform_type == "SVD"){
	     SpEOMatrixD ImY_2D    = SpEOMatrixD::Zero(ImY->get_NCh(),   ImY->get_sizeU()   *ImY->get_sizeV());
	     transform_SpEODataset_to_2D(ImY    ,ImY_2D    );
	     ImY_2D /= ImY_2D.mean();
	     JacobiSVD<SpEOMatrixD> svd(ImY_2D, ComputeThinU);
	     SpEOMatrixD SVD_mat_U = svd.matrixU();
	     SubspaceTransformMat = SVD_mat_U.leftCols(fSetting.subspace_dim);	
	     ImY_2D    = SpEOMatrixD::Zero(0,0);
	  }else if(fSetting.subspace_transform_type == "PCA"){
	     SpEOMatrixD ImY_2D    = SpEOMatrixD::Zero(ImY->get_NCh(),   ImY->get_sizeU()   *ImY->get_sizeV());
             transform_SpEODataset_to_2D(ImY    ,ImY_2D    );
             ImY_2D /= ImY_2D.mean();
             JacobiSVD<SpEOMatrixD> svd(ImY_2D, ComputeThinU);
	     SpEOMatrixD SVD_mat_U = svd.matrixU();
	     SpEOVectorD SVD_vec_S = svd.singularValues();
	     // reduce U and S to subspace
	     SpEOMatrixD SVD_mat_U_sub = SVD_mat_U.leftCols(fSetting.subspace_dim);
	     SpEOVectorD SVD_vec_S_sub = SVD_vec_S.head(fSetting.subspace_dim);
	     SpEOMatrixD SVD_mat_S_sub = SVD_vec_S_sub.asDiagonal();
	     
             SubspaceTransformMat = SVD_mat_U_sub*SVD_mat_S_sub;
             //[U,S,~] = svd(XH_2D');
	     ImY_2D    = SpEOMatrixD::Zero(0,0);
	  }else{
	     cerr << endl << "unknown subspace Transformation Type: " << fSetting.subspace_transform_type << endl;
	     exit(2); 
	  }

	  glPrms.NChY_subspace = SubspaceTransformMat.rows();
	  // declare, initialize and possibly load the initial value of the abundance matrix
	  fSetting.use_init_value_Eq1Unmixing = false;
	  SpEOMatrixD* AbundanceMat = new SpEOMatrixD[glPrms.NChY_subspace];	  
          for(int iChY_sub=0; iChY_sub<glPrms.NChY_subspace; iChY_sub++){
               AbundanceMat[iChY_sub] = SpEOMatrixD::Zero(ImX->get_sizeU(),ImX->get_sizeV());
      	  }
          MPI_Barrier(MPI_COMM_WORLD);

	//==================================================================//
	//  Low-pass filter and down-sample ImX to generate ImX_LR          //
	//==================================================================//
	  SpEODataset *ImX_LR;
	  ImX_LR = new SpEODataset(LR, imFlag_X_LR);
	  ImX_LR->copyMetaInfoFromDatasets(ImY, ImX, ImY,SpEOFloat);
	  lowPassFilter_and_downSample(ImX,ImX_LR,SpEOFloat,SpEOFloat,glPrms);

	//==================================================================//
	//  Calculate Gauss filter for low-pass filtering and down-sampling //
	//==================================================================//
	  // calculate Gauss filter for low-pass filtering
	  SpEOMatrixD gauss_filter;
	  int filter_size = 2*glPrms.fDS - glPrms.fDS%2;
	  create_Gaussian_filter(gauss_filter, filter_size);
	  SpEOMatrixD filter_coeff;
	  calc_filter_boundary_coeff(filter_coeff, gauss_filter, glPrms.fDS);

	//==================================================================//
	//              Read spectral response functions file               //
    //   or estimate SRFs from ImY and ImX_LR together with a band-wise //
    //   shift in ImX                                                   //
	//==================================================================//
	  SpEOMatrixD SRF(ImX->get_NCh(), ImY->get_NCh());
	  
	  if(!fSetting.use_estimated_SRFs){
		  bool normalize_SRFs_to_fulfill_SumToOne_constraint = false; // potentially beneficial when processing more realistic data (not synthesized ones)
		  read_SRF(&glPrms, &SRF, paths.fname_SRF,',', normalize_SRFs_to_fulfill_SumToOne_constraint);
		  MPI_Barrier(MPI_COMM_WORLD);
          }else{
                  if(my_rank == 0){
                     cout << "#############################################################" << endl 
                          << "##                                                         ##" << endl
                          << "##   estimate spectral response functions and ImX shift    ##" << endl
                          << "##                                                         ##" << endl
                          << "#############################################################" << endl;
		  }
		  double *ImX_shift = new double[ImX->get_NCh()];
		  SpEOMatrixD ImY_2D    = SpEOMatrixD::Zero(ImY->get_NCh(),   ImY->get_sizeU()   *ImY->get_sizeV());
		  SpEOMatrixD ImX_LR_2D = SpEOMatrixD::Zero(ImX_LR->get_NCh(),ImX_LR->get_sizeU()*ImX_LR->get_sizeV());
	
		  transform_SpEODataset_to_2D(ImY    ,ImY_2D    );
		  transform_SpEODataset_to_2D(ImX_LR ,ImX_LR_2D );
	          if(my_rank == 0) cout << "estimate SFRs and ImX shift (band-wise) simultaneously via non-negative Least Squares.." << endl;
		  bool est_successful = estimate_SRFs_and_ImX_shift(SRF, ImX_shift, ImY_2D, ImX_LR_2D, my_rank);
	 	  if(my_rank == 0) cout << ".. done! The estimation completed with status: successful=" << est_successful << endl << endl;
		  if(my_rank == 0) cout << "shift both ImX and ImX_LR band-wise..." << endl;
		  ImX->shift_image_bandwise(ImX_shift);
		  ImX_LR->shift_image_bandwise(ImX_shift);
	 	  if(my_rank==0) cout << ".. done!" << endl << endl;
	
	      // clean up
	      delete[] ImX_shift;
	      ImY_2D    = SpEOMatrixD::Zero(0,0);
	      ImX_LR_2D = SpEOMatrixD::Zero(0,0);
	}

	//==================================================================//
	//                  Calculate decision matrix C                     //
	//     according to the spectral grouping concept based on SRFs     //
	//==================================================================//
	  if (paths.fname_SRF != paths.fname_SRF_for_Spectral_Grouping){
		  if(my_rank==0){
			  cout << "The path names stored in fname_SRF and fname_SRF_for_Spectral_Grouping are not identical. " <<
					  "Therefore, there will be an additional SRF file loaded, supposedly containing a " <<
					  "modified version of the original SRFs, which decreases the uncertainty in the " <<
					  "calculation of the decision matrix in the spectral grouping step based on the SRFs." << endl;
		  }
		  SpEOMatrixD SRF_SG;
		  read_SRF(&glPrms, &SRF_SG, paths.fname_SRF_for_Spectral_Grouping,',', false);
		  MPI_Barrier(MPI_COMM_WORLD);
		  calcDecisionMat(&SRF_SG, &dSetting, &fSetting, &pSetting, &glPrms);
	  }else{
		  if(my_rank==0){
			  cout << "The path names stored in fname_SRF and fname_SRF_for_Spectral_Grouping are identical. " <<
					  "Therefore, the original SRFs are used for the calculation of the decision matrix in the spectral grouping step." << endl;
		  }
		  MPI_Barrier(MPI_COMM_WORLD);
		  calcDecisionMat(&SRF, &dSetting, &fSetting, &pSetting, &glPrms);
	  }
	  MPI_Barrier(MPI_COMM_WORLD);

	 //===================================================================//
	 //              Create MPI groups and sub-communicators              //
	 //===================================================================//
	  if(my_rank==0){
		  cout << "Generate MPI communicators .." << endl;
	  }

	  MPI_Group mpi_group_orig;
	  MPI_Comm_group(MPI_COMM_WORLD, &mpi_group_orig);

	  //################
	  //## comm_busy  ##
	  //################
	  MPI_Group group_busy;
#ifdef _OPENMP
	  int numProcBusy = glPrms.numPatchGroups;
#else
	  int numProcBusy = min(glPrms.numPatchGroups*pSetting.numProcPerPatch,glPrms.NP_sub*pSetting.numProcPerPatch);
#endif
	  int ranges_busy[1][3] = {{0,numProcBusy-1,1}};
	  int rangesIdl_busy[1][3] = {{numProcBusy,my_processes-1,1}};
	  if(my_rank < numProcBusy){
		  MPI_Group_range_incl(mpi_group_orig, 1, ranges_busy, &group_busy);
	  }else{
		  MPI_Group_range_incl(mpi_group_orig, 1, rangesIdl_busy, &group_busy);
	  }
	  MPI_Comm comm_busy;
	  MPI_Comm_create(MPI_COMM_WORLD, group_busy, &comm_busy);
	  
	  //################
	  //## comm_write ##
	  //################
	  MPI_Group group_write;
	  int numProcWrite = min(numProcBusy,pSetting.parWrNumProc);
	  int ranges_write[1][3] = {{0,numProcWrite-1,1}};
	  int rangesIdl_write[1][3] = {{numProcWrite,my_processes-1,1}};
	  if(my_rank < numProcWrite){
		  MPI_Group_range_incl(mpi_group_orig, 1, ranges_write, &group_write);
	  }else{
		  MPI_Group_range_incl(mpi_group_orig, 1, rangesIdl_write, &group_write);
	  }
	  MPI_Comm comm_write;
	  MPI_Comm_create(MPI_COMM_WORLD, group_write, &comm_write);

	//==================================================================//
	//                      Call image fusion method                    //
	//==================================================================//
	  if(my_rank < numProcBusy){
		  int final_evaluation = false;
		  if(my_rank==0){
			 //==================================================================//
			 //                   Assess initial image ImZ_init                  //
			 //==================================================================//
			  if(fSetting.evaluate){
			  	string dir_eval = paths.dir_out + "/" + "eval";
			  	mkdir(dir_eval.c_str(), 0777);
			  	chmod(dir_eval.c_str(), S_IRWXU | S_IRWXG | S_IRWXO);
			  }
			  if(fSetting.evaluate && fSetting.evaluate_ImZ_init){
				  report.file.open(report.fileName.c_str(),	fstream::in | fstream::out | fstream::app);
				  report.file << "evaluation of ImZ_init:" << "\n";
				  report.file.close();
				  cout << "evaluation of ImZ_init:" << "\n";
				  // load high resolution reference image
				  SpEODataset *ImZ_ref;
				  ImZ_ref =	new SpEODataset(HR, imFlag_Z_ref);
				  ImZ_ref->dataRead(paths.fname_ImZ_ref, &report);
				  // possibly correct Z_ref and ImY to selected spectral and spatial subset
				  cutRelevantInput(HR, imFlag_Z_ref, ImZ_ref, &dSetting, &glPrms,false);
				  cutRelevantInput(LR, imFlag_Y, ImY, &dSetting, &glPrms,false);
				  SpEOAssessmentMetrics assMetrics_HR;
				  // do the assessment
				  PanSharp_Assessment_Eval(ImZ_ref, ImZ_init, ImY, &assMetrics_HR, &fSetting, &glPrms);
				  // write assessment results in report
				  report.addEvaluation(&assMetrics_HR, &dSetting, &fSetting, &oSetting, &glPrms, HR_ASSESSMENT);
				  // save assessment results
				  int iterMain=0;
				  bool init_image_eval=true;
				  save_evalResults(&paths, &glPrms, &fSetting, &assMetrics_HR, iterMain, fSetting.iterMain, true,fSetting.doFullImOptWithoutPatRec, init_image_eval, final_evaluation);
				  delete ImZ_ref;
				  delete[] assMetrics_HR.RMSE_sep;
				  delete[] assMetrics_HR.PSNR_sep;
				  delete[] assMetrics_HR.CC_sep;
				  delete[] assMetrics_HR.ERGAS_sep;
				  delete[] assMetrics_HR.UIQI_sep;
				  delete[] assMetrics_HR.DD_sep;
				  delete[] assMetrics_HR.AG_orig_sep;
				  delete[] assMetrics_HR.AG_rec_sep;
				  for(int iChZ=0; iChZ<glPrms.NChZ; iChZ++){
					  delete[] assMetrics_HR.DLambda_mat[iChZ];
				  }
				  delete[] assMetrics_HR.DLambda_mat;
			  }
		  }

		  if(fSetting.doFullImOptWithoutPatRec){
			  int iChZ;
			  for(iChZ=0; iChZ<ImZ->get_NCh(); iChZ++){
				  ImZ->get_rasterBands()[iChZ]->bandDataMatD = ImZ_init->get_rasterBands()[iChZ]->bandDataMatD;
			  }
			  glPrms.timeMainLoop = 0;
			  glPrms.timeDictSelect = 0;
			  glPrms.timeDictSelect_avg = 0;
		  }

		  // for stopping criteria:
		  int stoppingCriteriaFulfilled = 0;
		  SpEOVectorD SAM_patRec_in_all_iters  = SpEOVectorD::Zero(fSetting.iterMain);
		  SpEOVectorD CC_patRec_in_all_iters   = SpEOVectorD::Zero(fSetting.iterMain);
		  SpEOVectorD UIQI_patRec_in_all_iters = SpEOVectorD::Zero(fSetting.iterMain);
		  SpEOVectorD RMSE_patRec_in_all_iters = SpEOVectorD::Zero(fSetting.iterMain);
		  SpEOVectorD PSNR_patRec_in_all_iters = SpEOVectorD::Zero(fSetting.iterMain);
		  SpEOVectorD SAM_fullImRec_in_all_iters  = SpEOVectorD::Zero(fSetting.iterMain);
		  SpEOVectorD CC_fullImRec_in_all_iters   = SpEOVectorD::Zero(fSetting.iterMain);
		  SpEOVectorD UIQI_fullImRec_in_all_iters = SpEOVectorD::Zero(fSetting.iterMain);
		  SpEOVectorD RMSE_fullImRec_in_all_iters = SpEOVectorD::Zero(fSetting.iterMain);
		  SpEOVectorD PSNR_fullImRec_in_all_iters = SpEOVectorD::Zero(fSetting.iterMain);

		  int iterMain;
		  for(iterMain=0; iterMain<fSetting.iterMain && stoppingCriteriaFulfilled==0; iterMain++){
			  if(my_rank==0){
				  cout << endl
					   << "##########################################" << endl
					   << "#                                        #" << endl
					   << "#            iteration "<<iterMain<<"                 #" << endl
					   << "#                                        #" << endl
					   << "##########################################" << endl << endl;
			  }
			  if (!fSetting.doFullImOptWithoutPatRec){

				  if(fSetting.fMethod == GroupedJSparseFI){
					  //*************************************************
					  // currently only experimental version for the WorldView-2 pan-sharpening case and manual spectral grouping
					  // -----------------------------------------------
					  // (1) calculate main groups using the Joint Sparsity Model (JSM) and the
					  //     original Pan image for the generation of the coupled dictionaries
					  //     (in the WV-2 case there are the two main channel groups 1~4 and 5)
					  // (2) calculate outer channel groups using reconstruction results from
					  //     neighboring channels (as Pan image) for dictionary generatation
					  //************************************************
					  // check compatibility of input to the J-SparseFI algorithm (currently called 'GroupedJSparseFI')
					  checkInputForJSpFI(glPrms, fSetting, pSetting);
					  // save original parameters that need to be adapted to the individual groups
					  // ---------------->
					  int NChY_orig          = glPrms.NChY;
					  int NChZ_orig          = glPrms.NChZ;
					  int Nc_vec_0_orig      = glPrms.Nc_vec[0];
					  int Nc_orig            = fSetting.Nc;
					  int chBundleFirst_orig = dSetting.chBundleFirst;
					  int chBundleLast_orig  = dSetting.chBundleLast;
					  SpEOMatrixD SRF_orig   = SRF;
					  // <----------------

					  SpEODataset *ImX_tmp, *ImX_LR_tmp, *ImY_tmp, *ImZ_tmp;
					  int firstBandY, lastBandY, panBand;
					  ImX_tmp    = new SpEODataset(HR, imFlag_X);
					  ImX_LR_tmp = new SpEODataset(LR, imFlag_X_LR);
					  ImY_tmp    = new SpEODataset(LR, imFlag_Y);
					  ImZ_tmp    = new SpEODataset(HR, imFlag_Z);

		//			  // calculate Gauss filter for low-pass filtering the 'pan' images
		//			  SpEOMatrixD gauss_filter;
		//			  int filter_size = 2*glPrms.fDS - glPrms.fDS%2;
		//			  create_Gaussian_filter(gauss_filter, filter_size);
		//			  SpEOMatrixD filter_coeff;
		//			  calc_filter_boundary_coeff(filter_coeff, gauss_filter, glPrms.fDS);

					  //## main group, channels 1~4 #################

					  firstBandY = 1;
					  lastBandY  = 4;
					  panBand    = -1; // use original panchromatic image for the generation of the coupled dictionaries
					  prepGroupCalcForJSpFI( firstBandY, lastBandY, panBand, glPrms, fSetting, dSetting, SRF_orig, SRF, filter_coeff, gauss_filter,
											 ImX,     ImX_LR,     ImY,     ImZ,
											 ImX_tmp, ImX_LR_tmp, ImY_tmp, ImZ_tmp);
					  JSparseFI_alg(SubspaceTransformMat,AbundanceMat,iterMain, fSetting.iterMain, &paths, &dSetting, &fSetting, &oSetting, &sSetting, &pSetting, &glPrms, ImX_tmp, ImX_LR_tmp,                      ImY_tmp, ImZ_tmp, &SRF, comm_busy, group_busy, ImZ_init, report);
					  ImZ->fill_band_data(ImZ_tmp,firstBandY,lastBandY);
					  // restore original parameters
					  glPrms.NChY=NChY_orig; glPrms.NChZ=NChZ_orig; glPrms.Nc_vec[0]=Nc_vec_0_orig;fSetting.Nc=Nc_orig; dSetting.chBundleFirst=chBundleFirst_orig; dSetting.chBundleLast=chBundleLast_orig; SRF=SRF_orig;

					  //## main group, channel 5 #################
					  firstBandY = 5;
					  lastBandY  = 5;
					  panBand    = -1; // use original panchromatic image for the generation of the coupled dictionaries
					  prepGroupCalcForJSpFI( firstBandY, lastBandY, panBand, glPrms, fSetting, dSetting, SRF_orig, SRF, filter_coeff, gauss_filter,
											 ImX,     ImX_LR,     ImY,     ImZ,
											 ImX_tmp, ImX_LR_tmp, ImY_tmp, ImZ_tmp);

					  JSparseFI_alg(SubspaceTransformMat,AbundanceMat,iterMain, fSetting.iterMain, &paths, &dSetting, &fSetting, &oSetting, &sSetting, &pSetting, &glPrms, ImX_tmp, ImX_LR_tmp, ImY_tmp, ImZ_tmp, &SRF, comm_busy, group_busy, ImZ_init, report);
					  ImZ->fill_band_data(ImZ_tmp,firstBandY,lastBandY);
					  // restore original parameters
					  glPrms.NChY=NChY_orig; glPrms.NChZ=NChZ_orig; glPrms.Nc_vec[0]=Nc_vec_0_orig;fSetting.Nc=Nc_orig; dSetting.chBundleFirst=chBundleFirst_orig; dSetting.chBundleLast=chBundleLast_orig; SRF=SRF_orig;

					  //## side group, channel 0 #################
					  firstBandY = 0;
					  lastBandY  = 0;
					  panBand    = 1; // use original panchromatic image for the generation of the coupled dictionaries
					  prepGroupCalcForJSpFI( firstBandY, lastBandY, panBand, glPrms, fSetting, dSetting, SRF_orig, SRF, filter_coeff, gauss_filter,
											 ImX,     ImX_LR,     ImY,     ImZ,
											 ImX_tmp, ImX_LR_tmp, ImY_tmp, ImZ_tmp);
					  JSparseFI_alg(SubspaceTransformMat,AbundanceMat,iterMain, fSetting.iterMain, &paths, &dSetting, &fSetting, &oSetting, &sSetting, &pSetting, &glPrms, ImX_tmp, ImX_LR_tmp, ImY_tmp, ImZ_tmp, &SRF, comm_busy, group_busy, ImZ_init, report);
					  ImZ->fill_band_data(ImZ_tmp,firstBandY,lastBandY);
					  // restore original parameters
					  glPrms.NChY=NChY_orig; glPrms.NChZ=NChZ_orig; glPrms.Nc_vec[0]=Nc_vec_0_orig;fSetting.Nc=Nc_orig; dSetting.chBundleFirst=chBundleFirst_orig; dSetting.chBundleLast=chBundleLast_orig; SRF=SRF_orig;

					  //## side group, channels 6~7 #################
					  firstBandY = 6;
					  lastBandY  = 7;
					  panBand    = 5; // use original panchromatic image for the generation of the coupled dictionaries
					  prepGroupCalcForJSpFI( firstBandY, lastBandY, panBand, glPrms, fSetting, dSetting, SRF_orig, SRF, filter_coeff, gauss_filter,
											 ImX,     ImX_LR,     ImY,     ImZ,
											 ImX_tmp, ImX_LR_tmp, ImY_tmp, ImZ_tmp);
					  JSparseFI_alg(SubspaceTransformMat,AbundanceMat,iterMain, fSetting.iterMain, &paths, &dSetting, &fSetting, &oSetting, &sSetting, &pSetting, &glPrms, ImX_tmp, ImX_LR_tmp, ImY_tmp, ImZ_tmp, &SRF, comm_busy, group_busy, ImZ_init, report);
					  ImZ->fill_band_data(ImZ_tmp,firstBandY,lastBandY);
					  // restore original parameters
					  glPrms.NChY=NChY_orig; glPrms.NChZ=NChZ_orig; glPrms.Nc_vec[0]=Nc_vec_0_orig;fSetting.Nc=Nc_orig; dSetting.chBundleFirst=chBundleFirst_orig; dSetting.chBundleLast=chBundleLast_orig; SRF=SRF_orig;

					  delete ImX_tmp;
					  delete ImX_LR_tmp;
					  delete ImY_tmp;
					  delete ImZ_tmp;
				  }else{
					  JSparseFI_alg(SubspaceTransformMat,AbundanceMat,iterMain, fSetting.iterMain, &paths, &dSetting, &fSetting, &oSetting, &sSetting, &pSetting, &glPrms,
									ImX, ImX_LR, ImY, ImZ, &SRF, comm_busy, group_busy, ImZ_init, report);
				  }
							   if(my_rank==0){
								 //==================================================================//
								 //                   Assess reconstruction results                  //
								 //==================================================================//
								  if(fSetting.evaluate){
									  if(iterMain==fSetting.iterMain-1 && !fSetting.LQ_post_opt_im
											  && !(    fSetting.set_neg_to_0==0
												   || (fSetting.set_neg_to_0==1  &&  fSetting.doFullImOptWithoutPatRec)
												   || (fSetting.set_neg_to_0==2  &&  !fSetting.LQ_post_opt_im)
												  )
										){
										  final_evaluation = true;
									  }
									  // load high resolution reconstructed image
									  if(oSetting.writeImageFile && pSetting.store_patches_tmp_on_drive){
										  ImZ->dataRead(paths.fname_ImZ, &report);
									  }
									  // load high resolution reference image
									  SpEODataset *ImZ_ref;
									  ImZ_ref =	new SpEODataset(HR, imFlag_Z_ref);
									  ImZ_ref->dataRead(paths.fname_ImZ_ref, &report);
									  // possibly correct Z_ref and ImY to selected spectral and spatial subset
									  cutRelevantInput(HR, imFlag_Z_ref, ImZ_ref, &dSetting, &glPrms,false);
									  cutRelevantInput(LR, imFlag_Y, ImY, &dSetting, &glPrms,false);
									  SpEOAssessmentMetrics assMetrics_HR;
									  // do the assessment
									  PanSharp_Assessment_Eval(ImZ_ref, ImZ, ImY, &assMetrics_HR, &fSetting, &glPrms);
									  // write assessment results in report
									  report.addEvaluation(&assMetrics_HR, &dSetting, &fSetting, &oSetting, &glPrms, HR_ASSESSMENT);
									  // save assessment results
									  bool init_image_eval=false;
									  save_evalResults(&paths, &glPrms, &fSetting, &assMetrics_HR, iterMain, fSetting.iterMain, true,fSetting.doFullImOptWithoutPatRec, init_image_eval, final_evaluation);

									  SAM_patRec_in_all_iters(iterMain)  = assMetrics_HR.SAM;
									  CC_patRec_in_all_iters(iterMain)   = assMetrics_HR.CC_mean;
									  UIQI_patRec_in_all_iters(iterMain) = assMetrics_HR.UIQI_mean;
									  RMSE_patRec_in_all_iters(iterMain) = assMetrics_HR.RMSE_mean;
									  PSNR_patRec_in_all_iters(iterMain) = assMetrics_HR.PSNR_mean;

									  delete ImZ_ref;
									  delete[] assMetrics_HR.RMSE_sep;
									  delete[] assMetrics_HR.PSNR_sep;
									  delete[] assMetrics_HR.CC_sep;
									  delete[] assMetrics_HR.ERGAS_sep;
									  delete[] assMetrics_HR.UIQI_sep;
									  delete[] assMetrics_HR.DD_sep;
									  delete[] assMetrics_HR.AG_orig_sep;
									  delete[] assMetrics_HR.AG_rec_sep;
									  for(int iChZ=0; iChZ<glPrms.NChZ; iChZ++){
										  delete[] assMetrics_HR.DLambda_mat[iChZ];
									  }
									  delete[] assMetrics_HR.DLambda_mat;
								  }
							  }
			  }
			  
              //========================================================#
			  //#   Equation (1): improve reconstructed high resolution #
			  //#   image spatially and spectrally using SRFs           #
			  //#   and PSFs ( new step: added in Sep/Oct 2015 )        #
			  //========================================================>
			  if (fSetting.LQ_post_opt_im){
				  full_im_optimization_LS(report, SubspaceTransformMat, AbundanceMat, &glPrms, &fSetting, &sSetting, &dSetting, ImX, ImY, ImZ, &SRF, comm_busy, group_busy);
                  for (int iChZ=0; iChZ<glPrms.NChZ; iChZ++){
                      int datasize = ImZ->rasterBands[iChZ]->bandDataMatD.rows()*ImZ->rasterBands[iChZ]->bandDataMatD.cols();
                      MPI_Bcast((ImZ->rasterBands[iChZ]->bandDataMatD.data()), datasize, MPI_DOUBLE, 0, comm_busy);//MPI_COMM_WORLD);
                  }
                  if(my_rank==0){
                      //==================================================================//
					  //                   Assess reconstruction results                  //
					  //==================================================================//
					   if(fSetting.evaluate){
					       if(iterMain==fSetting.iterMain-1
					     		  && !(    fSetting.set_neg_to_0==0
					     			   || (fSetting.set_neg_to_0==1  &&  fSetting.doFullImOptWithoutPatRec)
					     			   || (fSetting.set_neg_to_0==2  &&  !fSetting.LQ_post_opt_im)
					     			  )
					     		  ){
					     	  final_evaluation = true;
					       }
					       // load high resolution reconstructed image
					       if(oSetting.writeImageFile && pSetting.store_patches_tmp_on_drive){
					     	  ImZ->dataRead(paths.fname_ImZ, &report);
					       }
					       // load high resolution reference image
					       SpEODataset *ImZ_ref;
					       ImZ_ref =	new SpEODataset(HR, imFlag_Z_ref);
					       ImZ_ref->dataRead(paths.fname_ImZ_ref, &report);
					       // possibly correct Z_ref and ImY to selected spectral and spatial subset
					       cutRelevantInput(HR, imFlag_Z_ref, ImZ_ref, &dSetting, &glPrms,false);
					       cutRelevantInput(LR, imFlag_Y, ImY, &dSetting, &glPrms,false);
					       SpEOAssessmentMetrics assMetrics_HR;
					       // do the assessment
					       PanSharp_Assessment_Eval(ImZ_ref, ImZ, ImY, &assMetrics_HR, &fSetting, &glPrms);
					       // write assessment results in report
					       report.addEvaluation(&assMetrics_HR, &dSetting, &fSetting, &oSetting, &glPrms, HR_ASSESSMENT);
					       // save assessment results
					       bool init_image_eval=true;
					       if(fSetting.evaluate_ImZ_init || !fSetting.doFullImOptWithoutPatRec){
					     	  init_image_eval=false;
					       }
					       save_evalResults(&paths, &glPrms, &fSetting, &assMetrics_HR, iterMain, fSetting.iterMain, false,fSetting.doFullImOptWithoutPatRec, init_image_eval, final_evaluation);

					       SAM_fullImRec_in_all_iters(iterMain)  = assMetrics_HR.SAM;
					       CC_fullImRec_in_all_iters(iterMain)   = assMetrics_HR.CC_mean;
					       UIQI_fullImRec_in_all_iters(iterMain) = assMetrics_HR.UIQI_mean;
					       RMSE_fullImRec_in_all_iters(iterMain) = assMetrics_HR.RMSE_mean;
					       PSNR_fullImRec_in_all_iters(iterMain) = assMetrics_HR.PSNR_mean;

					       delete ImZ_ref;
					       delete[] assMetrics_HR.RMSE_sep;
					       delete[] assMetrics_HR.PSNR_sep;
					       delete[] assMetrics_HR.CC_sep;
					       delete[] assMetrics_HR.ERGAS_sep;
					       delete[] assMetrics_HR.UIQI_sep;
					       delete[] assMetrics_HR.DD_sep;
					       delete[] assMetrics_HR.AG_orig_sep;
					       delete[] assMetrics_HR.AG_rec_sep;
					       for(int iChZ=0; iChZ<glPrms.NChZ; iChZ++){
					     	  delete[] assMetrics_HR.DLambda_mat[iChZ];
					       }
					       delete[] assMetrics_HR.DLambda_mat;
					   }
                  }
			  }
			//==================================================================//
			//                  Add global parameters to report                 //
			//==================================================================//
			  if(my_rank==0 && iterMain==0){
				  report.addGlobalParams(&glPrms, &fSetting);
			  }

			//==================================================================//
			//              Write reconstructed high resolution data            //
			//==================================================================//
			  if(oSetting.writeImageFile){
				  if(iterMain==0){
					  paths.fname_ImZ = paths.dir_out+"/"+paths.fname_ImZ;
				  }
				  if(pSetting.store_patches_tmp_on_drive){
					  if(iterMain==0){
						  paths.fname_ImZ = paths.fname_ImZ+"_iter"+(char)iterMain;
					  }
					  MPI_Barrier(comm_busy);
					  double writeTSt = MPI_Wtime();
					  if(my_rank < numProcWrite){
						  ImZ->dataWriteParMPIIO(&report, &dSetting, &fSetting, &pSetting, &glPrms, &paths, comm_write);
					  }
					  if(my_rank==0){
						  ImZ->writeENVIHeader(&report, &dSetting, &fSetting, &glPrms, &paths);
					  }
					  MPI_Barrier(comm_busy);
					  glPrms.timeFileWrite = MPI_Wtime() - writeTSt;
				  }
			  }

			  MPI_Barrier(comm_busy);

			  if(my_rank==0){
				  //==================================================================//
				  //                     Save important information                   //
				  //==================================================================//
				  if(iterMain==0){
					  save_fusion_setup(&paths, &dSetting, &fSetting, &pSetting, &sSetting, &oSetting, &glPrms);
				  }
			  }
			  if( iterMain==fSetting.iterMain-1
			      && (    fSetting.set_neg_to_0==0
			          || (fSetting.set_neg_to_0==1  &&  fSetting.doFullImOptWithoutPatRec)
			          || (fSetting.set_neg_to_0==2  &&  !fSetting.LQ_post_opt_im)
			         )
			      && my_rank==0
			      && fSetting.evaluate
				){
				  cout << endl << endl
				       << "------------------------------------" << endl
					   << "| Set all negative values to zero  |" << endl
					   << "------------------------------------" << endl
					   << endl;
				  for(int iChY=0; iChY<ImZ->get_NCh(); iChY++){
					  ImZ->rasterBands[iChY]->bandDataMatD = ((ImZ->rasterBands[iChY]->bandDataMatD).array()>0).select(ImZ->rasterBands[iChY]->bandDataMatD,0);
				  }
				 //==================================================================//
				 //                   Assess reconstruction results                  //
				 //==================================================================//
				  final_evaluation = true;
				  // load high resolution reconstructed image
				  if(oSetting.writeImageFile && pSetting.store_patches_tmp_on_drive){
					  ImZ->dataRead(paths.fname_ImZ, &report);
				  }
				  // load high resolution reference image
				  SpEODataset *ImZ_ref;
				  ImZ_ref =	new SpEODataset(HR, imFlag_Z_ref);
				  ImZ_ref->dataRead(paths.fname_ImZ_ref, &report);
				  // possibly correct Z_ref and ImY to selected spectral and spatial subset
				  cutRelevantInput(HR, imFlag_Z_ref, ImZ_ref, &dSetting, &glPrms,false);
				  cutRelevantInput(LR, imFlag_Y, ImY, &dSetting, &glPrms,false);
				  SpEOAssessmentMetrics assMetrics_HR;
				  // do the assessment
				  PanSharp_Assessment_Eval(ImZ_ref, ImZ, ImY, &assMetrics_HR, &fSetting, &glPrms);
				  // write assessment results in report
				  report.addEvaluation(&assMetrics_HR, &dSetting, &fSetting, &oSetting, &glPrms, HR_ASSESSMENT);
				  // save assessment results
				  bool init_image_eval=true;
				  if(fSetting.evaluate_ImZ_init || !fSetting.doFullImOptWithoutPatRec){
					  init_image_eval=false;
				  }
				  save_evalResults(&paths, &glPrms, &fSetting, &assMetrics_HR, iterMain, fSetting.iterMain, false,fSetting.doFullImOptWithoutPatRec, init_image_eval, final_evaluation);
				  delete ImZ_ref;
				  delete[] assMetrics_HR.RMSE_sep;
				  delete[] assMetrics_HR.PSNR_sep;
				  delete[] assMetrics_HR.CC_sep;
				  delete[] assMetrics_HR.ERGAS_sep;
				  delete[] assMetrics_HR.UIQI_sep;
				  delete[] assMetrics_HR.DD_sep;
				  delete[] assMetrics_HR.AG_orig_sep;
				  delete[] assMetrics_HR.AG_rec_sep;
				  for(int iChZ=0; iChZ<glPrms.NChZ; iChZ++){
					  delete[] assMetrics_HR.DLambda_mat[iChZ];
				  }
				  delete[] assMetrics_HR.DLambda_mat;
			  }
			  if(oSetting.writeImageFile){
				  if(pSetting.store_patches_tmp_on_drive){
					  if (fSetting.LQ_post_opt_im){
						  if(my_rank < numProcWrite){
							  ImZ->dataWriteParMPIIO(&report, &dSetting, &fSetting, &pSetting, &glPrms, &paths, comm_write);
						  }
					  }
				  }else{
					  if(iterMain==fSetting.iterMain-1 || oSetting.writeImageFileAfterEveryIter){
						  string fname_ImZ_backup = paths.fname_ImZ;
						  if(oSetting.writeImageFileAfterEveryIter){
							  // modify image out filename intermediately for this iteration
							  stringstream fname_ImZ_tmp;
							  fname_ImZ_tmp << paths.fname_ImZ << "_iter" << iterMain << ".dat";
							  paths.fname_ImZ = fname_ImZ_tmp.str();
						  }
						  double writeTSt = MPI_Wtime();
						  if(my_rank < numProcWrite){
							  ImZ->dataWrite(&report, &dSetting, &fSetting, &pSetting, &glPrms, &paths, comm_write, ImZ);//ImZ_ref_tmp);
						  }
						  if(my_rank==0){
							  ImZ->writeENVIHeader(&report, &dSetting, &fSetting, &glPrms, &paths);
						  }
						  MPI_Barrier(comm_write);
						  glPrms.timeFileWrite = MPI_Wtime() - writeTSt;

						  if(oSetting.writeImageFileAfterEveryIter){
							  // restore original file name
							  paths.fname_ImZ = fname_ImZ_backup;
						  }
					  }
				  }
			  }

			//==================================================================//
			//                       Remove temporary files                     //
			//==================================================================//
			  if(my_rank==0){
				  if(pSetting.store_patches_tmp_on_drive && dSetting.delete_tmp_patch_folders){
					  cout << "delete tmp patch folders and files from disk.." << endl;
					  remove_dir(paths.dir_tmp_patches.c_str());
				  }
			  }
			  MPI_Barrier(comm_busy);
			  //==================================================================//
			  //      copy the results of ImZ to ImZ_init for next iteration      //
			  //==================================================================//
			  if(iterMain < fSetting.iterMain-1){
				  if(!fSetting.doFullImOptWithoutPatRec){
					  if(my_rank==0){
						  cout << "Copy the reconstruction results, ImZ, to ImZ_init for next iteration.." << endl;
					  }
					  for(int iChY=0; iChY<glPrms.NChY; iChY++){
						  ImZ_init->get_rasterBands()[iChY]->bandDataMatD = ImZ->get_rasterBands()[iChY]->bandDataMatD;
						  ImZ->get_rasterBands()[iChY]->bandDataMatD = SpEOMatrixD::Zero(ImZ->get_rasterBands()[iChY]->bandDataMatD.rows(),ImZ->get_rasterBands()[iChY]->bandDataMatD.cols());
					  }
				  }
			  }



			  if(my_rank==0){
				  // check stopping criteria
				  if(iterMain>1){
					  if(    ( SAM_patRec_in_all_iters(iterMain) > SAM_patRec_in_all_iters(iterMain-1)
							&& SAM_patRec_in_all_iters(iterMain) > SAM_patRec_in_all_iters(iterMain-2)
							&& CC_patRec_in_all_iters(iterMain) < CC_patRec_in_all_iters(iterMain-1)
							&& CC_patRec_in_all_iters(iterMain) < CC_patRec_in_all_iters(iterMain-2)
							&& UIQI_patRec_in_all_iters(iterMain) < UIQI_patRec_in_all_iters(iterMain-1)
							&& UIQI_patRec_in_all_iters(iterMain) < UIQI_patRec_in_all_iters(iterMain-2)
							&& RMSE_patRec_in_all_iters(iterMain) > RMSE_patRec_in_all_iters(iterMain-1)
							&& RMSE_patRec_in_all_iters(iterMain) > RMSE_patRec_in_all_iters(iterMain-2)
							&& PSNR_patRec_in_all_iters(iterMain) > PSNR_patRec_in_all_iters(iterMain-1)
							&& PSNR_patRec_in_all_iters(iterMain) > PSNR_patRec_in_all_iters(iterMain-2))
						 &&  ( SAM_fullImRec_in_all_iters(iterMain) > SAM_fullImRec_in_all_iters(iterMain-1)
						    && SAM_fullImRec_in_all_iters(iterMain) > SAM_fullImRec_in_all_iters(iterMain-2)
						    && CC_fullImRec_in_all_iters(iterMain) < CC_fullImRec_in_all_iters(iterMain-1)
						    && CC_fullImRec_in_all_iters(iterMain) < CC_fullImRec_in_all_iters(iterMain-2)
						    && UIQI_fullImRec_in_all_iters(iterMain) < UIQI_fullImRec_in_all_iters(iterMain-1)
						    && UIQI_fullImRec_in_all_iters(iterMain) < UIQI_fullImRec_in_all_iters(iterMain-2)
						    && RMSE_fullImRec_in_all_iters(iterMain) > RMSE_fullImRec_in_all_iters(iterMain-1)
						    && RMSE_fullImRec_in_all_iters(iterMain) > RMSE_fullImRec_in_all_iters(iterMain-2)
						    && PSNR_fullImRec_in_all_iters(iterMain) > PSNR_fullImRec_in_all_iters(iterMain-1)
						    && PSNR_fullImRec_in_all_iters(iterMain) > PSNR_fullImRec_in_all_iters(iterMain-2)  )
						){
						  stoppingCriteriaFulfilled = 1;
						  fSetting.iterMain = iterMain+1;

						  // code copied from the function save_fSetting(&paths, &fSetting):
						  string dir_fSetting = paths.dir_out + "/" + "fSetting";
						  cout << "write fusion settings to files in directory: " << endl << "     " << dir_fSetting << " .. ";
						  SpEOMatrixD tmp_mat = SpEOMatrixD::Constant(1,1,-99999);
						  string fname_tmp="";
						  tmp_mat(0,0) = fSetting.iterMain;
						  fname_tmp = dir_fSetting + "/" + "iterMain.csv";
						  write_Mat_to_CSV(&tmp_mat, fname_tmp.c_str());


						  // communicate to all processes in comm_busy
						  cout << "###################################################" << endl
							   << "# Stopping criteria for iterMain is fulfilled!!   #" << endl
							   << "# Loop terminates now!                            #" << endl
							   << "###################################################" << endl << endl;
					  }
				  }
			  }

			  MPI_Bcast(&stoppingCriteriaFulfilled, 1, MPI_INT, 0,comm_busy);
		  }
		  MPI_Barrier(comm_busy);
	  }

	//==================================================================//
	//                        Free MPI structures                       //
	//==================================================================//
	  MPI_Barrier(comm_busy);
	  MPI_Group_free(&group_busy);
	  MPI_Group_free(&mpi_group_orig);

	//==================================================================//
	//                          Clean up memory                         //
	//==================================================================//
//	  delete ImX;
	  delete ImX_LR;
//	  delete ImX_sim;
//	  delete ImX_sim_LR;
	  delete ImY;
          delete[] AbundanceMat;
	//==================================================================//
	//                        free communicators                        //
	//==================================================================//
	  MPI_Barrier(MPI_COMM_WORLD);
	  // clean up
	  MPI_Group_free(&group_write);
	  MPI_Comm_free(&comm_write);
	  MPI_Comm_free(&comm_busy);

	//==================================================================//
	//                         Finalize report                          //
	//==================================================================//
	  if(my_rank==0){
		  report.finalize(&glPrms);
		  save_fusion_setup(&paths, &dSetting, &fSetting, &pSetting, &sSetting, &oSetting, &glPrms);
	  }

	// clean up
	  delete ImX;
	  delete ImZ;

	  delete[] glPrms.Nm;
	  delete[] glPrms.Nm2;
	  for(int ii=0; ii<glPrms.NChX; ii++){
		  delete[] glPrms.km[ii];
	  }
	  for(int ii=0; ii<glPrms.Ng; ii++){
		  delete[] glPrms.km2[ii];
	  }
	  delete[] glPrms.km;
	  delete[] glPrms.km2;
	  delete[] glPrms.Nc_vec;
	  delete[] glPrms.P_lmd_idx_row;
	  for(int iG=0; iG<glPrms.Ng; iG++){
		  delete[] glPrms.P_lmd_idx_bl[iG];
	  }
	  delete[] glPrms.P_lmd_idx_bl;
	  delete[] glPrms.P_lmd_vecs;

#ifndef _OPENMP
	  delete[] glPrms.myChX;
	  delete[] glPrms.myBundle;
#endif

	//==================================================================//
	//                           Finalize MPI                           //
	//==================================================================//
	  MPI_Finalize();
	  return 0;
}
