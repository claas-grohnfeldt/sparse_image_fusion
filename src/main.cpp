/*#############################################################################
##                                                                           ##
##    File:        main.cpp                                                  ##
##                                                                           ##
##    Author:      Claas Grohnfeldt                                          ##
##    Contact:     claas.grohnfeldt@gmail.com                                ##
##                                                                           ##
##    Version:     3.0 (12/2018)                                             ##
##                                                                           ##
##    Last modification: December 23, 2018                                   ##
##                                                                           ##
##    Copyright:   Claas Grohnfeldt                                          ##
##                                                                           ##
##    Description: Main file for the "sparse image fusion" project software  ##
##                 suite, which is based on an implementation of the         ##
##                 J-SparseFI-HM (Jointly Sparse Fusion of Hyper- and Multi- ##
##                 spectral Imagery).                                        ##
##                                                                           ##
#############################################################################*/

#include "main.h"
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
      int my_rank; int my_processes;
      MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
      MPI_Comm_size(MPI_COMM_WORLD, &my_processes);

  	//==================================================================//
  	//                         Get user settings                        //
  	//==================================================================//
	  SpEOPaths paths;
      SpEODataIOSetting dSetting;
      SpEOFusionSetting fSetting;
      SpEOOutputSetting oSetting;
      SpEOSolverSetting sSetting;
      SpEOParallelSetting pSetting;

      getUserSettings(&paths, &dSetting, &fSetting, &oSetting, &sSetting, &pSetting, argc, argv);

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
	  ImZ_init->dataRead(paths.fname_ImZ_init, &report);

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

	//=======================================================================================//
	//  Read subspace transformation matrix file (e.g. endmember matrix in case of unmixing) //
	//=======================================================================================//
	  // declare, initialize and possibly load endmember matrix
	  SpEOMatrixD SubspaceTransformMat(ImY->get_NCh(),fSetting.subspace_dim);
	  if(fSetting.subspace_transform_type == "none"){
	     SubspaceTransformMat =  SpEOMatrixD::Identity(ImY->get_NCh(),ImY->get_NCh());
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
	//   Read spectral response functions file                          //
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
	 MPI_Barrier(MPI_COMM_WORLD);

	 //===================================================================//
	 //              Create MPI groups and sub-communicators              //
	 //===================================================================//
	  if(my_rank==0){
		  cout << "Generate MPI communicators .." << endl;
	  }

	  MPI_Group mpi_group_orig;
	  MPI_Comm_group(MPI_COMM_WORLD, &mpi_group_orig);

	  glPrms.numPatchGroups = my_processes / pSetting.numProcPerPatch;
	  
	  //################
	  //## comm_busy  ##
	  //################
	  MPI_Group group_busy;
#ifdef _OPENMP
	  int numProcBusy = glPrms.numPatchGroups;
#else
	  
	  int numProcBusy = pSetting.numProcPerPatch * min(glPrms.numPatchGroups, glPrms.NP_sub);
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
		  }

		  int iterMain;
		  for(iterMain=0; iterMain<fSetting.iterMain; iterMain++){
			  if(my_rank==0){
				  cout << endl
					   << "##########################################" << endl
					   << "#                                        #" << endl
					   << "#            iteration "<<iterMain<<"                 #" << endl
					   << "#                                        #" << endl
					   << "##########################################" << endl << endl;
			  }
			  if (!fSetting.doFullImOptWithoutPatRec){
					  JSparseFIHM_alg(SubspaceTransformMat,AbundanceMat,iterMain, fSetting.iterMain, &paths, &dSetting, &fSetting, &oSetting, &sSetting, &pSetting, &glPrms,
									ImX, ImX_LR, ImY, ImZ, &SRF, comm_busy, group_busy, ImZ_init, report);
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
										  ImZ->dataRead(paths.fname_ImZ_out, &report);
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
			  //#   and PSFs                                            #
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
					     	  ImZ->dataRead(paths.fname_ImZ_out, &report);
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
					  paths.fname_ImZ_out = paths.dir_out+"/"+paths.fname_ImZ_out;
				  }
				  if(pSetting.store_patches_tmp_on_drive){
					  if(iterMain==0){
						  paths.fname_ImZ_out = paths.fname_ImZ_out+"_iter"+(char)iterMain;
					  }
					  if(my_rank < numProcWrite){
						  ImZ->dataWriteParMPIIO(&report, &dSetting, &fSetting, &pSetting, &glPrms, &paths, comm_write);
					  }
					  if(my_rank==0){
						  ImZ->writeENVIHeader(&report, &dSetting, &fSetting, &glPrms, &paths);
					  }
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
					  ImZ->dataRead(paths.fname_ImZ_out, &report);
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
						  string fname_ImZ_out_backup = paths.fname_ImZ_out;
						  if(oSetting.writeImageFileAfterEveryIter){
							  // modify image out filename intermediately for this iteration
							  stringstream fname_ImZ_out_tmp;
							  fname_ImZ_out_tmp << paths.fname_ImZ_out << "_iter" << iterMain << ".dat";
							  paths.fname_ImZ_out = fname_ImZ_out_tmp.str();
						  }
						  if(my_rank < numProcWrite){
							  ImZ->dataWrite(&report, &dSetting, &fSetting, &pSetting, &glPrms, &paths, comm_write, ImZ);
						  }
						  if(my_rank==0){
							  ImZ->writeENVIHeader(&report, &dSetting, &fSetting, &glPrms, &paths);
						  }
						  if(oSetting.writeImageFileAfterEveryIter){
							  // restore original file name
							  paths.fname_ImZ_out = fname_ImZ_out_backup;
						  }
					  }
				  }
			  }

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
	  delete ImX_LR;
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

	//==================================================================//
	//                           Finalize MPI                           //
	//==================================================================//
	  MPI_Finalize();
	  return 0;
}
