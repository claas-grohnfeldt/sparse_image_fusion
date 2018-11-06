/*
 * eval_alg.cpp
 *
 *  Created on: Nov 23, 2013
 *      Author: Claas Grohnfeldt
 */

#include "eval_alg.h"


using namespace std;
using namespace Eigen;

double Average_Gradient(SpEODataset *ImZ, int kCh, bool is_ref){
	int sizeU = ImZ->get_sizeU();
	int sizeV = ImZ->get_sizeV();
	double AG;
	if (is_ref){
		SpEOMatrixF diffY = ImZ->get_rasterBands()[kCh]->get_bandDataMat()->block(0,0,sizeU-1,sizeV-1) //get_bandDataMat()->block(0,0,sizeU-1,sizeV-1)
						   -ImZ->get_rasterBands()[kCh]->get_bandDataMat()->block(1,0,sizeU-1,sizeV-1);
		SpEOMatrixF diffX = ImZ->get_rasterBands()[kCh]->get_bandDataMat()->block(0,0,sizeU-1,sizeV-1)
						   -ImZ->get_rasterBands()[kCh]->get_bandDataMat()->block(0,1,sizeU-1,sizeV-1);
		diffY = diffY.cwiseAbs2();
		diffX = diffX.cwiseAbs2();
		diffY = (diffY+diffX)/2.0;
		diffY = diffY.cwiseSqrt();
		AG = (double)diffY.mean();
	}else{
		SpEOMatrixD diffY = ImZ->get_rasterBands()[kCh]->bandDataMatD.block(0,0,sizeU-1,sizeV-1) //get_bandDataMat()->block(0,0,sizeU-1,sizeV-1)
							   -ImZ->get_rasterBands()[kCh]->bandDataMatD.block(1,0,sizeU-1,sizeV-1);
		SpEOMatrixD diffX = ImZ->get_rasterBands()[kCh]->bandDataMatD.block(0,0,sizeU-1,sizeV-1)
						   -ImZ->get_rasterBands()[kCh]->bandDataMatD.block(0,1,sizeU-1,sizeV-1);
		diffY = diffY.cwiseAbs2();
		diffX = diffX.cwiseAbs2();
		diffY = (diffY+diffX)/2.0;
		diffY = diffY.cwiseSqrt();
		AG = diffY.mean();
	}
	return AG;
}

double Spectral_Angle(SpEODataset *ImZ_ref, SpEODataset *ImZ){
	int NChY = ImZ_ref->get_NCh();
	double tmp;
	SpEOMatrixD SAM_mat(ImZ_ref->get_sizeU(),ImZ_ref->get_sizeV());
	SpEOVectorD spectrumVec_o(NChY);
	SpEOVectorD spectrumVec_r(NChY);

	double pi = acos(-1.0);
	for(int y=0; y<ImZ_ref->get_sizeU(); y++){
		for(int x=0; x<ImZ_ref->get_sizeV(); x++){
			for(int k=0; k<NChY; k++){
				spectrumVec_o(k) = (double)(ImZ_ref->get_rasterBands()[k]->get_bandDataMat()->coeff(y,x));
				spectrumVec_r(k) = ImZ->get_rasterBands()[k]->bandDataMatD.coeff(y,x);
			}
			tmp = spectrumVec_o.dot(spectrumVec_r) / (spectrumVec_o.norm()*spectrumVec_r.norm());
			if(tmp>1)       tmp =  1.0;
			else if(tmp<-1) tmp = -1.0;			
			// transform rad to deg
			SAM_mat(y,x)= acos(tmp);
			if(SAM_mat(y,x) != SAM_mat(y,x)){ // this means if(SAM_mat(y,x)==nan)
				//cout << "Warning: Bad SAM value in pixel (y,x)=("<<y<<","<<x<<") go corrected from SAM(y,x)=" << SAM_mat(y,x) << " to SAM(y,x)=0" << endl;
				SAM_mat(y,x) = 0;
			}
			SAM_mat(y,x) = SAM_mat(y,x)*180.0/pi;
		}
	}
	//SAM(isnan(SAM))=0;
	return SAM_mat.mean();

}




double Universal_Image_Quality_Index(SpEOMatrixD *ImZ_ref, int kCh, SpEOMatrixD *ImZ, int lCh){
	int numEl = ImZ_ref[kCh].rows()*ImZ_ref[kCh].cols();
	double mean_o = ImZ_ref[kCh].mean();
	double mean_r = ImZ[lCh].mean();
	SpEOMatrixD im_o_zrMn = ImZ_ref[kCh].array()-mean_o;
	SpEOMatrixD im_r_zrMn = ImZ[lCh].array()-mean_r;
	im_o_zrMn.resize(numEl,1);
	im_r_zrMn.resize(numEl,1);
	double stdDev_o = sqrt(im_o_zrMn.squaredNorm()/((double)(numEl-1)));
	double stdDev_r = sqrt(im_r_zrMn.squaredNorm()/((double)(numEl-1)));
	double stdDev_or = (im_o_zrMn.transpose()*im_r_zrMn)(0,0)/((double)(numEl-1));
	double UIQI = (     stdDev_or      *        2*mean_o*mean_r          *           2*stdDev_o*stdDev_r          ) /
				  ( stdDev_o*stdDev_r  *  (mean_o*mean_o+mean_r*mean_r)  *  (stdDev_o*stdDev_o+stdDev_r*stdDev_r) );
	return UIQI;
}

double Degree_of_Distorion(SpEODataset *ImZ_ref, SpEODataset *ImZ, int kCh){
	// This following line could lead to a runtime error. to be checked!
//	SpEOMatrixD imDiff = *(ImZ_ref->get_rasterBands()[kCh]->get_bandDataMat()->cast<double>()) - *(ImZ->get_rasterBands()[kCh]->bandDataMatD);
	SpEOMatrixD imDiff =  (ImZ_ref->get_rasterBands()[kCh]->get_bandDataMat()->cast<double>()) - (ImZ->get_rasterBands()[kCh]->bandDataMatD);
	double DD = imDiff.cwiseAbs().mean();
	return DD;
}
double Corr_Coef(SpEODataset *ImZ_ref, SpEODataset *ImZ, int kCh){
	int numEl = ImZ_ref->get_sizeU()*ImZ_ref->get_sizeV();
	double mean_o = (ImZ_ref->get_rasterBands()[kCh]->get_bandDataMat()->cast<double>()).mean();				     // = mean1
	double mean_r = ImZ->get_rasterBands()[kCh]->bandDataMatD.mean();				     // = mean2
	SpEOMatrixD im_o_zrMn = ImZ_ref->get_rasterBands()[kCh]->get_bandDataMat()->cast<double>().array()-mean_o; // = d1_prel
	SpEOMatrixD im_r_zrMn = ImZ->get_rasterBands()[kCh]->bandDataMatD.array()-mean_r;  // = d2_prel
	im_o_zrMn.resize(numEl,1);																	 // = d1
	im_r_zrMn.resize(numEl,1);																	 // = d2
	double CC = (im_o_zrMn.transpose()*im_r_zrMn)(0,0)/(sqrt((im_o_zrMn.transpose()*im_o_zrMn)(0,0))*sqrt((im_r_zrMn.transpose()*im_r_zrMn)(0,0)));
	return CC;
}

double Root_Mean_Square_Error(SpEODataset *ImZ_ref, SpEODataset *ImZ, int kCh){
	// potential bug (possibly the asterisks have to be removed)
//	SpEOMatrixD imDiff = *(ImZ_ref->get_rasterBands()[kCh]->get_bandDataMat()->cast<double>())
//			             - *(ImZ->get_rasterBands()[kCh]->bandDataMatD);

//	cout << "bp2.2.1" << endl;
	SpEOMatrixD tmp1 = ImZ_ref->get_rasterBands()[kCh]->get_bandDataMat()->cast<double>();
	SpEOMatrixD tmp2 = ImZ->get_rasterBands()[kCh]->bandDataMatD;

//	cout << "tmp1.rows() = " << tmp1.rows() << ", tmp1.cols() = " << tmp1.cols() << endl;
//	cout << "tmp2.rows() = " << tmp2.rows() << ", tmp2.cols() = " << tmp2.cols() << endl;

	SpEOMatrixD imDiff = tmp1 - tmp2;
//	SpEOMatrixD imDiff = ImZ_ref->get_rasterBands()[kCh]->get_bandDataMat()->cast<double>()
//				             - ImZ->get_rasterBands()[kCh]->bandDataMatD;

//	cout << "bp2.2.2" << endl;


	double RMSE = imDiff.norm()/sqrt((double)(ImZ_ref->get_sizeU()*ImZ_ref->get_sizeV()));
	return RMSE;
};

void PanSharp_Assessment_Eval(SpEODataset* ImZ_ref,
							  SpEODataset* ImZ,
							  SpEODataset* ImY,
							  SpEOAssessmentMetrics *assMetrics,
							  SpEOFusionSetting *fSetting,
							  SpEOGlobalParams *glPrms){

	cout << "\n"
		<< "###########################################################" << "\n"
		<< "##                                                       ##" << "\n"
		<< "##                  assessment evaluation                ##" << "\n"
		<< "##                                                       ##" << "\n"
		<< "###########################################################"
		<< "\n";

	int NChZ  = ImZ_ref->get_NCh();
	VectorXd RMSE_Sep(NChZ);
	VectorXd PSNR_Sep(NChZ);
	VectorXd Mean_Sep(NChZ);
	MatrixXd Ratio(NChZ,1);
	VectorXd CC_Sep(NChZ);
	ArrayXd ERGAS_Sep(NChZ);
	ArrayXd UIQI_Sep(NChZ);
	ArrayXd DD_Sep(NChZ);
	ArrayXd SAM_Sep(NChZ);
	MatrixXd DLambda_mat(NChZ,NChZ);
	ArrayXd AG_Sep_Org(NChZ);
	ArrayXd AG_Sep_Rec(NChZ);

	 //cout << endl << "bp1" << endl;
	SpEOMatrixD *ImZ_o_mat = new SpEOMatrixD[NChZ];
	SpEOMatrixD *ImZ_r_mat = new SpEOMatrixD[NChZ];
	SpEOMatrixD *ImY_mat    = new SpEOMatrixD[NChZ];
	for(int k=0; k<NChZ; k++){
		// cout << endl << "bp1.1, k="<<k << endl;
		ImZ_o_mat[k] = ImZ_ref->get_rasterBands()[k]->get_bandDataMat()->cast<double>();
		// cout << endl << "bp1.2, k="<<k << endl;
		ImZ_r_mat[k] = SpEOMatrixD::Map(ImZ->get_rasterBands()[k]->get_bandDataMatD()->data(), ImZ->get_rasterBands()[k]->get_bandDataMatD()->rows(), ImZ->get_rasterBands()[k]->get_bandDataMatD()->cols());
//		ImZ_r_mat[k] << ImZ->get_rasterBands()[k]->get_bandDataMatD()->data();
		// cout << endl << "bp1.3, k="<<k << endl;
		ImY_mat[k]    = ImY->get_rasterBands()[k]->get_bandDataMat()->cast<double>();
		// cout << endl << "bp1.4, k="<<k << endl;
	}
	//cout << endl << "bp2" << endl;
	assMetrics->RMSE_sep    = new double[NChZ];
	assMetrics->PSNR_sep    = new double[NChZ];
	assMetrics->CC_sep      = new double[NChZ];
	assMetrics->ERGAS_sep   = new double[NChZ];
	assMetrics->UIQI_sep    = new double[NChZ];
	assMetrics->DD_sep      = new double[NChZ];
	assMetrics->AG_orig_sep = new double[NChZ];
	assMetrics->AG_rec_sep  = new double[NChZ];
	assMetrics->DLambda_mat = new double*[NChZ];
	for(int iChZ=0; iChZ<NChZ; iChZ++){
		assMetrics->DLambda_mat[iChZ] = new double[NChZ];
	}




	//cout << endl << "bp2.1" << endl;
//	RMSE & PSNR
	for(int k=0; k<NChZ; k++){
		RMSE_Sep(k) = Root_Mean_Square_Error(ImZ_ref, ImZ, k);
		//cout << endl << "k=" << k << " bp2.2" << endl;
		Mean_Sep(k)=ImZ->get_rasterBands()[k]->bandDataMatD.mean();
		//cout << endl << "k=" << k << " bp2.3" << endl;
		double max_tmp = ImZ_ref->get_rasterBands()[k]->get_bandDataMat()->maxCoeff();
		//cout << endl << "k=" << k << " bp2.4" << endl;
		PSNR_Sep(k) = 20.0*log10(max_tmp/RMSE_Sep(k));
		//cout << endl << "k=" << k << " bp2.5" << endl;
		assMetrics->RMSE_sep[k]=RMSE_Sep(k);
		//cout << endl << "k=" << k << " bp2.6" << endl;
		assMetrics->PSNR_sep[k]=PSNR_Sep(k);
		//cout << endl << "k=" << k << " bp2.7" << endl;
	}
	assMetrics->RMSE_mean = RMSE_Sep.mean();
	assMetrics->PSNR_mean = PSNR_Sep.mean();
//	CC
	//cout << endl << "bp3" << endl;
	for(int k=0; k<NChZ; k++){
		CC_Sep(k) = Corr_Coef(ImZ_ref, ImZ, k);
		assMetrics->CC_sep[k]=CC_Sep(k);	}
	assMetrics->CC_mean = CC_Sep.mean();
//	ERGAS
	//cout << endl << "bp4" << endl;
	int fDS = glPrms->fDS;
	double tmp0;
	for(int kCh=0; kCh<NChZ; kCh++){
		tmp0 = RMSE_Sep(kCh)/Mean_Sep(kCh);
		assMetrics->ERGAS_sep[kCh] = 100.0*tmp0/((double)fDS);
	}
	//cout << endl << "bp5" << endl;
	double tmp1 = RMSE_Sep.cwiseQuotient(Mean_Sep).norm()/sqrt((double)NChZ);
	assMetrics->ERGAS_mean=100.0*tmp1/((double)fDS);
//	UIQI
	//cout << endl << "bp6" << endl;
	for(int k=0; k<NChZ; k++){
		UIQI_Sep(k) = Universal_Image_Quality_Index(ImZ_o_mat, k, ImZ_r_mat, k);
		assMetrics->UIQI_sep[k]=UIQI_Sep(k);
	}
	assMetrics->UIQI_mean = UIQI_Sep.mean();

	//cout << endl << "bp7" << endl;
//	D_lambda
	for(int k=0; k<NChZ; k++){
		for(int l=0; l<NChZ; l++){
			double UIQI_ImZ_kl = Universal_Image_Quality_Index(ImZ_r_mat, k, ImZ_r_mat, l);
		    double UIQI_ImY_kl = Universal_Image_Quality_Index(ImY_mat,    k, ImY_mat,    l);
		    if(UIQI_ImZ_kl<0.0){
		    	UIQI_ImZ_kl = 0.0;
		    }
		    if(UIQI_ImY_kl<0.0){
		    	UIQI_ImY_kl = 0.0;
		    }
			DLambda_mat(k,l) = abs(UIQI_ImZ_kl-UIQI_ImY_kl);
			assMetrics->DLambda_mat[k][l]=DLambda_mat(k,l);
		}
	}
	assMetrics->DLambda_mean = DLambda_mat.mean();

	//cout << endl << "bp8" << endl;
//	DD
	for(int k=0; k<NChZ; k++){
			DD_Sep(k) = Degree_of_Distorion(ImZ_ref, ImZ, k);
			assMetrics->DD_sep[k]=DD_Sep(k);
	}
	assMetrics->DD_mean =DD_Sep.mean();
//	SAM
	//cout << endl << "bp9" << endl;
	assMetrics->SAM=Spectral_Angle(ImZ_ref, ImZ);
	 //cout << endl << "bp10" << endl;
//	Average Gradient (AG)
	for(int k=0; k<NChZ; k++){
		AG_Sep_Org(k) = Average_Gradient(ImZ_ref, k, true);
		assMetrics->AG_orig_sep[k]=AG_Sep_Org(k);
		AG_Sep_Rec(k) = Average_Gradient(ImZ, k, false);
		assMetrics->AG_rec_sep[k]=AG_Sep_Rec(k);
	}
	 //cout << endl << "bp11" << endl;
	assMetrics->AG_orig_mean = AG_Sep_Org.mean();
	assMetrics->AG_rec_mean = AG_Sep_Rec.mean();

	delete[] ImZ_o_mat;
	delete[] ImZ_r_mat;
	delete[] ImY_mat;

	//cout << endl << "bp12" << endl;
	return;
};
