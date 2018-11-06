/*
 * eval_alg.h
 *
 *  Created on: Nov 23, 2013
 *      Author: Claas Grohnfeldt
 */

#ifndef EVAL_ALG_H_
#define EVAL_ALG_H_

#include "includes.h"

// project headers
#include "dataIO.h"


using namespace std;
using namespace Eigen;



double Average_Gradient(SpEODataset* MSHR, int kCh, bool is_ref);
double Universal_Image_Quality_Index(SpEOMatrixD* ImZ_ref, int kCh, SpEOMatrixD* ImZ, int lCh);
double Spectral_Angle(SpEOMatrixD *ImZ_ref, SpEOMatrixD *ImZ);
double Degree_of_Distorion(SpEODataset *ImZ_ref, SpEODataset *ImZ, int kCh);
double Corr_Coef(SpEODataset *ImZ_ref, SpEODataset *ImZ, int kCh);
double Root_Mean_Square_Error(SpEODataset *ImZ_ref, SpEODataset *ImZ, int kCh);
void   PanSharp_Assessment_Eval(SpEODataset* ImZ_ref, SpEODataset* ImZ, SpEODataset* ImY, SpEOAssessmentMetrics *assMetrics, SpEOFusionSetting *sSetting, SpEOGlobalParams *glPrms);


#endif /* EVAL_ALG_H_ */
