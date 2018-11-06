/*
 * JSparseFI_alg.h
 *
 *  Created on: Feb 20, 2014
 *      Author: Claas Grohnfeldt
 */

#ifndef JSPARSEFI_ALG_H_
#define JSPARSEFI_ALG_H_

#include "includes.h"

// project headers
#include "dataIO.h"
#include "JS.h"
#include "mpi_counter.h"
#include "auxFcts.h"

//void JSparseFI_alg(	int iterMain,
//					int numIterMain,
//					SpEOPaths *paths,
//					SpEODataIOSetting *dSet,
//					SpEOFusionSetting *fSet,
//					SpEOOutputSetting *oSet,
//					SpEOSolverSetting *sSet,
//					SpEOParallelSetting *pSet,
//					SpEOGlobalParams *glPrms,
//					SpEODataset *ImX,
//					SpEODataset *ImX_LR,
//					SpEODataset *ImX_sym,
//					SpEODataset *ImX_sym_LR,
//					SpEODataset *ImY,
//					SpEODataset *ImZ,
//					SpEOMatrixD *SRF,
//					MPI_Comm comm_busy,
//					MPI_Group group_busy,
//					SpEODataset *ImZ_init,
//					SpEOReport &report);
void JSparseFI_alg(	                SpEOMatrixD &EndmemberMat, 
					SpEOMatrixD* &AbundanceMat, 
                                        int iterMain,
					int numIterMain,
					SpEOPaths *paths,
					SpEODataIOSetting *dSet,
					SpEOFusionSetting *fSet,
					SpEOOutputSetting *oSet,
					SpEOSolverSetting *sSet,
					SpEOParallelSetting *pSet,
					SpEOGlobalParams *glPrms,
					SpEODataset *ImX,
					SpEODataset *ImX_LR,
//					SpEODataset *ImX_sym,
//					SpEODataset *ImX_sym_LR,
					SpEODataset *ImY,
					SpEODataset *ImZ,
					SpEOMatrixD *SRF,
					MPI_Comm comm_busy,
					MPI_Group group_busy,
					SpEODataset *ImZ_init,
					SpEOReport &report);

void full_im_optimization_LS(SpEOReport &report, SpEOMatrixD &EndmemberMat, SpEOMatrixD* &AbundanceMat, SpEOGlobalParams *glPrms, SpEOFusionSetting *fSetting, SpEOSolverSetting *sSetting, SpEODataIOSetting *dSetting, SpEODataset *ImX, SpEODataset *ImY, SpEODataset *ImZ, SpEOMatrixD *SRF, MPI_Comm comm_busy, MPI_Group group_busy);
 
#endif /* JSPARSEFI_ALG_H_ */
