/*
 * auxFcts.h
 *
 *  Created on: Apr 23, 2013
 *      Author: Claas Grohnfeldt
 */

#ifndef AUXFCTS_H_
#define AUXFCTS_H_

#include "includes.h"

// project headers
#include "dataIO.h"

/*****************************************************************************/
/*
  Function prototypes.
*/
/*****************************************************************************/

void calcGlobalParams(SpEOGlobalParams *glPrms, SpEOPaths *paths,
                      SpEODataIOSetting *dSetting, SpEOFusionSetting *fSetting,
                      SpEODataset *ImY, SpEODataset *ImX);
SpEOVectorI *kNearestIndices(int y0, int x0, int yLim, int xLim, int K);

// misc
double l1l2norm(SpEOMatrixD x);
SpEOMatrixD randn(
    int m, int N, double mean,
    double stddev);  // gives back matrix with standard normal distribution
SpEOVectorI randperm(
    int N);  // gives back random permutation of the numbers 0 to N-1
double spec_norm(SpEOMatrixD A);

#endif /* AUXFCTS_H_ */
