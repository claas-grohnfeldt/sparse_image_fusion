#ifndef FILTER_H
#define FILTER_H

#include "includes.h"

// project includes
#include "dataIO.h"

using namespace Eigen;
using namespace std;

// auxiliary filter operations
void create_Gaussian_filter(SpEOMatrixD &filter, int filter_size);
void calc_filter_boundary_coeff(SpEOMatrixD &filter_coeff, SpEOMatrixD filter, int fDS);

void apply_filter(SpEOMatrixD &image, SpEOMatrixD filter, int fDS);
void apply_filter_adjoint(SpEOMatrixD &image, SpEOMatrixD filter, int fDS);
void apply_filter_boundary_coeff(SpEOMatrixD &image, SpEOMatrixD filter_coeff, int filter_size, int fDS); // filter_size even <=> image sizes even
void select(SpEOMatrixD &LRp, SpEOMatrixD HRp, int fDS);
void distribute(SpEOMatrixD &HRp, SpEOMatrixD LRp, int fDS);

// filter (alias X*BS)
void standard_filter(SpEOMatrixD &LRp, SpEOMatrixD HRp, SpEOMatrixD filter, SpEOMatrixD filter_coeff, int fDS);
void fast_filter(SpEOMatrixD &LRp, SpEOMatrixD HRp, SpEOMatrixD filter, SpEOMatrixD filter_coeff, int fDS);

// adjoint filter (alias Y*(BS)^T)
void standard_adjoint_filter(SpEOMatrixD &HRp, SpEOMatrixD LRp, SpEOMatrixD filter, SpEOMatrixD filter_coeff, int fDS);

#endif // FILTER_H