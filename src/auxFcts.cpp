/*
 * auxFcts.cpp
 *
 *  Created on: Apr 23, 2013
 *      Author: Claas Grohnfeldt
 */

#include "auxFcts.h"

using namespace Eigen;
using namespace std;


void calcGlobalParams(SpEOGlobalParams *glPrms,
					  SpEOPaths *paths,
					  SpEODataIOSetting *dSetting,
					  SpEOFusionSetting *fSetting,
					  SpEODataset *ImY,
					  SpEODataset *ImX) {

	int my_rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

	if(my_rank==0){
		cout << "\n"
		<< "###########################################################\n"
		<< "##                                                       ##\n"
		<< "##              calculate global parameters              ##\n"
		<< "##                                                       ##\n"
		<< "###########################################################"
		<< "\n" << "\n";
	}

    int resultlen;
	char *myProcName = new char[MPI_MAX_PROCESSOR_NAME];
	MPI_Get_processor_name(myProcName, &resultlen);
	glPrms->myProcName.assign(myProcName,resultlen);

	glPrms->NChX = ImX->get_NCh();
	glPrms->NChY = ImY->get_NCh();
	glPrms->NChZ = dSetting->chBundleLast-dSetting->chBundleFirst+1;
	// check the necessary condition of dimH being an even multiple of dimL
	float tol = 0.0000001;
	if (((ImX->get_sizeU() / ImY->get_sizeU()) % 1) > tol
			|| ((ImX->get_sizeV() / ImY->get_sizeV()) % 1) > tol
			|| abs(
					ImX->get_sizeV() / ImY->get_sizeV()
							- ImX->get_sizeU() / ImY->get_sizeU()) > tol) {

		cerr << "\n\n ERROR: Dimension of pan image must be an even multiple of the dimension of the low-resolution multispectral image!" << "\n\n";
		exit(2);

		return;
	} else
	glPrms->fDS = ImX->get_sizeU() / ImY->get_sizeU();
	glPrms->sizeUL = ImY->get_sizeU();
	glPrms->sizeVL = ImY->get_sizeV();
	glPrms->sizeUH = ImX->get_sizeU();
	glPrms->sizeVH = ImX->get_sizeV();

	if(dSetting->uLFirst<0){
		if(my_rank==0){
			cout << "\n uLFirst got corrected from uLFirst=" << dSetting->uLFirst << " to uLFirst=0!" << endl;
		}
		dSetting->uLFirst = 0;
	}
	if(dSetting->vLFirst<0){
		if(my_rank==0){
			cout << "\n vLFirst got corrected from vLFirst=" << dSetting->vLFirst << " to vLFirst=0!" << endl;
		}
		dSetting->vLFirst = 0;
	}
	if(dSetting->uLLast>=ImY->get_sizeU()){
		if(my_rank==0){
			cout << "\n uLLast got corrected from uLLast=" << dSetting->uLLast << " to uLLast=" << ImY->get_sizeU()-1 << "=ImY->get_sizeU()-1!" << endl;
		}
		dSetting->uLLast = ImY->get_sizeU()-1;
	}
	if(dSetting->vLLast>=ImY->get_sizeV()){
		if(my_rank==0){
			cout << "\n vLLast got corrected from vLLast=" << dSetting->vLLast << " to vLLast=" << ImY->get_sizeV()-1 << "=ImY->get_sizeV()-1!" << endl;
		}
		dSetting->vLLast = ImY->get_sizeV()-1;
	}
	glPrms->sizeUL_red = dSetting->uLLast-dSetting->uLFirst+1;
	glPrms->sizeVL_red = dSetting->vLLast-dSetting->vLFirst+1;
	glPrms->sizeUH_red = glPrms->sizeUL_red*glPrms->fDS;
	glPrms->sizeVH_red = glPrms->sizeVL_red*glPrms->fDS;

	int a     = fSetting->patchsize-fSetting->overlap;
	glPrms->NPU   = ceil(double((ImY->get_sizeU() - fSetting->patchsize)) / a) + 1; // number of patches in a column
	glPrms->NPV   = ceil(double((ImY->get_sizeV() - fSetting->patchsize)) / a) + 1; // number of patches in a row
	glPrms->NP = glPrms->NPU*glPrms->NPV;

	glPrms->uPFirst = dSetting->uLFirst/a;          if(glPrms->uPFirst>=glPrms->NPU){ glPrms->uPFirst = glPrms->NPU-1; }
	glPrms->uPLast = dSetting->uLLast/a;            if(glPrms->uPLast>=glPrms->NPU) { glPrms->uPLast  = glPrms->NPU-1; }
	glPrms->vPFirst = dSetting->vLFirst/a;          if(glPrms->vPFirst>=glPrms->NPV){ glPrms->vPFirst = glPrms->NPV-1; }
	glPrms->vPLast = dSetting->vLLast/a;            if(glPrms->vPLast>=glPrms->NPV) { glPrms->vPLast  = glPrms->NPV-1; }

	if(glPrms->uPFirst<0 || glPrms->uPFirst >= glPrms->NPU){
		if(my_rank==0){ cout << "\n WARNING: uPFirst got corrected from uPFirst=" << glPrms->uPFirst << " to uPFirst=0" << endl;}
		glPrms->uPFirst = 0;
	}
	if(glPrms->uPLast<0 || glPrms->uPLast<glPrms->uPFirst || glPrms->uPLast >= glPrms->NPU){
		if(my_rank==0){ cout << "\n WARNING: uPLast got corrected from uPLast=" << glPrms->uPLast << " to uPLast=" << glPrms->NPU << "==NPU!" << endl;}
		glPrms->uPLast = glPrms->NPU-1;
	}
	if(glPrms->vPFirst<0 || glPrms->vPFirst >= glPrms->NPV){
		if(my_rank==0){ cout << "\n WARNING: vPFirst got corrected from vPFirst=" << glPrms->vPFirst << " to vPFirst=0" << endl;}
		glPrms->vPFirst = 0;
	}
	if(glPrms->vPLast<0 || glPrms->vPLast<glPrms->vPFirst || glPrms->vPLast >= glPrms->NPV){
		if(my_rank==0){ cout << "\n WARNING: vPLast got corrected from vPLast=" << glPrms->vPLast << " to vPLast=" << glPrms->NPV << "==NPV!" << endl;}
		glPrms->vPLast = glPrms->NPV-1;
	}

	glPrms->NPV_sub = glPrms->vPLast-glPrms->vPFirst+1;
	glPrms->NPU_sub = glPrms->uPLast-glPrms->uPFirst+1;
	glPrms->NP_sub = glPrms->NPU_sub*glPrms->NPV_sub;

	if (fSetting->NDP > glPrms->NP){
		if(my_rank==0){
			cout << "\n     WARNING: NDP got corrected from NDP=" << fSetting->NDP << " to NDP=" << glPrms->NP << "==NP!" << endl;
		}
		fSetting->NDP = glPrms->NP;
	}


	// in case of shortage in memory: possibly replace by operations:
	SpEOVectorI idxPUH = SpEOVectorI::Zero(glPrms->NPU, 1);
	SpEOVectorI idxPVH = SpEOVectorI::Zero(glPrms->NPV, 1);
	//SpEOVectorI idxPUL = SpEOVectorI::Zero(glPrms->NPU, 1);
	int uP;
	for(uP=0; uP<glPrms->NPU-1; uP++){
		idxPUH(uP) = glPrms->fDS*a*uP;
//		idxPUL(uP) = a*uP;
	}
	int pszH  = fSetting->patchsize*glPrms->fDS;
	idxPVH(glPrms->NPV-1) = ImX->get_sizeV()-pszH;
	idxPUH(glPrms->NPU-1) = ImX->get_sizeU()-pszH;
	int vP;
	for(vP=0; vP<glPrms->NPV-1; vP++){
		idxPVH(vP) = glPrms->fDS*a*vP;
//		idxPVL(vP) = a*vP;
	}
	glPrms->idxPUH = idxPUH;
	glPrms->idxPVH = idxPVH;

}

SpEOVectorI* kNearestIndices(int u0, int v0, int uLim, int vLim, int K) {
	SpEOVectorI* indices = new SpEOVectorI[K];
	//
	indices[0] = SpEOVectorI::Zero(2);
	indices[0](0) = u0;
	indices[0](1) = v0;
	int pntCnt = 0;
	// (du, dv) is a vector - direction in which we move right now
	int du = 1;
	int dv = 0;
	// length of current segment
	int segment_length = 1;

	// current position (u, v) and how much of current segment we passed
	int u = u0;
	int v = v0;
	int segment_passed = 0;
	int k = 1;
	while (k < K) {
		// make a step, add 'direction' vector (du, dv) to current position (u, v)
		u += du;
		v += dv;
		++segment_passed;
		//cout << u <<  " " << v << endl;
		if (segment_passed == segment_length) {
			// done with current segment
			segment_passed = 0;
			// 'rotate' directions
			int buffer = du;
			du = -dv;
			dv = buffer;
			// increase segment length if necessary
			if (dv == 0) {
				++segment_length;
			}
		}
		if (u >= 0 && u <= uLim && v >= 0 && v <= vLim) {
			indices[k] = SpEOVectorI::Zero(2);
			indices[k](0) = u;
			indices[k](1) = v;
			pntCnt++;
			k++;
		}
	}
	return indices;
}

double l1l2norm(SpEOMatrixD x) {
    return x.rowwise().norm().lpNorm<1>();
}

double rand_normal(double mean, double stddev)
{//Box muller method
    static double n2 = 0.0;
    static int n2_cached = 0;
    if (!n2_cached)
    {
        double x, y, r;
        do
        {
            x = 2.0*rand()/RAND_MAX - 1;
            y = 2.0*rand()/RAND_MAX - 1;

            r = x*x + y*y;
        }
        while (r == 0.0 || r > 1.0);
        {
            double d = sqrt(-2.0*log(r)/r);
            double n1 = x*d;
            n2 = y*d;
            double result = n1*stddev + mean;
            n2_cached = 1;
            return result;
        }
    }
    else
    {
        n2_cached = 0;
        return n2*stddev + mean;
    }
}

SpEOMatrixD randn(int m, int N, double mean, double stddev) {
  SpEOMatrixD mat = SpEOMatrixD::Zero(m,N);

  int i,j;
  for (i=0; i<m; i++) {
    for (j=0; j<N; j++) {
      //mat(i,j) = distribution(generator);
      mat(i,j) = rand_normal(mean, stddev);
    }
  }
  return mat;
}

SpEOVectorI randperm(int N) {
  int* vec = new int[N];
  int i;
  for (i = 0; i < N; i++) {
    vec[i] = i;
  }  
  random_shuffle(vec, vec + N);
  SpEOVectorI vec_out = SpEOVectorI::Zero(N);
  for (i=0; i < N; i++) {
    vec_out(i) = vec[i];
  }
  delete [] vec;
  return vec_out;
}

double spec_norm(SpEOMatrixD A) {
  //SpEOMatrixD S = A*A.transpose(); // use selfadjoint property of spectral norm and use the computationally more efficient way (smaller matrix)
  MatrixXd S = A*A.transpose(); // use selfadjoint property of spectral norm and use the computationally more efficient way (smaller matrix)

  // old Eigen Library version (eigen-eigen-2249f9c22fe8)
  // PowerIteration<SpEOMatrixD> e2 = PowerIteration<SpEOMatrixD>();
  // e2.setPrecision(1e-2);
  // e2.setMaxIter(1e2);
  // e2.compute(S,1);
  // SpEOVectorD esev2 = e2.eigenvalues().transpose();
  // return sqrt(esev2(0));

  // newer Eigen Library versions
  // EigenSolver<SpEOMatrixD> es;
  // es.compute(S, /*computeEigenvectors = */ false);
  // VectorXcd eivals = es.eigenvalues().transpose();
  // double esev2 = real(eivals(0));
  // return sqrt(esev2);

  // EigenSolver<MatrixXd> es;
  // es.compute(S, /*computeEigenvectors = */ false);
  EigenSolver<MatrixXd> es(S, /*computeEigenvectors = */ false);
  //cout << "  es.info() = " << endl << es.info() << endl;
  double largest_eival_real;
  if(!es.info()){
	  complex<double> largest_eival = es.eigenvalues()[0];
  	  largest_eival_real = real(largest_eival);
  	  //cout << " real(largest_eival) = " << largest_eival_real << endl;
  }else{
	  largest_eival_real = 0.1;
	  cout << "Warning: eigenvalue computation for spectral norm normalization did not succeed. Normalize by static value instead" << endl;
  }
  return sqrt(largest_eival_real);
}