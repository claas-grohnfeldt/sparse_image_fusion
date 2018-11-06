#ifndef JS_H
#define JS_H

#include "includes.h"

// project includes
#include "dataIO.h"
#include "auxFcts.h"
#include "filter.h"

using namespace Eigen;
using namespace std;

/**
 * Forward Backward Splitting (FBS) Solver
 * (including FISTA)
 */

// Forward Backward Splitting (FBS) Solver option struct and option enumerations
enum Stopping_criterion {
  FUNCTIONAL_DIFFERENCE,
  FIRST_ORDER_CONDITIONS
};
enum Prediction_rule {
  FISTA,
  FISTA_RESTART,
  STANDARD,
  STANDARD_AGGRESSIVE,
  LINE_SEARCH,
  LINE_SEARCH_SIMPLEX,
  LINE_SEARCH_ARMIJO
};
enum Backtracking_strategy {
  NO,
  DEC,
  DECPAR,
  DECEFF
};

enum Stepsize_update {
  INIT,
  LARGE,
  PSCL
};

struct FBSsolverOptions{
  FBSsolverOptions( int _Ndecomp, int _maxiter_in, int _maxiter_out, double _tol, Stopping_criterion _stop_crit, Prediction_rule _pred_rule, Backtracking_strategy _backtracking, Stepsize_update _stepsize_update, bool _stepsize_adaptive, bool _silent, double _gamma) : 
			  Ndecomp(_Ndecomp),
			  maxiter_in(_maxiter_in), 
			  maxiter_out(_maxiter_out),
			  tol(_tol),
			  stop_crit(_stop_crit),
			  pred_rule(_pred_rule),
			  backtracking(_backtracking),
			  stepsize_update(_stepsize_update),
			  stepsize_adaptive(_stepsize_adaptive),
			  silent(_silent),
			  gamma(_gamma) {}
  FBSsolverOptions() : Ndecomp(1), maxiter_in(1), maxiter_out(10000), tol(1e-2), stop_crit(FIRST_ORDER_CONDITIONS), pred_rule(STANDARD), backtracking(NO), stepsize_update(INIT), stepsize_adaptive(false), silent(false), gamma(1e2) {}

  int Ndecomp;
  int maxiter_in;
  int maxiter_out;
  double tol;
  Stopping_criterion stop_crit;
  Prediction_rule pred_rule;
  Backtracking_strategy backtracking;
  Stepsize_update stepsize_update;
  bool stepsize_adaptive;
  bool silent;
  double gamma;
};

// FBS Solver (parallel)
void FBSSolver(int &iter_need, double &rel_res, SpEOMatrixD &X, int N_all, SpEOMatrixD Phi, SpEOMatrixD y, double alpha, FBSsolverOptions opts, MPI_Comm comm_group, SpEOMatrixD &timestat, SpEOMatrixD &btstat);
// FBS Solver (sequential)
void FBSSolver(int &iter_need, double &rel_res, SpEOMatrixD &X, SpEOMatrixD Phi, SpEOMatrixD y, double alpha, FBSsolverOptions opts, SpEOMatrixD &timestat, SpEOMatrixD &btstat);


// Equation 1 Solver Classes:
class Eq1Op;
class Eq1MVar;
class Eq1SVar;
class Eq1SVar_sub;

class Eq1Op
{
public:
	Eq1Op(SpEOMatrixD* R, SpEOMatrixD gauss_filter, SpEOMatrixD filter_coeff, int fDS, double coeff1, double coeff2, double coeff3);
	
	SpEOMatrixD* R;
	SpEOMatrixD gauss_filter;
	SpEOMatrixD filter_coeff;
	int fDS;
	double coeff1;
	double coeff2;
	double coeff3;
	bool fastfilter;
};

class Eq1MVar
{
public:
	Eq1MVar(int N_upper, int N_center, int N_lower);
	void SubtractOperatorXvar(    Eq1Op op, SpEOMatrixD* var, MPI_Comm eq1_comm);
	void SubtractOperatorXvar_sub(Eq1Op op, SpEOMatrixD* var, MPI_Comm eq1_comm, SpEOMatrixD &EndmemberMat, SpEOVectorD &LAMBDA_ZY);
	void SetOperatorXvar(Eq1Op op, Eq1SVar*     var    , MPI_Comm eq1_comm);
	void SetOperatorXvar(Eq1Op op, Eq1SVar_sub* var_sub, MPI_Comm eq1_comm, SpEOMatrixD &EndmemberMat);
	void SetOperatorXvar(Eq1Op op, Eq1SVar_sub* var_sub, MPI_Comm eq1_comm, SpEOMatrixD &EndmemberMat, SpEOVectorD &LAMBDA_ZY);
	double NormSquared(MPI_Comm eq1_comm);
	~Eq1MVar();
	
	SpEOMatrixD* upper;
	SpEOMatrixD* center;
	SpEOMatrixD* lower;
	int N_upper;
	int N_center;
	int N_lower;
};

class Eq1SVar
{
public:
	Eq1SVar(int N);
	void SetAdjointOperatorXvar(Eq1Op op, Eq1MVar* m_var);
	double NormSquared(MPI_Comm eq1_comm);
	~Eq1SVar();
	
	SpEOMatrixD* var;
	int N;
};
class Eq1SVar_sub
{
public:
	Eq1SVar_sub(int N);
	void SetAdjointOperatorXvar(Eq1Op op, Eq1MVar* m_var, SpEOMatrixD &EndmemberMat, SpEOVectorD &LAMBDA_ZY);
	double NormSquared(MPI_Comm eq1_comm);
	~Eq1SVar_sub();
	
	SpEOMatrixD* var;
	int N;
};

// Equation 6 Solver Classes:
class Eq6Op;
class Eq6MVar;
class Eq6SVar;

class Eq6Op
{
public:
	Eq6Op(SpEOMatrixD* R, double coeff1, double coeff2);
	
	SpEOMatrixD* R;
	double coeff1;
	double coeff2;
};

class Eq6MVar
{
public:
	Eq6MVar();
	void SubtractOperatorXvar(Eq6Op op, SpEOVectorD* var);
	void SetOperatorXvar(Eq6Op op, Eq6SVar* var);
	double NormSquared();
	~Eq6MVar();
	
	SpEOVectorD upper;
	SpEOVectorD lower;
};

class Eq6SVar
{
public:
	Eq6SVar();
	void SetAdjointOperatorXvar(Eq6Op op, Eq6MVar* m_var);
	double NormSquared();
	~Eq6SVar();
	
	SpEOVectorD var;
};

// Equation 7 Solver Classes:
class Eq7Op;
class Eq7MVar;
class Eq7SVar;

class Eq7Op
{
public:
	Eq7Op(SpEOVectorD *dH, SpEOVectorD *dL, SpEOMatrixD* R, SpEOVectorD* P_lmd_vecs, SpEOMatrixI* P_lmd_idx_row, int **P_lmd_idx_bl, double *coeff1, double *coeff2, double coeff3);
	
	SpEOVectorD *dH;
	SpEOVectorD *dL;
	SpEOMatrixD* R;
	SpEOVectorD* P_lmd_vecs;
	SpEOMatrixI* P_lmd_idx_row;
	int **P_lmd_idx_bl;
	double *coeff1;
	double *coeff2;
	double coeff3;
};

class Eq7MVar
{
public:
	Eq7MVar(int N_g);
	void SubtractOperatorXvar(Eq7Op op, SpEOVectorD* var);
	void SetOperatorXvar(Eq7Op op, Eq7SVar* var);
	double NormSquared();
	~Eq7MVar();
	
	SpEOMatrixD* upper;
	SpEOMatrixD* center;
	SpEOMatrixD lower;
	int N_g;
};

class Eq7SVar
{
public:
	Eq7SVar(int N_g);
	void SetAdjointOperatorXvar(Eq7Op op, Eq7MVar* m_var);
	double NormSquared();
	~Eq7SVar();
	
	SpEOVectorD* var;
	int N_g;
};

// Equation 10 Solver Classes:
class Eq9Op;
class Eq9MVar;
class Eq9SVar;

class Eq9Op
{
public:
	Eq9Op(SpEOMatrixD *DH, SpEOMatrixD *DL, SpEOMatrixD* R, SpEOVectorD* P_lmd_vecs, SpEOMatrixI* P_lmd_idx_row, int **P_lmd_idx_bl, double *coeff1, double *coeff2, double coeff3);
	
	SpEOMatrixD *DH;
	SpEOMatrixD *DL;
	SpEOMatrixD* R;
	SpEOVectorD* P_lmd_vecs;
	SpEOMatrixI* P_lmd_idx_row;
	int **P_lmd_idx_bl;
	double *coeff1;
	double *coeff2;
	double coeff3;
};

class Eq9MVar
{
public:
	Eq9MVar(int N_g);
	void SubtractOperatorXvar(Eq9Op op, SpEOMatrixD* var);
	void SetOperatorXvar(Eq9Op op, Eq9SVar* var);
	double NormSquared();
	~Eq9MVar();
	
	SpEOMatrixD* upper;
	SpEOMatrixD* center;
	SpEOMatrixD lower;
	int N_g;
};

class Eq9SVar
{
public:
	Eq9SVar(int N_g);
	void SetAdjointOperatorXvar(Eq9Op op, Eq9MVar* m_var);
	double NormSquared();
	~Eq9SVar();
	
	SpEOMatrixD* var;
	int N_g;
};

// Adapted Least Squares solution with CGLS
void solve_equation3(int &iter, double &rel_res, SpEOMatrixD &Z, SpEOVectorD m, SpEOVectorD &m_out, SpEOMatrixD *DH, SpEOMatrixD *DL, SpEOVectorD* P_lmd_vecs, SpEOMatrixI* P_lmd_idx_row, int **P_lmd_idx_bl, SpEOMatrixD R, SpEOMatrixD X, SpEOMatrixD *Y, SpEOMatrixD *A, SpEOMatrixD* &A_out, SpEOMatrixD Z0, double lambda_X, double lambda_Y, int N_g, int maxiter, double tol_r, int fix_Alpha, bool fix_delta_m);
void solve_equation3_grad(SpEOMatrixD &Z, SpEOVectorD m, SpEOVectorD &m_out, SpEOMatrixD *DH, SpEOMatrixD *DL, SpEOVectorD* P_lmd_vecs, SpEOMatrixI* P_lmd_idx_row, int **P_lmd_idx_bl, SpEOMatrixD R, SpEOMatrixD X, SpEOMatrixD *Y, SpEOMatrixD *A, SpEOMatrixD* &A_out, SpEOMatrixD Z0, double lambda_X, double lambda_Y, int N_g, int maxiter, double tol_r, bool fix_Alpha, bool fix_delta_m);
  void solve_equation1(         int &iter, double &rel_res,                                                        SpEOMatrixD* &I_Z, SpEOMatrixD R, int filter_size, int fDS, int N_Y, int N_X, SpEOMatrixD* I_Z_tilde, SpEOMatrixD* I_X, SpEOMatrixD* I_Y,                     double lambda_X, double lambda_Y, int maxiter, double tol_r, MPI_Comm eq1_comm, int channels_per_core, bool SNR_normalization, bool balance_ImX_term_coef);
  void solve_equation1_unmixing(bool use_starting_value, int &iter, double &rel_res, SpEOMatrixD &EndmemberMat, SpEOMatrixD* &AbundanceMat, SpEOMatrixD* &I_Z, SpEOMatrixD R, int filter_size, int fDS, int N_Y, int N_X, SpEOMatrixD* I_Z_tilde, SpEOMatrixD* I_X, SpEOMatrixD* I_Y,                     double lambda_X, double lambda_Y, int maxiter, double tol_r, MPI_Comm eq1_comm, int channels_per_core, bool SNR_normalization, bool balance_ImX_term_coef, SpEOVectorD &LAMBDA_ZY);

// alias eq 6
void calcOptMeanDiffLS(int &iter, double &rel_res, SpEOVectorD &delta_m, SpEOVectorD* delta_m_0, SpEOMatrixD* R, SpEOVectorD* m, SpEOMatrixD* Z, SpEOMatrixD* X, double lambda_X, double lambda_Z, int maxiter, double tol_r);
// alias eq 6 (explicit solution)
SpEOMatrixD optMeanDiffLS_inverse(SpEOMatrixD* R, double lambda_X, double lambda_Z);
void calcOptMeanDiffLS_explicit(SpEOVectorD &delta_m, SpEOMatrixD* Inv, SpEOMatrixD* R, SpEOVectorD* m, SpEOMatrixD* Z, SpEOMatrixD* X, double lambda_X, double lambda_Z);
// alias eq 7
void calcOptCoeffCurrPatchLS(int &iter, double &rel_res, SpEOVectorD* &a, SpEOVectorD* a0, SpEOVectorD *dH, SpEOVectorD *dL, SpEOVectorD* m, SpEOVectorD* delta_m, SpEOMatrixD* Z, SpEOMatrixD* Y, SpEOMatrixD* X, SpEOMatrixD* R, SpEOVectorD* P_lmd_vecs, SpEOMatrixI* P_lmd_idx_row, int **P_lmd_idx_bl, double lambda_X, double lambda_Y, double lambda_Z, int N_g, int maxiter, double tol_r);
// alias eq 8
void calcOptCoeffResPFISTA(int* &iter, double* &rel_res, SpEOMatrixD* &A, SpEOMatrixD* A_0, SpEOMatrixD *DH, SpEOMatrixD *DL, SpEOVectorD* m, SpEOVectorD* delta_m, SpEOMatrixD* Z, SpEOMatrixD* Y, int **P_lmd_idx_bl, double lambda_A, double lambda_Y, double lambda_Z, int N_g, FBSsolverOptions opts, bool spectral_normalizer, bool &write_testset, int testnr);
// alias eq 10
void calcOptCoeffResLS(int &iter, double &rel_res, SpEOMatrixD* &A_res, SpEOMatrixD* A_res_0, SpEOMatrixD *DH, SpEOMatrixD *DL, SpEOVectorD* m, SpEOVectorD* delta_m, SpEOMatrixD* Z, SpEOMatrixD* Y, SpEOMatrixD* X, SpEOMatrixD* R, SpEOVectorD* P_lmd_vecs, SpEOMatrixI* P_lmd_idx_row, int **P_lmd_idx_bl, double lambda_X, double lambda_Y, double lambda_Z, int N_g, int maxiter, double tol_r);

void multiply_2D_with_3D_mat(SpEOMatrixD* &Mat_3D_output ,SpEOMatrixD &Mat_2D, SpEOMatrixD* &var_tmp, bool transpose_2D_mat);

#endif // JS_H
