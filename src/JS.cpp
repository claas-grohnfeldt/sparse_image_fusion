
#include "JS.h"

/**
 * Equation 1 Solver Class Functions
 **/
Eq1Op::Eq1Op(SpEOMatrixD* R_m, SpEOMatrixD gauss_filter_m,
             SpEOMatrixD filter_coeff_m, int fDS_m, double coeff1_m,
             double coeff2_m, double coeff3_m) {
  R = R_m;
  gauss_filter = gauss_filter_m;
  filter_coeff = filter_coeff_m;
  fDS = fDS_m;
  coeff1 = coeff1_m;
  coeff2 = coeff2_m;
  coeff3 = coeff3_m;
  fastfilter = true;
  return;
}

Eq1MVar::Eq1MVar(int N_upper_m, int N_center_m, int N_lower_m) {
  upper = new SpEOMatrixD[N_upper_m];
  center = new SpEOMatrixD[N_center_m];
  lower = new SpEOMatrixD[N_lower_m];
  N_upper = N_upper_m;
  N_lower = N_lower_m;
  N_center = N_center_m;
  return;
}

double Eq1MVar::NormSquared(MPI_Comm eq1_comm) {
  double norm_squared = 0.0;
  int channel_Y, channel_X;
  for (channel_Y = 0; channel_Y < N_upper; channel_Y++) {
    norm_squared += pow(lower[channel_Y].norm(), 2);
    norm_squared += pow(upper[channel_Y].norm(), 2);
  }
  MPI_Allreduce(MPI_IN_PLACE, &norm_squared, 1, MPI_DOUBLE, MPI_SUM, eq1_comm);
  for (channel_X = 0; channel_X < N_center; channel_X++) {
    norm_squared += pow(center[channel_X].norm(), 2);
  }
  return norm_squared;
}

Eq1MVar::~Eq1MVar() {
  delete[] upper;
  delete[] center;
  delete[] lower;
  return;
}

void Eq1MVar::SubtractOperatorXvar(Eq1Op op, SpEOMatrixD* m_var,
                                   MPI_Comm eq1_comm) {
  int channel_X, channel_Y;
  SpEOMatrixD* TMP = new SpEOMatrixD[N_center];
  for (channel_X = 0; channel_X < N_center; channel_X++) {
    TMP[channel_X] = SpEOMatrixD::Zero((center[channel_X]).rows(),
                                       (center[channel_X]).cols());
  }
  for (channel_Y = 0; channel_Y < N_upper; channel_Y++) {
    upper[channel_Y] -= op.coeff1 * (m_var[channel_Y]);

    // calculate R*var
    for (channel_X = 0; channel_X < N_center; channel_X++) {
      TMP[channel_X] +=
          op.coeff2 * ((*op.R)(channel_X, channel_Y)) * (m_var[channel_Y]);
    }
    // calculate var*BS
    SpEOMatrixD tmp_image_HR = m_var[channel_Y];
    SpEOMatrixD tmp_image_LR;
    if (op.fastfilter == true) {
      fast_filter(tmp_image_LR, tmp_image_HR, op.gauss_filter, op.filter_coeff,
                  op.fDS);
    } else {
      standard_filter(tmp_image_LR, tmp_image_HR, op.gauss_filter,
                      op.filter_coeff, op.fDS);
    }
    lower[channel_Y] -= op.coeff3 * tmp_image_LR;
  }
  for (channel_X = 0; channel_X < N_center; channel_X++) {
    MPI_Allreduce(MPI_IN_PLACE, TMP[channel_X].data(),
                  (center[channel_X]).rows() * (center[channel_X]).cols(),
                  MPI_DOUBLE, MPI_SUM, eq1_comm);
    center[channel_X] -= TMP[channel_X];
  }
  return;
}

void Eq1MVar::SetOperatorXvar(Eq1Op op, Eq1SVar* m_var, MPI_Comm eq1_comm) {
  int channel_X, channel_Y;
  for (channel_X = 0; channel_X < N_center; channel_X++) {
    center[channel_X] = SpEOMatrixD::Zero((center[channel_X]).rows(),
                                          (center[channel_X]).cols());
  }
  for (channel_Y = 0; channel_Y < N_upper; channel_Y++) {
    upper[channel_Y] = op.coeff1 * (m_var->var[channel_Y]);
    // calculate R*var
    for (channel_X = 0; channel_X < N_center; channel_X++) {
      center[channel_X] +=
          op.coeff2 * ((*op.R)(channel_X, channel_Y)) * (m_var->var[channel_Y]);
    }
    // calculate var*BS
    SpEOMatrixD tmp_image_HR = m_var->var[channel_Y];
    SpEOMatrixD tmp_image_LR;
    if (op.fastfilter == true) {
      fast_filter(tmp_image_LR, tmp_image_HR, op.gauss_filter, op.filter_coeff,
                  op.fDS);
    } else {
      standard_filter(tmp_image_LR, tmp_image_HR, op.gauss_filter,
                      op.filter_coeff, op.fDS);
    }
    lower[channel_Y] = op.coeff3 * tmp_image_LR;
  }
  for (channel_X = 0; channel_X < N_center; channel_X++) {
    MPI_Allreduce(MPI_IN_PLACE, center[channel_X].data(),
                  (center[channel_X]).rows() * (center[channel_X]).cols(),
                  MPI_DOUBLE, MPI_SUM, eq1_comm);
  }
  return;
}

void Eq1MVar::SubtractOperatorXvar_sub(Eq1Op op, SpEOMatrixD* m_var_sub,
                                       MPI_Comm eq1_comm,
                                       SpEOMatrixD& EndmemberMat,
                                       SpEOVectorD& LAMBDA_ZY) {
  int channel_X, channel_Y;

  //=========== SpEOMatrixD m_var = W*m_var_sub; ========>
  SpEOMatrixD* m_var = new SpEOMatrixD[EndmemberMat.rows()];
  bool transpose_2D_mat = false;

  multiply_2D_with_3D_mat(m_var, EndmemberMat, m_var_sub, transpose_2D_mat);

  SpEOMatrixD* TMP = new SpEOMatrixD[N_center];
  for (channel_X = 0; channel_X < N_center; channel_X++) {
    TMP[channel_X] = SpEOMatrixD::Zero((center[channel_X]).rows(),
                                       (center[channel_X]).cols());
  }
  for (channel_Y = 0; channel_Y < N_upper; channel_Y++) {
    upper[channel_Y] -= op.coeff1 * LAMBDA_ZY(channel_Y) * (m_var[channel_Y]);

    // calculate R*var
    for (channel_X = 0; channel_X < N_center; channel_X++) {
      TMP[channel_X] +=
          op.coeff2 * ((*op.R)(channel_X, channel_Y)) * (m_var[channel_Y]);
    }
    // calculate var*BS
    SpEOMatrixD tmp_image_HR = m_var[channel_Y];
    SpEOMatrixD tmp_image_LR;
    if (op.fastfilter == true) {
      fast_filter(tmp_image_LR, tmp_image_HR, op.gauss_filter, op.filter_coeff,
                  op.fDS);
    } else {
      standard_filter(tmp_image_LR, tmp_image_HR, op.gauss_filter,
                      op.filter_coeff, op.fDS);
    }
    lower[channel_Y] -= op.coeff3 * tmp_image_LR;
  }
  for (channel_X = 0; channel_X < N_center; channel_X++) {
    MPI_Allreduce(MPI_IN_PLACE, TMP[channel_X].data(),
                  (center[channel_X]).rows() * (center[channel_X]).cols(),
                  MPI_DOUBLE, MPI_SUM, eq1_comm);
    center[channel_X] -= TMP[channel_X];
  }

  delete[] m_var;
  return;
}

void Eq1MVar::SetOperatorXvar(Eq1Op op, Eq1SVar_sub* m_var_sub,
                              MPI_Comm eq1_comm, SpEOMatrixD& EndmemberMat,
                              SpEOVectorD& LAMBDA_ZY) {
  int channel_X, channel_Y;
  SpEOMatrixD* m_var = new SpEOMatrixD[N_upper];
  bool transpose_2D_mat = false;

  multiply_2D_with_3D_mat(m_var, EndmemberMat, m_var_sub->var,
                          transpose_2D_mat);

  for (channel_X = 0; channel_X < N_center; channel_X++) {
    center[channel_X] = SpEOMatrixD::Zero((center[channel_X]).rows(),
                                          (center[channel_X]).cols());
  }
  for (channel_Y = 0; channel_Y < N_upper; channel_Y++) {
    upper[channel_Y] = op.coeff1 * (m_var[channel_Y]) * LAMBDA_ZY(channel_Y);
    // calculate R*var
    for (channel_X = 0; channel_X < N_center; channel_X++) {
      center[channel_X] +=
          op.coeff2 * ((*op.R)(channel_X, channel_Y)) * (m_var[channel_Y]);
    }
    // calculate var*BS
    SpEOMatrixD tmp_image_HR = m_var[channel_Y];
    SpEOMatrixD tmp_image_LR;
    if (op.fastfilter == true) {
      fast_filter(tmp_image_LR, tmp_image_HR, op.gauss_filter, op.filter_coeff,
                  op.fDS);
    } else {
      standard_filter(tmp_image_LR, tmp_image_HR, op.gauss_filter,
                      op.filter_coeff, op.fDS);
    }
    lower[channel_Y] = op.coeff3 * tmp_image_LR;
  }
  for (channel_X = 0; channel_X < N_center; channel_X++) {
    MPI_Allreduce(MPI_IN_PLACE, center[channel_X].data(),
                  (center[channel_X]).rows() * (center[channel_X]).cols(),
                  MPI_DOUBLE, MPI_SUM, eq1_comm);
  }
  delete[] m_var;
  return;
}

//  class SVAR
Eq1SVar::Eq1SVar(int N_m) {
  var = new SpEOMatrixD[N_m];
  N = N_m;
  return;
}

Eq1SVar::~Eq1SVar() {
  delete[] var;
  return;
}

double Eq1SVar::NormSquared(MPI_Comm eq1_comm) {
  double norm_squared = 0.0;
  int channel_Y;
  for (channel_Y = 0; channel_Y < N; channel_Y++) {
    norm_squared += pow(var[channel_Y].norm(), 2);
  }
  MPI_Allreduce(MPI_IN_PLACE, &norm_squared, 1, MPI_DOUBLE, MPI_SUM, eq1_comm);
  return norm_squared;
}

void Eq1SVar::SetAdjointOperatorXvar(Eq1Op op, Eq1MVar* m_var) {
  int channel_X, channel_Y;
  for (channel_Y = 0; channel_Y < m_var->N_upper; channel_Y++) {
    var[channel_Y] = op.coeff1 * (m_var->upper[channel_Y]);
    for (channel_X = 0; channel_X < m_var->N_center; channel_X++) {
      var[channel_Y] += op.coeff2 * ((*op.R)(channel_X, channel_Y)) *
                        (m_var->center[channel_X]);
    }
    SpEOMatrixD tmp_image_LR = m_var->lower[channel_Y];
    SpEOMatrixD tmp_image_HR;
    standard_adjoint_filter(tmp_image_HR, tmp_image_LR, op.gauss_filter,
                            op.filter_coeff, op.fDS);
    var[channel_Y] += op.coeff3 * tmp_image_HR;
  }
  return;
}

// ================= Svar_sub =================================================
// >

Eq1SVar_sub::Eq1SVar_sub(int N_m) {
  var = new SpEOMatrixD[N_m];
  N = N_m;
  return;
}

Eq1SVar_sub::~Eq1SVar_sub() {
  delete[] var;
  return;
}

double Eq1SVar_sub::NormSquared(MPI_Comm eq1_comm) {
  double norm_squared = 0.0;
  int channel_Y;
  for (channel_Y = 0; channel_Y < N; channel_Y++) {
    norm_squared += pow(var[channel_Y].norm(), 2);
  }
  MPI_Allreduce(MPI_IN_PLACE, &norm_squared, 1, MPI_DOUBLE, MPI_SUM, eq1_comm);
  return norm_squared;
}

void Eq1SVar_sub::SetAdjointOperatorXvar(Eq1Op op, Eq1MVar* m_var,
                                         SpEOMatrixD& EndmemberMat,
                                         SpEOVectorD& LAMBDA_ZY) {
  int channel_X, channel_Y;
  SpEOMatrixD* var_tmp = new SpEOMatrixD[m_var->N_upper];
  for (channel_Y = 0; channel_Y < m_var->N_upper; channel_Y++) {
    var_tmp[channel_Y] =
        SpEOMatrixD::Zero(m_var->upper[0].rows(), m_var->upper[0].cols());
  }
  for (channel_Y = 0; channel_Y < m_var->N_upper; channel_Y++) {
    var_tmp[channel_Y] =
        op.coeff1 * (m_var->upper[channel_Y]) * LAMBDA_ZY(channel_Y);
    for (channel_X = 0; channel_X < m_var->N_center; channel_X++) {
      var_tmp[channel_Y] += op.coeff2 * ((*op.R)(channel_X, channel_Y)) *
                            (m_var->center[channel_X]);
    }
    SpEOMatrixD tmp_image_LR = m_var->lower[channel_Y];
    SpEOMatrixD tmp_image_HR;
    standard_adjoint_filter(tmp_image_HR, tmp_image_LR, op.gauss_filter,
                            op.filter_coeff, op.fDS);
    var_tmp[channel_Y] += op.coeff3 * tmp_image_HR;
  }
  bool transpose_2D_mat = true;
  multiply_2D_with_3D_mat(var, EndmemberMat, var_tmp, transpose_2D_mat);
  delete[] var_tmp;
  return;
}

//<<<============

void multiply_2D_with_3D_mat(SpEOMatrixD*& Mat_3D_output,
                             SpEOMatrixD& Mat_2D_input,
                             SpEOMatrixD*& Mat_3D_input,
                             bool transpose_2D_mat) {
  int channelU, channelV;
  SpEOMatrixD Mat_2D;
  if (transpose_2D_mat) {
    Mat_2D = Mat_2D_input.transpose();
  } else {
    Mat_2D = Mat_2D_input;
  }
  for (channelU = 0; channelU < Mat_2D.rows(); channelU++) {
    Mat_3D_output[channelU] =
        SpEOMatrixD::Zero(Mat_3D_input[0].rows(), Mat_3D_input[0].cols());
    for (channelV = 0; channelV < Mat_2D.cols(); channelV++) {
      Mat_3D_output[channelU] +=
          Mat_2D(channelU, channelV) * Mat_3D_input[channelV];
    }
  }
}

/**
 * Equation 6 Solver Class Functions
 **/
Eq6Op::Eq6Op(SpEOMatrixD* R_m, double coeff1_m, double coeff2_m) {
  R = R_m;
  coeff1 = coeff1_m;
  coeff2 = coeff2_m;
  return;
}

Eq6MVar::Eq6MVar() { return; }

double Eq6MVar::NormSquared() {
  return pow(lower.norm(), 2) + pow(upper.norm(), 2);
}

Eq6MVar::~Eq6MVar() { return; }

void Eq6MVar::SubtractOperatorXvar(Eq6Op op, SpEOVectorD* m_var) {
  upper -= op.coeff1 * (*m_var);
  // calculate R*var
  lower -= op.coeff2 * ((*op.R) * (*m_var));
  return;
}

void Eq6MVar::SetOperatorXvar(Eq6Op op, Eq6SVar* m_var) {
  upper = op.coeff1 * m_var->var;
  // calculate R*var
  lower = op.coeff2 * ((*op.R) * m_var->var);
  return;
}

Eq6SVar::Eq6SVar() { return; }

Eq6SVar::~Eq6SVar() { return; }

double Eq6SVar::NormSquared() { return pow(var.norm(), 2); }

void Eq6SVar::SetAdjointOperatorXvar(Eq6Op op, Eq6MVar* m_var) {
  var = op.coeff1 * m_var->upper +
        op.coeff2 * ((*op.R).transpose() * m_var->lower);
  return;
}

/**
 * Equation 7 Solver Class Functions
 **/
Eq7Op::Eq7Op(SpEOVectorD* dH_m, SpEOVectorD* dL_m, SpEOMatrixD* R_m,
             SpEOVectorD* P_lmd_vecs_m, SpEOMatrixI* P_lmd_idx_row_m,
             int** P_lmd_idx_bl_m, double* coeff1_m, double* coeff2_m,
             double coeff3_m) {
  dH = dH_m;
  dL = dL_m;
  R = R_m;
  P_lmd_vecs = P_lmd_vecs_m;
  P_lmd_idx_row = P_lmd_idx_row_m;
  P_lmd_idx_bl = P_lmd_idx_bl_m;
  coeff1 = coeff1_m;
  coeff2 = coeff2_m;
  coeff3 = coeff3_m;
  return;
}

Eq7MVar::Eq7MVar(int N_g_m) {
  upper = new SpEOMatrixD[N_g_m];
  center = new SpEOMatrixD[N_g_m];
  N_g = N_g_m;
  return;
}

double Eq7MVar::NormSquared() {
  double norm_squared = 0.0;
  int g;
  for (g = 0; g < N_g; g++) {
    norm_squared += pow(upper[g].norm(), 2);
    norm_squared += pow(center[g].norm(), 2);
  }
  norm_squared += pow(lower.norm(), 2);
  return norm_squared;
}

Eq7MVar::~Eq7MVar() {
  delete[] upper;
  delete[] center;
  return;
}

void Eq7MVar::SubtractOperatorXvar(Eq7Op op, SpEOVectorD* m_var) {
  int g, j, k;
  for (g = 0; g < N_g; g++) {
    upper[g] -= op.coeff1[g] * (op.dH[g] * m_var[g].transpose());
    center[g] -= op.coeff2[g] * (op.dL[g] * m_var[g].transpose());
  }
  int N_h = op.dH[0].size();
  int N_Y = op.R->cols();
  int sum_N_c = 0;
  for (g = 0; g < N_g; g++) {
    sum_N_c += op.P_lmd_vecs[g].size();
  }
  SpEOMatrixD tmpDA = SpEOMatrixD::Zero(sum_N_c, N_h);
  int running_idx = 0;
  for (g = 0; g < N_g; g++) {
    tmpDA.block(running_idx, 0, op.P_lmd_vecs[g].size(), N_h) =
        (op.dH[g] * m_var[g].transpose()).transpose();
    running_idx += op.P_lmd_vecs[g].size();
  }
  SpEOMatrixD tmpPDA = SpEOMatrixD::Zero(N_Y, N_h);
  for (j = 0; j < N_Y; j++) {
    for (k = 0; k < op.P_lmd_idx_row[j].cols(); k++) {
      int block = op.P_lmd_idx_row[j].coeff(0, k);
      int relidx = op.P_lmd_idx_row[j].coeff(1, k);
      tmpPDA.row(j) += op.P_lmd_vecs[block](relidx) *
                       tmpDA.row(op.P_lmd_idx_bl[block][1] + relidx);
    }
  }
  lower -= op.coeff3 * ((*op.R) * (tmpPDA));
  return;
}

void Eq7MVar::SetOperatorXvar(Eq7Op op, Eq7SVar* m_var) {
  int g, j, k;
  for (g = 0; g < N_g; g++) {
    upper[g] = op.coeff1[g] * (op.dH[g] * m_var->var[g].transpose());
    center[g] = op.coeff2[g] * (op.dL[g] * m_var->var[g].transpose());
  }
  int N_h = op.dH[0].size();
  int N_Y = op.R->cols();
  int sum_N_c = 0;
  for (g = 0; g < N_g; g++) {
    sum_N_c += op.P_lmd_vecs[g].size();
  }
  SpEOMatrixD tmpDA = SpEOMatrixD::Zero(sum_N_c, N_h);
  int running_idx = 0;
  for (g = 0; g < N_g; g++) {
    tmpDA.block(running_idx, 0, op.P_lmd_vecs[g].size(), N_h) =
        (op.dH[g] * m_var->var[g].transpose()).transpose();
    running_idx += op.P_lmd_vecs[g].size();
  }
  SpEOMatrixD tmpPDA = SpEOMatrixD::Zero(N_Y, N_h);
  for (j = 0; j < N_Y; j++) {
    for (k = 0; k < op.P_lmd_idx_row[j].cols(); k++) {
      int block = op.P_lmd_idx_row[j].coeff(0, k);
      int relidx = op.P_lmd_idx_row[j].coeff(1, k);
      tmpPDA.row(j) += op.P_lmd_vecs[block](relidx) *
                       tmpDA.row(op.P_lmd_idx_bl[block][1] + relidx);
    }
  }
  lower = op.coeff3 * ((*op.R) * (tmpPDA));
  return;
}

Eq7SVar::Eq7SVar(int N_g_m) {
  var = new SpEOVectorD[N_g_m];
  N_g = N_g_m;
  return;
}

Eq7SVar::~Eq7SVar() {
  delete[] var;
  return;
}

double Eq7SVar::NormSquared() {
  double norm_squared = 0.0;
  int g;
  for (g = 0; g < N_g; g++) {
    norm_squared += pow(var[g].norm(), 2);
  }
  return norm_squared;
}

void Eq7SVar::SetAdjointOperatorXvar(Eq7Op op, Eq7MVar* m_var) {
  int g, j;
  int N_h = op.dH[0].size();
  SpEOMatrixD tempR = (*op.R).transpose() * m_var->lower;
  for (g = 0; g < N_g; g++) {
    var[g] = (op.coeff1[g] * (op.dH[g].transpose() * m_var->upper[g]) +
              op.coeff2[g] * (op.dL[g].transpose() * m_var->center[g]))
                 .transpose();
    SpEOMatrixD temp = SpEOMatrixD::Zero(N_h, op.P_lmd_vecs[g].size());
    for (j = 0; j < op.P_lmd_vecs[g].size(); j++) {  // write temp = B^t*P(k)
      temp.col(j) = (op.P_lmd_vecs[g](j)) *
                    tempR.row(op.P_lmd_idx_bl[g][0] + j).transpose();
    }
    var[g] += op.coeff3 * (op.dH[g].transpose() * temp).transpose();
  }
  return;
}

/**
 * Equation 10 Solver Class Functions
 **/
Eq9Op::Eq9Op(SpEOMatrixD* DH_m, SpEOMatrixD* DL_m, SpEOMatrixD* R_m,
             SpEOVectorD* P_lmd_vecs_m, SpEOMatrixI* P_lmd_idx_row_m,
             int** P_lmd_idx_bl_m, double* coeff1_m, double* coeff2_m,
             double coeff3_m) {
  DH = DH_m;
  DL = DL_m;
  R = R_m;
  P_lmd_vecs = P_lmd_vecs_m;
  P_lmd_idx_row = P_lmd_idx_row_m;
  P_lmd_idx_bl = P_lmd_idx_bl_m;
  coeff1 = coeff1_m;
  coeff2 = coeff2_m;
  coeff3 = coeff3_m;
  return;
}

Eq9MVar::Eq9MVar(int N_g_m) {
  upper = new SpEOMatrixD[N_g_m];
  center = new SpEOMatrixD[N_g_m];
  N_g = N_g_m;
  return;
}

double Eq9MVar::NormSquared() {
  double norm_squared = 0.0;
  int g;
  for (g = 0; g < N_g; g++) {
    norm_squared += pow(upper[g].norm(), 2);
    norm_squared += pow(center[g].norm(), 2);
  }
  norm_squared += pow(lower.norm(), 2);
  return norm_squared;
}

Eq9MVar::~Eq9MVar() {
  delete[] upper;
  delete[] center;
  return;
}

void Eq9MVar::SubtractOperatorXvar(Eq9Op op, SpEOMatrixD* m_var) {
  int g, j, k;
  for (g = 0; g < N_g; g++) {
    int N_DP = op.DH[g].cols();
    upper[g] -= op.coeff1[g] *
                (op.DH[g].rightCols(N_DP - 1) * m_var[g].bottomRows(N_DP - 1));
    center[g] -= op.coeff2[g] *
                 (op.DL[g].rightCols(N_DP - 1) * m_var[g].bottomRows(N_DP - 1));
  }
  int N_h = op.DH[0].rows();
  int N_Y = op.R->cols();
  int sum_N_c = 0;
  for (g = 0; g < N_g; g++) {
    sum_N_c += op.P_lmd_vecs[g].size();
  }
  SpEOMatrixD tmpDA = SpEOMatrixD::Zero(sum_N_c, N_h);
  int running_idx = 0;
  for (g = 0; g < N_g; g++) {
    int N_DP = op.DH[g].cols();
    tmpDA.block(running_idx, 0, op.P_lmd_vecs[g].size(), N_h) =
        (op.DH[g].rightCols(N_DP - 1) * m_var[g].bottomRows(N_DP - 1))
            .transpose();
    running_idx += op.P_lmd_vecs[g].size();
  }
  SpEOMatrixD tmpPDA = SpEOMatrixD::Zero(N_Y, N_h);
  for (j = 0; j < N_Y; j++) {
    for (k = 0; k < op.P_lmd_idx_row[j].cols(); k++) {
      int block = op.P_lmd_idx_row[j].coeff(0, k);
      int relidx = op.P_lmd_idx_row[j].coeff(1, k);
      tmpPDA.row(j) += op.P_lmd_vecs[block](relidx) *
                       tmpDA.row(op.P_lmd_idx_bl[block][1] + relidx);
    }
  }
  lower -= op.coeff3 * ((*op.R) * (tmpPDA));
  return;
}

void Eq9MVar::SetOperatorXvar(Eq9Op op, Eq9SVar* m_var) {
  int g, j, k;
  for (g = 0; g < N_g; g++) {
    int N_DP = op.DH[g].cols();
    upper[g] = op.coeff1[g] * (op.DH[g].rightCols(N_DP - 1) * m_var->var[g]);
    center[g] = op.coeff2[g] * (op.DL[g].rightCols(N_DP - 1) * m_var->var[g]);
  }
  int N_h = op.DH[0].rows();
  int N_Y = op.R->cols();
  int sum_N_c = 0;
  for (g = 0; g < N_g; g++) {
    sum_N_c += op.P_lmd_vecs[g].size();
  }
  SpEOMatrixD tmpDA = SpEOMatrixD::Zero(sum_N_c, N_h);
  int running_idx = 0;
  for (g = 0; g < N_g; g++) {
    int N_DP = op.DH[g].cols();
    tmpDA.block(running_idx, 0, op.P_lmd_vecs[g].size(), N_h) =
        (op.DH[g].rightCols(N_DP - 1) * m_var->var[g]).transpose();
    running_idx += op.P_lmd_vecs[g].size();
  }
  SpEOMatrixD tmpPDA = SpEOMatrixD::Zero(N_Y, N_h);
  for (j = 0; j < N_Y; j++) {
    for (k = 0; k < op.P_lmd_idx_row[j].cols(); k++) {
      int block = op.P_lmd_idx_row[j].coeff(0, k);
      int relidx = op.P_lmd_idx_row[j].coeff(1, k);
      tmpPDA.row(j) += op.P_lmd_vecs[block](relidx) *
                       tmpDA.row(op.P_lmd_idx_bl[block][1] + relidx);
    }
  }
  lower = op.coeff3 * ((*op.R) * (tmpPDA));
  return;
}

Eq9SVar::Eq9SVar(int N_g_m) {
  var = new SpEOMatrixD[N_g_m];
  N_g = N_g_m;
  return;
}

Eq9SVar::~Eq9SVar() {
  delete[] var;
  return;
}

double Eq9SVar::NormSquared() {
  double norm_squared = 0.0;
  int g;
  for (g = 0; g < N_g; g++) {
    norm_squared += pow(var[g].norm(), 2);
  }
  return norm_squared;
}

void Eq9SVar::SetAdjointOperatorXvar(Eq9Op op, Eq9MVar* m_var) {
  int g, j;
  int N_h = op.DH[0].rows();
  SpEOMatrixD tempR = (*op.R).transpose() * m_var->lower;
  for (g = 0; g < N_g; g++) {
    int N_DP = op.DH[g].cols();
    var[g] = op.coeff1[g] *
                 (op.DH[g].rightCols(N_DP - 1).transpose() * m_var->upper[g]) +
             op.coeff2[g] *
                 (op.DL[g].rightCols(N_DP - 1).transpose() * m_var->center[g]);
    SpEOMatrixD temp = SpEOMatrixD::Zero(N_h, op.P_lmd_vecs[g].size());
    for (j = 0; j < op.P_lmd_vecs[g].size(); j++) {  // write temp = B^t*P(k)
      temp.col(j) = (op.P_lmd_vecs[g](j)) *
                    tempR.row(op.P_lmd_idx_bl[g][0] + j).transpose();
    }
    var[g] += op.coeff3 * (op.DH[g].rightCols(N_DP - 1).transpose() * temp);
  }
  return;
}

/**
 * FBS Solver
 **/
void softThresholdJS(double alpha, SpEOMatrixD& x) {
  int i;
  for (i = 0; i < x.rows(); i++) {
    double x_abs = x.row(i).norm();
    x.row(i) = max(0.0, 1.0 - alpha / x_abs) * x.row(i);
  }
  return;
}

void softThresholdJS(double alpha, SpEOMatrixD* A, SpEOMatrixD& x,
                     SpEOMatrixD& Ax) {
  int i;
  for (i = 0; i < x.rows(); i++) {
    double x_abs = x.row(i).norm();
    if (x_abs > alpha) {
      x.row(i) = (1.0 - alpha / x_abs) * x.row(i);
    } else {
      x.row(i) = SpEOMatrixD::Zero(1, x.cols());
    }
  }
  Ax = (*A) * x;
  return;
}

void softThresholdJS_Phit(double alpha, SpEOMatrixD* AtA, SpEOMatrixD& x,
                          SpEOMatrixD& AtAx) {
  int i;
  for (i = 0; i < x.rows(); i++) {
    double x_abs = x.row(i).norm();
    if (x_abs > alpha) {
      x.row(i) = (1.0 - alpha / x_abs) * x.row(i);
    } else {
      x.row(i) = SpEOMatrixD::Zero(1, x.cols());
    }
  }
  AtAx = (*AtA) * x;
  return;
}
double J_1D(double x, SpEOMatrixD Ad_h, SpEOMatrixD y_Ab, SpEOMatrixD d_h,
            SpEOMatrixD b_h, double tau) {
  return 0.5 * pow((x * Ad_h - y_Ab).norm(), 2) + tau * l1l2norm(x * d_h + b_h);
}

double J_1D(double x, double Ad_h_sq, double Ad_h_y_Ab, double y_Ab_sq,
            SpEOVectorD d_s, SpEOVectorD d_s_b_s, SpEOVectorD b_s, double tau) {
  return 0.5 * (Ad_h_sq * pow(x, 2) - Ad_h_y_Ab * x + y_Ab_sq) +
         tau *
             (pow(x, 2) * d_s + x * d_s_b_s + b_s).cwiseAbs().cwiseSqrt().sum();
}

double J_1D(double x, double Ad_h_sq, double Ad_h_y_Ab, double y_Ab_sq,
            SpEOVectorD d_s, SpEOVectorD d_s_b_s, SpEOVectorD b_s, double tau,
            MPI_Comm comm_group) {
  double commdata =
      (pow(x, 2) * d_s + x * d_s_b_s + b_s).cwiseAbs().cwiseSqrt().sum();
  MPI_Allreduce(MPI_IN_PLACE, &commdata, 1, MPI_DOUBLE, MPI_SUM, comm_group);
  return 0.5 * (Ad_h_sq * pow(x, 2) - Ad_h_y_Ab * x + y_Ab_sq) + tau * commdata;
}

double J(SpEOMatrixD* A, SpEOMatrixD* y, double mu, SpEOMatrixD* X) {
  return 0.5 * pow(((*A) * (*X) - (*y)).norm(), 2) + mu * l1l2norm(*X);
}

double J_AX(SpEOMatrixD* AX, SpEOMatrixD* y, double mu, SpEOMatrixD* X) {
  return 0.5 * pow((*AX - *y).norm(), 2) + mu * l1l2norm(*X);
}

double J_minsearch(SpEOMatrixD Ad_h, SpEOMatrixD y_Ab, SpEOMatrixD d_h,
                   SpEOMatrixD b_h, double tau, double x0,
                   double tol) {  // Nelder Mead Simplex in 1D
  double alpha = 1;
  double gamma = 2;
  double beta = 0.5;
  double sigma = 0.5;
  double x1 = x0 + 1.0;
  double f0 = J_1D(x0, Ad_h, y_Ab, d_h, b_h, tau);
  double f1 = J_1D(x1, Ad_h, y_Ab, d_h, b_h, tau);
  while (abs(x1 - x0) > tol) {
    if (f0 > f1) {  // sort
      double xtmp = x0;
      x0 = x1;
      x1 = xtmp;
      double ftmp = f0;
      f0 = f1;
      f1 = ftmp;
    }
    double r = (1 + alpha) * x0 - alpha * x1;
    double fr = J_1D(r, Ad_h, y_Ab, d_h, b_h, tau);
    if (fr < f0) {
      double e = (1 + gamma) * x0 - gamma * x1;
      double fe = J_1D(e, Ad_h, y_Ab, d_h, b_h, tau);
      if (fr < fe) {
        x1 = r;
        f1 = fr;
      } else {
        x1 = e;
        f1 = fe;
      }
    } else {
      double c;
      if (fr < f1) {
        c = beta * x0 + (1 - beta) * r;
      } else {
        c = beta * x0 + (1 - beta) * x1;
      }
      double fc = J_1D(c, Ad_h, y_Ab, d_h, b_h, tau);
      if (fc < f1) {
        x1 = c;
        f1 = fc;
      } else {
        x1 = sigma * x0 + (1 - sigma) * x1;
        f1 = J_1D(x1, Ad_h, y_Ab, d_h, b_h, tau);
      }
    }
  }
  return x0;
}

double J_minsearch(double Ad_h_sq, double Ad_h_y_Ab, double y_Ab_sq,
                   SpEOVectorD d_s, SpEOVectorD d_s_b_s, SpEOVectorD b_s,
                   double tau, double x0, double tol,
                   int& iter) {  // Nelder Mead Simplex in 1D
  double alpha = 1.0;
  double gamma = 2.0;
  double beta = 0.5;
  double sigma = 0.5;
  double x1 = x0 + 1.0;
  double f0 = J_1D(x0, Ad_h_sq, Ad_h_y_Ab, y_Ab_sq, d_s, d_s_b_s, b_s, tau);
  double f1 = J_1D(x1, Ad_h_sq, Ad_h_y_Ab, y_Ab_sq, d_s, d_s_b_s, b_s, tau);
  iter = 0;
  while (abs(x1 - x0) > tol) {
    iter++;
    if (f0 > f1) {  // sort
      double xtmp = x0;
      x0 = x1;
      x1 = xtmp;
      double ftmp = f0;
      f0 = f1;
      f1 = ftmp;
    }
    double r = (1.0 + alpha) * x0 - alpha * x1;
    double fr = J_1D(r, Ad_h_sq, Ad_h_y_Ab, y_Ab_sq, d_s, d_s_b_s, b_s, tau);
    if (fr < f0) {
      double e = (1.0 + gamma) * x0 - gamma * x1;
      double fe = J_1D(e, Ad_h_sq, Ad_h_y_Ab, y_Ab_sq, d_s, d_s_b_s, b_s, tau);
      if (fr < fe) {
        x1 = r;
        f1 = fr;
      } else {
        x1 = e;
        f1 = fe;
      }
    } else {
      double c;
      if (fr < f1) {
        c = beta * x0 + (1.0 - beta) * r;
      } else {
        c = beta * x0 + (1.0 - beta) * x1;
      }
      double fc = J_1D(c, Ad_h_sq, Ad_h_y_Ab, y_Ab_sq, d_s, d_s_b_s, b_s, tau);
      if (fc < f1) {
        x1 = c;
        f1 = fc;
      } else {
        x1 = sigma * x0 + (1.0 - sigma) * x1;
        f1 = J_1D(x1, Ad_h_sq, Ad_h_y_Ab, y_Ab_sq, d_s, d_s_b_s, b_s, tau);
      }
    }
  }
  return x0;
}

double J_minsearch(double Ad_h_sq, double Ad_h_y_Ab, double y_Ab_sq,
                   SpEOVectorD d_s, SpEOVectorD d_s_b_s, SpEOVectorD b_s,
                   double tau, double x0, double tol, int& iter,
                   MPI_Comm comm_group) {  // Nelder Mead Simplex in 1D
  double alpha = 1.0;
  double gamma = 2.0;
  double beta = 0.5;
  double sigma = 0.5;
  double x1 = x0 + 1.0;
  double f0 =
      J_1D(x0, Ad_h_sq, Ad_h_y_Ab, y_Ab_sq, d_s, d_s_b_s, b_s, tau, comm_group);
  double f1 =
      J_1D(x1, Ad_h_sq, Ad_h_y_Ab, y_Ab_sq, d_s, d_s_b_s, b_s, tau, comm_group);
  iter = 0;
  while (abs(x1 - x0) > tol) {
    iter++;
    if (f0 > f1) {  // sort
      double xtmp = x0;
      x0 = x1;
      x1 = xtmp;
      double ftmp = f0;
      f0 = f1;
      f1 = ftmp;
    }
    double r = (1.0 + alpha) * x0 - alpha * x1;
    double fr = J_1D(r, Ad_h_sq, Ad_h_y_Ab, y_Ab_sq, d_s, d_s_b_s, b_s, tau,
                     comm_group);
    if (fr < f0) {
      double e = (1.0 + gamma) * x0 - gamma * x1;
      double fe = J_1D(e, Ad_h_sq, Ad_h_y_Ab, y_Ab_sq, d_s, d_s_b_s, b_s, tau,
                       comm_group);
      if (fr < fe) {
        x1 = r;
        f1 = fr;
      } else {
        x1 = e;
        f1 = fe;
      }
    } else {
      double c;
      if (fr < f1) {
        c = beta * x0 + (1.0 - beta) * r;
      } else {
        c = beta * x0 + (1.0 - beta) * x1;
      }
      double fc = J_1D(c, Ad_h_sq, Ad_h_y_Ab, y_Ab_sq, d_s, d_s_b_s, b_s, tau,
                       comm_group);
      if (fc < f1) {
        x1 = c;
        f1 = fc;
      } else {
        x1 = sigma * x0 + (1.0 - sigma) * x1;
        f1 = J_1D(x1, Ad_h_sq, Ad_h_y_Ab, y_Ab_sq, d_s, d_s_b_s, b_s, tau,
                  comm_group);
      }
    }
  }
  return x0;
}

// parallel FBS
void FBSSolver(int& iter_need, double& rel_res, SpEOMatrixD& X, int N,
               SpEOMatrixD Phi, SpEOMatrixD y, double alpha,
               FBSsolverOptions opts, MPI_Comm comm_group,
               SpEOMatrixD& timestat, SpEOMatrixD& btstat) {
  double eta = 0.5;

  int i;

  int my_rank;
  int my_processes;
  MPI_Comm_rank(comm_group, &my_rank);
  MPI_Comm_size(comm_group, &my_processes);

  timestat = SpEOMatrixD::Zero(
      1,
      5);  // total, communication, backtracking, prediction, stopping criterion
  double time_tmp = MPI_Wtime();

  /* #################
   * global parameters
   * ################# */
  int L_max = opts.maxiter_in;
  if (opts.pred_rule == FISTA && L_max > 1) {
    if (my_rank == 0 && !opts.silent) {
      cout << "FBS Solver WARNING: Prediction Rule FISTA is only allowed with "
              "opts.maxiter_in == 1. Value is set to 1!"
           << endl;
    }
    L_max = 1;
  } else if (opts.stepsize_update == PSCL && L_max > 1) {
    if (my_rank == 0 && !opts.silent) {
      cout << "FBS Solver WARNING: Stepsize Update Rule PSCL is only allowed "
              "with opts.maxiter_in == 1. Value is set to 1!"
           << endl;
    }
    L_max = 1;
  }
  int P = opts.Ndecomp;
  if (P != my_processes) {
    if (my_rank == 0 && !opts.silent) {
      cout << "FBS Solver WARNING: P = " << P
           << " is not equal to number of processes. P is set to "
           << my_processes << "!" << endl;
    }
  }
  double stepsize_global = 1.0;
  int N_loc = Phi.cols();
  int m = Phi.rows();
  int d = y.cols();
  // only FISTA
  double t = 1.0;
  // only adaptive stepsize
  double stepsize_new = stepsize_global;

  SpEOMatrixD Phit = Phi.transpose();
  SpEOMatrixD Phity = Phit * y;
  SpEOMatrixD PhitPhi = Phit * Phi;
  SpEOMatrixD Phity_loc = Phity;
  SpEOMatrixD y_loc = y;

  SpEOMatrixD X_loc = SpEOMatrixD::Zero(N_loc, d);
  SpEOMatrixD X_loc_prev = SpEOMatrixD::Zero(N_loc, d);
  SpEOMatrixD Xtmp = SpEOMatrixD::Zero(N_loc, d);
  SpEOMatrixD Xtmp_prev = SpEOMatrixD::Zero(N_loc, d);
  SpEOMatrixD Xin = SpEOMatrixD::Zero(N_loc, d);

  SpEOMatrixD PhiX = SpEOMatrixD::Zero(m, d);
  SpEOMatrixD PhiX_prev = SpEOMatrixD::Zero(m, d);
  SpEOMatrixD PhiX_loc = SpEOMatrixD::Zero(m, d);
  SpEOMatrixD PhiX_loc_prev = SpEOMatrixD::Zero(m, d);
  SpEOMatrixD PhiXtmp = SpEOMatrixD::Zero(m, d);
  SpEOMatrixD PhiXtmp_loc = SpEOMatrixD::Zero(m, d);

  SpEOMatrixD PhitPhiX = SpEOMatrixD::Zero(N_loc, d);
  SpEOMatrixD PhitPhiX_prev = SpEOMatrixD::Zero(N_loc, d);
  SpEOMatrixD PhitPhiX_loc = SpEOMatrixD::Zero(N_loc, d);
  SpEOMatrixD PhitPhiX_loc_prev = SpEOMatrixD::Zero(N_loc, d);
  SpEOMatrixD PhitPhiXtmp = SpEOMatrixD::Zero(N_loc, d);
  SpEOMatrixD PhitPhiXtmp_prev = SpEOMatrixD::Zero(N_loc, d);
  SpEOMatrixD PhitPhiXtmp_loc = SpEOMatrixD::Zero(N_loc, d);
  SpEOMatrixD PhitPhiXtmp_loc_prev = SpEOMatrixD::Zero(N_loc, d);
  SpEOMatrixD PhitPhiXin = SpEOMatrixD::Zero(N_loc, d);
  SpEOMatrixD PhitPhiXin_res = SpEOMatrixD::Zero(N_loc, d);

  SpEOMatrixD Phity_PhitPhiX = SpEOMatrixD::Zero(N_loc, d);
  SpEOMatrixD Phity_PhitPhiX_prev = SpEOMatrixD::Zero(N_loc, d);

  // normalizer for stopping criterion
  double Phitynorm = pow(Phity.norm(), 2);
  MPI_Allreduce(MPI_IN_PLACE, &Phitynorm, 1, MPI_DOUBLE, MPI_SUM, comm_group);
  Phitynorm = sqrt(Phitynorm);

  double y_norm_sq = pow(y.norm(), 2);

  double Jnew = 0.5 * y_norm_sq;
  double Jold = 0.0;
  double weight_predict = 1.0;
  double f_hat = y_norm_sq;

  int iterations = 0;
  iter_need = 0;
  bool stop_crit = false;
  int L = 1;
  while (!stop_crit && (iterations < opts.maxiter_out)) {
    // update of X_loc
    X_loc_prev = X_loc;
    PhiX_loc_prev = PhiX_loc;
    PhitPhiX_loc_prev = PhitPhiX_loc;
    PhiX_prev = PhiX;

    Xin = Xtmp;
    PhitPhiXin = PhitPhiXtmp_loc;

    y_loc = y - PhiXtmp + PhiXtmp_loc;
    Phity_loc = Phity - PhitPhiXtmp + PhitPhiXtmp_loc;

    double y_loc_norm_sq_half = 0.5 * pow(y_loc.norm(), 2);

    timestat(0, 0) += MPI_Wtime() - time_tmp;
    time_tmp = MPI_Wtime();

    // update L
    if (L_max == 1 || iterations == 0 || opts.gamma * opts.tol < rel_res) {
      L = 1;
    } else {
      L = max(L,
              min(L_max, (int)(2 + log10((rel_res / opts.tol) / opts.gamma) *
                                       (L_max - 1.0) / (-log10(opts.gamma)))));
    }
    double stepsize;
    int iter_in;
    for (iter_in = 0; iter_in < L; iter_in++) {
      PhitPhiXin_res = Phity_loc - PhitPhiXin;

      // update stepsize
      if (opts.stepsize_update == INIT) {
        stepsize = stepsize_new;
      } else if (opts.stepsize_update == LARGE) {
        stepsize = 0.99 * pow(eta, -10);
      } else if (opts.stepsize_update == PSCL) {
        if (iterations == 0) {
          stepsize = min(1 + 1.665 * (1 - m / N), 1.999);
        } else {
          double tmp_sz1 = pow((Xtmp - Xtmp_prev).norm(), 2);
          double tmp_sz2 = 2 * (Xtmp - Xtmp_prev)
                                   .cwiseProduct(PhitPhiXtmp - PhitPhiXtmp_prev)
                                   .sum();
          stepsize = max(1.7 * tmp_sz1 / tmp_sz2, 1.0);
        }
      }

      /*
       * UPDATE Xnew and BACKTRACKING
       */
      if (opts.backtracking == DEC || opts.backtracking == DECEFF) {
        double f_new;
        bool step_cond = false;
        stepsize /= eta;
        while (!step_cond && stepsize > 0.01) {
          stepsize *= eta;
          X_loc = Xin + stepsize * (PhitPhiXin_res);
          if (opts.backtracking == DEC) {
            softThresholdJS(stepsize * alpha, &Phi, X_loc, PhiX_loc);
          } else {
            softThresholdJS_Phit(stepsize * alpha, &PhitPhi, X_loc,
                                 PhitPhiX_loc);
          }
          if (opts.backtracking == DEC) {
            f_new = 0.5 * pow((y_loc - PhiX_loc).norm(), 2);
          } else {
            f_new = y_loc_norm_sq_half +
                    X_loc.cwiseProduct(0.5 * PhitPhiX_loc - Phity_loc).sum();
          }
          step_cond =
              f_new < 0.5 * f_hat +
                          (X_loc - Xin)
                              .cwiseProduct(-PhitPhiXin_res +
                                            1 / (2 * stepsize) * (X_loc - Xin))
                              .sum();
        }
        f_hat = f_new;
        if (stepsize <= 0.01 && !step_cond && !opts.silent) {
          cout << "FBS Solver WARNING: backtracking stepsize limit reached in "
                  "iteration "
               << iterations << " on process with rank " << my_rank
               << ". Condition violation by "
               << 0.5 * pow((y - PhiX_loc).norm(), 2) - 0.5 * f_hat +
                      (X_loc - Xin).cwiseProduct(PhitPhiXin_res).sum() -
                      1 / (2 * stepsize) * pow((X_loc - Xin).norm(), 2)
               << "." << endl;
        }
      } else if (opts.backtracking == DECPAR) {
        bool step_cond = false;
        stepsize /= eta;
        while (!step_cond && stepsize > 0.01) {
          stepsize *= eta;
          X_loc = Xin + stepsize * (PhitPhiXin_res);
          softThresholdJS(stepsize * alpha, &Phi, X_loc, PhiX_loc);
          double tmp_stepcond =
              -(X_loc - Xin).cwiseProduct(PhitPhiXin_res).sum() +
              1 / (2 * stepsize) * pow((X_loc - Xin).norm(), 2);
          timestat(0, 2) += MPI_Wtime() - time_tmp;
          time_tmp = MPI_Wtime();
          MPI_Allreduce(PhiX_loc.data(), PhiX.data(), m * d, MPI_DOUBLE,
                        MPI_SUM, comm_group);
          MPI_Allreduce(MPI_IN_PLACE, &tmp_stepcond, 1, MPI_DOUBLE, MPI_SUM,
                        comm_group);
          timestat(0, 1) += MPI_Wtime() - time_tmp;
          time_tmp = MPI_Wtime();
          step_cond =
              0.5 * pow((y - PhiX).norm(), 2) < 0.5 * f_hat + tmp_stepcond;
        }
        if (stepsize <= 0.01 && !step_cond && !opts.silent) {
          cout << "FBS Solver WARNING: backtracking stepsize limit reached in "
                  "iteration "
               << iterations << " on process with rank " << my_rank << ". "
               << endl;
        }
      } else if (opts.backtracking == NO) {
        // Bregman step
        X_loc = Xin + stepsize * (PhitPhiXin_res);
        softThresholdJS_Phit(stepsize * alpha, &PhitPhi, X_loc, PhitPhiX_loc);
      }

      if (opts.backtracking != DECEFF && opts.backtracking != NO) {
        PhitPhiX_loc = Phit * (PhiX_loc);
      }
      if (iter_in < L - 1) {
        Xin = X_loc;
        PhitPhiXin = PhitPhiX_loc;
      }
    }

    /*
     * SHARE RESULTS -> COMMUNICATION
     */
    if (opts.backtracking == DECEFF || opts.backtracking == NO) {
      PhiX_loc = Phi * X_loc;
    }

    timestat(0, 2) += MPI_Wtime() - time_tmp;
    time_tmp = MPI_Wtime();
    if (opts.backtracking != DECPAR) {
      MPI_Allreduce(PhiX_loc.data(), PhiX.data(), m * d, MPI_DOUBLE, MPI_SUM,
                    comm_group);
    }
    timestat(0, 1) += MPI_Wtime() - time_tmp;
    time_tmp = MPI_Wtime();

    // PhitPhiXall = Phit*(PhiXall);
    PhitPhiX_prev = PhitPhiX;
    PhitPhiX = Phit * (PhiX);
    Phity_PhitPhiX_prev = Phity_PhitPhiX;
    Phity_PhitPhiX = Phity - PhitPhiX;

    timestat(0, 2) += MPI_Wtime() - time_tmp;
    time_tmp = MPI_Wtime();

    /*
     * CHECK STOPPING CRITERION
     */
    switch (opts.stop_crit) {
      case FUNCTIONAL_DIFFERENCE: {
        // update J
        Jold = Jnew;
        double buf_comm = 0.5 * X_loc.cwiseProduct(PhitPhiX - 2 * Phity).sum() +
                          alpha * l1l2norm(X_loc);
        MPI_Allreduce(MPI_IN_PLACE, &buf_comm, 1, MPI_DOUBLE, MPI_SUM,
                      comm_group);
        Jnew = 0.5 * pow(y.norm(), 2) + buf_comm;
        // Jnew =0.5*pow((y-PhiXall).norm(),2) + alpha * buf_l1l2norm;
        rel_res = abs(Jnew - Jold) / abs(Jold);
        stop_crit = rel_res < opts.tol;
        break;
      }
      case FIRST_ORDER_CONDITIONS: {
        // update stop_crit preparation
        double norm_squared = 0.0;
        for (i = 0; i < N_loc; i++) {
          double normXnew = X_loc.row(i).norm();
          if (normXnew > 1e-16) {
            norm_squared +=
                pow((-Phity_PhitPhiX.row(i) + (alpha / normXnew) * X_loc.row(i))
                        .norm(),
                    2);
          } else {
            norm_squared +=
                pow(max(0.0, Phity_PhitPhiX.row(i).norm() - alpha), 2);
          }
        }

        // update stepsize preparation
        SpEOMatrixD Dgrad_tmp = -Phity_PhitPhiX + Phity_PhitPhiX_prev;
        SpEOMatrixD DX_tmp = X_loc - X_loc_prev;

        double* buf = new double[4];
        double* buf_recv = new double[4];
        buf[0] = norm_squared;
        buf[1] = pow(DX_tmp.norm(), 2);
        buf[2] = DX_tmp.cwiseProduct(Dgrad_tmp).sum();
        buf[3] = pow(Dgrad_tmp.norm(), 2);

        // communication
        MPI_Allreduce(buf, buf_recv, 4, MPI_DOUBLE, MPI_SUM, comm_group);

        // update stop_crit final evaluation
        // rel_res = sqrt(buf_recv[0]*N*m)/y.norm();
        rel_res = sqrt(buf_recv[0]) / Phitynorm;
        stop_crit = rel_res < opts.tol;

        if (opts.stepsize_adaptive == true) {
          // update stepsize final evaluation
          double tau_s = buf_recv[1] / buf_recv[2];
          double tau_m = buf_recv[2] / buf_recv[3];
          if (2 * tau_m > tau_s) {
            stepsize_new = tau_m;
          } else {
            stepsize_new = tau_s - 0.5 * tau_m;
          }
          if (stepsize_new <= 0) {
            stepsize_new = stepsize_global;
          }
        }

        delete[] buf;
        delete[] buf_recv;
        break;
      }
      default: {
        cout << "FBS Solver ERROR: Stopping Rule not known";
        exit(EXIT_FAILURE);
        break;
      }
    }

    timestat(0, 4) += MPI_Wtime() - time_tmp;
    time_tmp = MPI_Wtime();
    /*
     * PREDICTION
     */
    switch (opts.pred_rule) {
      case FISTA:  // guaranteed to converge only for param.maxiter_in == 1
      case FISTA_RESTART: {
        double t_prev = t;
        if (opts.pred_rule == FISTA_RESTART) {
          double restartval =
              ((Xtmp - X_loc).cwiseProduct(X_loc - X_loc_prev)).sum();
          MPI_Allreduce(MPI_IN_PLACE, &restartval, 1, MPI_DOUBLE, MPI_SUM,
                        comm_group);
          if (restartval >= 0) {  // restart
            if (my_rank == 0 && !opts.silent) {
              cout << "FBS Solver HINT: FISTA was restarted in iteration "
                   << iterations << endl;
            }
            t_prev = 1;
          }
        }
        t = (1 + sqrt(1 + 4 * pow(t_prev, 2))) / 2;
        weight_predict = (1 + (t_prev - 1) / t);
        break;
      }
      case STANDARD: {
        weight_predict = 1.0 / P;  // standard FBS guaranteed to converge
        break;
      }
      case STANDARD_AGGRESSIVE: {
        weight_predict = 1.0;  // guaranteed to converge only for P == 1 or
                               // param.maxiter_in == 1
        break;
      }
      case LINE_SEARCH_SIMPLEX: {
        weight_predict = 1.0 / P;
        double weight_eps = 1e-2;
        SpEOMatrixD d_h = X_loc - Xtmp;
        SpEOMatrixD b_h = Xtmp;

        SpEOMatrixD PhitPhiX_PhitPhiXtmp = PhitPhiX - PhitPhiXtmp;
        SpEOMatrixD Phity_PhitPhiXtmp = Phity - PhitPhiXtmp;

        double Ad_h_sq = d_h.cwiseProduct(PhitPhiX_PhitPhiXtmp).sum();
        double y_Ab_sq = 0.0;
        double Ad_h_y_Ab = 2 * d_h.cwiseProduct(Phity_PhitPhiXtmp).sum();

        MPI_Allreduce(MPI_IN_PLACE, &Ad_h_sq, 1, MPI_DOUBLE, MPI_SUM,
                      comm_group);
        MPI_Allreduce(MPI_IN_PLACE, &Ad_h_y_Ab, 1, MPI_DOUBLE, MPI_SUM,
                      comm_group);

        // parallel:
        SpEOVectorD d_s = d_h.cwiseAbs2().rowwise().sum();
        SpEOVectorD b_s = b_h.cwiseAbs2().rowwise().sum();
        SpEOVectorD d_s_b_s = 2 * d_h.cwiseProduct(b_h).rowwise().sum();
        int iter_loc;

        double w = J_minsearch(Ad_h_sq, Ad_h_y_Ab, y_Ab_sq, d_s, d_s_b_s, b_s,
                               alpha, 0, weight_eps, iter_loc, comm_group);

        if (w > 0) {
          weight_predict = w;
        }
        break;
      }
      case LINE_SEARCH_ARMIJO: {
        weight_predict = 1.0 / P;
        double rho = 0.5;
        double sigma = 0.5;
        double eta_armijo = 2.0;
        SpEOMatrixD d_h = X_loc - Xtmp;
        SpEOMatrixD b_h = Xtmp;

        SpEOMatrixD PhitPhiX_PhitPhiXtmp = PhitPhiX - PhitPhiXtmp;
        SpEOMatrixD Phity_PhitPhiXtmp = Phity - PhitPhiXtmp;

        double Ad_h_sq = d_h.cwiseProduct(PhitPhiX_PhitPhiXtmp).sum();
        double y_Ab_sq = 0.0;
        double Ad_h_y_Ab = 2 * d_h.cwiseProduct(Phity_PhitPhiXtmp).sum();

        MPI_Allreduce(MPI_IN_PLACE, &Ad_h_sq, 1, MPI_DOUBLE, MPI_SUM,
                      comm_group);
        MPI_Allreduce(MPI_IN_PLACE, &Ad_h_y_Ab, 1, MPI_DOUBLE, MPI_SUM,
                      comm_group);

        double l1l2b_h = l1l2norm(b_h);
        double allred_result =
            2 * alpha * (l1l2norm(b_h + 1.0 / P * d_h) - l1l2b_h);
        MPI_Allreduce(MPI_IN_PLACE, &allred_result, 1, MPI_DOUBLE, MPI_SUM,
                      comm_group);
        double allsum = -Ad_h_y_Ab + P * allred_result;

        weight_predict = eta_armijo * rho;
        bool armijo_cond = false;
        int max_armijo = 10;
        int step_armijo = 0;
        while (armijo_cond == false && step_armijo < max_armijo) {
          double allred_result_loc =
              -2 * alpha * l1l2norm(b_h + weight_predict * d_h) +
              2 * alpha * l1l2b_h;
          MPI_Allreduce(MPI_IN_PLACE, &allred_result_loc, 1, MPI_DOUBLE,
                        MPI_SUM, comm_group);
          if (pow(weight_predict, 2) * Ad_h_sq - weight_predict * Ad_h_y_Ab <=
              sigma * weight_predict * allsum + allred_result_loc) {
            armijo_cond = true;
          } else {
            weight_predict *= rho;
          }
          if (weight_predict < 0.0001 / P) {
            if (my_rank == 0 && !opts.silent) {
              cout << "FBS Solver WARNING: ARMIJO stepsize very small. Value "
                      "is reset to 1/P!"
                   << endl;
            }
            weight_predict = 1.0 / P;
            armijo_cond == true;
          }
          step_armijo++;
        }
        break;
      }
      default: {
        cout << "FBS Solver ERROR: Prediction Rule not known";
        exit(EXIT_FAILURE);
        break;
      }
    }

    timestat(0, 3) += MPI_Wtime() - time_tmp;
    time_tmp = MPI_Wtime();
    // update Xtmp
    if (opts.pred_rule == FISTA || opts.pred_rule == FISTA_RESTART) {
      Xtmp = weight_predict * X_loc + (1 - weight_predict) * X_loc_prev;
      PhiXtmp = weight_predict * PhiX + (1 - weight_predict) * PhiX_prev;
      PhiXtmp_loc =
          weight_predict * PhiX_loc + (1 - weight_predict) * PhiX_loc_prev;
      PhitPhiXtmp_loc = weight_predict * PhitPhiX_loc +
                        (1 - weight_predict) * PhitPhiX_loc_prev;
      PhitPhiXtmp =
          weight_predict * PhitPhiX + (1 - weight_predict) * PhitPhiX_prev;
    } else {
      Xtmp_prev = Xtmp;  // only PSCL
      Xtmp = weight_predict * X_loc + (1 - weight_predict) * Xtmp_prev;

      PhiXtmp = weight_predict * PhiX + (1 - weight_predict) * PhiXtmp;
      PhiXtmp_loc =
          weight_predict * PhiX_loc + (1 - weight_predict) * PhiXtmp_loc;
      PhitPhiXtmp_loc = weight_predict * PhitPhiX_loc +
                        (1 - weight_predict) * PhitPhiXtmp_loc;
      PhitPhiXtmp_prev = PhitPhiXtmp;  // only PSCL
      PhitPhiXtmp =
          weight_predict * PhitPhiX + (1 - weight_predict) * PhitPhiXtmp_prev;
    }

    if (opts.backtracking == DEC || opts.backtracking == DECPAR ||
        opts.backtracking == DECEFF) {
      f_hat = pow((y - PhiXtmp).norm(), 2);
    }

    iterations++;
    iter_need += L;
  }
  // result
  X = Xtmp;
  timestat(0, 0) +=
      timestat(0, 1) + timestat(0, 2) + timestat(0, 3) + timestat(0, 4);
  return;
}

// sequential FBS
void FBSSolver(int& iter_need, double& rel_res, SpEOMatrixD& X, SpEOMatrixD Phi,
               SpEOMatrixD y, double alpha, FBSsolverOptions opts,
               SpEOMatrixD& timestat, SpEOMatrixD& btstat) {
  double eta = 0.5;

  int L = opts.maxiter_in;
  if ((opts.pred_rule == FISTA || opts.pred_rule == FISTA_RESTART) && L > 1) {
    if (!opts.silent) {
      cout << "FBS Solver WARNING: Prediction Rule FISTA is only allowed with "
              "maxiter_in == 1. Value is set to 1!"
           << endl;
    }
    L = 1;
  }

  int Ndecomp = opts.Ndecomp;
  if (opts.Ndecomp > 1) {
    cout << "FBS Solver ERROR: This is the sequential version of the solver. "
            "Ndecomp has to be set to 1"
         << endl;
    exit(EXIT_FAILURE);
  }

  timestat = SpEOMatrixD::Zero(
      1,
      5);  // total, communication, backtracking, prediction, stopping criterion
  double time_tmp = MPI_Wtime();

  double stepsize = 1.0;
  double threshold = stepsize * alpha;
  int max_backtracking_steps = 10;

  int i, iter_in;

  int N = Phi.cols();
  int m = Phi.rows();
  int d = y.cols();

  int total_iter_ls = 0;

  SpEOMatrixD Phit = Phi.transpose();

  // run the recovery algorithm

  int iterations = 0;
  iter_need = 0;
  SpEOMatrixD Phity = Phit * y;
  SpEOMatrixD Xold = SpEOMatrixD::Zero(N, d);
  SpEOMatrixD Xnew = SpEOMatrixD::Zero(N, d);
  SpEOMatrixD Xnew_old = SpEOMatrixD::Zero(N, d);
  SpEOMatrixD Xtmp = SpEOMatrixD::Zero(N, d);
  SpEOMatrixD Xin = SpEOMatrixD::Zero(N, d);
  SpEOMatrixD Xin_old = SpEOMatrixD::Zero(N, d);
  SpEOMatrixD PhiXtmp = SpEOMatrixD::Zero(m, d);
  SpEOMatrixD PhiXnew = SpEOMatrixD::Zero(m, d);
  SpEOMatrixD PhiXnew_old = SpEOMatrixD::Zero(m, d);
  SpEOMatrixD PhiXin = SpEOMatrixD::Zero(m, d);
  SpEOMatrixD PhiXin_old = SpEOMatrixD::Zero(m, d);
  SpEOMatrixD PhiXall = SpEOMatrixD::Zero(m, d);
  SpEOMatrixD PhitPhiXnew_res = SpEOMatrixD::Zero(N, d);
  SpEOMatrixD PhitPhiXtmp = SpEOMatrixD::Zero(N, d);
  SpEOMatrixD PhitPhiXnew = SpEOMatrixD::Zero(N, d);
  SpEOMatrixD PhitPhiXnew_old = SpEOMatrixD::Zero(N, d);
  SpEOMatrixD PhitPhiXold = SpEOMatrixD::Zero(N, d);
  SpEOMatrixD PhitPhiXin = SpEOMatrixD::Zero(N, d);
  SpEOMatrixD PhitPhiXin_old = SpEOMatrixD::Zero(N, d);
  SpEOMatrixD PhitPhiXin_res = SpEOMatrixD::Zero(N, d);

  // only armijo
  SpEOMatrixD Xtmp_old = SpEOMatrixD::Zero(N, d);
  SpEOMatrixD PhitPhiXtmp_all_old = SpEOMatrixD::Zero(N, d);

  // only for backtracking:
  SpEOMatrixD Xnew_prev = SpEOMatrixD::Zero(N, d);
  SpEOMatrixD PhiXnew_prev = SpEOMatrixD::Zero(m, d);

  // normalization for stopping criterion
  double Phitynorm = Phity.norm();

  double Jnew = pow(y.norm(), 2);
  double Jold = 0.0;
  double told = 1.0;
  double tnew = 1.0;
  double weight_predict = 1.0;

  double f_hat;

  bool stop_crit = false;
  while (!stop_crit && (iterations < opts.maxiter_out)) {
    Xin = Xtmp;
    PhiXin = PhiXtmp;
    PhitPhiXin = PhitPhiXtmp;

    timestat(0, 0) += MPI_Wtime() - time_tmp;
    time_tmp = MPI_Wtime();

    double stepsize_loc = stepsize;

    // update L
    if (iterations == 0 || opts.gamma * opts.tol < rel_res) {
      L = 1;
    } else {
      L = max(L, min(opts.maxiter_in,
                     (int)(2 + log10((rel_res / opts.tol) / opts.gamma) *
                                   (opts.maxiter_in - 1.0) /
                                   (-log10(opts.gamma)))));
    }

    double stepsize_init, stepsize_valid;
    for (iter_in = 0; iter_in < L; iter_in++) {
      PhitPhiXin_res = Phity - PhitPhiXin;

      // update stepsize
      stepsize_init = stepsize_loc;
      if (opts.stepsize_update == INIT) {
        stepsize_init = stepsize_loc;
      } else if (opts.stepsize_update == LARGE) {
        stepsize_init = 0.99 * pow(eta, -3);
      } else if (opts.stepsize_update == PSCL) {
        if (iter_in == 0) {
          if (iterations == 0) {
            stepsize_init = min(1 + 1.665 * (1 - m / N), 1.999);
          } else {
            double tmp_sz1 = pow((Xtmp - Xtmp_old).norm(), 2);
            double tmp_sz2 =
                2 * (Xtmp - Xtmp_old)
                        .cwiseProduct(PhitPhiXtmp - PhitPhiXtmp_all_old)
                        .sum();
            Xtmp_old = Xtmp;
            PhitPhiXtmp_all_old = PhitPhiXtmp;
            stepsize_init = max(1.7 * tmp_sz1 / tmp_sz2, 1.0);
          }
        } else {
          double tmp_sz1 = pow((Xin - Xin_old).norm(), 2);
          double tmp_sz2 =
              2 *
              (Xin - Xin_old).cwiseProduct(PhitPhiXin - PhitPhiXin_old).sum();
          stepsize_init = max(1.7 * tmp_sz1 / tmp_sz2, 1.0);
        }
      }
      stepsize_valid = stepsize_init;

      /*
       * UPDATE Xnew and BACKTRACKING
       */
      if (opts.backtracking == DEC || opts.backtracking == DECPAR ||
          opts.backtracking == DECEFF) {
        f_hat = pow((y - PhiXin).norm(), 2);
        stepsize_valid /= eta;
        bool step_cond = false;
        while (!step_cond && stepsize_valid > 0.01) {
          stepsize_valid *= eta;
          threshold = stepsize_valid * alpha;
          Xnew = Xin + stepsize_valid * (PhitPhiXin_res);
          softThresholdJS(threshold, &Phi, Xnew, PhiXnew);
          step_cond =
              0.5 * pow((y - PhiXnew).norm(), 2) <
              0.5 * f_hat - (Xnew - Xin).cwiseProduct(PhitPhiXin_res).sum() +
                  1 / (2 * stepsize_valid) * pow((Xnew - Xin).norm(), 2);
        }
        if (stepsize_valid <= 0.01 && !step_cond && !opts.silent) {
          cout << "FBS Solver WARNING: backtracking stepsize limit reached in "
                  "iteration "
               << iterations << ". Condition violation by "
               << 0.5 * pow((y - PhiXnew).norm(), 2) - 0.5 * f_hat +
                      (Xnew - Xin).cwiseProduct(PhitPhiXin_res).sum() -
                      1 / (2 * stepsize_valid) * pow((Xnew - Xin).norm(), 2)
               << "." << endl;
        }
      } else {
        // Bregman step
        Xnew = Xin + stepsize_valid * (PhitPhiXin_res);
        softThresholdJS(stepsize_valid * alpha, &Phi, Xnew, PhiXnew);
      }

      PhitPhiXnew = Phit * (PhiXnew);
      if (iter_in < L - 1) {
        Xin_old = Xin;
        PhiXin_old = PhiXin;
        PhitPhiXin_old = PhitPhiXin;
        Xin = Xnew;
        PhiXin = PhiXnew;
        PhitPhiXin = PhitPhiXnew;
      }
    }

    timestat(0, 2) += MPI_Wtime() - time_tmp;
    time_tmp = MPI_Wtime();

    /*
     * CHECK STOPPING CRITERION
     */
    switch (opts.stop_crit) {
      case FUNCTIONAL_DIFFERENCE: {
        // update J
        Jold = Jnew;
        Jnew = J_AX(&PhiXnew, &y, alpha, &Xnew);
        rel_res = abs(Jnew - Jold) / abs(Jold);
        stop_crit = rel_res < opts.tol;
        break;
      }
      case FIRST_ORDER_CONDITIONS: {
        PhitPhiXnew_res = Phity - PhitPhiXnew;
        double norm_squared = 0.0;
        for (i = 0; i < N; i++) {
          double normXnew = Xnew.row(i).norm();
          if (normXnew > 1e-16) {
            norm_squared +=
                pow((-PhitPhiXnew_res.row(i) + (alpha / normXnew) * Xnew.row(i))
                        .norm(),
                    2);
          } else {
            norm_squared +=
                pow(max(0.0, PhitPhiXnew_res.row(i).norm() - alpha), 2);
          }
        }
        rel_res = sqrt(norm_squared) / Phitynorm;
        stop_crit = rel_res < opts.tol;

        if (opts.stepsize_adaptive == true) {
          double stepsize_new = stepsize;
          SpEOMatrixD Dgrad = PhitPhiXnew - PhitPhiXold;
          SpEOMatrixD DX = Xnew - Xold;
          Xold = Xnew;
          PhitPhiXold = PhitPhiXnew;
          double DXnorm_sq = pow(DX.norm(), 2);
          double DX_Dgrad = DX.cwiseProduct(Dgrad).sum();
          double Dgradnorm_sq = pow(Dgrad.norm(), 2);
          double tau_s = DXnorm_sq / DX_Dgrad;
          double tau_m = max(DX_Dgrad / Dgradnorm_sq, 0.0);
          if (2 * tau_m > tau_s) {
            stepsize_new = tau_m;
          } else {
            stepsize_new = tau_s - 0.5 * tau_m;
          }
          if (stepsize_new > 0) {
            stepsize = stepsize_new;
          }
        }

        break;
      }
      default: {
        cout << "FBS Solver ERROR: Stopping Rule not known";
        exit(EXIT_FAILURE);
        break;
      }
    }
    timestat(0, 4) += MPI_Wtime() - time_tmp;
    time_tmp = MPI_Wtime();
    /*
     * PREDICTION
     */
    switch (opts.pred_rule) {
      case FISTA:  // guaranteed to converge only for param.maxiter_in == 1
      {
        told = tnew;
        tnew = (1 + sqrt(1 + 4 * pow(told, 2))) / 2;
        weight_predict = (1 + (told - 1) / tnew);
        break;
      }
      case FISTA_RESTART:  // guaranteed to converge only for param.maxiter_in
                           // == 1
      {
        told = tnew;
        if (((Xtmp - Xnew).cwiseProduct(Xnew - Xnew_old)).sum() >=
            0) {  // restart
          told = 1;
        }
        tnew = (1 + sqrt(1 + 4 * pow(told, 2))) / 2;
        weight_predict = (1 + (told - 1) / tnew);
        break;
      }
      case STANDARD: {
        weight_predict = 1.0 / Ndecomp;  // standard FBS guaranteed to converge
        break;
      }
      case STANDARD_AGGRESSIVE: {
        weight_predict = 1.0;  // guaranteed to converge only for Ndecomp == 1
                               // or param.maxiter_in == 1
        break;
      }
      case LINE_SEARCH_SIMPLEX: {
        weight_predict = 1.0 / Ndecomp;
        double weight_eps = 1e-2;

        SpEOMatrixD d_h = Xnew - Xtmp;
        SpEOMatrixD b_h = Xtmp;
        SpEOMatrixD Ad_h = PhiXnew - PhiXtmp;
        SpEOMatrixD y_Ab = y - PhiXtmp;
        double Ad_h_sq = pow(Ad_h.norm(), 2);
        double y_Ab_sq = pow(y_Ab.norm(), 2);
        double Ad_h_y_Ab = 2 * Ad_h.cwiseProduct(y_Ab).sum();
        SpEOVectorD d_s = d_h.cwiseAbs2().rowwise().sum();
        SpEOVectorD b_s = b_h.cwiseAbs2().rowwise().sum();
        SpEOVectorD d_s_b_s = 2 * d_h.cwiseProduct(b_h).rowwise().sum();
        int iter_loc;
        double w = J_minsearch(Ad_h_sq, Ad_h_y_Ab, y_Ab_sq, d_s, d_s_b_s, b_s,
                               alpha, 0, weight_eps, iter_loc);
        if (w > 0) {
          weight_predict = w;
        }
        total_iter_ls += iter_loc;
        break;
      }
      case LINE_SEARCH_ARMIJO: {
        weight_predict = 1.0 / Ndecomp;
        double rho = 0.5;
        double sigma = 0.5;
        double eta_armijo = 2.0;
        SpEOMatrixD d_h = Xnew - Xtmp;
        SpEOMatrixD b_h = Xtmp;
        SpEOMatrixD Ad_h = PhiXnew - PhiXtmp;
        SpEOMatrixD y_Ab = y - PhiXtmp;
        double Ad_h_sq = pow(Ad_h.norm(), 2);
        double y_Ab_sq = pow(y_Ab.norm(), 2);
        double Ad_h_y_Ab = 2 * Ad_h.cwiseProduct(y_Ab).sum();

        double l1l2b_h = l1l2norm(b_h);
        double allred_result =
            2 * alpha * (l1l2norm(b_h + 1.0 / Ndecomp * d_h) - l1l2b_h);
        double allsum = -Ad_h_y_Ab + Ndecomp * allred_result;

        weight_predict = eta_armijo * rho;
        bool armijo_cond = false;
        while (armijo_cond == false) {
          double allred_result_loc =
              -2 * alpha * l1l2norm(b_h + weight_predict * d_h) +
              2 * alpha * l1l2b_h;
          if (pow(weight_predict, 2) * Ad_h_sq - weight_predict * Ad_h_y_Ab <=
              sigma * weight_predict * allsum + allred_result_loc) {
            armijo_cond = true;
          } else {
            weight_predict *= rho;
          }
          if (weight_predict < 0.0001 / Ndecomp) {
            cout << "Problem in ARMIJO " << Ndecomp << " " << L << endl;
            weight_predict = 1.0 / Ndecomp;
            armijo_cond == true;
          }
        }
        break;
      }
      default: {
        cout << "FBS Solver ERROR: Stopping Rule not known";
        exit(EXIT_FAILURE);
        break;
      }
    }
    timestat(0, 3) += MPI_Wtime() - time_tmp;
    time_tmp = MPI_Wtime();
    // update Xtmp
    if (opts.pred_rule == FISTA || opts.pred_rule == FISTA_RESTART) {
      Xtmp = weight_predict * Xnew + (1 - weight_predict) * Xnew_old;
      PhiXtmp = weight_predict * PhiXnew + (1 - weight_predict) * PhiXnew_old;
      PhitPhiXtmp =
          weight_predict * PhitPhiXnew + (1 - weight_predict) * PhitPhiXnew_old;
      Xnew_old = Xnew;
      PhiXnew_old = PhiXnew;
      PhitPhiXnew_old = PhitPhiXnew;
    } else {
      Xtmp = weight_predict * Xnew + (1 - weight_predict) * Xtmp;
      PhiXtmp = weight_predict * PhiXnew + (1 - weight_predict) * PhiXtmp;
      PhitPhiXtmp =
          weight_predict * PhitPhiXnew + (1 - weight_predict) * PhitPhiXtmp;
    }

    iterations++;
    iter_need += L;
  }
  // result
  X = Xtmp;

  timestat(0, 0) +=
      timestat(0, 1) + timestat(0, 2) + timestat(0, 3) + timestat(0, 4);
  return;
}

void solve_equation3(int& iter, double& rel_res, SpEOMatrixD& Z, SpEOVectorD m,
                     SpEOVectorD& m_out, SpEOMatrixD* DH, SpEOMatrixD* DL,
                     SpEOVectorD* P_lmd_vecs, SpEOMatrixI* P_lmd_idx_row,
                     int** P_lmd_idx_bl, SpEOMatrixD R, SpEOMatrixD X,
                     SpEOMatrixD* Y, SpEOMatrixD* A, SpEOMatrixD*& A_out,
                     SpEOMatrixD Z0, double lambda_X, double lambda_Y, int N_g,
                     int maxiter, double tol_r, int fix_Alpha,
                     bool fix_delta_m) {
  bool use_starting_value = 1;

  int j, k;  //, iter;
  // compute global stuff
  int N_Y = m.size();
  int N_h = X.cols();
  int N_X = X.rows();
  int N_l = Y[0].rows();
  int sum_omega = 0;
  for (k = 0; k < N_g; k++) {
    sum_omega += DH[k].cols();
  }
  int sum_N_c = 0;
  for (k = 0; k < N_g; k++) {
    sum_N_c += A[k].cols();
  }
  double coeff1 = sqrt(0.5 / (N_Y * N_h));
  double coeff2 = sqrt(0.5 * lambda_X / (N_X * N_h));
  double coeff3 = sqrt(0.5 * lambda_Y / (N_l * sum_N_c));
  SpEOVectorD ones = SpEOVectorD::Ones(N_h);

  // compute initial right-hand side
  SpEOMatrixD B_upper = coeff1 * (m * ones.transpose());
  if (fix_Alpha == 1) {
    SpEOMatrixD tmpDA = SpEOMatrixD::Zero(sum_N_c, N_h);
    int running_idx = 0;
    for (k = 0; k < N_g; k++) {
      tmpDA.block(running_idx, 0, A[k].cols(), N_h) =
          (DH[k] * A[k]).transpose();
      running_idx += A[k].cols();
    }
    SpEOMatrixD tmpPDA = SpEOMatrixD::Zero(N_Y, N_h);
    for (j = 0; j < N_Y; j++) {
      for (k = 0; k < P_lmd_idx_row[j].cols(); k++) {
        int block = P_lmd_idx_row[j].coeff(0, k);
        int relidx = P_lmd_idx_row[j].coeff(1, k);
        tmpPDA.row(j) += P_lmd_vecs[block](relidx) *
                         tmpDA.row(P_lmd_idx_bl[block][1] + relidx);
      }
    }
    B_upper += coeff1 * tmpPDA;
  }
  SpEOMatrixD B_lower = coeff2 * X;
  SpEOMatrixD* B_lower2 = new SpEOMatrixD[N_g];
  for (k = 0; k < N_g; k++) {
    B_lower2[k] = coeff3 * Y[k];
    if (fix_Alpha == 1) {
      B_lower2[k] -= coeff3 * (DL[k] * A[k]);
    }
  }

  // compute initial values for CGLS
  SpEOMatrixD x_upper = SpEOMatrixD::Zero(N_Y, N_h);  // alias Z^p
  if (use_starting_value) {
    x_upper = Z0;
  }
  SpEOMatrixD* x_center = new SpEOMatrixD[N_g];  // alias A_i(Omega_i)'s
  for (k = 0; k < N_g; k++) {
    x_center[k] = SpEOMatrixD::Zero(A[k].rows(), A[k].cols());
    if (use_starting_value) {
      x_center[k] = A[k];
    }
  }
  SpEOVectorD x_lower = SpEOVectorD::Zero(N_Y);  // alias m

  // compute d= right-hand side
  SpEOMatrixD d_upper = B_upper;
  SpEOMatrixD d_lower = B_lower;
  SpEOMatrixD* d_lower2 = new SpEOMatrixD[N_g];
  for (k = 0; k < N_g; k++) {
    d_lower2[k] = B_lower2[k];
  }
  SpEOMatrixD d_upper_tmp = d_upper;
  SpEOMatrixD d_lower_tmp = d_lower;
  SpEOMatrixD* d_lower2_tmp = new SpEOMatrixD[N_g];
  for (k = 0; k < N_g; k++) {
    d_lower2_tmp[k] = d_lower2[k];
  }

  if (use_starting_value) {  // d -= operator * x (only take the A's into
                             // account since delta_m_0 = 0, Z_0 = 0)
    d_upper -= coeff1 * Z0;
    if (fix_Alpha == 0 || fix_Alpha == 2) {
      SpEOMatrixD tmpDA = SpEOMatrixD::Zero(sum_N_c, N_h);
      int running_idx = 0;
      for (k = 0; k < N_g; k++) {
        tmpDA.block(running_idx, 0, A[k].cols(), N_h) =
            (DH[k] * A[k]).transpose();
        running_idx += A[k].cols();
      }
      SpEOMatrixD tmpPDA = SpEOMatrixD::Zero(N_Y, N_h);
      for (j = 0; j < N_Y; j++) {
        for (k = 0; k < P_lmd_idx_row[j].cols(); k++) {
          int block = P_lmd_idx_row[j].coeff(0, k);
          int relidx = P_lmd_idx_row[j].coeff(1, k);
          tmpPDA.row(j) += P_lmd_vecs[block](relidx) *
                           tmpDA.row(P_lmd_idx_bl[block][1] + relidx);
        }
      }
      d_upper -= coeff1 * (-tmpPDA);
    }
    d_lower -= coeff2 * (R * Z0);
    if (fix_Alpha == 0 || fix_Alpha == 2) {
      for (k = 0; k < N_g; k++) {
        d_lower2[k] -= coeff3 * (DL[k] * A[k]);
      }
    }
  }

  // compute r = operator.transposed()*d;
  SpEOMatrixD r_upper = coeff1 * d_upper + coeff2 * (R.transpose() * d_lower);
  SpEOMatrixD* r_center = new SpEOMatrixD[N_g];
  for (k = 0; k < N_g; k++) {
    SpEOMatrixD temp = SpEOMatrixD::Zero(N_h, A[k].cols());
    for (j = 0; j < A[k].cols(); j++) {  // write temp = B^t*P(k)
      temp.col(j) =
          (P_lmd_vecs[k](j)) * d_upper.row(P_lmd_idx_bl[k][0] + j).transpose();
    }
    r_center[k] = -coeff1 * (DH[k].transpose() * temp) +
                  coeff3 * (DL[k].transpose() * d_lower2[k]);
    if (fix_Alpha == 2) {
      r_center[k].row(0) = SpEOMatrixD::Zero(1, A[k].cols());
    }
  }
  SpEOVectorD r_lower = -coeff1 * d_upper.rowwise().sum();

  // compute r = operator.transposed()*d;
  SpEOMatrixD r_upper_tmp =
      coeff1 * d_upper_tmp + coeff2 * (R.transpose() * d_lower_tmp);
  SpEOMatrixD* r_center_tmp = new SpEOMatrixD[N_g];
  for (k = 0; k < N_g; k++) {
    SpEOMatrixD temp = SpEOMatrixD::Zero(N_h, A[k].cols());
    for (j = 0; j < A[k].cols(); j++) {  // write temp = B^t*P(k)
      temp.col(j) = (P_lmd_vecs[k](j)) *
                    d_upper_tmp.row(P_lmd_idx_bl[k][0] + j).transpose();
    }
    r_center_tmp[k] = -coeff1 * (DH[k].transpose() * temp) +
                      coeff3 * (DL[k].transpose() * d_lower2_tmp[k]);
  }
  SpEOVectorD r_lower_tmp = -coeff1 * d_upper_tmp.rowwise().sum();

  // compute p = r
  SpEOMatrixD p_upper = r_upper;
  SpEOMatrixD* p_center = new SpEOMatrixD[N_g];
  for (k = 0; k < N_g; k++) {
    p_center[k] = r_center[k];
    if (fix_Alpha == 2) {
      p_center[k].row(0) = SpEOMatrixD(1, A[k].cols());
    }
  }
  SpEOVectorD p_lower = r_lower;

  // compute t = operator*p;
  SpEOMatrixD t_upper = coeff1 * (p_upper);
  if (fix_Alpha == 0 || fix_Alpha == 2) {
    SpEOMatrixD tmpDA = SpEOMatrixD::Zero(sum_N_c, N_h);
    int running_idx = 0;
    for (k = 0; k < N_g; k++) {
      tmpDA.block(running_idx, 0, A[k].cols(), N_h) =
          (DH[k] * p_center[k]).transpose();
      running_idx += A[k].cols();
    }
    SpEOMatrixD tmpPDA = SpEOMatrixD::Zero(N_Y, N_h);
    for (j = 0; j < N_Y; j++) {
      for (k = 0; k < P_lmd_idx_row[j].cols(); k++) {
        int block = P_lmd_idx_row[j].coeff(0, k);
        int relidx = P_lmd_idx_row[j].coeff(1, k);
        tmpPDA.row(j) += P_lmd_vecs[block](relidx) *
                         tmpDA.row(P_lmd_idx_bl[block][1] + relidx);
      }
    }
    t_upper += coeff1 * (-tmpPDA);
  }
  if (!fix_delta_m) {
    t_upper += coeff1 * (-p_lower * ones.transpose());
  }
  SpEOMatrixD t_lower = coeff2 * (R * p_upper);
  SpEOMatrixD* t_lower2 = new SpEOMatrixD[N_g];
  for (k = 0; k < N_g; k++) {
    t_lower2[k] = SpEOMatrixD::Zero(DL[k].rows(), A[k].cols());
    if (fix_Alpha == 0 || fix_Alpha == 2) {
      t_lower2[k] = coeff3 * (DL[k] * p_center[k]);
    }
  }

  double alpha, beta;
  double r_norm_squared = pow(r_upper.norm(), 2);
  double norm_rc = 0;
  if (fix_Alpha == 0 || fix_Alpha == 2) {
    for (k = 0; k < N_g; k++) {
      norm_rc += pow(r_center[k].norm(), 2);
    }
    r_norm_squared += norm_rc;
  }
  if (!fix_delta_m) {
    r_norm_squared += pow(r_lower.norm(), 2);
  }
  double r0_norm = pow(r_upper_tmp.norm(), 2);
  double norm_rc_tmp = 0;
  if (fix_Alpha == 0 || fix_Alpha == 2) {
    for (k = 0; k < N_g; k++) {
      norm_rc_tmp += pow(r_center_tmp[k].norm(), 2);
    }
    r0_norm += norm_rc_tmp;
  }
  if (!fix_delta_m) {
    r0_norm += pow(r_lower_tmp.norm(), 2);
  }
  r0_norm = sqrt(r0_norm);

  iter = 0;
  while (iter < maxiter && sqrt(r_norm_squared) / r0_norm > tol_r) {
    double t_lower2_sq = 0;
    for (k = 0; k < N_g; k++) {
      t_lower2_sq += pow(t_lower2[k].norm(), 2);
    }
    alpha = r_norm_squared /
            (pow(t_lower.norm(), 2) + pow(t_upper.norm(), 2) + t_lower2_sq);

    x_upper += alpha * p_upper;
    if (fix_Alpha == 0 || fix_Alpha == 2) {
      for (k = 0; k < N_g; k++) {
        x_center[k] += alpha * p_center[k];
      }
    }
    if (!fix_delta_m) {
      x_lower += alpha * p_lower;
    }

    d_upper -= alpha * t_upper;
    d_lower -= alpha * t_lower;
    for (k = 0; k < N_g; k++) {
      d_lower2[k] -= alpha * t_lower2[k];
    }

    r_upper = coeff1 * d_upper + coeff2 * (R.transpose() * d_lower);
    if (fix_Alpha == 0 || fix_Alpha == 2) {
      for (k = 0; k < N_g; k++) {
        SpEOMatrixD temp = SpEOMatrixD::Zero(N_h, A[k].cols());
        for (j = 0; j < A[k].cols(); j++) {  // write temp = B^t*P(k)
          temp.col(j) = (P_lmd_vecs[k](j)) *
                        d_upper.row(P_lmd_idx_bl[k][0] + j).transpose();
        }
        r_center[k] = -coeff1 * (DH[k].transpose() * temp) +
                      coeff3 * (DL[k].transpose() * d_lower2[k]);
        if (fix_Alpha == 2) {
          r_center[k].row(0) = SpEOMatrixD::Zero(1, A[k].cols());
        }
      }
    }
    if (!fix_delta_m) {
      r_lower = -coeff1 * d_upper.rowwise().sum();
    }

    double r_norm_squared_old = r_norm_squared;
    r_norm_squared = pow(r_upper.norm(), 2);
    norm_rc = 0;
    if (fix_Alpha == 0 || fix_Alpha == 2) {
      for (k = 0; k < N_g; k++) {
        norm_rc += pow(r_center[k].norm(), 2);
      }
      r_norm_squared += norm_rc;
    }
    if (!fix_delta_m) {
      r_norm_squared += pow(r_lower.norm(), 2);
    }
    beta = r_norm_squared / r_norm_squared_old;

    p_upper = r_upper + beta * p_upper;
    if (fix_Alpha == 0 || fix_Alpha == 2) {
      for (k = 0; k < N_g; k++) {
        p_center[k] = r_center[k] + beta * p_center[k];
        if (fix_Alpha == 2) {
          p_center[k].row(0) = SpEOMatrixD::Zero(1, A[k].cols());
        }
      }
    }
    if (!fix_delta_m) {
      p_lower = r_lower + beta * p_lower;
    }

    t_upper = coeff1 * (p_upper);
    if (fix_Alpha == 0 || fix_Alpha == 2) {
      SpEOMatrixD tmpDA = SpEOMatrixD::Zero(sum_N_c, N_h);
      int running_idx = 0;
      for (k = 0; k < N_g; k++) {
        tmpDA.block(running_idx, 0, A[k].cols(), N_h) =
            (DH[k] * p_center[k]).transpose();
        running_idx += A[k].cols();
      }
      SpEOMatrixD tmpPDA = SpEOMatrixD::Zero(N_Y, N_h);
      for (j = 0; j < N_Y; j++) {
        for (k = 0; k < P_lmd_idx_row[j].cols(); k++) {
          int block = P_lmd_idx_row[j].coeff(0, k);
          int relidx = P_lmd_idx_row[j].coeff(1, k);
          tmpPDA.row(j) += P_lmd_vecs[block](relidx) *
                           tmpDA.row(P_lmd_idx_bl[block][1] + relidx);
        }
      }
      t_upper += coeff1 * (-tmpPDA);
    }

    if (!fix_delta_m) {
      t_upper += coeff1 * (-p_lower * ones.transpose());
    }
    t_lower = coeff2 * (R * p_upper);
    if (fix_Alpha == 0 || fix_Alpha == 2) {
      for (k = 0; k < N_g; k++) {
        t_lower2[k] = coeff3 * (DL[k] * p_center[k]);
      }
    }

    iter++;
  }

  Z = x_upper;
  A_out = new SpEOMatrixD[N_g];
  for (k = 0; k < N_g; k++) {
    A_out[k] = x_center[k];
  }
  m_out = x_lower;

  rel_res = sqrt(r_norm_squared) / r0_norm;

  return;
}
void solve_equation1(int& iter, double& rel_res, SpEOMatrixD*& I_Z,
                     SpEOMatrixD R, int filter_size, int fDS, int N_Y, int N_X,
                     SpEOMatrixD* I_Z_tilde, SpEOMatrixD* I_X, SpEOMatrixD* I_Y,
                     double lambda_X, double lambda_Y, int maxiter,
                     double tol_r, MPI_Comm eq1_comm, int channels_per_core,
                     bool SNR_normalization, bool balance_ImX_term_coef) {
  bool use_starting_value = false;

  // compute global stuff

  int my_local_rank;
  MPI_Comm_rank(eq1_comm, &my_local_rank);

  int my_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

  int N_l = I_Y[0].cols() * I_Y[0].rows();
  int N_h = I_X[0].cols() * I_X[0].rows();
  int HR_cols = I_Z_tilde[0].cols();
  int HR_rows = I_Z_tilde[0].rows();

  int channel_Y, channel_X;

  double coeff1 = sqrt(0.5 / (N_Y * N_h));
  double coeff2 = sqrt(0.5 * lambda_X / (N_X * N_h));
  double coeff3 = sqrt(0.5 * lambda_Y / (N_Y * N_l));

  if (SNR_normalization) {
    coeff1 = sqrt(0.5);
    coeff2 = sqrt(0.5 * lambda_X);
    coeff3 = sqrt(0.5 * lambda_Y);
    if (balance_ImX_term_coef) {
      coeff2 *= sqrt((double)N_Y / (double)N_X);
    }
  }
  // define Gaussian filter
  SpEOMatrixD gauss_filter;
  create_Gaussian_filter(gauss_filter, filter_size);
  SpEOMatrixD filter_coeff;
  calc_filter_boundary_coeff(filter_coeff, gauss_filter, fDS);

  // define operator
  Eq1Op op = Eq1Op(&R, gauss_filter, filter_coeff, fDS, coeff1, coeff2, coeff3);

  // compute initial right-hand side
  Eq1MVar b = Eq1MVar(channels_per_core, N_X, channels_per_core);
  Eq1MVar d = Eq1MVar(channels_per_core, N_X, channels_per_core);
  for (channel_Y = 0; channel_Y < channels_per_core; channel_Y++) {
    b.upper[channel_Y] = coeff1 * I_Z_tilde[channel_Y];
    b.lower[channel_Y] = coeff3 * I_Y[channel_Y];
    // d = right-hand side
    d.upper[channel_Y] = b.upper[channel_Y];
    d.lower[channel_Y] = b.lower[channel_Y];
  }
  for (channel_X = 0; channel_X < N_X; channel_X++) {
    b.center[channel_X] = coeff2 * I_X[channel_X];
    // d = right-hand side
    d.center[channel_X] = b.center[channel_X];
  }

  // compute initial values for CGLS
  Eq1SVar x = Eq1SVar(channels_per_core);
  for (channel_Y = 0; channel_Y < channels_per_core; channel_Y++) {
    if (use_starting_value) {
      // x = I_Z_tilde //= I_Z_0
      x.var[channel_Y] = I_Z_tilde[channel_Y];
    } else {
      // x = 0
      x.var[channel_Y] = SpEOMatrixD::Zero(HR_rows, HR_cols);
    }
  }

  Eq1SVar r = Eq1SVar(channels_per_core);
  // compute norm(operator.transposed()*d);
  r.SetAdjointOperatorXvar(op, &d);

  double r0_norm = sqrt(r.NormSquared(eq1_comm));
  cout << endl << endl << "### r0_norm = " << r0_norm << endl << endl;

  if (use_starting_value) {
    // d =  d - operator * x
    d.SubtractOperatorXvar(op, I_Z_tilde, eq1_comm);
  }

  // compute r = operator.transposed()*d;
  r.SetAdjointOperatorXvar(op, &d);

  // p = r
  Eq1SVar p = Eq1SVar(channels_per_core);
  for (channel_Y = 0; channel_Y < channels_per_core; channel_Y++) {
    p.var[channel_Y] = r.var[channel_Y];
  }

  // t = operator * p;
  Eq1MVar t = Eq1MVar(channels_per_core, N_X, channels_per_core);
  for (channel_Y = 0; channel_Y < channels_per_core; channel_Y++) {
    t.upper[channel_Y] = SpEOMatrixD::Zero((I_Z_tilde[channel_Y]).rows(),
                                           (I_Z_tilde[channel_Y]).cols());
    t.lower[channel_Y] =
        SpEOMatrixD::Zero((I_Y[channel_Y]).rows(), (I_Y[channel_Y]).cols());
  }
  for (channel_X = 0; channel_X < N_X; channel_X++) {
    t.center[channel_X] =
        SpEOMatrixD::Zero((I_X[channel_X]).rows(), (I_X[channel_X]).cols());
  }
  t.SetOperatorXvar(op, &p, eq1_comm);

  double alpha, beta;
  double r_norm_squared = r.NormSquared(eq1_comm);

  if (my_rank == 0) {
    cout << "resudial in all iterations:" << endl;
  }
  iter = 0;
  while (iter < maxiter && sqrt(r_norm_squared) / r0_norm > tol_r) {
    if (my_rank == 0) {
      cout << sqrt(r_norm_squared) / r0_norm << endl;
    }
    // alpha = (norm(r)/norm(t))^2
    alpha = r_norm_squared / t.NormSquared(eq1_comm);

    // x = x + alpha * p ; d = d - alpha * t
    for (channel_Y = 0; channel_Y < channels_per_core; channel_Y++) {
      x.var[channel_Y] += alpha * p.var[channel_Y];

      d.upper[channel_Y] -= alpha * t.upper[channel_Y];
      d.lower[channel_Y] -= alpha * t.lower[channel_Y];
    }
    for (channel_X = 0; channel_X < N_X; channel_X++) {
      d.center[channel_X] -= alpha * t.center[channel_X];
    }

    // r = operator^T * d
    r.SetAdjointOperatorXvar(op, &d);

    // beta = (norm(r)/norm(r_old))^2
    double r_norm_squared_old = r_norm_squared;
    r_norm_squared = r.NormSquared(eq1_comm);
    beta = r_norm_squared / r_norm_squared_old;

    // p = r + beta * p
    for (channel_Y = 0; channel_Y < channels_per_core; channel_Y++) {
      p.var[channel_Y] = r.var[channel_Y] + beta * p.var[channel_Y];
    }

    // t = operator * p
    t.SetOperatorXvar(op, &p, eq1_comm);

    iter++;
  }

  for (channel_Y = 0; channel_Y < channels_per_core; channel_Y++) {
    I_Z[channel_Y] = x.var[channel_Y];
  }
  rel_res = sqrt(r_norm_squared) / r0_norm;

  return;
}

void solve_equation1_unmixing(
    bool use_starting_value, int& iter, double& rel_res,
    SpEOMatrixD& EndmemberMat, SpEOMatrixD*& AbundanceMat, SpEOMatrixD*& I_Z,
    SpEOMatrixD R, int filter_size, int fDS, int N_Y, int N_X,
    SpEOMatrixD* I_Z_tilde, SpEOMatrixD* I_X, SpEOMatrixD* I_Y, double lambda_X,
    double lambda_Y, int maxiter, double tol_r, MPI_Comm eq1_comm,
    int channels_per_core, bool SNR_normalization, bool balance_ImX_term_coef,
    SpEOVectorD& LAMBDA_ZY) {
  // compute global stuff

  int my_local_rank;
  MPI_Comm_rank(eq1_comm, &my_local_rank);

  int my_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

  int N_Y_sub = EndmemberMat.cols();

  int N_l = I_Y[0].cols() * I_Y[0].rows();
  int N_h = I_X[0].cols() * I_X[0].rows();
  int HR_cols = I_Z_tilde[0].cols();
  int HR_rows = I_Z_tilde[0].rows();

  int channel_Y, channel_X;
  double coeff1 = sqrt(0.5 / (N_Y * N_h));
  double coeff2 = sqrt(0.5 * lambda_X / (N_X * N_h));
  double coeff3 = sqrt(0.5 * lambda_Y / (N_Y * N_l));
  if (SNR_normalization) {
    coeff1 = sqrt(0.5);
    coeff2 = sqrt(0.5 * lambda_X);
    coeff3 = sqrt(0.5 * lambda_Y);
    if (balance_ImX_term_coef) {
      coeff2 *= sqrt(((double)N_Y) / ((double)N_X));
    }
  }

  // define Gaussian filter
  SpEOMatrixD gauss_filter;
  create_Gaussian_filter(gauss_filter, filter_size);
  SpEOMatrixD filter_coeff;
  calc_filter_boundary_coeff(filter_coeff, gauss_filter, fDS);

  // define operator
  Eq1Op op = Eq1Op(&R, gauss_filter, filter_coeff, fDS, coeff1, coeff2, coeff3);

  // compute initial right-hand side
  Eq1MVar b = Eq1MVar(N_Y, N_X, N_Y);
  Eq1MVar d = Eq1MVar(N_Y, N_X, N_Y);
  for (channel_Y = 0; channel_Y < N_Y; channel_Y++) {
    b.upper[channel_Y] = coeff1 * I_Z_tilde[channel_Y];
    b.lower[channel_Y] = coeff3 * I_Y[channel_Y];
    // d = right-hand side
    d.upper[channel_Y] = b.upper[channel_Y];
    d.lower[channel_Y] = b.lower[channel_Y];
  }
  for (channel_X = 0; channel_X < N_X; channel_X++) {
    b.center[channel_X] = coeff2 * I_X[channel_X];
    // d = right-hand side
    d.center[channel_X] = b.center[channel_X];
  }

  // compute initial values for CGLS
  Eq1SVar x = Eq1SVar(N_Y_sub);
  for (channel_Y = 0; channel_Y < N_Y_sub; channel_Y++) {
    if (use_starting_value) {
      // x = I_Z_tilde //= I_Z_0
      x.var[channel_Y] = AbundanceMat[channel_Y];
    } else {
      // x = 0
      x.var[channel_Y] = SpEOMatrixD::Zero(HR_rows, HR_cols);
    }
  }

  Eq1SVar_sub r = Eq1SVar_sub(N_Y_sub);
  // compute norm(operator.transposed()*d);
  r.SetAdjointOperatorXvar(op, &d, EndmemberMat, LAMBDA_ZY);

  double r0_norm = sqrt(r.NormSquared(eq1_comm));
  cout << endl << endl << "### r0_norm = " << r0_norm << endl << endl;
  if (use_starting_value) {
    // d =  d - operator * x
    d.SubtractOperatorXvar_sub(op, AbundanceMat, eq1_comm, EndmemberMat,
                               LAMBDA_ZY);
  }

  r.SetAdjointOperatorXvar(op, &d, EndmemberMat, LAMBDA_ZY);

  // p = r
  Eq1SVar_sub p = Eq1SVar_sub(N_Y_sub);
  for (channel_Y = 0; channel_Y < N_Y_sub; channel_Y++) {
    p.var[channel_Y] = r.var[channel_Y];
  }

  // t = operator * p;
  Eq1MVar t = Eq1MVar(N_Y, N_X, N_Y);
  for (channel_Y = 0; channel_Y < N_Y; channel_Y++) {
    t.upper[channel_Y] = SpEOMatrixD::Zero((I_Z_tilde[channel_Y]).rows(),
                                           (I_Z_tilde[channel_Y]).cols());
    t.lower[channel_Y] =
        SpEOMatrixD::Zero((I_Y[channel_Y]).rows(), (I_Y[channel_Y]).cols());
  }
  for (channel_X = 0; channel_X < N_X; channel_X++) {
    t.center[channel_X] =
        SpEOMatrixD::Zero((I_X[channel_X]).rows(), (I_X[channel_X]).cols());
  }
  t.SetOperatorXvar(op, &p, eq1_comm, EndmemberMat, LAMBDA_ZY);

  double alpha, beta;
  double r_norm_squared = r.NormSquared(eq1_comm);

  iter = 0;
  if (my_rank == 0) {
    cout << "residual in all iterations: " << endl;
  }
  while (iter < maxiter && sqrt(r_norm_squared) / r0_norm > tol_r) {
    if (my_rank == 0) {
      cout << sqrt(r_norm_squared) / r0_norm << endl;
    }
    // alpha = (norm(r)/norm(t))^2
    alpha = r_norm_squared / t.NormSquared(eq1_comm);

    // x = x + alpha * p ; d = d - alpha * t
    for (channel_Y = 0; channel_Y < N_Y_sub; channel_Y++) {
      x.var[channel_Y] += alpha * p.var[channel_Y];
    }
    for (channel_Y = 0; channel_Y < N_Y; channel_Y++) {
      d.upper[channel_Y] -= alpha * t.upper[channel_Y];
      d.lower[channel_Y] -= alpha * t.lower[channel_Y];
    }
    for (channel_X = 0; channel_X < N_X; channel_X++) {
      d.center[channel_X] -= alpha * t.center[channel_X];
    }

    // r = operator^T * d
    r.SetAdjointOperatorXvar(op, &d, EndmemberMat, LAMBDA_ZY);

    // beta = (norm(r)/norm(r_old))^2
    double r_norm_squared_old = r_norm_squared;
    r_norm_squared = r.NormSquared(eq1_comm);
    beta = r_norm_squared / r_norm_squared_old;
    for (channel_Y = 0; channel_Y < N_Y_sub; channel_Y++) {
      p.var[channel_Y] = r.var[channel_Y] + beta * p.var[channel_Y];
    }
    t.SetOperatorXvar(op, &p, eq1_comm, EndmemberMat, LAMBDA_ZY);

    iter++;
  }
  for (channel_Y = 0; channel_Y < N_Y_sub; channel_Y++) {
    AbundanceMat[channel_Y] = x.var[channel_Y];
  }
  bool transpose_2D_mat = false;
  multiply_2D_with_3D_mat(I_Z, EndmemberMat, x.var, transpose_2D_mat);
  rel_res = sqrt(r_norm_squared) / r0_norm;

  return;
}

void solve_equation3_grad(SpEOMatrixD& Z, SpEOVectorD m, SpEOVectorD& m_out,
                          SpEOMatrixD* DH, SpEOMatrixD* DL,
                          SpEOVectorD* P_lmd_vecs, SpEOMatrixI* P_lmd_idx_row,
                          int** P_lmd_idx_bl, SpEOMatrixD R, SpEOMatrixD X,
                          SpEOMatrixD* Y, SpEOMatrixD* A, SpEOMatrixD*& A_out,
                          SpEOMatrixD Z0, double lambda_X, double lambda_Y,
                          int N_g, int maxiter, double tol_r, bool fix_Alpha,
                          bool fix_delta_m) {
  bool use_starting_value = 1;

  int j, k;
  // compute global stuff
  int N_Y = m.size();
  int N_h = X.cols();
  int N_X = X.rows();
  int N_l = Y[0].rows();
  int sum_omega = 0;
  for (k = 0; k < N_g; k++) {
    sum_omega += DH[k].cols();
  }
  int sum_N_c = 0;
  for (k = 0; k < N_g; k++) {
    sum_N_c += A[k].cols();
  }
  double coeff1 = sqrt(0.5 / (N_Y * N_h));
  double coeff2 = sqrt(0.5 * lambda_X / (N_X * N_h));
  double coeff3 = sqrt(0.5 * lambda_Y / (N_l * sum_N_c));
  SpEOVectorD ones = SpEOVectorD::Ones(N_h);

  // compute initial right-hand side
  SpEOMatrixD B_upper = coeff1 * (m * ones.transpose());
  if (fix_Alpha) {
    SpEOMatrixD tmpDA = SpEOMatrixD::Zero(sum_N_c, N_h);
    int running_idx = 0;
    for (k = 0; k < N_g; k++) {
      tmpDA.block(running_idx, 0, A[k].cols(), N_h) =
          (DH[k] * A[k]).transpose();
      running_idx += A[k].cols();
    }
    SpEOMatrixD tmpPDA = SpEOMatrixD::Zero(N_Y, N_h);
    for (j = 0; j < N_Y; j++) {
      for (k = 0; k < P_lmd_idx_row[j].cols(); k++) {
        int block = P_lmd_idx_row[j].coeff(0, k);
        int relidx = P_lmd_idx_row[j].coeff(1, k);
        tmpPDA.row(j) += P_lmd_vecs[block](relidx) *
                         tmpDA.row(P_lmd_idx_bl[block][1] + relidx);
      }
    }
    B_upper += coeff1 * tmpPDA;
  }
  SpEOMatrixD B_lower = coeff2 * X;
  SpEOMatrixD* B_lower2 = new SpEOMatrixD[N_g];
  for (k = 0; k < N_g; k++) {
    B_lower2[k] = coeff3 * Y[k];
    if (fix_Alpha) {
      B_lower2[k] -= coeff3 * (DL[k] * A[k]);
    }
  }

  // compute r_resp = At*b (respective value for diff)
  double r_resp = 0.0;
  r_resp +=
      pow((coeff1 * B_upper + coeff2 * (R.transpose() * B_lower)).norm(), 2);
  if (!fix_Alpha) {
    for (k = 0; k < N_g; k++) {
      SpEOMatrixD temp = SpEOMatrixD::Zero(N_h, A[k].cols());
      for (j = 0; j < A[k].cols(); j++) {  // write temp = B^t*P(k)
        temp.col(j) = (P_lmd_vecs[k](j)) *
                      B_upper.row(P_lmd_idx_bl[k][0] + j).transpose();
      }
      r_resp += pow((-coeff1 * (DH[k].transpose() * temp) +
                     coeff3 * (DL[k].transpose() * B_lower2[k]))
                        .norm(),
                    2);
    }
  }
  if (!fix_delta_m) {
    r_resp += pow((-coeff1 * B_upper.rowwise().sum()).norm(), 2);
  }
  r_resp = sqrt(r_resp);

  // compute initial values for gradient method
  SpEOMatrixD x_upper = SpEOMatrixD::Zero(N_Y, N_h);  // alias Z^p
  if (use_starting_value) {
    x_upper = Z0;
  }
  SpEOMatrixD* x_center = new SpEOMatrixD[N_g];  // alias A_i(Omega_i)'s
  for (k = 0; k < N_g; k++) {
    x_center[k] = SpEOMatrixD::Zero(A[k].rows(), A[k].cols());
    if (use_starting_value) {
      x_center[k] = A[k];
    }
  }
  SpEOVectorD x_lower = SpEOVectorD::Zero(N_Y);  // alias m

  SpEOMatrixD d_upper;
  SpEOMatrixD d_lower;
  SpEOMatrixD* d_lower2 = new SpEOMatrixD[N_g];
  SpEOMatrixD r_upper;
  SpEOMatrixD* r_center = new SpEOMatrixD[N_g];
  SpEOVectorD r_lower;

  // compute stepsize
  double stepsize = 0.0;
  double frob_norm_op = 0.0;
  frob_norm_op += coeff1 * x_upper.rows() * x_upper.cols();
  if (!fix_Alpha) {
    double tmpnorm = 0.0;
    for (k = 0; k < N_g; k++) {
      tmpnorm += pow(A[k].norm() * DH[k].norm() * P_lmd_vecs[k].norm(), 2);
    }
    frob_norm_op += coeff1 * tmpnorm;
  }
  if (!fix_delta_m) {
    frob_norm_op += coeff1 * N_h * N_Y;
  }
  frob_norm_op += coeff2 * pow(R.norm(), 2) * x_upper.cols();
  if (!fix_Alpha) {
    for (k = 0; k < N_g; k++) {
      frob_norm_op += coeff3 * pow(DL[k].norm(), 2);
    }
  }
  stepsize = 2e12 / frob_norm_op;

  // loop
  int iter = 0;
  double diff = (tol_r + 1.0) * r_resp;
  while (iter < maxiter && diff / r_resp > tol_r) {
    // compute d=Ax-b
    d_upper = coeff1 * (x_upper)-B_upper;
    if (!fix_Alpha) {
      SpEOMatrixD tmpDA = SpEOMatrixD::Zero(sum_N_c, N_h);
      int running_idx = 0;
      for (k = 0; k < N_g; k++) {
        tmpDA.block(running_idx, 0, A[k].cols(), N_h) =
            (DH[k] * x_center[k]).transpose();
        running_idx += A[k].cols();
      }
      SpEOMatrixD tmpPDA = SpEOMatrixD::Zero(N_Y, N_h);
      for (j = 0; j < N_Y; j++) {
        for (k = 0; k < P_lmd_idx_row[j].cols(); k++) {
          int block = P_lmd_idx_row[j].coeff(0, k);
          int relidx = P_lmd_idx_row[j].coeff(1, k);
          tmpPDA.row(j) += P_lmd_vecs[block](relidx) *
                           tmpDA.row(P_lmd_idx_bl[block][1] + relidx);
        }
      }
      d_upper += coeff1 * (-tmpPDA);
    }
    if (!fix_delta_m) {
      d_upper += coeff1 * (-x_lower * ones.transpose());
    }
    d_lower = coeff2 * (R * x_upper) - B_lower;
    if (!fix_Alpha) {
      for (k = 0; k < N_g; k++) {
        d_lower2[k] = coeff3 * (DL[k] * x_center[k]) - B_lower2[k];
      }
    }

    // compute r = At*d
    r_upper = coeff1 * d_upper + coeff2 * (R.transpose() * d_lower);
    if (!fix_Alpha) {
      for (k = 0; k < N_g; k++) {
        SpEOMatrixD temp = SpEOMatrixD::Zero(N_h, A[k].cols());
        for (j = 0; j < A[k].cols(); j++) {  // write temp = B^t*P(k)
          temp.col(j) = (P_lmd_vecs[k](j)) *
                        d_upper.row(P_lmd_idx_bl[k][0] + j).transpose();
        }
        r_center[k] = -coeff1 * (DH[k].transpose() * temp) +
                      coeff3 * (DL[k].transpose() * d_lower2[k]);
      }
    }
    if (!fix_delta_m) {
      r_lower = -coeff1 * d_upper.rowwise().sum();
    }

    // compute x_new = x_old - r, diff = norm(r);
    diff = 0.0;
    x_upper -= stepsize * r_upper;
    diff += pow(r_upper.norm(), 2);
    for (k = 0; k < N_g; k++) {
      x_center[k] -= stepsize * r_center[k];
      diff += pow(r_center[k].norm(), 2);
    }
    x_lower -= stepsize * r_lower;
    diff += pow(r_lower.norm(), 2);
    diff = sqrt(diff);

    iter++;
  }

  Z = x_upper;
  A_out = new SpEOMatrixD[N_g];
  for (k = 0; k < N_g; k++) {
    A_out[k] = x_center[k];
  }
  m_out = x_lower;

  cout << "The solution of equation (3) required " << iter
       << " CGLS-iterations. The final relative residual of CGLS: "
       << diff / r_resp << endl;

  return;
}

void calcOptMeanDiffLS(int& iter, double& rel_res, SpEOVectorD& delta_m,
                       SpEOVectorD* delta_m_0, SpEOMatrixD* R, SpEOVectorD* m,
                       SpEOMatrixD* Z, SpEOMatrixD* X, double lambda_X,
                       double lambda_Z, int maxiter, double tol_r) {
  bool use_starting_value = true;

  // compute global stuff
  int N_Y = Z->rows();
  int N_X = X->rows();
  int N_h = Z->cols();
  double coeff1 = sqrt(lambda_Z / N_Y);
  double coeff2 = sqrt(lambda_X / N_X);

  // define operator
  Eq6Op op = Eq6Op(R, coeff1, coeff2);

  // compute initial right-hand side
  Eq6MVar b = Eq6MVar();
  Eq6MVar d = Eq6MVar();
  b.upper = coeff1 * ((1.0 / N_h) * Z->rowwise().sum() - (*m));
  b.lower = coeff2 * ((1.0 / N_h) * X->rowwise().sum() - (*R) * (*m));
  // d = right-hand side
  d.upper = b.upper;
  d.lower = b.lower;

  // compute initial values for CGLS
  Eq6SVar x = Eq6SVar();
  if (use_starting_value) {
    // x = delta_m_0
    x.var = *delta_m_0;
  } else {
    // x = 0
    x.var = SpEOVectorD::Zero(N_Y);
  }

  Eq6SVar r = Eq6SVar();
  // compute norm(operator.transposed()*d);
  r.SetAdjointOperatorXvar(op, &d);
  double r0_norm = sqrt(r.NormSquared());

  if (use_starting_value) {
    // d =  d - operator * x
    d.SubtractOperatorXvar(op, delta_m_0);
  }

  // compute r = operator.transposed()*d;
  r.SetAdjointOperatorXvar(op, &d);

  // p = r
  Eq6SVar p = Eq6SVar();
  p.var = r.var;

  // t = operator * p;
  Eq6MVar t = Eq6MVar();
  t.upper = SpEOVectorD::Zero(N_Y);
  t.lower = SpEOVectorD::Zero(N_X);
  t.SetOperatorXvar(op, &p);

  // CG loop
  double alpha, beta;
  double r_norm_squared = r.NormSquared();

  iter = 0;
  while (iter < maxiter && sqrt(r_norm_squared) / r0_norm > tol_r) {
    // alpha = (norm(r)/norm(t))^2
    alpha = r_norm_squared / t.NormSquared();

    // x = x + alpha * p ; d = d - alpha * t
    x.var += alpha * p.var;

    d.upper -= alpha * t.upper;
    d.lower -= alpha * t.lower;

    // r = operator^T * d
    r.SetAdjointOperatorXvar(op, &d);

    // beta = (norm(r)/norm(r_old))^2
    double r_norm_squared_old = r_norm_squared;
    r_norm_squared = r.NormSquared();
    beta = r_norm_squared / r_norm_squared_old;

    // p = r + beta * p
    p.var = r.var + beta * p.var;

    // t = operator * p
    t.SetOperatorXvar(op, &p);

    iter++;
  }

  delta_m = x.var;

  rel_res = sqrt(r_norm_squared) / r0_norm;
  return;
}

void calcOptCoeffResLS(int& iter, double& rel_res, SpEOMatrixD*& A_res,
                       SpEOMatrixD* A_res_0, SpEOMatrixD* DH, SpEOMatrixD* DL,
                       SpEOVectorD* m, SpEOVectorD* delta_m, SpEOMatrixD* Z,
                       SpEOMatrixD* Y, SpEOMatrixD* X, SpEOMatrixD* R,
                       SpEOVectorD* P_lmd_vecs, SpEOMatrixI* P_lmd_idx_row,
                       int** P_lmd_idx_bl, double lambda_X, double lambda_Y,
                       double lambda_Z, int N_g, int maxiter, double tol_r) {
  bool use_starting_value = true;

  int g, j, k;

  // compute global stuff
  int N_l = Y->cols();
  int N_h = X->cols();
  int N_Y = Y->rows();
  int N_X = X->rows();
  int sum_N_c = 0;
  for (g = 0; g < N_g; g++) {
    sum_N_c += P_lmd_vecs[g].size();
  }

  double* coeff1 = new double[N_g];
  double* coeff2 = new double[N_g];
  for (g = 0; g < N_g; g++) {
    coeff1[g] = sqrt(lambda_Z / (2 * N_h * N_g * A_res_0[g].cols()));
    coeff2[g] = sqrt(lambda_Y / (2 * N_l * N_g * A_res_0[g].cols()));
  }
  double coeff3 = sqrt(0.5 * lambda_X / (N_h * N_X));

  SpEOVectorD onesNh = SpEOVectorD::Ones(N_h);
  SpEOVectorD onesNl = SpEOVectorD::Ones(N_l);

  // define operator
  Eq9Op op = Eq9Op(DH, DL, R, P_lmd_vecs, P_lmd_idx_row, P_lmd_idx_bl, coeff1,
                   coeff2, coeff3);

  // compute initial right-hand side
  Eq9MVar b = Eq9MVar(N_g);
  Eq9MVar d = Eq9MVar(N_g);
  SpEOMatrixD tmp_upper =
      ((*Z) - ((*m) + (*delta_m)) * onesNh.transpose()).transpose();
  SpEOMatrixD tmp_center = ((*Y) - ((*m) * onesNl.transpose())).transpose();
  for (g = 0; g < N_g; g++) {
    b.upper[g] = coeff1[g] *
                 (tmp_upper.middleCols(P_lmd_idx_bl[g][0], A_res_0[g].cols()) -
                  DH[g].col(0) * A_res_0[g].row(0));
    b.center[g] = coeff2[g] * (tmp_center.middleCols(P_lmd_idx_bl[g][0],
                                                     A_res_0[g].cols()) -
                               DL[g].col(0) * A_res_0[g].row(0));
    // d = right-hand side
    d.upper[g] = b.upper[g];
    d.center[g] = b.center[g];
  }
  SpEOMatrixD tmpDA = SpEOMatrixD::Zero(sum_N_c, N_h);
  int running_idx = 0;
  for (g = 0; g < N_g; g++) {
    tmpDA.block(running_idx, 0, A_res_0[g].cols(), N_h) =
        (DH[g].col(0) * A_res_0[g].row(0)).transpose();
    running_idx += A_res_0[g].cols();
  }
  SpEOMatrixD tmpPDA = SpEOMatrixD::Zero(N_Y, N_h);
  for (j = 0; j < N_Y; j++) {
    for (k = 0; k < P_lmd_idx_row[j].cols(); k++) {
      int block = P_lmd_idx_row[j].coeff(0, k);
      int relidx = P_lmd_idx_row[j].coeff(1, k);
      tmpPDA.row(j) += P_lmd_vecs[block](relidx) *
                       tmpDA.row(P_lmd_idx_bl[block][1] + relidx);
    }
  }
  b.lower = coeff3 *
            ((*X) - (*R) * (tmpPDA + ((*m) + (*delta_m)) * onesNh.transpose()));
  d.lower = b.lower;

  // compute initial values for CGLS
  Eq9SVar x = Eq9SVar(N_g);
  for (g = 0; g < N_g; g++) {
    if (use_starting_value) {
      // x = I_Z_tilde // = I_Z_0
      x.var[g] = A_res_0[g].bottomRows(A_res_0[g].rows() - 1);
    } else {
      // x = 0
      x.var[g] = SpEOMatrixD::Zero(A_res_0[g].rows() - 1, A_res_0[g].cols());
    }
  }

  Eq9SVar r = Eq9SVar(N_g);
  // compute norm(operator.transposed()*d);
  r.SetAdjointOperatorXvar(op, &d);
  double r0_norm = sqrt(r.NormSquared());

  if (use_starting_value) {
    // d =  d - operator * x
    d.SubtractOperatorXvar(op, A_res_0);
  }

  // compute r = operator.transposed()*d;
  r.SetAdjointOperatorXvar(op, &d);

  // p = r
  Eq9SVar p = Eq9SVar(N_Y);
  for (g = 0; g < N_g; g++) {
    p.var[g] = r.var[g];
  }

  // t = operator * p;
  Eq9MVar t = Eq9MVar(N_g);
  for (g = 0; g < N_g; g++) {
    t.upper[g] = SpEOMatrixD::Zero(N_h, A_res_0[g].cols());
    t.center[g] = SpEOMatrixD::Zero(N_l, A_res_0[g].cols());
  }
  t.lower = SpEOMatrixD::Zero(N_X, N_h);
  t.SetOperatorXvar(op, &p);

  // CG loop
  double alpha, beta;
  double r_norm_squared = r.NormSquared();

  iter = 0;
  while (iter < maxiter && sqrt(r_norm_squared) / r0_norm > tol_r) {
    // alpha = (norm(r)/norm(t))^2
    alpha = r_norm_squared / t.NormSquared();

    // x = x + alpha * p ; d = d - alpha * t
    for (g = 0; g < N_g; g++) {
      x.var[g] += alpha * p.var[g];

      d.upper[g] -= alpha * t.upper[g];
      d.center[g] -= alpha * t.center[g];
    }
    d.lower -= alpha * t.lower;

    // r = operator^T * d
    r.SetAdjointOperatorXvar(op, &d);

    // beta = (norm(r)/norm(r_old))^2
    double r_norm_squared_old = r_norm_squared;
    r_norm_squared = r.NormSquared();
    beta = r_norm_squared / r_norm_squared_old;

    // p = r + beta * p
    for (g = 0; g < N_g; g++) {
      p.var[g] = r.var[g] + beta * p.var[g];
    }

    // t = operator * p
    t.SetOperatorXvar(op, &p);

    iter++;
  }

  for (g = 0; g < N_g; g++) {
    A_res[g] = SpEOMatrixD::Zero(A_res_0[g].rows(), A_res_0[g].cols());
    A_res[g].topRows(1) = A_res_0[g].topRows(1);
    A_res[g].bottomRows(A_res_0[g].rows() - 1) = x.var[g];
  }
  rel_res = sqrt(r_norm_squared) / r0_norm;

  return;
}

void calcOptCoeffCurrPatchLS(int& iter, double& rel_res, SpEOVectorD*& a,
                             SpEOVectorD* a0, SpEOVectorD* dH, SpEOVectorD* dL,
                             SpEOVectorD* m, SpEOVectorD* delta_m,
                             SpEOMatrixD* Z, SpEOMatrixD* Y, SpEOMatrixD* X,
                             SpEOMatrixD* R, SpEOVectorD* P_lmd_vecs,
                             SpEOMatrixI* P_lmd_idx_row, int** P_lmd_idx_bl,
                             double lambda_X, double lambda_Y, double lambda_Z,
                             int N_g, int maxiter, double tol_r) {
  bool use_starting_value = true;

  int g, j, k;

  // compute global stuff
  int N_l = Y->cols();
  int N_h = X->cols();
  int N_Y = Y->rows();
  int N_X = X->rows();
  int sum_N_c = 0;
  for (g = 0; g < N_g; g++) {
    sum_N_c += P_lmd_vecs[g].size();
  }

  double* coeff1 = new double[N_g];
  double* coeff2 = new double[N_g];
  for (g = 0; g < N_g; g++) {
    coeff1[g] = sqrt(lambda_Z / (2 * N_h * N_g * a0[g].size()));
    coeff2[g] = sqrt(lambda_Y / (2 * N_l * N_g * a0[g].size()));
  }
  double coeff3 = sqrt(0.5 * lambda_X / (N_h * N_X));

  SpEOVectorD onesNh = SpEOVectorD::Ones(N_h);
  SpEOVectorD onesNl = SpEOVectorD::Ones(N_l);

  // define operator
  Eq7Op op = Eq7Op(dH, dL, R, P_lmd_vecs, P_lmd_idx_row, P_lmd_idx_bl, coeff1,
                   coeff2, coeff3);

  // compute initial right-hand side
  Eq7MVar b = Eq7MVar(N_g);
  Eq7MVar d = Eq7MVar(N_g);
  SpEOMatrixD tmp_upper =
      ((*Z) - ((*m) + (*delta_m)) * onesNh.transpose()).transpose();
  SpEOMatrixD tmp_center = ((*Y) - ((*m) * onesNl.transpose())).transpose();
  for (g = 0; g < N_g; g++) {
    b.upper[g] =
        coeff1[g] * tmp_upper.middleCols(P_lmd_idx_bl[g][0], a0[g].size());
    b.center[g] =
        coeff2[g] * tmp_center.middleCols(P_lmd_idx_bl[g][0], a0[g].size());
    // d = right-hand side
    d.upper[g] = b.upper[g];
    d.center[g] = b.center[g];
  }
  b.lower = coeff3 * ((*X) - (*R) * (((*m) + (*delta_m)) * onesNh.transpose()));
  d.lower = b.lower;

  // compute initial values for CGLS
  Eq7SVar x = Eq7SVar(N_g);
  for (g = 0; g < N_g; g++) {
    if (use_starting_value) {
      // x = I_Z_tilde//= I_Z_0
      x.var[g] = a0[g];
    } else {
      // x = 0
      x.var[g] = SpEOVectorD::Zero(a0[g].size());
    }
  }

  Eq7SVar r = Eq7SVar(N_g);
  // compute norm(operator.transposed()*d);
  r.SetAdjointOperatorXvar(op, &d);
  double r0_norm = sqrt(r.NormSquared());

  if (use_starting_value) {
    // d =  d - operator * x
    d.SubtractOperatorXvar(op, a0);
  }
  // compute r = operator.transposed()*d;
  r.SetAdjointOperatorXvar(op, &d);

  // p = r
  Eq7SVar p = Eq7SVar(N_Y);
  for (g = 0; g < N_g; g++) {
    p.var[g] = r.var[g];
  }

  // t = operator * p;
  Eq7MVar t = Eq7MVar(N_g);
  for (g = 0; g < N_g; g++) {
    t.upper[g] = SpEOMatrixD::Zero(N_h, a0[g].size());
    t.center[g] = SpEOMatrixD::Zero(N_l, a0[g].size());
  }
  t.lower = SpEOMatrixD::Zero(N_X, N_h);
  t.SetOperatorXvar(op, &p);

  // CG loop
  double alpha, beta;
  double r_norm_squared = r.NormSquared();

  iter = 0;
  while (iter < maxiter && sqrt(r_norm_squared) / r0_norm > tol_r) {
    // alpha = (norm(r)/norm(t))^2
    alpha = r_norm_squared / t.NormSquared();

    // x = x + alpha * p ; d = d - alpha * t
    for (g = 0; g < N_g; g++) {
      x.var[g] += alpha * p.var[g];

      d.upper[g] -= alpha * t.upper[g];
      d.center[g] -= alpha * t.center[g];
    }
    d.lower -= alpha * t.lower;

    // r = operator^T * d
    r.SetAdjointOperatorXvar(op, &d);

    // beta = (norm(r)/norm(r_old))^2
    double r_norm_squared_old = r_norm_squared;
    r_norm_squared = r.NormSquared();
    beta = r_norm_squared / r_norm_squared_old;

    // p = r + beta * p
    for (g = 0; g < N_g; g++) {
      p.var[g] = r.var[g] + beta * p.var[g];
    }

    // t = operator * p
    t.SetOperatorXvar(op, &p);

    iter++;
  }

  for (g = 0; g < N_g; g++) {
    a[g] = x.var[g];
  }
  rel_res = sqrt(r_norm_squared) / r0_norm;
  return;
}

void calcOptCoeffResPFISTA(int*& iter, double*& rel_res, SpEOMatrixD*& A,
                           SpEOMatrixD* A_0, SpEOMatrixD* DH, SpEOMatrixD* DL,
                           SpEOVectorD* m, SpEOVectorD* delta_m, SpEOMatrixD* Z,
                           SpEOMatrixD* Y, int** P_lmd_idx_bl, double lambda_A,
                           double lambda_Y, double lambda_Z, int N_g,
                           FBSsolverOptions opts, bool spectral_normalizer,
                           bool& write_testset, int testnr) {
  int g;
  int N_h = DH[0].rows();
  int N_l = DL[0].rows();
  int N_DP = DH[0].cols();

  SpEOVectorD onesNh = SpEOVectorD::Ones(N_h);
  SpEOVectorD onesNl = SpEOVectorD::Ones(N_l);

  SpEOMatrixD tmp_upper =
      ((*Z) - ((*m) + (*delta_m)) * onesNh.transpose()).transpose();
  SpEOMatrixD tmp_center = ((*Y) - ((*m)) * onesNl.transpose()).transpose();
  for (g = 0; g < N_g; g++) {
    double N_c = A_0[g].cols();
    double coeff0 = lambda_A / sqrt(N_c);
    double coeff1 = sqrt(
        lambda_Z /
        (N_h * N_c));  // exclude 1/2 in order to fulfill standard formulation
    double coeff2 = sqrt(
        lambda_Y /
        (N_l * N_c));  // exclude 1/2 in order to fulfill standard formulation
    SpEOMatrixD y = SpEOMatrixD::Zero(N_h + N_l, N_c);
    y.topRows(N_h) = coeff1 * (tmp_upper.middleCols(P_lmd_idx_bl[g][0], N_c) -
                               DH[g].col(0) * A_0[g].row(0));
    y.bottomRows(N_l) =
        coeff2 * (tmp_center.middleCols(P_lmd_idx_bl[g][0], N_c) -
                  DL[g].col(0) * A_0[g].row(0));
    SpEOMatrixD D = SpEOMatrixD::Zero(N_h + N_l, N_DP - 1);
    D.topRows(N_h) = coeff1 * DH[g].rightCols(N_DP - 1);
    D.bottomRows(N_l) = coeff2 * DL[g].rightCols(N_DP - 1);
    double normalizer = 1.0;
    if (spectral_normalizer) {  // compute spectral norm
      normalizer =
          1.5 * spec_norm(D);  // use the factor 1.5 in order to make sure that
                               // the normalizer is large enough, since we only
                               // compute an approximation of the spectral norm
    } else {                   // compute frobenius norm
      normalizer = D.norm();
    }
    double dictFac = 1.0 / normalizer;
    SpEOMatrixD A_out;

    // write testset for solver
    string dirout = "/gpfs/work/pr45ne/di72kaq/solvertestsets/datasetJSpFIHMH";
    SpEOMatrixD timestat;
    SpEOMatrixD btstat;
    FBSSolver(iter[g], rel_res[g], A_out, dictFac * D, dictFac * y,
              pow(dictFac, 2) * coeff0, opts, timestat, btstat);
    A[g].bottomRows(A_0[g].rows() - 1) = A_out;
    A[g].row(0) = A_0[g].row(0);
  }
  return;
}

SpEOMatrixD optMeanDiffLS_inverse(SpEOMatrixD* R, double lambda_X,
                                  double lambda_Z) {
  int N_X = R->rows();
  int N_Y = R->cols();

  return ((lambda_Z / N_Y) * SpEOMatrixD::Identity(N_Y, N_Y) +
          (lambda_X / N_X) * R->transpose() * (*R))
      .inverse();
}

void calcOptMeanDiffLS_explicit(SpEOVectorD& delta_m, SpEOMatrixD* Inv,
                                SpEOMatrixD* R, SpEOVectorD* m, SpEOMatrixD* Z,
                                SpEOMatrixD* X, double lambda_X,
                                double lambda_Z) {
  int N_X = R->rows();
  int N_Y = R->cols();
  int N_h = Z->cols();
  SpEOVectorD RHS =
      (lambda_Z / N_Y) * ((1.0 / N_h) * Z->rowwise().sum() - (*m)) +
      (lambda_X / N_X) *
          (R->transpose() * ((1.0 / N_h) * X->rowwise().sum() - (*R) * (*m)));
  delta_m = (*Inv) * RHS;
  return;
}
