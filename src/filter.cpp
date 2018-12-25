#include "filter.h"

void create_Gaussian_filter(SpEOMatrixD &filter, int filter_size) {
  int i, j;

  double kappa = 2.7725887 / pow((filter_size + 1) / 2, 2);

  filter = SpEOMatrixD::Zero(filter_size, filter_size);
  if (filter_size % 2 == 0) {
    int radius = filter_size / 2;
    for (i = 0; i < radius; i++) {
      for (j = i; j < radius; j++) {
        filter(radius + i, radius + j) =
            exp(-(pow(i + 0.5, 2) + pow(j + 0.5, 2)) * kappa);
        filter(radius + j, radius + i) = filter(radius + i, radius + j);
        filter(radius - i - 1, radius + j) = filter(radius + i, radius + j);
        filter(radius - j - 1, radius + i) = filter(radius + i, radius + j);
        filter(radius + i, radius - j - 1) = filter(radius + i, radius + j);
        filter(radius + j, radius - i - 1) = filter(radius + i, radius + j);
        filter(radius - i - 1, radius - j - 1) = filter(radius + i, radius + j);
        filter(radius - j - 1, radius - i - 1) = filter(radius + i, radius + j);
      }
    }
  } else {
    int radius = (filter_size - 1) / 2;
    for (i = 0; i < radius + 1; i++) {
      for (j = i; j < radius + 1; j++) {
        filter(radius + i, radius + j) = exp(-(pow(i, 2) + pow(j, 2)) * kappa);
        filter(radius + j, radius + i) = filter(radius + i, radius + j);
        filter(radius - i, radius + j) = filter(radius + i, radius + j);
        filter(radius - j, radius + i) = filter(radius + i, radius + j);
        filter(radius + i, radius - j) = filter(radius + i, radius + j);
        filter(radius + j, radius - i) = filter(radius + i, radius + j);
        filter(radius - i, radius - j) = filter(radius + i, radius + j);
        filter(radius - j, radius - i) = filter(radius + i, radius + j);
      }
    }
  }
  filter = 1.0 / filter.sum() * filter;

  return;
}

void apply_filter(SpEOMatrixD &HRp, SpEOMatrixD filter,
                  int fDS) {  // only square filters
  int i, j;

  int filter_size = filter.rows();
  int N_i = HRp.rows();
  int N_j = HRp.cols();

  SpEOMatrixD HRp_tmp = HRp;

  if (filter_size % 2 == 0) {
    int radius = filter_size / 2;
    int start_index;
    if (fDS % 4 == 0) {
      start_index = -1;
    } else {
      start_index = 0;
    }
    for (i = start_index; i < N_i; i = i + 2) {
      int lm = max(-i, -radius + 1);
      int rm = min(radius, N_i - (i + 1));
      for (j = start_index; j < N_j; j = j + 2) {
        int tm = max(-j, -radius + 1);
        int bm = min(radius, N_j - (j + 1));
        double filteredpatch_sum =
            (HRp_tmp.block(i + lm, j + tm, rm - lm + 1, bm - tm + 1)
                 .cwiseProduct(filter.block(radius - 1 + lm, radius - 1 + tm,
                                            rm - lm + 1, bm - tm + 1)))
                .sum();
        if (i >= 0 && j >= 0) HRp(i, j) = filteredpatch_sum;
        if (i + 1 < N_i && j >= 0) HRp(i + 1, j) = filteredpatch_sum;
        if (j + 1 < N_j && i >= 0) HRp(i, j + 1) = filteredpatch_sum;
        if (i + 1 < N_i && j + 1 < N_j) HRp(i + 1, j + 1) = filteredpatch_sum;
      }
    }
  } else {
    int radius = (filter_size - 1) / 2;
    for (i = 0; i < N_i; i++) {
      int lm = max(-i, -radius);
      int rm = min(radius, N_i - (i + 1));
      for (j = 0; j < N_j; j++) {
        int tm = max(-j, -radius);
        int bm = min(radius, N_j - (j + 1));
        HRp(i, j) = (HRp_tmp.block(i + lm, j + tm, rm - lm + 1, bm - tm + 1)
                         .cwiseProduct(filter.block(radius + lm, radius + tm,
                                                    rm - lm + 1, bm - tm + 1)))
                        .sum();
      }
    }
  }

  return;
}

void apply_filter_adjoint(SpEOMatrixD &HRp, SpEOMatrixD filter,
                          int fDS) {  // only square filters
  int i, j;

  int filter_size = filter.rows();
  int N_i = HRp.rows();
  int N_j = HRp.cols();

  if (filter_size % 2 == 0) {
    SpEOMatrixD HRp_tmp = SpEOMatrixD::Zero(N_i, N_j);

    int radius = filter_size / 2;
    int start_index;
    if (fDS % 4 == 0) {
      start_index = -1;
    } else {
      start_index = 0;
    }
    for (i = start_index; i < N_i; i = i + 2) {
      int lm = max(-i, -radius + 1);
      int rm = min(radius, N_i - (i + 1));
      for (j = start_index; j < N_j; j = j + 2) {
        int tm = max(-j, -radius + 1);
        int bm = min(radius, N_j - (j + 1));
        if (i >= 0 && j >= 0)
          HRp_tmp.block(i + lm, j + tm, rm - lm + 1, bm - tm + 1) +=
              HRp(i, j) * filter.block(radius - 1 + lm, radius - 1 + tm,
                                       rm - lm + 1, bm - tm + 1);
        if (i + 1 < N_i && j >= 0)
          HRp_tmp.block(i + lm, j + tm, rm - lm + 1, bm - tm + 1) +=
              HRp(i + 1, j) * filter.block(radius - 1 + lm, radius - 1 + tm,
                                           rm - lm + 1, bm - tm + 1);
        if (j + 1 < N_j && i >= 0)
          HRp_tmp.block(i + lm, j + tm, rm - lm + 1, bm - tm + 1) +=
              HRp(i, j + 1) * filter.block(radius - 1 + lm, radius - 1 + tm,
                                           rm - lm + 1, bm - tm + 1);
        if (i + 1 < N_i && j + 1 < N_j)
          HRp_tmp.block(i + lm, j + tm, rm - lm + 1, bm - tm + 1) +=
              HRp(i + 1, j + 1) * filter.block(radius - 1 + lm, radius - 1 + tm,
                                               rm - lm + 1, bm - tm + 1);
      }
    }
    HRp = HRp_tmp;
  } else {
    SpEOMatrixD HRp_tmp = HRp;

    int radius = (filter_size - 1) / 2;
    for (i = 0; i < N_i; i++) {
      int lm = max(-i, -radius);
      int rm = min(radius, N_i - (i + 1));
      for (j = 0; j < N_j; j++) {
        int tm = max(-j, -radius);
        int bm = min(radius, N_j - (j + 1));
        HRp(i, j) = (HRp_tmp.block(i + lm, j + tm, rm - lm + 1, bm - tm + 1)
                         .cwiseProduct(filter.block(radius + lm, radius + tm,
                                                    rm - lm + 1, bm - tm + 1)))
                        .sum();
      }
    }
  }

  return;
}

void select(SpEOMatrixD &LRp, SpEOMatrixD HRp, int fDS) {
  int i, j;

  int LR_rows = HRp.rows() / fDS;
  int LR_cols = HRp.cols() / fDS;
  LRp = SpEOMatrixD::Zero(LR_rows, LR_cols);

  for (i = 0; i < LR_rows; i++) {
    for (j = 0; j < LR_cols; j++) {
      LRp(i, j) = HRp(
          i * fDS + fDS / 2,
          j * fDS + fDS / 2);  // lower right corner of center 2x2 square if fDS
                               // even -> ok for Gaussian filter and assumption
                               // that dimensions of HRp are multiple of fDS
    }
  }

  return;
}

void distribute(SpEOMatrixD &HRp, SpEOMatrixD LRp,
                int fDS) {  // adjoint operator of "select"
  int i, j;

  int HR_rows = LRp.rows() * fDS;
  int HR_cols = LRp.cols() * fDS;
  int LR_rows = LRp.rows();
  int LR_cols = LRp.cols();

  HRp = SpEOMatrixD::Zero(HR_rows, HR_cols);

  for (i = 0; i < LR_rows; i++) {
    for (j = 0; j < LR_cols; j++) {
      HRp(i * fDS + fDS / 2, j * fDS + fDS / 2) = LRp(i, j);
    }
  }
  return;
}

void calc_filter_boundary_coeff(SpEOMatrixD &filter_coeff, SpEOMatrixD filter,
                                int fDS) {  // assume filter.sum() == 1
  int i, j;

  int filter_size = filter.rows();
  int table_size = 0;
  if (filter_size % 2 == 0) {
    table_size = filter_size / 2;
    if (fDS % 4 == 0) {
      table_size++;
    }
  } else {
    table_size = (filter_size + 1) / 2;
  }
  filter_coeff = SpEOMatrixD::Zero(table_size, table_size);

  for (i = 0; i < table_size; i++) {
    for (j = i; j < table_size; j++) {
      filter_coeff(i, j) =
          1.0 / filter.block(i, j, filter_size - i, filter_size - j).sum();
      filter_coeff(j, i) = filter_coeff(i, j);
    }
  }
  return;
}

void apply_filter_boundary_coeff(
    SpEOMatrixD &HRp, SpEOMatrixD filter_coeff, int filter_size,
    int fDS) {  // filter_size even <=> image sizes even
  int i, j;

  int HR_rows = HRp.rows();
  int HR_cols = HRp.cols();

  // can be computed more efficiently -> later
  if (filter_size % 2 == 0) {
    int max_overlap = filter_size / 2 - 1;
    int start_index;
    if (fDS % 4 == 0) {
      start_index = -1;
    } else {
      start_index = 0;
    }
    for (i = start_index; i < HR_rows; i = i + 2) {
      int i_overlap =
          max(0, max(max_overlap - i, max_overlap - (HR_rows - (i + 2))));
      for (j = start_index; j < HR_cols; j = j + 2) {
        int j_overlap =
            max(0, max(max_overlap - j, max_overlap - (HR_cols - (j + 2))));
        if (i_overlap > 0 || j_overlap > 0) {
          if (i >= 0 && j >= 0) HRp(i, j) *= filter_coeff(i_overlap, j_overlap);
          if (i + 1 < HR_rows && j >= 0)
            HRp(i + 1, j) *= filter_coeff(i_overlap, j_overlap);
          if (j + 1 < HR_cols && i >= 0)
            HRp(i, j + 1) *= filter_coeff(i_overlap, j_overlap);
          if (i + 1 < HR_rows && j + 1 < HR_cols)
            HRp(i + 1, j + 1) *= filter_coeff(i_overlap, j_overlap);
        }
      }
    }
  } else {
    int max_overlap = (filter_size - 1) / 2;
    for (i = 0; i < HR_rows; i++) {
      int i_overlap =
          max(0, max(max_overlap - i, max_overlap - (HR_rows - (i + 1))));
      for (j = 0; j < HR_cols; j++) {
        int j_overlap =
            max(0, max(max_overlap - j, max_overlap - (HR_cols - (j + 1))));
        if (i_overlap > 0 || j_overlap > 0) {
          HRp(i, j) *= filter_coeff(i_overlap, j_overlap);
        }
      }
    }
  }

  return;
}

void standard_filter(SpEOMatrixD &LRp, SpEOMatrixD HRp, SpEOMatrixD filter,
                     SpEOMatrixD filter_coeff, int fDS) {
  int filter_size = filter.rows();
  SpEOMatrixD tmp_image_HR = HRp;
  apply_filter(tmp_image_HR, filter, fDS);
  apply_filter_boundary_coeff(tmp_image_HR, filter_coeff, filter_size, fDS);
  select(LRp, tmp_image_HR, fDS);
  return;
}

void standard_adjoint_filter(SpEOMatrixD &HRp, SpEOMatrixD LRp,
                             SpEOMatrixD filter, SpEOMatrixD filter_coeff,
                             int fDS) {
  int filter_size = filter.rows();
  distribute(HRp, LRp, fDS);
  apply_filter_boundary_coeff(HRp, filter_coeff, filter_size, fDS);
  apply_filter_adjoint(HRp, filter, fDS);
  return;
}

void fast_filter(SpEOMatrixD &LRp, SpEOMatrixD HRp, SpEOMatrixD filter,
                 SpEOMatrixD filter_coeff, int fDS) {
  int i_LR, i_HR, j_LR, j_HR;

  int filter_size = filter.rows();
  int LR_rows = HRp.rows() / fDS;
  int LR_cols = HRp.cols() / fDS;
  int HR_rows = HRp.rows();
  int HR_cols = HRp.cols();

  LRp = SpEOMatrixD::Zero(LR_rows, LR_cols);
  if (filter_size % 2 == 0) {
    int radius = filter_size / 2;
    int max_overlap = filter_size / 2 - 1;
    for (i_LR = 0; i_LR < LR_rows; i_LR++) {
      i_HR = fDS * i_LR + fDS / 2 - 1;
      int lm = max(-i_HR, -radius + 1);
      int rm = min(radius, HR_rows - (i_HR + 1));
      int i_overlap =
          max(0, max(max_overlap - i_HR, max_overlap - (HR_rows - (i_HR + 2))));
      for (j_LR = 0; j_LR < LR_cols; j_LR++) {
        j_HR = fDS * j_LR + fDS / 2 - 1;
        int tm = max(-j_HR, -radius + 1);
        int bm = min(radius, HR_cols - (j_HR + 1));
        int j_overlap = max(
            0, max(max_overlap - j_HR, max_overlap - (HR_cols - (j_HR + 2))));
        LRp(i_LR, j_LR) =
            filter_coeff(i_overlap, j_overlap) *
            (HRp.block(i_HR + lm, j_HR + tm, rm - lm + 1, bm - tm + 1)
                 .cwiseProduct(filter.block(radius - 1 + lm, radius - 1 + tm,
                                            rm - lm + 1, bm - tm + 1)))
                .sum();
      }
    }
  } else {
    int radius = (filter_size - 1) / 2;
    int max_overlap = (filter_size - 1) / 2;
    for (i_LR = 0; i_LR < LR_rows; i_LR++) {
      i_HR = fDS * i_LR + fDS / 2;
      int lm = max(-i_HR, -radius);
      int rm = min(radius, HR_rows - (i_HR + 1));
      int i_overlap =
          max(0, max(max_overlap - i_HR, max_overlap - (HR_rows - (i_HR + 1))));
      for (j_LR = 0; j_LR < LR_cols; j_LR++) {
        j_HR = fDS * j_LR + fDS / 2;
        int tm = max(-j_HR, -radius);
        int bm = min(radius, HR_cols - (j_HR + 1));
        int j_overlap = max(
            0, max(max_overlap - j_HR, max_overlap - (HR_cols - (j_HR + 1))));
        LRp(i_LR, j_LR) =
            filter_coeff(i_overlap, j_overlap) *
            (HRp.block(i_HR + lm, j_HR + tm, rm - lm + 1, bm - tm + 1)
                 .cwiseProduct(filter.block(radius + lm, radius + tm,
                                            rm - lm + 1, bm - tm + 1)))
                .sum();
      }
    }
  }
  return;
}