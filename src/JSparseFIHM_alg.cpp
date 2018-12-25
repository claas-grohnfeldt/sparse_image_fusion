/*
 * JSparseFIHM_alg.cpp
 *
 *  Created on: Feb 19, 2014
 *      Author: Claas Grohnfeldt
 */

#include "JSparseFIHM_alg.h"

void JSparseFIHM_alg(SpEOMatrixD &EndmemberMat, SpEOMatrixD *&AbundanceMat,
                     int iterMain, int numIterMain, SpEOPaths *paths,
                     SpEODataIOSetting *dSet, SpEOFusionSetting *fSet,
                     SpEOOutputSetting *oSet, SpEOSolverSetting *sSet,
                     SpEOParallelSetting *pSet, SpEOGlobalParams *glPrms,
                     SpEODataset *ImX, SpEODataset *ImX_LR, SpEODataset *ImY,
                     SpEODataset *ImZ, SpEOMatrixD *SRF, MPI_Comm comm_busy,
                     MPI_Group group_busy, SpEODataset *ImZ_init,
                     SpEOReport &report) {
  int my_rank;
  int my_processes;

  MPI_Comm_rank(comm_busy, &my_rank);
  MPI_Comm_size(comm_busy, &my_processes);

#ifndef _OPENMP
  MPI_Barrier(comm_busy);
  MPI_Comm comm_patch;
  MPI_Group group_patch;
  MPI_Comm_group(comm_busy, &group_patch);

  int ranges[1][3] = {
      {(my_rank / pSet->numProcPerPatch) * pSet->numProcPerPatch,
       (my_rank / pSet->numProcPerPatch + 1) * pSet->numProcPerPatch - 1, 1}};
  MPI_Group_range_incl(group_busy, 1, ranges, &group_patch);
  MPI_Comm_create(comm_busy, group_patch, &comm_patch);
  int my_patch_rank;
  int my_patch_processes;
  MPI_Comm_rank(comm_patch, &my_patch_rank);
  MPI_Comm_size(comm_patch, &my_patch_processes);
#endif
  bool print_CSG_output = false;
  bool print_optimization_input_output = false;
  bool print_other_stuff_during_patch_rec = false;
  bool print_for_debugging = false;

  int testnr = 0;
  bool write_testset = false;
  if (my_rank == 0) {
    cout << "\n"
         << "###########################################################\n"
         << "##                                                       ##\n"
         << "##     algorithm for patch-wise image fusion started     ##\n"
         << "##                                                       ##\n"
         << "###########################################################"
         << "\n"
         << "\n";
  }
  check_for_inf_or_nan(my_rank, (*SRF), " ", -123, "(*SRF)");

  if (my_rank == 0) {
    cout << "     step 1: Declarations and initializations.. ";
  }

  // running variables
  int jP, iP, uP, vP, iG;

  iP = glPrms->numPatchGroups;

  // frequently used variables
  int pszL = fSet->patchsize;
  short fDS = glPrms->fDS;
  int NChX = ImX->get_NCh();
  int NChY = ImY->get_NCh();
  int pszH = pszL * fDS;
  int pszL2 = pszL * pszL;
  int pszH2 = pszH * pszH;
  int a = pszL - fSet->overlap;
  int NP = glPrms->NP;
  int NPU = glPrms->NPU;  // number of patches in a column
  int NPV = glPrms->NPV;  // number of patches in a row

  int NPU_sub = glPrms->NPU_sub;  // uPLast-uPFirst+1;
  int NPV_sub = glPrms->NPV_sub;  // vPLast-vPFirst+1;
  int NP_sub = NPU_sub * NPV_sub;

  int NDP = fSet->NDP;

  int uPFirst = glPrms->uPFirst;
  int vPFirst = glPrms->vPFirst;
  // PARAMETER TO BE SET FOR NORM PROPAGATION -> depending on the test outcome,
  // the Frobenius can be removed
  bool spectral = fSet->matrixNorm;  // frobenius = 0, spectral norm = 1
  // form index matrices for
  // - the link between the total patch index and the horizontal and vertical
  //   patch indices (idxM, idxME)
  // - the upper left corner of all patches (idxPUH, idxPUL, idxPVH, idxPVL)

  int col_idx = 0, row_idx;
  SpEOVectorI idxPUL = SpEOVectorI::Zero(NPU, 1);
  SpEOVectorI idxPVL = SpEOVectorI::Zero(NPV, 1);
  SpEOVectorI idxPUH = glPrms->idxPUH;
  SpEOVectorI idxPVH = glPrms->idxPVH;

  for (uP = 0; uP < NPU - 1; uP++) {
    idxPUL(uP) = a * uP;
  }
  idxPUL(NPU - 1) = ImY->get_sizeU() - pszL;
  for (vP = 0; vP < NPV - 1; vP++) {
    idxPVL(vP) = a * vP;
  }
  idxPVL(NPV - 1) = ImY->get_sizeV() - pszL;
  // dynamic dictionaries for dynamically changing (selected) possibly nearest
  // NDP patches
  SpEOMatrixD zHR = SpEOMatrixD::Zero(pszH2, NChY);
  string iPStr;

  // Dictionary Selection Initialization
  SpEOMatrixF patchComp = SpEOMatrixF::Zero(NP, 3);
  SpEOMatrixD patchCompD = SpEOMatrixD::Zero(NP, 3);

  for (iP = 0; iP < NP; iP++) {
    uP = iP / NPV;
    vP = iP % NPV;
    patchComp(iP, 0) = uP;
    patchComp(iP, 1) = vP;
  }
  patchComp = SpEOMatrixF::Zero(0, 0);
  double mytimeDictSum = 0;
  double mytimeOptSum = 0;

  double mytimeTMP1 = 0;
  double mytimeTMP2 = 0;
  double mytimeTMP3 = 0;
  double mytimeTMP4 = 0;

  SpEOMatrixI incomplPatchList;
  bool *patchExists;
  patchExists = new bool[NP];
  for (iP = 0; iP < NP; iP++) {
    patchExists[iP] = false;
  }

  int numRemPatches = NP_sub;
  int *remainingPatches;
  remainingPatches = new int[1];
  // clean up
  delete[] patchExists;

  //============================================================================//
  //	         Sparse high-resolution image reconstruction                  //
  //	             (local-non-local processing modul)                       //
  //============================================================================//
  if (my_rank == 0) {
    cout << endl << "     step 2: Sparse Reconstruction.. " << endl;
    cout << "         number of patches = NPU x NPV = " << NPU << " x " << NPV
         << " = " << glPrms->NP << endl;
    cout << "     [my_rank] (iP, uP, vP) =" << endl;
  }

  // parallelization of the patches
  int groupN = pSet->numProcGrp;  // number of processors per patch

  // define subcommunicators (maybe already in the beginning)
  MPI_Comm group_comm;
  MPI_Comm_split(comm_busy, my_rank / groupN, my_rank % groupN,
                 &group_comm);  // create group communicators

  // get your own (local) rank
  int my_local_rank;
  MPI_Comm_rank(group_comm, &my_local_rank);

  // work accounting variables
  int numberProcessGroups = my_processes / groupN;
  if (my_processes % groupN > 0) {
    numberProcessGroups++;
  }

#ifndef _OPENMP
  jP = my_rank / pSet->numProcPerPatch;
#else
  jP = my_rank;
#endif
  int pLast_sub = NP_sub - 1;

  // additional work stealing
  int lastFixedPatch = pLast_sub;
  if (pSet->workStealingTurns >= 0) {
    lastFixedPatch -= (pLast_sub + 1) % glPrms->numPatchGroups +
                      pSet->workStealingTurns * glPrms->numPatchGroups;
    lastFixedPatch = max(lastFixedPatch, glPrms->numPatchGroups - 1);
  }

  struct mpi_counter_t *c;
  c = create_counter(0, lastFixedPatch + 1, comm_busy);

  // Compute NP indices for all patches, for use in Non-Local Dictionary
  // Training
  while (jP <= pLast_sub) {  //      as long as there are patches to work on
    uP = uPFirst + jP / NPV_sub;
    vP = vPFirst + jP % NPV_sub;
    iP = uP * NPV + vP;
    if (print_for_debugging && (my_rank == 0))
      cout << "[" << my_rank << "] bp_0" << endl;
#ifndef _OPENMP
    if (my_rank % pSet->numProcPerPatch == 0) {
      if (my_rank == 0) {
        cout << "     [" << my_rank << "] (iP, uP, vP) = (" << jP << " (of "
             << pLast_sub << "), " << jP / NPV_sub << ", " << jP % NPV_sub
             << ")" << endl;  //, and (iP,uP,vP)_total=(" << iP << "," << uP <<
                              //"," << vP << ")" << endl;
        report.file.open(report.fileName.c_str(),
                         fstream::in | fstream::out | fstream::app);

        report.file << "     [" << my_rank << "] (iP, uP, vP) = (" << jP
                    << " (of " << pLast_sub << "), " << jP / NPV_sub << ", "
                    << jP % NPV_sub
                    << ")\n";  //, and (iP,uP,vP)_total=(" << iP << "," << uP <<
                               //"," << vP << ")" << "\n";
        report.file.close();
      }
    }
#else
    // if(!fSet->LQ_post_opt){
    if (my_rank == 0) {
      cout << endl
           << "     "
              "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
              "++++++++++++++++"
           << endl
           << "     [" << my_rank << "] (iP,uP,vP)_local=(" << jP << "(of "
           << pLast_sub << ")," << jP / NPV_sub << "," << jP % NPV_sub
           << "), and (iP,uP,vP)_total=(" << iP << "," << uP << "," << vP << ")"
           << endl
           << "     "
              "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
              "++++++++++++++++"
           << endl
           << endl;
    }
    //}
#endif
    //###########################################################################
    //
    //          new model (from version 1.0 -> 2.0) for calculating the patch
    //          patZ beginning ...
    //
    //###########################################################################
    NDP = fSet->NDP;
    //+++++++++++++++++++++
    //+  extract patches  +
    //+++++++++++++++++++++
    if (print_for_debugging && (my_rank == 0))
      cout << "[" << my_rank << "] bp-01: extract patches" << endl;
    // extract patch patY from ImY and calculate mean values in all channels
    // individually
    SpEOMatrixD patY = SpEOMatrixD::Zero(NChY, pszL2);
    SpEOVectorD m_Y = SpEOVectorD::Zero(NChY);
    int iChY;
    for (iChY = 0; iChY < NChY; iChY++) {
      SpEOMatrixD patchBand =
          (ImY->get_rasterBands()[iChY]->get_bandDataMat()->block(
               idxPUL.coeff(uP), idxPVL.coeff(vP), pszL, pszL))
              .cast<double>();
      patY.row(iChY) = SpEOVectorD::Map(patchBand.data(), pszL2);
      m_Y(iChY) = patchBand.mean();
    }
    check_for_inf_or_nan(my_rank, m_Y, " ", -123, "m_Y");
    // extract patch patX from ImX
    SpEOMatrixD patX = SpEOMatrixD::Zero(NChX, pszH2);
    int iChX;
    for (iChX = 0; iChX < NChX; iChX++) {
      SpEOMatrixD patchBand =
          (ImX->get_rasterBands()[iChX]->get_bandDataMat()->block(
               idxPUH.coeff(uP), idxPVH.coeff(vP), pszH, pszH))
              .cast<double>();
      patX.row(iChX) = SpEOVectorD::Map(patchBand.data(), pszH2);
    }
    check_for_inf_or_nan(my_rank, patX, " ", -123, "patX");
    // extract patch patX_LR from ImX_LR
    SpEOMatrixD patX_LR = SpEOMatrixD::Zero(NChX, pszL2);
    for (iChX = 0; iChX < NChX; iChX++) {
      SpEOMatrixD patchBand =
          (ImX_LR->get_rasterBands()[iChX]->get_bandDataMat()->block(
               idxPUL.coeff(uP), idxPVL.coeff(vP), pszL, pszL))
              .cast<double>();
      patX_LR.row(iChX) = SpEOVectorD::Map(patchBand.data(), pszL2);
    }
    check_for_inf_or_nan(my_rank, patX_LR, " ", -123, "patX_LR");

    //+++++++++++++++++++++
    //+  extract windows  +
    //+++++++++++++++++++++
    if (print_for_debugging && (my_rank == 0))
      cout << "[" << my_rank << "] bp-02: extract windows" << endl;
    int winSizeL = fSet->winSize;  // must have the same sign as pszL in order
                                   // to have both centers matched
    int winSizeH = winSizeL * fDS;

    int idxWUL =
        max(0, idxPUL.coeff(uP) - (int)(0.5 * (double)(winSizeL - pszL)));
    int idxWVL =
        max(0, idxPVL.coeff(vP) - (int)(0.5 * (double)(winSizeL - pszL)));
    int winSizeUL =
        min((int)ImY->get_sizeU() - 1,
            idxPUL.coeff(uP) - (int)(0.5 * (double)(winSizeL - pszL)) +
                winSizeL - 1) -
        idxWUL + 1;
    int winSizeVL =
        min((int)ImY->get_sizeV() - 1,
            idxPVL.coeff(vP) - (int)(0.5 * (double)(winSizeL - pszL)) +
                winSizeL - 1) -
        idxWVL + 1;

    int idxWUH = fDS * idxWUL;
    int idxWVH = fDS * idxWVL;
    int winSizeUH = fDS * winSizeUL;
    int winSizeVH = fDS * winSizeVL;

    // extract window winY from ImY in each channel individually
    if (print_for_debugging && (my_rank == 0))
      cout << "[" << my_rank << "] bp-02-01" << endl;
    SpEOMatrixD winY = SpEOMatrixD::Zero(NChY, winSizeUL * winSizeVL);
    for (iChY = 0; iChY < NChY; iChY++) {
      SpEOMatrixD windowBand =
          (ImY->get_rasterBands()[iChY]->get_bandDataMat()->block(
               idxWUL, idxWVL, winSizeUL, winSizeVL))
              .cast<double>();
      SpEOVectorD winY_row =
          SpEOVectorD::Map(windowBand.data(), winSizeUL * winSizeVL);
      double winY_row_mean = winY_row.mean();
      winY_row.array() -= winY_row_mean;
      double winY_row_norm = winY_row.norm();
      if (winY_row_norm > 1e-8) {
        winY_row /= winY_row_norm;
      }
      winY.row(iChY) = winY_row;
    }
    check_for_inf_or_nan(my_rank, winY, " ", -123, "winY");
    // extract window winX from ImX
    if (print_for_debugging && (my_rank == 0))
      cout << "[" << my_rank << "] bp-02-02" << endl;
    SpEOMatrixD winX = SpEOMatrixD::Zero(NChX, winSizeUH * winSizeVH);
    for (iChX = 0; iChX < NChX; iChX++) {
      SpEOMatrixD windowBand =
          (ImX->get_rasterBands()[iChX]->get_bandDataMat()->block(
               idxWUH, idxWVH, winSizeUH, winSizeVH))
              .cast<double>();
      winX.row(iChX) =
          SpEOVectorD::Map(windowBand.data(), winSizeUH * winSizeVH);
      winX.row(iChX).array() -= winX.row(iChX).mean();
      double winX_row_norm = winX.row(iChX).norm();
      if (winX_row_norm > 1e-8) {
        winX.row(iChX) /= winX_row_norm;
      }
    }
    check_for_inf_or_nan(my_rank, winX, " ", -123, "winX");
    // extract window winX_LR from ImX_LR
    if (print_for_debugging && (my_rank == 0))
      cout << "[" << my_rank << "] bp-02-03" << endl;
    SpEOMatrixD winX_LR = SpEOMatrixD::Zero(NChX, winSizeUL * winSizeVL);
    for (iChX = 0; iChX < NChX; iChX++) {
      SpEOMatrixD windowBand =
          (ImX_LR->get_rasterBands()[iChX]->get_bandDataMat()->block(
               idxWUL, idxWVL, winSizeUL, winSizeVL))
              .cast<double>();
      winX_LR.row(iChX) =
          SpEOVectorD::Map(windowBand.data(), winSizeUL * winSizeVL);
      winX_LR.row(iChX).array() -= winX_LR.row(iChX).mean();
      double winX_LR_row_norm = winX_LR.row(iChX).norm();
      if (winX_LR_row_norm > 1e-8) {
        winX_LR.row(iChX) /= winX_LR_row_norm;
      }
    }
    check_for_inf_or_nan(my_rank, winX_LR, " ", -123, "winX_LR");
    //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++//
    //+              correlation based spectral grouping            +//
    //+                          - START -                          +//
    //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++//
    if (print_for_debugging && (my_rank == 0))
      cout << "[" << my_rank << "] bp-03: Correlation based spectral grouping"
           << endl;

    //++++++++++++++++++++++++++++++++++++++++++++++
    //+   calculate Ng, Nc_vec and idxChY via      +
    //+   correlation based spectral grouping      +
    //++++++++++++++++++++++++++++++++++++++++++++++
    int Ng;
    int Nc_vec[NChY];
    int idxChY[NChY];

    CSG_corr_based_spectral_grouping(Ng, idxChY, Nc_vec, winY, fSet->theta,
                                     fSet->Nc_max, my_rank);
    check_for_inf_or_nan(my_rank, Ng, " ", -123, "Ng");
    for (int ig = 0; ig < Ng; ig++) {
      check_for_inf_or_nan(my_rank, idxChY[ig], "ig=", ig, ": idxChY[ig]");
      check_for_inf_or_nan(my_rank, Nc_vec[ig], "ig=", ig, ": Nc_vec[ig]");
    }

    if (print_CSG_output && my_rank == 0) {
      cout << "         CGS outcome:" << endl
           << "             Ng=" << Ng << endl;
      int ig;
      for (ig = 0; ig < Ng - 1; ig++) {
        cout << "             (idxChY[ig=" << ig << "], Nc_vec[ig=" << ig
             << "], o_fw) = (" << idxChY[ig] << ", " << Nc_vec[ig] << ", "
             << idxChY[ig] + Nc_vec[ig] - idxChY[ig + 1] << ")" << endl;
      }
      ig = Ng - 1;
      cout << "             (idxChY[ig=" << ig << "], Nc_vec[ig=" << ig
           << "], o_fw) = (" << idxChY[ig] << ", " << Nc_vec[ig] << ", " << 0
           << ")" << endl;
    }

    //++++++++++++++++++++++++++++++++++++++++++++++
    //+                                            +
    //+  calculate sum_Nc_vec, P_lmd_vecs_loc,     +
    //+  P_lmd_idx_bl, P_lmd_idx_row, avrgBnds     +
    //+                                            +
    //++++++++++++++++++++++++++++++++++++++++++++++
    if (print_for_debugging && (my_rank == 0))
      cout << "[" << my_rank << "] bp-03-01" << endl;
    // initialize the vectors (containing the non-trivial entries along the
    // diagonal in the corresponding blocks) in P_lmd
    SpEOVectorD *P_lmd_vecs_loc = new SpEOVectorD[Ng];
    int **P_lmd_idx_bl_loc = new int *[Ng];
    for (int ig = 0; ig < Ng; ig++) {
      P_lmd_idx_bl_loc[ig] = new int[2];
    }
    SpEOMatrixI *
        P_lmd_idx_row_loc;  // for every row (each corresponding to one HS
                            // channel iChY) these matrices contain the relevant
                            // corresponding block indexes (first row) and the
                            // corresponding entry indexes relative to the
                            // corresponding block's origin each starting at 0.
    P_lmd_idx_row_loc = new SpEOMatrixI[glPrms->NChZ];

    calc_P_matrices(P_lmd_vecs_loc, P_lmd_idx_bl_loc, P_lmd_idx_row_loc, Ng,
                    Nc_vec, glPrms->NChZ, idxChY, my_rank);

    for (int ig = 0; ig < Ng; ig++) {
      check_for_inf_or_nan(my_rank, P_lmd_vecs_loc[ig], "ig=", ig,
                           ": P_lmd_vecs_loc[ig]");
      check_for_inf_or_nan(my_rank, P_lmd_idx_row_loc[ig], "ig=", ig,
                           ": P_lmd_idx_row_loc[ig]");
      check_for_inf_or_nan(my_rank, P_lmd_idx_bl_loc[ig][0], "ig=", ig,
                           ": P_lmd_idx_bl_loc[ig][0]");
      check_for_inf_or_nan(my_rank, P_lmd_idx_bl_loc[ig][1], "ig=", ig,
                           ": P_lmd_idx_bl_loc[ig][1]");
    }

    //++++++++++++++++++++++++++++++++++++++++++++++
    //+                                            +
    //+      calculate ImX_sim                     +
    //+                                            +
    //++++++++++++++++++++++++++++++++++++++++++++++
    if (print_for_debugging && (my_rank == 0))
      cout << "[" << my_rank << "] bp-04: Calculate ImX_sim" << endl;
    //==================================================================//
    //      Calculate a preliminary (temporary) version of ImX_sim      //
    //==================================================================//
    if (print_for_debugging && (my_rank == 0))
      cout << "[" << my_rank
           << "] bp-04-01: Calculate a preliminary (temporary) version of "
              "ImX_sim"
           << endl;
    SpEODataset *ImX_sim;  //, *ImX_sim_LR;
    ImX_sim = new SpEODataset(HR, imFlag_X_sim);

    if (print_other_stuff_during_patch_rec && my_rank == 0) {
      cout << "         calc_ImX_sim .." << endl;
    }
    int sim_mode = fSet->ImX_sim_mode;
    calc_ImX_sim(ImX_sim, ImX, ImX_LR, ImY, fSet, dSet, glPrms, &winX_LR, &winY,
                 &idxPUL, &idxPVL, uP, vP, true, Ng, Nc_vec, idxChY, sim_mode,
                 comm_busy);

    for (int ig = 0; ig < Ng; ig++) {
      check_for_inf_or_nan(
          my_rank, ImX_sim->get_rasterBands()[ig]->bandDataMatD, "ig=", ig,
          ": ImX_sim->get_rasterBands()[ig]->bandDataMatD");
    }
    //===================================================================//
    //  Low-pass filter and down-sample ImX_sim to generate ImX_sim_LR   //
    //===================================================================//
    if (print_for_debugging && (my_rank == 0))
      cout << "[" << my_rank
           << "] bp-04-02: Low-pass filter and down-sample ImX_sim to generate "
              "ImX_sim_LR"
           << endl;
    if (print_other_stuff_during_patch_rec && my_rank == 0) {
      cout << "         Declare and initialize ImX_sim_LR dataset .." << endl;
    }
    SpEODataset *ImX_sim_LR;
    ImX_sim_LR = new SpEODataset(LR, imFlag_X_sim_LR);
    ImX_sim_LR->copyMetaInfoFromDatasets(ImY, ImX_sim, ImY, SpEODouble);
    if (print_other_stuff_during_patch_rec && my_rank == 0) {
      cout << "         Calculate ImX_sim_LR by low pass filtering and "
              "down-sampling ImX_sim ..";
    }
    lowPassFilter_and_downSample(ImX_sim, ImX_sim_LR, SpEODouble, SpEODouble,
                                 *glPrms);
    if (print_other_stuff_during_patch_rec && my_rank == 0) {
      cout << " done!" << endl;
    }
    for (int ig = 0; ig < Ng; ig++) {
      check_for_inf_or_nan(
          my_rank, ImX_sim_LR->get_rasterBands()[ig]->bandDataMatD, "ig=", ig,
          ": ImX_sim_LR->get_rasterBands()[ig]->bandDataMatD");
    }
    //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++//
    //+                                                             +//
    //+              correlation based spectral grouping            +//
    //+                           - END -                           +//
    //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++//
    int iG, iChZ, iC, ig;

    SpEOMatrixD patZ = SpEOMatrixD::Zero(NChY, pszH2);
    // ******************************************************************************
    // generate coupled local dictionaries
    // ******************************************************************************
    if (print_for_debugging && (my_rank == 0))
      cout << "[" << my_rank << "] bp-05: generate coupled local dictionaries"
           << endl;
    SpEOVectorD *alpha = new SpEOVectorD[Ng];
    SpEOVectorD *alpha_init = new SpEOVectorD[Ng];
    SpEOVectorD *dictHR = new SpEOVectorD[Ng];
    SpEOVectorD *dictLR = new SpEOVectorD[Ng];
    for (ig = 0; ig < Ng; ig++) {
      alpha[ig] = SpEOVectorD::Zero(Nc_vec[ig]);
      alpha_init[ig] = SpEOVectorD::Zero(Nc_vec[ig]);
      SpEOMatrixD patchBand_HR =
          (ImX_sim->get_rasterBands()[ig]->bandDataMatD.block(
              idxPUH.coeff(uP), idxPVH.coeff(vP), pszH, pszH));
      SpEOMatrixD patchBand_LR =
          (ImX_sim_LR->get_rasterBands()[ig]->bandDataMatD.block(
              idxPUL.coeff(uP), idxPVL.coeff(vP), pszL, pszL));
      if (fSet->substrMean) {
        patchBand_LR.array() -= patchBand_LR.mean();
      }
      if (fSet->nrmlDicts) {
        double patchBand_norm = patchBand_HR.norm();
        if (fSet->use_LRnorm_for_dic_normalization) {
          patchBand_norm = patchBand_LR.norm();
        }
        if (patchBand_norm < 1e-8) {
          patchBand_HR.block(0, 0, fDS, fDS) = SpEOMatrixD::Ones(fDS, fDS);
          patchBand_LR(0, 0) = 1;
        }
      }
      if (fSet->substrMean) {
        patchBand_HR.array() -= patchBand_HR.mean();
        patchBand_LR.array() -= patchBand_LR.mean();
      }
      if (fSet->nrmlDicts) {
        double norm_tmp = patchBand_HR.norm();
        if (fSet->use_LRnorm_for_dic_normalization) {
          norm_tmp = patchBand_LR.norm();
        }
        patchBand_LR /= norm_tmp;
        patchBand_HR /= norm_tmp;
      }
      dictHR[ig] = SpEOVectorD::Map(patchBand_HR.data(), pszH2);
      dictLR[ig] = SpEOVectorD::Map(patchBand_LR.data(), pszL2);

      check_for_inf_or_nan(my_rank, dictHR[ig], "ig=", ig, ": dictHR[ig]");
      check_for_inf_or_nan(my_rank, dictLR[ig], "ig=", ig, ": dictLR[ig]");
    }
    if (print_for_debugging && (my_rank == 0))
      cout << "[" << my_rank << "] bp-05-1: Init Alpha, DictHR and DictLR"
           << endl;
    SpEOMatrixD *Alpha = new SpEOMatrixD[Ng];
    SpEOMatrixD *Alpha_init = new SpEOMatrixD[Ng];
    SpEOMatrixD *DictHR = new SpEOMatrixD[Ng];
    SpEOMatrixD *DictLR = new SpEOMatrixD[Ng];

    int NDP_min = NDP;
    if (print_for_debugging && (my_rank == 0))
      cout << "[" << my_rank << "] bp-05-2: Fill DictHR and DictLR" << endl;
    for (ig = 0; ig < Ng; ig++) {
      Alpha[ig] = SpEOMatrixD::Zero(NDP_min, Nc_vec[ig]);
      Alpha_init[ig] = SpEOMatrixD::Zero(NDP_min, Nc_vec[ig]);
      DictHR[ig] = SpEOMatrixD::Zero(pszH2, NDP_min);
      DictLR[ig] = SpEOMatrixD::Zero(pszL2, NDP_min);
      SpEOMatrixD patchCompSubsetDouble = SpEOMatrixD::Zero(NDP_min, 3);
      // Generate patchCompSubset matrix of selected patch indices ordered
      // according to respective comparison metric:
      int NDP_tmp = NDP_min;
      dictSelectFunc(&patchCompSubsetDouble, ImX_sim_LR, ImX_sim, ImY, fSet,
                     glPrms, NP, iP, uP, vP, &idxPUH, &idxPVH, &idxPUL, &idxPVL,
                     SRF, &patchCompD, my_rank, ig, NDP_tmp);
      if (NDP_tmp < NDP_min) {
        NDP_min = NDP_tmp;
      }
      // Save dictionary patch coordinates AND sorting metric (column 3)
      if (oSet->saveDicts && iP >= oSet->pFirstDict && iP <= oSet->pLastDict) {
        // write local patch dictionary coordinates to .csv file
        char buf[paths->dir_out.length() + 60];
        sprintf(buf,
                "%s/dictCoords/iP%07d_chX%03d_uP%05d_vP%05d_bundle%03d.csv",
                paths->dir_out.c_str(), iP, ig, uP, vP, ig);
        write_Mat_to_CSV(&patchCompSubsetDouble, buf);
      }
      SpEOMatrixI patchCompSubset = patchCompSubsetDouble.cast<int>();
      int iNP;
      for (iNP = 0; iNP < NDP_min; iNP++) {
        SpEOMatrixD patchHR_tmp =
            (ImX_sim->get_rasterBands()[ig]->bandDataMatD.block(
                 idxPUH.coeff(patchCompSubset(iNP, 0)),
                 idxPVH.coeff(patchCompSubset(iNP, 1)), pszH, pszH))
                .cast<double>();
        SpEOMatrixD patchLR_tmp =
            (ImX_sim_LR->get_rasterBands()[ig]->bandDataMatD.block(
                 idxPUL.coeff(patchCompSubset(iNP, 0)),
                 idxPVL.coeff(patchCompSubset(iNP, 1)), pszL, pszL))
                .cast<double>();
        if (fSet->substrMean) {
          patchLR_tmp.array() -= patchLR_tmp.mean();
          patchHR_tmp.array() -= patchHR_tmp.mean();
        }
        if (fSet->nrmlDicts) {
          double norm_tmp = patchHR_tmp.norm();
          if (fSet->use_LRnorm_for_dic_normalization) {
            norm_tmp = patchLR_tmp.norm();
          }
          if (norm_tmp < 1e-8) {
            cerr << endl
                 << endl
                 << "======================================>" << endl
                 << "=======         WARNING        ========" << endl
                 << "[my_rank=" << my_rank << "] ig= " << ig << ", iNP=" << iNP
                 << ", norm of ImX_sim_LR patch is close to zero!!!: "
                    "patchLR_tmp.norm()="
                 << patchLR_tmp.norm() << endl
                 << "<======================================" << endl
                 << endl;

            patchHR_tmp =
                SpEOMatrixD::Zero(patchHR_tmp.rows(), patchHR_tmp.cols());
            patchHR_tmp.block(0, 0, fDS, fDS) = SpEOMatrixD::Ones(fDS, fDS);
            patchHR_tmp.array() -= patchHR_tmp.mean();
            patchLR_tmp = SpEOVectorD::Zero(patchLR_tmp.rows());
            patchLR_tmp(0) = 1;
            patchLR_tmp.array() -= patchLR_tmp.mean();
            norm_tmp = patchHR_tmp.norm();
            if (fSet->use_LRnorm_for_dic_normalization) {
              norm_tmp = patchLR_tmp.norm();
            }
          }
          patchLR_tmp /= norm_tmp;
          patchHR_tmp /= norm_tmp;
        }
        DictHR[ig].col(iNP) = SpEOVectorD::Map(patchHR_tmp.data(), pszH2);
        DictLR[ig].col(iNP) = SpEOVectorD::Map(patchLR_tmp.data(), pszL2);
      }
      check_for_inf_or_nan(my_rank, DictHR[ig], "ig=", ig, ": DictHR[ig]");
      check_for_inf_or_nan(my_rank, DictLR[ig], "ig=", ig, ": DictLR[ig]");
    }
    if (print_for_debugging && (my_rank == 0))
      cout << "[" << my_rank << "] bp-05-3: check sizes" << endl;
    if (NDP_min < NDP) {
      // ensure that the dictionaries in all groups corresponding to this
      // current patch have the same size (number of atoms)
      NDP = NDP_min;
      for (ig = 0; ig < Ng; ig++) {
        Alpha[ig] = Alpha[ig].block(0, 0, NDP, Nc_vec[ig]);
        Alpha_init[ig] = Alpha_init[ig].block(0, 0, NDP, Nc_vec[ig]);
        DictHR[ig] = DictHR[ig].block(0, 0, pszH2, NDP);
        DictLR[ig] = DictLR[ig].block(0, 0, pszL2, NDP);
      }
    }

    double lambda_X_ABC = fSet->lambdaX_ABC;
    double lambda_Y_ABC = fSet->lambdaY_ABC;
    double lambda_Z_ABC;
    if (iterMain == 0) {
      lambda_Z_ABC = fSet->lambdaZ_ABC_in_1st_iter;
    } else {
      lambda_Z_ABC = fSet->lambdaZ_ABC;
    }
    double lambda_A_ABC = fSet->lambda;
    int optMeanDiffLS_maxiter = 5000;
    double optMeanDiffLS_tol_r = 1e-12;

    // extract patch Z_init from ImZ_init
    if (print_for_debugging && (my_rank == 0))
      cout << "[" << my_rank << "] bp-05-04: extract patch Z_init from ImZ_init"
           << endl;
    SpEOMatrixD Z_init = SpEOMatrixD::Zero(NChY, pszH2);
    for (iChY = 0; iChY < NChY; iChY++) {
      SpEOMatrixD patchBand =
          (ImZ_init->get_rasterBands()[iChY]->bandDataMatD.block(
              idxPUH.coeff(uP), idxPVH.coeff(vP), pszH, pszH));
      Z_init.row(iChY) = SpEOVectorD::Map(patchBand.data(), pszH2);
    }
    check_for_inf_or_nan(my_rank, Z_init, " ", -123, "Z_init");
    SpEOVectorD delta_m_Z_0 = SpEOVectorD::Zero(NChY);
    SpEOVectorD delta_m_Z;
    int optMeanDiffLS_iter;
    double optMeanDiffLS_rel_res;
    // calculate delta_m_Z via least squares
    if (print_for_debugging && (my_rank == 0))
      cout << "[" << my_rank << "] bp-05-05: calculate delta_m_Z" << endl;
    int g, j, k;
    int N_l = patY.cols();
    int N_h = patX.cols();
    int N_Y = SRF->cols();
    int N_X = SRF->rows();
    int N_g = Ng;
    int sum_N_c = 0;
    for (g = 0; g < N_g; g++) {
      sum_N_c += P_lmd_vecs_loc[g].rows();
    }
    double calcOptMeanDiffLS_functional_before =
        (lambda_Z_ABC / NChY) *
            pow((m_Y + delta_m_Z_0 - (1.0 / N_h) * Z_init.rowwise().sum())
                    .norm(),
                2) +
        (lambda_X_ABC / NChX) * pow(((*SRF) * (m_Y + delta_m_Z_0) -
                                     (1.0 / N_h) * patX.rowwise().sum())
                                        .norm(),
                                    2);
    check_for_inf_or_nan(my_rank, calcOptMeanDiffLS_functional_before, " ",
                         -123, "calcOptMeanDiffLS_functional_before");

    calcOptMeanDiffLS(optMeanDiffLS_iter, optMeanDiffLS_rel_res, delta_m_Z,
                      &delta_m_Z_0, SRF, &m_Y, &Z_init, &patX, lambda_X_ABC,
                      lambda_Z_ABC, optMeanDiffLS_maxiter, optMeanDiffLS_tol_r);
    check_for_inf_or_nan(my_rank, optMeanDiffLS_rel_res, " ", -123,
                         "optMeanDiffLS_rel_res");
    check_for_inf_or_nan(my_rank, delta_m_Z, " ", -123, "delta_m_Z");
    double calcOptMeanDiffLS_functional_after =
        (lambda_Z_ABC / NChY) *
            pow((m_Y + delta_m_Z - (1.0 / N_h) * Z_init.rowwise().sum()).norm(),
                2) +
        (lambda_X_ABC / NChX) * pow(((*SRF) * (m_Y + delta_m_Z) -
                                     (1.0 / N_h) * patX.rowwise().sum())
                                        .norm(),
                                    2);
    check_for_inf_or_nan(my_rank, calcOptMeanDiffLS_functional_after, " ", -123,
                         "calcOptMeanDiffLS_functional_after");

    if (print_optimization_input_output && my_rank == 0) {
      cout << "         "
              "----------------------------------------------------------------"
              "--------------------------"
           << endl;
      cout << "         [" << my_rank
           << "] calcOptMeanDiffLS (Eq.6):  iter=" << optMeanDiffLS_iter
           << ", rel_res=" << optMeanDiffLS_rel_res
           << ", functional (before , after)=( "
           << calcOptMeanDiffLS_functional_before << " , "
           << calcOptMeanDiffLS_functional_after << " )" << endl;
      cout << "         "
              "----------------------------------------------------------------"
              "--------------------------"
           << endl;
      cout
          << "                   iChY  | m_Y(iChY) | m_Z(iChY) | delta_m_Z"
          << endl
          << "                   ----------------------------------------------"
          << endl;
      for (iChY = 0; iChY < NChY; iChY++) {
        cout << fixed << setprecision(0) << "                   ";
        cout.fill(' ');
        cout.width(4);
        cout << iChY;
        cout << "  |   ";
        cout.fill(' ');
        cout.width(6);
        cout << m_Y(iChY) << "  |  ";
        cout.fill(' ');
        cout.width(6);
        cout << m_Y(iChY) + delta_m_Z(iChY) << "   |";
        cout.fill(' ');
        cout.width(4);
        cout << delta_m_Z(iChY) << endl;
      }
    }

    // ******************************************************************************
    // calcOptCoeffCurrPatchLS:
    // calculate magnitudes of coefficients corresponding to the current patch
    // under reconstruction in each channel group separately
    // (i.e. alpha which is the first row in Alpha)
    // ******************************************************************************
    if (print_for_debugging && (my_rank == 0))
      cout << "[" << my_rank << "] bp-06: calcOptCoeffCurrPatchLS" << endl;
    int optCoeffCurrPatchLS_maxiter = 5000;
    double optCoeffCurrPatchLS_tol_r = 1e-12;
    int optCoeffCurrPatchLS_iter;
    double optCoeffCurrPatchLS_rel_res;
    // calculate alpha via least squares
    // calculate Functional Value before solver call
    if (print_for_debugging && (my_rank == 0))
      cout << "[" << my_rank
           << "] bp-06-01: calculate Functional Value before solver call"
           << endl;
    SpEOVectorD onesNh = SpEOVectorD::Ones(N_h);
    SpEOVectorD onesNl = SpEOVectorD::Ones(N_l);
    double fct_value = 0.0;
    SpEOMatrixD tmp_upper =
        (Z_init - (m_Y + delta_m_Z) * onesNh.transpose()).transpose();
    SpEOMatrixD tmp_center = (patY - (m_Y * onesNl.transpose())).transpose();
    for (g = 0; g < N_g; g++) {
      fct_value += lambda_Z_ABC / (2 * N_h * N_g * alpha_init[g].size()) *
                   pow((dictHR[g] * alpha_init[g].transpose() -
                        tmp_upper.middleCols(P_lmd_idx_bl_loc[g][0],
                                             alpha_init[g].size()))
                           .norm(),
                       2);
      fct_value += lambda_Y_ABC / (2 * N_l * N_g * alpha_init[g].size()) *
                   pow((dictLR[g] * alpha_init[g].transpose() -
                        tmp_center.middleCols(P_lmd_idx_bl_loc[g][0],
                                              alpha_init[g].size()))
                           .norm(),
                       2);
    }
    SpEOMatrixD tmpDA = SpEOMatrixD::Zero(sum_N_c, N_h);
    int running_idx = 0;
    for (g = 0; g < N_g; g++) {
      tmpDA.block(running_idx, 0, P_lmd_vecs_loc[g].size(), N_h) =
          (dictHR[g] * alpha_init[g].transpose()).transpose();
      running_idx += P_lmd_vecs_loc[g].size();
    }
    SpEOMatrixD tmpPDA = SpEOMatrixD::Zero(N_Y, N_h);
    fct_value +=
        0.5 * lambda_X_ABC / (N_h * N_X) *
        pow((patX - (*SRF) * (tmpPDA + (m_Y + delta_m_Z) * onesNh.transpose()))
                .norm(),
            2);
    double calcOptCoeffCurrPatchLS_functional_before = fct_value;
    check_for_inf_or_nan(my_rank, calcOptCoeffCurrPatchLS_functional_before,
                         " ", -123,
                         "calcOptCoeffCurrPatchLS_functional_before");
    // call solver
    if (print_for_debugging && (my_rank == 0))
      cout << "[" << my_rank << "] bp-06-02: call solver" << endl;
    calcOptCoeffCurrPatchLS(
        optCoeffCurrPatchLS_iter, optCoeffCurrPatchLS_rel_res, alpha,
        alpha_init, dictHR, dictLR, &m_Y, &delta_m_Z, &Z_init, &patY, &patX,
        SRF, P_lmd_vecs_loc, P_lmd_idx_row_loc, P_lmd_idx_bl_loc, lambda_X_ABC,
        lambda_Y_ABC, lambda_Z_ABC, Ng, optCoeffCurrPatchLS_maxiter,
        optCoeffCurrPatchLS_tol_r);

    check_for_inf_or_nan(my_rank, optCoeffCurrPatchLS_rel_res, " ", -123,
                         "optCoeffCurrPatchLS_rel_res");
    for (int ig = 0; ig < N_g; ig++) {
      check_for_inf_or_nan(my_rank, alpha[ig], "ig=", ig, ": alpha[ig]");
    }
    // calculate Functional Value after solver call
    if (print_for_debugging && (my_rank == 0))
      cout << "[" << my_rank
           << "] bp-06-03: calculate Functional Value after solver call"
           << endl;
    fct_value = 0.0;
    tmp_upper = (Z_init - (m_Y + delta_m_Z) * onesNh.transpose()).transpose();
    tmp_center = (patY - (m_Y * onesNl.transpose())).transpose();
    for (g = 0; g < N_g; g++) {
      fct_value += lambda_Z_ABC / (2 * N_h * N_g * alpha_init[g].size()) *
                   pow((dictHR[g] * alpha[g].transpose() -
                        tmp_upper.middleCols(P_lmd_idx_bl_loc[g][0],
                                             alpha_init[g].size()))
                           .norm(),
                       2);
      fct_value += lambda_Y_ABC / (2 * N_l * N_g * alpha_init[g].size()) *
                   pow((dictLR[g] * alpha[g].transpose() -
                        tmp_center.middleCols(P_lmd_idx_bl_loc[g][0],
                                              alpha_init[g].size()))
                           .norm(),
                       2);
    }
    tmpDA = SpEOMatrixD::Zero(sum_N_c, N_h);
    running_idx = 0;
    for (g = 0; g < N_g; g++) {
      tmpDA.block(running_idx, 0, P_lmd_vecs_loc[g].size(), N_h) =
          (dictHR[g] * alpha[g].transpose()).transpose();
      running_idx += P_lmd_vecs_loc[g].size();
    }
    tmpPDA = SpEOMatrixD::Zero(N_Y, N_h);
    for (j = 0; j < N_Y; j++) {
      for (k = 0; k < P_lmd_idx_row_loc[j].cols(); k++) {
        int block = P_lmd_idx_row_loc[j].coeff(0, k);
        int relidx = P_lmd_idx_row_loc[j].coeff(1, k);
        tmpPDA.row(j) += P_lmd_vecs_loc[block](relidx) *
                         tmpDA.row(P_lmd_idx_bl_loc[block][1] + relidx);
      }
    }
    fct_value +=
        0.5 * lambda_X_ABC / (N_h * N_X) *
        pow((patX - (*SRF) * (tmpPDA + (m_Y + delta_m_Z) * onesNh.transpose()))
                .norm(),
            2);
    // output
    if (print_for_debugging && (my_rank == 0))
      cout << "[" << my_rank << "] bp-06-04: print output" << endl;
    double calcOptCoeffCurrPatchLS_functional_after = fct_value;
    check_for_inf_or_nan(my_rank, calcOptCoeffCurrPatchLS_functional_after, " ",
                         -123, "calcOptCoeffCurrPatchLS_functional_after");
    if (print_optimization_input_output && my_rank == 0) {
      cout << "         "
              "----------------------------------------------------------------"
              "--------------------------"
           << endl;
      cout << "         [" << my_rank
           << "] optCoeffCurrPatchLS (Eq.7):  iter=" << optCoeffCurrPatchLS_iter
           << ",   rel_res=" << optCoeffCurrPatchLS_rel_res
           << ",   functional (before , after) = ( "
           << calcOptCoeffCurrPatchLS_functional_before << " , "
           << calcOptCoeffCurrPatchLS_functional_after << " )" << endl;
      cout << "         "
              "----------------------------------------------------------------"
              "--------------------------"
           << endl;

      for (int ig = 0; ig < Ng; ig++) {
        cout << "           alpha[ig=" << ig << "] = " << alpha[ig].transpose()
             << endl;
      }
    }

    // copy the vector alpha to the first row of Alpha and Alpha_init
    if (print_for_debugging && (my_rank == 0))
      cout << "[" << my_rank
           << "] bp-06-05: copy the vector alpha to the first row of Alpha and "
              "Alpha_init"
           << endl;
    for (int ig = 0; ig < Ng; ig++) {
      Alpha[ig].row(0) = SpEOVectorD::Map(alpha[ig].data(), Alpha[ig].cols());
      Alpha_init[ig].row(0) = Alpha[ig].row(0);
    }

    // ******************************************************************************
    // calcOptCoeffResPFISTA:
    // for each spectral group:
    // do a preliminary joint sparse reconstruction of the residuals to find
    // the column-sparse coefficient matrices Alpha_res
    // ******************************************************************************

    if (print_for_debugging && (my_rank == 0))
      cout << "[" << my_rank << "] bp-07: calcOptCoeffResPFISTA" << endl;
    FBSsolverOptions opts = FBSsolverOptions(
        1,      // decomposition parameter
        1,      // inner iterations (1 for standard FISTA)
        10000,  // max. outer iterations (former optCoeffResPFISTA_maxiter)
        1e-3,   // respective tolerance (dep. on stopping crit.) (former
                // optCoeffResPFISTA_tol_r)
        FIRST_ORDER_CONDITIONS,  // stopping criterion
        FISTA,                   // prediction step rule
        DECPAR,                  // backtracking strategy NO, INC, or DEC
        INIT, false, true, 1e2);
    if (print_for_debugging && (my_rank == 0))
      cout << "[" << my_rank << "] bp-07-01: after fnct call" << endl;
    int spectral_normalizer = fSet->matrixNorm;
    int *optCoeffResPFISTA_iter = new int[Ng];
    double *optCoeffResPFISTA_rel_res = new double[Ng];
    fct_value = 0.0;
    SpEOMatrixD tmp_upper1 =
        (Z_init - (m_Y + delta_m_Z) * onesNh.transpose()).transpose();
    SpEOMatrixD tmp_center1 = (patY - (m_Y)*onesNl.transpose()).transpose();
    check_for_inf_or_nan(my_rank, tmp_upper1, " ", -123, "tmp_upper1");
    check_for_inf_or_nan(my_rank, tmp_center1, " ", -123, "tmp_center1");
    for (int ig = 0; ig < N_g; ig++) {
      fct_value +=
          lambda_A_ABC / sqrt(Alpha_init[ig].cols()) * l1l2norm(Alpha_init[ig]);
      check_for_inf_or_nan(
          my_rank, fct_value, "ig=", ig,
          ": fct_value1 (for calcOptCoeffResPFISTA_functional_before)");
      fct_value += lambda_Z_ABC / (2 * N_h * Alpha_init[ig].cols()) *
                   pow((DictHR[ig] * Alpha_init[ig] -
                        tmp_upper1.middleCols(P_lmd_idx_bl_loc[ig][0],
                                              Alpha_init[ig].cols()))
                           .norm(),
                       2);
      check_for_inf_or_nan(
          my_rank, fct_value, "ig=", ig,
          ": fct_value2 (for calcOptCoeffResPFISTA_functional_before)");
      fct_value += lambda_Y_ABC / (2 * N_l * Alpha_init[ig].cols()) *
                   pow((DictLR[ig] * Alpha_init[ig] -
                        tmp_center1.middleCols(P_lmd_idx_bl_loc[ig][0],
                                               Alpha_init[ig].cols()))
                           .norm(),
                       2);
      check_for_inf_or_nan(
          my_rank, fct_value, "ig=", ig,
          ": fct_value3 (for calcOptCoeffResPFISTA_functional_before)");
    }
    double calcOptCoeffResPFISTA_functional_before = fct_value;
    check_for_inf_or_nan(my_rank, calcOptCoeffResPFISTA_functional_before, " ",
                         -123, "calcOptCoeffResPFISTA_functional_before");
    if (print_for_debugging && (my_rank == 0))
      cout << "[" << my_rank << "] bp-07-02: calling calcOptCoeffResPFISTA"
           << endl;
    calcOptCoeffResPFISTA(optCoeffResPFISTA_iter, optCoeffResPFISTA_rel_res,
                          Alpha, Alpha_init, DictHR, DictLR, &m_Y, &delta_m_Z,
                          &Z_init, &patY, P_lmd_idx_bl_loc, lambda_A_ABC,
                          lambda_Y_ABC, lambda_Z_ABC, Ng, opts,
                          spectral_normalizer, write_testset, testnr);
    fct_value = 0.0;
    int iter_mean = 0;
    double rel_res_mean = 0.0;
    for (g = 0; g < N_g; g++) {
      fct_value += lambda_A_ABC / sqrt(Alpha[g].cols()) * l1l2norm(Alpha[g]);
      fct_value +=
          lambda_Z_ABC / (2 * N_h * Alpha[g].cols()) *
          pow((DictHR[g] * Alpha[g] -
               tmp_upper1.middleCols(P_lmd_idx_bl_loc[g][0], Alpha[g].cols()))
                  .norm(),
              2);
      fct_value +=
          lambda_Y_ABC / (2 * N_l * Alpha[g].cols()) *
          pow((DictLR[g] * Alpha[g] -
               tmp_center1.middleCols(P_lmd_idx_bl_loc[g][0], Alpha[g].cols()))
                  .norm(),
              2);
      iter_mean += optCoeffResPFISTA_iter[g] / N_g;
      rel_res_mean += optCoeffResPFISTA_rel_res[g] / N_g;
    }
    double calcOptCoeffResPFISTA_functional_after = fct_value;
    if (print_optimization_input_output && my_rank == 0) {
      cout << "         "
              "----------------------------------------------------------------"
              "--------------------------"
           << endl;
      cout << "         [" << my_rank
           << "] calcOptCoeffResPFISTA (Eq.8):   iter(mean)=" << iter_mean
           << ",   rel_res(mean)=" << rel_res_mean
           << ",   functional sum (before , after) = ( "
           << calcOptCoeffResPFISTA_functional_before << " , "
           << calcOptCoeffResPFISTA_functional_after << " )" << endl;
      cout << "                                                            "
              "iter = ";
      for (g = 0; g < N_g; g++) {
        cout << optCoeffResPFISTA_iter[g] << " ";
      }
      cout << endl;
      cout << "                                                            "
              "rel_res = ";
      for (g = 0; g < N_g; g++) {
        cout << optCoeffResPFISTA_rel_res[g] << " ";
      }
      cout << endl;
      cout << "         "
              "----------------------------------------------------------------"
              "--------------------------"
           << endl;
    }
    for (int ig = 0; ig < N_g; ig++) {
      check_for_inf_or_nan(my_rank, optCoeffResPFISTA_rel_res[ig], "ig=", ig,
                           ": optCoeffResPFISTA_rel_res[ig]");
      check_for_inf_or_nan(my_rank, Alpha[ig], "ig=", ig, ": Alpha[ig]");
    }
    delete[] optCoeffResPFISTA_iter;
    delete[] optCoeffResPFISTA_rel_res;

    // ******************************************************************************
    // for each spectral group:
    // identify Omega = support(Alpha), i.e. find the
    // non-trivial rows of the jointly sparse coefficients in Alpha
    // ******************************************************************************
    if (print_for_debugging && (my_rank == 0))
      cout
          << "[" << my_rank
          << "] bp-08: for each spectral group: identify Omega = support(Alpha)"
          << endl;
    double tol_support = 1e-13;
    SpEOVectorI *support = new SpEOVectorI[Ng];
    for (int ig = 0; ig < Ng; ig++) {
      int support_tmp[NDP_min];
      int support_num = 1;
      // Assure that the first atom will be contained in the support
      support_tmp[0] = 0;
      // Calculate the remaining support from the remaining atoms 2,...,NDP_min
      for (int idp = 1; idp < NDP_min; idp++) {
        if (Alpha[ig].row(idp).norm() > tol_support) {
          support_tmp[support_num] = idp;
          support_num++;
        }
      }
      support[ig] = SpEOVectorI::Zero(support_num);
      for (int iSup = 0; iSup < support_num; iSup++) {
        support[ig](iSup) = support_tmp[iSup];
      }
      check_for_inf_or_nan(my_rank, support[ig], "ig=", ig, ": support[ig]");
    }
    if (print_optimization_input_output && my_rank == 0) {
      for (int ig = 0; ig < Ng; ig++) {
        cout << "           [" << my_rank << "] support[ig=" << ig
             << "] =  " << support[ig].transpose() << endl;
      }
    }

    // ******************************************************************************
    // for each spectral group:
    // reduce Alpha, Alpha_init, DictHR and DictLR to
    // their entries on the support of Alpha
    // ******************************************************************************/
    if (print_for_debugging && (my_rank == 0))
      cout << "[" << my_rank
           << "] bp-09: for each spectral group: reduce Alpha, Alpha_init, "
              "DictHR and DictLR"
           << endl;
    SpEOMatrixD *Alpha_red = new SpEOMatrixD[Ng];
    SpEOMatrixD *Alpha_red_init = new SpEOMatrixD[Ng];
    SpEOMatrixD *DictHR_red = new SpEOMatrixD[Ng];
    SpEOMatrixD *DictLR_red = new SpEOMatrixD[Ng];
    int supp_numel_sum = 0;
    for (int ig = 0; ig < Ng; ig++) {
      int supp_numel = support[ig].rows();
      supp_numel_sum += supp_numel;
      Alpha_red[ig] = SpEOMatrixD::Zero(supp_numel, Nc_vec[ig]);
      Alpha_red_init[ig] = SpEOMatrixD::Zero(supp_numel, Nc_vec[ig]);
      DictHR_red[ig] = SpEOMatrixD::Zero(pszH2, supp_numel);
      DictLR_red[ig] = SpEOMatrixD::Zero(pszL2, supp_numel);
      for (int iSup = 0; iSup < supp_numel; iSup++) {
        Alpha_red[ig].row(iSup) = Alpha[ig].row(support[ig].coeff(iSup));
        Alpha_red_init[ig].row(iSup) = Alpha[ig].row(support[ig].coeff(iSup));
        DictHR_red[ig].col(iSup) = DictHR[ig].col(support[ig].coeff(iSup));
        DictLR_red[ig].col(iSup) = DictLR[ig].col(support[ig].coeff(iSup));
      }
      check_for_inf_or_nan(my_rank, Alpha_red[ig], "ig=", ig,
                           ": Alpha_red[ig]");
      check_for_inf_or_nan(my_rank, Alpha_red_init[ig], "ig=", ig,
                           ": Alpha_red_init[ig]");
      check_for_inf_or_nan(my_rank, DictHR_red[ig], "ig=", ig,
                           ": DictHR_red[ig]");
      check_for_inf_or_nan(my_rank, DictLR_red[ig], "ig=", ig,
                           ": DictLR_red[ig]");
    }

    // ******************************************************************************
    // calcOptCoeffResLS:
    // Calculate the Alpha_res (rows 2 to end of Alpha) by re-estimating the
    // magnitudes of Alpha_res_init on the support Omega via least squares
    // ******************************************************************************
    if (print_for_debugging && (my_rank == 0))
      cout << "[" << my_rank << "] bp-10: calcOptCoeffResLS" << endl;
    int optCoeffResLS_maxiter = 5000;
    double optCoeffResLS_tol_r = 1e-12;
    int optCoeffResLS_iter;
    double optCoeffResLS_rel_res;
    // calculate Functional Value before solver call
    fct_value = 0.0;
    tmp_upper = (Z_init - (m_Y + delta_m_Z) * onesNh.transpose()).transpose();
    tmp_center = (patY - (m_Y * onesNl.transpose())).transpose();
    for (g = 0; g < N_g; g++) {
      fct_value += lambda_Z_ABC / (2 * N_h * N_g * Alpha_red_init[g].cols()) *
                   pow((DictHR_red[g] * Alpha_red_init[g] -
                        tmp_upper.middleCols(P_lmd_idx_bl_loc[g][0],
                                             Alpha_red_init[g].cols()))
                           .norm(),
                       2);
      fct_value += lambda_Y_ABC / (2 * N_l * N_g * Alpha_red_init[g].cols()) *
                   pow((DictLR_red[g] * Alpha_red_init[g] -
                        tmp_center.middleCols(P_lmd_idx_bl_loc[g][0],
                                              Alpha_red_init[g].cols()))
                           .norm(),
                       2);
    }
    tmpDA = SpEOMatrixD::Zero(sum_N_c, N_h);
    running_idx = 0;
    for (g = 0; g < N_g; g++) {
      tmpDA.block(running_idx, 0, P_lmd_vecs_loc[g].size(), N_h) =
          (DictHR_red[g] * Alpha_red_init[g]).transpose();
      running_idx += P_lmd_vecs_loc[g].size();
    }
    tmpPDA = SpEOMatrixD::Zero(N_Y, N_h);
    for (j = 0; j < N_Y; j++) {
      for (k = 0; k < P_lmd_idx_row_loc[j].cols(); k++) {
        int block = P_lmd_idx_row_loc[j].coeff(0, k);
        int relidx = P_lmd_idx_row_loc[j].coeff(1, k);
        tmpPDA.row(j) += P_lmd_vecs_loc[block](relidx) *
                         tmpDA.row(P_lmd_idx_bl_loc[block][1] + relidx);
      }
    }
    fct_value +=
        0.5 * lambda_X_ABC / (N_h * N_X) *
        pow((patX - (*SRF) * (tmpPDA + (m_Y + delta_m_Z) * onesNh.transpose()))
                .norm(),
            2);
    double calcOptCoeffResLS_functional_before = fct_value;
    check_for_inf_or_nan(my_rank, calcOptCoeffResLS_functional_before, " ",
                         -123, "calcOptCoeffResLS_functional_before");
    check_for_inf_or_nan(my_rank, tmpPDA, " ", -123, "tmpPDA");
    // call solver
    if (print_for_debugging && (my_rank == 0))
      cout << "[" << my_rank << "] bp-10-01: call solver" << endl;
    if (supp_numel_sum >
        N_g) {  // otherwise it's a trivial case that every dictionary (the
                // dictionary in any group) contains only the current patch ->
                // this step can be skipped
      calcOptCoeffResLS(optCoeffResLS_iter, optCoeffResLS_rel_res, Alpha_red,
                        Alpha_red_init, DictHR_red, DictLR_red, &m_Y,
                        &delta_m_Z, &Z_init, &patY, &patX, SRF, P_lmd_vecs_loc,
                        P_lmd_idx_row_loc, P_lmd_idx_bl_loc, lambda_X_ABC,
                        lambda_Y_ABC, lambda_Z_ABC, Ng, optCoeffResLS_maxiter,
                        optCoeffResLS_tol_r);
    } else {
      cout << "         [" << my_rank << "] WWWW Warning:  (iP,uP,vP)=(" << iP
           << "," << uP << "," << vP
           << "): After FISTA, Alpha[iG] contains only current patch for all "
              "iG=0,...,Ng-1 => Eq. 10 is skipped!"
           << endl;
      optCoeffResLS_iter = 0;
      optCoeffResLS_rel_res = -1;
    }
    check_for_inf_or_nan(my_rank, optCoeffResLS_rel_res, " ", -123,
                         "optCoeffResLS_rel_res");
    // calculate Functional Value after solver call
    if (print_for_debugging && (my_rank == 0))
      cout << "[" << my_rank << "] bp-10-02: calc functional value" << endl;
    fct_value = 0.0;
    tmp_upper = (Z_init - (m_Y + delta_m_Z) * onesNh.transpose()).transpose();
    tmp_center = (patY - (m_Y * onesNl.transpose())).transpose();
    for (g = 0; g < N_g; g++) {
      fct_value += lambda_Z_ABC / (2 * N_h * N_g * Alpha_red_init[g].cols()) *
                   pow((DictHR_red[g] * Alpha_red[g] -
                        tmp_upper.middleCols(P_lmd_idx_bl_loc[g][0],
                                             Alpha_red_init[g].cols()))
                           .norm(),
                       2);
      fct_value += lambda_Y_ABC / (2 * N_l * N_g * Alpha_red_init[g].cols()) *
                   pow((DictLR_red[g] * Alpha_red[g] -
                        tmp_center.middleCols(P_lmd_idx_bl_loc[g][0],
                                              Alpha_red_init[g].cols()))
                           .norm(),
                       2);
    }
    tmpDA = SpEOMatrixD::Zero(sum_N_c, N_h);
    running_idx = 0;
    for (g = 0; g < N_g; g++) {
      tmpDA.block(running_idx, 0, P_lmd_vecs_loc[g].size(), N_h) =
          (DictHR_red[g] * Alpha_red[g]).transpose();
      running_idx += P_lmd_vecs_loc[g].size();
    }
    tmpPDA = SpEOMatrixD::Zero(N_Y, N_h);
    for (j = 0; j < N_Y; j++) {
      for (k = 0; k < P_lmd_idx_row_loc[j].cols(); k++) {
        int block = P_lmd_idx_row_loc[j].coeff(0, k);
        int relidx = P_lmd_idx_row_loc[j].coeff(1, k);
        tmpPDA.row(j) += P_lmd_vecs_loc[block](relidx) *
                         tmpDA.row(P_lmd_idx_bl_loc[block][1] + relidx);
      }
    }
    fct_value +=
        0.5 * lambda_X_ABC / (N_h * N_X) *
        pow((patX - (*SRF) * (tmpPDA + (m_Y + delta_m_Z) * onesNh.transpose()))
                .norm(),
            2);
    // output
    if (print_for_debugging && (my_rank == 0))
      cout << "[" << my_rank << "] bp-10-03: output" << endl;
    double calcOptCoeffResLS_functional_after = fct_value;
    if (print_optimization_input_output && my_rank == 0) {
      cout << "         "
              "----------------------------------------------------------------"
              "--------------------------"
           << endl;
      cout << "         [" << my_rank
           << "] calcOptCoeffResLS (Eq.10):  iter=" << optCoeffResLS_iter
           << ",   rel_res=" << optCoeffResLS_rel_res
           << ",   functional (before , after) = ( "
           << calcOptCoeffResLS_functional_before << " , "
           << calcOptCoeffResLS_functional_after << " )" << endl;
      cout << "         "
              "----------------------------------------------------------------"
              "--------------------------"
           << endl;

      for (int ig = 0; ig < Ng; ig += max(1, (Ng / 3) - 2)) {
        if (Nc_vec[ig] < 13) {
          cout << "         Alpha_red_init[ig=" << ig << "]:" << endl;
          cout << Alpha_red_init[ig] << endl;
          cout << "         Alpha_red     [ig=" << ig << "]:" << endl
               << Alpha_red[ig] << endl;
        }
      }
    }
    check_for_inf_or_nan(my_rank, calcOptCoeffResLS_functional_after, " ", -123,
                         "calcOptCoeffResLS_functional_after");
    // ******************************************************************************
    // Calculate the new patch patZ
    // ******************************************************************************
    if (print_for_debugging && (my_rank == 0))
      cout << "[" << my_rank << "] bp-11: calculate the new patch patZ" << endl;
    patZ = tmpPDA + (m_Y + delta_m_Z) * onesNh.transpose();
    check_for_inf_or_nan(my_rank, patZ, " ", -123, "patZ");

    zHR = patZ.transpose();
    if (fSet->set_neg_to_0 == 1 || fSet->set_neg_to_0 == 3) {
      zHR = (zHR.array() > 0).select(zHR, 0);
    }

    check_for_inf_or_nan(my_rank, zHR, " ", -123, "zHR");

    delete[] support;
    delete[] alpha;
    delete[] alpha_init;
    delete[] Alpha;
    delete[] Alpha_init;
    delete[] Alpha_red;
    delete[] Alpha_red_init;
    delete[] dictHR;
    delete[] dictLR;
    delete[] DictHR;
    delete[] DictLR;
    delete[] DictHR_red;
    delete[] DictLR_red;
    for (ig = 0; ig < Ng; ig++) {
      delete[] P_lmd_idx_bl_loc[ig];
    }
    delete[] P_lmd_idx_bl_loc;
    delete[] P_lmd_idx_row_loc;
    delete[] P_lmd_vecs_loc;
    delete ImX_sim;
    delete ImX_sim_LR;

    //###########################################################################
    //
    //        ... new model for calculating the patch patZ ends here
    //
    //###########################################################################

    if (pSet->store_patches_tmp_on_drive) {
      char buf0[paths->dir_tmp_patches.length() + 7];
      sprintf(buf0, "%s/%06d", paths->dir_tmp_patches.c_str(), iP);
      mkdir(buf0, 0777);

      char buf[paths->dir_tmp_patches.length() + 42];
      sprintf(buf, "%s/%06d/patch_u%04d_v%04d_iP%06d.csv",
              paths->dir_tmp_patches.c_str(), iP, uP, vP, iP);
      SpEOMatrixF zHRfloat = zHR.cast<float>();
      write_Mat_to_CSV(&zHRfloat, buf);
    } else {
      bool boundary_patch = false;
      int pszHU_tmp = pszH;
      // special treatment of patches that exceed the bottom image boundary
      // (only relevant of the reconstructed image is a spatial subset of the
      // input data)
      if (ImZ->get_sizeU() + dSet->uLFirst * fDS < glPrms->sizeUH &&
          idxPUH.coeff(uP) + pszH >
              ImZ->rasterBands[0]->bandDataMatD.rows() + dSet->uLFirst * fDS) {
        boundary_patch = true;
        pszHU_tmp = ImZ->rasterBands[0]->bandDataMatD.rows() -
                    (idxPUH.coeff(uP) - dSet->uLFirst * fDS);
      }
      int pszHV_tmp = pszH;
      // special treatment of patches that exceed the bottom image boundary
      // (only relevant of the reconstructed image is a spatial subset of the
      // input data)
      if (ImZ->get_sizeV() + dSet->vLFirst * fDS < glPrms->sizeVH &&
          idxPVH.coeff(vP) + pszH >
              ImZ->rasterBands[0]->bandDataMatD.cols() + dSet->vLFirst * fDS) {
        boundary_patch = true;
        pszHV_tmp = ImZ->rasterBands[0]->bandDataMatD.cols() -
                    (idxPVH.coeff(vP) - dSet->vLFirst * fDS);
      }
      SpEOMatrixD patchHR, patchHR_red;
      SpEOVectorD patchVec;

      for (int iChY = 0; iChY < NChY; iChY++) {
        patchVec = zHR.col(iChY);
        patchHR = SpEOMatrixD::Map(patchVec.data(), pszH, pszH);

        int start_uH_local_imZ = idxPUH.coeff(uP) - dSet->uLFirst * fDS;
        int start_vH_local_imZ = idxPVH.coeff(vP) - dSet->vLFirst * fDS;

        if (boundary_patch) {
          patchHR_red = SpEOMatrixD::Zero(pszHU_tmp, pszHV_tmp);
          patchHR_red = patchHR.block(0, 0, pszHU_tmp, pszHV_tmp);
          ImZ->rasterBands[iChY]->bandDataMatD.block(
              start_uH_local_imZ, start_vH_local_imZ, pszHU_tmp, pszHV_tmp) +=
              patchHR_red;
        } else {
          ImZ->rasterBands[iChY]->bandDataMatD.block(
              start_uH_local_imZ, start_vH_local_imZ, pszH, pszH) += patchHR;
        }
      }
    }
    zHR = SpEOMatrixD::Zero(pszH2, NChY);

    if (jP <= lastFixedPatch - glPrms->numPatchGroups ||
        pSet->workStealingTurns < 0) {
      jP += glPrms->numPatchGroups;
    } else {
#ifndef _OPENMP
      if (my_rank % pSet->numProcPerPatch == 0) {
        jP = increment_counter(c, 1) -
             1;  // the root process of each group gets a patch number to work
                 // on, by the current counter
      }
      int jP_tmp[1];
      if (my_rank % pSet->numProcPerPatch == 0) {
        jP_tmp[0] = jP;
      }
      MPI_Bcast(
          &jP_tmp, 1, MPI_INT, 0,
          comm_patch);  // broadcast the patch to work on to the group members
      jP = jP_tmp[0];
#else
      jP = increment_counter(c, 1) - 1;
#endif
    }
    if (print_for_debugging && (my_rank == 0))
      cout << "[" << my_rank << "] bp_3" << endl;
  }
  if (print_for_debugging && (my_rank == 0))
    cout << "[" << my_rank << "] bp_4" << endl;
  if (print_for_debugging && (my_rank == 0))
    cout << "[" << my_rank << "] bp_5" << endl;
  MPI_Barrier(comm_busy);
  if (print_for_debugging && (my_rank == 0))
    cout << "[" << my_rank << "] bp_6" << endl;
  delete_counter(&c);  // the barrier is necessary here (otherwise the counter
                       // is deleted, while others might access it)
  if (print_for_debugging && (my_rank == 0))
    cout << "[" << my_rank << "] bp_7" << endl;

    // TIME STATISTICS
    // obtain additional idling time for proc != 0
#ifndef _OPENMP
  MPI_Barrier(comm_busy);
#endif

  if (my_rank == 0) {
    cout << endl
         << "########################################################### "
         << endl
         << "#        algorithm (patch-wise fusion) completed          # "
         << endl
         << "########################################################### "
         << endl
         << endl;
  }

  if (print_for_debugging && (my_rank == 0))
    cout << "[" << my_rank << "] bp_9" << endl;
  // communicate all patches and take the average on the overlapping regions
  if (!pSet->store_patches_tmp_on_drive) {
    // take the average on the overlapping regions
    int pszHU_tmp, pszHV_tmp;
    SpEOMatrixD avg_mat =
        SpEOMatrixD::Zero(ImZ->rasterBands[0]->bandDataMatD.rows(),
                          ImZ->rasterBands[0]->bandDataMatD.cols());
    for (jP = 0; jP <= pLast_sub; jP++) {
      uP = uPFirst + jP / NPV_sub;
      vP = vPFirst + jP % NPV_sub;
      pszHU_tmp = pszH;
      // special treatment of patches that exceed the bottom image boundary
      // (only relevant of the reconstructed image is a spatial subset of the
      // input data)
      if (dSet->uLFirst * fDS + ImZ->get_sizeU() < glPrms->sizeUH &&
          idxPUH.coeff(uP) + pszH >
              dSet->uLFirst * fDS + ImZ->rasterBands[0]->bandDataMatD.rows()) {
        pszHU_tmp = ImZ->rasterBands[0]->bandDataMatD.rows() -
                    idxPUH.coeff(uP) + dSet->uLFirst * fDS;
      }
      pszHV_tmp = pszH;
      // special treatment of patches that exceed the right image boundary (only
      // relevant of the reconstructed image is a spatial subset of the input
      // data)
      if (dSet->vLFirst * fDS + ImZ->get_sizeV() < glPrms->sizeVH &&
          idxPVH.coeff(vP) + pszH >
              dSet->vLFirst * fDS + ImZ->rasterBands[0]->bandDataMatD.cols()) {
        pszHV_tmp = ImZ->rasterBands[0]->bandDataMatD.cols() -
                    idxPVH.coeff(vP) + dSet->vLFirst * fDS;
      }
      avg_mat.block(idxPUH.coeff(uP) - dSet->uLFirst * fDS,
                    idxPVH.coeff(vP) - dSet->vLFirst * fDS, pszHU_tmp,
                    pszHV_tmp) += SpEOMatrixD::Ones(pszHU_tmp, pszHV_tmp);
    }
    // communicate all results band by band
    for (int iChY = 0; iChY < NChY; iChY++) {
      double *buffer = new double[ImZ->get_sizeU() * ImZ->get_sizeV()];
      MPI_Reduce(ImZ->rasterBands[iChY]->bandDataMatD.data(), buffer,
                 ImZ->get_sizeU() * ImZ->get_sizeV(), MPI_DOUBLE, MPI_SUM, 0,
                 comm_busy);
      ImZ->rasterBands[iChY]->bandDataMatD =
          SpEOMatrixD::Map(buffer, ImZ->get_sizeU(), ImZ->get_sizeV());
      delete[] buffer;
      int datasize = ImZ->rasterBands[iChY]->bandDataMatD.rows() *
                     ImZ->rasterBands[iChY]->bandDataMatD.cols();
      MPI_Bcast((ImZ->rasterBands[iChY]->bandDataMatD.data()), datasize,
                MPI_DOUBLE, 0, comm_busy);
    }
    MPI_Barrier(comm_busy);
    for (int iChY = 0; iChY < NChY; iChY++) {
      ImZ->rasterBands[iChY]->bandDataMatD =
          (ImZ->rasterBands[iChY]->bandDataMatD).cwiseQuotient(avg_mat);
    }
  }

#ifndef _OPENMP
  if (print_for_debugging && (my_rank == 0))
    cout << "[" << my_rank << "] bp_20" << endl;
  MPI_Comm_free(&comm_patch);
  MPI_Group_free(&group_patch);
  MPI_Comm_free(&group_comm);
  if (print_for_debugging && (my_rank == 0))
    cout << "[" << my_rank << "] bp_21" << endl;
#endif
}

void full_im_optimization_LS(
    SpEOReport &report, SpEOMatrixD &EndmemberMat, SpEOMatrixD *&AbundanceMat,
    SpEOGlobalParams *glPrms, SpEOFusionSetting *fSetting,
    SpEOSolverSetting *sSetting, SpEODataIOSetting *dSetting, SpEODataset *ImX,
    SpEODataset *ImY, SpEODataset *ImZ, SpEOMatrixD *SRF, MPI_Comm comm_busy,
    MPI_Group group_busy) {
  int my_rank, my_processes;
  MPI_Comm_rank(comm_busy, &my_rank);
  MPI_Comm_size(comm_busy, &my_processes);
  int iChX, iChY;
  int iter;        // number of iterations needed
  double rel_res;  // relative residual after optimizations
  int filter_size = 2 * glPrms->fDS - glPrms->fDS % 2;
  int fDS = (int)glPrms->fDS;
  int N_Y = glPrms->NChY;
  int N_X = glPrms->NChX;
  double lambda_X = fSetting->lambdaX_im;
  double lambda_Y = fSetting->lambdaY_im;
  int maxiter = sSetting->maxiter_CGLS_im;
  double tol_r = sSetting->tol_r_CGLS_im;

  SpEOMatrixD *I_Z = new SpEOMatrixD[glPrms->NChY];  // final image

  SpEOMatrixD *I_Z_tilde =
      new SpEOMatrixD[glPrms->NChY];  // preliminary reconstruction result
  SpEOMatrixD *I_Y =
      new SpEOMatrixD[glPrms->NChY];  // low resolution input image
  for (iChY = 0; iChY < glPrms->NChY; iChY++) {
    I_Z[iChY] = SpEOMatrixD::Zero(glPrms->sizeUH_red, glPrms->sizeVH_red);
    I_Z_tilde[iChY] = ImZ->get_rasterBands()[iChY]->bandDataMatD;
    I_Y[iChY] =
        (ImY->get_rasterBands()[iChY]->get_bandDataMat())->cast<double>();
  }
  // possibly correct ImX to selected spectral and spatial subset
  cutRelevantInput(HR, imFlag_X, ImX, dSetting, glPrms, false);
  SpEOMatrixD *I_X = new SpEOMatrixD[glPrms->NChX];
  for (iChX = 0; iChX < glPrms->NChX; iChX++) {
    I_X[iChX] =
        (ImX->get_rasterBands()[iChX]->get_bandDataMat())->cast<double>();
  }

  //******************************************************
  // Calculate the residual of the term (R*Z - X)
  //******************************************************
  SpEOMatrixD *RZ_error = new SpEOMatrixD[glPrms->NChX];
  for (iChX = 0; iChX < glPrms->NChX; iChX++) {
    RZ_error[iChX] = SpEOMatrixD::Zero(I_X[iChX].rows(), I_X[iChX].cols());
    for (iChY = 0; iChY < glPrms->NChY; iChY++) {
      RZ_error[iChX] += SRF->coeff(iChX, iChY) * I_Z_tilde[iChY];
    }
    RZ_error[iChX] -=
        (ImX->get_rasterBands()[iChX]->get_bandDataMat()->cast<double>());
  }
  SpEOMatrixD *RZ_error_rel = new SpEOMatrixD[glPrms->NChX];
  for (iChX = 0; iChX < glPrms->NChX; iChX++) {
    RZ_error_rel[iChX] = RZ_error[iChX].cwiseQuotient(
        (ImX->get_rasterBands()[iChX]->get_bandDataMat())->cast<double>());
  }

  //**************************************************
  // Check input of LS solver for inf and NaN values
  //**************************************************
  if (contains_inf_or_nan(*SRF)) {
    cerr << "[" << my_rank << "] ERROR: SRF contains 'inf' of 'NaN' values!! "
         << endl
         << endl;
    cout << "SRF = " << endl << SRF << endl << endl;
  }
  for (iChX = 0; iChX < N_X; iChX++) {
    if (contains_inf_or_nan(I_X[iChX])) {
      cerr << "[" << my_rank << "] ERROR: I_X[iChX=" << iChX
           << "] contains 'inf' of 'NaN' values!! " << endl
           << endl;
      cout << "I_X[iChX] = " << endl << I_X[iChX] << endl << endl;
    }
  }
  for (iChY = 0; iChY < N_Y; iChY++) {
    if (contains_inf_or_nan(I_Y[iChY])) {
      cerr << "[" << my_rank << "] ERROR: I_Y[iChY=" << iChY
           << "] contains 'inf' of 'NaN' values!! " << endl
           << endl;
      cout << "I_Y[iChY] = " << endl << I_Y[iChY] << endl << endl;
    }
    if (contains_inf_or_nan(I_Z_tilde[iChY])) {
      cerr << "[" << my_rank << "] ERROR: I_Z_tilde[iChY=" << iChY
           << "] contains 'inf' of 'NaN' values!! " << endl
           << endl;
      cout << "I_Z_tilde[iChY] = " << endl << I_Z_tilde[iChY] << endl << endl;
    }
    if (contains_inf_or_nan(I_Z_tilde[iChY])) {
      cerr << "[" << my_rank << "] ERROR: I_Z_tilde[iChY=" << iChY
           << "] contains 'inf' of 'NaN' values!! " << endl
           << endl;
      cout << "I_Z_tilde[iChY] = " << endl << I_Z_tilde[iChY] << endl << endl;
    }
  }
  // MPI_Barrier(MPI_COMM_WORLD);
  iChX = 0;
  if (my_rank == 0) {
    cout << "########## Full Image Optimization ######### =======>" << endl
         << "lamndaY_im = " << fSetting->lambdaY_im << endl
         << "lambdaX_im = " << fSetting->lambdaX_im << endl;
  }
  // calculate normalization factor (without the corresponding lambda value) of
  // this term
  int N_h = I_X[0].cols() * I_X[0].rows();
  double coeff_RZ = 1. / ((double)(2 * N_X * N_h));
  double RZ_res_tmp = 0;
  for (iChX = 0; iChX < N_X; iChX++) {
    RZ_res_tmp += (RZ_error[iChX].norm()) * (RZ_error[iChX].norm());
  }
  double RZ_res = sqrt(RZ_res_tmp);
  RZ_res *= sqrt(coeff_RZ);

  //******************************************************
  // Calculate the residual of the term (Z*B*S - Y)
  //******************************************************
  SpEOMatrixD gauss_filter;
  create_Gaussian_filter(gauss_filter, filter_size);
  SpEOMatrixD filter_coeff;
  calc_filter_boundary_coeff(filter_coeff, gauss_filter, fDS);

  double ZBS_res_tmp = 0;
  for (iChY = 0; iChY < N_Y; iChY++) {
    SpEOMatrixD tmp_image_HR_fast = I_Z_tilde[iChY];
    SpEOMatrixD tmp_image_LR_fast;
    fast_filter(tmp_image_LR_fast, tmp_image_HR_fast, gauss_filter,
                filter_coeff, fDS);  // fast method (intuitive)
    SpEOMatrixD ZBS_error =
        tmp_image_LR_fast -
        (ImY->get_rasterBands()[iChY]->get_bandDataMat())->cast<double>();
    SpEOMatrixD ZBS_error_rel = ZBS_error.cwiseQuotient(
        (ImY->get_rasterBands()[iChY]->get_bandDataMat())->cast<double>());
    ZBS_res_tmp += ZBS_error.norm() * ZBS_error.norm();
  }
  double ZBS_res = sqrt(ZBS_res_tmp);
  // calculate normalization factor (without the corresponding lambda value) of
  // this term
  int N_l = I_Y[0].cols() * I_Y[0].rows();
  double coeff_ZBS = 1. / ((double)(2 * N_Y * N_l));
  ZBS_res *= sqrt(coeff_ZBS);
  if (my_rank == 0) {
    cout << "## Residuals before L2-optimization: " << endl
         << "coeff1*|I_Z - I_Z_int|_2 = " << 0.0 << endl
         << "coeff2*|R*I_Z - I_X|_F^2   = " << pow(RZ_res, 2) << endl
         << "coeff3*|I_Z*B*S - I_Y|_F^2 = " << pow(ZBS_res, 2) << endl
         << "overall residual of entire functional (where coeff2 and coeff3 "
            "are multiplied by their corresp. lambda): res = "
         << lambda_X * pow(RZ_res, 2) + lambda_Y * pow(ZBS_res, 2) << endl;
  }
  // = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
  //									call least squares
  //solver
  if (my_rank == 0) {
    cout << "start LS-optimization..";
  }
  double tLS_opt_start = MPI_Wtime();
  MPI_Group eq1_group;
  int numProcEq1 = N_Y;
  if (numProcEq1 > my_processes || fSetting->fullImOptOnSubspace) {
    numProcEq1 = 1;
    if (my_rank == 0) {
      cout << "Not enough processes for parallelization of Eq. 1" << endl;
    }
  }
  int ranges_eq1[1][3] = {{0, numProcEq1 - 1, 1}};
  int rangesIdl_eq1[1][3] = {{numProcEq1, my_processes - 1, 1}};
  if (my_rank < numProcEq1) {
    MPI_Group_range_incl(group_busy, 1, ranges_eq1, &eq1_group);
  } else {
    MPI_Group_range_incl(group_busy, 1, rangesIdl_eq1, &eq1_group);
  }
  MPI_Comm eq1_comm;
  MPI_Comm_create(comm_busy, eq1_group, &eq1_comm);
  bool SNR_normalization = fSetting->SNR_normalization;
  bool balance_ImX_term_coef = fSetting->balance_ImX_term_coef;
  if (my_rank < numProcEq1) {
    if (numProcEq1 == 1) {
      if (!fSetting->fullImOptOnSubspace) {
        // calculate LAMBDA_X LAMBDA_Y and LAMBDA_Z  (can be oursource to
        // main.cpp)
        SpEOVectorD LAMBDA_X(N_X);
        SpEOVectorD LAMBDA_Y(N_Y);
        SpEOVectorD LAMBDA_Z(N_Y);
        SpEOVectorD SNR_X_db =
            35.0 *
            SpEOVectorD::Ones(N_X);  // estimated band-wise from the image ImX
        SpEOVectorD SNR_Y_db =
            35.0 *
            SpEOVectorD::Ones(N_Y);  // estimated band-wise from the image ImY
        SpEOVectorD SNR_Z_db =
            35.0 * SpEOVectorD::Ones(
                       N_Y);  // estimated band-wise from the image ImZ_init
        if (SNR_normalization) {
          for (iChX = 0; iChX < N_X; iChX++) {
            LAMBDA_X(iChX) =
                sqrt(pow(10.0, 0.1 * SNR_X_db(iChX))) / (I_X[iChX].norm());
          }
          for (iChY = 0; iChY < N_Y; iChY++) {
            LAMBDA_Y(iChY) =
                sqrt(pow(10.0, 0.1 * SNR_Y_db(iChY))) / (I_Y[iChY].norm());
          }
          for (iChY = 0; iChY < N_Y; iChY++) {
            LAMBDA_Z(iChY) = sqrt(pow(10.0, 0.1 * SNR_Z_db(iChY))) /
                             (I_Z_tilde[iChY].norm());
          }

          // Z_tilde = LAMBDA_Z * Z
          // Y_tilde = LAMBDA_Y * Y
          for (iChY = 0; iChY < N_Y; iChY++) {
            I_Z_tilde[iChY] *= LAMBDA_Z(iChY);
            I_Y[iChY] *= LAMBDA_Y(iChY);
          }
          // X_tilde = LAMBDA_X * X
          // R_tilde = LAMBDA_X * R * (LAMBDA_Y^-1)
          for (iChX = 0; iChX < N_X; iChX++) {
            I_X[iChX] *= LAMBDA_X(iChX);
            SRF->row(iChX) *= LAMBDA_X(iChX);
          }
          for (iChY = 0; iChY < N_Y; iChY++) {
            SRF->col(iChY) /= LAMBDA_Y(iChY);
          }
        }
        // ###################
        solve_equation1(iter, rel_res, I_Z, *SRF, filter_size, fDS, N_Y, N_X,
                        I_Z_tilde, I_X, I_Y, lambda_X, lambda_Y, maxiter, tol_r,
                        eq1_comm, N_Y, SNR_normalization,
                        balance_ImX_term_coef);
        // ###################
        if (SNR_normalization) {
          // Z = (LAMBDA_Y^-1) * Z_tilde
          // Y = (LAMBDA_Y^-1) * Y_tilde
          for (iChY = 0; iChY < N_Y; iChY++) {
            I_Z_tilde[iChY] /= LAMBDA_Z(iChY);
            I_Z[iChY] /= LAMBDA_Z(iChY);
            I_Y[iChY] /= LAMBDA_Y(iChY);
          }
          // X = (LAMBDA_X^-1) * X_tilde
          // R = (LAMBDA_X^-1) * R_tilde * LAMBDA_Y
          for (iChX = 0; iChX < N_X; iChX++) {
            I_X[iChX] /= LAMBDA_X(iChX);
            SRF->row(iChX) /= LAMBDA_X(iChX);
          }
          for (iChY = 0; iChY < N_Y; iChY++) {
            SRF->col(iChY) *= LAMBDA_Y(iChY);
          }
        }
      } else {
        // calculate LAMBDA_X LAMBDA_Y and LAMBDA_Z   TBD: calculation of
        // LAMBDA_X and LAMBDA_Y can be oursource to main.cpp
        SpEOVectorD LAMBDA_X(N_X);
        SpEOVectorD LAMBDA_Y(N_Y);
        SpEOVectorD LAMBDA_Z(N_Y);
        SpEOVectorD LAMBDA_ZY(N_Y);
        SpEOVectorD SNR_X_db =
            35.0 *
            SpEOVectorD::Ones(N_X);  // estimated band-wise from the image ImX
        SpEOVectorD SNR_Y_db =
            35.0 *
            SpEOVectorD::Ones(N_Y);  // estimated band-wise from the image ImY
        SpEOVectorD SNR_Z_db =
            35.0 * SpEOVectorD::Ones(
                       N_Y);  // estimated band-wise from the image ImZ_init

        for (iChX = 0; iChX < N_X; iChX++) {
          LAMBDA_X(iChX) =
              sqrt(pow(10.0, 0.1 * SNR_X_db(iChX))) / (I_X[iChX].norm());
        }
        for (iChY = 0; iChY < N_Y; iChY++) {
          LAMBDA_Y(iChY) =
              sqrt(pow(10.0, 0.1 * SNR_Y_db(iChY))) / (I_Y[iChY].norm());
          LAMBDA_Z(iChY) =
              sqrt(pow(10.0, 0.1 * SNR_Z_db(iChY))) / (I_Z_tilde[iChY].norm());
          LAMBDA_ZY(iChY) = LAMBDA_Z(iChY) / LAMBDA_Y(iChY);
        }

        if (SNR_normalization) {
          // T_tilde = LAMBDA_Y * T
          // Z_tilde = LAMBDA_Z * Z
          // Y_tilde = LAMBDA_Y * Y
          for (iChY = 0; iChY < N_Y; iChY++) {
            EndmemberMat.row(iChY) *= LAMBDA_Y(iChY);
            I_Z_tilde[iChY] *= LAMBDA_Z(iChY);
            I_Y[iChY] *= LAMBDA_Y(iChY);
          }
          // X_tilde = LAMBDA_X * X
          // R_tilde = LAMBDA_X * R * (LAMBDA_Y^-1)
          for (iChX = 0; iChX < N_X; iChX++) {
            I_X[iChX] *= LAMBDA_X(iChX);
            SRF->row(iChX) *= LAMBDA_X(iChX);
          }
          for (iChY = 0; iChY < N_Y; iChY++) {
            SRF->col(iChY) /= LAMBDA_Y(iChY);
          }
        }
        solve_equation1_unmixing(
            fSetting->use_init_value_Eq1Unmixing, iter, rel_res, EndmemberMat,
            AbundanceMat, I_Z, *SRF, filter_size, fDS, N_Y, N_X, I_Z_tilde, I_X,
            I_Y, lambda_X, lambda_Y, maxiter, tol_r, eq1_comm, N_Y,
            SNR_normalization, balance_ImX_term_coef, LAMBDA_ZY);
        fSetting->use_init_value_Eq1Unmixing =
            true;  // for any iteration after the first one, use the last result
                   // of AbundanceMat as initial value for next iteration

        if (SNR_normalization) {
          // T = (LAMBDA_Y^-1) * T_tilde
          // Z = (LAMBDA_Y^-1) * Z_tilde
          // Y = (LAMBDA_Y^-1) * Y_tilde
          for (iChY = 0; iChY < N_Y; iChY++) {
            EndmemberMat.row(iChY) /= LAMBDA_Y(iChY);
            I_Z_tilde[iChY] /= LAMBDA_Z(iChY);
            I_Z[iChY] /= LAMBDA_Y(iChY);
            I_Y[iChY] /= LAMBDA_Y(iChY);
          }
          // X = (LAMBDA_X^-1) * X_tilde
          // R = (LAMBDA_X^-1) * R_tilde * LAMBDA_Y
          for (iChX = 0; iChX < N_X; iChX++) {
            I_X[iChX] /= LAMBDA_X(iChX);
            SRF->row(iChX) /= LAMBDA_X(iChX);
          }
          for (iChY = 0; iChY < N_Y; iChY++) {
            SRF->col(iChY) *= LAMBDA_Y(iChY);
          }
        }
      }
    } else {
      SpEOMatrixD *I_Ztmp = new SpEOMatrixD;
      SpEOMatrixD SRF_tmp = SpEOMatrixD(SRF->rows(), 1);
      SRF_tmp = SRF->col(my_rank);
      solve_equation1(iter, rel_res, I_Ztmp, SRF_tmp, filter_size, fDS, N_Y,
                      N_X, &(I_Z_tilde[my_rank]), I_X, &(I_Y[my_rank]),
                      lambda_X, lambda_Y, maxiter, tol_r, eq1_comm, 1,
                      SNR_normalization, balance_ImX_term_coef);
      I_Z[my_rank] = *(I_Ztmp);
      for (iChY = 0; iChY < N_Y; iChY++) {
        MPI_Allreduce(MPI_IN_PLACE, I_Z[iChY].data(),
                      glPrms->sizeUH_red * glPrms->sizeVH_red, MPI_DOUBLE,
                      MPI_SUM, eq1_comm);
      }
      delete I_Ztmp;
    }
  }

  MPI_Barrier(comm_busy);
  MPI_Comm_free(&eq1_comm);
  MPI_Group_free(&eq1_group);
  if (my_rank == 0) {
    cout << ".. LS-optimization finished!" << endl;
    cout << "took " << MPI_Wtime() - tLS_opt_start << " sec. " << endl;
  }
  if (my_rank == 0) {
    cout << "Solver output: iter=" << iter << ", rel_res=" << rel_res << endl;
    report.file.open(report.fileName.c_str(),
                     fstream::in | fstream::out | fstream::app);
    report.file << "\n Eq. 1 Solver output: iter=" << iter
                << ", rel_res=" << rel_res << "\n\n";
    report.file.close();
  }

  //*******************************
  //                             //
  //  re-calculate residuals     //
  //                             //
  //******************************
  //******************************************************
  // Calculate the residual of the term (R*Z - X)
  //******************************************************
  // residual IZ_res
  double IZ_res_tmp = 0;
  for (iChY = 0; iChY < N_Y; iChY++) {
    IZ_res_tmp += (I_Z[iChY] - I_Z_tilde[iChY]).norm() *
                  (I_Z[iChY] - I_Z_tilde[iChY]).norm();
  }
  double IZ_res = sqrt(IZ_res_tmp);
  double coeff_IZ = 1. / ((double)(2 * N_Y * N_h));
  IZ_res *= sqrt(coeff_IZ);
  //******************************************************
  // Calculate the residual of the term (R*Z - X)
  //******************************************************
  for (iChX = 0; iChX < glPrms->NChX; iChX++) {
    RZ_error[iChX] = SpEOMatrixD::Zero(I_X[iChX].rows(), I_X[iChX].cols());
    for (iChY = 0; iChY < glPrms->NChY; iChY++) {
      RZ_error[iChX] += SRF->coeff(iChX, iChY) * I_Z[iChY];
    }
    RZ_error[iChX] -=
        (ImX->get_rasterBands()[iChX]->get_bandDataMat()->cast<double>());
  }
  for (iChX = 0; iChX < glPrms->NChX; iChX++) {
    RZ_error_rel[iChX] = RZ_error[iChX].cwiseQuotient(
        (ImX->get_rasterBands()[iChX]->get_bandDataMat())->cast<double>());
  }
  RZ_res_tmp = 0;
  for (iChX = 0; iChX < N_X; iChX++) {
    RZ_res_tmp += (RZ_error[iChX].norm()) * (RZ_error[iChX].norm());
  }
  RZ_res = sqrt(RZ_res_tmp);
  RZ_res *= sqrt(coeff_RZ);
  delete[] RZ_error;
  delete[] RZ_error_rel;
  //******************************************************
  // Calculate the residual of the term (Z*B*S - Y)
  //******************************************************
  ZBS_res_tmp = 0;
  for (iChY = 0; iChY < N_Y; iChY++) {
    SpEOMatrixD tmp_image_HR_fast = I_Z[iChY];
    SpEOMatrixD tmp_image_LR_fast;
    fast_filter(tmp_image_LR_fast, tmp_image_HR_fast, gauss_filter,
                filter_coeff, fDS);  // fast method (intuitive)

    SpEOMatrixD ZBS_error =
        tmp_image_LR_fast -
        (ImY->get_rasterBands()[iChY]->get_bandDataMat())->cast<double>();
    SpEOMatrixD ZBS_error_rel = ZBS_error.cwiseQuotient(
        (ImY->get_rasterBands()[iChY]->get_bandDataMat())->cast<double>());

    ZBS_res_tmp += ZBS_error.norm() * ZBS_error.norm();
  }
  ZBS_res = sqrt(ZBS_res_tmp);
  ZBS_res *= sqrt(coeff_ZBS);

  if (my_rank == 0) {
    cout << "## Residuals after L2-optimization: " << endl
         << "coeff1*|I_Z - I_Z_int|_F^2 = " << pow(IZ_res, 2) << endl
         << "coeff2*|R*I_Z - I_X|_F^2   = " << pow(RZ_res, 2) << endl
         << "coeff3*|I_Z*B*S - I_Y|_F^2 = " << pow(ZBS_res, 2) << endl
         << "overall residual of entire functional (where coeff2 and coeff3 "
            "are multiplied by their corresp. lambda): res = "
         << pow(IZ_res, 2) + lambda_X * pow(RZ_res, 2) +
                lambda_Y * pow(ZBS_res, 2)
         << endl
         << "############################################ <=======" << endl;
  }

  //**************************************************
  // Check output of LS solver for inf and NaN values
  //**************************************************
  for (iChY = 0; iChY < N_Y; iChY++) {
    if (contains_inf_or_nan(I_Z[iChY])) {
      cerr << "[" << my_rank << "] ERROR: I_Z[iChY=" << iChY
           << "] contains 'inf' of 'NaN' values!! " << endl
           << endl;
      cout << "I_Z[iChY] = " << endl << I_Z[iChY] << endl << endl;
    }
  }
  // copy I_Z to ImZ while setting all non-negative values to 0
  if (fSetting->set_neg_to_0 == 2 || fSetting->set_neg_to_0 == 3) {
    for (iChY = 0; iChY < glPrms->NChY; iChY++) {
      ImZ->rasterBands[iChY]->bandDataMatD =
          ((I_Z[iChY]).array() > 0).select(I_Z[iChY], 0);
    }
  } else {  // copy I_Z to ImZ without setting all non-negative values to 0
    for (iChY = 0; iChY < glPrms->NChY; iChY++) {
      ImZ->rasterBands[iChY]->bandDataMatD = I_Z[iChY];
    }
  }
  delete[] I_Z;
  delete[] I_Z_tilde;
  delete[] I_Y;
  delete[] I_X;
}
