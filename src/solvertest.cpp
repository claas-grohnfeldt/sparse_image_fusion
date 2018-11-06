/*#############################################################################
##                                                                           ##
##    File:        solvertest.cpp                                            ##
##                                                                           ##
##    Author:      Dipl.-Math.techn. Steffen PETER (since 10/2015)           ##
##                                                                           ##
##    Contact / corresponding author:                                        ##
##                 Steffen Peter                                             ##
##                 Technical University Munich                               ##
##                 Institute:  Mathematics                                   ##
##                 Department: Applied Numerical Analysis                    ##
##                 Location:   Garching b. München, Germany                  ##
##                 E-mail:     Steffen.Peter@ma.tum.de                       ##
##                                                                           ##
##    Version:     0.1 (since 10/2015)                                       ##
##                                                                           ##
##    Last modification: October 20, 2015                                    ##
##                                                                           ##
##    Security Classification: Confidential                                  ##
##                                                                           ##
##    Copyright:   [2013] - [2015] DLR. All Rights Reserved.                 ##
##                                                                           ##
##                 Notice: All information contained herein is, and remains  ##
##                 the property of DLR and its suppliers, if any. The        ##
##                 intellectual and technical concepts contained herein are  ##
##                 proprietary to DLR and its suppliers and may be covered   ##
##                 by German and Foreign Patents, patents in process, and    ##
##                 are protected by trade secret or copyright law.           ##
##                 Dissemination of this information or reproduction of this ##
##                 material is strictly forbidden unless prior written       ##
##                 permission is obtained from DLR.                          ##
##                                                                           ##
##    Description: FBS Solver testing by means of standard gaussian matrices ##
##                                                                           ##
##    Input:       ---                                                       ##
##                                                                           ##
##    Output:      ---                                                       ##
##                                                                           ##
#############################################################################*/

#include "solvertest.h"

using namespace std;

void write_var_to_matlab(std::ofstream &in_file_matlab, SpEOMatrixD A, string varname) {
  in_file_matlab << varname << " = [";
  for (int i = 0; i < A.rows(); i++) {
    for (int j = 0; j < A.cols(); j++) {
      in_file_matlab << A(i,j);
      if (j < A.cols()-1) {
	in_file_matlab << ",";
      }
    }
    if (i < A.rows()-1) {
      in_file_matlab << ";";
    }
    else {
      in_file_matlab << "];"<< endl << endl;
    }
  }
  return;
}

void write_var_to_matlab(std::ofstream &in_file_matlab, SpEOVectorD A, string varname) {
  in_file_matlab << varname << " = [";
  for (int i = 0; i < A.size(); i++) {
    in_file_matlab << A(i);
    if (i < A.size()-1) {
      in_file_matlab << ";";
    }
    else {
      in_file_matlab << "];"<< endl << endl;
    }
  }
  return;
}

void write_var_to_matlab(std::ofstream &in_file_matlab, SpEOVectorI A, string varname) {
  in_file_matlab << varname << " = [";
  for (int i = 0; i < A.size(); i++) {
    in_file_matlab << A(i);
    if (i < A.size()-1) {
      in_file_matlab << ";";
    }
    else {
      in_file_matlab << "];"<< endl << endl;
    }
  }
  return;
}

int main(int argc, char* argv[]) {

  //==================================================================//
  //                        MPI initialization                        //
  //==================================================================//
#ifdef _OPENMP
  int iprovided;
  MPI_Init_thread(&argc,&argv,MPI_THREAD_FUNNELED,&iprovided);
#else
  MPI_Init(&argc,&argv);
#endif
  //MPI_Status recv_status;
  int my_rank; int my_processes;
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &my_processes);

  //==================================================================//
  //                        READ ARGUMENTS / SET GLOBAL VARS          //
  //==================================================================//
  
  string curTime = get_current_time();
      
  // user arguments
  int m = atoi(argv[1]);
  int N = atoi(argv[2]);
  int k = atoi(argv[3]);
  int d = atoi(argv[4]);

  double sigma_y = (double) atof(argv[5]);
  double sigma_x = (double) atof(argv[6]);

  double alpha = (double) atof(argv[7]);

  double switching_gamma = (double) atof(argv[8]);
  int maxiterIn = atoi(argv[9]);
  int maxiterOut = atoi(argv[10]);
  double stopCritTol = (double) atof(argv[11]);
  
  int testN = atoi(argv[12]);
 
  // global variables
  srand(time(NULL));
  bool singletest = false;
  bool write_to_matlab = false;
  bool write_output_matlab = false;
  bool read_example = false;
  SpEOVectorI Ndecomp_def = SpEOVectorI::LinSpaced(8,(int) log2(my_processes)-7,(int) log2(my_processes));
  SpEOVectorI maxiter_in_def = SpEOVectorI::LinSpaced(maxiterIn,1,maxiterIn);
  
  int methods = 3;
  int methodsS = 2;
  Stopping_criterion stop_rule = FIRST_ORDER_CONDITIONS;
  
  FBSsolverOptions* opts = new FBSsolverOptions[methods];
  opts[0] = FBSsolverOptions(
			    1,  // decomposition parameter
			    2, // max. inner iterations (1 for standard FISTA)
			    maxiterOut, // max. outer iterations (former optCoeffResPFISTA_maxiter)
			    stopCritTol, // respective tolerance (dep. on stopping crit.) (former optCoeffResPFISTA_tol_r)
			    stop_rule, // stopping criterion
			    LINE_SEARCH_SIMPLEX, // prediction step rule
			    NO, // backtracking strategy NO, INC, or DEC
			    INIT,
			    false, //adaptive
			    true,
			    switching_gamma); 
  opts[1] = FBSsolverOptions(
    			    3,  // decomposition parameter
			    1, // max. inner iterations (1 for standard FISTA)
			    maxiterOut, // max. outer iterations (former optCoeffResPFISTA_maxiter)
			    stopCritTol, // respective tolerance (dep. on stopping crit.) (former optCoeffResPFISTA_tol_r)
			    stop_rule, // stopping criterion
			    LINE_SEARCH_SIMPLEX, // prediction step rule
			    DECEFF, // backtracking strategy NO, INC, or DEC // DECNO
			    INIT,
			    true,//adaptive
			    true,
			    0.1*switching_gamma);
  opts[2] = FBSsolverOptions(
    			    3,  // decomposition parameter
			    2, // max. inner iterations (1 for standard FISTA)
			    maxiterOut, // max. outer iterations (former optCoeffResPFISTA_maxiter)
			    stopCritTol, // respective tolerance (dep. on stopping crit.) (former optCoeffResPFISTA_tol_r)
			    stop_rule, // stopping criterion
			    LINE_SEARCH_SIMPLEX, // prediction step rule
			    NO, // backtracking strategy NO, INC, or DEC // NO
			    NORM,
			    false,//adaptive
			    true,
			    switching_gamma);
  FBSsolverOptions* optsS = new FBSsolverOptions[methodsS];
  optsS[0] = FBSsolverOptions(
    			    1,  // decomposition parameter
			    1, // max. inner iterations (1 for standard FISTA)
			    5*maxiterOut, // max. outer iterations
			    stopCritTol, // respective tolerance (dep. on stopping crit.) (former optCoeffResPFISTA_tol_r)
			    stop_rule, // stopping criterion
			    FISTA, // prediction step rule
			    NO, // backtracking strategy NO, INC, or DEC
			    INIT,
			    false,//adaptive
			    true,
			    switching_gamma);
  optsS[1] = FBSsolverOptions(
    			    1,  // decomposition parameter
			    1, // max. inner iterations (1 for standard FISTA)
			    maxiterOut, // max. outer iterations
			    stopCritTol, // respective tolerance (dep. on stopping crit.) (former optCoeffResPFISTA_tol_r)
			    stop_rule, // stopping criterion
			    LINE_SEARCH_ARMIJO, // prediction step rule
			    NO, // backtracking strategy NO, INC, or DEC
			    PSCL,
			    false,//adaptive
			    true,
			    switching_gamma);

  // repeat process
  if (singletest) {
    testN = 1;
    Ndecomp_def = SpEOVectorI::LinSpaced(1,0,0);
    maxiter_in_def = SpEOVectorI::LinSpaced(1,1,1);
  }

  SpEOMatrixD* result_iter_data = new SpEOMatrixD[methods];
  SpEOMatrixD* result_time_total = new SpEOMatrixD[methods];
  SpEOMatrixD* result_time_commu = new SpEOMatrixD[methods];
  SpEOMatrixD* result_time_backt = new SpEOMatrixD[methods];
  SpEOMatrixD* result_time_predi = new SpEOMatrixD[methods];
  SpEOMatrixD* result_time_stopp = new SpEOMatrixD[methods];
  int meth;
  for (meth = 0; meth < methods; meth++){
    result_iter_data[meth] = SpEOMatrixD::Zero(Ndecomp_def.size(),maxiter_in_def.size());
    result_time_total[meth] = SpEOMatrixD::Zero(Ndecomp_def.size(),maxiter_in_def.size());
    result_time_commu[meth] = SpEOMatrixD::Zero(Ndecomp_def.size(),maxiter_in_def.size());
    result_time_backt[meth] = SpEOMatrixD::Zero(Ndecomp_def.size(),maxiter_in_def.size());
    result_time_predi[meth] = SpEOMatrixD::Zero(Ndecomp_def.size(),maxiter_in_def.size());
    result_time_stopp[meth] = SpEOMatrixD::Zero(Ndecomp_def.size(),maxiter_in_def.size());
  }
  SpEOVectorD* result_iter_data_single = new SpEOVectorD[methodsS];
  SpEOVectorD* result_time_total_single = new SpEOVectorD[methodsS];
  SpEOVectorD* result_time_commu_single = new SpEOVectorD[methodsS];
  SpEOVectorD* result_time_backt_single = new SpEOVectorD[methodsS];
  SpEOVectorD* result_time_predi_single = new SpEOVectorD[methodsS];
  SpEOVectorD* result_time_stopp_single = new SpEOVectorD[methodsS];
  for (meth = 0; meth < methodsS; meth++){
    result_iter_data_single[meth] = SpEOVectorD::Zero(Ndecomp_def.size());
    result_time_total_single[meth] = SpEOVectorD::Zero(Ndecomp_def.size());
    result_time_commu_single[meth] = SpEOVectorD::Zero(Ndecomp_def.size());
    result_time_backt_single[meth] = SpEOVectorD::Zero(Ndecomp_def.size());
    result_time_predi_single[meth] = SpEOVectorD::Zero(Ndecomp_def.size());
    result_time_stopp_single[meth] = SpEOVectorD::Zero(Ndecomp_def.size());
  }
  
  
  if (my_rank == 0) {
    cout << "+++++++ INITIALIZATION DONE " << endl;
  }

  // test runs
  int testnr;
  for (testnr = 1; testnr <= testN; testnr++) {
    if (my_rank == 0) {
      cout << "++++++++++++ DO TEST NR  " << testnr << endl;
    }
    // data test set creation (with broadcasting)
    SpEOMatrixD A, y;
    if (read_example) {
      int stat_CSV_read_A = read_CSV(&A, "./A.csv", ',', 0);
      int stat_CSV_read_y = read_CSV(&y, "./y.csv", ',', 0);
      N = A.cols();
      m = A.rows();
      d = y.cols();
    } 
    else {
      double * bufA_recv = new double[m*N];
      double * bufy_recv = new double[m*d];
      if (my_rank == 0) {
	// data test set creation
	A = randn(m,N,0,1);
	A = (1.0/spec_norm(A))*A;  // normalize
	SpEOVectorI I = randperm(N); SpEOVectorI S = I.head(k); // generate k random locations for the support T of x
	SpEOMatrixD xr = randn(k,d,0,1); //generate signal
	SpEOMatrixD xsol = SpEOMatrixD::Zero(N,d); 
	int k_i;
	for (k_i = 0; k_i < k; k_i++) {
	  xsol.row(S(k_i)) = xr.row(k_i);
	}
	SpEOMatrixD w = sigma_y*randn(m,d,0,1);
	SpEOMatrixD wx = sigma_x*randn(N,d,0,1);
	SpEOMatrixD x_pert = 1*(xsol + wx);
	y = A*(x_pert) + 1*w; // generate measurements
	*(bufA_recv) = *(A.data());
	*(bufy_recv) = *(y.data());
      }
      else {
	A = SpEOMatrixD::Zero(m,N);
	y = SpEOMatrixD::Zero(m,d);
      }
      MPI_Allreduce(y.data(),bufy_recv,m*d,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
      MPI_Allreduce(A.data(),bufA_recv,m*N,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
      A = SpEOMatrixD::Map(bufA_recv,m,N);
      y = SpEOMatrixD::Map(bufy_recv,m,d);
      delete [] bufA_recv;
      delete [] bufy_recv;
    }
    
    if (write_to_matlab) {    
      string fname_matlab;
      fname_matlab = "./matlab_input_A_y.m";
      std::ofstream in_file_matlab(fname_matlab.c_str());
      if (in_file_matlab.is_open()) {
	write_var_to_matlab(in_file_matlab, A, "A");
	write_var_to_matlab(in_file_matlab, y, "y");
      }
      else {
	cout << endl << "WARNING: in_matlab.m could not be written!" << endl;
      }
      in_file_matlab.close();
      chmod(fname_matlab.c_str(), S_IRWXU | S_IRWXG | S_IRWXO);
    }
    
    // testrun on all parameter pairs
    int i_Ndecomp, i_maxiter_in;
    for (i_Ndecomp = 0; i_Ndecomp < Ndecomp_def.size(); i_Ndecomp++) {
      int Ndecomp = pow(2,Ndecomp_def(i_Ndecomp));
      for (i_maxiter_in = 0; i_maxiter_in < maxiter_in_def.size(); i_maxiter_in++) {
	int maxiter_in = maxiter_in_def(i_maxiter_in);
	
	// output
	/*string fname_matlab;
	if (singletest) {
	  fname_matlab = "./matlab_input_backtrack.m";
	}
	else {
	  string nd; stringstream c_nd; c_nd << Ndecomp; nd = c_nd.str();
	  string mi; stringstream c_mi; c_mi << maxiter_in; mi = c_mi.str();
	  fname_matlab = "./matlab_input_backtrack_nd" + nd +  "_mi" +  mi + ".m";
	}
	std::ofstream in_file_matlab(fname_matlab.c_str());
	bool canwrite = false;
	if (in_file_matlab.is_open()) {
	  canwrite = true;
	  if (my_rank == 0) {
	    in_file_matlab << "bstMethod = cell(1, " << methods << ");" << endl << endl;
	  }
	}
	else {
	  cout << endl << "WARNING: in_matlab.m could not be written!" << endl;
	}*/
	
	
	// sparse reconstruction
	//double alpha_normalized = alpha*y.norm()/sqrt(N*m);
	double alpha_normalized = alpha*((A.transpose()*y).rowwise().norm().maxCoeff());
	
	for (meth = 0; meth < methods; meth++) {
	  
	  if (!singletest) {
	    opts[meth].Ndecomp = Ndecomp;
	    opts[meth].maxiter_in = maxiter_in;
	  }
	  int Ndecomp_here = opts[meth].Ndecomp;
	  
	  // create decomposition (uniform)	 
	  MPI_Group mpi_group_orig;
	  MPI_Comm_group(MPI_COMM_WORLD, &mpi_group_orig);
	  MPI_Group group_busy;
	  int numProcBusy = Ndecomp_here;
	  if (Ndecomp_here > my_processes) {
	    cout << "ERROR: NOT ENOUGH PROCESSES AVAILABLE" << endl;
	    exit(EXIT_FAILURE);
	  }
	  int ranges_busy[1][3] = {{0,numProcBusy-1,1}};
	  int rangesIdl_busy[1][3] = {{numProcBusy,my_processes-1,1}};
	  if(my_rank < numProcBusy){
		  MPI_Group_range_incl(mpi_group_orig, 1, ranges_busy, &group_busy);
	  }else{
		  MPI_Group_range_incl(mpi_group_orig, 1, rangesIdl_busy, &group_busy);
	  }
	  MPI_Comm comm_busy;
	  MPI_Comm_create(MPI_COMM_WORLD, group_busy, &comm_busy);
	  
	  if (my_rank < numProcBusy) {
	    // obtain your part of the data
	    int colstart, colN;
	    if (my_rank < N%Ndecomp_here) {
	      colstart = my_rank*(N/Ndecomp_here+1);
	      colN = N/Ndecomp_here+1;
	    }
	    else {
	      colstart = (N%Ndecomp_here)*(N/Ndecomp_here+1)+(my_rank-(N%Ndecomp_here))*(N/Ndecomp_here);
	      colN = N/Ndecomp_here;
	    }
	    SpEOMatrixD A_loc = SpEOMatrixD::Zero(m, colN);
	    int li;
	    for (li = 0; li < colN; li++) {
	      A_loc.col(li) = A.col(colstart + li);
	    }
	    /*SpEOMatrixD A_loc;
	    if (my_rank < N%Ndecomp_here) {
	      A_loc = A.middleCols(my_rank*(N/Ndecomp_here+1),N/Ndecomp_here+1);
	    }
	    else {
	      A_loc = A.middleCols((N%Ndecomp_here)*(N/Ndecomp_here+1)+(my_rank-(N%Ndecomp_here))*(N/Ndecomp_here), N/Ndecomp_here);
	    }*/
	    MPI_Barrier(comm_busy);
	    
	    // call solver
	    int iter;
	    double rel_res;
	    SpEOMatrixD x_rec;
	    SpEOMatrixD timestat;
	    SpEOMatrixD btstat;
	    if (my_rank == 0) {
	      cout << "+++++++ CALL SOLVER for Method " << meth << " with P = " << Ndecomp << ", L_max = " << maxiter_in << endl;
	    }
	    MPI_Barrier(comm_busy);
	    double time000 = MPI_Wtime();
	    if (Ndecomp_here == 1) {
	      FBSSolver(iter, rel_res, x_rec, A, y, alpha_normalized, opts[meth], timestat, btstat);
	    }
	    else {
	      FBSSolver(iter, rel_res, x_rec, N, A_loc, y, alpha_normalized, opts[meth], comm_busy, timestat, btstat);
	    }
	    MPI_Barrier(comm_busy);
	    result_time_total[meth](i_Ndecomp,i_maxiter_in) += (MPI_Wtime()-time000)/testN;
	    //result_time_total[meth](i_Ndecomp,i_maxiter_in) += timestat(0,0)/testN;
	    result_time_commu[meth](i_Ndecomp,i_maxiter_in) += timestat(0,1)/testN;
	    result_time_backt[meth](i_Ndecomp,i_maxiter_in) += timestat(0,2)/testN;
	    result_time_predi[meth](i_Ndecomp,i_maxiter_in) += timestat(0,3)/testN;
	    result_time_stopp[meth](i_Ndecomp,i_maxiter_in) += timestat(0,4)/testN;
	    result_iter_data[meth](i_Ndecomp,i_maxiter_in) += (1.0*iter)/testN;
	    /*if (canwrite) {
	      string Result;
	      stringstream convert;
	      convert << meth+1;
	      Result = convert.str();
	      MPI_Allreduce(MPI_IN_PLACE,btstat.data(),btstat.rows()*btstat.cols(),MPI_DOUBLE,MPI_SUM,comm_busy);
	      if (my_rank == 0) {
		write_var_to_matlab(in_file_matlab, btstat, "bstMETHOD{" + Result + "}");
	      }
	    }*/
	  }
	  MPI_Barrier(MPI_COMM_WORLD);
	  MPI_Comm_free(&comm_busy);
	  MPI_Group_free(&group_busy);
	  MPI_Group_free(&mpi_group_orig);	  
	}
	 
	if (maxiter_in == 1) {
	  for (meth = 0; meth < methodsS; meth++) {
	    if (!singletest) {
	      optsS[meth].Ndecomp = Ndecomp;
	    }
	    int Ndecomp_here = optsS[meth].Ndecomp;	    
	    // create decomposition (uniform)	  
	    MPI_Group mpi_group_orig;
	    MPI_Comm_group(MPI_COMM_WORLD, &mpi_group_orig);
	    MPI_Group group_busy;
	    int numProcBusy = Ndecomp_here;
	    if (Ndecomp_here > my_processes) {
	      cout << "ERROR: NOT ENOUGH PROCESSES AVAILABLE" << endl;
	      exit(EXIT_FAILURE);
	    }
	    int ranges_busy[1][3] = {{0,numProcBusy-1,1}};
	    int rangesIdl_busy[1][3] = {{numProcBusy,my_processes-1,1}};
	    if(my_rank < numProcBusy){
		    MPI_Group_range_incl(mpi_group_orig, 1, ranges_busy, &group_busy);
	    }else{
		    MPI_Group_range_incl(mpi_group_orig, 1, rangesIdl_busy, &group_busy);
	    }
	    MPI_Comm comm_busy;
	    MPI_Comm_create(MPI_COMM_WORLD, group_busy, &comm_busy);
	    if (my_rank < numProcBusy) {
	      int colstart, colN;
	      if (my_rank < N%Ndecomp_here) {
		colstart = my_rank*(N/Ndecomp_here+1);
		colN = N/Ndecomp_here+1;
	      }
	      else {
		colstart = (N%Ndecomp_here)*(N/Ndecomp_here+1)+(my_rank-(N%Ndecomp_here))*(N/Ndecomp_here);
		colN = N/Ndecomp_here;
	      }
	      SpEOMatrixD A_loc = SpEOMatrixD::Zero(m, colN);
	      int li;
	      for (li = 0; li < colN; li++) {
		A_loc.col(li) = A.col(colstart + li);
	      }
	      /*SpEOMatrixD A_loc;
	      if (my_rank < N%Ndecomp_here) {
		A_loc = A.middleCols(my_rank*(N/Ndecomp_here+1),N/Ndecomp_here+1);
	      }
	      else {
		A_loc = A.middleCols((N%Ndecomp_here)*(N/Ndecomp_here+1)+(my_rank-(N%Ndecomp_here))*(N/Ndecomp_here), N/Ndecomp_here);
	      }*/
	      int iter;
	      double rel_res;
	      SpEOMatrixD x_rec;
	      SpEOMatrixD timestat;
	      SpEOMatrixD btstat;
	      if (my_rank == 0) {
		cout << "+++++++ CALL SOLVER for Method Single " << meth << " with P = " << Ndecomp <<  endl;
	      }
	      MPI_Barrier(comm_busy);
	      double time000 = MPI_Wtime();
	      if (Ndecomp_here == 1) {
		FBSSolver(iter, rel_res, x_rec, A, y, alpha_normalized, optsS[meth], timestat, btstat);
	      }
	      else {
		FBSSolver(iter, rel_res, x_rec, N, A_loc, y, alpha_normalized, optsS[meth], comm_busy, timestat, btstat);
	      }
	      result_time_total_single[meth](i_Ndecomp) += (MPI_Wtime()-time000)/testN;
	      //result_time_total_fista(i_Ndecomp) += timestat(0,0)/testN;
	      result_time_commu_single[meth](i_Ndecomp) += timestat(0,1)/testN;
	      result_time_backt_single[meth](i_Ndecomp) += timestat(0,2)/testN;
	      result_time_predi_single[meth](i_Ndecomp) += timestat(0,3)/testN;
	      result_time_stopp_single[meth](i_Ndecomp) += timestat(0,4)/testN;
	      result_iter_data_single[meth](i_Ndecomp) += (1.0*iter)/testN;
	      /*if (canwrite) {
		MPI_Allreduce(MPI_IN_PLACE,btstat.data(),btstat.rows()*btstat.cols(),MPI_DOUBLE,MPI_SUM,comm_busy);
		if (my_rank == 0) {
		  write_var_to_matlab(in_file_matlab, btstat, "bstFISTA");
		}
	      }*/
	    }
	    //cout << "++++++++++++ FISTA DONE " << endl;
	    
	    // output close
	    /*MPI_Barrier(comm_busy);
	    in_file_matlab.close();
	    chmod(fname_matlab.c_str(), S_IRWXU | S_IRWXG | S_IRWXO);*/
	    
	    // free structures
	    MPI_Barrier(MPI_COMM_WORLD);
	    MPI_Comm_free(&comm_busy);
	    MPI_Group_free(&group_busy);
	    MPI_Group_free(&mpi_group_orig);
	  }
	}
      }
    }
  }
  MPI_Barrier(MPI_COMM_WORLD);
  
  // Write Report
  if (my_rank == 0) {
    // open output file
    string mydate = curTime.substr(0, 6);
    string mytime = curTime.substr(7, 6);
    string dir_out = "./recResults/res_at_" + mydate + "_" + mytime;
    mkdir(dir_out.c_str(), 0777);
    chmod(dir_out.c_str(), S_IRWXU | S_IRWXG | S_IRWXO);  
    string reportFileName = dir_out + "/report.txt";
    string reportFileNameM = dir_out + "/report.m";
    ofstream reportFile, reportFileM;
    reportFile.open(reportFileName.c_str(), fstream::in | fstream::out | fstream::app);
    reportFileM.open(reportFileNameM.c_str(), fstream::in | fstream::out | fstream::app); 
    
    // write parameters
    reportFile << "Global Parameters " << endl;
    reportFile << "---------------------------------------------------------------------" << endl;
    reportFile << "m = " << m << endl;
    reportFile << "N = " << N << endl;
    reportFile << "k = " << k << endl;
    reportFile << "d = " << d << endl;
    reportFile << "sigma_y = " << sigma_y << endl;
    reportFile << "sigma_x = " << sigma_x << endl;
    reportFile << "alpha = " << alpha << endl;
    reportFile << "gamma (switching parameter) (glob.) = " << switching_gamma << endl;
    reportFile << "maximum inner iterations = " << maxiterIn << endl;
    reportFile << "maximum outer iterations = " << maxiterOut << endl;
    reportFile << "stopping criterion = " << stop_rule << endl;
    reportFile << "stopping criterion tolerance = " << stopCritTol << endl;
    reportFile << "number of test samples = " << testN << endl;
    reportFile << "single test? = " << singletest << endl;
    reportFile << "external example? = " << read_example << endl;
    reportFile << endl << endl;
    
    reportFile << "Tested methods " << endl;
    reportFile << "---------------------------------------------------------------------" << endl;
    for (meth = 0; meth < methods; meth++){
      reportFile << "Method " << meth << endl;
      if (singletest) {
	reportFile << "P = " << opts[meth].Ndecomp << ", maximum inner iterations = " << opts[meth].maxiter_in << endl;
      }
      reportFile << "Prediction Rule = " << opts[meth].pred_rule << ", backtracking strategy = " << opts[meth].backtracking << ", inner stepsize strategy = " << opts[meth].stepsize_update << ", adaptive stepsize? = " << opts[meth].stepsize_adaptive << endl;
      reportFile << "gamma (switching parameter) = " << opts[meth].gamma << endl;
    }
    for (meth = 0; meth < methodsS; meth++){
      reportFile << "Method SINGLE " << meth << endl;
      if (singletest) {
	reportFile << "P = " << optsS[meth].Ndecomp << ", maximum inner iterations = " << optsS[meth].maxiter_in << endl;
      }
      reportFile << "Prediction Rule = " << optsS[meth].pred_rule << ", backtracking strategy = " << optsS[meth].backtracking << ", inner stepsize strategy = " << optsS[meth].stepsize_update << ", adaptive stepsize? = " << optsS[meth].stepsize_adaptive << endl;
      reportFile << "gamma (switching parameter) = " << optsS[meth].gamma << endl;
      reportFile << endl << endl;
    }
    
    reportFile << "Tested parameters " << endl;
    reportFile << "---------------------------------------------------------------------" << endl;
    reportFile << "P = 2^ ... " << Ndecomp_def.transpose() << endl;
    reportFile << "maximum inner iterations = " << maxiter_in_def.transpose() << endl;
    reportFile << endl << endl;
    
    reportFile << "Results " << endl;
    reportFile << "---------------------------------------------------------------------" << endl;
    for (meth = 0; meth < methods; meth++) {
      reportFile << "### Result Method " << meth << endl;
      reportFile << "Iterations:" << endl;
      reportFile << result_iter_data[meth] << endl << endl;
      reportFile << "Total Time:" << endl;
      reportFile << result_time_total[meth] << endl << endl;
      /*cout << 100*result_time_commu[meth].cwiseQuotient(result_time_total[meth]) << endl << endl;
      cout << 100*result_time_backt[meth].cwiseQuotient(result_time_total[meth]) << endl << endl;
      cout << 100*result_time_predi[meth].cwiseQuotient(result_time_total[meth]) << endl << endl;
      cout << 100*result_time_stopp[meth].cwiseQuotient(result_time_total[meth]) << endl << endl;*/
      reportFile << "Communication Time:" << endl;
      reportFile << result_time_commu[meth] << endl << endl;
      reportFile << "Working Time:" << endl;
      reportFile << result_time_backt[meth] << endl << endl;
      reportFile << "Prediction Time:" << endl;
      reportFile << result_time_predi[meth] << endl << endl;
      reportFile << "Stopping Criterion Time:" << endl;
      reportFile << result_time_stopp[meth] << endl << endl;
    }
    for (meth = 0; meth < methodsS; meth++) {
      reportFile << "### Result Method Single" << meth << endl;
      reportFile << "Iterations:" << endl;
      reportFile << result_iter_data_single[meth] << endl << endl;
      reportFile << "Total Time:" << endl;
      reportFile << result_time_total_single[meth] << endl << endl;
      /*cout << 100*result_time_commu_fista.cwiseQuotient(result_time_total_fista) << endl << endl;
      cout << 100*result_time_backt_fista.cwiseQuotient(result_time_total_fista) << endl << endl;
      cout << 100*result_time_predi_fista.cwiseQuotient(result_time_total_fista) << endl << endl;
      cout << 100*result_time_stopp_fista.cwiseQuotient(result_time_total_fista) << endl << endl;*/
      reportFile << "Communication Time:" << endl;
      reportFile << result_time_commu_single[meth] << endl << endl;
      reportFile << "Working Time:" << endl;
      reportFile << result_time_backt_single[meth] << endl << endl;
      reportFile << "Prediction Time:" << endl;
      reportFile << result_time_predi_single[meth] << endl << endl;
      reportFile << "Stopping Criterion Time:" << endl;
      reportFile << result_time_stopp_single[meth] << endl << endl;
    }
    reportFile.close();
    
    // MATLAB
    reportFileM << "m = " << m << "; " << endl;
    reportFileM << "N = " << N << "; " << endl;
    reportFileM << "k = " << k << "; " << endl;
    reportFileM << "d = " << d << "; " << endl;
    reportFileM << "sigma_y = " << sigma_y << "; " << endl;
    reportFileM << "sigma_x = " << sigma_x << "; " << endl;
    reportFileM << "alpha = " << alpha << "; " << endl;
    reportFileM << "gamma_switching_parameter_glob = " << switching_gamma << "; " << endl;
    reportFileM << "maximum_inner_iterations = " << maxiterIn << "; " << endl;
    reportFileM << "maximum_outer_iterations = " << maxiterOut << "; " << endl;
    reportFileM << "stopping_criterion = " << stop_rule << "; " << endl;
    reportFileM << "stopping_criterion_tolerance = " << stopCritTol << "; " << endl;
    reportFileM << "number_of_test_samples = " << testN << "; " << endl;
    reportFileM << "single_test = " << singletest << "; " << endl;
    reportFileM << "external_example = " << read_example << "; " << endl;
    
    for (meth = 0; meth < methods; meth++){
      if (singletest) {
	reportFileM << "P{" << meth+1 << "} = " << opts[meth].Ndecomp << "; maximum_inner_iterations{" << meth+1 << "} = " << opts[meth].maxiter_in << "; " << endl;
      }
      reportFileM << "Prediction_Rule{" << meth+1 << "} = " << opts[meth].pred_rule << "; backtracking_strategy{" << meth+1 << "} = " << opts[meth].backtracking << "; inner_stepsize_strategy{" << meth+1 << "} = " << opts[meth].stepsize_update << "; adaptive_stepsize{" << meth+1 << "} = " << opts[meth].stepsize_adaptive << "; " << endl;
      reportFileM << "gamma_switching_parameter{" << meth+1 << "} = " << opts[meth].gamma << "; " << endl;
    }
    for (meth = 0; meth < methodsS; meth++){
      if (singletest) {
	reportFileM << "P{" << meth+1 << "} = " << optsS[meth].Ndecomp << "; maximum_inner_iterations{" << meth+1 << "} = " << optsS[meth].maxiter_in << "; " << endl;
      }
      reportFileM << "Prediction_Rule{" << meth+1 << "} = " << optsS[meth].pred_rule << "; backtracking_strategy{" << meth+1 << "} = " << optsS[meth].backtracking << "; inner_stepsize_strategy{" << meth+1 << "} = " << optsS[meth].stepsize_update << "; adaptive_stepsize{" << meth+1 << "} = " << optsS[meth].stepsize_adaptive << "; " << endl;
      reportFileM << "gamma_switching_parameter{" << meth+1 << "} = " << optsS[meth].gamma << "; " << endl;
    }
    
    write_var_to_matlab(reportFileM, Ndecomp_def, "P_array");
    write_var_to_matlab(reportFileM, maxiter_in_def, "L_array");
    
    for (meth = 0; meth < methods; meth++) {
      stringstream ss;
      ss << meth+1;
      string methStr = ss.str();
      write_var_to_matlab(reportFileM, result_iter_data[meth], "result_iter_data{"+ methStr + "}");
      write_var_to_matlab(reportFileM, result_time_total[meth], "result_time_total{"+methStr + "}");
      write_var_to_matlab(reportFileM, result_time_commu[meth], "result_time_commu{"+methStr + "}");
      write_var_to_matlab(reportFileM, result_time_backt[meth], "result_time_backt{"+methStr + "}");
      write_var_to_matlab(reportFileM, result_time_predi[meth], "result_time_predi{"+methStr + "}");
      write_var_to_matlab(reportFileM, result_time_stopp[meth], "result_time_stopp{"+methStr + "}");
    }
    for (meth = 0; meth < methods; meth++) {
      stringstream ss;
      ss << methods + meth + 1;
      string methStr = ss.str();
      write_var_to_matlab(reportFileM, result_iter_data_single[meth], "result_iter_data{"+methStr + "}");
      write_var_to_matlab(reportFileM, result_time_total_single[meth], "result_time_total{"+methStr + "}");
      write_var_to_matlab(reportFileM, result_time_commu_single[meth], "result_time_commu{"+methStr + "}");
      write_var_to_matlab(reportFileM, result_time_backt_single[meth], "result_time_backt{"+methStr + "}");
      write_var_to_matlab(reportFileM, result_time_predi_single[meth], "result_time_predi{"+methStr + "}");
      write_var_to_matlab(reportFileM, result_time_stopp_single[meth], "result_time_stopp{"+methStr + "}");
    }
    reportFileM.close();
  }



  //==================================================================//
  //                        Free MPI structures                       //
  //==================================================================//
    /*MPI_Barrier(comm_busy);
    MPI_Group_free(&group_busy);
    MPI_Group_free(&mpi_group_orig);*/

  //==================================================================//
  //                          Clean up memory                         //
  //==================================================================//
  delete [] opts;
  delete [] optsS;
  delete [] result_iter_data;
  delete [] result_time_total;
  delete [] result_time_commu;
  delete [] result_time_backt;
  delete [] result_time_predi;
  delete [] result_time_stopp;
  delete [] result_iter_data_single;
  delete [] result_time_total_single;
  delete [] result_time_commu_single;
  delete [] result_time_backt_single;
  delete [] result_time_predi_single;
  delete [] result_time_stopp_single;

    
  //==================================================================//
  //                        free communicators                        //
  //==================================================================//
    /*MPI_Barrier(MPI_COMM_WORLD);
    // clean up
    MPI_Comm_free(&comm_busy);*/

  //==================================================================//
  //                           Finalize MPI                           //
  //==================================================================//

    MPI_Finalize();
    return 0;
}