#include <mpi.h>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <iomanip>
#include <time.h>
#include <string>
#include <math.h>
#include "mpi_counter.h"

#ifdef _OPENMP
#include <omp.h>
#endif


using namespace std;

int main(int argc, char* argv[]) {

	/*===================================================================*
	 *                        MPI initialization                         *
	 *===================================================================*/
#ifdef _OPENMP
      int iprovided;
      MPI_Init_thread(&argc,&argv,MPI_THREAD_FUNNELED,&iprovided);
#else
      MPI_Init(&argc,&argv);
#endif
      int time0 = MPI_Wtime();
      
      int p_groups = 2;
      int subprocesses = 2;
      int pLast_sub = 100;
      int workStealingTurns = 2;
      int subtasks = 20;
      int my_global_rank; int my_global_processes;
      MPI_Comm_rank(MPI_COMM_WORLD, &my_global_rank);
      MPI_Comm_size(MPI_COMM_WORLD, &my_global_processes);
      
      MPI_Barrier(MPI_COMM_WORLD);
      /*===================================================================*
	*              Create MPI groups and sub-communicators              *
	*===================================================================*/
	MPI_Group mpi_group_orig;
	MPI_Comm_group(MPI_COMM_WORLD, &mpi_group_orig);

	//################
	//## comm_busy  ##
	//################
	MPI_Group group_busy;
#ifdef _OPENMP
	int numProcBusy = p_groups;
#else
	int numProcBusy = p_groups*subprocesses;
#endif
	int ranges_busy[1][3] = {{0,numProcBusy-1,1}};
	int rangesIdl_busy[1][3] = {{numProcBusy,my_global_processes-1,1}};
	if(my_global_rank < numProcBusy){
		MPI_Group_range_incl(mpi_group_orig, 1, ranges_busy, &group_busy);
	}else{
		MPI_Group_range_incl(mpi_group_orig, 1, rangesIdl_busy, &group_busy);
	}
	MPI_Comm comm_busy;
	MPI_Comm_create(MPI_COMM_WORLD, group_busy, &comm_busy);
	  
	  

	/*===================================================================*
	 *                      DO something                    *
	 *===================================================================*/
	if(my_global_rank < numProcBusy){

		int my_rank;
		int my_processes;
		MPI_Comm_rank(comm_busy, &my_rank);
		MPI_Comm_size(comm_busy, &my_processes);
		
#ifndef _OPENMP
		// Definition of groups of tasks which work on the same process (similarly to threads in OpenMP)d
		MPI_Barrier(comm_busy);
		MPI_Comm comm_patch;
		MPI_Group group_patch;
		MPI_Comm_group(comm_busy, &group_patch);
		int ranges[1][3] = {{(my_rank/subprocesses)*subprocesses,(my_rank/subprocesses+1)*subprocesses-1,1}};
		MPI_Group_range_incl(group_busy, 1, ranges, &group_patch);
		MPI_Comm_create(comm_busy, group_patch, &comm_patch);
		
		int my_patch_rank;
		int my_patch_processes;
		MPI_Comm_rank(comm_patch, &my_patch_rank);
		MPI_Comm_size(comm_patch, &my_patch_processes);
#endif

		int jP = 0;
#ifndef _OPENMP
		jP = my_rank / subprocesses;
#else
		jP = my_rank;
#endif

		// additional work stealing
		int lastFixedPatch = pLast_sub;
		if (workStealingTurns >= 0){
			lastFixedPatch -= (pLast_sub+1)%p_groups + workStealingTurns*p_groups;
			lastFixedPatch = max(lastFixedPatch, p_groups-1);
		}

		struct mpi_counter_t *c;
		c = create_counter(0, lastFixedPatch+1, comm_busy);

		while(jP<=pLast_sub){//      as long as there are patches to work on
#ifndef _OPENMP
			if(my_rank%subprocesses==0){
				cout << "[" << my_rank << "] iP_local=" << jP  << endl;
			}
#else
			cout << "[" << my_rank << "] iP_local=" << jP << endl;
#endif

			
			
#ifndef _OPENMP
			for (int ipp=0; ipp < subtasks / subprocesses; ipp++){
#else
#pragma omp parallel for schedule(dynamic)
			for (int ipo=0; ipo < subtasks; ipo++){
#endif
				int c = 1+1;
			}
			

			if (jP <= lastFixedPatch - p_groups || workStealingTurns < 0) {
				jP += p_groups;
			}
			else {
#ifndef _OPENMP
				if (my_rank%subprocesses == 0) {
					jP = increment_counter(c, 1) - 1; // the root process of each group gets a patch number to work on, by the current counter
				}
				int jP_tmp[1];
				if (my_rank%subprocesses == 0) {
					jP_tmp[0] = jP;
				}
				MPI_Bcast(&jP_tmp, 1, MPI_INT, 0, comm_patch); // broadcast the patch to work on to the group members
				jP = jP_tmp[0];
#else
				jP = increment_counter(c, 1) - 1;
#endif
			}
		}

		MPI_Barrier(comm_busy);delete_counter(&c); // the barrier is necessary here (otherwise the counter is deleted, while others might access it)

#ifndef _OPENMP
		MPI_Comm_free(&comm_patch);
		MPI_Group_free(&group_patch);
#endif

	}

      /*===================================================================*
	*                        Free MPI structures                        *
	*===================================================================*/
	MPI_Barrier(comm_busy);
	MPI_Group_free(&group_busy);
	MPI_Group_free(&mpi_group_orig);
	MPI_Comm_free(&comm_busy);
	
	cout << MPI_Wtime()-time0 << endl;

	MPI_Finalize();
	return 0;
}
