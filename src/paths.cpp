/*
 * paths.cpp
 *
 *  Created on: Mar 10, 2014
 *      Author: Claas Grohnfeldt
 */

#include "paths.h"

void getPaths(SpEOPaths *paths, SpEODataIOSetting *dSetting, SpEOParallelSetting *pSetting, int argc, char **argv){

	int my_rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    // outsourced as program arguments
    paths->fname_ImX           = argv[36]; // paths->dir_in + "/" + "slink_to_ImX.dat";
    paths->fname_ImY           = argv[37]; // paths->dir_in + "/" + "slink_to_ImY.dat";
    paths->fname_ImZ_init      = argv[38]; // paths->dir_in + "/" + "slink_to_ImZ_init_rec.dat";
    paths->fname_ImZ_ref       = argv[39]; // paths->dir_in + "/" + "slink_to_ImZ_ref.dat";
    paths->fname_SRF           = argv[40]; // paths->dir_in + "/" + "slink_to_SRF.csv";
    
    paths->dir_out             = argv[41];
    paths->dir_tmp             = argv[42];
    paths->fname_ImZ_out       = "JSparesFIHM_fusion_result";
    
     // get current time and broadcast to all processes

     char buf[15]="";
     if(my_rank==0){
             time_t t = time(0);
             struct tm  tstruct;
             tstruct = *localtime(&t);
             strftime(buf, sizeof(buf), "%y%m%d_%H%M%S", &tstruct);
     }
     MPI_Bcast(&buf, sizeof(buf), MPI_CHAR, 0,MPI_COMM_WORLD);
     string mystringstr(buf);
     if(pSetting->store_patches_tmp_on_drive){
             paths->dir_tmp_patches = paths->dir_tmp + "/patches/" + mystringstr;// + "_" + dSetting->jobID;
             if(my_rank==0){
                     string tmpp(paths->dir_tmp);
                     tmpp += "/patches";
                     mkdir(paths->dir_tmp.c_str(), 0777);
                     chmod(paths->dir_tmp.c_str(), S_IRWXU | S_IRWXG | S_IRWXO);
                     mkdir(tmpp.c_str(), 0777);
                     chmod(tmpp.c_str(), S_IRWXU | S_IRWXG | S_IRWXO);
                     mkdir(paths->dir_tmp_patches.c_str(), 0777);
                     chmod(paths->dir_tmp_patches.c_str(), S_IRWXU | S_IRWXG | S_IRWXO);
             }
     }

    // obsolete:
	// for(int i=0; i<dSetting->dir_tmp_patches_additional_num; i++){
	// 	if(i==0){
	// 		paths->dir_tmp_patches_additional = new string[dSetting->dir_tmp_patches_additional_num];
	// 	}
	// 	paths->dir_tmp_patches_additional[i] = argv[argc-dSetting->dir_tmp_patches_additional_num+i];
	// }
	// if(dSetting->contUnfinishedRec){
	// 	paths->PathToIncompletePatchSetCSV = "<path_to_csv_file>"; // formally set as program argument argv[37];
	// }
}
