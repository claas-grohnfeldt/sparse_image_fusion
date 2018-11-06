/*
 * paths.cpp
 *
 *  Created on: Mar 10, 2014
 *      Author: Claas Grohnfeldt
 */

#include "paths.h"

void getPaths(SpEOPaths *paths, SpEODataIOSetting *dSetting, SpEOParallelSetting *pSetting, int argc, char **argv){

	paths->dataSetID_str = argv[3];
	for(int i=0; i<dSetting->dir_tmp_patches_additional_num; i++){
		if(i==0){
			paths->dir_tmp_patches_additional = new string[dSetting->dir_tmp_patches_additional_num];
		}
		paths->dir_tmp_patches_additional[i] = argv[argc-dSetting->dir_tmp_patches_additional_num+i];
	}

	if(dSetting->contUnfinishedRec){
		paths->PathToIncompletePatchSetCSV = argv[35];
	}
	string maindir_path;
	switch(dSetting->platformID){
        case 0:{ // generic
                   maindir_path = "data";
		           paths->dir_out = "results";
                   paths->dir_tmp = "tmp";
                   break;
        }case 1:{ //SuperMUC - pr45ne - ga39yoz2
                maindir_path = "/gpfs/work/pr45ne/ga39yoz2/data/links_for_fusion";
		paths->dir_out = "/gpfs/work/pr45ne/ga39yoz2/recResults";
                paths->dir_tmp = "/gpfs/scratch/pr45ne/ga39yoz2/JSparseFI/tmp";
		break;
	    }case 2:{ // CG local
		    maindir_path = "./data";
		    paths->dir_out = "results";
		    paths->dir_tmp = "tmp";
                break;
           }default:{
		cout << "UNKNOWN platformID!" << endl << endl;
		cerr << "UNKNOWN platformID!" << endl << endl;
	  }
	}
	// define dataset paths and set ID according to the following encryption:
		/*
		 * 1st digit: ID of Machine the date is stored on:
		 *                 1 = SuperMUC - pr45ne - ga39yoz2
		 *                 2 = CG local
		 *                 3 = ...
		 *                 ...
		 * 2nd digit: ID of sensor of original data from which the dataset was synthesized:
		 *                 0 = none: dataset is original (was not synthesized - no high resolution reference image is available)
		 *                 1 = HySpex
		 *                 2 = HyMap
		 *                 3 = Aviris
		 *                 4 = WorldView2
		 *                 5 = Quickbird
		 *                 6 = IKONOS
		 *                 7 = UltraCam
		 *                 8 = HYDICE
		 *                 9 = ...
		 *                 ...
		 * 3rd digit: ID of low resolution sensor (corresponding to ImY):
		 *                 1 = HySpex
		 *                 2 = HyMap
		 *                 3 = Aviris
		 *                 4 = WorldView2
		 *                 5 = Quickbird
		 *                 6 = IKONOS
		 *                 7 = UltraCam
		 *                 8 = HYDICE
		 *                 9 = ...
		 *
		 * 4th digit: ID of high resolution sensor (corresponding to ImX):
		 *                 1 = HySpex
		 *                 2 = HyMap
		 *                 3 = Aviris
		 *                 4 = WorldView2
		 *                 5 = Quickbird
		 *                 6 = IKONOS
		 *                 7 = UltraCam
		 *                 8 = HYDICE
		 *                 9 = ...
		 *
		 * 5th digit: kind of input image pair:
		 *                 1 = multispectral-panchromatic
		 *                 2 = hyperspectral-multispectral
		 *                 3 = hyperspectral-panchromatic
		 *
		 * 6th digit: filter kernel:
	     *                 1 = gauss
	     *                 2 = box-car
	     *
		 * 7th digit: Sight ID relative to sensor if multiple scenes are available:
		 *  e.g. for HySpex  0 = Munich Allianz arena (2012)
		 *                   1 = Munich Olympic Park (2012)
		 *                   2 = Munich Olympic Park - polished (2012)
		 *                   3 = Kaufbeuern (?)
		 *                   4 = Oberpfaffenhofen (2015)
		 *                 ...
		 * 8th digit: size ID (relative to digits above):
		 *   e.g. :        0 = 1500 x 1500
		 *                 2 =  800 x  800
		 *                 3 =  400 x  400
		 *                 4 =  360 x  520
		 *                 5 = ...
		 *
		 * 9th & 10th digit: down-sampling factor / resolution ratio ID:
		 *                 01 = 1
		 *                 02 = 2
		 *                 04 = 2
		 *                 06 = 6
		 *                 10 = 10
		 *                 etc.
		 *
		 * 11th & 12th digit: SNR ID:
		 *                 00 <= (SNR=inf)
		 *                 10 <= (SNR=10)
		 *                 20 <= (SNR=20)
		 *                 30 <= (SNR=30)
		 *                 40 <= (SNR=40)
		 *
		 * 13th digit: redundant digit for non-regular datasets
		 *
		 */

	int my_rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);


    ////*****************************************************************************************************////
    ////                                                                                                     ////
    ////                                      Platform independent                                           ////
    ////                                                                                                     ////
    ////*****************************************************************************************************////

    // Hyperspectral-Multispectral (J-SparseFI-HM)
    if(      paths->dataSetID_str == "11119211105350"){ paths->dir_in = maindir_path + "/" + "HS_MS"  + "/" + "11119211105350_2013IEEEGRSSDFC_Sentinel2_Univ"   + "/" + "InputData" + "/" + "links";
    }else if(paths->dataSetID_str == "109211103990"){   paths->dir_in = maindir_path + "/" + "HS_MS"  + "/" + "109211103990_EnMAP_Sentinel2"                    + "/" + "InputData" + "/" + "links";
    }else if(paths->dataSetID_str == "3315211304350"){  paths->dir_in = maindir_path + "/" + "HS_MS"  + "/" + "3315211304350_Aviris_IndianPines_WV3_VNIR_SWIR"  + "/" + "InputData" + "/" + "links";
    }else if(paths->dataSetID_str == "3315212405350"){  paths->dir_in = maindir_path + "/" + "HS_MS"  + "/" + "3315212405350_Aviris_Cuprite_sc03_WV3_VNIR_SWIR" + "/" + "InputData" + "/" + "links";
    }else if(paths->dataSetID_str == "335213304350"){   paths->dir_in = maindir_path + "/" + "HS_MS"  + "/" + "335213304350_Aviris_Moffett_Field"               + "/" + "InputData" + "/" + "links";
    }else if(paths->dataSetID_str == "665211108350"){   paths->dir_in = maindir_path + "/" + "HS_MS"  + "/" + "665211108350_ROSIS_Pavia_Univeristy"             + "/" + "InputData" + "/" + "links";
    }else if(paths->dataSetID_str == "774212106350"){   paths->dir_in = maindir_path + "/" + "HS_MS"  + "/" + "774212106350_Headwall_Chikusei_nonUrban"         + "/" + "InputData" + "/" + "links";
    }else if(paths->dataSetID_str == "885211404350"){   paths->dir_in = maindir_path + "/" + "HS_MS"  + "/" + "885211404350_HYDICE_WashDC_Mall"                 + "/" + "InputData" + "/" + "links";
    }
    // Multispectral-Panchromatic (Pan-sharpening: SparseFI & J-SparseFI)
    else if(paths->dataSetID_str == "155111203350"){         paths->dir_in = maindir_path + "/" + "MS_PAN"  + "/" + "155111203350_HySpex_Olymp_3600x1200"       + "/" + "InputData" + "/" + "links";
    }else if(paths->dataSetID_str == "2444101104000"){  paths->dir_in = maindir_path + "/" + "MS_PAN" + "/" + "2444101104000_WV2_REAL_scene"                    + "/" + "InputData" + "/" + "links";
    }else{ 
    	if(my_rank==0){
    		cout << endl << "ERROR: Undefined dataset: " << paths->dataSetID_str << "! " << endl << endl;
    	}
    	MPI_Barrier(MPI_COMM_WORLD);
    	cerr << endl << "ERROR: Undefined dataset: " << paths->dataSetID_str << "! " << endl << endl;
    	exit(2);
    }
    
    paths->fname_ImZ_ref                   = paths->dir_in + "/" + "slink_to_ImZ_ref.dat";
    paths->fname_ImZ_init_rec              = paths->dir_in + "/" + "slink_to_ImZ_init_rec.dat";
    paths->fname_ImZ_init_ImY_US           = paths->dir_in + "/" + "slink_to_ImZ_init_ImY_US.dat";
    paths->fname_ImY                       = paths->dir_in + "/" + "slink_to_ImY.dat";
    paths->fname_ImX                       = paths->dir_in + "/" + "slink_to_ImX.dat";
    paths->fname_SRF                       = paths->dir_in + "/" + "slink_to_SRF.csv";
    paths->fname_SRF_for_Spectral_Grouping = paths->dir_in + "/" + "slink_to_SRF_for_Spectral_Grouping.csv";
    paths->fname_SubspaceTransformMat      = paths->dir_in + "/" + "slink_to_HySure_output" + "/" + "EndmemberMat.csv";
    paths->fname_DictLR                    = paths->dir_in + "/" + "slink_to_DictLR.csv";
    paths->fname_DictHR                    = paths->dir_in + "/" + "slink_to_DictHR.csv";
    paths->fname_ImZ                       = paths->dataSetID_str + "_rec";





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
             paths->dir_tmp_patches = paths->dir_tmp + "/patches/" + mystringstr + "_" + dSetting->jobID;
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
}
