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
	    case 1:{ //SuperMUC
                maindir_path = "/gpfs/work/pr45ne/ga39yoz2/data/links_for_fusion";
		paths->dir_out = "/gpfs/work/pr45ne/ga39yoz2/recResults";
                paths->dir_tmp = "/gpfs/scratch/pr45ne/ga39yoz2/JSparseFI/tmp";
		break;
	   }case 2:{ // CG local
		maindir_path = "/mnt/ssd2/data/links_for_fusion"; //"/mnt/NTFS/data/links_for_fusion";
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
		 *                 1 = SuperMUC
		 *                 2 = CG-PC
		 *                 3 = SP-PC
		 *                 4 = ...
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

//	cout << endl << endl << "[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[ paths->dataSetID_str = " << paths->dataSetID_str << "]]]]]]]]]]]]]]]] " << endl << endl;
//	else if(paths->dataSetID_str == "2665211108350"){
//		cout << "++++++++++ ++++++++++ ++++++++ hhh1" << endl << endl << endl;
//	}else else if(paths->dataSetID_str == "2665211108990"){
//		cout << "++++++++++ ++++++++++ ++++++++ hhh2" << endl << endl << endl;
//	}else else if(paths->dataSetID_str == "22109211103990"){
//		cout << "++++++++++ ++++++++++ ++++++++ hhh3" << endl << endl << endl;
//	}else else if(paths->dataSetID_str == "211119211105350"){
//		cout << "++++++++++ ++++++++++ ++++++++ hhh4" << endl << endl << endl;
//	}
//	switch(paths->dataSetID){







        ////*****************************************************************************************************////
        ////*****************************************************************************************************////
        ////                                                                                                     ////
        ////                                                                                                     ////
        ////                                      Platform independent                                           ////
        ////                                                                                                     ////
        ////                                                                                                     ////
        ////*****************************************************************************************************////
        ////*****************************************************************************************************////

	//####################################################################################################################################################################################################################
	//#  orig. sensor    - LR sensor       - HR sensor       - fusion type - filter kernel - scene            - size ID       - fDS         - SNR          - redundant digit (>0 for non-regular datasets #
	//#  3 (Aviris)      - 3 (Aviris)      - 5 (Quickbird)   - 2 (HS-MS)   - 1 (gauss)     - 3 (Moffett Field)- 3 (360x360)   - 04 (fDS=4)  - 35 (SNR=35db)- 0                                            #
	//####################################################################################################################################################################################################################
        if(paths->dataSetID_str == "11119211105350"){
		paths->dir_in        = maindir_path + "/" + "11119211105350_2013IEEEGRSSDFC_Sentinel2_Univ" + "/" + "InputData" + "/" + "links";
	}else if(paths->dataSetID_str == "11119212105350"){
                  paths->dir_in        = maindir_path + "/" + "11119212105350_2013IEEEGRSSDFC_Sentinel2_cloud" + "/" + "InputData" + "/" + "links";
        }else if(paths->dataSetID_str == "109211103990"){
                paths->dir_in        = maindir_path + "/" + "109211103990_EnMAP_Sentinel2" + "/" + "InputData" + "/" + "links";
        }else if(paths->dataSetID_str == "3313211304350"){
                paths->dir_in        = maindir_path + "/" + "3313211304350_Aviris_IndianPines_WV3_VNIR" + "/" + "InputData" + "/" + "links";
        }else if(paths->dataSetID_str == "3313212405350"){
                paths->dir_in        = maindir_path + "/" + "3313212405350_Aviris_Cuprite_sc03_WV3_VNIR" + "/" + "InputData" + "/" + "links";
        }else if(paths->dataSetID_str == "3314211304350"){
                paths->dir_in        = maindir_path + "/" + "3314211304350_Aviris_IndianPines_WV3_SWIR" + "/" + "InputData" + "/" + "links";
        }else if(paths->dataSetID_str == "3314212405350"){
                paths->dir_in        = maindir_path + "/" + "3314212405350_Aviris_Cuprite_sc03_WV3_SWIR" + "/" + "InputData" + "/" + "links";
        }else if(paths->dataSetID_str == "3315211304350"){
                paths->dir_in        = maindir_path + "/" + "3315211304350_Aviris_IndianPines_WV3_VNIR_SWIR" + "/" + "InputData" + "/" + "links";
        }else if(paths->dataSetID_str == "3315212405350"){
                paths->dir_in        = maindir_path + "/" + "3315212405350_Aviris_Cuprite_sc03_WV3_VNIR_SWIR" + "/" + "InputData" + "/" + "links";
        }else if(paths->dataSetID_str == "334211304350"){
                paths->dir_in        = maindir_path + "/" + "334211304350_Aviris_IndianPines" + "/" + "InputData" + "/" + "links";
        }else if(paths->dataSetID_str == "334212404350"){
                paths->dir_in        = maindir_path + "/" + "334212404350_Aviris_Cuprite_sc03_fDS4" + "/" + "InputData" + "/" + "links";
        }else if(paths->dataSetID_str == "334212405350"){
                paths->dir_in        = maindir_path + "/" + "334212405350_Aviris_Cuprite_sc03_fDS5" + "/" + "InputData" + "/" + "links";
        }else if(paths->dataSetID_str == "335213304350"){
                paths->dir_in        = maindir_path + "/" + "335213304350_Aviris_Moffett_Field" + "/" + "InputData" + "/" + "links_SNR35";
        }else if(paths->dataSetID_str == "335213304990"){
                paths->dir_in        = maindir_path + "/" + "335213304350_Aviris_Moffett_Field" + "/" + "InputData" + "/" + "links_SNRinf";
        }else if(paths->dataSetID_str == "665211108350"){
                paths->dir_in        = maindir_path + "/" + "665211108350_ROSIS_Pavia_Univeristy" + "/" + "InputData" + "/" + "links";
        }else if(paths->dataSetID_str == "665211108990"){
                paths->dir_in        = maindir_path + "/" + "665211108990_ROSIS_Pavia_University_SNRinf" + "/" + "InputData" + "/" + "links";
        }else if(paths->dataSetID_str == "774211106350"){
                paths->dir_in        = maindir_path + "/" + "774211106350_Headwall_Chikusei_urban" + "/" + "InputData" + "/" + "links";
        }else if(paths->dataSetID_str == "774212106350"){
                paths->dir_in        = maindir_path + "/" + "774212106350_Headwall_Chikusei_nonUrban" + "/" + "InputData" + "/" + "links";
        }else if(paths->dataSetID_str == "885211404350"){
                paths->dir_in        = maindir_path + "/" + "885211404350_HYDICE_WashDC_Mall" + "/" + "InputData" + "/" + "links";
        }else if(paths->dataSetID_str == "155111203350"){
                paths->dir_in        = maindir_path + "/" + "155111203350_HySpex_Olymp_3600x1200" + "/" + "InputData" + "/" + "links";
        }else if(paths->dataSetID_str == "2444101104000"){
                paths->dir_in        = maindir_path + "/" + "2444101104000_WV2_REAL_scene" + "/" + "InputData" + "/" + "links";
        }else if(paths->dataSetID_str == "2444102104000"){
                paths->dir_in        = maindir_path + "/" + "2444102104000_WV2_REAL_HongKong_from_Naoto" + "/" + "InputData" + "/" + "links";
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
	//paths->fname_SRF_estimated             = paths->dir_in + "/" + "slink_to_SRF_estimated.csv";
	//paths->fname_ImX_shifted               = paths->dir_in + "/" + "slink_to_ImX_shifted.dat";
	paths->fname_SRF_for_Spectral_Grouping = paths->dir_in + "/" + "slink_to_SRF_for_Spectral_Grouping.csv";
	paths->fname_SubspaceTransformMat      = paths->dir_in + "/" + "slink_to_HySure_output" + "/" + "EndmemberMat.csv";
	paths->fname_DictLR                    = paths->dir_in + "/" + "slink_to_DictLR.csv";
	paths->fname_DictHR                    = paths->dir_in + "/" + "slink_to_DictHR.csv";
	
	paths->fname_ImZ                       = paths->dataSetID_str + "_rec";




//########## DO NOT MODIFY BEYOND THIS LINE ##########
  
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




        //####################################################################################################################################################################################################################
	//        //#  orig. sensor    - LR sensor       - HR sensor       - fusion type - filter kernel - scene            - size ID       - fDS         - SNR          - redundant digit (>0 for non-regular datasets #
	//                //#  3 (Aviris)      - 3 (Aviris)      - 5 (Quickbird)   - 2 (HS-MS)   - 1 (gauss)     - 3 (Moffett Field)- 3 (360x360)   - 04 (fDS=4)  - 35 (SNR=35db)- 0                                            #
	//                        //####################################################################################################################################################################################################################i
	//
	//                        //#################################################################################################################################################################################################################
	//                        //      //#  platform     - orig. sensor    - LR sensor       - HR sensor       - fusion type - filter kernel - scene         - size ID       - fDS         - SNR           - redundant digit (>0 for non-regular datasets #
	//                        //      //#  1 (SuperMUC) - 1 (HySpex)      - 5 (QuickBird)   - 5 (QuickBird)   - 1 (MS-Pan)  - 1 (gauss)     - 1 (MUC_Oly)   - 2 (3600x1200) - 03 (fDS=3) - 35 (SNR=35db) - 0                                            #
	//                        //      //#################################################################################################################################################################################################################
	//
	//                        1551112033t
	//                        5
	//
	//                        334211304350
	//
	//
	//	*                 1 = HySpex
	//	*                 2 = HyMap
	//	*                 3 = Aviris
	//	*                 4 = WorldView2
	//	*                 5 = Quickbird
	//	*                 6 = IKONOS
	//	*                 7 = UltraCam
	//	*                 8 = HYDICE
	




//        ////*****************************************************************************************************////
//        ////*****************************************************************************************************////
//        ////                                                                                                     ////
//        ////                                                                                                     ////
//        ////                                      Platform: SuperMUC - CG                                        ////
//        ////                                                                                                     ////
//        ////                                                                                                     ////
//        ////*****************************************************************************************************////
//        ////*****************************************************************************************************////
//
//	/*******************************************
//	 *  small test data sets                   *
//	 *******************************************/
//
//	if(paths->dataSetID_str == "1114211004350"){
//	//#################################################################################################################################################################################################################
//	//#  platform     - orig. sensor    - LR sensor       - HR sensor       - fusion type - filter kernel - scene         - size ID       - fDS         - SNR          - redundant digit (>0 for non-regular datasets #
//	//#  1 (SuperMUC) - 1 (HySpex)      - 1 (HySpex)      - 4 (WorldView-2) - 2 (HS-MS)   - 1 (gauss)     - 1 (MUC_Oly)   - 0 (60x80)     - 04 (fDS=4)  - 35 (SNR=35)  - 0                                            #
//	//#################################################################################################################################################################################################################
//		paths->dir_in         = "/gpfs/work/pr45ne/ga39yoz2/data/HySpex_MUC/Olymp/posished_2012_2m/gauss/60x80";
//		paths->fname_ImZ_ref         = paths->dir_in + "/" + "HSHR"                      + "/" + "HySpex_OLY_60x80_HSHR.dat";
//		paths->fname_ImZ_init_rec    = paths->dir_in + "/" + "../FusionResults"             + "/" + "090459_ID2114211004350_psz3_itr1_NwCf1_ImXSm1_lX10_lY10_lZ1_lZ0In1stIter1_fllImOpt1_lXIm10_lYIm10" + "/" + "2114211004350_rec";
//		paths->fname_ImZ_init_ImY_US = paths->dir_in + "/" + "HSLR_fDS4_US_via_bilinear" + "/" + "HySpex_VNIR_60x80_HSLR_fDS4_US_via_bilinear_SNR35.dat";
//		paths->fname_ImY      = paths->dir_in + "/" + "HSLR_fDS4"     + "/" + "HySpex_VNIR_60x80_HSLR_fDS4_SNR35.dat";
//		paths->fname_ImX      = paths->dir_in + "/" + "WV2_MSHR"      + "/" + "HySpex_VNIR_60x80_WV2_MSHR_SNR35.dat";
//		paths->dir_out        = "/gpfs/work/pr45ne/ga39yoz2/recResults";
//		paths->fname_ImZ      = "1114211004350_rec";
//		paths->fname_SRF                          = "/gpfs/work/pr45ne/ga39yoz2/data/SRF/high_precision/SRF_WV2_MS_gridded_to_HySpex_VNIR_centers.csv";
//		paths->fname_SRF_for_Spectral_Grouping    = "/gpfs/work/pr45ne/ga39yoz2/data/SRF/high_precision/SRF_WV2_MS_gridded_to_HySpex_VNIR_centers.csv";
//		paths->dir_tmp        = "/gpfs/scratch/pr45ne/ga39yoz2/JSparseFI/tmp";
//		//break;
//	}
//	else if(paths->dataSetID_str == "1114211108000"){
//	//#################################################################################################################################################################################################################
//	//#  platform     - orig. sensor    - LR sensor       - HR sensor       - fusion type - filter kernel - scene         - size ID       - fDS         - SNR          - redundant digit (>0 for non-regular datasets #
//	//#  1 (SuperMUC) - 1 (HySpex)      - 1 (HySpex)      - 4 (WorldView-2) - 2 (HS-MS)   - 1 (gauss)     - 1 (MUC_Oly)   - 1 (240x120)   - 08 (fDS=8)  - 00 (SNR=inf) - 0                                            #
//	//#################################################################################################################################################################################################################
//		paths->dir_in         = "/gpfs/work/pr45ne/ga39yoz2/data/HySpex_MUC/Olymp/posished_2012_2m/gauss/240x120";
//		paths->fname_ImZ_ref  = paths->dir_in + "/" + "HSHR"          + "/" + "HySpex_2m_240x120_HSHR.dat";
//		paths->fname_ImZ_init_rec    = paths->dir_in + "/" + "../FusionResults"             + "/" + "2114211108000_rec_BayesianSparse";
//		paths->fname_ImZ_init_ImY_US = paths->dir_in + "/" + "HSLR_fDS8_US_via_bilinear" + "/" + "HySpex_240x120_HSLR_fDS8_US_via_bilinear.dat";
//		paths->fname_ImY      = paths->dir_in + "/" + "HSLR_fDS8"     + "/" + "HySpex_240x120_HSLR_fDS8.dat";
//		paths->fname_ImX      = paths->dir_in + "/" + "WV2_MSHR"      + "/" + "HySpex_240x120_WV2_MSHR.dat";
//		paths->dir_out        = "/gpfs/work/pr45ne/ga39yoz2/recResults";
//		paths->fname_ImZ      = "1114211108000_rec";
//		paths->fname_SRF                          = "/gpfs/work/pr45ne/ga39yoz2/data/SRF/high_precision/SRF_WV2_MS_gridded_to_HySpex_VNIR_centers.csv";
//		paths->fname_SRF_for_Spectral_Grouping    = "/gpfs/work/pr45ne/ga39yoz2/data/SRF/high_precision/SRF_WV2_MS_gridded_to_HySpex_VNIR_centers.csv";
//		paths->dir_tmp        = "/gpfs/scratch/pr45ne/ga39yoz2/JSparseFI/tmp";
//		//break;
//	}
//	else if(paths->dataSetID_str == "1115211105000"){
//	//#################################################################################################################################################################################################################
//	//#  platform     - orig. sensor    - LR sensor       - HR sensor       - fusion type - filter kernel - scene         - size ID       - fDS         - SNR          - redundant digit (>0 for non-regular datasets #
//	//#  1 (SuperMUC) - 1 (HySpex)      - 1 (HySpex)      - 5 (Quickbird)   - 2 (HS-MS)   - 1 (gauss)     - 1 (MUC_Oly)   - 1 (240x120)   - 05 (fDS=5)  - 00 (SNR=inf) - 0                                            #
//	//#################################################################################################################################################################################################################
//		paths->dir_in         = "/gpfs/work/pr45ne/ga39yoz2/data/HySpex_MUC/Olymp/posished_2012_2m/gauss/240x120";
//		paths->fname_ImZ_ref = paths->dir_in + "/" + "HSHR"          + "/" + "HySpex_2m_240x120_HSHR.dat";
//		paths->fname_ImZ_init_rec    = paths->dir_in + "/" + "../FusionResults"             + "/" + "2115211105000_rec_BayesianSparse";
//		paths->fname_ImZ_init_ImY_US = paths->dir_in + "/" + "HSLR_fDS5_US_via_bilinear" + "/" + "HySpex_240x120_HSLR_fDS5_US_via_bilinear.dat";
//		paths->fname_ImY      = paths->dir_in + "/" + "HSLR_fDS5"     + "/" + "HySpex_240x120_HSLR_fDS5.dat";
//		paths->fname_ImX      = paths->dir_in + "/" + "QB_MSHR"      + "/" + "HySpex_240x120_QB_MSHR.dat";
//		paths->dir_out        = "/gpfs/work/pr45ne/ga39yoz2/recResults";
//		paths->fname_ImZ      = "1115211105000_rec";
//		paths->fname_SRF                          = "/gpfs/work/pr45ne/ga39yoz2/data/SRF/high_precision/SRF_QB_MS_gridded_to_HySpex_VNIR_centers.csv";
//		paths->fname_SRF_for_Spectral_Grouping    = "/gpfs/work/pr45ne/ga39yoz2/data/SRF/high_precision/SRF_QB_MS_gridded_to_HySpex_VNIR_centers_modForSpectralGrouping.csv";
//		paths->dir_tmp        = "/gpfs/scratch/pr45ne/ga39yoz2/JSparseFI/tmp";
//		//break;
//	}
//	else if(paths->dataSetID_str == "1144111104000"){
//	//#################################################################################################################################################################################################################
//	//#  platform     - orig. sensor    - LR sensor       - HR sensor       - fusion type - filter kernel - scene         - size ID       - fDS         - SNR          - redundant digit (>0 for non-regular datasets #
//	//#  1 (SuperMUC) - 1 (HySpex)      - 4 (WorldView-2) - 4 (WorldView-2) - 1 (MS-Pan)  - 1 (gauss)     - 1 (MUC_Oly)   - 1 (240x120)   - 04 (fDS=4)  - 00 (SNR=inf) - 0                                            #
//	//#################################################################################################################################################################################################################
//		paths->dir_in         = "/gpfs/work/pr45ne/ga39yoz2/data/HySpex_MUC/Olymp/posished_2012_2m/gauss/240x120";
//		paths->fname_ImZ_ref = paths->dir_in + "/" + "WV2_MSHR"       + "/" + "HySpex_240x120_WV2_MSHR.dat";
//		paths->fname_ImZ_init_rec    = paths->dir_in + "/" + "WV2_MSHR"       + "/" + "HySpex_240x120_WV2_MSHR.dat";
//		paths->fname_ImZ_init_ImY_US = paths->dir_in + "/" + "WV2_MSHR"       + "/" + "HySpex_240x120_WV2_MSHR.dat";
//		paths->fname_ImZ_init = paths->dir_in + "/" + "WV2_MSHR"       + "/" + "HySpex_240x120_WV2_MSHR.dat";
//		paths->fname_ImY      = paths->dir_in + "/" + "WV2_MSLR_fDS4"  + "/" + "HySpex_240x120_WV2_MSLR_fDS4.dat";
//		paths->fname_ImX      = paths->dir_in + "/" + "WV2_PanHR"      + "/" + "HySpex_240x120_WV2_PanHR.dat";
//		paths->dir_out        = "/gpfs/work/pr45ne/ga39yoz2/recResults";
//		paths->fname_ImZ      = "1144111104000_rec";
//		paths->fname_SRF                          = "/gpfs/work/pr45ne/ga39yoz2/data/SRF/high_precision/old/SRF_WV2_Pan_gridded_to_WV2_MS_centers_new.csv";
//		paths->fname_SRF_for_Spectral_Grouping    = "/gpfs/work/pr45ne/ga39yoz2/data/SRF/high_precision/old/SRF_WV2_Pan_gridded_to_WV2_MS_centers_new.csv";
//		paths->dir_tmp        = "/gpfs/scratch/pr45ne/ga39yoz2/JSparseFI/tmp";
//		//break;
//	}
//
//	/*******************************************
//	 *  standard datasets                      *
//	 *******************************************/
//	//####################################################################################################################################################################################################################
//	//#  platform     - orig. sensor    - LR sensor       - HR sensor       - fusion type - filter kernel - scene            - size ID       - fDS         - SNR          - redundant digit (>0 for non-regular datasets #
//	//#  1 (SuperMUC) - 1 (HySpex)      - 1 (HySpex)      - 4 (WorldView-2) - 2 (HS-MS)   - 1 (gauss)     - 1 (MUC_Oly)      - 5 (540x540)   - 10 (fDS=10) - 35 (SNR=35db)- 0                                            #
//	//####################################################################################################################################################################################################################
//	else if(paths->dataSetID_str == "1114211510350"){
//		paths->dir_in        = "/gpfs/work/pr45ne/ga39yoz2/data/HySpex/MUC/Olymp/polished/gauss/540x540/InputData";
////		paths->fname_ImZ_ref         = paths->dir_in + "/" + "HSHR_ref_raw"               + "/" + "Muenchen_28_VNIR_FOVx2_raw_rad_atm_polish_geo_540x540_raw.dat";
//		paths->fname_ImZ_ref         = paths->dir_in + "/" + "HSHR_ref_denoised"          + "/" + "HySpex_VNIR_540x540_HSHR_ref_denoised.dat";
////		paths->fname_ImZ_init_rec    = paths->dir_in + "/" + "../FusionResults"              + "/" + "2114211510350_HySpex_fDS10" + "/" + "BayesianSparse" + "/" + "BEST_ID2114211510350_nbSub10"                              + "/" + "2114211510350_rec_.dat";
////		paths->fname_ImZ_init_rec    = paths->dir_in + "/" + "../FusionResults"              + "/" + "2114211510350_HySpex_fDS10" + "/" + "CNMF"           + "/" + "2114211510350_CNMF_mod0PInvBackStep0_M40_I2_2_I1_400_BEST" + "/" + "2114211510350_rec";
////		paths->fname_ImZ_init_rec    = paths->dir_in + "/" + "../FusionResults"              + "/" + "2114211510350_HySpex_fDS10" + "/" + "MAPSMM"         + "/" + "2114211510350_sz540x540_fDS10_snr800_nkt4_npure4_BEST"     + "/" + "2114211510350_HSHR_rec_MAPSMM.dat";
////		paths->fname_ImZ_init_rec    = paths->dir_in + "/" + "../FusionResults"              + "/" + "2114211510350_HySpex_fDS10" + "/" + "JSpFIHM"        + "/" + "after_one_patRec" + "/" + "ID1114211510350_itr1_init1_lXIm10_lYIm10_BEST" + "/" + "1114211510350_rec";
//		paths->fname_ImZ_init_rec    = paths->dir_in + "/" + "HSHR_ref_denoised"          + "/" + "HySpex_VNIR_540x540_HSHR_ref_denoised.dat";
//		paths->fname_ImZ_init_ImY_US = paths->dir_in + "/" + "HSLR_fDS10_US_via_bilinear" + "/" + "HySpex_VNIR_540x540_HSLR_fDS10_US_via_bilinear_SNR35.dat";
//		paths->fname_ImY             = paths->dir_in + "/" + "HSLR_fDS10"                 + "/" + "HySpex_VNIR_540x540_HSLR_fDS10_SNR35.dat";
//		paths->fname_ImX             = paths->dir_in + "/" + "WV2_MSHR"                   + "/" + "HySpex_VNIR_540x540_WV2_MSHR_SNR35.dat";
//		paths->dir_out               = "/gpfs/work/pr45ne/ga39yoz2/recResults";
//		paths->fname_ImZ             = "1114211510350_rec";
//
//		// ===> SRFs (apriori and estimate) and bandwise offset estimation of ImX
//		paths->fname_SRF             = paths->dir_in + "/" + "SRFs"     + "/" + ".csv";
//		// via estimation on the LR level
//		paths->fname_SRF_estimated   = paths->dir_in + "/" + "SRFs"     + "/" + ".csv";
//		paths->fname_ImX_shifted     = paths->dir_in + "/" + "" + "/" + ".dat";
//		// via estimation on the HR level
////		paths->fname_SRF_estimated   = paths->dir_in + "/" + "SRFs"     + "/" + ".csv";
////		paths->fname_ImX_shifted     = paths->dir_in + "/" + "" + "/" + ".dat";
//		// <===
//
////		paths->fname_SRF                       = paths->dir_in + "/" + "SRFs"     + "/" + "SRF_WV2_MS_gridded_to_HySpex_VNIR_centers.csv";
////		paths->fname_SRF_estimated             = paths->dir_in + "/" + "SRFs"     + "/" + "SRF_WV2_MS_gridded_to_HySpex_VNIR_centers_est_LR_A0_B0_G0.csv";
//
//		paths->fname_SRF_for_Spectral_Grouping = paths->dir_in + "/" + "SRFs"     + "/" + "SRF_WV2_MS_gridded_to_HySpex_VNIR_centers.csv";
//		paths->dir_tmp = "/gpfs/scratch/pr45ne/ga39yoz2/JSparseFI/tmp";
//		//break;
//	}
//	//####################################################################################################################################################################################################################
//	//#  platform     - orig. sensor    - LR sensor       - HR sensor       - fusion type - filter kernel - scene            - size ID       - fDS         - SNR          - redundant digit (>0 for non-regular datasets #
//	//#  1 (SuperMUC) - 3 (Aviris)      - 3 (Aviris)      - 4 (WorldView-2) - 2 (HS-MS)   - 1 (gauss)     - 1 (Ind.Pines)    - 3 (360x360)   - 04 (fDS=4)  - 35 (SNR=35db)- 0                                            #
//	//####################################################################################################################################################################################################################
//	else if(paths->dataSetID_str == "1334211304350"){
//		paths->dir_in        = "/gpfs/work/pr45ne/ga39yoz2/data/Aviris/Indian_Pine/220Band_AVIRIS_12June_1992_Indian_Pine_Test_Site_3/supporting/aviris_hyperspectral_data/NS-line_360x360/waterAbsBands_removed/InputData";
////		paths->fname_ImZ_ref         = paths->dir_in + "/" + "HSHR_ref_raw"              + "/" + "19920612_AVIRIS_IndianPine_NS-line_360x360_BBR.dat";
//		paths->fname_ImZ_ref         = paths->dir_in + "/" + "HSHR_ref_denoised"         + "/" + "Aviris_IndianPines_360x360_HSHR_ref_denoised.dat";
////		paths->fname_ImZ_init_rec    = paths->dir_in + "/" + "../FusionResults"             + "/" + "2334211304350_Aviris_IndianPines" + "/" + "BayesianSparse" + "/" + "BEST_ID2334211304350_nbSub35"                              + "/" + "2334211304350_rec_";
////		paths->fname_ImZ_init_rec    = paths->dir_in + "/" + "../FusionResults"             + "/" + "2334211304350_Aviris_IndianPines" + "/" + "CNMF"           + "/" + "2334211304350_CNMF_mod0PInvBackStep0_M20_I2_1_I1_400_BEST" + "/" + "2334211304350_rec";
////		paths->fname_ImZ_init_rec    = paths->dir_in + "/" + "../FusionResults"             + "/" + "2334211304350_Aviris_IndianPines" + "/" + "MAPSMM"         + "/" + "2334211304350_sz360x360_fDS4_snr800_nkt4_npure4_BEST"      + "/" + "2334211304350_HSHR_rec_MAPSMM.dat";
////		paths->fname_ImZ_init_rec    = paths->dir_in + "/" + "../FusionResults"             + "/" + "2334211304350_Aviris_IndianPines" + "/" + "JSpFIHM"        + "/" + "after_one_patRec" + "/" + "ID1334211304350_itr1_init1_lXIm31d6_lYIm1_BEST" + "/" + "1334211304350_rec";
//		paths->fname_ImZ_init_rec    = paths->dir_in + "/" + "HSHR_ref_denoised"         + "/" + "Aviris_IndianPines_360x360_HSHR_ref_denoised.dat";
//		paths->fname_ImZ_init_ImY_US = paths->dir_in + "/" + "HSLR_fDS4_US_via_bilinear" + "/" + "Aviris_IndianPines_360x360_HSLR_fDS4_US_via_bilinear_SNR35.dat";
//		paths->fname_ImY             = paths->dir_in + "/" + "HSLR_fDS4"                 + "/" + "Aviris_IndianPines_360x360_HSLR_fDS4_SNR35.dat";
//		paths->fname_ImX             = paths->dir_in + "/" + "WV2_MSHR"                  + "/" + "Aviris_IndianPines_360x360_WV2_MSHR_SNR35.dat";
//		paths->dir_out               = "/gpfs/work/pr45ne/ga39yoz2/recResults";
//		paths->fname_ImZ             = "1334211304350_rec";
//		// ===> SRFs (apriori and estimate) and bandwise offset estimation of ImX
//		paths->fname_SRF             = paths->dir_in + "/" + "SRFs"     + "/" + "SRF_WV2_MS_gridded_to_Aviris_IndianPines_centers.csv";
//		// via estimation on the LR level
//		paths->fname_SRF_estimated   = paths->dir_in + "/" + "SRFs"     + "/" + "SRF_WV2_MS_gridded_to_Aviris_IndianPines_centers_est_with_offset_LR_A0_B0_G0.csv";
//		paths->fname_ImX_shifted     = paths->dir_in + "/" + "WV2_MSHR_bandwiseShifted_via_SRF_est" + "/" + "Aviris_IndianPines_360x360_WV2_MSHR_SNR35_bandwise_shifted_via_LR_est.dat";
//		// via estimation on the HR level
////		paths->fname_SRF_estimated   = paths->dir_in + "/" + "SRFs"     + "/" + "SRF_WV2_MS_gridded_to_Aviris_IndianPines_centers_est_with_offset_HR_ref_A0_B0_G0.csv";
////		paths->fname_ImX_shifted     = paths->dir_in + "/" + "WV2_MSHR_bandwiseShifted_via_SRF_est" + "/" + "Aviris_IndianPines_360x360_WV2_MSHR_SNR35_bandwise_shifted_via_HR_ref_est.dat";
//		// <===
//		paths->fname_SRF_for_Spectral_Grouping = paths->dir_in + "/" + "SRFs"     + "/" + "SRF_WV2_MS_gridded_to_Aviris_IndianPines_centers_mod.csv";
//		paths->dir_tmp = "/gpfs/scratch/pr45ne/ga39yoz2/JSparseFI/tmp";
//		//break;
//	}
//	//####################################################################################################################################################################################################################
//	//#  platform     - orig. sensor    - LR sensor       - HR sensor       - fusion type - filter kernel - scene            - size ID       - fDS         - SNR          - redundant digit (>0 for non-regular datasets #
//	//#  1 (SuperMUC) - 3 (Aviris)      - 3 (Aviris)      - 4 (WorldView-2) - 2 (HS-MS)   - 1 (gauss)     - 2 (Cuprite Sc03) - 4 (420x360)   - 04 (fDS=4)  - 35 (SNR=35db)- 0                                            #
//	//####################################################################################################################################################################################################################
//	else if(paths->dataSetID_str == "1334212404350"){
//		paths->dir_in        = "/gpfs/work/pr45ne/ga39yoz2/data/Aviris/cuprite_from_aviris_website_f970619t01p02r02c_rfl/sc03_a_rfl_420x360/waterAbsBands_removed/InputData";
////		paths->fname_ImZ_ref         = paths->dir_in + "/" + "HSHR_ref_raw"              + "/" + "f970619t01p02_r02_sc03_a_rfl_420x360.dat";
//		paths->fname_ImZ_ref         = paths->dir_in + "/" + "HSHR_ref_denoised"         + "/" + "Aviris_Cuprite_sc03_420x360_HSHR_ref_denoised.dat";
////		paths->fname_ImZ_init_rec    = paths->dir_in + "/" + "../FusionResults"             + "/" + "2334212404350_Aviris_Cuprite_sc03" + "/" + "BayesianSparse" + "/" + "BEST_ID2334212404350_nbSub52"                              + "/" + "2334212404350_rec_.dat";
////		paths->fname_ImZ_init_rec    = paths->dir_in + "/" + "../FusionResults"             + "/" + "2334212404350_Aviris_Cuprite_sc03" + "/" + "CNMF"           + "/" + "2334212404350_CNMF_mod0PInvBackStep0_M20_I2_1_I1_400_BEST" + "/" + "2334212404350_rec";
////		paths->fname_ImZ_init_rec    = paths->dir_in + "/" + "../FusionResults"             + "/" + "2334212404350_Aviris_Cuprite_sc03" + "/" + "MAPSMM"         + "/" + "2334212404350_sz420x360_fDS4_snr800_nkt6_npure5_BEST"      + "/" + "2334212404350_HSHR_rec_MAPSMM.dat";
////		paths->fname_ImZ_init_rec    = paths->dir_in + "/" + "../FusionResults"             + "/" + "2334212404350_Aviris_Cuprite_sc03" + "/" + "JSpFIHM"        + "/" + "after_one_patRec" + "/" + "ID1334212404350_itr1_init1_lXIm0d1_lYIm1d7_BEST" + "/" + "1334212404350_rec";
//		paths->fname_ImZ_init_rec    = paths->dir_in + "/" + "HSHR_ref_denoised"         + "/" + "Aviris_Cuprite_sc03_420x360_HSHR_ref_denoised.dat";
//		paths->fname_ImZ_init_ImY_US = paths->dir_in + "/" + "HSLR_fDS4_US_via_bilinear" + "/" + "Aviris_Cuprite_sc03_420x360_HSLR_fDS4_US_via_bilinear_SNR35.dat";
//		paths->fname_ImY             = paths->dir_in + "/" + "HSLR_fDS4"                 + "/" + "Aviris_Cuprite_sc03_420x360_HSLR_fDS4_SNR35.dat";
//		paths->fname_ImX             = paths->dir_in + "/" + "WV2_MSHR"                  + "/" + "Aviris_Cuprite_sc03_420x360_WV2_MSHR_SNR35.dat";
//		paths->dir_out               = "/gpfs/work/pr45ne/ga39yoz2/recResults";
//		paths->fname_ImZ             = "1334212404350_rec";
//		// ===> SRFs (apriori and estimate) and bandwise offset estimation of ImX
//		paths->fname_SRF             = paths->dir_in + "/" + "SRFs"     + "/" + "SRF_WV2_MS_gridded_to_Aviris_Cuprite_sc03_centers.csv";
//		// via estimation on the LR level
//		paths->fname_SRF_estimated   = paths->dir_in + "/" + "SRFs"     + "/" + "SRF_WV2_MS_gridded_to_Aviris_Cuprite_sc03_centers_est_with_offset_LR_A0_B0_G0.csv";
//		paths->fname_ImX_shifted     = paths->dir_in + "/" + "WV2_MSHR_bandwiseShifted_via_SRF_est" + "/" + "Aviris_Cuprite_sc03_420x360_WV2_MSHR_SNR35_bandwise_shifted_via_LR_est.dat";
//		// via estimation on the HR level
////		paths->fname_SRF_estimated   = paths->dir_in + "/" + "SRFs"     + "/" + "SRF_WV2_MS_gridded_to_Aviris_Cuprite_sc03_centers_est_with_offset_HR_ref_A0_B0_G0.csv";
////		paths->fname_ImX_shifted     = paths->dir_in + "/" + "WV2_MSHR_bandwiseShifted_via_SRF_est" + "/" + "Aviris_Cuprite_sc03_420x360_WV2_MSHR_SNR35_bandwise_shifted_via_HR_ref_est.dat";
//		// <===
//		paths->fname_SRF_for_Spectral_Grouping = paths->dir_in + "/" + "SRFs"     + "/" + "SRF_WV2_MS_gridded_to_Aviris_Cuprite_sc03_centers_mod.csv";
//		paths->dir_tmp = "/gpfs/scratch/pr45ne/ga39yoz2/JSparseFI/tmp";
//		//break;
//	}
//	//####################################################################################################################################################################################################################
//	//#  platform     - orig. sensor    - LR sensor       - HR sensor       - fusion type - filter kernel - scene            - size ID       - fDS         - SNR          - redundant digit (>0 for non-regular datasets #
//	//#  1 (SuperMUC) - 3 (Aviris)      - 3 (Aviris)      - 4 (WorldView-2) - 2 (HS-MS)   - 1 (gauss)     - 2 (Cuprite Sc03) - 4 (420x360)   - 05(fDS=5)  - 35 (SNR=35db)- 0                                            #
//	//####################################################################################################################################################################################################################
//	else if(paths->dataSetID_str == "1334212405350"){
//		paths->dir_in        = "/gpfs/work/pr45ne/ga39yoz2/data/Aviris/cuprite_from_aviris_website_f970619t01p02r02c_rfl/sc03_a_rfl_420x360/waterAbsBands_removed_fDS5/InputData";
////		paths->fname_ImZ_ref         = paths->dir_in + "/" + "HSHR_ref_raw"              + "/" + "f970619t01p02_r02_sc03_a_rfl_420x360.dat";
//		paths->fname_ImZ_ref         = paths->dir_in + "/" + "HSHR_ref_denoised"         + "/" + "Aviris_Cuprite_sc03_420x360_HSHR_ref_denoised.dat";
////		paths->fname_ImZ_init_rec    = paths->dir_in + "/" + "../FusionResults"             + "/" + "2334212404350_Aviris_Cuprite_sc03" + "/" + "BayesianSparse" + "/" + "BEST_ID2334212404350_nbSub52"                              + "/" + "2334212404350_rec_.dat";
////		paths->fname_ImZ_init_rec    = paths->dir_in + "/" + "../FusionResults"             + "/" + "2334212404350_Aviris_Cuprite_sc03" + "/" + "CNMF"           + "/" + "2334212404350_CNMF_mod0PInvBackStep0_M20_I2_1_I1_400_BEST" + "/" + "2334212404350_rec";
////		paths->fname_ImZ_init_rec    = paths->dir_in + "/" + "../FusionResults"             + "/" + "2334212404350_Aviris_Cuprite_sc03" + "/" + "MAPSMM"         + "/" + "2334212404350_sz420x360_fDS4_snr800_nkt6_npure5_BEST"      + "/" + "2334212404350_HSHR_rec_MAPSMM.dat";
////		paths->fname_ImZ_init_rec    = paths->dir_in + "/" + "../FusionResults"             + "/" + "2334212404350_Aviris_Cuprite_sc03" + "/" + "JSpFIHM"        + "/" + "after_one_patRec" + "/" + "ID1334212404350_itr1_init1_lXIm0d1_lYIm1d7_BEST" + "/" + "1334212404350_rec";
//		paths->fname_ImZ_init_rec    = paths->dir_in + "/" + "HSHR_ref_denoised"         + "/" + "Aviris_Cuprite_sc03_420x360_HSHR_ref_denoised.dat";
//		paths->fname_ImZ_init_ImY_US = paths->dir_in + "/" + "HSLR_fDS5_US_via_bilinear" + "/" + "Aviris_Cuprite_sc03_420x360_HSLR_fDS5_US_via_bilinear_SNR35.dat";
//		paths->fname_ImY             = paths->dir_in + "/" + "HSLR_fDS5"                 + "/" + "Aviris_Cuprite_sc03_420x360_HSLR_fDS5_SNR35.dat";
//		paths->fname_ImX             = paths->dir_in + "/" + "WV2_MSHR"                  + "/" + "Aviris_Cuprite_sc03_420x360_WV2_MSHR_SNR35.dat";
//		paths->dir_out               = "/gpfs/work/pr45ne/ga39yoz2/recResults";
//		paths->fname_ImZ             = "1334212405350_rec";
//		// ===> SRFs (apriori and estimate) and bandwise offset estimation of ImX
//		paths->fname_SRF             = paths->dir_in + "/" + "SRFs"     + "/" + "SRF_WV2_MS_gridded_to_Aviris_Cuprite_sc03_centers.csv";
//		// via estimation on the LR level
//		paths->fname_SRF_estimated   = paths->dir_in + "/" + "SRFs"     + "/" + "SRF_WV2_MS_gridded_to_Aviris_Cuprite_sc03_centers_est_with_offset_LR_A0_B0_G0.csv";
//		paths->fname_ImX_shifted     = paths->dir_in + "/" + "WV2_MSHR_bandwiseShifted_via_SRF_est" + "/" + "Aviris_Cuprite_sc03_420x360_WV2_MSHR_SNR35_bandwise_shifted_via_LR_est.dat";
//		// via estimation on the HR level
////		paths->fname_SRF_estimated   = paths->dir_in + "/" + "SRFs"     + "/" + "SRF_WV2_MS_gridded_to_Aviris_Cuprite_sc03_centers_est_with_offset_HR_ref_A0_B0_G0.csv";
////		paths->fname_ImX_shifted     = paths->dir_in + "/" + "WV2_MSHR_bandwiseShifted_via_SRF_est" + "/" + "Aviris_Cuprite_sc03_420x360_WV2_MSHR_SNR35_bandwise_shifted_via_HR_ref_est.dat";
//		// <===
//		paths->fname_SRF_for_Spectral_Grouping = paths->dir_in + "/" + "SRFs"     + "/" + "SRF_WV2_MS_gridded_to_Aviris_Cuprite_sc03_centers_mod.csv";
//		paths->dir_tmp = "/gpfs/scratch/pr45ne/ga39yoz2/JSparseFI/tmp";
//		//break;
//	}
//	//####################################################################################################################################################################################################################
//	//#  platform     - orig. sensor    - LR sensor       - HR sensor       - fusion type - filter kernel - scene            - size ID       - fDS         - SNR          - redundant digit (>0 for non-regular datasets #
//	//#  1 (SuperMUC) - 3 (Aviris)      - 3 (Aviris)      - 5 (Quickbird)   - 2 (HS-MS)   - 1 (gauss)     - 3 (Moffett Field)- 3 (360x360)   - 04 (fDS=4)  - 35 (SNR=35db)- 0                                            #
//	//####################################################################################################################################################################################################################
//	else if(paths->dataSetID_str == "1335213304350"){
//		paths->dir_in        = "/gpfs/work/pr45ne/ga39yoz2/data/Aviris/Moffet_Field_f080611t01p00r07rdn_c/sc01_ort_img_360x360/waterAbsBands_removed/InputData";
////		paths->fname_ImZ_ref         = paths->dir_in + "/" + "HSHR_ref_raw"              + "/" + "f080611t01p00r07rdn_c_sc01_ort_img_360x360_BBR.dat";
//		paths->fname_ImZ_ref         = paths->dir_in + "/" + "HSHR_ref_denoised"         + "/" + "Aviris_Moffett_Field_360x360_HSHR_ref_denoised.dat";
////		paths->fname_ImZ_init_rec    = paths->dir_in + "/" + "../FusionResults"             + "/" + "2335213304350_Aviris_Moffett_Field" + "/" + "BayesianSparse" + "/" + "BEST_ID2335213304350_nbSub13"                              + "/" + "2335213304350_rec_.dat";
////		paths->fname_ImZ_init_rec    = paths->dir_in + "/" + "../FusionResults"             + "/" + "2335213304350_Aviris_Moffett_Field" + "/" + "CNMF"           + "/" + "2335213304350_CNMF_mod0PInvBackStep0_M30_I2_1_I1_400_BEST" + "/" + "2335213304350_rec";
////		paths->fname_ImZ_init_rec    = paths->dir_in + "/" + "../FusionResults"             + "/" + "2335213304350_Aviris_Moffett_Field" + "/" + "MAPSMM"         + "/" + "2335213304350_sz360x360_fDS4_snr800_nkt7_npure5_BEST"      + "/" + "2335213304350_HSHR_rec_MAPSMM.dat";
////		paths->fname_ImZ_init_rec    = paths->dir_in + "/" + "../FusionResults"             + "/" + "2335213304350_Aviris_Moffett_Field" + "/" + "JSpFIHM"        + "/" + "after_one_patRec" + "/" + "ID1335213304350_itr1_init1_lXIm10_lYIm3d16_BEST" + "/" +  "1335213304350_rec";
//		paths->fname_ImZ_init_rec    = paths->dir_in + "/" + "HSHR_ref_denoised"         + "/" + "Aviris_Moffett_Field_360x360_HSHR_ref_denoised.dat";
//		paths->fname_ImZ_init_ImY_US = paths->dir_in + "/" + "HSLR_fDS4_US_via_bilinear" + "/" + "Aviris_Moffett_Field_360x360_HSLR_fDS4_US_via_bilinear_SNR35.dat";
//		paths->fname_ImY             = paths->dir_in + "/" + "HSLR_fDS4"                 + "/" + "Aviris_Moffett_Field_360x360_HSLR_fDS4_SNR35.dat";
//		paths->fname_ImX             = paths->dir_in + "/" + "QB_MSHR"                   + "/" + "Aviris_Moffett_Field_360x360_QB_MSHR_SNR35.dat";
//		paths->dir_out               = "/gpfs/work/pr45ne/ga39yoz2/recResults";
//		paths->fname_ImZ             = "1335213304350_rec";
//		// ===> SRFs (apriori and estimate) and bandwise offset estimation of ImX
//		paths->fname_SRF             = paths->dir_in + "/" + "SRFs"     + "/" + "SRF_QB_MS_gridded_to_Aviris_Moffett_Field_centers.csv";
//		// via estimation on the LR level
//		paths->fname_SRF_estimated   = paths->dir_in + "/" + "SRFs"     + "/" + "SRF_QB_MS_gridded_to_Aviris_Moffett_Field_centers_est_with_offset_LR_A0_B0_G0.csv";
//		paths->fname_ImX_shifted     = paths->dir_in + "/" + "QB_MSHR_bandwiseShifted_via_SRF_est" + "/" + "Aviris_Moffett_Field_360x360_QB_MSHR_SNR35_bandwise_shifted_via_LR_est.dat";
//		// via estimation on the HR level
////		paths->fname_SRF_estimated   = paths->dir_in + "/" + "SRFs"     + "/" + "SRF_QB_MS_gridded_to_Aviris_Moffett_Field_centers_est_with_offset_HR_ref_A0_B0_G0.csv";
////		paths->fname_ImX_shifted     = paths->dir_in + "/" + "QB_MSHR_bandwiseShifted_via_SRF_est" + "/" + "Aviris_Moffett_Field_360x360_QB_MSHR_SNR35_bandwise_shifted_via_HR_ref_est.dat";
//		// <===
//		paths->fname_SRF_for_Spectral_Grouping = paths->dir_in + "/" + "SRFs"     + "/" + "SRF_QB_MS_gridded_to_Aviris_Moffett_Field_centers_mod.csv";
//		paths->dir_tmp = "/gpfs/scratch/pr45ne/ga39yoz2/JSparseFI/tmp";
//		//break;
//	}
//	//####################################################################################################################################################################################################################
//	//#  platform     - orig. sensor    - LR sensor       - HR sensor       - fusion type - filter kernel - scene            - size ID       - fDS         - SNR          - redundant digit (>0 for non-regular datasets #
//	//#  1 (SuperMUC) - 8 (HYDICE)      - 8 (HYDICE)      - 5 (Quickbird)   - 2 (HS-MS)   - 1 (gauss)     - 1 (Wash.DC Mall) - 4 (420x300)   - 04 (fDS=4)  - 35 (SNR=35db)- 0                                            #
//	//####################################################################################################################################################################################################################
//	else if(paths->dataSetID_str == "1885211404350"){
//		paths->dir_in        = "/gpfs/work/pr45ne/ga39yoz2/data/HYDICE/191band_HYDICE_image_Washington_DC_Mall/420x300/InputData";
////		paths->fname_ImZ_ref         = paths->dir_in + "/" + "HSHR_ref_raw"              + "/" + "HYDICE_191band_Washington_DC_Mall_420x300.dat";
//		paths->fname_ImZ_ref         = paths->dir_in + "/" + "HSHR_ref_denoised"          + "/" + "HYDICE_WashDC_Mall_420x300_HSHR_ref_denoised.dat";
////		paths->fname_ImZ_init_rec    = paths->dir_in + "/" + "../FusionResults"             + "/" + "2885211404350_HYDICE_WashDC_Mall" + "/" + "BayesianSparse" + "/" + "BEST_ID2885211404350_nbSub7"                               + "/" + "2885211404350_rec_.dat";
////		paths->fname_ImZ_init_rec    = paths->dir_in + "/" + "../FusionResults"             + "/" + "2885211404350_HYDICE_WashDC_Mall" + "/" + "CNMF"           + "/" + "2885211404350_CNMF_mod0PInvBackStep0_M30_I2_5_I1_400_BEST" + "/" + "2885211404350_rec";
////		paths->fname_ImZ_init_rec    = paths->dir_in + "/" + "../FusionResults"             + "/" + "2885211404350_HYDICE_WashDC_Mall" + "/" + "MAPSMM"         + "/" + "2885211404350_sz420x300_fDS4_snr800_nkt4_npure5_BEST"      + "/" + "2885211404350_HSHR_rec_MAPSMM.dat";
////		paths->fname_ImZ_init_rec    = paths->dir_in + "/" + "../FusionResults"             + "/" + "2885211404350_HYDICE_WashDC_Mall" + "/" + "JSpFIHM"        + "/" + "after_one_patRec" + "/" + "ID1885211404350_itr1_init1_lXIm10_lYIm1d7_BEST" + "/" + "1885211404350_rec";
//		paths->fname_ImZ_init_rec    = paths->dir_in + "/" + "HSHR_ref_denoised"          + "/" + "HYDICE_WashDC_Mall_420x300_HSHR_ref_denoised.dat";
//		paths->fname_ImZ_init_ImY_US = paths->dir_in + "/" + "HSLR_fDS4_US_via_bilinear" + "/" + "HYDICE_WashDC_Mall_420x300_HSLR_fDS4_US_via_bilinear_SNR35.dat";
//		paths->fname_ImY             = paths->dir_in + "/" + "HSLR_fDS4"                 + "/" + "HYDICE_WashDC_Mall_420x300_HSLR_fDS4_SNR35.dat";
//		paths->fname_ImX             = paths->dir_in + "/" + "QB_MSHR"                   + "/" + "HYDICE_WashDC_Mall_420x300_QB_MSHR_SNR35.dat";
//		paths->dir_out               = "/gpfs/work/pr45ne/ga39yoz2/recResults";
//		paths->fname_ImZ             = "1885211404350_rec";
//		// ===> SRFs (apriori and estimate) and bandwise offset estimation of ImX
//		paths->fname_SRF             = paths->dir_in + "/" + "SRFs"     + "/" + "SRF_QB_MS_gridded_to_HYDICE_WashDC_Mall_centers.csv";
//		// via estimation on the LR level
//		paths->fname_SRF_estimated   = paths->dir_in + "/" + "SRFs"     + "/" + "SRF_QB_MS_gridded_to_HYDICE_WashDC_Mall_centers_est_with_offset_LR_A0_B0_G0.csv";
//		paths->fname_ImX_shifted     = paths->dir_in + "/" + "QB_MSHR_bandwiseShifted_via_SRF_est" + "/" + "HYDICE_WashDC_Mall_420x300_QB_MSHR_SNR35_bandwise_shifted_via_LR_est.dat";
//		// via estimation on the HR level
////		paths->fname_SRF_estimated   = paths->dir_in + "/" + "SRFs"     + "/" + "SRF_QB_MS_gridded_to_HYDICE_WashDC_Mall_centers_est_with_offset_HR_ref_A0_B0_G0.csv";
////		paths->fname_ImX_shifted     = paths->dir_in + "/" + "QB_MSHR_bandwiseShifted_via_SRF_est" + "/" + "HYDICE_WashDC_Mall_420x300_QB_MSHR_SNR35_bandwise_shifted_via_HR_ref_est.dat";
//		// <===
//		paths->fname_SRF_for_Spectral_Grouping = paths->dir_in + "/" + "SRFs"     + "/" + "SRF_QB_MS_gridded_to_HYDICE_WashDC_Mall_centers_mod.csv";
//		paths->dir_tmp = "/gpfs/scratch/pr45ne/ga39yoz2/JSparseFI/tmp";
//		//break;
//	}
//	//####################################################################################################################################################################################################################
//	//#  platform     - orig. sensor    - LR sensor       - HR sensor       - fusion type - filter kernel - scene            - size ID       - fDS         - SNR          - redundant digit (>0 for non-regular datasets #
//	//#  1 (SuperMUC) - 7 (Headwall)    - 7 (Headwall)    - 4 (WorldView-2) - 2 (HS-MS)   - 1 (gauss)     - 1 (Chikusei)     - 1 (540x420)   - 06 (fDS=6)  - 35 (SNR=35db)- 0                                            #
//	//####################################################################################################################################################################################################################
//	else if(paths->dataSetID_str == "1774211106350"){
//		paths->dir_in        = "/gpfs/work/pr45ne/ga39yoz2/data/Headwall_Chikusei/InputData";
////		paths->fname_ImZ_ref         = paths->dir_in + "/" + "HSHR_ref_raw"              + "/" + "20140729_L1_atm_bcor_mosaic_polish_540x420_HSHR_ref_raw.dat";
//		paths->fname_ImZ_ref         = paths->dir_in + "/" + "HSHR_ref_denoised"          + "/" + "Headwall_Chikusei_540x420_HSHR_ref_denoised.dat";
////		paths->fname_ImZ_init_rec    = paths->dir_in + "/" + "../FusionResults"             + "/" + ".dat";
//		paths->fname_ImZ_init_rec    = paths->dir_in + "/" + "HSHR_ref_denoised"          + "/" + "Headwall_Chikusei_540x420_HSHR_ref_denoised.dat";
//		paths->fname_ImZ_init_ImY_US = paths->dir_in + "/" + "HSLR_fDS6_US_via_bilinear" + "/" + "Headwall_Chikusei_540x420_HSLR_fDS6_US_via_bilinear_SNR35.dat";
//		paths->fname_ImY             = paths->dir_in + "/" + "HSLR_fDS6"                 + "/" + "Headwall_Chikusei_540x420_HSLR_fDS6_SNR35.dat";
//		paths->fname_ImX             = paths->dir_in + "/" + "WV2_MSHR"                  + "/" + "Headwall_Chikusei_540x420_WV2_MSHR_SNR35.dat";
//		paths->dir_out               = "/gpfs/work/pr45ne/ga39yoz2/recResults";
//		paths->fname_ImZ             = "1774211106350_rec";
//		// ===> SRFs (apriori and estimate) and bandwise offset estimation of ImX
//		paths->fname_SRF             = paths->dir_in + "/" + "SRFs"     + "/" + "SRF_WV2_MS_gridded_to_Headwall_Chikusei_centers.csv";
//		// via estimation on the LR level
//		paths->fname_SRF_estimated   = paths->dir_in + "/" + "SRFs"     + "/" + "SRF_WV2_MS_gridded_to_Headwall_Chikusei_centers_est_with_offset_LR_A0_B0_G0.csv";
//		paths->fname_ImX_shifted     = paths->dir_in + "/" + "WV2_MSHR_bandwiseShifted_via_SRF_est" + "/" + "Headwall_Chikusei_540x420_WV2_MSHR_SNR35_bandwise_shifted_via_LR_est.dat";
//		// via estimation on the HR level
////		paths->fname_SRF_estimated   = paths->dir_in + "/" + "SRFs"     + "/" + "SRF_WV2_MS_gridded_to_Headwall_Chikusei_centers_est_with_offset_HR_ref_A0_B0_G0.csv";
////		paths->fname_ImX_shifted     = paths->dir_in + "/" + "WV2_MSHR_bandwiseShifted_via_SRF_est" + "/" + "Headwall_Chikusei_540x420_WV2_MSHR_SNR35_bandwise_shifted_via_HR_ref_est.dat";
//		// <===
//		paths->fname_SRF_for_Spectral_Grouping = paths->dir_in + "/" + "SRFs"     + "/" + "SRF_WV2_MS_gridded_to_Headwall_Chikusei_centers.csv";
//		paths->dir_tmp = "/gpfs/scratch/pr45ne/ga39yoz2/JSparseFI/tmp";
//		//break;
//	}
//	//####################################################################################################################################################################################################################
//	//#  platform     - orig. sensor    - LR sensor       - HR sensor       - fusion type - filter kernel - scene            - size ID       - fDS         - SNR          - redundant digit (>0 for non-regular datasets #
//	//#  1 (SuperMUC) - 7 (Headwall)    - 7 (Headwall)    - 4 (WorldView-2) - 2 (HS-MS)   - 1 (gauss)     - 2 (Chikusei n.u.)- 1 (540x420)   - 06 (fDS=6)  - 35 (SNR=35db)- 0                                            #
//	//####################################################################################################################################################################################################################
//	else if(paths->dataSetID_str == "1774212106350"){
//		paths->dir_in        = "/gpfs/work/pr45ne/ga39yoz2/data/Headwall_Chikusei/non_urban/InputData";
////		paths->fname_ImZ_ref         = paths->dir_in + "/" + "HSHR_ref_raw"              + "/" + "20140729_L1_atm_bcor_mosaic_polish_non_urban_540x420_HSHR_ref_raw.dat";
//		paths->fname_ImZ_ref         = paths->dir_in + "/" + "HSHR_ref_denoised"          + "/" + "Headwall_Chikusei_540x420_HSHR_ref_denoised.dat";
////		paths->fname_ImZ_init_rec    = paths->dir_in + "/" + "../FusionResults"             + "/" + ".dat";
//		paths->fname_ImZ_init_rec    = paths->dir_in + "/" + "HSHR_ref_denoised"          + "/" + "Headwall_Chikusei_540x420_HSHR_ref_denoised.dat";
//		paths->fname_ImZ_init_ImY_US = paths->dir_in + "/" + "HSLR_fDS6_US_via_bilinear" + "/" + "Headwall_Chikusei_540x420_HSLR_fDS6_US_via_bilinear_SNR35.dat";
//		paths->fname_ImY             = paths->dir_in + "/" + "HSLR_fDS6"                 + "/" + "Headwall_Chikusei_540x420_HSLR_fDS6_SNR35.dat";
//		paths->fname_ImX             = paths->dir_in + "/" + "WV2_MSHR"                  + "/" + "Headwall_Chikusei_540x420_WV2_MSHR_SNR35.dat";
//		paths->dir_out               = "/gpfs/work/pr45ne/ga39yoz2/recResults";
//		paths->fname_ImZ             = "1774212106350_rec";
//		// ===> SRFs (apriori and estimate) and bandwise offset estimation of ImX
//		paths->fname_SRF             = paths->dir_in + "/" + "SRFs"     + "/" + "SRF_WV2_MS_gridded_to_Headwall_Chikusei_centers.csv";
//		// via estimation on the LR level
//		paths->fname_SRF_estimated   = paths->dir_in + "/" + "SRFs"     + "/" + "SRF_WV2_MS_gridded_to_Headwall_Chikusei_centers_est_with_offset_LR_A0_B0_G0.csv";
//		paths->fname_ImX_shifted     = paths->dir_in + "/" + "WV2_MSHR_bandwiseShifted_via_SRF_est" + "/" + "Headwall_Chikusei_540x420_WV2_MSHR_SNR35_bandwise_shifted_via_LR_est.dat";
//		// via estimation on the HR level
////		paths->fname_SRF_estimated   = paths->dir_in + "/" + "SRFs"     + "/" + "SRF_WV2_MS_gridded_to_Headwall_Chikusei_centers_est_with_offset_HR_ref_A0_B0_G0.csv";
////		paths->fname_ImX_shifted     = paths->dir_in + "/" + "WV2_MSHR_bandwiseShifted_via_SRF_est" + "/" + "Headwall_Chikusei_540x420_WV2_MSHR_SNR35_bandwise_shifted_via_HR_ref_est.dat";
//		// <===
//		paths->fname_SRF_for_Spectral_Grouping = paths->dir_in + "/" + "SRFs" + "/" + "SRF_WV2_MS_gridded_to_Headwall_Chikusei_centers.csv";
//		paths->dir_tmp = "/gpfs/scratch/pr45ne/ga39yoz2/JSparseFI/tmp";
//		//break;
//	}
//	//####################################################################################################################################################################################################################
//	//#  platform     - orig. sensor    - LR sensor       - HR sensor       - fusion type - filter kernel - scene            - size ID       - fDS         - SNR          - redundant digit (>0 for non-regular datasets #
//	//#  1 (SuperMUC) - 2 (HyMap)       - 10 (EnMAP)      - 9 (Sentinel2)   - 2 (HS-MS)   - 1 (gauss)     - 1 (Rodalquilar)  - 1 (261x867)   - 03 (fDS=3)  - 99 (SNR=inf) - 0                                            #
//	//####################################################################################################################################################################################################################
//	else if(paths->dataSetID_str == "12109211103990"){
//		paths->dir_in        = "/gpfs/work/pr45ne/ga39yoz2/data/HyMap_EnMAP_Sentinel2_from_KS/Fusion/full_image/InputData";
//		paths->fname_ImZ_ref         = paths->dir_in + "/" + "HSHR_ref_raw"                  + "/" + "EnMAP_Sentinel2_rodalquilar_reference_BBR.dat";
////		paths->fname_ImZ_init_rec    = paths->dir_in + "/" + "../FusionResults"                 + "/" + ".dat";
//		paths->fname_ImZ_init_rec    = paths->dir_in + "/" + "HSHR_ref_raw"                  + "/" + "EnMAP_Sentinel2_rodalquilar_reference_BBR.dat";
//		paths->fname_ImZ_init_ImY_US = paths->dir_in + "/" + "HSLR_fDS3_raw_US_via_bilinear" + "/" + "EnMAP_rodalquilar_inpaint_BBR_US_via_bilinear.dat";
//		paths->fname_ImY             = paths->dir_in + "/" + "HSLR_fDS3_raw"                 + "/" + "EnMAP_rodalquilar_inpaint_BBR.dat";
//		paths->fname_ImX             = paths->dir_in + "/" + "Sentinel210m_MSHR_raw"         + "/" + "Sentinel2_rodalquilar_inpaint_10m_bands.dat";
//		paths->dir_out               = "/gpfs/work/pr45ne/ga39yoz2/recResults";
//		paths->fname_ImZ             = "12109211103990_rec";
//		// ===> SRFs (apriori and estimate) and bandwise offset estimation of ImX
//		paths->fname_SRF             = paths->dir_in + "/" + "SRFs"     + "/" + "SRF_Sentinel210m_MS_gridded_to_EnMAP_centers.csv";
//		// via estimation on the LR level
//		paths->fname_SRF_estimated   = paths->dir_in + "/" + "SRFs"     + "/" + "SRF_Sentinel210m_MS_gridded_to_EnMAP_centers_est_with_offset_LR_A0_B0_G0.csv";
//		paths->fname_ImX_shifted     = paths->dir_in + "/" + "Sentinel210m_MSHR_raw_bandwiseShifted_via_SRF_est" + "/" + "Sentinel2_rodalquilar_inpaint_10m_bands_bandwise_shifted_via_LR_est.dat";
//		// via estimation on the HR level
////		paths->fname_SRF_estimated   = paths->dir_in + "/" + "SRFs"     + "/" + "SRF_Sentinel210m_MS_gridded_to_EnMAP_centers_est_with_offset_HR_ref_A0_B0_G0.csv";
////		paths->fname_ImX_shifted     = paths->dir_in + "/" + "Sentinel210m_MSHR_raw_bandwiseShifted_via_SRF_est" + "/" + "Sentinel2_rodalquilar_inpaint_10m_bands_bandwise_shifted_via_HR_ref_est.dat";
//		// <===
//		paths->fname_SRF_for_Spectral_Grouping = paths->dir_in + "/" + "SRFs" + "/" + "SRF_Sentinel210m_MS_gridded_to_EnMAP_centers.csv";
//		paths->dir_tmp = "/gpfs/scratch/pr45ne/ga39yoz2/JSparseFI/tmp";
//		//break;
//	}
//	//####################################################################################################################################################################################################################
//	//#  platform     - orig. sensor    - LR sensor       - HR sensor       - fusion type - filter kernel - scene            - size ID       - fDS         - SNR          - redundant digit (>0 for non-regular datasets #
//	//#  1 (SuperMUC) - 11 (GRSSDFC2013)- 11 (GRSSDFC2013)- 9 (Sentinel2)   - 2 (HS-MS)   - 1 (gauss)     - 1 (Housten Univ.)- 1 (320x540)   - 05 (fDS=5)  - 35 (SNR=inf) - 0                                            #
//	//####################################################################################################################################################################################################################
//	else if(paths->dataSetID_str == "111119211105350"){
//		paths->dir_in        = "/gpfs/work/pr45ne/ga39yoz2/data/2013IEEEGRSSDFC/Fusion/InputData";
////		paths->fname_ImZ_ref         = paths->dir_in + "/" + "HSHR_ref_raw"              + "/" + "2013IEEEGRSSDFC_320x540_HSHR_ref_raw.dat";
//		paths->fname_ImZ_ref         = paths->dir_in + "/" + "HSHR_ref_denoised"         + "/" + "2013IEEEGRSSDFC_320x540_HSHR_ref_denoised.dat";
////		paths->fname_ImZ_init_rec    = paths->dir_in + "/" + "../FusionResults"             + "/" + ".dat";
//		paths->fname_ImZ_init_rec    = paths->dir_in + "/" + "HSHR_ref_denoised"         + "/" + "2013IEEEGRSSDFC_320x540_HSHR_ref_denoised.dat";
//		paths->fname_ImZ_init_ImY_US = paths->dir_in + "/" + "HSLR_fDS5_US_via_bilinear" + "/" + "2013IEEEGRSSDFC_320x540_HSLR_fDS5_US_via_bilinear_SNR35.dat";
//		paths->fname_ImY             = paths->dir_in + "/" + "HSLR_fDS5"                 + "/" + "2013IEEEGRSSDFC_320x540_HSLR_fDS5_SNR35.dat";
//		paths->fname_ImX             = paths->dir_in + "/" + "Sentinel210m_MSHR"         + "/" + "2013IEEEGRSSDFC_320x540_Sentinel210m_MSHR_SNR35.dat";
//		paths->dir_out               = "/gpfs/work/pr45ne/ga39yoz2/recResults";
//		paths->fname_ImZ             = "111119211105350_rec";
//		// ===> SRFs (apriori and estimate) and bandwise offset estimation of ImX
//		paths->fname_SRF             = paths->dir_in + "/" + "SRFs"     + "/" + "SRF_Sentinel210m_MS_gridded_to_2013IEEEGRSSDFC_centers.csv";
//		// via estimation on the LR level
//		paths->fname_SRF_estimated   = paths->dir_in + "/" + "SRFs"     + "/" + "SRF_Sentinel210m_MS_gridded_to_2013IEEEGRSSDFC_centers_est_with_offset_LR_A0_B0_G0.csv";
//		paths->fname_ImX_shifted     = paths->dir_in + "/" + "Sentinel210m_MSHR_bandwiseShifted_via_SRF_est" + "/" + "2013IEEEGRSSDFC_320x540_Sentinel210m_MSHR_SNR35_bandwise_shifted_via_LR_est.dat";
//		// via estimation on the HR level
////		paths->fname_SRF_estimated   = paths->dir_in + "/" + "SRFs"     + "/" + "SRF_Sentinel210m_MS_gridded_to_2013IEEEGRSSDFC_centers_est_with_offset_HR_ref_A0_B0_G0.csv";
////		paths->fname_ImX_shifted     = paths->dir_in + "/" + "Sentinel210m_MSHR_bandwiseShifted_via_SRF_est" + "/" + "2013IEEEGRSSDFC_320x540_Sentinel210m_MSHR_SNR35_bandwise_shifted_via_HR_ref_est.dat";
//		// <===
//		paths->fname_SRF_for_Spectral_Grouping = paths->dir_in + "/" + "SRFs" + "/" + "SRF_Sentinel210m_MS_gridded_to_2013IEEEGRSSDFC_centers.csv";
//		paths->dir_tmp = "/gpfs/scratch/pr45ne/ga39yoz2/JSparseFI/tmp";
//		//break;
//	}
//	//####################################################################################################################################################################################################################
//	//#  platform     - orig. sensor    - LR sensor       - HR sensor       - fusion type - filter kernel - scene            - size ID       - fDS         - SNR          - redundant digit (>0 for non-regular datasets #
//	//#  1 (SuperMUC) - 6 (ROSIS)       - 6 (ROSIS)       - 5 (Quickbird)   - 2 (HS-MS)   - 1 (gauss)     - 1 (Pavia Univ.)  - 1 (560x320)   - 08 (fDS=8)  - 35 (SNR=35db)- 0                                            #
//	//####################################################################################################################################################################################################################
//	else if(paths->dataSetID_str == "1665211108350"){
//		paths->dir_in        = "/gpfs/work/pr45ne/ga39yoz2/data/ROSIS/Pavia_University/InputData";
////		paths->fname_ImZ_ref         = paths->dir_in + "/" + "HSHR_ref_raw"              + "/" + "ROSIS_Pavia_University_560x320_HSHR_ref_raw.dat";
//		paths->fname_ImZ_ref         = paths->dir_in + "/" + "HSHR_ref_denoised"          + "/" + "ROSIS_Pavia_University_560x320_HSHR_ref_denoised.dat";
////		paths->fname_ImZ_init_rec    = paths->dir_in + "/" + "../FusionResults"             + "/" + ".dat";
//		paths->fname_ImZ_init_rec    = paths->dir_in + "/" + "HSHR_ref_denoised"          + "/" + "ROSIS_Pavia_University_560x320_HSHR_ref_denoised.dat";
//		paths->fname_ImZ_init_ImY_US = paths->dir_in + "/" + "HSLR_fDS8_US_via_bilinear" + "/" + "ROSIS_Pavia_University_560x320_HSLR_fDS8_US_via_bilinear_SNR35.dat";
//		paths->fname_ImY             = paths->dir_in + "/" + "HSLR_fDS8"                 + "/" + "ROSIS_Pavia_University_560x320_HSLR_fDS8_SNR35.dat";
//		paths->fname_ImX             = paths->dir_in + "/" + "QB_MSHR"                   + "/" + "ROSIS_Pavia_University_560x320_QB_MSHR_SNR35.dat";
//		paths->dir_out               = "/gpfs/work/pr45ne/ga39yoz2/recResults";
//		paths->fname_ImZ             = "1665211108350_rec";
//		// ===> SRFs (apriori and estimate) and bandwise offset estimation of ImX
//		paths->fname_SRF             = paths->dir_in + "/" + "SRFs"     + "/" + "SRF_QB_MS_gridded_to_ROSIS_Pavia_University_centers.csv";
//		// via estimation on the LR level
//		paths->fname_SRF_estimated   = paths->dir_in + "/" + "SRFs"     + "/" + "SRF_QB_MS_gridded_to_ROSIS_Pavia_University_centers_SNR35_est_with_offset_LR_A0_B0_G0.csv";
//		paths->fname_ImX_shifted     = paths->dir_in + "/" + "QB_MSHR_bandwiseShifted_via_SRF_est" + "/" + "ROSIS_Pavia_University_560x320_QB_MSHR_SNR35_bandwise_shifted_via_LR_est.dat";
//		// via estimation on the HR level
////		paths->fname_SRF_estimated   = paths->dir_in + "/" + "SRFs"     + "/" + "SRF_QB_MS_gridded_to_ROSIS_Pavia_University_centers_est_SNR35_with_offset_HR_ref_A0_B0_G0.csv";
////		paths->fname_ImX_shifted     = paths->dir_in + "/" + "QB_MSHR_bandwiseShifted_via_SRF_est" + "/" + "ROSIS_Pavia_University_560x320_QB_MSHR_SNR35_bandwise_shifted_via_HR_ref_est.dat";
//		// <===
//		paths->fname_SRF_for_Spectral_Grouping = paths->dir_in + "/" + "SRFs"     + "/" + "SRF_QB_MS_gridded_to_ROSIS_Pavia_University_centers.csv";
//		paths->dir_tmp = "/gpfs/scratch/pr45ne/ga39yoz2/JSparseFI/tmp";
//		//break;
//	}
//	//####################################################################################################################################################################################################################
//	//#  platform     - orig. sensor    - LR sensor       - HR sensor       - fusion type - filter kernel - scene            - size ID       - fDS         - SNR          - redundant digit (>0 for non-regular datasets #
//	//#  1 (SuperMUC) - 6 (ROSIS)       - 6 (ROSIS)       - 5 (Quickbird)   - 2 (HS-MS)   - 1 (gauss)     - 1 (Pavia Univ.)  - 1 (560x320)   - 08 (fDS=8)  - 99 (SNR=inf) - 0                                            #
//	//####################################################################################################################################################################################################################
//	else if(paths->dataSetID_str == "1665211108990"){
//		paths->dir_in        = "/gpfs/work/pr45ne/ga39yoz2/data/ROSIS/Pavia_University/InputData";
////		paths->fname_ImZ_ref         = paths->dir_in + "/" + "HSHR_ref_raw"              + "/" + "ROSIS_Pavia_University_560x320_HSHR_ref_raw.dat";
//		paths->fname_ImZ_ref         = paths->dir_in + "/" + "HSHR_ref_denoised"          + "/" + "ROSIS_Pavia_University_560x320_HSHR_ref_denoised.dat";
////		paths->fname_ImZ_init_rec    = paths->dir_in + "/" + "../FusionResults"             + "/" + ".dat";
//		paths->fname_ImZ_init_rec    = paths->dir_in + "/" + "HSHR_ref_denoised"          + "/" + "ROSIS_Pavia_University_560x320_HSHR_ref_denoised.dat";
//		paths->fname_ImZ_init_ImY_US = paths->dir_in + "/" + "HSLR_fDS8_US_via_bilinear" + "/" + "ROSIS_Pavia_University_560x320_HSLR_fDS8_US_via_bilinear_SNRinf.dat";
//		paths->fname_ImY             = paths->dir_in + "/" + "HSLR_fDS8"                 + "/" + "ROSIS_Pavia_University_560x320_HSLR_fDS8_SNRinf.dat";
//		paths->fname_ImX             = paths->dir_in + "/" + "QB_MSHR"                   + "/" + "ROSIS_Pavia_University_560x320_QB_MSHR_SNRinf.dat";
//		paths->dir_out               = "/gpfs/work/pr45ne/ga39yoz2/recResults";
//		paths->fname_ImZ             = "1665211108990_rec";
//		// ===> SRFs (apriori and estimate) and bandwise offset estimation of ImX
//		paths->fname_SRF             = paths->dir_in + "/" + "SRFs"     + "/" + "SRF_QB_MS_gridded_to_ROSIS_Pavia_University_centers.csv";
//		// via estimation on the LR level
//		paths->fname_SRF_estimated   = paths->dir_in + "/" + "SRFs"     + "/" + "SRF_QB_MS_gridded_to_ROSIS_Pavia_University_centers_SNRinf_est_with_offset_LR_A0_B0_G0.csv";
//		paths->fname_ImX_shifted     = paths->dir_in + "/" + "QB_MSHR_bandwiseShifted_via_SRF_est" + "/" + "ROSIS_Pavia_University_560x320_QB_MSHR_SNRinf_bandwise_shifted_via_LR_est.dat";
//		// via estimation on the HR level
////		paths->fname_SRF_estimated   = paths->dir_in + "/" + "SRFs"     + "/" + "SRF_QB_MS_gridded_to_ROSIS_Pavia_University_centers_SNRinf_est_with_offset_HR_ref_A0_B0_G0.csv";
////		paths->fname_ImX_shifted     = paths->dir_in + "/" + "QB_MSHR_bandwiseShifted_via_SRF_est" + "/" + "ROSIS_Pavia_University_560x320_QB_MSHR_SNRinf_bandwise_shifted_via_HR_ref_est.dat";
//		// <===
//		paths->fname_SRF_for_Spectral_Grouping = paths->dir_in + "/" + "SRFs"     + "/" + "SRF_QB_MS_gridded_to_ROSIS_Pavia_University_centers.csv";
//
//		//break;
//	}
//
//////	//####################################################################################################################################################################################################################
//////	//#  platform     - orig. sensor    - LR sensor       - HR sensor       - fusion type - filter kernel - scene            - size ID       - fDS         - SNR          - redundant digit (>0 for non-regular datasets #
//////	//#  1 (SuperMUC) - 1 (HySpex)      - 1 (HySpex)      - 4 (WorldView-2) - 2 (HS-MS)   - 1 (gauss)     - 1 (MUC_Oly)      - 5 (540x540)   - 04 (fDS=4)  - 35 (SNR=35db)- 0                                            #
//////	//####################################################################################################################################################################################################################
//////	else if(paths->dataSetID_str == "1114211504350"){
//////		paths->dir_in        = "/gpfs/work/pr45ne/ga39yoz2/data/HySpex_MUC/Olymp/posished_2012_2m/gauss/540x540";
//////		paths->fname_ImZ_ref         = paths->dir_in + "/" + "HSHR"                      + "/" + "Muenchen_28_VNIR_FOVx2_raw_rad_atm_polish_geo_540x540.dat";
//////		paths->fname_ImZ_init_rec    = paths->dir_in + "/" + "../FusionResults"             + "/" + ".dat"; // TO BE DONE
//////		paths->fname_ImZ_init_ImY_US = paths->dir_in + "/" + "HSLR_fDS4_US_via_bilinear" + "/" + "HySpex_VNIR_540x540_HSLR_fDS4_US_via_bilinear_SNR35.dat";
//////		paths->fname_ImY             = paths->dir_in + "/" + "HSLR_fDS4"                 + "/" + "HySpex_VNIR_540x540_HSLR_fDS4_SNR35.dat";
//////		paths->fname_ImX             = paths->dir_in + "/" + "WV2_MSHR"                  + "/" + "HySpex_VNIR_540x540_WV2_MSHR_SNR35.dat";
//////		paths->dir_out               = "/gpfs/work/pr45ne/ga39yoz2/recResults";
//////		paths->fname_ImZ             = "1114211504350_rec";
//////		paths->fname_SRF                          = "/gpfs/work/pr45ne/ga39yoz2/data/SRF/high_precision/SRF_WV2_MS_gridded_to_HySpex_VNIR_centers.csv";
//////		paths->fname_SRF_estimated                = "/gpfs/work/pr45ne/ga39yoz2/data/SRF/high_precision/estimated/SRF_WV2_MS_gridded_to_HySpex_VNIR_centers_est_LR_A0_B0_G0.csv";
//////		paths->fname_SRF_for_Spectral_Grouping    = "/gpfs/work/pr45ne/ga39yoz2/data/SRF/high_precision/SRF_WV2_MS_gridded_to_HySpex_VNIR_centers.csv";
//////		paths->dir_tmp       = "/gpfs/scratch/pr45ne/ga39yoz2/JSparseFI/tmp";
//////		//break;
//////	}
//////	//####################################################################################################################################################################################################################
//////	//#  platform     - orig. sensor    - LR sensor       - HR sensor       - fusion type - filter kernel - scene            - size ID       - fDS         - SNR          - redundant digit (>0 for non-regular datasets #
//////	//#  1 (SuperMUC) - 1 (HySpex)      - 1 (HySpex)      - 4 (WorldView-2) - 2 (HS-MS)   - 1 (gauss)     - 1 (MUC_Oly)      - 5 (540x540)   - 10 (fDS=10) - 35 (SNR=35db)- 0                                            #
//////	//####################################################################################################################################################################################################################
//////	else if(paths->dataSetID_str == "1114211510350"){
//////		paths->dir_in        = "/gpfs/work/pr45ne/ga39yoz2/data/HySpex_MUC/Olymp/posished_2012_2m/gauss/540x540";
//////		paths->fname_ImZ_ref         = paths->dir_in + "/" + "HSHR"                       + "/" + "Muenchen_28_VNIR_FOVx2_raw_rad_atm_polish_geo_540x540.dat";
////////		paths->fname_ImZ_init_rec    = paths->dir_in + "/" + "../FusionResults"             + "/" + "2114211510350_HySpex_fDS10" + "/" + "BayesianSparse" + "/" + "BEST_ID2114211510350_nbSub10"                              + "/" + "2114211510350_rec_.dat";
////////		paths->fname_ImZ_init_rec    = paths->dir_in + "/" + "../FusionResults"             + "/" + "2114211510350_HySpex_fDS10" + "/" + "CNMF"           + "/" + "2114211510350_CNMF_mod0PInvBackStep0_M40_I2_2_I1_400_BEST" + "/" + "2114211510350_rec";
////////		paths->fname_ImZ_init_rec    = paths->dir_in + "/" + "../FusionResults"             + "/" + "2114211510350_HySpex_fDS10" + "/" + "MAPSMM"         + "/" + "2114211510350_sz540x540_fDS10_snr800_nkt4_npure4_BEST"     + "/" + "2114211510350_HSHR_rec_MAPSMM.dat";
//////		paths->fname_ImZ_init_rec    = paths->dir_in + "/" + "../FusionResults"             + "/" + "2114211510350_HySpex_fDS10" + "/" + "JSpFIHM"        + "/" + "after_one_patRec" + "/" + "ID1114211510350_itr1_init1_lXIm10_lYIm10_BEST" + "/" + "1114211510350_rec";
//////		paths->fname_ImZ_init_ImY_US = paths->dir_in + "/" + "HSLR_fDS10_US_via_bilinear" + "/" + "HySpex_VNIR_540x540_HSLR_fDS10_US_via_bilinear_SNR35.dat";
//////		paths->fname_ImY             = paths->dir_in + "/" + "HSLR_fDS10"                 + "/" + "HySpex_VNIR_540x540_HSLR_fDS10_SNR35.dat";
//////		paths->fname_ImX             = paths->dir_in + "/" + "WV2_MSHR"                   + "/" + "HySpex_VNIR_540x540_WV2_MSHR_SNR35.dat";
//////		paths->dir_out               = "/gpfs/work/pr45ne/ga39yoz2/recResults";
//////		paths->fname_ImZ             = "1114211510350_rec";
//////		paths->fname_SRF                          = "/gpfs/work/pr45ne/ga39yoz2/data/SRF/high_precision/SRF_WV2_MS_gridded_to_HySpex_VNIR_centers.csv";
//////		paths->fname_SRF_estimated                = "/gpfs/work/pr45ne/ga39yoz2/data/SRF/high_precision/estimated/SRF_WV2_MS_gridded_to_HySpex_VNIR_centers_est_LR_A0_B0_G0.csv";
//////		paths->fname_SRF_for_Spectral_Grouping    = "/gpfs/work/pr45ne/ga39yoz2/data/SRF/high_precision/SRF_WV2_MS_gridded_to_HySpex_VNIR_centers.csv";
//////		paths->dir_tmp       = "/gpfs/scratch/pr45ne/ga39yoz2/JSparseFI/tmp";
//////		//break;
//////	}
////	//####################################################################################################################################################################################################################
////		//#  platform     - orig. sensor    - LR sensor       - HR sensor       - fusion type - filter kernel - scene            - size ID       - fDS         - SNR          - redundant digit (>0 for non-regular datasets #
////		//#  1 (SuperMUC) - 1 (HySpex)      - 1 (HySpex)      - 4 (WorldView-2) - 2 (HS-MS)   - 1 (gauss)     - 1 (MUC_Oly)      - 5 (540x540)   - 10 (fDS=10) - 35 (SNR=35db)- 0                                            #
////		//####################################################################################################################################################################################################################
////		else if(paths->dataSetID_str == "1114211510350"){
////			paths->dir_in        = "/gpfs/work/pr45ne/ga39yoz2/data/HySpex_MUC/Olymp/posished_2012_2m/gauss/540x540/InputData";
////	//		paths->fname_ImZ_ref         = paths->dir_in + "/" + "HSHR_ref_raw"               + "/" + "Muenchen_28_VNIR_FOVx2_raw_rad_atm_polish_geo_540x540_raw.dat";
////			paths->fname_ImZ_ref         = paths->dir_in + "/" + "HSHR_ref_denoised"          + "/" + "HySpex_VNIR_540x540_HSHR_ref_denoised.dat";
////	//		paths->fname_ImZ_init_rec    = paths->dir_in + "/" + "../FusionResults"              + "/" + "2114211510350_HySpex_fDS10" + "/" + "BayesianSparse" + "/" + "BEST_ID2114211510350_nbSub10"                              + "/" + "2114211510350_rec_.dat";
////	//		paths->fname_ImZ_init_rec    = paths->dir_in + "/" + "../FusionResults"              + "/" + "2114211510350_HySpex_fDS10" + "/" + "CNMF"           + "/" + "2114211510350_CNMF_mod0PInvBackStep0_M40_I2_2_I1_400_BEST" + "/" + "2114211510350_rec";
////	//		paths->fname_ImZ_init_rec    = paths->dir_in + "/" + "../FusionResults"              + "/" + "2114211510350_HySpex_fDS10" + "/" + "MAPSMM"         + "/" + "2114211510350_sz540x540_fDS10_snr800_nkt4_npure4_BEST"     + "/" + "2114211510350_HSHR_rec_MAPSMM.dat";
////	//		paths->fname_ImZ_init_rec    = paths->dir_in + "/" + "../FusionResults"              + "/" + "2114211510350_HySpex_fDS10" + "/" + "JSpFIHM"        + "/" + "after_one_patRec" + "/" + "ID1114211510350_itr1_init1_lXIm10_lYIm10_BEST" + "/" + "1114211510350_rec";
////			paths->fname_ImZ_init_rec    = paths->dir_in + "/" + "HSHR_ref_denoised"          + "/" + "HySpex_VNIR_540x540_HSHR_ref_denoised.dat";
////			paths->fname_ImZ_init_ImY_US = paths->dir_in + "/" + "HSLR_fDS10_US_via_bilinear" + "/" + "HySpex_VNIR_540x540_HSLR_fDS10_US_via_bilinear_SNR35.dat";
////			paths->fname_ImY             = paths->dir_in + "/" + "HSLR_fDS10"                 + "/" + "HySpex_VNIR_540x540_HSLR_fDS10_SNR35.dat";
////			paths->fname_ImX             = paths->dir_in + "/" + "WV2_MSHR"                   + "/" + "HySpex_VNIR_540x540_WV2_MSHR_SNR35.dat";
////			paths->dir_out               = "/gpfs/work/pr45ne/ga39yoz2/recResults";
////			paths->fname_ImZ             = "1114211510350_rec";
////			paths->fname_SRF                       = paths->dir_in + "/" + "SRFs"     + "/" + "SRF_WV2_MS_gridded_to_HySpex_VNIR_centers.csv";
////			paths->fname_SRF_estimated             = paths->dir_in + "/" + "SRFs"     + "/" + "SRF_WV2_MS_gridded_to_HySpex_VNIR_centers_est_LR_A0_B0_G0.csv";
////			paths->fname_SRF_for_Spectral_Grouping = paths->dir_in + "/" + "SRFs"     + "/" + "SRF_WV2_MS_gridded_to_HySpex_VNIR_centers.csv";
////			paths->dir_tmp       = "/gpfs/scratch/pr45ne/ga39yoz2/JSparseFI/tmp";
////			//break;
////		}
////
////
//////	//####################################################################################################################################################################################################################
//////	//#  platform     - orig. sensor    - LR sensor       - HR sensor       - fusion type - filter kernel - scene            - size ID       - fDS         - SNR          - redundant digit (>0 for non-regular datasets #
//////	//#  1 (SuperMUC) - 1 (HySpex)      - 1 (HySpex)      - 4 (WorldView-2) - 2 (HS-MS)   - 1 (gauss)     - 1 (MUC_Oly)      - 5 (540x540)   - 15 (fDS=15) - 35 (SNR=35db)- 0                                            #
//////	//####################################################################################################################################################################################################################
//////	else if(paths->dataSetID_str == "1114211515350"){
//////		paths->dir_in        = "/gpfs/work/pr45ne/ga39yoz2/data/HySpex_MUC/Olymp/posished_2012_2m/gauss/540x540";
//////		paths->fname_ImZ_ref         = paths->dir_in + "/" + "HSHR"                       + "/" + "Muenchen_28_VNIR_FOVx2_raw_rad_atm_polish_geo_540x540.dat";
//////		paths->fname_ImZ_init_rec    = paths->dir_in + "/" + "../FusionResults"             + "/" + ".dat"; // TO BE DONE
//////		paths->fname_ImZ_init_ImY_US = paths->dir_in + "/" + "HSLR_fDS15_US_via_bilinear" + "/" + "HySpex_VNIR_540x540_HSLR_fDS15_US_via_bilinear_SNR35.dat";
//////		paths->fname_ImY             = paths->dir_in + "/" + "HSLR_fDS15"                 + "/" + "HySpex_VNIR_540x540_HSLR_fDS15_SNR35.dat";
//////		paths->fname_ImX             = paths->dir_in + "/" + "WV2_MSHR"                   + "/" + "HySpex_VNIR_540x540_WV2_MSHR_SNR35.dat";
//////		paths->dir_out               = "/gpfs/work/pr45ne/ga39yoz2/recResults";
//////		paths->fname_ImZ             = "1114211515350_rec";
//////		paths->fname_SRF                          = "/gpfs/work/pr45ne/ga39yoz2/data/SRF/high_precision/SRF_WV2_MS_gridded_to_HySpex_VNIR_centers.csv";
//////		paths->fname_SRF_estimated                = "/gpfs/work/pr45ne/ga39yoz2/data/SRF/high_precision/estimated/SRF_WV2_MS_gridded_to_HySpex_VNIR_centers_est_LR_A0_B0_G0.csv";
//////		paths->fname_SRF_for_Spectral_Grouping    = "/gpfs/work/pr45ne/ga39yoz2/data/SRF/high_precision/SRF_WV2_MS_gridded_to_HySpex_VNIR_centers.csv";
//////		paths->dir_tmp       = "/gpfs/scratch/pr45ne/ga39yoz2/JSparseFI/tmp";
//////		//break;
//////	}
////	//####################################################################################################################################################################################################################
////	//#  platform     - orig. sensor    - LR sensor       - HR sensor       - fusion type - filter kernel - scene            - size ID       - fDS         - SNR          - redundant digit (>0 for non-regular datasets #
////	//#  1 (SuperMUC) - 3 (Aviris)      - 3 (Aviris)      - 4 (WorldView-2) - 2 (HS-MS)   - 1 (gauss)     - 1 (Ind.Pines)    - 3 (360x360)   - 04 (fDS=4)  - 35 (SNR=35db)- 0                                            #
////	//####################################################################################################################################################################################################################
////	else if(paths->dataSetID_str == "1334211304350"){
////		paths->dir_in        = "/gpfs/work/pr45ne/ga39yoz2/data/Aviris/Indian_Pine/220Band_AVIRIS_12June_1992_Indian_Pine_Test_Site_3/supporting/aviris_hyperspectral_data/NS-line_360x360/waterAbsBands_removed/InputData";
//////		paths->fname_ImZ_ref         = paths->dir_in + "/" + "HSHR_ref_raw"              + "/" + "19920612_AVIRIS_IndianPine_NS-line_360x360_BBR.dat";
////		paths->fname_ImZ_ref         = paths->dir_in + "/" + "HSHR_ref_denoised"         + "/" + "Aviris_IndianPines_360x360_HSHR_ref_denoised.dat";
//////		paths->fname_ImZ_init_rec    = paths->dir_in + "/" + "../FusionResults"             + "/" + "2334211304350_Aviris_IndianPines" + "/" + "BayesianSparse" + "/" + "BEST_ID2334211304350_nbSub35"                              + "/" + "2334211304350_rec_";
//////		paths->fname_ImZ_init_rec    = paths->dir_in + "/" + "../FusionResults"             + "/" + "2334211304350_Aviris_IndianPines" + "/" + "CNMF"           + "/" + "2334211304350_CNMF_mod0PInvBackStep0_M20_I2_1_I1_400_BEST" + "/" + "2334211304350_rec";
//////		paths->fname_ImZ_init_rec    = paths->dir_in + "/" + "../FusionResults"             + "/" + "2334211304350_Aviris_IndianPines" + "/" + "MAPSMM"         + "/" + "2334211304350_sz360x360_fDS4_snr800_nkt4_npure4_BEST"      + "/" + "2334211304350_HSHR_rec_MAPSMM.dat";
//////		paths->fname_ImZ_init_rec    = paths->dir_in + "/" + "../FusionResults"             + "/" + "2334211304350_Aviris_IndianPines" + "/" + "JSpFIHM"        + "/" + "after_one_patRec" + "/" + "ID1334211304350_itr1_init1_lXIm31d6_lYIm1_BEST" + "/" + "1334211304350_rec";
////		paths->fname_ImZ_init_rec    = paths->dir_in + "/" + "HSHR_ref_denoised"         + "/" + "Aviris_IndianPines_360x360_HSHR_ref_denoised.dat";
////		paths->fname_ImZ_init_ImY_US = paths->dir_in + "/" + "HSLR_fDS4_US_via_bilinear" + "/" + "Aviris_IndianPines_360x360_HSLR_fDS4_US_via_bilinear_SNR35.dat";
////		paths->fname_ImY             = paths->dir_in + "/" + "HSLR_fDS4"                 + "/" + "Aviris_IndianPines_360x360_HSLR_fDS4_SNR35.dat";
////		paths->fname_ImX             = paths->dir_in + "/" + "WV2_MSHR"                  + "/" + "Aviris_IndianPines_360x360_WV2_MSHR_SNR35.dat";
////		paths->dir_out               = "/gpfs/work/pr45ne/ga39yoz2/recResults";
////		paths->fname_ImZ             = "1334211304350_rec";
////		paths->fname_SRF                       = paths->dir_in + "/" + "SRFs"     + "/" + "SRF_WV2_MS_gridded_to_Aviris_IndianPines_centers.csv";
////		paths->fname_SRF_estimated             = paths->dir_in + "/" + "SRFs"     + "/" + "SRF_WV2_MS_gridded_to_Aviris_IndianPines_centers_est_LR_A0_B0_G0.csv";
////		paths->fname_SRF_for_Spectral_Grouping = paths->dir_in + "/" + "SRFs"     + "/" + "SRF_WV2_MS_gridded_to_Aviris_IndianPines_centers_mod.csv";
////		paths->dir_tmp       = "/gpfs/scratch/pr45ne/ga39yoz2/JSparseFI/tmp";
////		//break;
////	}
////	//####################################################################################################################################################################################################################
////	//#  platform     - orig. sensor    - LR sensor       - HR sensor       - fusion type - filter kernel - scene            - size ID       - fDS         - SNR          - redundant digit (>0 for non-regular datasets #
////	//#  1 (SuperMUC) - 3 (Aviris)      - 3 (Aviris)      - 4 (WorldView-2) - 2 (HS-MS)   - 1 (gauss)     - 2 (Cuprite Sc03) - 4 (420x360)   - 04 (fDS=4)  - 35 (SNR=35db)- 0                                            #
////	//####################################################################################################################################################################################################################
////	else if(paths->dataSetID_str == "1334212404350"){
////		paths->dir_in        = "/gpfs/work/pr45ne/ga39yoz2/data/Aviris/cuprite_from_aviris_website_f970619t01p02r02c_rfl/sc03_a_rfl_420x360/waterAbsBands_removed/InputData";
//////		paths->fname_ImZ_ref         = paths->dir_in + "/" + "HSHR_ref_raw"              + "/" + "f970619t01p02_r02_sc03_a_rfl_420x360.dat";
////		paths->fname_ImZ_ref         = paths->dir_in + "/" + "HSHR_ref_denoised"         + "/" + "Aviris_Cuprite_sc03_420x360_HSHR_ref_denoised.dat";
//////		paths->fname_ImZ_init_rec    = paths->dir_in + "/" + "../FusionResults"             + "/" + "2334212404350_Aviris_Cuprite_sc03" + "/" + "BayesianSparse" + "/" + "BEST_ID2334212404350_nbSub52"                              + "/" + "2334212404350_rec_.dat";
//////		paths->fname_ImZ_init_rec    = paths->dir_in + "/" + "../FusionResults"             + "/" + "2334212404350_Aviris_Cuprite_sc03" + "/" + "CNMF"           + "/" + "2334212404350_CNMF_mod0PInvBackStep0_M20_I2_1_I1_400_BEST" + "/" + "2334212404350_rec";
//////		paths->fname_ImZ_init_rec    = paths->dir_in + "/" + "../FusionResults"             + "/" + "2334212404350_Aviris_Cuprite_sc03" + "/" + "MAPSMM"         + "/" + "2334212404350_sz420x360_fDS4_snr800_nkt6_npure5_BEST"      + "/" + "2334212404350_HSHR_rec_MAPSMM.dat";
//////		paths->fname_ImZ_init_rec    = paths->dir_in + "/" + "../FusionResults"             + "/" + "2334212404350_Aviris_Cuprite_sc03" + "/" + "JSpFIHM"        + "/" + "after_one_patRec" + "/" + "ID1334212404350_itr1_init1_lXIm0d1_lYIm1d7_BEST" + "/" + "1334212404350_rec";
////		paths->fname_ImZ_init_rec    = paths->dir_in + "/" + "HSHR_ref_denoised"         + "/" + "Aviris_Cuprite_sc03_420x360_HSHR_ref_denoised.dat";
////		paths->fname_ImZ_init_ImY_US = paths->dir_in + "/" + "HSLR_fDS4_US_via_bilinear" + "/" + "Aviris_Cuprite_sc03_420x360_HSLR_fDS4_US_via_bilinear_SNR35.dat";
////		paths->fname_ImY             = paths->dir_in + "/" + "HSLR_fDS4"                 + "/" + "Aviris_Cuprite_sc03_420x360_HSLR_fDS4_SNR35.dat";
////		paths->fname_ImX             = paths->dir_in + "/" + "WV2_MSHR"                  + "/" + "Aviris_Cuprite_sc03_420x360_WV2_MSHR_SNR35.dat";
////		paths->dir_out               = "/gpfs/work/pr45ne/ga39yoz2/recResults";
////		paths->fname_ImZ             = "1334212404350_rec";
////		paths->fname_SRF                       = paths->dir_in + "/" + "SRFs"     + "/" + "SRF_WV2_MS_gridded_to_Aviris_Cuprite_sc03_centers.csv";
////		paths->fname_SRF_estimated             = paths->dir_in + "/" + "SRFs"     + "/" + "SRF_WV2_MS_gridded_to_Aviris_Cuprite_sc03_centers_est_without_offset_LR_A0_B0_G0.csv";
////		paths->fname_SRF_for_Spectral_Grouping = paths->dir_in + "/" + "SRFs"     + "/" + "SRF_WV2_MS_gridded_to_Aviris_Cuprite_sc03_centers_mod.csv";
////		paths->dir_tmp       = "/gpfs/scratch/pr45ne/ga39yoz2/JSparseFI/tmp";
////		//break;
////	}
////	//####################################################################################################################################################################################################################
////	//#  platform     - orig. sensor    - LR sensor       - HR sensor       - fusion type - filter kernel - scene            - size ID       - fDS         - SNR          - redundant digit (>0 for non-regular datasets #
////	//#  1 (SuperMUC)  - 3 (Aviris)      - 3 (Aviris)      - 5 (Quickbird)   - 2 (HS-MS)   - 1 (gauss)     - 3 (Moffett Field)- 3 (360x360)   - 04 (fDS=4)  - 35 (SNR=35db)- 0                                            #
////	//####################################################################################################################################################################################################################
////	else if(paths->dataSetID_str == "1335213304350"){
////		paths->dir_in        = "/gpfs/work/pr45ne/ga39yoz2/data/Aviris/Moffet_Field_f080611t01p00r07rdn_c/sc01_ort_img_360x360/waterAbsBands_removed/InputData";
//////		paths->fname_ImZ_ref         = paths->dir_in + "/" + "HSHR_ref_raw"              + "/" + "f080611t01p00r07rdn_c_sc01_ort_img_360x360_BBR.dat";
////		paths->fname_ImZ_ref         = paths->dir_in + "/" + "HSHR_ref_denoised"         + "/" + "Aviris_Moffett_Field_360x360_HSHR_ref_denoised.dat";
//////		paths->fname_ImZ_init_rec    = paths->dir_in + "/" + "../FusionResults"             + "/" + "2335213304350_Aviris_Moffett_Field" + "/" + "BayesianSparse" + "/" + "BEST_ID2335213304350_nbSub13"                              + "/" + "2335213304350_rec_.dat";
//////		paths->fname_ImZ_init_rec    = paths->dir_in + "/" + "../FusionResults"             + "/" + "2335213304350_Aviris_Moffett_Field" + "/" + "CNMF"           + "/" + "2335213304350_CNMF_mod0PInvBackStep0_M30_I2_1_I1_400_BEST" + "/" + "2335213304350_rec";
//////		paths->fname_ImZ_init_rec    = paths->dir_in + "/" + "../FusionResults"             + "/" + "2335213304350_Aviris_Moffett_Field" + "/" + "MAPSMM"         + "/" + "2335213304350_sz360x360_fDS4_snr800_nkt7_npure5_BEST"      + "/" + "2335213304350_HSHR_rec_MAPSMM.dat";
//////		paths->fname_ImZ_init_rec    = paths->dir_in + "/" + "../FusionResults"             + "/" + "2335213304350_Aviris_Moffett_Field" + "/" + "JSpFIHM"        + "/" + "after_one_patRec" + "/" + "ID1335213304350_itr1_init1_lXIm10_lYIm3d16_BEST" + "/" +  "1335213304350_rec";
////		paths->fname_ImZ_init_rec    = paths->dir_in + "/" + "HSHR_ref_denoised"         + "/" + "Aviris_Moffett_Field_360x360_HSHR_ref_denoised.dat";
////		paths->fname_ImZ_init_ImY_US = paths->dir_in + "/" + "HSLR_fDS4_US_via_bilinear" + "/" + "Aviris_Moffett_Field_360x360_HSLR_fDS4_US_via_bilinear_SNR35.dat";
////		paths->fname_ImY             = paths->dir_in + "/" + "HSLR_fDS4"                 + "/" + "Aviris_Moffett_Field_360x360_HSLR_fDS4_SNR35.dat";
////		paths->fname_ImX             = paths->dir_in + "/" + "QB_MSHR"                   + "/" + "Aviris_Moffett_Field_360x360_QB_MSHR_SNR35.dat";
////		paths->dir_out               = "/gpfs/work/pr45ne/ga39yoz2/recResults";
////		paths->fname_ImZ             = "1335213304350_rec";
////		paths->fname_SRF                       = paths->dir_in + "/" + "SRFs"     + "/" + "SRF_QB_MS_gridded_to_Aviris_Moffett_Field_centers.csv";
////		paths->fname_SRF_estimated             = paths->dir_in + "/" + "SRFs"     + "/" + "SRF_QB_MS_gridded_to_Aviris_Moffett_Field_centers_est_without_offset_LR_A0_B0_G0.csv";
////		paths->fname_SRF_for_Spectral_Grouping = paths->dir_in + "/" + "SRFs"     + "/" + "SRF_QB_MS_gridded_to_Aviris_Moffett_Field_centers_mod.csv";
////		paths->dir_tmp       = "/gpfs/scratch/pr45ne/ga39yoz2/JSparseFI/tmp";
////		//break;
////	}
////	//####################################################################################################################################################################################################################
////	//#  platform     - orig. sensor    - LR sensor       - HR sensor       - fusion type - filter kernel - scene            - size ID       - fDS         - SNR          - redundant digit (>0 for non-regular datasets #
////	//#  1 (SuperMUC)  - 8 (HYDICE)      - 8 (HYDICE)      - 5 (Quickbird)   - 2 (HS-MS)   - 1 (gauss)     - 1 (Wash.DC Mall) - 4 (420x300)   - 04 (fDS=4)  - 35 (SNR=35db)- 0                                            #
////	//####################################################################################################################################################################################################################
////	else if(paths->dataSetID_str == "1885211404350"){
////		paths->dir_in        = "/gpfs/work/pr45ne/ga39yoz2/data/HYDICE/191band_HYDICE_image_Washington_DC_Mall/420x300/InputData";
//////		paths->fname_ImZ_ref         = paths->dir_in + "/" + "HSHR_ref_raw"              + "/" + "HYDICE_191band_Washington_DC_Mall_420x300.dat";
////		paths->fname_ImZ_ref         = paths->dir_in + "/" + "HSHR_ref_denoised"          + "/" + "HYDICE_WashDC_Mall_420x300_HSHR_ref_denoised.dat";
//////		paths->fname_ImZ_init_rec    = paths->dir_in + "/" + "../FusionResults"             + "/" + "2885211404350_HYDICE_WashDC_Mall" + "/" + "BayesianSparse" + "/" + "BEST_ID2885211404350_nbSub7"                               + "/" + "2885211404350_rec_.dat";
//////		paths->fname_ImZ_init_rec    = paths->dir_in + "/" + "../FusionResults"             + "/" + "2885211404350_HYDICE_WashDC_Mall" + "/" + "CNMF"           + "/" + "2885211404350_CNMF_mod0PInvBackStep0_M30_I2_5_I1_400_BEST" + "/" + "2885211404350_rec";
//////		paths->fname_ImZ_init_rec    = paths->dir_in + "/" + "../FusionResults"             + "/" + "2885211404350_HYDICE_WashDC_Mall" + "/" + "MAPSMM"         + "/" + "2885211404350_sz420x300_fDS4_snr800_nkt4_npure5_BEST"      + "/" + "2885211404350_HSHR_rec_MAPSMM.dat";
//////		paths->fname_ImZ_init_rec    = paths->dir_in + "/" + "../FusionResults"             + "/" + "2885211404350_HYDICE_WashDC_Mall" + "/" + "JSpFIHM"        + "/" + "after_one_patRec" + "/" + "ID1885211404350_itr1_init1_lXIm10_lYIm1d7_BEST" + "/" + "1885211404350_rec";
////		paths->fname_ImZ_init_rec    = paths->dir_in + "/" + "HSHR_ref_denoised"          + "/" + "HYDICE_WashDC_Mall_420x300_HSHR_ref_denoised.dat";
////		paths->fname_ImZ_init_ImY_US = paths->dir_in + "/" + "HSLR_fDS4_US_via_bilinear" + "/" + "HYDICE_WashDC_Mall_420x300_HSLR_fDS4_US_via_bilinear_SNR35.dat";
////		paths->fname_ImY             = paths->dir_in + "/" + "HSLR_fDS4"                 + "/" + "HYDICE_WashDC_Mall_420x300_HSLR_fDS4_SNR35.dat";
////		paths->fname_ImX             = paths->dir_in + "/" + "QB_MSHR"                   + "/" + "HYDICE_WashDC_Mall_420x300_QB_MSHR_SNR35.dat";
////		paths->dir_out               = "/gpfs/work/pr45ne/ga39yoz2/recResults";
////		paths->fname_ImZ             = "1885211404350_rec";
////		paths->fname_SRF                       = paths->dir_in + "/" + "SRFs"     + "/" + "SRF_QB_MS_gridded_to_HYDICE_WashDC_Mall_centers.csv";
////		paths->fname_SRF_estimated             = paths->dir_in + "/" + "SRFs"     + "/" + "SRF_QB_MS_gridded_to_HYDICE_WashDC_Mall_centers_est_without_offset_LR_A0_B0_G0.csv";
////		paths->fname_SRF_for_Spectral_Grouping = paths->dir_in + "/" + "SRFs"     + "/" + "SRF_QB_MS_gridded_to_HYDICE_WashDC_Mall_centers_mod.csv";
////		paths->dir_tmp       = "/gpfs/scratch/pr45ne/ga39yoz2/JSparseFI/tmp";
////		//break;
////	}
////	//####################################################################################################################################################################################################################
////	//#  platform     - orig. sensor    - LR sensor       - HR sensor       - fusion type - filter kernel - scene            - size ID       - fDS         - SNR          - redundant digit (>0 for non-regular datasets #
////	//#  1 (SuperMUC) - 7 (Headwall)    - 7 (Headwall)    - 4 (WorldView-2) - 2 (HS-MS)   - 1 (gauss)     - 1 (Chikusei)     - 1 (540x420)   - 06 (fDS=6)  - 35 (SNR=35db)- 0                                            #
////	//####################################################################################################################################################################################################################
////	else if(paths->dataSetID_str == "1774211106350"){
////		paths->dir_in        = "/gpfs/work/pr45ne/ga39yoz2/data/Headwall_Chikusei/InputData";
//////		paths->fname_ImZ_ref         = paths->dir_in + "/" + "HSHR_ref_raw"              + "/" + "20140729_L1_atm_bcor_mosaic_polish_540x420_HSHR_ref_raw.dat";
////		paths->fname_ImZ_ref         = paths->dir_in + "/" + "HSHR_ref_denoised"          + "/" + "Headwall_Chikusei_540x420_HSHR_ref_denoised.dat";
//////		paths->fname_ImZ_init_rec    = paths->dir_in + "/" + "../FusionResults"             + "/" + ".dat";
////		paths->fname_ImZ_init_rec    = paths->dir_in + "/" + "HSHR_ref_denoised"          + "/" + "Headwall_Chikusei_540x420_HSHR_ref_denoised.dat";
////		paths->fname_ImZ_init_ImY_US = paths->dir_in + "/" + "HSLR_fDS6_US_via_bilinear" + "/" + "Headwall_Chikusei_540x420_HSLR_fDS6_US_via_bilinear_SNR35.dat";
////		paths->fname_ImY             = paths->dir_in + "/" + "HSLR_fDS6"                 + "/" + "Headwall_Chikusei_540x420_HSLR_fDS6_SNR35.dat";
////		paths->fname_ImX             = paths->dir_in + "/" + "WV2_MSHR"                  + "/" + "Headwall_Chikusei_540x420_WV2_MSHR_SNR35.dat";
////		paths->dir_out               = "/gpfs/work/pr45ne/ga39yoz2/recResults";
////		paths->fname_ImZ             = "1774211106350_rec";
////		paths->fname_SRF                       = paths->dir_in + "/" + "SRFs"     + "/" + "SRF_WV2_MS_gridded_to_Headwall_Chikusei_centers.csv";
////		paths->fname_SRF_estimated             = paths->dir_in + "/" + "SRFs"     + "/" + "SRF_WV2_MS_gridded_to_Headwall_Chikusei_centers_est_without_offset_LR_A0_B0_G0.csv";
////		paths->fname_SRF_for_Spectral_Grouping = paths->dir_in + "/" + "SRFs"     + "/" + "SRF_WV2_MS_gridded_to_Headwall_Chikusei_centers.csv";
////		paths->dir_tmp       = "/gpfs/scratch/pr45ne/ga39yoz2/JSparseFI/tmp";
////		//break;
////	}
////	//####################################################################################################################################################################################################################
////	//#  platform     - orig. sensor    - LR sensor       - HR sensor       - fusion type - filter kernel - scene            - size ID       - fDS         - SNR          - redundant digit (>0 for non-regular datasets #
////	//#  1 (SuperMUC) - 6 (ROSIS)       - 6 (ROSIS)       - 5 (Quickbird)   - 2 (HS-MS)   - 1 (gauss)     - 1 (Pavia Univ.)  - 1 (560x320)   - 08 (fDS=8)  - 35 (SNR=35db)- 0                                            #
////	//####################################################################################################################################################################################################################
////	else if(paths->dataSetID_str == "1665211108350"){
////		paths->dir_in        = "/gpfs/work/pr45ne/ga39yoz2/data/ROSIS/Pavia_University/InputData";
//////		paths->fname_ImZ_ref         = paths->dir_in + "/" + "HSHR_ref_raw"              + "/" + "ROSIS_Pavia_University_560x320_HSHR_ref_raw.dat";
////		paths->fname_ImZ_ref         = paths->dir_in + "/" + "HSHR_ref_denoised"          + "/" + "ROSIS_Pavia_University_560x320_HSHR_ref_denoised.dat";
//////		paths->fname_ImZ_init_rec    = paths->dir_in + "/" + "../FusionResults"             + "/" + ".dat";
////		paths->fname_ImZ_init_rec    = paths->dir_in + "/" + "HSHR_ref_denoised"          + "/" + "ROSIS_Pavia_University_560x320_HSHR_ref_denoised.dat";
////		paths->fname_ImZ_init_ImY_US = paths->dir_in + "/" + "HSLR_fDS8_US_via_bilinear" + "/" + "ROSIS_Pavia_University_560x320_HSLR_fDS8_US_via_bilinear_SNR35.dat";
////		paths->fname_ImY             = paths->dir_in + "/" + "HSLR_fDS8"                 + "/" + "ROSIS_Pavia_University_560x320_HSLR_fDS8_SNR35.dat";
////		paths->fname_ImX             = paths->dir_in + "/" + "QB_MSHR"                   + "/" + "ROSIS_Pavia_University_560x320_QB_MSHR_SNR35.dat";
////		paths->dir_out               = "/gpfs/work/pr45ne/ga39yoz2/recResults";
////		paths->fname_ImZ             = "1665211108350_rec";
////		paths->fname_SRF                       = paths->dir_in + "/" + "SRFs"     + "/" + "SRF_QB_MS_gridded_to_ROSIS_Pavia_University_centers.csv";
////		paths->fname_SRF_estimated             = paths->dir_in + "/" + "SRFs"     + "/" + "SRF_QB_MS_gridded_to_ROSIS_Pavia_University_centers_est_without_offset_LR_A0_B0_G0.csv";
////		paths->fname_SRF_for_Spectral_Grouping = paths->dir_in + "/" + "SRFs"     + "/" + "SRF_QB_MS_gridded_to_ROSIS_Pavia_University_centers.csv";
////		paths->dir_tmp       = "/gpfs/scratch/pr45ne/ga39yoz2/JSparseFI/tmp";
////		//break;
////	}
////
////	//####################################################################################################################################################################################################################
////	//#  platform     - orig. sensor    - LR sensor       - HR sensor       - fusion type - filter kernel - scene            - size ID       - fDS         - SNR          - redundant digit (>0 for non-regular datasets #
////	//#  1 (SuperMUC) - 6 (ROSIS)       - 6 (ROSIS)       - 5 (Quickbird)   - 2 (HS-MS)   - 1 (gauss)     - 1 (Pavia Univ.)  - 1 (560x320)   - 08 (fDS=8)  - 99 (SNR=inf) - 0                                            #
////	//####################################################################################################################################################################################################################
////	else if(paths->dataSetID_str == "1665211108990"){
////		paths->dir_in        = "/gpfs/work/pr45ne/ga39yoz2/data/ROSIS/Pavia_University/InputData";
//////		paths->fname_ImZ_ref         = paths->dir_in + "/" + "HSHR_ref_raw"              + "/" + "ROSIS_Pavia_University_560x320_HSHR_ref_raw.dat";
////		paths->fname_ImZ_ref         = paths->dir_in + "/" + "HSHR_ref_denoised"          + "/" + "ROSIS_Pavia_University_560x320_HSHR_ref_denoised.dat";
//////		paths->fname_ImZ_init_rec    = paths->dir_in + "/" + "../FusionResults"             + "/" + ".dat";
////		paths->fname_ImZ_init_rec    = paths->dir_in + "/" + "HSHR_ref_denoised"          + "/" + "ROSIS_Pavia_University_560x320_HSHR_ref_denoised.dat";
////		paths->fname_ImZ_init_ImY_US = paths->dir_in + "/" + "HSLR_fDS8_US_via_bilinear" + "/" + "ROSIS_Pavia_University_560x320_HSLR_fDS8_US_via_bilinear_SNRinf.dat";
////		paths->fname_ImY             = paths->dir_in + "/" + "HSLR_fDS8"                 + "/" + "ROSIS_Pavia_University_560x320_HSLR_fDS8_SNRinf.dat";
////		paths->fname_ImX             = paths->dir_in + "/" + "QB_MSHR"                   + "/" + "ROSIS_Pavia_University_560x320_QB_MSHR_SNRinf.dat";
////		paths->dir_out               = "/gpfs/work/pr45ne/ga39yoz2/recResults";
////		paths->fname_ImZ             = "1665211108990_rec";
////		paths->fname_SRF                       = paths->dir_in + "/" + "SRFs"     + "/" + "SRF_QB_MS_gridded_to_ROSIS_Pavia_University_centers.csv";
////		paths->fname_SRF_estimated             = paths->dir_in + "/" + "SRFs"     + "/" + "SRF_QB_MS_gridded_to_ROSIS_Pavia_University_centers_est_without_offset_LR_A0_B0_G0.csv";
////		paths->fname_SRF_for_Spectral_Grouping = paths->dir_in + "/" + "SRFs"     + "/" + "SRF_QB_MS_gridded_to_ROSIS_Pavia_University_centers.csv";
////		paths->dir_tmp       = "/gpfs/scratch/pr45ne/ga39yoz2/JSparseFI/tmp";
////		//break;
////	}
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//	/*******************************************
//	 *  datasets used in JSpFI TGRS paper      *
//	 *******************************************/
//	else if(paths->dataSetID_str == "1144111210000"){
//	//#################################################################################################################################################################################################################
//	//#  platform     - orig. sensor    - LR sensor       - HR sensor       - fusion type - filter kernel - scene         - size ID       - fDS         - SNR          - redundant digit (>0 for non-regular datasets #
//	//#  1 (SuperMUC) - 1 (HySpex)      - 4 (WorldView-2) - 4 (WorldView-2) - 1 (MS-Pan)  - 1 (gauss)     - 1 (MUC_Oly)   - 2 (3600x1200) - 10 (fDS=10) - 00 (SNR=inf) - 0                                            #
//	//#################################################################################################################################################################################################################
//		paths->dir_in         = "/gpfs/work/pr45ne/ga39yoz2/data/HySpex_MUC/Olymp/gauss/3600x1200";
//		paths->fname_ImZ_ref  = paths->dir_in + "/" + "WV2_MSHR"        + "/" + "HySpex_3600x1200_WV2_MSHR.dat";
//		paths->fname_ImZ_init = paths->dir_in + "/" + "WV2_MSHR"        + "/" + "HySpex_3600x1200_WV2_MSHR.dat";
//		paths->fname_ImY      = paths->dir_in + "/" + "WV2_MSLR_fDS10"  + "/" + "HySpex_3600x1200_WV2_MSLR_fDS10.dat";
//		paths->fname_ImX      = paths->dir_in + "/" + "WV2_PanHR"       + "/" + "HySpex_3600x1200_WV2_PanHR.dat";
//		paths->dir_out        = "/gpfs/work/pr45ne/ga39yoz2/recResults";
//		paths->fname_ImZ      = "1144111210000_rec";
//		paths->fname_SRF                          = "/gpfs/work/pr45ne/ga39yoz2/data/SRF/high_precision/SRF_WV2_Pan_gridded_to_WV2_MS_centers_new.csv";
//		paths->fname_SRF_for_Spectral_Grouping    = "/gpfs/work/pr45ne/ga39yoz2/data/SRF/high_precision/SRF_WV2_Pan_gridded_to_WV2_MS_centers_new.csv";
//		paths->dir_tmp        = "/gpfs/scratch/pr45ne/ga39yoz2/JSparseFI/tmp";
//		//break;
//	}
//	//------> varies SNRs
//			else if(paths->dataSetID_str == "1144111210100"){
//			//#################################################################################################################################################################################################################
//			//#  platform     - orig. sensor    - LR sensor       - HR sensor       - fusion type - filter kernel - scene         - size ID       - fDS         - SNR          - redundant digit (>0 for non-regular datasets #
//			//#  1 (SuperMUC) - 1 (HySpex)      - 4 (WorldView-2) - 4 (WorldView-2) - 1 (MS-Pan)  - 1 (gauss)     - 1 (MUC_Oly)   - 2 (3600x1200) - 10 (fDS=10) - 10 (SNR=10)  - 0                                            #
//			//#################################################################################################################################################################################################################
//				paths->dir_in         = "/gpfs/work/pr45ne/ga39yoz2/data/HySpex_MUC/Olymp/gauss/3600x1200";
//				paths->fname_ImZ_ref = paths->dir_in + "/" + "WV2_MSHR"        + "/" + "HySpex_3600x1200_WV2_MSHR.dat";
//				paths->fname_ImZ_init = paths->dir_in + "/" + "WV2_MSHR"        + "/" + "HySpex_3600x1200_WV2_MSHR.dat";
//				paths->fname_ImY      = paths->dir_in + "/" + "WV2_MSLR_fDS10"  + "/" + "HySpex_3600x1200_WV2_MSLR_fDS10_SNR10.dat";
//				paths->fname_ImX      = paths->dir_in + "/" + "WV2_PanHR"       + "/" + "HySpex_3600x1200_WV2_PanHR.dat";
//				paths->dir_out        = "/gpfs/work/pr45ne/ga39yoz2/recResults";
//				paths->fname_ImZ      = "1144111210100_rec";
//				paths->fname_SRF                          = "/gpfs/work/pr45ne/ga39yoz2/data/SRF/high_precision/SRF_WV2_Pan_gridded_to_WV2_MS_centers_new.csv";
//				paths->fname_SRF_for_Spectral_Grouping    = "/gpfs/work/pr45ne/ga39yoz2/data/SRF/high_precision/SRF_WV2_Pan_gridded_to_WV2_MS_centers_new.csv";
//				paths->dir_tmp        = "/gpfs/scratch/pr45ne/ga39yoz2/JSparseFI/tmp";
//				//break;
//			}
//			else if(paths->dataSetID_str == "1144111210150"){
//			//#################################################################################################################################################################################################################
//			//#  platform     - orig. sensor    - LR sensor       - HR sensor       - fusion type - filter kernel - scene         - size ID       - fDS         - SNR          - redundant digit (>0 for non-regular datasets #
//			//#  1 (SuperMUC) - 1 (HySpex)      - 4 (WorldView-2) - 4 (WorldView-2) - 1 (MS-Pan)  - 1 (gauss)     - 1 (MUC_Oly)   - 2 (3600x1200) - 10 (fDS=10) - 15 (SNR=15)  - 0                                            #
//			//#################################################################################################################################################################################################################
//				paths->dir_in         = "/gpfs/work/pr45ne/ga39yoz2/data/HySpex_MUC/Olymp/l11vnir/gauss/3600x1200";
//				paths->fname_ImZ_ref = paths->dir_in + "/" + "WV2_MSHR"        + "/" + "HySpex_3600x1200_WV2_MSHR.dat";
//				paths->fname_ImZ_init = paths->dir_in + "/" + "WV2_MSHR"        + "/" + "HySpex_3600x1200_WV2_MSHR.dat";
//				paths->fname_ImY      = paths->dir_in + "/" + "WV2_MSLR_fDS10"  + "/" + "HySpex_3600x1200_WV2_MSLR_fDS10_SNR15.dat";
//				paths->fname_ImX      = paths->dir_in + "/" + "WV2_PanHR"       + "/" + "HySpex_3600x1200_WV2_PanHR.dat";
//				paths->dir_out        = "/gpfs/work/pr45ne/ga39yoz2/recResults";
//				paths->fname_ImZ      = "1144111210150_rec";
//				paths->fname_SRF                          = "/gpfs/work/pr45ne/ga39yoz2/data/SRF/high_precision/SRF_WV2_Pan_gridded_to_WV2_MS_centers_new.csv";
//				paths->fname_SRF_for_Spectral_Grouping    = "/gpfs/work/pr45ne/ga39yoz2/data/SRF/high_precision/SRF_WV2_Pan_gridded_to_WV2_MS_centers_new.csv";
//				paths->dir_tmp        = "/gpfs/scratch/pr45ne/ga39yoz2/JSparseFI/tmp";
//				//break;
//			}
//			else if(paths->dataSetID_str == "1144111210200"){
//			//#################################################################################################################################################################################################################
//			//#  platform     - orig. sensor    - LR sensor       - HR sensor       - fusion type - filter kernel - scene         - size ID       - fDS         - SNR          - redundant digit (>0 for non-regular datasets #
//			//#  1 (SuperMUC) - 1 (HySpex)      - 4 (WorldView-2) - 4 (WorldView-2) - 1 (MS-Pan)  - 1 (gauss)     - 1 (MUC_Oly)   - 2 (3600x1200) - 10 (fDS=10) - 20 (SNR=20)  - 0                                            #
//			//#################################################################################################################################################################################################################
//				paths->dir_in         = "/gpfs/work/pr45ne/ga39yoz2/data/HySpex_MUC/Olymp/gauss/3600x1200";
//				paths->fname_ImZ_ref = paths->dir_in + "/" + "WV2_MSHR"        + "/" + "HySpex_3600x1200_WV2_MSHR.dat";
//				paths->fname_ImZ_init = paths->dir_in + "/" + "WV2_MSHR"        + "/" + "HySpex_3600x1200_WV2_MSHR.dat";
//				paths->fname_ImY      = paths->dir_in + "/" + "WV2_MSLR_fDS10"  + "/" + "HySpex_3600x1200_WV2_MSLR_fDS10_SNR20.dat";
//				paths->fname_ImX      = paths->dir_in + "/" + "WV2_PanHR"       + "/" + "HySpex_3600x1200_WV2_PanHR.dat";
//				paths->dir_out        = "/gpfs/work/pr45ne/ga39yoz2/recResults";
//				paths->fname_ImZ      = "1144111210200_rec";
//				paths->fname_SRF                          = "/gpfs/work/pr45ne/ga39yoz2/data/SRF/high_precision/SRF_WV2_Pan_gridded_to_WV2_MS_centers_new.csv";
//				paths->fname_SRF_for_Spectral_Grouping    = "/gpfs/work/pr45ne/ga39yoz2/data/SRF/high_precision/SRF_WV2_Pan_gridded_to_WV2_MS_centers_new.csv";
//				paths->dir_tmp        = "/gpfs/scratch/pr45ne/ga39yoz2/JSparseFI/tmp";
//				//break;
//			}
//			else if(paths->dataSetID_str == "1144111210300"){
//			//#################################################################################################################################################################################################################
//			//#  platform     - orig. sensor    - LR sensor       - HR sensor       - fusion type - filter kernel - scene         - size ID       - fDS         - SNR          - redundant digit (>0 for non-regular datasets #
//			//#  1 (SuperMUC) - 1 (HySpex)      - 4 (WorldView-2) - 4 (WorldView-2) - 1 (MS-Pan)  - 1 (gauss)     - 1 (MUC_Oly)   - 2 (3600x1200) - 10 (fDS=10) - 30 (SNR=30)  - 0                                            #
//			//#################################################################################################################################################################################################################
//				paths->dir_in         = "/gpfs/work/pr45ne/ga39yoz2/data/HySpex_MUC/Olymp/gauss/3600x1200";
//				paths->fname_ImZ_ref = paths->dir_in + "/" + "WV2_MSHR"        + "/" + "HySpex_3600x1200_WV2_MSHR.dat";
//				paths->fname_ImZ_init = paths->dir_in + "/" + "WV2_MSHR"        + "/" + "HySpex_3600x1200_WV2_MSHR.dat";
//				paths->fname_ImY      = paths->dir_in + "/" + "WV2_MSLR_fDS10"  + "/" + "HySpex_3600x1200_WV2_MSLR_fDS10_SNR30.dat";
//				paths->fname_ImX      = paths->dir_in + "/" + "WV2_PanHR"       + "/" + "HySpex_3600x1200_WV2_PanHR.dat";
//				paths->dir_out        = "/gpfs/work/pr45ne/ga39yoz2/recResults";
//				paths->fname_ImZ      = "1144111210300_rec";
//				paths->fname_SRF                          = "/gpfs/work/pr45ne/ga39yoz2/data/SRF/high_precision/SRF_WV2_Pan_gridded_to_WV2_MS_centers_new.csv";
//				paths->fname_SRF_for_Spectral_Grouping    = "/gpfs/work/pr45ne/ga39yoz2/data/SRF/high_precision/SRF_WV2_Pan_gridded_to_WV2_MS_centers_new.csv";
//				paths->dir_tmp        = "/gpfs/scratch/pr45ne/ga39yoz2/JSparseFI/tmp";
//				//break;
//			}
//			else if(paths->dataSetID_str == "1144111210400"){
//			//#################################################################################################################################################################################################################
//			//#  platform     - orig. sensor    - LR sensor       - HR sensor       - fusion type - filter kernel - scene         - size ID       - fDS         - SNR          - redundant digit (>0 for non-regular datasets #
//			//#  1 (SuperMUC) - 1 (HySpex)      - 4 (WorldView-2) - 4 (WorldView-2) - 1 (MS-Pan)  - 1 (gauss)     - 1 (MUC_Oly)   - 2 (3600x1200) - 10 (fDS=10) - 40 (SNR=40)  - 0                                            #
//			//#################################################################################################################################################################################################################
//				paths->dir_in         = "/gpfs/work/pr45ne/ga39yoz2/data/HySpex_MUC/Olymp/gauss/3600x1200";
//				paths->fname_ImZ_ref = paths->dir_in + "/" + "WV2_MSHR"        + "/" + "HySpex_3600x1200_WV2_MSHR.dat";
//				paths->fname_ImZ_init = paths->dir_in + "/" + "WV2_MSHR"        + "/" + "HySpex_3600x1200_WV2_MSHR.dat";
//				paths->fname_ImY      = paths->dir_in + "/" + "WV2_MSLR_fDS10"  + "/" + "HySpex_3600x1200_WV2_MSLR_fDS10_SNR40.dat";
//				paths->fname_ImX      = paths->dir_in + "/" + "WV2_PanHR"       + "/" + "HySpex_3600x1200_WV2_PanHR.dat";
//				paths->dir_out        = "/gpfs/work/pr45ne/ga39yoz2/recResults";
//				paths->fname_ImZ      = "1144111210400_rec";
//				paths->fname_SRF                          = "/gpfs/work/pr45ne/ga39yoz2/data/SRF/high_precision/SRF_WV2_Pan_gridded_to_WV2_MS_centers_new.csv";
//				paths->fname_SRF_for_Spectral_Grouping    = "/gpfs/work/pr45ne/ga39yoz2/data/SRF/high_precision/SRF_WV2_Pan_gridded_to_WV2_MS_centers_new.csv";
//				paths->dir_tmp        = "/gpfs/scratch/pr45ne/ga39yoz2/JSparseFI/tmp";
//				//break;
//			}
//
//	//------> for reconstruction of band 1, with band 2 previously reconstructed via JSM of bands 2~5
//			else if(paths->dataSetID_str == "1144111210002"){
//			//#################################################################################################################################################################################################################
//			//#  platform     - orig. sensor    - LR sensor       - HR sensor       - fusion type - filter kernel - scene         - size ID       - fDS         - SNR          - redundant digit (>0 for non-regular datasets #
//			//#  1 (SuperMUC) - 1 (HySpex)      - 4 (WorldView-2) - 4 (WorldView-2) - 1 (MS-Pan)  - 1 (gauss)     - 1 (MUC_Oly)   - 2 (3600x1200) - 10 (fDS=10) - 00 (SNR=inf) - 2  (with Pan-band replaced by rec. band 2)   #
//			//#################################################################################################################################################################################################################
//			paths->dir_in         = "/gpfs/work/pr45ne/ga39yoz2/data/HySpex_MUC/Olymp/gauss/3600x1200";
//			paths->fname_ImZ_ref = paths->dir_in + "/" + "WV2_MSHR"        + "/" + "HySpex_3600x1200_WV2_MSHR.dat";
//			paths->fname_ImZ_init = paths->dir_in + "/" + "WV2_MSHR"        + "/" + "HySpex_3600x1200_WV2_MSHR.dat";
//			paths->fname_ImY      = paths->dir_in + "/" + "WV2_MSLR_fDS10"  + "/" + "HySpex_3600x1200_WV2_MSLR_fDS10.dat";
//
//			string dir_in_Pan     = paths->dir_in + "/" + "rec" + "/" + "b2_from_fDS10_as_pan";
//			paths->fname_ImX      = dir_in_Pan + "/" + "WV2_PanHR"       + "/" + "1144111210000_rec_band2.dat";
//
//			paths->dir_out        = "/gpfs/work/pr45ne/ga39yoz2/recResults";
//			paths->fname_ImZ      = "1144111210002_rec";
//			paths->fname_SRF                          = "/gpfs/work/pr45ne/ga39yoz2/data/SRF/high_precision/SRF_WV2_Pan_gridded_to_WV2_MS_centers_new.csv";
//			paths->fname_SRF_for_Spectral_Grouping    = "/gpfs/work/pr45ne/ga39yoz2/data/SRF/high_precision/SRF_WV2_Pan_gridded_to_WV2_MS_centers_new.csv";
//			paths->dir_tmp        = "/gpfs/scratch/pr45ne/ga39yoz2/JSparseFI/tmp";
//			//break;
//			}
//	//------> for reconstruction of bands 7~8, with band 6 previously reconstructed via SparseFI
//			else if(paths->dataSetID_str == "1144111210006"){
//			//#################################################################################################################################################################################################################
//			//#  platform     - orig. sensor    - LR sensor       - HR sensor       - fusion type - filter kernel - scene         - size ID       - fDS         - SNR          - redundant digit (>0 for non-regular datasets #
//			//#  1 (SuperMUC) - 1 (HySpex)      - 4 (WorldView-2) - 4 (WorldView-2) - 1 (MS-Pan)  - 1 (gauss)     - 1 (MUC_Oly)   - 2 (3600x1200) - 10 (fDS=10) - 00 (SNR=inf) - 6  (with Pan-band replaced by rec. band 6)   #
//			//#################################################################################################################################################################################################################
//			paths->dir_in         = "/gpfs/work/pr45ne/ga39yoz2/data/HySpex_MUC/Olymp/gauss/3600x1200";
//			paths->fname_ImZ_ref = paths->dir_in + "/" + "WV2_MSHR"        + "/" + "HySpex_3600x1200_WV2_MSHR.dat";
//			paths->fname_ImZ_init = paths->dir_in + "/" + "WV2_MSHR"        + "/" + "HySpex_3600x1200_WV2_MSHR.dat";
//			paths->fname_ImY      = paths->dir_in + "/" + "WV2_MSLR_fDS10"  + "/" + "HySpex_3600x1200_WV2_MSLR_fDS10.dat";
//
//			string dir_in_Pan     = paths->dir_in + "/" + "rec" + "/" + "b6_from_fDS10_as_pan";
//			paths->fname_ImX      = dir_in_Pan + "/" + "WV2_PanHR"       + "/" + "1144111210000_rec_band6.dat";
//
//			paths->dir_out        = "/gpfs/work/pr45ne/ga39yoz2/recResults";
//			paths->fname_ImZ      = "1144111210006_rec";
//			paths->fname_SRF                          = "/gpfs/work/pr45ne/ga39yoz2/data/SRF/high_precision/SRF_WV2_Pan_gridded_to_WV2_MS_centers_new.csv";
//			paths->fname_SRF_for_Spectral_Grouping    = "/gpfs/work/pr45ne/ga39yoz2/data/SRF/high_precision/SRF_WV2_Pan_gridded_to_WV2_MS_centers_new.csv";
//			paths->dir_tmp        = "/gpfs/scratch/pr45ne/ga39yoz2/JSparseFI/tmp";
//			//break;
//			}
//	else if(paths->dataSetID_str == "1144111204000"){
//	//#################################################################################################################################################################################################################
//	//#  platform     - orig. sensor    - LR sensor       - HR sensor       - fusion type - filter kernel - scene         - size ID       - fDS         - SNR          - redundant digit (>0 for non-regular datasets #
//	//#  1 (SuperMUC) - 1 (HySpex)      - 4 (WorldView-2) - 4 (WorldView-2) - 1 (MS-Pan)  - 1 (gauss)     - 1 (MUC_Oly)   - 2 (3600x1200) - 04 (fDS=4)  - 00 (SNR=inf) - 0                                            #
//	//#################################################################################################################################################################################################################
//		paths->dir_in         = "/gpfs/work/pr45ne/ga39yoz2/data/HySpex_MUC/Olymp/gauss/3600x1200";
//		paths->fname_ImZ_ref = paths->dir_in + "/" + "WV2_MSHR"        + "/" + "HySpex_3600x1200_WV2_MSHR.dat";
//		paths->fname_ImZ_init = paths->dir_in + "/" + "WV2_MSHR"        + "/" + "HySpex_3600x1200_WV2_MSHR.dat";
//		paths->fname_ImY      = paths->dir_in + "/" + "WV2_MSLR_fDS4"   + "/" + "HySpex_3600x1200_WV2_MSLR_fDS4.dat";
//		paths->fname_ImX      = paths->dir_in + "/" + "WV2_PanHR"       + "/" + "HySpex_3600x1200_WV2_PanHR.dat";
//		paths->dir_out        = "/gpfs/work/pr45ne/ga39yoz2/recResults";
//		paths->fname_ImZ      = "1144111104000_rec";
//		paths->fname_SRF                          = "/gpfs/work/pr45ne/ga39yoz2/data/SRF/high_precision/SRF_WV2_Pan_gridded_to_WV2_MS_centers_new.csv";
//		paths->fname_SRF_for_Spectral_Grouping    = "/gpfs/work/pr45ne/ga39yoz2/data/SRF/high_precision/SRF_WV2_Pan_gridded_to_WV2_MS_centers_new.csv";
//		paths->dir_tmp        = "/gpfs/scratch/pr45ne/ga39yoz2/JSparseFI/tmp";
//		//break;
//	}
//	//------> for reconstruction of band 1, with band 2 previously reconstructed via JSM of bands 2~5
//		else if(paths->dataSetID_str == "1144111204002"){
//		//#################################################################################################################################################################################################################
//		//#  platform     - orig. sensor    - LR sensor       - HR sensor       - fusion type - filter kernel - scene         - size ID       - fDS         - SNR          - redundant digit (>0 for non-regular datasets #
//		//#  1 (SuperMUC) - 1 (HySpex)      - 4 (WorldView-2) - 4 (WorldView-2) - 1 (MS-Pan)  - 1 (gauss)     - 1 (MUC_Oly)   - 2 (3600x1200) - 04 (fDS=4)  - 00 (SNR=inf) - 2  (with Pan-band replaced by rec. band 2)   #
//		//#################################################################################################################################################################################################################
//			paths->dir_in         = "/gpfs/work/pr45ne/ga39yoz2/data/HySpex_MUC/Olymp/gauss/3600x1200";
//			paths->fname_ImZ_ref = paths->dir_in + "/" + "WV2_MSHR"        + "/" + "HySpex_3600x1200_WV2_MSHR.dat";
//			paths->fname_ImZ_init = paths->dir_in + "/" + "WV2_MSHR"        + "/" + "HySpex_3600x1200_WV2_MSHR.dat";
//			paths->fname_ImY      = paths->dir_in + "/" + "WV2_MSLR_fDS4"   + "/" + "HySpex_3600x1200_WV2_MSLR_fDS4.dat";
//
//			string dir_in_Pan     = paths->dir_in + "/" + "rec" + "/" + "b2_from_fDS4_as_pan";
//			paths->fname_ImX      = dir_in_Pan + "/" + "WV2_PanHR"       + "/" + "1144111204000_rec_band2.dat";
//
//			paths->dir_out        = "/gpfs/work/pr45ne/ga39yoz2/recResults";
//			paths->fname_ImZ      = "1144111104002_rec";
//			paths->fname_SRF                          = "/gpfs/work/pr45ne/ga39yoz2/data/SRF/high_precision/SRF_WV2_Pan_gridded_to_WV2_MS_centers_new.csv";
//			paths->fname_SRF_for_Spectral_Grouping    = "/gpfs/work/pr45ne/ga39yoz2/data/SRF/high_precision/SRF_WV2_Pan_gridded_to_WV2_MS_centers_new.csv";
//			paths->dir_tmp        = "/gpfs/scratch/pr45ne/ga39yoz2/JSparseFI/tmp";
//			//break;
//		}
//	//------> for reconstruction of band 1, with band 6 previously reconstructed via SparseFI
//		else if(paths->dataSetID_str == "1144111204006"){
//		//#################################################################################################################################################################################################################
//		//#  platform     - orig. sensor    - LR sensor       - HR sensor       - fusion type - filter kernel - scene         - size ID       - fDS         - SNR          - redundant digit (>0 for non-regular datasets #
//		//#  1 (SuperMUC) - 1 (HySpex)      - 4 (WorldView-2) - 4 (WorldView-2) - 1 (MS-Pan)  - 1 (gauss)     - 1 (MUC_Oly)   - 2 (3600x1200) - 04 (fDS=4)  - 00 (SNR=inf) - 6  (with Pan-band replaced by rec. band 6)   #
//		//#################################################################################################################################################################################################################
//			paths->dir_in         = "/gpfs/work/pr45ne/ga39yoz2/data/HySpex_MUC/Olymp/gauss/3600x1200";
//			paths->fname_ImZ_ref = paths->dir_in + "/" + "WV2_MSHR"        + "/" + "HySpex_3600x1200_WV2_MSHR.dat";
//			paths->fname_ImZ_init = paths->dir_in + "/" + "WV2_MSHR"        + "/" + "HySpex_3600x1200_WV2_MSHR.dat";
//			paths->fname_ImY      = paths->dir_in + "/" + "WV2_MSLR_fDS4"   + "/" + "HySpex_3600x1200_WV2_MSLR_fDS4.dat";
//
//			string dir_in_Pan     = paths->dir_in + "/" + "rec" + "/" + "b6_from_fDS4_as_pan";
//			paths->fname_ImX      = dir_in_Pan + "/" + "WV2_PanHR"       + "/" + "1144111204000_rec_band6.dat";
//
//			paths->dir_out        = "/gpfs/work/pr45ne/ga39yoz2/recResults";
//			paths->fname_ImZ      = "1144111104006_rec";
//			paths->fname_SRF                          = "/gpfs/work/pr45ne/ga39yoz2/data/SRF/high_precision/SRF_WV2_Pan_gridded_to_WV2_MS_centers_new.csv";
//			paths->fname_SRF_for_Spectral_Grouping    = "/gpfs/work/pr45ne/ga39yoz2/data/SRF/high_precision/SRF_WV2_Pan_gridded_to_WV2_MS_centers_new.csv";
//			paths->dir_tmp        = "/gpfs/scratch/pr45ne/ga39yoz2/JSparseFI/tmp";
//			//break;
//		}
//
//	else if(paths->dataSetID_str == "1444111104000"){
//	//#################################################################################################################################################################################################################
//	//#  platform     - orig. sensor    - LR sensor       - HR sensor       - fusion type - filter kernel - scene         - size ID       - fDS         - SNR          - redundant digit (>0 for non-regular datasets #
//	//#  1 (SuperMUC) - 4 (WorldView-2) - 4 (WorldView-2) - 4 (WorldView-2) - 1 (MS-Pan)  - 1 (gauss)     - 1 (MUC)       - 1 (960x1320)  - 04 (fDS=4)  - 00 (SNR=inf) - 0                                            #
//	//#################################################################################################################################################################################################################
//		paths->dir_in         = "/gpfs/work/pr45ne/ga39yoz2/data/WorldView2/MUC/gauss/960x1320";
//		paths->fname_ImZ_ref = paths->dir_in + "/" + "MSHR"           + "/" + "WV2_MUC_960x1320_MSHR_ref.dat";
//		paths->fname_ImZ_init = paths->dir_in + "/" + "MSHR"           + "/" + "WV2_MUC_960x1320_MSHR_ref.dat";
//		paths->fname_ImY      = paths->dir_in + "/" + "MSLR_fDS4"      + "/" + "WV2_960x1320_MSLR_fDS4.dat";
//		paths->fname_ImX      = paths->dir_in + "/" + "WV2_PanHR"      + "/" + "WV2_960x1320_WV2_PanHR.dat";
//		paths->dir_out        = "/gpfs/work/pr45ne/ga39yoz2/recResults";
//		paths->fname_ImZ      = "1444111104000_rec";
//		paths->fname_SRF                          = "/gpfs/work/pr45ne/ga39yoz2/data/SRF/high_precision/SRF_WV2_Pan_gridded_to_WV2_MS_centers_new.csv";
//		paths->fname_SRF_for_Spectral_Grouping    = "/gpfs/work/pr45ne/ga39yoz2/data/SRF/high_precision/SRF_WV2_Pan_gridded_to_WV2_MS_centers_new.csv";
//		paths->dir_tmp        = "/gpfs/scratch/pr45ne/ga39yoz2/JSparseFI/tmp";
//		//break;
//	}
//	//------> for reconstruction of band 1, with band 2 previously reconstructed via JSM of bands 2~5
//		else if(paths->dataSetID_str == "1444111104002"){
//		//#################################################################################################################################################################################################################
//		//#  platform     - orig. sensor    - LR sensor       - HR sensor       - fusion type - filter kernel - scene         - size ID       - fDS         - SNR          - redundant digit (>0 for non-regular datasets #
//		//#  1 (SuperMUC) - 4 (WorldView-2) - 4 (WorldView-2) - 4 (WorldView-2) - 1 (MS-Pan)  - 1 (gauss)     - 1 (MUC)       - 1 (960x1320)  - 04 (fDS=4)  - 00 (SNR=inf) - 2  (with Pan-band replaced by rec. band 2)   #
//		//#################################################################################################################################################################################################################
//			paths->dir_in         = "/gpfs/work/pr45ne/ga39yoz2/data/WorldView2/MUC/gauss/960x1320";
//			paths->fname_ImZ_ref = paths->dir_in + "/" + "MSHR"           + "/" + "WV2_MUC_960x1320_MSHR_ref.dat";
//			paths->fname_ImZ_init = paths->dir_in + "/" + "MSHR"           + "/" + "WV2_MUC_960x1320_MSHR_ref.dat";
//			paths->fname_ImY      = paths->dir_in + "/" + "MSLR_fDS4"      + "/" + "WV2_960x1320_MSLR_fDS4.dat";
//
//			string dir_in_Pan     = paths->dir_in + "/" + "rec" + "/" + "b2_as_pan";
//			paths->fname_ImX      = dir_in_Pan + "/" + "WV2_PanHR"       + "/" + "1444111104000_rec_band2.dat";
//
//			paths->dir_out        = "/gpfs/work/pr45ne/ga39yoz2/recResults";
//			paths->fname_ImZ      = "1444111104002_rec";
//			paths->fname_SRF                          = "/gpfs/work/pr45ne/ga39yoz2/data/SRF/high_precision/SRF_WV2_Pan_gridded_to_WV2_MS_centers_new.csv";
//			paths->fname_SRF_for_Spectral_Grouping    = "/gpfs/work/pr45ne/ga39yoz2/data/SRF/high_precision/SRF_WV2_Pan_gridded_to_WV2_MS_centers_new.csv";
//			paths->dir_tmp        = "/gpfs/scratch/pr45ne/ga39yoz2/JSparseFI/tmp";
//			//break;
//		}
//	//------> for reconstruction of band 1, with band 6 previously reconstructed via SparseFI
//		else if(paths->dataSetID_str == "1444111104006"){
//		//#################################################################################################################################################################################################################
//		//#  platform     - orig. sensor    - LR sensor       - HR sensor       - fusion type - filter kernel - scene         - size ID       - fDS         - SNR          - redundant digit (>0 for non-regular datasets #
//		//#  1 (SuperMUC) - 4 (WorldView-2) - 4 (WorldView-2) - 4 (WorldView-2) - 1 (MS-Pan)  - 1 (gauss)     - 1 (MUC)       - 1 (960x1320)  - 04 (fDS=4)  - 00 (SNR=inf) - 6  (with Pan-band replaced by rec. band 6)   #
//		//#################################################################################################################################################################################################################
//			paths->dir_in         = "/gpfs/work/pr45ne/ga39yoz2/data/WorldView2/MUC/gauss/960x1320";
//			paths->fname_ImZ_ref = paths->dir_in + "/" + "MSHR"           + "/" + "WV2_MUC_960x1320_MSHR_ref.dat";
//			paths->fname_ImZ_init = paths->dir_in + "/" + "MSHR"           + "/" + "WV2_MUC_960x1320_MSHR_ref.dat";
//			paths->fname_ImY      = paths->dir_in + "/" + "MSLR_fDS4"      + "/" + "WV2_960x1320_MSLR_fDS4.dat";
//
//			string dir_in_Pan     = paths->dir_in + "/" + "rec" + "/" + "b6_as_pan";
//			paths->fname_ImX      = dir_in_Pan + "/" + "WV2_PanHR"       + "/" + "1444111104000_rec_band6.dat";
//
//			paths->dir_out        = "/gpfs/work/pr45ne/ga39yoz2/recResults";
//			paths->fname_ImZ      = "1444111104006_rec";
//			paths->fname_SRF                          = "/gpfs/work/pr45ne/ga39yoz2/data/SRF/high_precision/SRF_WV2_Pan_gridded_to_WV2_MS_centers_new.csv";
//			paths->fname_SRF_for_Spectral_Grouping    = "/gpfs/work/pr45ne/ga39yoz2/data/SRF/high_precision/SRF_WV2_Pan_gridded_to_WV2_MS_centers_new.csv";
//			paths->dir_tmp        = "/gpfs/scratch/pr45ne/ga39yoz2/JSparseFI/tmp";
//			//break;
//		}
//
//        ////*****************************************************************************************************////
//        ////*****************************************************************************************************////
//        ////                                                                                                     ////
//        ////                                                                                                     ////
//        ////                                          Platform: PC CG                                            ////
//        ////                                                                                                     ////
//        ////                                                                                                     ////
//        ////*****************************************************************************************************////
//        ////*****************************************************************************************************////
//	
//// 	//####################################################################################################################################################################################################################
////         //#  platform     - orig. sensor    - LR sensor       - HR sensor       - fusion type - filter kernel - scene            - size ID       - fDS         - SNR          - redundant digit (>0 for non-regular datasets #
////         //#  2 (CG-PC)    - 7 (Headwall)    - 7 (Headwall)    - 4 (WorldView-2) - 2 (HS-MS)   - 1 (gauss)     - 2 (Chikusei n.u.)- 1 (540x420)   - 06 (fDS=6)  - 35 (SNR=35db)- 0                                            #
////         //####################################################################################################################################################################################################################
////         else if(paths->dataSetID_str == "2774212106350"){
////                 paths->dir_in        = "/data/Projects/JSparseFI/Data/Headwall_Chikusei/non_urban/InputData";
////                 paths->fname_ImZ_ref         = paths->dir_in + "/" + "HSHR_ref_denoised"          + "/" + "Headwall_Chikusei_540x420_HSHR_ref_denoised.dat";
////                 paths->fname_ImZ_init_rec    = paths->fname_ImZ_ref; // paths->dir_in + "/" + "../FusionResults"             + "/" + ".dat";
////                 paths->fname_ImZ_init_ImY_US = paths->dir_in + "/" + "HSLR_fDS6_US_via_bilinear" + "/" + "Headwall_Chikusei_540x420_HSLR_fDS6_US_via_bilinear_SNR35.dat";
////                 paths->fname_ImY             = paths->dir_in + "/" + "HSLR_fDS6"                 + "/" + "Headwall_Chikusei_540x420_HSLR_fDS6_SNR35.dat";
////                 paths->fname_ImX             = paths->dir_in + "/" + "WV2_MSHR"                  + "/" + "Headwall_Chikusei_540x420_WV2_MSHR_SNR35.dat";
////                 paths->dir_out               = "recResults";
////                 paths->fname_ImZ             = "2774212106350_rec";
////                 // ===> SRFs (apriori and estimate) and bandwise offset estimation of ImX
////                 paths->fname_SRF             = paths->dir_in + "/" + "SRFs"     + "/" + "SRF_WV2_MS_gridded_to_Headwall_Chikusei_centers.csv";
////                 // via estimation on the LR level
////                 paths->fname_SRF_estimated   = paths->dir_in + "/" + "SRFs"     + "/" + "SRF_WV2_MS_gridded_to_Headwall_Chikusei_centers_est_with_offset_LR_A0_B0_G0.csv";
////                 paths->fname_ImX_shifted     = paths->dir_in + "/" + "WV2_MSHR_bandwiseShifted_via_SRF_est" + "/" + "Headwall_Chikusei_540x420_WV2_MSHR_SNR35_bandwise_shifted_via_LR_est.dat";
////                 // via estimation on the HR level
////   //            paths->fname_SRF_estimated   = paths->dir_in + "/" + "SRFs"     + "/" + "SRF_WV2_MS_gridded_to_Headwall_Chikusei_centers_est_with_offset_HR_ref_A0_B0_G0.csv";
////   //            paths->fname_ImX_shifted     = paths->dir_in + "/" + "WV2_MSHR_bandwiseShifted_via_SRF_est" + "/" + "Headwall_Chikusei_540x420_WV2_MSHR_SNR35_bandwise_shifted_via_HR_ref_est.dat";
////                 // <===
////                 paths->fname_SRF_for_Spectral_Grouping = paths->dir_in + "/" + "SRFs" + "/" + "SRF_WV2_MS_gridded_to_Headwall_Chikusei_centers.csv";
////                 paths->dir_tmp       = "tmp";
////                 //break;
////         }
//
//
//
//
//
//
//	
//	
//
//	//####################################################################################################################################################################################################################
//	//#  platform     - orig. sensor    - LR sensor       - HR sensor       - fusion type - filter kernel - scene            - size ID       - fDS         - SNR          - redundant digit (>0 for non-regular datasets #
//	//#  2 (CG-PC)    - 1 (HySpex)      - 1 (HySpex)      - 4 (WorldView-2) - 2 (HS-MS)   - 1 (gauss)     - 1 (MUC_Oly)      - 5 (540x540)   - 10 (fDS=10) - 35 (SNR=35db)- 0                                            #
//	//####################################################################################################################################################################################################################
//	else if(paths->dataSetID_str == "2114211510350"){
//		paths->dir_in        = "/data/Projects/JSparseFI/Data/HySpex/MUC/Olymp/polished/gauss/540x540/InputData";
//		paths->fname_ImZ_ref         = paths->dir_in + "/" + "HSHR_ref_denoised"          + "/" + "HySpex_VNIR_540x540_HSHR_ref_denoised.dat";
////		paths->fname_ImZ_init_rec    = paths->dir_in + "/" + "../FusionResults"           + "/" + "";
//		paths->fname_ImZ_init_rec    = paths->dir_in + "/" + "HSHR_ref_denoised"          + "/" + "HySpex_VNIR_540x540_HSHR_ref_denoised.dat";
//		paths->fname_ImZ_init_ImY_US = paths->dir_in + "/" + "HSLR_fDS10_US_via_bilinear" + "/" + "HySpex_VNIR_540x540_HSLR_fDS10_US_via_bilinear_SNR35.dat";
//		paths->fname_ImY             = paths->dir_in + "/" + "HSLR_fDS10"                 + "/" + "HySpex_VNIR_540x540_HSLR_fDS10_SNR35.dat";
//		paths->fname_ImX             = paths->dir_in + "/" + "WV2_MSHR"                   + "/" + "HySpex_VNIR_540x540_WV2_MSHR_SNR35.dat";
//		paths->dir_out               = "recResults";
//		paths->fname_ImZ             = "2114211510350_rec";
//		// ===> SRFs (apriori and estimate) and bandwise offset estimation of ImX
//		paths->fname_SRF             = paths->dir_in + "/" + "SRFs"     + "/" + ".csv";
//		// via estimation on the LR level
//		paths->fname_SRF_estimated   = paths->dir_in + "/" + "SRFs"     + "/" + ".csv";
//		paths->fname_ImX_shifted     = paths->dir_in + "/" + "" + "/" + ".dat";
//		// via estimation on the HR level
////		paths->fname_SRF_estimated   = paths->dir_in + "/" + "SRFs"     + "/" + ".csv";
////		paths->fname_ImX_shifted     = paths->dir_in + "/" + "" + "/" + ".dat";
//		paths->fname_SRF_for_Spectral_Grouping = paths->dir_in + "/" + "SRFs"     + "/" + "SRF_WV2_MS_gridded_to_HySpex_VNIR_centers.csv";
//		// <===
//		paths->dir_tmp       = "tmp";
//		//break;
//	}
//	//####################################################################################################################################################################################################################
//	//#  platform     - orig. sensor    - LR sensor       - HR sensor       - fusion type - filter kernel - scene            - size ID       - fDS         - SNR          - redundant digit (>0 for non-regular datasets #
//	//#  2 (CG-PC)    - 3 (Aviris)      - 3 (Aviris)      - 4 (WorldView-2) - 2 (HS-MS)   - 1 (gauss)     - 1 (Ind.Pines)    - 3 (360x360)   - 04 (fDS=4)  - 35 (SNR=35db)- 0                                            #
//	//####################################################################################################################################################################################################################
//	else if(paths->dataSetID_str == "2334211304350"){
//		paths->dir_in        = "/data/Projects/JSparseFI/Data/Aviris/Indian_Pine/220Band_AVIRIS_12June_1992_Indian_Pine_Test_Site_3/supporting/aviris_hyperspectral_data/NS-line_360x360/waterAbsBands_removed/InputData";
//		paths->fname_ImZ_ref         = paths->dir_in + "/" + "HSHR_ref_denoised"         + "/" + "Aviris_IndianPines_360x360_HSHR_ref_denoised.dat";
////		paths->fname_ImZ_init_rec    = paths->dir_in + "/" + "../FusionResults"          + "/" + "";
//		paths->fname_ImZ_init_rec    = paths->dir_in + "/" + "HSHR_ref_denoised"         + "/" + "Aviris_IndianPines_360x360_HSHR_ref_denoised.dat";
//		paths->fname_ImZ_init_ImY_US = paths->dir_in + "/" + "HSLR_fDS4_US_via_bilinear" + "/" + "Aviris_IndianPines_360x360_HSLR_fDS4_US_via_bilinear_SNR35.dat";
//		paths->fname_ImY             = paths->dir_in + "/" + "HSLR_fDS4"                 + "/" + "Aviris_IndianPines_360x360_HSLR_fDS4_SNR35.dat";
//		paths->fname_ImX             = paths->dir_in + "/" + "WV2_MSHR"                  + "/" + "Aviris_IndianPines_360x360_WV2_MSHR_SNR35.dat";
//		paths->dir_out               = "recResults";
//		paths->fname_ImZ             = "2334211304350_rec";
//		// ===> SRFs (apriori and estimate) and bandwise offset estimation of ImX
//		paths->fname_SRF             = paths->dir_in + "/" + "SRFs"     + "/" + "SRF_WV2_MS_gridded_to_Aviris_IndianPines_centers.csv";
//		// via estimation on the LR level
//		paths->fname_SRF_estimated   = paths->dir_in + "/" + "SRFs"     + "/" + "SRF_WV2_MS_gridded_to_Aviris_IndianPines_centers_est_with_offset_LR_A0_B0_G0.csv";
//		paths->fname_ImX_shifted     = paths->dir_in + "/" + "WV2_MSHR_bandwiseShifted_via_SRF_est" + "/" + "Aviris_IndianPines_360x360_WV2_MSHR_SNR35_bandwise_shifted_via_LR_est.dat";
//		// via estimation on the HR level
////		paths->fname_SRF_estimated   = paths->dir_in + "/" + "SRFs"     + "/" + "SRF_WV2_MS_gridded_to_Aviris_IndianPines_centers_est_with_offset_HR_ref_A0_B0_G0.csv";
////		paths->fname_ImX_shifted     = paths->dir_in + "/" + "WV2_MSHR_bandwiseShifted_via_SRF_est" + "/" + "Aviris_IndianPines_360x360_WV2_MSHR_SNR35_bandwise_shifted_via_HR_ref_est.dat";
//		// <===
//		paths->fname_SRF_for_Spectral_Grouping = paths->dir_in + "/" + "SRFs"     + "/" + "SRF_WV2_MS_gridded_to_Aviris_IndianPines_centers_mod.csv";
//		paths->dir_tmp       = "tmp";
//		//break;
//	}
//	//####################################################################################################################################################################################################################
//	//#  platform     - orig. sensor    - LR sensor       - HR sensor       - fusion type - filter kernel - scene            - size ID       - fDS         - SNR          - redundant digit (>0 for non-regular datasets #
//	//#  2 (CG-PC)    - 3 (Aviris)      - 3 (Aviris)      - 4 (WorldView-2) - 2 (HS-MS)   - 1 (gauss)     - 2 (Cuprite Sc03) - 4 (420x360)   - 04 (fDS=4)  - 35 (SNR=35db)- 0                                            #
//	//####################################################################################################################################################################################################################
//	else if(paths->dataSetID_str == "2334212404350"){
//		paths->dir_in        = "/mnt/NTFS/data/Aviris/cuprite_from_aviris_website_f970619t01p02r02c_rfl/sc03_a_rfl_420x360/waterAbsBands_removed";
//		paths->fname_ImZ_ref         = paths->dir_in + "/" + "HSHR_ref_denoised"         + "/" + "Aviris_Cuprite_sc03_420x360_HSHR_ref_denoised.dat";
////		paths->fname_ImZ_init_rec    = paths->dir_in + "/" + "../FusionResults"             + "/" + "";
//		paths->fname_ImZ_init_rec    = paths->dir_in + "/" + "HSHR_ref_denoised"         + "/" + "Aviris_Cuprite_sc03_420x360_HSHR_ref_denoised.dat";
//		paths->fname_ImZ_init_ImY_US = paths->dir_in + "/" + "HSLR_fDS4_US_via_bilinear" + "/" + "Aviris_Cuprite_sc03_420x360_HSLR_fDS4_US_via_bilinear_SNR35.dat";
//		paths->fname_ImY             = paths->dir_in + "/" + "HSLR_fDS4"                 + "/" + "Aviris_Cuprite_sc03_420x360_HSLR_fDS4_SNR35.dat";
//		paths->fname_ImX             = paths->dir_in + "/" + "WV2_MSHR"                  + "/" + "Aviris_Cuprite_sc03_420x360_WV2_MSHR_SNR35.dat";
//		paths->dir_out               = "recResults";
//		paths->fname_ImZ             = "2334212404350_rec";
//		// ===> SRFs (apriori and estimate) and bandwise offset estimation of ImX
//		paths->fname_SRF             = paths->dir_in + "/" + "SRFs"     + "/" + "SRF_WV2_MS_gridded_to_Aviris_Cuprite_sc03_centers.csv";
//		// via estimation on the LR level
//		paths->fname_SRF_estimated   = paths->dir_in + "/" + "SRFs"     + "/" + "SRF_WV2_MS_gridded_to_Aviris_Cuprite_sc03_centers_est_with_offset_LR_A0_B0_G0.csv";
//		paths->fname_ImX_shifted     = paths->dir_in + "/" + "WV2_MSHR_bandwiseShifted_via_SRF_est" + "/" + "Aviris_Cuprite_sc03_420x360_WV2_MSHR_SNR35_bandwise_shifted_via_LR_est.dat";
//		// via estimation on the HR level
////		paths->fname_SRF_estimated   = paths->dir_in + "/" + "SRFs"     + "/" + "SRF_WV2_MS_gridded_to_Aviris_Cuprite_sc03_centers_est_with_offset_HR_ref_A0_B0_G0.csv";
////		paths->fname_ImX_shifted     = paths->dir_in + "/" + "WV2_MSHR_bandwiseShifted_via_SRF_est" + "/" + "Aviris_Cuprite_sc03_420x360_WV2_MSHR_SNR35_bandwise_shifted_via_HR_ref_est.dat";
//		// <===
//		paths->fname_SRF_for_Spectral_Grouping = paths->dir_in + "/" + "SRFs"     + "/" + "SRF_WV2_MS_gridded_to_Aviris_Cuprite_sc03_centers_mod.csv";
//		paths->dir_tmp       = "tmp";
//		//break;
//	}
//	//####################################################################################################################################################################################################################
//	//#  platform     - orig. sensor    - LR sensor       - HR sensor       - fusion type - filter kernel - scene            - size ID       - fDS         - SNR          - redundant digit (>0 for non-regular datasets #
//	//#  2 (CG-PC)    - 3 (Aviris)      - 3 (Aviris)      - 4 (WorldView-2) - 2 (HS-MS)   - 1 (gauss)     - 2 (Cuprite Sc03) - 4 (420x360)   - 05 (fDS=5)  - 35 (SNR=35db)- 0                                            #
//	//####################################################################################################################################################################################################################
//	else if(paths->dataSetID_str == "2334212405350"){
////		paths->dir_in        = "/data/Projects/JSparseFI/Data/Aviris/cuprite_from_aviris_website_f970619t01p02r02c_rfl/sc03_a_rfl_420x360/waterAbsBands_removed_fDS5/InputData";
//		paths->dir_in        = "/mnt/NTFS/data/Aviris/cuprite_from_aviris_website_f970619t01p02r02c_rfl/sc03_a_rfl_420x360/waterAbsBands_removed_fDS5/InputData";
//		paths->fname_ImZ_ref         = paths->dir_in + "/" + "HSHR_ref_denoised"         + "/" + "Aviris_Cuprite_sc03_420x360_HSHR_ref_denoised.dat";
////		paths->fname_ImZ_init_rec    = paths->dir_in + "/" + "../FusionResults"             + "/" + "";
//		paths->fname_ImZ_init_rec    = paths->dir_in + "/" + "HSHR_ref_denoised"         + "/" + "Aviris_Cuprite_sc03_420x360_HSHR_ref_denoised.dat";
//		paths->fname_ImZ_init_ImY_US = paths->dir_in + "/" + "HSLR_fDS5_US_via_bilinear" + "/" + "Aviris_Cuprite_sc03_420x360_HSLR_fDS5_US_via_bilinear_SNR35.dat";
//		paths->fname_ImY             = paths->dir_in + "/" + "HSLR_fDS5"                 + "/" + "Aviris_Cuprite_sc03_420x360_HSLR_fDS5_SNR35.dat";
//		paths->fname_ImX             = paths->dir_in + "/" + "WV2_MSHR"                  + "/" + "Aviris_Cuprite_sc03_420x360_WV2_MSHR_SNR35.dat";
//		paths->dir_out               = "Results";
//		paths->fname_ImZ             = "2334212405350_rec";
//		// ===> SRFs (apriori and estimate) and bandwise offset estimation of ImX
//		paths->fname_SRF             = paths->dir_in + "/" + "SRFs"     + "/" + "SRF_WV2_MS_gridded_to_Aviris_Cuprite_sc03_centers.csv";
//		// via estimation on the LR level
//		paths->fname_SRF_estimated   = paths->dir_in + "/" + "SRFs"     + "/" + "SRF_WV2_MS_gridded_to_Aviris_Cuprite_sc03_centers_est_with_offset_LR_A0_B0_G0.csv";
//		paths->fname_ImX_shifted     = paths->dir_in + "/" + "WV2_MSHR_bandwiseShifted_via_SRF_est" + "/" + "Aviris_Cuprite_sc03_420x360_WV2_MSHR_SNR35_bandwise_shifted_via_LR_est.dat";
//		// via estimation on the HR level
////		paths->fname_SRF_estimated   = paths->dir_in + "/" + "SRFs"     + "/" + "SRF_WV2_MS_gridded_to_Aviris_Cuprite_sc03_centers_est_with_offset_HR_ref_A0_B0_G0.csv";
////		paths->fname_ImX_shifted     = paths->dir_in + "/" + "WV2_MSHR_bandwiseShifted_via_SRF_est" + "/" + "Aviris_Cuprite_sc03_420x360_WV2_MSHR_SNR35_bandwise_shifted_via_HR_ref_est.dat";
//		// <===
//		paths->fname_SRF_for_Spectral_Grouping = paths->dir_in + "/" + "SRFs"     + "/" + "SRF_WV2_MS_gridded_to_Aviris_Cuprite_sc03_centers_mod.csv";
//		paths->dir_tmp       = "tmp";
//		//break;
//	}
//	//####################################################################################################################################################################################################################
//	//#  platform     - orig. sensor    - LR sensor       - HR sensor       - fusion type - filter kernel - scene            - size ID       - fDS         - SNR          - redundant digit (>0 for non-regular datasets #
//	//#  2 (CG-PC)    - 3 (Aviris)      - 3 (Aviris)      - 5 (Quickbird)   - 2 (HS-MS)   - 1 (gauss)     - 3 (Moffett Field)- 3 (360x360)   - 04 (fDS=4)  - 35 (SNR=35db)- 0                                            #
//	//####################################################################################################################################################################################################################
//	else if(paths->dataSetID_str == "2335213304350"){
//	        //paths->dataSetID_str = SSTR(dSetting->platformID) + paths->dataSetID_str;
//		paths->dir_in        = maindir_path + "/" + "2335213304350_Aviris_Moffett_Field" + "/" + "InputData";
//		paths->fname_ImZ_ref         = paths->dir_in + "/" + "links" + "/" + "slink_to_ImZ_ref.dat";// "HSHR_ref_denoised"         + "/" + "Aviris_Moffett_Field_360x360_HSHR_ref_denoised.dat";
//		paths->fname_ImZ_init_rec    = paths->dir_in + "/" + "links" + "/" + "slink_to_ImZ_ref.dat";// "HSHR_ref_denoised"         + "/" + "Aviris_Moffett_Field_360x360_HSHR_ref_denoised.dat";
//		paths->fname_ImZ_init_ImY_US = paths->dir_in + "/" + "links" + "/" + "slink_to_ImZ_init_ImY_US.dat";// "HSLR_fDS4_US_via_bilinear" + "/" + "Aviris_Moffett_Field_360x360_HSLR_fDS4_US_via_bilinear_SNR35.dat";
//		paths->fname_ImY             = paths->dir_in + "/" + "links" + "/" + "slink_to_ImY.dat";// "HSLR_fDS4"                 + "/" + "Aviris_Moffett_Field_360x360_HSLR_fDS4_SNR35.dat";
//		paths->fname_ImX             = paths->dir_in + "/" + "links" + "/" + "slink_to_ImX.dat";// "QB_MSHR"                   + "/" + "Aviris_Moffett_Field_360x360_QB_MSHR_SNR35.dat";
//		paths->fname_ImZ             = paths->dataSetID_str + "_rec";
//		// ===> SRFs (apriori and estimate) and bandwise offset estimation of ImX
//		paths->fname_SRF             = paths->dir_in + "/" + "links" + "/" + "slink_to_SRF.csv";// "SRFs"     + "/" + "SRF_QB_MS_gridded_to_Aviris_Moffett_Field_centers.csv";
//		// via estimation on the LR level
//		paths->fname_SRF_estimated   = paths->dir_in + "/" + "links" + "/" + "slink_to_SRF_estimated.csv";// "SRFs"     + "/" + "SRF_QB_MS_gridded_to_Aviris_Moffett_Field_centers_est_with_offset_LR_A0_B0_G0.csv";
//		paths->fname_ImX_shifted     = paths->dir_in + "/" + "links" + "/" + "slink_to_ImX_shifted.dat";// "QB_MSHR_bandwiseShifted_via_SRF_est" + "/" + "Aviris_Moffett_Field_360x360_QB_MSHR_SNR35_bandwise_shifted_via_LR_est.dat";
//		// via estimation on the HR level
////		paths->fname_SRF_estimated   = paths->dir_in + "/" + "SRFs"     + "/" + "SRF_QB_MS_gridded_to_Aviris_Moffett_Field_centers_est_with_offset_HR_ref_A0_B0_G0.csv";
////		paths->fname_ImX_shifted     = paths->dir_in + "/" + "QB_MSHR_bandwiseShifted_via_SRF_est" + "/" + "Aviris_Moffett_Field_360x360_QB_MSHR_SNR35_bandwise_shifted_via_HR_ref_est.dat";
//		// <===
//		//paths->fname_SRF_for_Spectral_Grouping = paths->dir_in + "/" + "SRFs"     + "/" + "SRF_QB_MS_gridded_to_Aviris_Moffett_Field_centers_mod.csv";
//		paths->fname_SRF_for_Spectral_Grouping = paths->dir_in + "/" + "links" + "/" + "slink_to_SRF_for_Spectral_Grouping.csv";//"/" + "SRFs"     + "/" + "fname_SRF_for_Spectral_Grouping.csv"; //"SRF_QB_MS_gridded_to_Aviris_Moffett_Field_centers.csv";
//		//paths->dir_tmp       = "tmp";
//		//break;
//	}
//	//####################################################################################################################################################################################################################
//	//#  platform     - orig. sensor    - LR sensor       - HR sensor       - fusion type - filter kernel - scene            - size ID       - fDS         - SNR          - redundant digit (>0 for non-regular datasets #
//	//#  2 (CG-PC)    - 8 (HYDICE)      - 8 (HYDICE)      - 5 (Quickbird)   - 2 (HS-MS)   - 1 (gauss)     - 1 (Wash.DC Mall) - 4 (420x300)   - 04 (fDS=4)  - 35 (SNR=35db)- 0                                            #
//	//####################################################################################################################################################################################################################
//	else if(paths->dataSetID_str == "2885211404350"){
//		paths->dir_in        = "/data/Projects/JSparseFI/Data/HYDICE/191band_HYDICE_image_Washington_DC_Mall/420x300/InputData";
////		paths->fname_ImZ_ref         = paths->dir_in + "/" + "HSHR_ref_raw"              + "/" + "HYDICE_191band_Washington_DC_Mall_420x300.dat";
//		paths->fname_ImZ_ref         = paths->dir_in + "/" + "HSHR_ref_denoised"          + "/" + "HYDICE_WashDC_Mall_420x300_HSHR_ref_denoised.dat";
////		paths->fname_ImZ_init_rec    = paths->dir_in + "/" + "../FusionResults"             + "/" + "2885211404350_HYDICE_WashDC_Mall" + "/" + "BayesianSparse" + "/" + "BEST_ID2885211404350_nbSub7"                               + "/" + "2885211404350_rec_.dat";
//		paths->fname_ImZ_init_rec    = paths->dir_in + "/" + "HSHR_ref_denoised"          + "/" + "HYDICE_WashDC_Mall_420x300_HSHR_ref_denoised.dat";
//		paths->fname_ImZ_init_ImY_US = paths->dir_in + "/" + "HSLR_fDS4_US_via_bilinear" + "/" + "HYDICE_WashDC_Mall_420x300_HSLR_fDS4_US_via_bilinear_SNR35.dat";
//		paths->fname_ImY             = paths->dir_in + "/" + "HSLR_fDS4"                 + "/" + "HYDICE_WashDC_Mall_420x300_HSLR_fDS4_SNR35.dat";
//		paths->fname_ImX             = paths->dir_in + "/" + "QB_MSHR"                   + "/" + "HYDICE_WashDC_Mall_420x300_QB_MSHR_SNR35.dat";
//		paths->dir_out               = "recResults";
//		paths->fname_ImZ             = "2885211404350_rec";
//		// ===> SRFs (apriori and estimate) and bandwise offset estimation of ImX
//		paths->fname_SRF             = paths->dir_in + "/" + "SRFs"     + "/" + "SRF_QB_MS_gridded_to_HYDICE_WashDC_Mall_centers.csv";
//		// via estimation on the LR level
//		paths->fname_SRF_estimated   = paths->dir_in + "/" + "SRFs"     + "/" + "SRF_QB_MS_gridded_to_HYDICE_WashDC_Mall_centers_est_with_offset_LR_A0_B0_G0.csv";
//		paths->fname_ImX_shifted     = paths->dir_in + "/" + "QB_MSHR_bandwiseShifted_via_SRF_est" + "/" + "HYDICE_WashDC_Mall_420x300_QB_MSHR_SNR35_bandwise_shifted_via_LR_est.dat";
//		// via estimation on the HR level
////		paths->fname_SRF_estimated   = paths->dir_in + "/" + "SRFs"     + "/" + "SRF_QB_MS_gridded_to_HYDICE_WashDC_Mall_centers_est_with_offset_HR_ref_A0_B0_G0.csv";
////		paths->fname_ImX_shifted     = paths->dir_in + "/" + "QB_MSHR_bandwiseShifted_via_SRF_est" + "/" + "HYDICE_WashDC_Mall_420x300_QB_MSHR_SNR35_bandwise_shifted_via_HR_ref_est.dat";
//		// <===
//		paths->fname_SRF_for_Spectral_Grouping = paths->dir_in + "/" + "SRFs"     + "/" + "SRF_QB_MS_gridded_to_HYDICE_WashDC_Mall_centers_mod.csv";
//		paths->dir_tmp       = "tmp";
//		//break;
//	}
//	//####################################################################################################################################################################################################################
//	//#  platform     - orig. sensor    - LR sensor       - HR sensor       - fusion type - filter kernel - scene            - size ID       - fDS         - SNR          - redundant digit (>0 for non-regular datasets #
//	//#  2 (CG-PC)    - 7 (Headwall)    - 7 (Headwall)    - 4 (WorldView-2) - 2 (HS-MS)   - 1 (gauss)     - 1 (Chikusei)     - 1 (540x420)   - 06 (fDS=6)  - 35 (SNR=35db)- 0                                            #
//	//####################################################################################################################################################################################################################
//	else if(paths->dataSetID_str == "2774211106350"){
//		paths->dir_in        = "/data/Projects/JSparseFI/Data/Headwall_Chikusei/InputData";
////		paths->fname_ImZ_ref         = paths->dir_in + "/" + "HSHR_ref_raw"              + "/" + "20140729_L1_atm_bcor_mosaic_polish_540x420_HSHR_ref_raw.dat";
//		paths->fname_ImZ_ref         = paths->dir_in + "/" + "HSHR_ref_denoised"          + "/" + "Headwall_Chikusei_540x420_HSHR_ref_denoised.dat";
////		paths->fname_ImZ_init_rec    = paths->dir_in + "/" + "../FusionResults"             + "/" + ".dat";
//		paths->fname_ImZ_init_rec    = paths->dir_in + "/" + "HSHR_ref_denoised"          + "/" + "Headwall_Chikusei_540x420_HSHR_ref_denoised.dat";
//		paths->fname_ImZ_init_ImY_US = paths->dir_in + "/" + "HSLR_fDS6_US_via_bilinear" + "/" + "Headwall_Chikusei_540x420_HSLR_fDS6_US_via_bilinear_SNR35.dat";
//		paths->fname_ImY             = paths->dir_in + "/" + "HSLR_fDS6"                 + "/" + "Headwall_Chikusei_540x420_HSLR_fDS6_SNR35.dat";
//		paths->fname_ImX             = paths->dir_in + "/" + "WV2_MSHR"                  + "/" + "Headwall_Chikusei_540x420_WV2_MSHR_SNR35.dat";
//		paths->dir_out               = "recResults";
//		paths->fname_ImZ             = "2774211106350_rec";
//		// ===> SRFs (apriori and estimate) and bandwise offset estimation of ImX
//		paths->fname_SRF             = paths->dir_in + "/" + "SRFs"     + "/" + "SRF_WV2_MS_gridded_to_Headwall_Chikusei_centers.csv";
//		// via estimation on the LR level
//		paths->fname_SRF_estimated   = paths->dir_in + "/" + "SRFs"     + "/" + "SRF_WV2_MS_gridded_to_Headwall_Chikusei_centers_est_with_offset_LR_A0_B0_G0.csv";
//		paths->fname_ImX_shifted     = paths->dir_in + "/" + "WV2_MSHR_bandwiseShifted_via_SRF_est" + "/" + "Headwall_Chikusei_540x420_WV2_MSHR_SNR35_bandwise_shifted_via_LR_est.dat";
//		// via estimation on the HR level
////		paths->fname_SRF_estimated   = paths->dir_in + "/" + "SRFs"     + "/" + "SRF_WV2_MS_gridded_to_Headwall_Chikusei_centers_est_with_offset_HR_ref_A0_B0_G0.csv";
////		paths->fname_ImX_shifted     = paths->dir_in + "/" + "WV2_MSHR_bandwiseShifted_via_SRF_est" + "/" + "Headwall_Chikusei_540x420_WV2_MSHR_SNR35_bandwise_shifted_via_HR_ref_est.dat";
//		// <===
//		paths->fname_SRF_for_Spectral_Grouping = paths->dir_in + "/" + "SRFs"     + "/" + "SRF_WV2_MS_gridded_to_Headwall_Chikusei_centers.csv";
//		paths->dir_tmp       = "tmp";
//		//break;
//	}
//	//####################################################################################################################################################################################################################
//	//#  platform     - orig. sensor    - LR sensor       - HR sensor       - fusion type - filter kernel - scene            - size ID       - fDS         - SNR          - redundant digit (>0 for non-regular datasets #
//	//#  2 (CG-PC)    - 7 (Headwall)    - 7 (Headwall)    - 4 (WorldView-2) - 2 (HS-MS)   - 1 (gauss)     - 2 (Chikusei n.u.)- 1 (540x420)   - 06 (fDS=6)  - 35 (SNR=35db)- 0                                            #
//	//####################################################################################################################################################################################################################
//	else if(paths->dataSetID_str == "2774212106350"){
//		paths->dir_in        = "/data/Projects/JSparseFI/Data/Headwall_Chikusei/non_urban/InputData";
////		paths->fname_ImZ_ref         = paths->dir_in + "/" + "HSHR_ref_raw"              + "/" + "20140729_L1_atm_bcor_mosaic_polish_non_urban_540x420_HSHR_ref_raw.dat";
//		paths->fname_ImZ_ref         = paths->dir_in + "/" + "HSHR_ref_denoised"          + "/" + "Headwall_Chikusei_540x420_HSHR_ref_denoised.dat";
////		paths->fname_ImZ_init_rec    = paths->dir_in + "/" + "../FusionResults"             + "/" + ".dat";
//		paths->fname_ImZ_init_rec    = paths->dir_in + "/" + "HSHR_ref_denoised"          + "/" + "Headwall_Chikusei_540x420_HSHR_ref_denoised.dat";
//		paths->fname_ImZ_init_ImY_US = paths->dir_in + "/" + "HSLR_fDS6_US_via_bilinear" + "/" + "Headwall_Chikusei_540x420_HSLR_fDS6_US_via_bilinear_SNR35.dat";
//		paths->fname_ImY             = paths->dir_in + "/" + "HSLR_fDS6"                 + "/" + "Headwall_Chikusei_540x420_HSLR_fDS6_SNR35.dat";
//		paths->fname_ImX             = paths->dir_in + "/" + "WV2_MSHR"                  + "/" + "Headwall_Chikusei_540x420_WV2_MSHR_SNR35.dat";
//		paths->dir_out               = "recResults";
//		paths->fname_ImZ             = "2774212106350_rec";
//		// ===> SRFs (apriori and estimate) and bandwise offset estimation of ImX
//		paths->fname_SRF             = paths->dir_in + "/" + "SRFs"     + "/" + "SRF_WV2_MS_gridded_to_Headwall_Chikusei_centers.csv";
//		// via estimation on the LR level
//		paths->fname_SRF_estimated   = paths->dir_in + "/" + "SRFs"     + "/" + "SRF_WV2_MS_gridded_to_Headwall_Chikusei_centers_est_with_offset_LR_A0_B0_G0.csv";
//		paths->fname_ImX_shifted     = paths->dir_in + "/" + "WV2_MSHR_bandwiseShifted_via_SRF_est" + "/" + "Headwall_Chikusei_540x420_WV2_MSHR_SNR35_bandwise_shifted_via_LR_est.dat";
//		// via estimation on the HR level
////		paths->fname_SRF_estimated   = paths->dir_in + "/" + "SRFs"     + "/" + "SRF_WV2_MS_gridded_to_Headwall_Chikusei_centers_est_with_offset_HR_ref_A0_B0_G0.csv";
////		paths->fname_ImX_shifted     = paths->dir_in + "/" + "WV2_MSHR_bandwiseShifted_via_SRF_est" + "/" + "Headwall_Chikusei_540x420_WV2_MSHR_SNR35_bandwise_shifted_via_HR_ref_est.dat";
//		// <===
//		paths->fname_SRF_for_Spectral_Grouping = paths->dir_in + "/" + "SRFs" + "/" + "SRF_WV2_MS_gridded_to_Headwall_Chikusei_centers.csv";
//		paths->dir_tmp       = "tmp";
//		//break;
//	}
//	//####################################################################################################################################################################################################################
//	//#  platform     - orig. sensor    - LR sensor       - HR sensor       - fusion type - filter kernel - scene            - size ID       - fDS         - SNR          - redundant digit (>0 for non-regular datasets #
//	//#  2 (CG-PC)    - 2 (HyMap)       - 10 (EnMAP)      - 9 (Sentinel210m)- 2 (HS-MS)   - 1 (gauss)     - 1 (Rodalquilar)  - 1 (261x867)   - 03 (fDS=3)  - 99 (SNR=inf) - 0                                            #
//	//####################################################################################################################################################################################################################
//	else if(paths->dataSetID_str == "22109211103990"){
//		paths->dir_in        = "/data/Projects/JSparseFI/Data/HyMap_EnMAP_Sentinel2_from_KS/Fusion/full_image/InputData";
//		paths->fname_ImZ_ref         = paths->dir_in + "/" + "HSHR_ref_raw"                  + "/" + "EnMAP_Sentinel2_rodalquilar_reference_BBR.dat";
////		paths->fname_ImZ_init_rec    = paths->dir_in + "/" + "../FusionResults"                 + "/" + ".dat";
//		paths->fname_ImZ_init_rec    = paths->dir_in + "/" + "HSHR_ref_raw"                  + "/" + "EnMAP_Sentinel2_rodalquilar_reference_BBR.dat";
//		paths->fname_ImZ_init_ImY_US = paths->dir_in + "/" + "HSLR_fDS3_raw_US_via_bilinear" + "/" + "EnMAP_rodalquilar_inpaint_BBR_US_via_bilinear.dat";
//		paths->fname_ImY             = paths->dir_in + "/" + "HSLR_fDS3_raw"                 + "/" + "EnMAP_rodalquilar_inpaint_BBR.dat";
//		paths->fname_ImX             = paths->dir_in + "/" + "Sentinel210m_MSHR_raw"         + "/" + "Sentinel2_rodalquilar_inpaint_10m_bands.dat";
//		paths->dir_out               = "recResults";
//		paths->fname_ImZ             = "22109211103990_rec";
//		// ===> SRFs (apriori and estimate) and bandwise offset estimation of ImX
//		paths->fname_SRF             = paths->dir_in + "/" + "SRFs"     + "/" + "SRF_Sentinel210m_MS_gridded_to_EnMAP_centers.csv";
//		// via estimation on the LR level
//		paths->fname_SRF_estimated   = paths->dir_in + "/" + "SRFs"     + "/" + "SRF_Sentinel210m_MS_gridded_to_EnMAP_centers_est_with_offset_LR_A0_B0_G0.csv";
//		paths->fname_ImX_shifted     = paths->dir_in + "/" + "Sentinel210m_MSHR_raw_bandwiseShifted_via_SRF_est" + "/" + "Sentinel2_rodalquilar_inpaint_10m_bands_bandwise_shifted_via_LR_est.dat";
//		// via estimation on the HR level
////		paths->fname_SRF_estimated   = paths->dir_in + "/" + "SRFs"     + "/" + "SRF_Sentinel210m_MS_gridded_to_EnMAP_centers_est_with_offset_HR_ref_A0_B0_G0.csv";
////		paths->fname_ImX_shifted     = paths->dir_in + "/" + "Sentinel210m_MSHR_raw_bandwiseShifted_via_SRF_est" + "/" + "Sentinel2_rodalquilar_inpaint_10m_bands_bandwise_shifted_via_HR_ref_est.dat";
//		// <===
//		paths->fname_SRF_for_Spectral_Grouping = paths->dir_in + "/" + "SRFs" + "/" + "SRF_Sentinel210m_MS_gridded_to_EnMAP_centers.csv";
//		paths->dir_tmp       = "tmp";
//		//break;
//	}
//	//####################################################################################################################################################################################################################
//	//#  platform     - orig. sensor    - LR sensor       - HR sensor       - fusion type - filter kernel - scene            - size ID       - fDS         - SNR          - redundant digit (>0 for non-regular datasets #
//	//#  2 (CG-PC)    - 11 (GRSSDFC2013)- 11 (GRSSDFC2013)- 9 (Sentinel210m)- 2 (HS-MS)   - 1 (gauss)     - 1 (Housten Univ.)- 1 (320x540)   - 05 (fDS=5)  - 35 (SNR=inf) - 0                                            #
//	//####################################################################################################################################################################################################################
//	else if(paths->dataSetID_str == "211119211105350"){
//		paths->dir_in        = "/data/Projects/JSparseFI/Data/2013IEEEGRSSDFC/Fusion/InputData";
////		paths->fname_ImZ_ref         = paths->dir_in + "/" + "HSHR_ref_raw"              + "/" + "2013IEEEGRSSDFC_320x540_HSHR_ref_raw.dat";
//		paths->fname_ImZ_ref         = paths->dir_in + "/" + "HSHR_ref_denoised"         + "/" + "2013IEEEGRSSDFC_320x540_HSHR_ref_denoised.dat";
////		paths->fname_ImZ_init_rec    = paths->dir_in + "/" + "../FusionResults"             + "/" + ".dat";
//		paths->fname_ImZ_init_rec    = paths->dir_in + "/" + "HSHR_ref_denoised"         + "/" + "2013IEEEGRSSDFC_320x540_HSHR_ref_denoised.dat";
//		paths->fname_ImZ_init_ImY_US = paths->dir_in + "/" + "HSLR_fDS5_US_via_bilinear" + "/" + "2013IEEEGRSSDFC_320x540_HSLR_fDS5_US_via_bilinear_SNR35.dat";
//		paths->fname_ImY             = paths->dir_in + "/" + "HSLR_fDS5"                 + "/" + "2013IEEEGRSSDFC_320x540_HSLR_fDS5_SNR35.dat";
//		paths->fname_ImX             = paths->dir_in + "/" + "Sentinel210m_MSHR"         + "/" + "2013IEEEGRSSDFC_320x540_Sentinel210m_MSHR_SNR35.dat";
//		paths->dir_out               = "recResults";
//		paths->fname_ImZ             = "211119211105350_rec";
//		// ===> SRFs (apriori and estimate) and bandwise offset estimation of ImX
//		paths->fname_SRF             = paths->dir_in + "/" + "SRFs"     + "/" + "SRF_Sentinel210m_MS_gridded_to_2013IEEEGRSSDFC_centers.csv";
//		// via estimation on the LR level
//		paths->fname_SRF_estimated   = paths->dir_in + "/" + "SRFs"     + "/" + "SRF_Sentinel210m_MS_gridded_to_2013IEEEGRSSDFC_centers_est_with_offset_LR_A0_B0_G0.csv";
//		paths->fname_ImX_shifted     = paths->dir_in + "/" + "Sentinel210m_MSHR_bandwiseShifted_via_SRF_est" + "/" + "2013IEEEGRSSDFC_320x540_Sentinel210m_MSHR_SNR35_bandwise_shifted_via_LR_est.dat";
//		// via estimation on the HR level
////		paths->fname_SRF_estimated   = paths->dir_in + "/" + "SRFs"     + "/" + "SRF_Sentinel210m_MS_gridded_to_2013IEEEGRSSDFC_centers_est_with_offset_HR_ref_A0_B0_G0.csv";
////		paths->fname_ImX_shifted     = paths->dir_in + "/" + "Sentinel210m_MSHR_bandwiseShifted_via_SRF_est" + "/" + "2013IEEEGRSSDFC_320x540_Sentinel210m_MSHR_SNR35_bandwise_shifted_via_HR_ref_est.dat";
//		// <===
//		paths->fname_SRF_for_Spectral_Grouping = paths->dir_in + "/" + "SRFs" + "/" + "SRF_Sentinel210m_MS_gridded_to_2013IEEEGRSSDFC_centers.csv";
//		paths->dir_tmp       = "tmp";
//		//break;
//	}
//	//####################################################################################################################################################################################################################
//	//#  platform     - orig. sensor    - LR sensor       - HR sensor       - fusion type - filter kernel - scene            - size ID       - fDS         - SNR          - redundant digit (>0 for non-regular datasets #
//	//#  2 (CG-PC)    - 6 (ROSIS)       - 6 (ROSIS)       - 5 (Quickbird)   - 2 (HS-MS)   - 1 (gauss)     - 1 (Pavia Univ.)  - 1 (560x320)   - 08 (fDS=8)  - 35 (SNR=35db)- 0                                            #
//	//####################################################################################################################################################################################################################
//	else if(paths->dataSetID_str == "2665211108350"){
//		paths->dir_in        = "/data/Projects/JSparseFI/Data/ROSIS/Pavia_University/InputData";
////		paths->fname_ImZ_ref         = paths->dir_in + "/" + "HSHR_ref_raw"              + "/" + "ROSIS_Pavia_University_560x320_HSHR_ref_raw.dat";
//		paths->fname_ImZ_ref         = paths->dir_in + "/" + "HSHR_ref_denoised"          + "/" + "ROSIS_Pavia_University_560x320_HSHR_ref_denoised.dat";
////		paths->fname_ImZ_init_rec    = paths->dir_in + "/" + "../FusionResults"             + "/" + ".dat";
//		paths->fname_ImZ_init_rec    = paths->dir_in + "/" + "HSHR_ref_denoised"          + "/" + "ROSIS_Pavia_University_560x320_HSHR_ref_denoised.dat";
//		paths->fname_ImZ_init_ImY_US = paths->dir_in + "/" + "HSLR_fDS8_US_via_bilinear" + "/" + "ROSIS_Pavia_University_560x320_HSLR_fDS8_US_via_bilinear_SNR35.dat";
//		paths->fname_ImY             = paths->dir_in + "/" + "HSLR_fDS8"                 + "/" + "ROSIS_Pavia_University_560x320_HSLR_fDS8_SNR35.dat";
//		paths->fname_ImX             = paths->dir_in + "/" + "QB_MSHR"                   + "/" + "ROSIS_Pavia_University_560x320_QB_MSHR_SNR35.dat";
//		paths->dir_out               = "recResults";
//		paths->fname_ImZ             = "2665211108350_rec";
//		// ===> SRFs (apriori and estimate) and bandwise offset estimation of ImX
//		paths->fname_SRF             = paths->dir_in + "/" + "SRFs"     + "/" + "SRF_QB_MS_gridded_to_ROSIS_Pavia_University_centers.csv";
//		// via estimation on the LR level
//		paths->fname_SRF_estimated   = paths->dir_in + "/" + "SRFs"     + "/" + "SRF_QB_MS_gridded_to_ROSIS_Pavia_University_centers_SNR35_est_with_offset_LR_A0_B0_G0.csv";
//		paths->fname_ImX_shifted     = paths->dir_in + "/" + "QB_MSHR_bandwiseShifted_via_SRF_est" + "/" + "ROSIS_Pavia_University_560x320_QB_MSHR_SNR35_bandwise_shifted_via_LR_est.dat";
//		// via estimation on the HR level
////		paths->fname_SRF_estimated   = paths->dir_in + "/" + "SRFs"     + "/" + "SRF_QB_MS_gridded_to_ROSIS_Pavia_University_centers_est_SNR35_with_offset_HR_ref_A0_B0_G0.csv";
////		paths->fname_ImX_shifted     = paths->dir_in + "/" + "QB_MSHR_bandwiseShifted_via_SRF_est" + "/" + "ROSIS_Pavia_University_560x320_QB_MSHR_SNR35_bandwise_shifted_via_HR_ref_est.dat";
//		// <===
//		paths->fname_SRF_for_Spectral_Grouping = paths->dir_in + "/" + "SRFs"     + "/" + "SRF_QB_MS_gridded_to_ROSIS_Pavia_University_centers.csv";
//		paths->dir_tmp       = "tmp";
//		//break;
//	}
//	//####################################################################################################################################################################################################################
//	//#  platform     - orig. sensor    - LR sensor       - HR sensor       - fusion type - filter kernel - scene            - size ID       - fDS         - SNR          - redundant digit (>0 for non-regular datasets #
//	//#  2 (CG-PC)    - 6 (ROSIS)       - 6 (ROSIS)       - 5 (Quickbird)   - 2 (HS-MS)   - 1 (gauss)     - 1 (Pavia Univ.)  - 1 (560x320)   - 08 (fDS=8)  - 99 (SNR=inf) - 0                                            #
//	//####################################################################################################################################################################################################################
//	else if(paths->dataSetID_str == "2665211108990"){
//		paths->dir_in        = "/data/Projects/JSparseFI/Data/ROSIS/Pavia_University/InputData";
////		paths->fname_ImZ_ref         = paths->dir_in + "/" + "HSHR_ref_raw"              + "/" + "ROSIS_Pavia_University_560x320_HSHR_ref_raw.dat";
//		paths->fname_ImZ_ref         = paths->dir_in + "/" + "HSHR_ref_denoised"          + "/" + "ROSIS_Pavia_University_560x320_HSHR_ref_denoised.dat";
////		paths->fname_ImZ_init_rec    = paths->dir_in + "/" + "../FusionResults"             + "/" + ".dat";
//		paths->fname_ImZ_init_rec    = paths->dir_in + "/" + "HSHR_ref_denoised"          + "/" + "ROSIS_Pavia_University_560x320_HSHR_ref_denoised.dat";
//		paths->fname_ImZ_init_ImY_US = paths->dir_in + "/" + "HSLR_fDS8_US_via_bilinear" + "/" + "ROSIS_Pavia_University_560x320_HSLR_fDS8_US_via_bilinear_SNRinf.dat";
//		paths->fname_ImY             = paths->dir_in + "/" + "HSLR_fDS8"                 + "/" + "ROSIS_Pavia_University_560x320_HSLR_fDS8_SNRinf.dat";
//		paths->fname_ImX             = paths->dir_in + "/" + "QB_MSHR"                   + "/" + "ROSIS_Pavia_University_560x320_QB_MSHR_SNRinf.dat";
//		paths->dir_out               = "recResults";
//		paths->fname_ImZ             = "2665211108990_rec";
//		// ===> SRFs (apriori and estimate) and bandwise offset estimation of ImX
//		paths->fname_SRF             = paths->dir_in + "/" + "SRFs"     + "/" + "SRF_QB_MS_gridded_to_ROSIS_Pavia_University_centers.csv";
//		// via estimation on the LR level
//		paths->fname_SRF_estimated   = paths->dir_in + "/" + "SRFs"     + "/" + "SRF_QB_MS_gridded_to_ROSIS_Pavia_University_centers_SNRinf_est_with_offset_LR_A0_B0_G0.csv";
//		paths->fname_ImX_shifted     = paths->dir_in + "/" + "QB_MSHR_bandwiseShifted_via_SRF_est" + "/" + "ROSIS_Pavia_University_560x320_QB_MSHR_SNRinf_bandwise_shifted_via_LR_est.dat";
//		// via estimation on the HR level
////		paths->fname_SRF_estimated   = paths->dir_in + "/" + "SRFs"     + "/" + "SRF_QB_MS_gridded_to_ROSIS_Pavia_University_centers_SNRinf_est_with_offset_HR_ref_A0_B0_G0.csv";
////		paths->fname_ImX_shifted     = paths->dir_in + "/" + "QB_MSHR_bandwiseShifted_via_SRF_est" + "/" + "ROSIS_Pavia_University_560x320_QB_MSHR_SNRinf_bandwise_shifted_via_HR_ref_est.dat";
//		// <===
//		paths->fname_SRF_for_Spectral_Grouping = paths->dir_in + "/" + "SRFs"     + "/" + "SRF_QB_MS_gridded_to_ROSIS_Pavia_University_centers.csv";
//		paths->dir_tmp       = "tmp";
//		//break;
//	}
//
//
//
//
//
//
//
//
//
//
//
//
//
//	/*******************************************
//	 *  datasets used in JSpFI TGRS paper      *
//	 *******************************************/
//	else if(paths->dataSetID_str == "2144111210000"){
//	//#################################################################################################################################################################################################################
//	//#  platform     - orig. sensor    - LR sensor       - HR sensor       - fusion type - filter kernel - scene         - size ID       - fDS         - SNR          - redundant digit (>0 for non-regular datasets #
//	//#  2 (CG-PC)    - 1 (HySpex)      - 4 (WorldView-2) - 4 (WorldView-2) - 1 (MS-Pan)  - 1 (gauss)     - 1 (MUC_Oly)   - 2 (3600x1200) - 10 (fDS=10) - 00 (SNR=inf) - 0                                            #
//	//#################################################################################################################################################################################################################
//		paths->dir_in        = "/data/Projects/JSparseFI/Data/HySpex/MUC/Olymp/l11vnir/gauss/3600x1200";
//		paths->fname_ImZ_ref= paths->dir_in + "/" + "WV2_MSHR"        + "/" + "HySpex_3600x1200_WV2_MSHR.dat";
//		paths->fname_ImZ_init= paths->dir_in + "/" + "WV2_MSHR"        + "/" + "HySpex_3600x1200_WV2_MSHR.dat";
//		paths->fname_ImY     = paths->dir_in + "/" + "WV2_MSLR_fDS10"  + "/" + "HySpex_3600x1200_WV2_MSLR_fDS10.dat";
//		paths->fname_ImX     = paths->dir_in + "/" + "WV2_PanHR"       + "/" + "HySpex_3600x1200_WV2_PanHR.dat";
//		paths->dir_out       = "recResults";
//		paths->fname_ImZ     = "2144111210000_rec";
//		paths->fname_SRF                          = "/data/Projects/JSparseFI/Data/SRF/for_cpp_routine/SRF_WV2_Pan_gridded_to_WV2_MS_centers_new.csv";
//		paths->fname_SRF_for_Spectral_Grouping    = "/data/Projects/JSparseFI/Data/SRF/for_cpp_routine/SRF_WV2_Pan_gridded_to_WV2_MS_centers_new.csv";
//		paths->dir_tmp       = "tmp";
//		//break;
//	}
//	//------> varies SNRs
//			else if(paths->dataSetID_str == "2144111210100"){
//			//#################################################################################################################################################################################################################
//			//#  platform     - orig. sensor    - LR sensor       - HR sensor       - fusion type - filter kernel - scene         - size ID       - fDS         - SNR          - redundant digit (>0 for non-regular datasets #
//			//#  2 (CG-PC)    - 1 (HySpex)      - 4 (WorldView-2) - 4 (WorldView-2) - 1 (MS-Pan)  - 1 (gauss)     - 1 (MUC_Oly)   - 2 (3600x1200) - 10 (fDS=10) - 10 (SNR=10)  - 0                                            #
//			//#################################################################################################################################################################################################################
//				paths->dir_in        = "/data/Projects/JSparseFI/Data/HySpex/MUC/Olymp/l11vnir/gauss/3600x1200";
//				paths->fname_ImZ_ref= paths->dir_in + "/" + "WV2_MSHR"        + "/" + "HySpex_3600x1200_WV2_MSHR.dat";
//				paths->fname_ImZ_init= paths->dir_in + "/" + "WV2_MSHR"        + "/" + "HySpex_3600x1200_WV2_MSHR.dat";
//				paths->fname_ImY     = paths->dir_in + "/" + "WV2_MSLR_fDS10"  + "/" + "HySpex_3600x1200_WV2_MSLR_fDS10_SNR10.dat";
//				paths->fname_ImX     = paths->dir_in + "/" + "WV2_PanHR"       + "/" + "HySpex_3600x1200_WV2_PanHR.dat";
//				paths->dir_out       = "recResults";
//				paths->fname_ImZ     = "2144111210100_rec";
//				paths->fname_SRF                          = "/data/Projects/JSparseFI/Data/SRF/for_cpp_routine/SRF_WV2_Pan_gridded_to_WV2_MS_centers_new.csv";
//				paths->fname_SRF_for_Spectral_Grouping    = "/data/Projects/JSparseFI/Data/SRF/for_cpp_routine/SRF_WV2_Pan_gridded_to_WV2_MS_centers_new.csv";
//				paths->dir_tmp       = "tmp";
//				//break;
//			}
//			else if(paths->dataSetID_str == "2144111210150"){
//			//#################################################################################################################################################################################################################
//			//#  platform     - orig. sensor    - LR sensor       - HR sensor       - fusion type - filter kernel - scene         - size ID       - fDS         - SNR          - redundant digit (>0 for non-regular datasets #
//			//#  2 (CG-PC)    - 1 (HySpex)      - 4 (WorldView-2) - 4 (WorldView-2) - 1 (MS-Pan)  - 1 (gauss)     - 1 (MUC_Oly)   - 2 (3600x1200) - 10 (fDS=10) - 15 (SNR=15)  - 0                                            #
//			//#################################################################################################################################################################################################################
//				paths->dir_in        = "/data/Projects/JSparseFI/Data/HySpex/MUC/Olymp/l11vnir/gauss/3600x1200";
//				paths->fname_ImZ_ref= paths->dir_in + "/" + "WV2_MSHR"        + "/" + "HySpex_3600x1200_WV2_MSHR.dat";
//				paths->fname_ImZ_init= paths->dir_in + "/" + "WV2_MSHR"        + "/" + "HySpex_3600x1200_WV2_MSHR.dat";
//				paths->fname_ImY     = paths->dir_in + "/" + "WV2_MSLR_fDS10"  + "/" + "HySpex_3600x1200_WV2_MSLR_fDS10_SNR15.dat";
//				paths->fname_ImX     = paths->dir_in + "/" + "WV2_PanHR"       + "/" + "HySpex_3600x1200_WV2_PanHR.dat";
//				paths->dir_out       = "recResults";
//				paths->fname_ImZ     = "2144111210150_rec";
//				paths->fname_SRF                          = "/data/Projects/JSparseFI/Data/SRF/for_cpp_routine/SRF_WV2_Pan_gridded_to_WV2_MS_centers_new.csv";
//				paths->fname_SRF_for_Spectral_Grouping    = "/data/Projects/JSparseFI/Data/SRF/for_cpp_routine/SRF_WV2_Pan_gridded_to_WV2_MS_centers_new.csv";
//				paths->dir_tmp       = "tmp";
//				//break;
//			}
//			else if(paths->dataSetID_str == "2144111210200"){
//			//#################################################################################################################################################################################################################
//			//#  platform     - orig. sensor    - LR sensor       - HR sensor       - fusion type - filter kernel - scene         - size ID       - fDS         - SNR          - redundant digit (>0 for non-regular datasets #
//			//#  2 (CG-PC)    - 1 (HySpex)      - 4 (WorldView-2) - 4 (WorldView-2) - 1 (MS-Pan)  - 1 (gauss)     - 1 (MUC_Oly)   - 2 (3600x1200) - 10 (fDS=10) - 20 (SNR=20)  - 0                                            #
//			//#################################################################################################################################################################################################################
//				paths->dir_in        = "/data/Projects/JSparseFI/Data/HySpex/MUC/Olymp/l11vnir/gauss/3600x1200";
//				paths->fname_ImZ_ref= paths->dir_in + "/" + "WV2_MSHR"        + "/" + "HySpex_3600x1200_WV2_MSHR.dat";
//				paths->fname_ImZ_init= paths->dir_in + "/" + "WV2_MSHR"        + "/" + "HySpex_3600x1200_WV2_MSHR.dat";
//				paths->fname_ImY     = paths->dir_in + "/" + "WV2_MSLR_fDS10"  + "/" + "HySpex_3600x1200_WV2_MSLR_fDS10_SNR20.dat";
//				paths->fname_ImX     = paths->dir_in + "/" + "WV2_PanHR"       + "/" + "HySpex_3600x1200_WV2_PanHR.dat";
//				paths->dir_out       = "recResults";
//				paths->fname_ImZ     = "2144111210200_rec";
//				paths->fname_SRF                          = "/data/Projects/JSparseFI/Data/SRF/for_cpp_routine/SRF_WV2_Pan_gridded_to_WV2_MS_centers_new.csv";
//				paths->fname_SRF_for_Spectral_Grouping    = "/data/Projects/JSparseFI/Data/SRF/for_cpp_routine/SRF_WV2_Pan_gridded_to_WV2_MS_centers_new.csv";
//				paths->dir_tmp       = "tmp";
//				//break;
//			}
//			else if(paths->dataSetID_str == "2144111210300"){
//			//#################################################################################################################################################################################################################
//			//#  platform     - orig. sensor    - LR sensor       - HR sensor       - fusion type - filter kernel - scene         - size ID       - fDS         - SNR          - redundant digit (>0 for non-regular datasets #
//			//#  2 (CG-PC)    - 1 (HySpex)      - 4 (WorldView-2) - 4 (WorldView-2) - 1 (MS-Pan)  - 1 (gauss)     - 1 (MUC_Oly)   - 2 (3600x1200) - 10 (fDS=10) - 30 (SNR=30)  - 0                                            #
//			//#################################################################################################################################################################################################################
//				paths->dir_in        = "/data/Projects/JSparseFI/Data/HySpex/MUC/Olymp/l11vnir/gauss/3600x1200";
//				paths->fname_ImZ_ref= paths->dir_in + "/" + "WV2_MSHR"        + "/" + "HySpex_3600x1200_WV2_MSHR.dat";
//				paths->fname_ImZ_init= paths->dir_in + "/" + "WV2_MSHR"        + "/" + "HySpex_3600x1200_WV2_MSHR.dat";
//				paths->fname_ImY     = paths->dir_in + "/" + "WV2_MSLR_fDS10"  + "/" + "HySpex_3600x1200_WV2_MSLR_fDS10_SNR30.dat";
//				paths->fname_ImX     = paths->dir_in + "/" + "WV2_PanHR"       + "/" + "HySpex_3600x1200_WV2_PanHR.dat";
//				paths->dir_out       = "recResults";
//				paths->fname_ImZ     = "2144111210300_rec";
//				paths->fname_SRF                          = "/data/Projects/JSparseFI/Data/SRF/for_cpp_routine/SRF_WV2_Pan_gridded_to_WV2_MS_centers_new.csv";
//				paths->fname_SRF_for_Spectral_Grouping    = "/data/Projects/JSparseFI/Data/SRF/for_cpp_routine/SRF_WV2_Pan_gridded_to_WV2_MS_centers_new.csv";
//				paths->dir_tmp       = "tmp";
//				//break;
//			}
//			else if(paths->dataSetID_str == "2144111210400"){
//			//#################################################################################################################################################################################################################
//			//#  platform     - orig. sensor    - LR sensor       - HR sensor       - fusion type - filter kernel - scene         - size ID       - fDS         - SNR          - redundant digit (>0 for non-regular datasets #
//			//#  2 (CG-PC)    - 1 (HySpex)      - 4 (WorldView-2) - 4 (WorldView-2) - 1 (MS-Pan)  - 1 (gauss)     - 1 (MUC_Oly)   - 2 (3600x1200) - 10 (fDS=10) - 40 (SNR=40)  - 0                                            #
//			//#################################################################################################################################################################################################################
//				paths->dir_in        = "/data/Projects/JSparseFI/Data/HySpex/MUC/Olymp/l11vnir/gauss/3600x1200";
//				paths->fname_ImZ_ref= paths->dir_in + "/" + "WV2_MSHR"        + "/" + "HySpex_3600x1200_WV2_MSHR.dat";
//				paths->fname_ImZ_init= paths->dir_in + "/" + "WV2_MSHR"        + "/" + "HySpex_3600x1200_WV2_MSHR.dat";
//				paths->fname_ImY     = paths->dir_in + "/" + "WV2_MSLR_fDS10"  + "/" + "HySpex_3600x1200_WV2_MSLR_fDS10_SNR40.dat";
//				paths->fname_ImX     = paths->dir_in + "/" + "WV2_PanHR"       + "/" + "HySpex_3600x1200_WV2_PanHR.dat";
//				paths->dir_out       = "recResults";
//				paths->fname_ImZ     = "2144111210400_rec";
//				paths->fname_SRF                          = "/data/Projects/JSparseFI/Data/SRF/for_cpp_routine/SRF_WV2_Pan_gridded_to_WV2_MS_centers_new.csv";
//				paths->fname_SRF_for_Spectral_Grouping    = "/data/Projects/JSparseFI/Data/SRF/for_cpp_routine/SRF_WV2_Pan_gridded_to_WV2_MS_centers_new.csv";
//				paths->dir_tmp       = "tmp";
//				//break;
//			}
//
//	//------> for reconstruction of band 1, with band 2 previously reconstructed via JSM of bands 2~5
//			else if(paths->dataSetID_str == "2144111210002"){
//			//#################################################################################################################################################################################################################
//			//#  platform     - orig. sensor    - LR sensor       - HR sensor       - fusion type - filter kernel - scene         - size ID       - fDS         - SNR          - redundant digit (>0 for non-regular datasets #
//			//#  2 (CG-PC)    - 1 (HySpex)      - 4 (WorldView-2) - 4 (WorldView-2) - 1 (MS-Pan)  - 1 (gauss)     - 1 (MUC_Oly)   - 2 (3600x1200) - 10 (fDS=10) - 00 (SNR=inf) - 2  (with Pan-band replaced by rec. band 2)   #
//			//#################################################################################################################################################################################################################
//				paths->dir_in        = "/data/Projects/JSparseFI/Data/HySpex/MUC/Olymp/l11vnir/gauss/3600x1200";
//				paths->fname_ImZ_ref= paths->dir_in + "/" + "WV2_MSHR"        + "/" + "HySpex_3600x1200_WV2_MSHR.dat";
//				paths->fname_ImZ_init= paths->dir_in + "/" + "WV2_MSHR"        + "/" + "HySpex_3600x1200_WV2_MSHR.dat";
//				paths->fname_ImY     = paths->dir_in + "/" + "WV2_MSLR_fDS10"  + "/" + "HySpex_3600x1200_WV2_MSLR_fDS10.dat";
//
//				string dir_in_Pan = "/data/Papers_Posters/Papers/141201_XX_CG_JSpFI_TGRS/1408_2nd_subm/data/fusion_results/new_gauss/HySpexWV2_fDS10/J_SparseFI/b1_via_SparseFI_using_b2_from_JSM/b2_input/ID1144111210000_psz5_ovl4_NDP316_lA1/band2_extracted";
//				paths->fname_ImX     = dir_in_Pan + "/" + "PanHR"       + "/" + "1144111210000_rec_band2.dat";
//
//				paths->dir_out       = "recResults";
//				paths->fname_ImZ     = "2144111210002_rec";
//				paths->fname_SRF                          = "/data/Projects/JSparseFI/Data/SRF/for_cpp_routine/SRF_WV2_Pan_gridded_to_WV2_MS_centers_new.csv";
//				paths->fname_SRF_for_Spectral_Grouping    = "/data/Projects/JSparseFI/Data/SRF/for_cpp_routine/SRF_WV2_Pan_gridded_to_WV2_MS_centers_new.csv";
//				paths->dir_tmp       = "tmp";
//				//break;
//			}
//	//------> for reconstruction of bands 7~8, with band 6 previously reconstructed via SparseFI
//			else if(paths->dataSetID_str == "2144111210006"){
//			//#################################################################################################################################################################################################################
//			//#  platform     - orig. sensor    - LR sensor       - HR sensor       - fusion type - filter kernel - scene         - size ID       - fDS         - SNR          - redundant digit (>0 for non-regular datasets #
//			//#  2 (CG-PC)    - 1 (HySpex)      - 4 (WorldView-2) - 4 (WorldView-2) - 1 (MS-Pan)  - 1 (gauss)     - 1 (MUC_Oly)   - 2 (3600x1200) - 10 (fDS=10) - 00 (SNR=inf) - 6  (with Pan-band replaced by rec. band 6)   #
//			//#################################################################################################################################################################################################################
//				paths->dir_in        = "/data/Projects/JSparseFI/Data/HySpex/MUC/Olymp/l11vnir/gauss/3600x1200";
//				paths->fname_ImZ_ref= paths->dir_in + "/" + "WV2_MSHR"        + "/" + "HySpex_3600x1200_WV2_MSHR.dat";
//				paths->fname_ImZ_init= paths->dir_in + "/" + "WV2_MSHR"        + "/" + "HySpex_3600x1200_WV2_MSHR.dat";
//				paths->fname_ImY     = paths->dir_in + "/" + "WV2_MSLR_fDS10"  + "/" + "HySpex_3600x1200_WV2_MSLR_fDS10.dat";
//
//				string dir_in_Pan = "/data/Papers_Posters/Papers/141201_XX_CG_JSpFI_TGRS/1408_2nd_subm/data/fusion_results/new_gauss/HySpexWV2_fDS10/J_SparseFI/b7to8_via_JSM_using_b6_from_SparseFI/b6_input/134231_ID1144111210000_psz5_ovl4_NDP126_lA1/band6_extracted";
//				paths->fname_ImX     = dir_in_Pan + "/" + "PanHR"       + "/" + "1144111210000_rec_band6.dat";
//
//				paths->dir_out       = "recResults";
//				paths->fname_ImZ     = "2144111210006_rec";
//				paths->fname_SRF                          = "/data/Projects/JSparseFI/Data/SRF/for_cpp_routine/SRF_WV2_Pan_gridded_to_WV2_MS_centers_new.csv";
//				paths->fname_SRF_for_Spectral_Grouping    = "/data/Projects/JSparseFI/Data/SRF/for_cpp_routine/SRF_WV2_Pan_gridded_to_WV2_MS_centers_new.csv";
//				paths->dir_tmp       = "tmp";
//				//break;
//			}
//
//	else if(paths->dataSetID_str == "2144111204000"){
//	//#################################################################################################################################################################################################################
//	//#  platform     - orig. sensor    - LR sensor       - HR sensor       - fusion type - filter kernel - scene         - size ID       - fDS         - SNR          - redundant digit (>0 for non-regular datasets #
//	//#  2 (CG-PC)    - 1 (HySpex)      - 4 (WorldView-2) - 4 (WorldView-2) - 1 (MS-Pan)  - 1 (gauss)     - 1 (MUC_Oly)   - 2 (3600x1200) - 04 (fDS=4)  - 00 (SNR=inf) - 0                                            #
//	//#################################################################################################################################################################################################################
//		paths->dir_in        = "/data/Projects/JSparseFI/Data/HySpex/MUC/Olymp/l11vnir/gauss/3600x1200";
//		paths->fname_ImZ_ref= paths->dir_in + "/" + "WV2_MSHR"        + "/" + "HySpex_3600x1200_WV2_MSHR.dat";
//		paths->fname_ImZ_init= paths->dir_in + "/" + "WV2_MSHR"        + "/" + "HySpex_3600x1200_WV2_MSHR.dat";
//		paths->fname_ImY     = paths->dir_in + "/" + "WV2_MSLR_fDS4"   + "/" + "HySpex_3600x1200_WV2_MSLR_fDS4.dat";
//		paths->fname_ImX     = paths->dir_in + "/" + "WV2_PanHR"       + "/" + "HySpex_3600x1200_WV2_PanHR.dat";
//		paths->dir_out       = "recResults";
//		paths->fname_ImZ     = "2144111104000_rec";
//		paths->fname_SRF                          = "/data/Projects/JSparseFI/Data/SRF/for_cpp_routine/SRF_WV2_Pan_gridded_to_WV2_MS_centers_new.csv";
//		paths->fname_SRF_for_Spectral_Grouping    = "/data/Projects/JSparseFI/Data/SRF/for_cpp_routine/SRF_WV2_Pan_gridded_to_WV2_MS_centers_new.csv";
//		paths->dir_tmp       = "tmp";
//		//break;
//	}
//	//------> for reconstruction of band 1, with band 2 previously reconstructed via JSM of bands 2~5
//			else if(paths->dataSetID_str == "2144111204002"){
//			//#################################################################################################################################################################################################################
//			//#  platform     - orig. sensor    - LR sensor       - HR sensor       - fusion type - filter kernel - scene         - size ID       - fDS         - SNR          - redundant digit (>0 for non-regular datasets #
//			//#  2 (CG-PC)    - 1 (HySpex)      - 4 (WorldView-2) - 4 (WorldView-2) - 1 (MS-Pan)  - 1 (gauss)     - 1 (MUC_Oly)   - 2 (3600x1200) - 04 (fDS=4)  - 00 (SNR=inf) - 2  (with Pan-band replaced by rec. band 2)   #
//			//#################################################################################################################################################################################################################
//				paths->dir_in        = "/data/Projects/JSparseFI/Data/HySpex/MUC/Olymp/l11vnir/gauss/3600x1200";
//				paths->fname_ImZ_ref= paths->dir_in + "/" + "WV2_MSHR"        + "/" + "HySpex_3600x1200_WV2_MSHR.dat";
//				paths->fname_ImZ_init= paths->dir_in + "/" + "WV2_MSHR"        + "/" + "HySpex_3600x1200_WV2_MSHR.dat";
//				paths->fname_ImY     = paths->dir_in + "/" + "WV2_MSLR_fDS4"  + "/" + "HySpex_3600x1200_WV2_MSLR_fDS4.dat";
//
//				string dir_in_Pan = "/data/Papers_Posters/Papers/141201_XX_CG_JSpFI_TGRS/1408_2nd_subm/data/fusion_results/new_gauss/HySpexWV2_fDS4/J_SparseFI/b1_via_SparseFI_using_b2_from_JSM/b2_input/ID1144111204000_psz5_ovl4_NDP251_lA1/band2_extracted";
//				paths->fname_ImX     = dir_in_Pan + "/" + "PanHR"       + "/" + "1144111204000_rec_band2.dat";
//
//				paths->dir_out       = "recResults";
//				paths->fname_ImZ     = "2144111204002_rec";
//				paths->fname_SRF                          = "/data/Projects/JSparseFI/Data/SRF/for_cpp_routine/SRF_WV2_Pan_gridded_to_WV2_MS_centers_new.csv";
//				paths->fname_SRF_for_Spectral_Grouping    = "/data/Projects/JSparseFI/Data/SRF/for_cpp_routine/SRF_WV2_Pan_gridded_to_WV2_MS_centers_new.csv";
//				paths->dir_tmp       = "tmp";
//				//break;
//			}
//	//------> for reconstruction of bands 7~8, with band 6 previously reconstructed via SparseFI
//			else if(paths->dataSetID_str == "2144111204006"){
//			//#################################################################################################################################################################################################################
//			//#  platform     - orig. sensor    - LR sensor       - HR sensor       - fusion type - filter kernel - scene         - size ID       - fDS         - SNR          - redundant digit (>0 for non-regular datasets #
//			//#  2 (CG-PC)    - 1 (HySpex)      - 4 (WorldView-2) - 4 (WorldView-2) - 1 (MS-Pan)  - 1 (gauss)     - 1 (MUC_Oly)   - 2 (3600x1200) - 04 (fDS=4)  - 00 (SNR=inf) - 6  (with Pan-band replaced by rec. band 6)   #
//			//#################################################################################################################################################################################################################
//				paths->dir_in        = "/data/Projects/JSparseFI/Data/HySpex/MUC/Olymp/l11vnir/gauss/3600x1200";
//				paths->fname_ImZ_ref= paths->dir_in + "/" + "WV2_MSHR"        + "/" + "HySpex_3600x1200_WV2_MSHR.dat";
//				paths->fname_ImZ_init= paths->dir_in + "/" + "WV2_MSHR"        + "/" + "HySpex_3600x1200_WV2_MSHR.dat";
//				paths->fname_ImY     = paths->dir_in + "/" + "WV2_MSLR_fDS4"  + "/" + "HySpex_3600x1200_WV2_MSLR_fDS4.dat";
//
//				string dir_in_Pan = "/data/Papers_Posters/Papers/141201_XX_CG_JSpFI_TGRS/1408_2nd_subm/data/fusion_results/new_gauss/HySpexWV2_fDS4/J_SparseFI/b7to8_via_JSM_using_b6_from_SparseFI/b6_input/120922_ID1144111204000_psz5_ovl4_NDP126_lA1/band6_extracted";
//				paths->fname_ImX     = dir_in_Pan + "/" + "PanHR"       + "/" + "1144111204000_rec_band6.dat";
//
//				paths->dir_out       = "recResults";
//				paths->fname_ImZ     = "2144111204006_rec";
//				paths->fname_SRF                          = "/data/Projects/JSparseFI/Data/SRF/for_cpp_routine/SRF_WV2_Pan_gridded_to_WV2_MS_centers_new.csv";
//				paths->fname_SRF_for_Spectral_Grouping    = "/data/Projects/JSparseFI/Data/SRF/for_cpp_routine/SRF_WV2_Pan_gridded_to_WV2_MS_centers_new.csv";
//				paths->dir_tmp       = "tmp";
//				//break;
//			}
//
//	else if(paths->dataSetID_str == "2444111104000"){
//	//#################################################################################################################################################################################################################
//	//#  platform     - orig. sensor    - LR sensor       - HR sensor       - fusion type - filter kernel - scene         - size ID       - fDS         - SNR          - redundant digit (>0 for non-regular datasets #
//	//#  2 (CG-PC)    - 4 (WorldView-2) - 4 (WorldView-2) - 4 (WorldView-2) - 1 (MS-Pan)  - 1 (gauss)     - 1 (MUC)       - 1 (960x1320)  - 04 (fDS=4)  - 00 (SNR=inf) - 0                                            #
//	//#################################################################################################################################################################################################################
//		paths->dir_in        = "/data/Projects/JSparseFI/Data/WorldView2/MUC/gauss/960x1320";
//		paths->fname_ImZ_ref= paths->dir_in + "/" + "MSHR"           + "/" + "WV2_MUC_960x1320_MSHR_ref.dat";
//		paths->fname_ImZ_init= paths->dir_in + "/" + "MSHR"           + "/" + "WV2_MUC_960x1320_MSHR_ref.dat";
//		paths->fname_ImY     = paths->dir_in + "/" + "MSLR_fDS4"      + "/" + "WV2_960x1320_MSLR_fDS4.dat";
//		paths->fname_ImX     = paths->dir_in + "/" + "WV2_PanHR"      + "/" + "WV2_960x1320_WV2_PanHR.dat";
//		paths->dir_out       = "recResults";
//		paths->fname_ImZ     = "2444111104000_rec";
//		paths->fname_SRF                          = "/data/Projects/JSparseFI/Data/SRF/for_cpp_routine/SRF_WV2_Pan_gridded_to_WV2_MS_centers_new.csv";
//		paths->fname_SRF_for_Spectral_Grouping    = "/data/Projects/JSparseFI/Data/SRF/for_cpp_routine/SRF_WV2_Pan_gridded_to_WV2_MS_centers_new.csv";
//		paths->dir_tmp       = "tmp";
//		//break;
//	}
//	//------> for reconstruction of band 1, with band 2 previously reconstructed via JSM of bands 2~5
//		else if(paths->dataSetID_str == "2444111104002"){
//		//#################################################################################################################################################################################################################
//		//#  platform     - orig. sensor    - LR sensor       - HR sensor       - fusion type - filter kernel - scene         - size ID       - fDS         - SNR          - redundant digit (>0 for non-regular datasets #
//		//#  2 (CG-PC)    - 4 (WorldView-2) - 4 (WorldView-2) - 4 (WorldView-2) - 1 (MS-Pan)  - 1 (gauss)     - 1 (MUC)       - 1 (960x1320)  - 04 (fDS=4)  - 00 (SNR=inf) - 2 (with Pan-band replaced by rec. band 2)    #
//		//#################################################################################################################################################################################################################
//			paths->dir_in        = "/data/Projects/JSparseFI/Data/WorldView2/MUC/gauss/960x1320";
//			paths->fname_ImZ_ref= paths->dir_in + "/" + "MSHR"           + "/" + "WV2_MUC_960x1320_MSHR_ref.dat";
//			paths->fname_ImZ_init= paths->dir_in + "/" + "MSHR"           + "/" + "WV2_MUC_960x1320_MSHR_ref.dat";
//			paths->fname_ImY     = paths->dir_in + "/" + "MSLR_fDS4"      + "/" + "WV2_960x1320_MSLR_fDS4.dat";
//
//			string dir_in_Pan = "/data/Papers_Posters/Papers/141201_XX_CG_JSpFI_TGRS/1408_2nd_subm/data/fusion_results/new_gauss/WV2WV2_fDS4/J_SparseFI/b1_via_SparseFI_using_b2_from_JSM/b2_input/ID1444111104000_psz5_ovl4_NDP316_lA1/band2_extracted";
//			paths->fname_ImX     = dir_in_Pan + "/" + "PanHR"       + "/" + "1444111104000_rec_band2.dat";
//
//			paths->dir_out       = "recResults";
//			paths->fname_ImZ     = "2444111104002_rec";
//			paths->fname_SRF                          = "/data/Projects/JSparseFI/Data/SRF/for_cpp_routine/SRF_WV2_Pan_gridded_to_WV2_MS_centers_new.csv";
//			paths->fname_SRF_for_Spectral_Grouping    = "/data/Projects/JSparseFI/Data/SRF/for_cpp_routine/SRF_WV2_Pan_gridded_to_WV2_MS_centers_new.csv";
//			paths->dir_tmp       = "tmp";
//			//break;
//		}
//	//------> for reconstruction of band 1, with band 2 previously reconstructed via JSM of bands 2~5
//		else if(paths->dataSetID_str == "2444111104006"){
//		//#################################################################################################################################################################################################################
//		//#  platform     - orig. sensor    - LR sensor       - HR sensor       - fusion type - filter kernel - scene         - size ID       - fDS         - SNR          - redundant digit (>0 for non-regular datasets #
//		//#  2 (CG-PC)    - 4 (WorldView-2) - 4 (WorldView-2) - 4 (WorldView-2) - 1 (MS-Pan)  - 1 (gauss)     - 1 (MUC)       - 1 (960x1320)  - 04 (fDS=4)  - 00 (SNR=inf) - 6 (with Pan-band replaced by rec. band 6)    #
//		//#################################################################################################################################################################################################################
//			paths->dir_in        = "/data/Projects/JSparseFI/Data/WorldView2/MUC/gauss/960x1320";
//			paths->fname_ImZ_ref= paths->dir_in + "/" + "MSHR"           + "/" + "WV2_MUC_960x1320_MSHR_ref.dat";
//			paths->fname_ImZ_init= paths->dir_in + "/" + "MSHR"           + "/" + "WV2_MUC_960x1320_MSHR_ref.dat";
//			paths->fname_ImY     = paths->dir_in + "/" + "MSLR_fDS4"      + "/" + "WV2_960x1320_MSLR_fDS4.dat";
//
//			string dir_in_Pan = "/data/Papers_Posters/Papers/141201_XX_CG_JSpFI_TGRS/1408_2nd_subm/data/fusion_results/new_gauss/WV2WV2_fDS4/J_SparseFI/b7to8_via_JSM_using_b6_from_SparseFI/b6_input/060643_ID1444111104000_psz5_ovl4_NDP100_lA1/band6_extracted";
//			paths->fname_ImX     = dir_in_Pan + "/" + "PanHR"       + "/" + "1444111104000_rec_band6.dat";
//			paths->dir_out       = "recResults";
//			paths->fname_ImZ     = "2444111104006_rec";
//			paths->fname_SRF                          = "/data/Projects/JSparseFI/Data/SRF/for_cpp_routine/SRF_WV2_Pan_gridded_to_WV2_MS_centers_new.csv";
//			paths->fname_SRF_for_Spectral_Grouping    = "/data/Projects/JSparseFI/Data/SRF/for_cpp_routine/SRF_WV2_Pan_gridded_to_WV2_MS_centers_new.csv";
//			paths->dir_tmp       = "tmp";
//			//break;
//		}
//
//
//
//
//
//
//
//
//    ////*****************************************************************************************************////
//    ////*****************************************************************************************************////
//    ////                                                                                                     ////
//    ////                                                                                                     ////
//    ////                                          Platform: SP PC                                            ////
//    ////                                                                                                     ////
//    ////                                                                                                     ////
//    ////*****************************************************************************************************////
//    ////*****************************************************************************************************////
//
//	else if(paths->dataSetID_str == "3114211108000"){
//	//######################################################################################################################################################################################################################
//	//#  platform     - orig. sensor    - LR sensor       - HR sensor       - fusion type - filter kernel - scene         - size ID       - fDS         - SNR          - redundant digit (>0 for non-regular datasets #
//	//#  3 (SP-PC)    - 1 (HySpex)      - 1 (HySpex)      - 4 (WorldView-2) - 2 (HS-MS)   - 1 (gauss)     - 1 (MUC_Oly)   - 1 (240x120)   - 08 (fDS=8)  - 00 (SNR=inf) - 0                                            #
//	//######################################################################################################################################################################################################################
//		paths->dir_in        = "/media/psf/Home/Programming/JSparseFI/JSpFIHM_input_sample";//***local_path***/data/240x120
//		paths->fname_ImZ_ref= paths->dir_in + "/" + "data/240x120" + "/" + "HSHR"          + "/" + "HySpex_2m_240x120_HSHR.dat";
//		paths->fname_ImZ_init= paths->dir_in + "/" + "data/240x120" + "/" + "HSHR"          + "/" + "HySpex_2m_240x120_HSHR.dat";
//		paths->fname_ImY     = paths->dir_in + "/" + "data/240x120" + "/" + "HSLR_fDS8"     + "/" + "HySpex_240x120_HSLR_fDS8.dat";
//		paths->fname_ImX     = paths->dir_in + "/" + "data/240x120" + "/" + "WV2_MSHR"      + "/" + "HySpex_240x120_WV2_MSHR.dat";
//		paths->dir_out       = paths->dir_in + "/" + "recResults";
//		paths->fname_ImZ     = "2114211108000_rec";
//		paths->fname_SRF                         = paths->dir_in + "/" + "SRF" + "/" + "SRF_WV2_MS_gridded_to_HySpex_centers_new.csv";
//		paths->fname_SRF_for_Spectral_Grouping   = paths->dir_in + "/" + "SRF" + "/" + "SRF_WV2_MS_gridded_to_HySpex_centers_new.csv";
//		paths->dir_tmp        = paths->dir_in + "/" + "tmp";
//		//break;
//	}
//	else if(paths->dataSetID_str == "3115211105000"){
//	//######################################################################################################################################################################################################################
//	//#  platform     - orig. sensor    - LR sensor       - HR sensor       - fusion type - filter kernel - scene         - size ID       - fDS         - SNR          - redundant digit (>0 for non-regular datasets #
//	//#  3 (SP-PC)    - 1 (HySpex)      - 1 (HySpex)      - 5 (Quickbird)   - 2 (HS-MS)   - 1 (gauss)     - 1 (MUC_Oly)   - 1 (240x120)   - 05 (fDS=5)  - 00 (SNR=inf) - 0                                            #
//	//######################################################################################################################################################################################################################
//		paths->dir_in        = "/media/psf/Home/Programming/JSparseFI/JSpFIHM_input_sample";//***local_path***/data/240x120
//		paths->fname_ImZ_ref= paths->dir_in + "/" + "data/240x120" + "/" + "HSHR"          + "/" + "HySpex_2m_240x120_HSHR.dat";
//		paths->fname_ImZ_init= paths->dir_in + "/" + "data/240x120" + "/" + "HSHR"          + "/" + "HySpex_2m_240x120_HSHR.dat";
//		paths->fname_ImY     = paths->dir_in + "/" + "data/240x120" + "/" + "HSLR_fDS5"     + "/" + "HySpex_240x120_HSLR_fDS5.dat";
//		paths->fname_ImX     = paths->dir_in + "/" + "data/240x120" + "/" + "QB_MSHR"      + "/" + "HySpex_240x120_QB_MSHR.dat";
//		paths->dir_out       = paths->dir_in + "/" + "recResults";
//		paths->fname_ImZ     = "3115211105000_rec";
//		paths->fname_SRF                         = paths->dir_in + "/" + "SRF" + "/" + "SRF_QB_MS_gridded_to_HySpex_centers_new.csv";
//		paths->fname_SRF_for_Spectral_Grouping   = paths->dir_in + "/" + "SRF" + "/" + "SRF_QB_MS_gridded_to_HySpex_centers_modForSpectralGrouping.csv";
//		paths->dir_tmp       = paths->dir_in + "/" + "tmp";
//		//break;
//	}
//	else if(paths->dataSetID_str == "3144111104000"){
//	//######################################################################################################################################################################################################################
//	//#  platform     - orig. sensor    - LR sensor       - HR sensor       - fusion type - filter kernel - scene         - size ID       - fDS         - SNR          - redundant digit (>0 for non-regular datasets #
//	//#  3 (SP-PC)    - 1 (HySpex)      - 4 (WorldView-2) - 4 (WorldView-2) - 1 (MS-Pan)  - 1 (gauss)     - 1 (MUC_Oly)   - 1 (240x120)   - 04 (fDS=4)  - 00 (SNR=inf) - 0                                            #
//	//######################################################################################################################################################################################################################
//		paths->dir_in        = "/media/psf/Home/Programming/JSparseFI/JSpFIHM_input_sample";//***local_path***/data/240x120
//		paths->fname_ImZ_ref= paths->dir_in + "/" + "data/240x120" "/" + "WV2_MSHR"       + "/" + "HySpex_240x120_WV2_MSHR.dat";
//		paths->fname_ImZ_init= paths->dir_in + "/" + "data/240x120" "/" + "WV2_MSHR"       + "/" + "HySpex_240x120_WV2_MSHR.dat";
//		paths->fname_ImY     = paths->dir_in + "/" + "data/240x120" "/" + "WV2_MSLR_fDS4"  + "/" + "HySpex_240x120_WV2_MSLR_fDS4.dat";
//		paths->fname_ImX     = paths->dir_in + "/" + "data/240x120" "/" + "WV2_PanHR"      + "/" + "HySpex_240x120_WV2_PanHR.dat";
//		paths->dir_out       = paths->dir_in + "/" + "recResults";
//		paths->fname_ImZ     = "3144111104000_rec";
//		paths->fname_SRF                         = paths->dir_in + "/" + "SRF" + "/" + "SRF_WV2_Pan_gridded_to_WV2_MS_centers_new.csv";
//		paths->fname_SRF_for_Spectral_Grouping   = paths->dir_in + "/" + "SRF" + "/" + "SRF_WV2_Pan_gridded_to_WV2_MS_centers_new.csv";
//		paths->dir_tmp       = paths->dir_in + "/" + "tmp";
//		//break;
//	}
//
//
//
//
//
//
//
//
//
//    ////**********************************************************************************************************************************************************************************************************////
//    ////**********************************************************************************************************************************************************************************************************////
//    ////**********************************************************************************************************************************************************************************************************////
//    ////**********************************************************************************************************************************************************************************************************////
//
//
//
//
//	else{
//		if(my_rank==0){
//			cerr << endl
////					<< "ERROR: Undefined dataset: " << paths->dataSetID << "! " << endl << endl
//					<< "ERROR: Undefined dataset: " << paths->dataSetID_str << "! " << endl << endl
//					<< "Troubleshooting: Use an integer number encrypted like the following:" << endl << endl
//					<< "201220: <- [  2    -    0     -       2       -      1      -     2     -    2   -    0    -    0   ]" << endl
//					<< "           [CG-PC  -  HySpex  -  WorldView-2  -  MUC_Olymp  -  960x420  -  fDS6  -  MS-Pan - SNR=inf]" << endl << endl
//					<< "1st digit: ID of Machine the date is stored on:" << endl
//					<< "                1 = SuperMUC" << endl
//					<< "                2 = CG-PC" << endl
//					<< "                3 = SP-PC" << endl
//					<< "                4 = XZ-PC"
//					<< "                ..." << endl
//					<< "2nd digit: Original Sensor ID:" << endl
//					<< "                0 = HySpex" << endl
//					<< "                1 = HyMap" << endl
//					<< "                2 = WorldView2" << endl
//					<< "                ..." << endl
//					<< "3rd digit: (possibly) Synthesized Sensor ID:" << endl
//					<< "                1 = QuickBird" << endl
//					<< "                2 = WordView-2" << endl
//					<< "                ..." << endl
//					<< "4th digit: Sight ID:" << endl
//					<< "                0 = Munich Allianz arena" << endl
//					<< "                1 = Munich Olympic Park" << endl
//					<< "                2 = Kaufbeuern" << endl
//					<< "                ..." << endl
//					<< "5th digit: size ID:" << endl
//					<< "  e.g. for dataset 200..: (200)0(...) = 1500 x 1500" << endl
//					<< "                          (200)1(...) = 1500 x 1500 - 94HS-6MS bands only" << endl
//					<< "                          (200)2(...) =  800 x  800" << endl
//					<< "                          (200)3(...) =  400 x  400" << endl
//					<< "                          (200)4(...) =  360 x  520" << endl
//					<< "                          (200)5(...) =  300 x  450" << endl
//					<< "                          (200)6(...) =   80 x   80" << endl
//					<< "   and for dataset 201..: (201)0(...) = 3600 x 1200" << endl
//					<< "                          (201)1(...) = 960  x  960" << endl
//					<< "                          (201)2(...) = 960  x  420" << endl
//					<< "                          (201)3(...) = 198  x  390" << endl
//					<< "6th digit: down-sampling factor / resolution ratio ID:" << endl
//					<< "                0 = 2"  << endl
//					<< "                1 = 4"  << endl
//					<< "                2 = 6"  << endl
//					<< "                3 = 10" << endl
//					<< "                4 = 15" << endl
//					<< "7th digit: kind of input image pair:" << endl
//					<< "                0 = multispectral-panchromatic"  << endl
//					<< "                1 = hyperspectral-multispectral" << endl
//					<< "8th digit: SNR ID for low resolution image:" << endl
//					<< "                0 = no noise (SNR = inf)"  << endl
//					<< "                1 = SNR/10 (SNR=10)" << endl
//					<< "                2 = SNR/10 (SNR=20)" << endl
//					<< "                3 = SNR/10 (SNR=30)" << endl
//					<< "                4 = SNR/10 (SNR=40)" << endl
//					<< endl;
//			exit(2);
//		}else{
//			exit(2);
//		}
//		//break;
//	}
//	}
//	//########## DO NOT MODIFY BEYOND THIS LINE ##########
//
//	// get current time and broadcast to all processes
//
//	char buf[15]="";
//	if(my_rank==0){
//		time_t t = time(0);
//		struct tm  tstruct;
//		tstruct = *localtime(&t);
//		strftime(buf, sizeof(buf), "%y%m%d_%H%M%S", &tstruct);
//	}
//	MPI_Bcast(&buf, sizeof(buf), MPI_CHAR, 0,MPI_COMM_WORLD);
//	string mystringstr(buf);
//	if(pSetting->store_patches_tmp_on_drive){
//		paths->dir_tmp_patches = paths->dir_tmp + "/patches/" + mystringstr + "_" + dSetting->jobID;
//		if(my_rank==0){
//			string tmpp(paths->dir_tmp);
//			tmpp += "/patches";
//			mkdir(paths->dir_tmp.c_str(), 0777);
//			chmod(paths->dir_tmp.c_str(), S_IRWXU | S_IRWXG | S_IRWXO);
//			mkdir(tmpp.c_str(), 0777);
//			chmod(tmpp.c_str(), S_IRWXU | S_IRWXG | S_IRWXO);
//			mkdir(paths->dir_tmp_patches.c_str(), 0777);
//			chmod(paths->dir_tmp_patches.c_str(), S_IRWXU | S_IRWXG | S_IRWXO);
//		}
//	}
//}
