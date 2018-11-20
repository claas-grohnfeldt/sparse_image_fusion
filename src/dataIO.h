/*
 * dataIO.h
 *
 *  Created on: Apr 25, 2013
 *      Author: Claas Grohnfeldt
 */

#ifndef DATAIO_H_
#define DATAIO_H_

/* -------------------------------------------------------------------- */
/*      Predeclare various classes before pulling in jSparseFI.h, the   */
/*      public declarations.                                            */
/* -------------------------------------------------------------------- */

class SpEODataset;
class SpEORasterBand;
class SpEODictionary;

struct SpEOReport;
struct SpEOFusionSetting;
struct SpEODataIOSetting;
struct SpEOOutputSetting;
struct SpEOSolverSetting;
struct SpEOParallelSetting;
struct SpEOPaths;
struct SpEOGlobalParams;
//struct SpEOGlobalParams;
struct SpEOParam;
struct SpEOAssessmentMetricsStr_Sep;
struct SpEOAssessmentMetricsStr;
struct SpEOAssessmentMetrics;

/* -------------------------------------------------------------------- */
/*      Pull in the public declarations.                                */
/* -------------------------------------------------------------------- */
#include "includes.h"
// ############# NNLS 160630 >>>>>>>>
#include "nnls.h"
// ######################### <<<<<<<<

using namespace Eigen;
using namespace std;

/*! Pixel data types */
typedef enum {
    /*! Unknown or unspecified type */ 	SpEODT_Unknown = 0,
    /*! Eight bit unsigned integer */		SpEODT_Byte = 1,
    /*! Sixteen bit unsigned integer */	SpEODT_UInt16 = 2,
    /*! Sixteen bit signed integer */		SpEODT_Int16 = 3,
    /*! Thirty two bit unsigned integer */	SpEODT_UInt32 = 4,
    /*! Thirty two bit signed integer */	SpEODT_Int32 = 5,
    /*! Thirty two bit floating point */	SpEODT_Float32 = 6,
    /*! Sixty four bit floating point */	SpEODT_Float64 = 7,
    /*! Complex Int16 */			SpEODT_CInt16 = 8,
    /*! Complex Int32 */			SpEODT_CInt32 = 9,
    /*! Complex Float32 */			SpEODT_CFloat32 = 10,
    /*! Complex Float64 */			SpEODT_CFloat64 = 11,
    SpEODT_TypeCount = 12
} SpEODataType;

// Use this typedef for easy sorting of tables using the std::sort() function
typedef std::pair<int, float> argsort_pair;

/*---------------------------------------------------------------------
 *        types for 16 and 32 bits integers, etc...
 *--------------------------------------------------------------------*/
#if UINT_MAX == 65535
	typedef long            SpEOInt32;
	typedef unsigned long   SpEOUInt32;
#else
	typedef int             SpEOInt32;
	typedef unsigned int    SpEOUInt32;
#endif

typedef short           SpEOInt16;
typedef unsigned short  SpEOUInt16;
typedef unsigned char   SpEOByte;

/* -------------------------------------------------------------------- */
/*      Significant constants.                                          */
/* -------------------------------------------------------------------- */

//currently there are two sparse image fusion methods implemented
enum SpEOFusionMethod{
	SparseFI = 0,
	JSparseFI = 1,
	JSparseFIHM = 2,
	LeastSquares = 3,
	GroupedJSparseFI = 4
};

enum SpEOImFlag{
	imFlag_X        = 0,
	imFlag_X_LR     = 1,
	imFlag_Y        = 2,
	imFlag_Z        = 3,
	imFlag_Z_LR     = 4,
	imFlag_Z_ref    = 5,
	imFlag_X_sim    = 6,
	imFlag_X_sim_LR = 7,
	imFlag_Z_init   = 8

};

enum SpEOSparseSolver{
	JPFISTA = 0,
	SPAMS = 1
};

enum SpEODataFormat{
	SpEODouble = 0,
	SpEOFloat = 1
};

enum SpEOResolution{
	LR = 0, /* low resolution */
	HR = 1  /* high resolution */
};



// My Eigen library matrix types
typedef Matrix<double, Dynamic, Dynamic, Dynamic, RowMajor> SpEOMatrix3DD;
typedef Matrix<double, Dynamic, Dynamic, RowMajor> SpEOMatrixD;
typedef Matrix<double, Dynamic, 1> SpEOVectorD;

typedef Matrix<float, Dynamic, Dynamic, Dynamic, RowMajor> SpEOMatrix3DF;
typedef Matrix<float, Dynamic, Dynamic, RowMajor> SpEOMatrixF;
typedef Matrix<float, Dynamic, 1> SpEOVectorF;

typedef Matrix<int, Dynamic, Dynamic, Dynamic, RowMajor> SpEOMatrix3DI;
typedef Matrix<int, Dynamic, Dynamic, RowMajor> SpEOMatrixI;
typedef Matrix<int, Dynamic, 1> SpEOVectorI;
//////////////////////////////////////////////////////////////////////////////

#define LR_ASSESSMENT   	0
#define HR_ASSESSMENT   	1
#define NUM_OF_ASSESSMENT   10 // RMSE, CC, UIQI, ERGAS, DD, SAM, DLambda, AG_orig, AG_rec, PSNR


#define SSTR( x ) static_cast< std::ostringstream & >( ( std::ostringstream() << std::dec << x ) ).str()


class SpEODataset
{
	 //friend void ChangePrivate( Point & );
	 friend void setMetaInfo(SpEODataset *ImZ, SpEODataset *ImY, SpEODataset *ImX, SpEODataIOSetting *dSet, SpEOGlobalParams *glPrms);
	 friend void prepGroupCalcForJSpFI( int firstBandY, int lastBandY, int panBand,
			                            SpEOGlobalParams &glPrmsglPrms, SpEOFusionSetting &fSetting, SpEODataIOSetting &dSetting,
			                            SpEOMatrixD &SRF_orig,SpEOMatrixD &SRF, SpEOMatrixD &filter_coeff, SpEOMatrixD &gauss_filter,
			                            SpEODataset *ImX    , SpEODataset *ImX_LR,     SpEODataset *ImY,     SpEODataset *ImZ,
			                            SpEODataset *ImX_tmp, SpEODataset *ImX_LR_tmp, SpEODataset *ImY_tmp, SpEODataset *ImZ_tmp );
	 friend void cutRelevantInput(SpEODataset *ImX, SpEODataset *ImY, SpEODataIOSetting *dSet);
	 friend void cutRelevantInput(SpEOResolution eRes, SpEOImFlag imFlag, SpEODataset *ImZ, SpEODataIOSetting *dSet, SpEOGlobalParams *glPrms, bool cutOffLRBoundaryPixels);
	 friend void lowPassFilter_and_downSample(SpEODataset *Im,SpEODataset *Im_LR, SpEODataFormat input_format, SpEODataFormat output_format, SpEOGlobalParams &glPrms);
private:
	std::string fname;

protected:
	//const char *byteOrder;
//	const char *imFormat;
	string imFormat;
//	const char *imFormatLong;
	string imFormatLong;
//	const char *projectionRef;
	string projectionRef;

//	char **papszMetadata;

	SpEODataFormat dataFormat; // SpEODouble or SpEOFloat
	
	long sizeU;
	long sizeV;
	short NCh;
	//unsigned short bitDepth; // e.g. bitDepth == 8 if the range if pixel values equals [0,255]

	SpEOResolution eResolution;
	SpEOImFlag imFlag;
//	SpEORasterBand **rasterBands;

public:
	// moved here from protected part, because the post-optimization needs to overwrite the data of ImZ
	SpEORasterBand **rasterBands;

    SpEORasterBand** get_rasterBands(void){return rasterBands;};
	void set_fname(std::string fnm){fname = fnm;};
	void set_fname(const char *fnm){fname = fnm;};
//	void set_metadata(char **metadata){papszMetadata = metadata;};
//	char **get_metadata(void){return papszMetadata;};//
	string get_imFormat(void){return imFormat;};//
	string get_imFormatLong(void){return imFormatLong;};//
	string get_projectionRef(void){return projectionRef;};//
	double *geoTransform; /* geoTransform[0] == top left x
		 *                               geoTransform[1] == w-e pixel resolution
		 *           					 geoTransform[2] == rotation, 0 if image is "north up"
		 *           					 geoTransform[3] == top left y
		 *           					 geoTransform[4] == rotation, 0 if image is "north up"
		 *           					 geoTransform[5] == n-s pixel resolution                 */
	double *get_geoTransform(void){return geoTransform;};
	const char *get_fname(void){return fname.c_str();};//
	void dataRead(std::string fnm, SpEOReport *report);
//	void dataWrite(const char *new_fname,SpEOReport *report);
	void dataWrite(SpEOReport *report, SpEODataIOSetting *dSet, SpEOFusionSetting *fSet, SpEOParallelSetting *pSet, SpEOGlobalParams *glPrms, SpEOPaths *paths, MPI_Comm comm_write, SpEODataset *ImZ_ref_tmp);

	void writeENVIHeader(SpEOReport *report, SpEODataIOSetting *dSet, SpEOFusionSetting *fSetting, SpEOGlobalParams *glPrms, SpEOPaths *paths);
	void dataWriteParMPIIO(SpEOReport *report, SpEODataIOSetting *dSetting, SpEOFusionSetting *fSetting, SpEOParallelSetting *pSetting, SpEOGlobalParams *glPrms, SpEOPaths *paths, MPI_Comm comm_write);
	void read_patch(SpEOPaths *paths, int iL, int &iP, int uP, int &pCnt, int &vP, int pszH, int uPH, int vPH, int fDS, SpEOMatrixF &myReadMat, SpEOMatrixF &patch, SpEOVectorF &pVecHR, SpEOMatrixF *area_sum, SpEODataIOSetting *dSet, int my_rank);
	long get_sizeV(void){return sizeV;};
	long get_sizeU(void){return sizeU;};
	short get_NCh(void){return NCh;};
	SpEOResolution get_eResolution(void){return eResolution;};
	SpEOImFlag get_imFlag(void){return imFlag;};
	void writePatch(SpEOVectorF vec, int psz, int posU, int posV, int band);
	void replBandByCwiseCuotient(SpEOMatrixF *rasterMat, int band);
	void cutNegCoeff(); // experimental, cuts off the negative values -> check model
	void copy_from_dataset( int firstBand, int lastBand, SpEODataset *Im);
	void copyMetaInfoFromDatasets(SpEODataset *sourceImSpatial, SpEODataset *sourceImSpectral, SpEODataset *sourceImFormat, SpEODataFormat dataFormat);
	void fill_band_data(SpEODataset *sourceIm, int firstBand, int lastBand);
        void shift_image_bandwise(double *&Im_shift);
	//void fillArrayFromMat();
//	SpEODataset(std::string fnm, SpEOResolution eRes)
//	{
//		eResolution = eRes;
//		dataRead(fnm);
//	};
	SpEODataset(SpEOResolution eRes, SpEOImFlag myImFlag){
		eResolution = eRes;
		imFlag = myImFlag;
		geoTransform = new double[6];
	};
	~SpEODataset();
};

/* ******************************************************************** */
/*                           SpEORasterBand                             */
/* ******************************************************************** */
class SpEORasterBand
{
private:
	friend class SpEODataset;
	friend void setMetaInfo(SpEODataset *ImZ, SpEODataset *ImY, SpEODataset *ImX, SpEODataIOSetting *dSet, SpEOGlobalParams *glPrms);
	friend void prepGroupCalcForJSpFI( int firstBandY, int lastBandY, int panBand,
										SpEOGlobalParams &glPrmsglPrms, SpEOFusionSetting &fSetting, SpEODataIOSetting &dSetting,
										SpEOMatrixD &SRF_orig,SpEOMatrixD &SRF, SpEOMatrixD &filter_coeff, SpEOMatrixD &gauss_filter,
										SpEODataset *ImX    , SpEODataset *ImX_LR,     SpEODataset *ImY,     SpEODataset *ImZ,
										SpEODataset *ImX_tmp, SpEODataset *ImX_LR_tmp, SpEODataset *ImY_tmp, SpEODataset *ImZ_tmp );
	friend void cutRelevantInput(SpEODataset *ImX, SpEODataset *ImY, SpEODataIOSetting *dSet);
	friend void cutRelevantInput(SpEOResolution eRes, SpEOImFlag imFlag, SpEODataset *Im, SpEODataIOSetting *dSet, SpEOGlobalParams *glPrms, bool cutOffLRBoundaryPixels);
	friend void lowPassFilter_and_downSample(SpEODataset *Im,SpEODataset *Im_LR, SpEODataFormat input_format, SpEODataFormat output_format, SpEOGlobalParams &glPrms);
	  protected:
	  //public:
	    //SpEODataset  *poDS;
	    //SpEODataType bandDataType;
	    int          nBand; /* 1 based */
	    double	     minMaxVal[2];
	    //int          nRasterVSize;
	    //int          nRasterUSize;
	    //SpEODataType eDataType;
	    string dataType;
//	    const char *dataType;
	    string colorInterp;
//	    const char *colorInterp;

	    //float *bandData; // this is where the data of the entire band is stored as 1D array. Maybe later, this will be replaced by using rasterBlocks only.
//	    SpEOMatrixF bandDataMat; // this is also where the data of the entire band is stored as 1D array. Maybe later, this will be replaced by using rasterBlocks only.
	    SpEOMatrixF bandDataMat; // this is also where the data of the entire band is stored as 1D array. Maybe later, this will be replaced by using rasterBlocks only.

public: 
	    // moved here from private part, because the content of ImZ needs to be overwritten after post-optimization
	    SpEOMatrixD bandDataMatD; // only used for ImZ, i.e. for the reconstructed high resolution image
	    //float *get_bandData(void){return bandData;};
	    SpEOMatrixF *get_bandDataMat(void){return &bandDataMat;};
	    SpEOMatrixD *get_bandDataMatD(void){return &bandDataMatD;};
	    SpEODataType get_bandDataType(void){return bandDataType;};

	    SpEODataType bandDataType;
	    int get_nBand(void){return nBand;}; /* 1 based */
	    double *get_minMaxVal(void){return minMaxVal;};
	    string get_dataType(void){return dataType;};
	    string get_colorInterp(void){return colorInterp;};
//	    SpEORasterBand(){};
//		~SpEORasterBand(){
//			//delete [] bandData;
////			delete const_cast<char*> (dataType);
////			delete [] dataType;
//			//delete const_cast<char*> (colorInterp);
////			delete [] dataType;
////			delete [] colorInterp;
//		};

	    ~SpEORasterBand();
};


std::string get_current_time();



/*********** Class declarations *************/
class SpEODictionary
{
public:
	int NCh;
	SpEOMatrixF* dictMat;
	SpEODataset* correspDataset;
//public:
  //SpEODictionary(SpEODataset *corrDS);
  SpEODictionary(SpEODataset *corrDS, SpEOFusionSetting *fSetting); // for pan images
  SpEODictionary(SpEODataset *corrDS){correspDataset = corrDS;}; // for pan images
  SpEODictionary(){}; // for resulting MS HR image
  ~SpEODictionary();
};

class SpEOReport
{
public:
	ofstream file;
	string fileName;
	string curTime;
	double curTimeSec;
	void initialize(SpEOPaths *paths, SpEODataIOSetting *dSetting, SpEOFusionSetting *fSetting, SpEOOutputSetting *oSetting, SpEOSolverSetting *sSetting, SpEOParallelSetting *pSetting, int argc, char **argv);
	void addEvaluation(SpEOAssessmentMetrics  *assessm_res_HR, SpEODataIOSetting *dSetting, SpEOFusionSetting *fSetting, SpEOOutputSetting *oSetting, SpEOGlobalParams *glPrms, int type);
	void addGlobalParams(SpEOGlobalParams *glPrms, SpEOFusionSetting *fSetting);
	void finalize(SpEOGlobalParams *glPrms);
};

struct SpEOFusionSetting
{
public:
	SpEOFusionMethod fMethod;//
	bool two_step_estimation;
	bool nrmlIm;			//
	bool nrmlDicts;			//
	bool substrMean;		//
	bool ImZ_ref_avlbl;	//
	short Nc;				// Number of HS channels per group, where each group of channels has a jointly sparse representation in at least one MS dictionary
	short No;				//
	double tol_SRF;			//
	int patchsize;			// e.g. patchsize = 5 results in patches of size 5 by 5, measured at low-resolution scale
	int winSize;            // size of window around patch: Must have the same sign as patchsize in order to have both centers matched; Used for correlation calculations (int)
	int overlap;			// number of pixels at overlapping regions of patches, measured at low resolution scale
	double lambda;			// regularization parameter trading the weighting of the l_2,1 norm term in the joint sparsity optimization problem
	double lambdaX;			// regularization parameter trading the relative weighting of the high resolution input patch xHR (double)
	double lambdaY;			// regularization parameter trading the relative weighting of the low resolution input patch yLR (double)
	double lambdaX_im;		// regularization parameter trading the relative weighting of the high resolution input image I_X (double)
	double lambdaY_im;		// regularization parameter trading the relative weighting of the low resolution input image I_Y (double)
	//
//	short *rgbChannelNo;	//
	int NDP;				// number of nearest patches = number of atoms in dynamic dictionaries
//	short fDS;				// down-sampling factor = resolution ration between HR and LR input image
//	int maxIter;			// maximum number of iterations. To be set to one, because it is the only valid number at the moment.
//	int pFirst;             // for tests and to write coefficients of specific patches, the user can specify a subset of patches he wants to be reconstructed. ...
//	int pLast;              // ... Therefore, he specifies the limits pFirst and pLast.
	bool evaluate;           // 1: evaluation will be done immetiately after the reconstruction during the same job; 0: evaluation will have to be done in a separate after processing step
	bool evaluate_ImZ_init;  // evaluate the initial image InZ_init
	int dictselect; 		// Select between NNP(0) and Nonlocal(1) dictionary selection
	//  TMP >>>>>>
	//  37: matrix (dictionary) normalization norm: 0: spectral norm, 1: Frobenious norm (bool)
	bool matrixNorm;
	bool addMeanPixelwise;
	//  <<<<<< TMP
	bool LQ_post_opt;       // flag used to decide whether or not the least square post-minimization (of any individual patch) is activated (bool)
	bool LQ_post_opt_im;    // flag used to decide whether or not the least square post-minimization (of the final image) is activated (bool)
	bool useSimulatedImXforDictLearn; // for all spectral groups: flag used to decide if the high resolution monochromatic image which is used
	                                  // for dictionary selection/learning is created as a correlation-based linear combination of the available
	                                  // channels in ImX (if 'true'), or if one of the available channels in ImX (selection according to SRF integrals)
	                                  // is used as it is (if 'false').
	bool useNewMethodForCalculatingZ;
	double lambdaX_ABC;
	double lambdaY_ABC;
	double lambdaZ_ABC;
	double lambdaZ_ABC_in_1st_iter;
	int ImZ_init_type;
	int iterMain;  // number of iterations/alternations between patchwise calculation of ImZ and the post-refinement of ImZ via full image optimization.
	bool doFullImOptWithoutPatRec;

	double Nc_max;  // maximum size of spectral group above groups will be double-checked and perhaps split into subgoups (int)
	double theta;  // minimum cross-correlation within spectral goups (double)

	int set_neg_to_0;   // = 0  => set negative values to zero only at the very end, before writing the final image
						// = 1  => set negative values to zero only after patch reconstruction
						// = 2  => set negative values to zero only after full image optimization
						// = 3  => set negative values to zero both after patch reconstruction and after full image optimization

	bool use_estimated_SRFs; // use estimated SRFs instead of apriori given ones
	bool fullImOptOnSubspace; // 0: do full image optimization on the full image space (not on a subspace); 1: in the full image optimization step, operate on a subspace (optimize for A, where  ImZ = E*A and E is known) (only effenctive if LQ_post_opt_im is set to 1 above)  (bool)
	bool use_init_value_Eq1Unmixing;

	string subspace_transform_type;
        int subspace_dim;
	int ImX_sim_mode;
	bool SNR_normalization;       // calc. coeff. in full image opt. eq. via SNR calc. of ImX and ImY
	bool balance_ImX_term_coef;   // set the coeff. of |R*Z-X| in full image opt. eq. to NChY/NChX [only relevant if SNR_normalization==1]
	bool use_LRnorm_for_dic_normalization; // 
	bool load_DictHR_and_DictLR;
};

struct SpEODataIOSetting
{
public:
	string jobName;
	string jobID;
	int chBundleFirst;//
	int chBundleLast;
	int uLFirst;                        // low resolution vertival pixel coordinate from where the reconstruction should start (default: 0) (int)
	int uLLast;                         // low resolution vertival pixel coordinate from where the reconstruction should end (default: 9999999) (int)
	int vLFirst;                        // low resolution horizontal pixel coordinate from where the reconstruction should start (default: 0) (int)
	int vLLast;                         // low resolution horizontal pixel coordinate from where the reconstruction should end (default: 9999999) (int)
	bool delete_tmp_patch_folders;      // delete tmp patch folders after the fusion process?
	int dir_tmp_patches_additional_num; // number of additional directories that contain previously processed tmp patches needed to reconstruct full image
	bool imageConstructionOnly;         // Construct the image exclusively from previously processed TMP patches. No patch processing will be done, which makes most of the parameters set above inaffective (bool)
	                                    //   TMP patches must be located in the directories given as last arguments in this list
	bool contUnfinishedRec;             // continue unfinished reconstruction
	int platformID;                     // e.g. =1 for SuperMUC or =2 for CG local
        bool saveAsDouble;
};

struct SpEOOutputSetting
{
public:
	short prec;         // precision: length of the strings of the metric/evaluation numbers
	bool saveAlphas;    // save coefficient vectors (alphas)
	int pFirstAlpha;    // first patch coefficient vector to be saved [applicable ONLY IF (saveAlpha==TRUE)]
	int pLastAlpha;     // last patch coefficient vector to be saved [applicable ONLY IF (saveAlpha==TRUE)]
	bool saveDicts;    // save coefficient vectors (alphas)
	int pFirstDict;    // first patch coefficient vector to be saved [applicable ONLY IF (saveAlpha==TRUE)]
	int pLastDict;     // last patch coefficient vector to be saved [applicable ONLY IF (saveAlpha==TRUE)]
	bool writeImageFile; // write fused image in file (1: create file and write resulting image in file; 0: to not write image in file (useful for analyses only)) (bool)
	bool writeImageFileAfterEveryIter; // write all intermediate image fusion resulta (after every iteration) (1: create file and write resulting image in file; 0: to not write image in file (useful for analyses only)) (bool)
};

struct SpEOSolverSetting
{
public:
	SpEOSparseSolver solver;//
	//int maxiter_in; 	// inner loop (8 is fine)
	int maxiter_out; 	// outer loop
	double tol; 		// tolerance
	//double weight;		//
	// for Least squares post processing on the patch level
	int maxiter_CGLS;   // maximum number of iterations in the GS step to solve the least squares problem
	double tol_r_CGLS;  // error tolerance
	int fix_Alpha;      // decides whether or not the coefficients get updates via least squares: 0: not fixed, 1: fixed, 2: only the first row is fixed.
	bool fix_delta_m;   // decides whether or not the mean values of Z are set to the same mean values as Y (i.e. either delta_m remains the initial zero vector or it gets updated via least squares)
	// for Least squares post processing on the final image level
	int maxiter_CGLS_im; // maximum number of iterations in the GS step to solve the least squares problem on the final image level (int)
	double tol_r_CGLS_im; // error tolerance
};

struct SpEOParallelSetting
{
public:
	int numProcTot;        // total number of (MPI) processes
	int numProcGrp;        // number of (MPI) processes in each group, i.e. number of processes jointly working on one patch
	int numProcPerPatch;
	int workStealingTurns; // switch between two patch reconstruction parallelization strategies.
//	bool parallelWriting;  // data will be written to the resulting file in parallel using MPI IO. This requires almost no memory!!! 0: only root process writes. This version requires twice as much memory an number of the sice
	bool store_patches_tmp_on_drive; // store the reconstructed patches temporarily on hard drive to avoid having the entire resonstructed image in the memory. (bool)
	int parWrNumProc;      // maximum number of processors to be used for writing final image in parallel. [only affective if parWr==1] (int)
};

struct SpEOPaths
{
public:
//	long long int dataSetID;
	std::string dataSetID_str;
	std::string dir_in;
	std::string dir_out;
	std::string fname_ImX;
	std::string fname_ImX_shifted;
//	std::string fname_ImX_LR;
	std::string fname_ImY;
	std::string fname_ImZ_ref;
	std::string fname_ImZ_init;
//	std::string fname_ImZ_init_ref;
	std::string fname_ImZ_init_rec;
	std::string fname_ImZ_init_ImY_US;
	std::string fname_ImZ;
	std::string fname_SRF; // Name of .csv file containing APRIORI GIVEN spectral response functions of multispectral sensor with values rastered to the centers of the spectral response function of the hyperspectral sensor
	std::string fname_SRF_estimated; // Name of .csv file containing the ESTIMATED spectral response functions of multispectral sensor with values rastered to the centers of the spectral response function of the hyperspectral sensor
	std::string fname_SRF_for_Spectral_Grouping; // Name of .csv file containing (possibly modified) spectral response functions of multispectral sensor with values rastered to the centers of the spectral response function of the hyperspectral sensor
//	const char *fname_SRF; // Name of .csv file containing spectral response functions of multispectral sensor with values rastered to the centers of the spectral response function of the hyperspectral sensor
	std::string dir_tmp;
	std::string dir_tmp_patches;
	int dir_tmp_patches_additional_num;
	std::string *dir_tmp_patches_additional; // additional directories that contain previously processed tmp patches needed to reconstruct full image
	std::string PathToIncompletePatchSetCSV; // path to .csv file that contains the list of linear patch IDs (iP numbers) which remain to be processed in oder to finish incomplete reconstruction
	                                         // [applicable ONLY IF (contUnfinishedRec==TRUE)] (string)
        std::string fname_SubspaceTransformMat; // e.g. Endmember matrix as VCA output
	std::string fname_DictHR; // pre-calculated high resolution dictionary
	std::string fname_DictLR; // pre-calculated low resolution dictionary
};

struct SpEOGlobalParams
{
public:
	int NChX;
	int NChY;
	int NChY_subspace;
	int NChX_orig;
	int NChY_orig;
	int NChZ;
	int NPU;                    // number of patches in vertical direction
	int NPV;                    // number of patches in horizontal direction
	int NPU_sub;				// number of patches in vertical direction
	int NPV_sub;				// number of patches in horizontal direction
	int uPFirst;				// vertical index of first patch the reconstructed image region
	int uPLast;					// vertical index of last patch the reconstructed image region
	int vPFirst;				// horizontal index of first patch the reconstructed image region
	int vPLast;					// horizontal index of last patch the reconstructed image region
	short fDS;					// down-sampling factor = resolution ration between HR and LR input image
	int NP;                     // total number of patches (in one band) (= numPatchesY * numPatchesX)
//	int usedPatches;			//
	SpEOVectorI idxPUH;         //
	SpEOVectorI idxPVH;         //
	int* idxChY;
	SpEOMatrixD ***SparseCoeffs;// originally: SpEOMatrixD *SparseCoeffs;
	SpEOMatrixD decMat_C;		// [Only relevant for J-SparseFI-HM] Decision matrix. depends on (1) spectral response function, (2) Nc, and (3) No
	SpEOVectorD decMat_C_IT_ChX;   // Index Table (ChX) for non-zero entries of decMat_C
	SpEOVectorD decMat_C_IT_Bundle;// Index Table (Bundle) for non-zero entries of decMat_C
	int Ng;						// [Only relevant for J-SparseFI-HM] Number of groups of hyperspectral channels that are jointly sharpened
	int *Nc_vec;     			// [Only relevant for J-SparseFI-HM] Vector of length Ng; Numbers of HS channels in channel groups (in case of Spectral Grouping,
	                            // all values in Nc_vec are identical and equal to the Nc value specified by the user in fSettings)
	int sum_Nc_vec;     		// sum of all values stored in Nc_vec
	SpEOMatrixI* P_lmd_idx_row; // for every row (each corresponding to one HS channel iChY) these matrices contain the relevant corresponding block indexes (first row) and the corresponding entry indexes relative to the corresponding block's origin each starting at 0.
	int **P_lmd_idx_bl;          // For every block (1,...,Ng) this array contains the absolute (total) row- and column index in the actual sparse matrix P_lmd. Row indexes are stored in the first row and the corresponding column indexes are stored in the second row.
	SpEOVectorD* P_lmd_vecs;    // vectors containing the non-trivial entries along the diagonal in the corresponding blocks (= spectral groups) in P_lmd
//	int * Nc_start_pos;         // vector containing the index of the first HS channel in each group
//	                            // (in case of Spectral Grouping, all groups (possibly except for the last group) are equally spaced, i.e. Nc_start_pos is linearly increasing)
//	SpEOMatrixF P_lmd;          // Matrix to conduct (possibly weighted) spectral averaging of the hyperspectral spectral bands that correspond to multiple spectral groups
//	SpEOVectorI *P_lmd_idx;     // NChZ vectors of indexes which indicate the non-trivial column entries in row j=1,...,NChZ in matrix P_lmd
	int *Nm;					/* [Only relevant for J-SparseFI-HM] Number of non-zero entries of C in row m=1,...,NChX (NChX = number of bands of X)
								 *   Nm[1] , ... , Nm[NChX] */
	int *Nm2;					/* [Only relevant for J-SparseFI-HM] Number of non-zero entries of C in bundle b=1,...,Ng (Ng = number of groups of HS bands)
								 *   Nm2[1] , ... , Nm2[Ng] */
	int **km;					/* [Only relevant for J-SparseFI-HM] Column indices of non-zero entries of C in row m=1,...,NChX (NChX = number of bands of X)
								 *   km[1][1]  , ... , km[1][Nm[1]]
								 *      ...    , ... ,    ...
								 *   km[NChX][1] , ... , km[NChX][Nm[NChX]] */
	int **km2;					/* [Only relevant for J-SparseFI-HM] Column indices of non-zero entries of C in row m=1,...,NChX (NChX = number of bands of X)
								 *   km2[1][1]  , ... , km2[1][Nm2[1]]
								 *      ...     , ... ,    ...
								 *   km2[Ng][1] , ... , km2[Ng][Nm2[Ng]] */
	double timeMainLoop;		//
	double timeFileWrite;
	double timeTotal;
	double timeDictSelect;
	double timeDictSelect_avg;
	// int numProcPerPatch;     // moved to SpEOParallelSetting
	int numProbPerPatch;
#ifndef _OPENMP
	int *myChX;					// varies for different processes
	int *myBundle;
	int myNumProbPerPatch;
#endif
	int numPatchGroups;
	int NP_sub;
	int sizeUL;
	int sizeVL;
	int sizeUH;
	int sizeVH;
	int sizeUL_red;
	int sizeVL_red;
	int sizeUH_red;
	int sizeVH_red;
	string myProcName;          // name of the processor/node the current process my_rank belongs to

//	int myJobInPatch;

	// more to come..
};


struct SpEOAssessmentMetrics
       {
       public:
          double RMSE_mean;
          double *RMSE_sep;
          double PSNR_mean;
          double *PSNR_sep;
          double CC_mean;
          double *CC_sep;
          double ERGAS_mean;
          double *ERGAS_sep;
          double UIQI_mean;
          double *UIQI_sep;
          double DD_mean;
          double *DD_sep;
          double SAM;
          double DLambda_mean;
          double **DLambda_mat;
          double AG_orig_mean;
          double *AG_orig_sep;
          double AG_rec_mean;
          double *AG_rec_sep;
       };

struct SpEOAssessmentMetricsStr
       {
       public:
		string RMSE_str;
		string PSNR_str;
		string CC_str;
		string ERGAS_str;
		string UIQI_str;
		string DD_str;
		string SAM_str;
		string DLambda_str;
		string AG_orig_str;
		string AG_rec_str;
       };

struct SpEOAssessmentMetricsStr_Sep
       {
       public:
		string RMSE_sep_str;
		string PSNR_sep_str;
		string CC_sep_str;
		string ERGAS_sep_str;
		string UIQI_sep_str;
		string DD_sep_str;
		string SAM_sep_str;
		string DLambda_sep_str;
		string AG_orig_sep_str;
		string AG_rec_sep_str;
       };

struct SpEOParam
    {
    public:
        int NCh; // number of channels
        int prec; // user-defined number precision
//        int maxIter; // number of iterations of the J-SparseFI algorithm
    };



// for CSV parser, needed to read spectral response functions
class CSVRow
{
    public:
		std::string line;
        std::string const& operator[](std::size_t index) const{return m_data[index];}
        std::size_t size() const{return m_data.size();}
        void readNextRow(std::istream& str, char delimiter);
    private:
        std::vector<std::string> m_data;
};


void cutRelevantInput(SpEODataset *ImX, SpEODataset *ImY, SpEODataIOSetting *dSet);
void cutRelevantInput(SpEODataset *ImZ, SpEOImFlag imFlag, SpEODataIOSetting *dSet, SpEOGlobalParams *glPrms, bool cutOffLRBoundaryPixels);

//void read_CSV(char *CSVHead, SpEOMatrixF *Mat, const char *fname_CSV, char delimiter, int skipLns);
int read_CSV(SpEOMatrixD *Mat, const char *fname_CSV, char delimiter, int skipLns);
int read_CSV(SpEOMatrixF *Mat, const char *fname_CSV, char delimiter, int skipLns);
int read_CSV(SpEOMatrixI *Mat, const char *fname_CSV, char delimiter, int skipLns);
//void read_SRF(SpEOGlobalParams *glPrms, SpEOMatrixF *Mat, string fname_CSV, char delimiter, bool normalize_SRF);
void read_SRF(SpEOGlobalParams *glPrms, SpEOMatrixD *Mat, string fname_CSV, char delimiter, bool normalize_SRF);
//void read_SRF(SpEOMatrixF *Mat, const char *fname_CSV, char delimiter);
void write_Mat_to_CSV(SpEOMatrixF *Mat, const char *fname_CSV);
void write_Mat_to_CSV(SpEOMatrixD *Mat, const char *fname_CSV);

void calcDecisionMat(SpEOMatrixD *SRF, SpEODataIOSetting *dSetting, SpEOFusionSetting *fSetting, SpEOParallelSetting *pSet, SpEOGlobalParams *glPrms);

//void createImFromTMPFiles(SpEOPaths *paths, SpEOFusionSetting *fSetting, SpEOGlobalParams *glPrms, SpEODataset *ImZ);

void printOneAsTwoDimArray(int nVSz, int nUSz, float *oneDimArray,int numRows,int numCols, int *rows,int *cols);
void convertOneToTwoDimArrayAndPrint(int nVSz, int nUSz, float *oneDimArray);
string createStr(int number, string Str);
string double2str(double num, int prec);
string int2str(int number);
//void createTable(SpEOAssessmentMetricsStr *assessmStr, SpEOAssessmentMetricsStr_Sep *assessmStrSep, SpEOFusionSetting *fSetting, SpEOGlobalParams *glPrms, int *maxDigits, int type, ofstream& reportFile, string reportFileName);
void createTable(SpEOAssessmentMetricsStr *assessmStr, SpEOAssessmentMetricsStr_Sep *assessmStrSep, SpEODataIOSetting *dSetting, SpEOFusionSetting *fSetting, SpEOGlobalParams *glPrms, int *maxDigits, int type, SpEOReport *report);
int remove_dir(const char *path);
void save_dSetting(SpEOPaths *paths, SpEODataIOSetting *fSetting);
void save_fSetting(SpEOPaths *paths, SpEOFusionSetting *fSetting);
void save_oSetting(SpEOPaths *paths, SpEOOutputSetting *oSetting);
void save_sSetting(SpEOPaths *paths, SpEOSolverSetting *sSetting);
void save_pSetting(SpEOPaths *paths, SpEOParallelSetting *pSetting);
void save_glPrms(SpEOPaths *paths, SpEOGlobalParams *glPrms);
void save_fusion_setup(SpEOPaths *paths, SpEODataIOSetting *dSetting, SpEOFusionSetting *fSetting, SpEOParallelSetting *pSetting, SpEOSolverSetting *sSetting, SpEOOutputSetting *oSetting, SpEOGlobalParams *glPrms);
void save_evalResults(SpEOPaths *paths, SpEOGlobalParams *glPrms, SpEOFusionSetting *fSetting, SpEOAssessmentMetrics *assMetrics_HR, int iterMain, int numIterMain, bool beforeFullImOpt, bool doFullImOptWithoutPatRec, bool init_image_eval, bool final_evaluation);
void sample_error(int error, char *string);

/////////////////////// dictionary selection  //////////////////////
void dictSelectFunc(SpEOMatrixF *patchCompSubsetfloat, SpEODataset *ImX_LR, SpEODataset *ImX, SpEODataset *ImY, SpEOFusionSetting *fSet, SpEOGlobalParams *glPrms, int NP, int iP, int uP, int vP, SpEOVectorI *idxPUH, SpEOVectorI *idxPVH, SpEOVectorI *idxPUL, SpEOVectorI *idxPVL, SpEOMatrixD *SRF, SpEOMatrixF *patchComp, int my_rank, int ipp, int &NDP);
void dictSelectFunc(SpEOMatrixD *patchCompSubsetfloat, SpEODataset *ImX_LR, SpEODataset *ImX, SpEODataset *ImY, SpEOFusionSetting *fSet, SpEOGlobalParams *glPrms, int NP, int iP, int uP, int vP, SpEOVectorI *idxPUH, SpEOVectorI *idxPVH, SpEOVectorI *idxPUL, SpEOVectorI *idxPVL, SpEOMatrixD *SRF, SpEOMatrixD *patchComp, int my_rank, int ipp, int &NDP);
// Use this function for easy sorting of tables using the std::sort() function
bool argsort_comp(const argsort_pair& left, const argsort_pair& right);
////////////////////////////////////////////////////////////////////

void checkInputForJSpFI(SpEOGlobalParams &glPrms, SpEOFusionSetting &fSetting, SpEOParallelSetting &pSetting);

//template<typename Derived>
//bool is_finite(const Eigen::MatrixBase<Derived>& x);
//template<typename Derived2>
//bool is_nan(const Eigen::MatrixBase<Derived2>& x);
bool is_inf_or_nan(double x);
bool is_inf_or_nan(float x);
bool is_inf_or_nan(int x);
bool contains_inf_or_nan(SpEOMatrixD x);
bool contains_inf_or_nan(SpEOMatrixF x);
bool contains_inf_or_nan(SpEOMatrixI x);
bool contains_inf_or_nan(SpEOVectorD x);
bool contains_inf_or_nan(SpEOVectorF x);
bool contains_inf_or_nan(SpEOVectorI x);
bool contains_nan(SpEOMatrixD x);
bool contains_nan(SpEOVectorD x);
void check_for_inf_or_nan(int my_rank, double x,      const char *numDescr, int number, const char *descr);
void check_for_inf_or_nan(int my_rank, float x,       const char *numDescr, int number, const char *descr);
void check_for_inf_or_nan(int my_rank, int x,         const char *numDescr, int number, const char *descr);
void check_for_inf_or_nan(int my_rank, SpEOMatrixD x, const char *numDescr, int number, const char *descr);
void check_for_inf_or_nan(int my_rank, SpEOMatrixF x, const char *numDescr, int number, const char *descr);
void check_for_inf_or_nan(int my_rank, SpEOMatrixI x, const char *numDescr, int number, const char *descr);
void check_for_inf_or_nan(int my_rank, SpEOVectorD x, const char *numDescr, int number, const char *descr);
void check_for_inf_or_nan(int my_rank, SpEOVectorF x, const char *numDescr, int number, const char *descr);
void check_for_inf_or_nan(int my_rank, SpEOVectorI x, const char *numDescr, int number, const char *descr);

double calc_Corr_Coeff(SpEOMatrixD mat1, SpEOMatrixD mat2);

void lowPassFilter_and_downSample(SpEODataset *Im,SpEODataset *Im_LR, SpEODataFormat input_format, SpEODataFormat output_format, SpEOGlobalParams &glPrms);
void calc_ImX_sim(SpEODataset *ImX_sim, SpEODataset *ImX, SpEODataset *ImX_LR, SpEODataset *ImY,
		          SpEOFusionSetting *fSetting, SpEODataIOSetting *dSetting, SpEOGlobalParams *glPrms,
		          SpEOMatrixD *patX_LR, SpEOMatrixD *patY, SpEOVectorI *idxPUL, SpEOVectorI *idxPVL, int uP, int vP,
		          bool localCalculation, int Ng, int *Nc_vec, int *idxChY, int sim_mode, MPI_Comm comm_busy);
void calc_P_matrices(SpEOVectorD* P_lmd_vecs_loc, int **P_lmd_idx_bl_loc, SpEOMatrixI* P_lmd_idx_row_loc, int Ng, int *Nc_vec, int NChZ, int *idxChY, int my_rank);
void CSG_corr_based_spectral_grouping(int &Ng, int *idxChY, int *Nc_vec, SpEOMatrixD &winY, double theta, int Nc_max, int my_rank);
void threshold_CCmat(SpEOMatrixD &CCmat, double theta, double q);
void cleanUp_CCmat(SpEOMatrixD &CCmat, double q);
int calc_Ng(SpEOMatrixD &CCmat, double q);
void calc_idxChY_and_Nc_vec(int *idxChY, int *Nc_vec, int Ng, SpEOMatrixD &CCmat, double q, int channelOffset);
void calc_refined_CCmins(double &CCmin1, double &CCmin2, double &CCmin3, SpEOMatrixD &CCmat_grp, double CCmin, int Nc_max);
void grp_sngl_bnd_grps_with_high_corr_nbrs(int **Nc_vec_grp, int **idxChY_grp, int *Ng_grp, int iG_lrg, SpEOMatrixD &CCmat);


// ############# NNLS 160630 >>>>>>>>
//template <typename MatrixType>
//bool testNNLS(const MatrixType &A,
//                const Matrix<typename MatrixType::Scalar, MatrixType::RowsAtCompileTime, 1> &b,
//                const Matrix<typename MatrixType::Scalar, MatrixType::ColsAtCompileTime, 1> &x,
//                double eps);
//bool case_1();
//bool case_2();
//bool case_3();
//bool case_4();
//bool case_5();
//bool case_6();
int NNLS_testing_CG();
// ########################## <<<<<<<<<<
void transform_SpEODataset_to_2D(SpEODataset *Im ,SpEOMatrixD &Im_2D);
bool estimate_SRFs_and_ImX_shift(SpEOMatrixD &SRF, double *&ImX_shift, SpEOMatrixD &ImY_2D, SpEOMatrixD &ImX_LR_2D, int my_rank);







#endif /* DATAIO_H_ */
