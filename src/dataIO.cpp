/*
 * dataIO.cpp
 *
 *  Created on: Apr 25, 2013
 *      Author: claas Grohnfeldt
 */

// source headers
#include "dataIO.h"

#include "filter.h"

// is there a more elegant way to use ENVIDataset??? ---------->
//#include "envidataset.cpp"
// <------------------------------------------------------------
using namespace std;
using namespace Eigen;

//SpEODictionary::SpEODictionary(SpEODataset *corrDS)
SpEODictionary::SpEODictionary(SpEODataset *corrDS,
		SpEOFusionSetting *fSetting) {
	cout << "Create dictionary corresponding to dataset:" << corrDS->get_fname()
			<< std::endl;
	this->correspDataset = corrDS;
}

void printOneAsTwoDimArray(int nVSz, int nUSz, float *oneDimArray, int numRows,
		int numCols, int *rows, int *cols) {
	/*  e.g.
	 int numRows=6;
	 int numCols=6;
	 int rows[6] = {36,37,38,39,40,41};
	 int cols[6] = {326,327,328,329,330,331}; */

	cout << endl << "Print 1D array as 2D array:" << endl;
	cout << "\tiV \t 0 \t 1 \t 2 \t  ...  \t";
	for (int col = 0; col < numCols; col++)
		cout << " " << cols[col] << " \t";
	cout << "...  \t" << nVSz - 3 << "\t" << nVSz - 2 << "\t" << nVSz - 1
			<< endl;
	cout
			<< "iU      ------------------------------------------------------------------------- "
			<< endl;
	for (int iU = 0; iU <= 3; iU++) {
		cout << iU << " \t|\t";
		for (int iV = 0; iV < 3; iV++) {
			cout << oneDimArray[iU * nVSz + iV] << "\t";
		}
		cout << "  ...  \t";
		for (int iV = cols[0]; iV <= cols[numCols - 1]; iV++) {
			cout << (int) oneDimArray[iU * nVSz + iV] << "\t";
		}
		cout << "  ...  \t";
		for (int iV = nVSz - 3; iV < nVSz; iV++) {
			cout << oneDimArray[iU * nVSz + iV] << "\t";
		}
		cout << "|" << endl;
	}
	cout
			<< "... \t|\t... \t... \t... \t... \t... \t... \t... \t  ... \t... \t... \t... \t|"
			<< endl;
	for (int iU = rows[0]; iU <= rows[numRows - 1]; iU++) {
		cout << iU << " \t|\t";
		for (int iV = 0; iV < 3; iV++) {
			cout << oneDimArray[iU * nVSz + iV] << "\t";
		}
		cout << "  ...  \t";
		for (int iV = cols[0]; iV <= cols[numCols - 1]; iV++) {
			cout << (int) oneDimArray[iU * nVSz + iV] << "\t";
		}
		cout << "  ...  \t";
		for (int iV = nVSz - 3; iV < nVSz; iV++) {
			cout << oneDimArray[iU * nVSz + iV] << "\t";
		}
		cout << "|" << endl;
	}
	cout
			<< "... \t|\t... \t... \t... \t... \t... \t... \t... \t  ... \t... \t... \t... \t|"
			<< endl;
	for (int iU = nUSz - 3; iU < nUSz; iU++) {
		cout << iU << " \t|\t";
		for (int iV = 0; iV < 3; iV++) {
			cout << oneDimArray[iU * nVSz + iV] << "\t";
		}
		cout << "  ...  \t";
		for (int iV = cols[0]; iV <= cols[numCols - 1]; iV++) {
			cout << (int) oneDimArray[iU * nVSz + iV] << "\t";
		}
		cout << "  ...  \t";
		for (int iV = nVSz - 3; iV < nVSz; iV++) {
			cout << oneDimArray[iU * nVSz + iV] << "\t";
		}
		cout << "|" << endl;
	}
	cout
			<< "        ------------------------------------------------------------------------- "
			<< endl;
}

void convertOneToTwoDimArrayAndPrint(int nVSz, int nUSz, float *oneDimArray) {
	// Convert 1D Array to 2D array:
	cout << endl << "Convert 1D array to 2D array ..." << endl;
	float twoDimArray[nUSz][nVSz];
	for (int iU = 0; iU < nUSz; iU++) {
		for (int iV = 0; iV < nVSz; iV++) {
			twoDimArray[iU][iV] = oneDimArray[iU * nVSz + iV];
		}
	}
	// Print 2D array:
	cout << "Print 2D array:" << endl;
	cout << endl << "\tiV \t 0 \t 1 \t 2 \t 3 \t  ...  \t" << nVSz - 3 << "\t"
			<< nVSz - 2 << "\t" << nVSz - 1 << endl;
	cout
			<< "iU      ------------------------------------------------------------------------- "
			<< endl;
	for (int iU = 0; iU < 3; iU++) {
		cout << iU << " \t|\t";
		for (int iV = 0; iV < 4; iV++) {
			cout << twoDimArray[iU][iV] << "\t";
		}
		cout << "  ...  \t";
		for (int iV = nVSz - 3; iV < nVSz; iV++) {
			cout << twoDimArray[iU][iV] << "\t";
		}
		cout << "|" << endl;
	}
	cout << "... \t|\t... \t... \t... \t... \t  ... \t... \t... \t... \t|"
			<< endl;
	for (int iU = nUSz - 2; iU < nUSz; iU++) {
		cout << iU << " \t|\t";
		for (int iV = 0; iV < 4; iV++) {
			cout << twoDimArray[iU][iV] << "\t";
		}
		cout << "  ...  \t"; // << endl << "   ..." << endl;
		for (int iV = nVSz - 3; iV < nVSz; iV++) {
			cout << twoDimArray[iU][iV] << "\t";
		}
		cout << "|" << endl;
	}
	cout
			<< "        ------------------------------------------------------------------------- "
			<< endl;
}

//void SpEODataset::set_metadata(char **metadata){
//  //this->papszMetadata = new *char[?]; ???????? to be done
//  this->papszMetadata = metadata;
//}
//char **SpEODataset::get_metadata(void){
//  return this->papszMetadata;
//}

void SpEODataset::dataRead(string fnm, SpEOReport *report) {
	//cout << "bp8.1" << endl;
	int my_rank;
	//cout << "bp8.2" << endl;
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
	//cout << "bp8.3" << endl;
//	string fnm(paths->fname_ImX);
	this->set_fname(fnm);
	//cout << "bp8.4" << endl;
	if (my_rank == 0) {
		cout << "read data from file: " << this->get_fname() << ".. " << endl;
		report->file.open(report->fileName.c_str(),
				fstream::in | fstream::out | fstream::app);
		report->file << "*==================================*" << "\n"
					 << "*  Read image: ";
		switch(this->get_imFlag()){
		case imFlag_X:        report->file <<    "ImX                 *" << "\n"; break;
		case imFlag_X_LR:     report->file <<    "ImX_LR              *" << "\n"; break;
		case imFlag_X_sim:    report->file <<    "ImX_sim             *" << "\n"; break;
		case imFlag_X_sim_LR: report->file <<    "ImX_sim_LR          *" << "\n"; break;
		case imFlag_Y:        report->file <<    "ImY                 *" << "\n"; break;
		case imFlag_Z:        report->file <<    "ImZ                 *" << "\n"; break;
		case imFlag_Z_ref:    report->file <<    "ImZ_ref             *" << "\n"; break;
		case imFlag_Z_init:   report->file <<    "ImZ_init            *" << "\n"; break;
		default: report->file << "\n\n ERROR: program terminated unsuccessfully due to try to read an image with an UNKNOWN IMAGE FLAG!" << "\n";
						 cerr << "\n\n ERROR: program terminated unsuccessfully due to try to read an image with an UNKNOWN IMAGE FLAG!" << "\n";
						 exit(2);
						 break;
		}
				report->file << "*==================================*" << "\n"
					 << " - file name: "
					 << this->get_fname() << "\n";
		//cout << "bp8.5" << endl;
		//cout << "read data from file:" << this->get_fname() << ".. ";
		//cout << "bp8.6" << endl;
	}

	//cout << "bp8.7" << endl;
	GDALDataset *poDataset;
	//cout << "bp8.8" << endl;
	GDALAllRegister();
	//cout << "bp8.9" << endl;
	poDataset = (GDALDataset *) GDALOpen(this->get_fname(), GA_ReadOnly);
	//cout << "bp8.10" << endl;

	//=============== NEWLY ADDED 13.01.2014 =============>
	/*
	 //poDataset->getMetadataItem()


	 // for envi: ------------>
	 //GDALDataset  *poDatasetENVI;
	 //	ENVIDataset *poDatasetENVI;
	 //GDALDataset *poDatasetENVI;
	 //	for(int i=0; i<3; i++){
	 //	  cout << "poDatasetENVI->papszHeader[" << i << "] = " << poDatasetENVI->papszHeader[i];
	 //	}
	 //cout << "poDatasetENVI->getMetadataItem(""wavelength_units"") = " << poDatasetENVI->getMetadataItem("wavelength_units");
	 cout << "poDataset->getMetadataItem(""wavelength_units"") = " << poDataset->getMetadataItem("wavelength_units");


	 GDALOpenInfo *my_poOpenInfo = new GDALOpenInfo( this->get_fname(), GA_ReadOnly, NULL);//char **papszSiblingFiles=NULL);
	 //GDALOpenInfo          ( const char *pszFile,     GDALAccess eAccessIn, char **papszSiblingFiles=NULL)
	 //GDALDatasetH my_GDALDatasetH = GDALOpen ( const char *pszFilename, GDALAccess eAccess)
	 // return: A GDALDatasetH handle or NULL on failure. For C++ applications this handle can be cast to a GDALDataset *.

	 GDALDataset  *poDatasetNew = (GDALDataset *) ENVIDataset::Open(my_poOpenInfo);
	 //poDatasetENVI = ENVIDataset::Open( GDALOpenInfo * poOpenInfo )

	 //<-------------------

	 char **papszMetadata;
	 GDALDriver *poDriver;
	 poDriver = GetGDALDriverManager()->GetDriverByName(poDataset->GetDriver()->GetDescription());
	 papszMetadata = poDriver->GetMetadata();
	 const char *myWL_unit = poDriver->GetMetadataItem("wavelength_units");
	 for(int ii=0; ii<7; ii++){
	 cout << "papszMetadata[" << ii << "] = " << papszMetadata[ii] << endl;
	 }
	 set_metadata(poDriver->GetMetadata()); // currently just forwards the address. TBD: forward content ??
	 //cout << "poDriver->GetMetadataItem(""wavelength_units"") = " << myWL_unit << endl;

	 //	char  **papszHeader
	 //	char	**papszBandNames =
	 //	char	**papszWL =
	 //	const char *pszWLUnits = NULL;
	 //	int nWLCount = CSLCount(papszWL);
	 //	->SetMetadataItem("wavelength_units", pszWLUnits);

	 //<============== NEWLY ADDED 13.01.2014 ==============
	 */
	this->sizeU = poDataset->GetRasterYSize();
	//cout << "bp8.11" << endl;
	this->sizeV = poDataset->GetRasterXSize();
	//cout << "bp8.12" << endl;
	this->NCh = poDataset->GetRasterCount();
	//cout << "bp8.13" << endl;
	this->imFormat = poDataset->GetDriver()->GetDescription();
	//cout << "bp8.14" << endl;
	this->imFormatLong = poDataset->GetDriver()->GetMetadataItem(
			GDAL_DMD_LONGNAME);
	//cout << "bp8.15" << endl;
	this->projectionRef = poDataset->GetProjectionRef();
	//cout << "bp8.16" << endl;

	if (my_rank == 0) {
		report->file << " - Information about the (GDAL-)dataset:\n"
				<< "     - Driver:\t" << this->imFormat << "/"
				<< this->imFormatLong << "\n" << "     - Image Size:\t"
				<< this->sizeV << "(VSize) x " << this->sizeU << "(USize) x "
				<< this->NCh << "(RasterCount)\n";
//		if (this->projectionRef != NULL)
		if (!this->projectionRef.empty())
			report->file << "     - Projection:\t`" << this->projectionRef
					<< "\n";
		else
			report->file << "     - No projection!\n";
		if (poDataset->GetGeoTransform(this->geoTransform) == CE_None) {
			report->file << "     - Origin =\t(top left x, top left y) = ("
					<< this->geoTransform[0] << "," << this->geoTransform[3]
					<< ")\n"
					<< "     - Pixel Size =\t(w-e pixel resolution, n-s pixel resolution) = ("
					<< this->geoTransform[1] << "," << this->geoTransform[5]
					<< ")\n";
		}
		report->file << " - Fetch and reading raster band data: \n";
	}
	//cout << "bp8.17" << endl;
	this->rasterBands = new SpEORasterBand*[this->NCh];
	int bGotMin, bGotMax;

	//cout << "bp8.18" << endl;
	for (int iCh = 1; iCh < this->NCh + 1; iCh++) {
		this->rasterBands[iCh - 1] = new SpEORasterBand();
//		SpEORasterBand *mySpEORasterBand;
//		this->rasterBands[iCh - 1] = mySpEORasterBand;

		// ====================================
		// fetching a raster band
		// ====================================
		if (my_rank == 0) {
			report->file << "     - Raster band no. " << iCh << ":  ";
		}

	//cout << "bp8.19" << endl;
		GDALRasterBand *poBand;
	//cout << "bp8.20" << endl;
		poBand = poDataset->GetRasterBand(iCh);
		

		this->rasterBands[iCh - 1]->bandDataType =
				static_cast<SpEODataType>((int) (poBand->GetRasterDataType()));

	//cout << "bp8.21" << endl;
		this->rasterBands[iCh - 1]->nBand = iCh;
		//poBand->GetBlockSize(&(this->rasterBands[iCh - 1]->nBlockVSize),
		//		&(this->rasterBands[iCh - 1]->nBlockUSize));
		this->rasterBands[iCh - 1]->minMaxVal[0] = poBand->GetMinimum(&bGotMin);
		this->rasterBands[iCh - 1]->minMaxVal[1] = poBand->GetMaximum(&bGotMax);
	//cout << "bp8.22" << endl;
		if (!(bGotMin && bGotMax)) {
			GDALComputeRasterMinMax(poBand, TRUE,	this->rasterBands[iCh - 1]->minMaxVal);
//			GDALComputeRasterMinMax((GDALRasterBandH) poBand, TRUE,	this->rasterBands[iCh - 1]->minMaxVal);
		}
	//cout << "bp8.23" << endl;
		this->rasterBands[iCh - 1]->dataType = GDALGetDataTypeName(
				poBand->GetRasterDataType());
		this->rasterBands[iCh - 1]->colorInterp =
				GDALGetColorInterpretationName(
						poBand->GetColorInterpretation());

	//cout << "bp8.24" << endl;
		if (my_rank == 0) {
//			report->file << "Block = "
//					<< this->rasterBands[iCh - 1]->nBlockVSize << "x"
//					<< this->rasterBands[iCh - 1]->nBlockUSize << ",  ";
			report->file << "Type = " << this->rasterBands[iCh - 1]->dataType
					<< ",  ";
			report->file << "ColorInterp = "
					<< this->rasterBands[iCh - 1]->colorInterp << ",  ";
			report->file << "(Min, Max) = ("
					<< this->rasterBands[iCh - 1]->minMaxVal[0] << ", "
					<< this->rasterBands[iCh - 1]->minMaxVal[1] << ")\n";

	//cout << "bp8.25" << endl;
			if (poBand->GetOverviewCount() > 0)
				report->file << "Band has " << poBand->GetOverviewCount()
						<< " overviews.\n";
	//cout << "bp8.26" << endl;
			if (poBand->GetColorTable() != NULL)
				report->file << "Band has a color table with "
						<< poBand->GetColorTable()->GetColorEntryCount()
						<< "entries.\n";
	//cout << "bp8.27" << endl;
		}
		// ====================================
		// reading raster data 
		// ====================================
		/*
		 * There are a few ways to read raster data, but the most common is via the
		 * GDALRasterBand::RasterIO() method. This method will automatically take care
		 * of data type conversion, up/down sampling and windowing. The following code
		 * will read the first scanline of data into a similarly sized buffer, converting
		 * it to floating point as part of the operation.
		 */

		// Eigen Matrix, which makes data handling easier. Perhaps, the band data array won't be needed in the future.

	//cout << "bp8.28" << endl;
		if(this->get_imFlag()==imFlag_Z || this->get_imFlag()==imFlag_Z_init){
	//cout << "bp8.29" << endl;
			double *bandData = (double *) CPLMalloc( sizeof(double) * this->sizeV * this->sizeU);
	//cout << "bp8.30" << endl;
			poBand->RasterIO(GF_Read, 0, 0, this->sizeV, this->sizeU, bandData, this->sizeV, this->sizeU, GDT_Float64, 0, 0);
	//cout << "bp8.31" << endl;
			this->rasterBands[iCh - 1]->bandDataMatD = SpEOMatrixD::Zero(this->sizeU, this->sizeV);
			this->rasterBands[iCh - 1]->bandDataMatD = SpEOMatrixD::Map(bandData, this->sizeU, this->sizeV);
	//cout << "bp8.32" << endl;
			//		this->rasterBands[iCh - 1]->bandDataMat = Map<SpEOMatrixF>(bandData, this->sizeU, this->sizeV);
			//delete this->rasterBands[iCh-1]->bandData;
	//cout << "bp8.33" << endl;
			CPLFree(bandData);
	//cout << "bp8.34" << endl;
			bandData = NULL;
	//cout << "bp8.35" << endl;
//			free(bandData);

			if(contains_inf_or_nan(this->rasterBands[iCh-1]->bandDataMatD)){
				cerr << endl << endl
					 << " EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE" << endl
					 << " E" << endl
					 << " E    NaN or INF   found in     this->rasterBands[iCh-1=" << iCh-1 << "]->bandDataMatD" << endl
					 << " E" << endl
					 << " EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE"
					 << endl << endl << endl;
				exit(2);
			}
	//cout << "bp8.36" << endl;
		}else{
	//cout << "bp8.37" << endl;
			float *bandData = (float *) CPLMalloc( sizeof(float) * this->sizeV * this->sizeU);
	//cout << "bp8.38" << endl;
			poBand->RasterIO(GF_Read, 0, 0, this->sizeV, this->sizeU, bandData, this->sizeV, this->sizeU, GDT_Float32, 0, 0);
	//cout << "bp8.39" << endl;
			this->rasterBands[iCh - 1]->bandDataMat = SpEOMatrixF::Zero(this->sizeU, this->sizeV);
			this->rasterBands[iCh - 1]->bandDataMat = SpEOMatrixF::Map(bandData, this->sizeU, this->sizeV);
	//cout << "bp8.40" << endl;
			//		this->rasterBands[iCh - 1]->bandDataMat = Map<SpEOMatrixF>(bandData, this->sizeU, this->sizeV);
			//delete this->rasterBands[iCh-1]->bandData;

			CPLFree(bandData);
	//cout << "bp8.41" << endl;
			bandData = NULL;
	//cout << "bp8.42" << endl;
//			free(bandData);


			if(contains_inf_or_nan(this->rasterBands[iCh-1]->bandDataMat)){
				cerr << endl << endl
					 << " EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE" << endl
					 << " E" << endl
					 << " E    NaN or INF   found in     this->rasterBands[iCh-1=" << iCh-1 << "]->bandDataMat" << endl
					 << " E" << endl
					 << " EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE"
					 << endl << endl << endl;
				exit(2);
			}
		}
	//cout << "bp8.43" << endl;

		//



//		this->rasterBands[iCh - 1]->bandData = (float *) CPLMalloc(
//				sizeof(float) * this->sizeV * this->sizeU);
//		poBand->RasterIO(GF_Read, 0, 0, this->sizeV, this->sizeU,
//				this->rasterBands[iCh - 1]->bandData, this->sizeV, this->sizeU,
//				GDT_Float32, 0, 0);
//
//		// Eigen Matrix, which makes data handling easier. Perhaps, the band data array won't be needed in the future.
//		this->rasterBands[iCh - 1]->bandDataMat = Map<SpEOMatrixF>(
//				this->rasterBands[iCh - 1]->bandData, this->sizeU, this->sizeV);
//		//delete this->rasterBands[iCh-1]->bandData;
//		CPLFree(this->rasterBands[iCh-1]->bandData);
//		this->rasterBands[iCh-1]->bandData = NULL;
//
//
//
//
//		// newly commented:
		//GDALClose(poBand);
	//cout << "bp8.44" << endl;
	}
	/**********************************************************************************************/
	//cout << "bp8.45" << endl;
	GDALClose(poDataset);
	//cout << "bp8.46" << endl;
	if (my_rank == 0) {
	//cout << "bp8.47" << endl;
		report->file << "\n";
	//cout << "bp8.48" << endl;
		report->file.close();
	//cout << "bp8.49" << endl;
		cout << "done!" << endl;
	//cout << "bp8.50" << endl;
	}
	//cout << "bp8.51" << endl;
}
//#############################################
void SpEODataset::dataWrite(SpEOReport *report, SpEODataIOSetting *dSet, SpEOFusionSetting *fSet, SpEOParallelSetting *pSet, SpEOGlobalParams *glPrms, SpEOPaths *paths, MPI_Comm comm_write, SpEODataset *ImZ) {

	int my_rank; int my_processes;
	MPI_Comm_rank(comm_write, &my_rank);
	MPI_Comm_size(comm_write, &my_processes);

	int u,v,iL, uL, vL, uH, vH, nL, fDS, iChZ;//, uPH, vPH;
	fDS = glPrms->fDS;
	nL = glPrms->sizeUL*glPrms->sizeVL;
	//a = pszL-fSet->overlap;

//	SpEOMatrixF myReadMat, patch;
//	SpEOVectorF pVecHR;
//	SpEOMatrixF *area_sum = new SpEOMatrixF[this->NCh];
	SpEOMatrixD *area_sum = new SpEOMatrixD[this->NCh];

	if(my_rank==0){
		cout << "\n"
				<< "###########################################################" << "\n"
				<< "##    Write image ImZ (from memory) in parallel using..  ##" << "\n"
				<< "##        - MPI I/O                                      ##" << "\n"
				<< "##        - ENVI header file                             ##" << "\n"
				<< "##        - band interleaved by pixel (BIP)              ##" << "\n"
				<< "##        - 16 bit unsigned integer                      ##" << "\n"
				<< "###########################################################"
				<< "\n" << "\n";

		report->file.open(report->fileName.c_str(),
				fstream::in | fstream::out | fstream::app);
		report->file << "\n"
				<< "###########################################################" << "\n"
				<< "##    Write image ImZ (from memory) in parallel using..  ##" << "\n"
				<< "##        - MPI I/O                                      ##" << "\n"
				<< "##        - ENVI header file                             ##" << "\n"
				<< "##        - band interleaved by pixel (BIP)              ##" << "\n"
				<< "##        - 16 bit unsigned integer                      ##" << "\n"
				<< "###########################################################"
				<< "\n" << "\n";
	}
	// ##################### DOUBLE 64bit
	if(dSet->saveAsDouble){
		// MPI IO declarations and initializations
		double *buf_spectrum_WR;
		buf_spectrum_WR = (double *)malloc(this->NCh*sizeof(double));
		MPI_File fh;
		MPI_Offset disp, offset;
		MPI_Datatype etype, ftype, buftype;
		MPI_Status status;
		int result;
		etype = MPI_DOUBLE;
		ftype = MPI_DOUBLE;
		//char *datarep = (char*)"native";
		buftype = MPI_DOUBLE;

		// communicator group comm_write opens the data file for writing only (and creating, if necessary) */
		result = MPI_File_open(comm_write, (char*)paths->fname_ImZ.c_str(), MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &fh);
//		result = MPI_File_open(comm_write, "recResults/151022/1510_testing_123456/093649_ID2115211105000_lXIm1_lYIm1000/2115211105000_rec", MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &fh);

		if(result != MPI_SUCCESS){
			sample_error(result, (char*)"MPI_File_open");
		}
		disp=0;
		result = MPI_File_set_view(fh, disp, etype, ftype, (char*)"native", MPI_INFO_NULL);
		//int MPI_File_set_view(MPI_File fh, MPI_Offset disp, MPI_Datatype etype, MPI_Datatype, filetype, char *datarep, MPI_Info info)!


		if(result != MPI_SUCCESS){
			sample_error(result, (char*)"MPI_File_set_view");
		}
//		if(my_rank < pSet->parWrNumProc){
//			iL = my_rank;
//		while(iL<=fSet->pLast){ // as long as there are LR pixels to work on
		iL = dSet->uLFirst*glPrms->sizeVL+dSet->vLFirst + my_rank;
		nL = dSet->uLLast*glPrms->sizeVL+dSet->vLLast +1;

		while(iL<nL){ // as long as there are LR pixels to work on

			//cout << "["<<my_rank<<"] iL(nL) = " << iL << "("<< nL<< ")" << endl;
			// calculate the coordinates of the patch
			uL = iL / glPrms->sizeVL;
			vL = iL % glPrms->sizeVL;
			uH = uL*fDS;
			vH = vL*fDS;

			int uH_sub = uH - dSet->uLFirst*fDS;
			int vH_sub = vH - dSet->vLFirst*fDS;

			if(uL>=dSet->uLFirst && vL>=dSet->vLFirst && uL<=dSet->uLLast && vL<=dSet->vLLast){

				//if(my_rank==0){cout << "["  << my_rank << "] write LR pixel: uL=" << uL << ", vL=" << vL << ", iL=" << iL << endl;
				// initialization
				for(iChZ=0; iChZ<this->NCh; iChZ++){
					//area_sum[iChZ] = SpEOMatrixF::Zero(fDS, fDS);
					area_sum[iChZ] = ImZ->get_rasterBands()[iChZ]->bandDataMatD.block(uH_sub, vH_sub, fDS, fDS);
				}
				// write area_sum to image file
				if(dSet->uLFirst==0 && dSet->uLLast==glPrms->sizeUL-1 && dSet->vLFirst==0 && dSet->vLLast==glPrms->sizeVL-1){
					for(u=0; u<fDS; u++){// all HR lines/rows in the current area
						for(v=0; v<fDS; v++){// all HR samples/columns in the current area
							// Allocate and initialize a buffer (buf_spectrum_WR) containing this->NCh
							// unsigned shorts, where the unsigned short at location iChZ will set to (unsigned int)area_sum[iChZ](u,v).
							for(iChZ=0;iChZ<this->NCh;iChZ++){
								buf_spectrum_WR[iChZ] = (double)area_sum[iChZ](u,v);
							}
							// specify the offset relative to the displacement position. In compare to displacement, the offset value is given NOT IN BYTES, but IN MULTIPLES OF etype (here etype==MPI_UNSIGNED_SHORT)
							offset = this->NCh * ( (uL*this->sizeV*fDS)+this->sizeV*u + (vL*fDS) + v);

							// Write the buffer (buf_spectrum_WR)

							result = MPI_File_write_at(fh, offset, buf_spectrum_WR, this->NCh, buftype, &status);
							if(result != MPI_SUCCESS){
								sample_error(result, (char*)"MPI_File_write_at");
							}
						}
					}
				}else{
					for(u=0; u<fDS; u++){// all HR lines/rows in the current area
						for(v=0; v<fDS; v++){// all HR samples/columns in the current area
							// Allocate and initialize a buffer (buf_spectrum_WR) containing this->NCh
							// unsigned shorts, where the unsigned short at location iChZ will set to (unsigned int)area_sum[iChZ](u,v).
							for(iChZ=0;iChZ<this->NCh;iChZ++){
								buf_spectrum_WR[iChZ] = (double)area_sum[iChZ](u,v);
							}
							// specify the offset relative to the displacement position. In compare to displacement, the offset value is given NOT IN BYTES, but IN MULTIPLES OF etype (here etype==MPI_UNSIGNED_SHORT)
							//offset = this->NCh * ( (uL*this->sizeV*fDS)+this->sizeV*u + (vL*fDS) + v);
//							offset =   this->NCh*( (uL*fDS+u)*this->sizeV + vL*fDS+v);
							//offset = this->NCh * ( (((uL-dSet->uLFirst)*fDS+u)*this->sizeV) * (dSet->vLLast-dSet->vLFirst+1)*fDS + (vL-dSet->vLFirst)*fDS+v );
							offset = this->NCh * ( (((uL-dSet->uLFirst)*fDS+u)*glPrms->sizeVH_red) + (vL-dSet->vLFirst)*fDS+v );
							// Write the buffer (buf_spectrum_WR)
							result = MPI_File_write_at(fh, offset, buf_spectrum_WR, this->NCh, buftype, &status);
							if(result != MPI_SUCCESS){
								sample_error(result, (char*)"MPI_File_write_at");
							}
						}
					}
				}
				for(iChZ=0; iChZ<this->NCh; iChZ++){
//					area_sum[iChZ] = SpEOMatrixF::Zero(fDS,fDS);
					area_sum[iChZ] = SpEOMatrixD::Zero(fDS,fDS);
				}
			}
			iL += my_processes;
		}
			MPI_Barrier(comm_write);
		result = MPI_File_close(&fh);
		if(result != MPI_SUCCESS){
			sample_error(result, (char*)"MPI_File_close");
		}
		if(my_rank==0){
			report->file.close();
		}
		delete[] area_sum;
		free(buf_spectrum_WR);
	}else{// ############### UNSIGNED INT 16bit
		// MPI IO declarations and initializations
		unsigned short *buf_spectrum_WR;
		buf_spectrum_WR = (unsigned short *)malloc(this->NCh*sizeof(unsigned short));
		MPI_File fh;
		MPI_Offset disp, offset;
		MPI_Datatype etype, ftype, buftype;
		MPI_Status status;
		int result;
		etype = MPI_UNSIGNED_SHORT;
		ftype = MPI_UNSIGNED_SHORT;
		//char *datarep = (char*)"native";
		buftype = MPI_UNSIGNED_SHORT;

		// communicator group comm_write opens the data file for writing only (and creating, if necessary) */
		result = MPI_File_open(comm_write, (char*)paths->fname_ImZ.c_str(), MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &fh);
//		result = MPI_File_open(comm_write, "recResults/151022/1510_testing_123456/093649_ID2115211105000_lXIm1_lYIm1000/2115211105000_rec", MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &fh);

		if(result != MPI_SUCCESS){
			sample_error(result, (char*)"MPI_File_open");
		}
		disp=0;
		result = MPI_File_set_view(fh, disp, etype, ftype, (char*)"native", MPI_INFO_NULL);
		//int MPI_File_set_view(MPI_File fh, MPI_Offset disp, MPI_Datatype etype, MPI_Datatype, filetype, char *datarep, MPI_Info info)!


		if(result != MPI_SUCCESS){
			sample_error(result, (char*)"MPI_File_set_view");
		}
//		if(my_rank < pSet->parWrNumProc){
//			iL = my_rank;
//		while(iL<=fSet->pLast){ // as long as there are LR pixels to work on
		iL = dSet->uLFirst*glPrms->sizeVL+dSet->vLFirst + my_rank;
		nL = dSet->uLLast*glPrms->sizeVL+dSet->vLLast +1;

		while(iL<nL){ // as long as there are LR pixels to work on

			//cout << "["<<my_rank<<"] iL(nL) = " << iL << "("<< nL<< ")" << endl;
			// calculate the coordinates of the patch
			uL = iL / glPrms->sizeVL;
			vL = iL % glPrms->sizeVL;
			uH = uL*fDS;
			vH = vL*fDS;

			int uH_sub = uH - dSet->uLFirst*fDS;
			int vH_sub = vH - dSet->vLFirst*fDS;

			if(uL>=dSet->uLFirst && vL>=dSet->vLFirst && uL<=dSet->uLLast && vL<=dSet->vLLast){

				//if(my_rank==0){cout << "["  << my_rank << "] write LR pixel: uL=" << uL << ", vL=" << vL << ", iL=" << iL << endl;
				// initialization
				for(iChZ=0; iChZ<this->NCh; iChZ++){
					//area_sum[iChZ] = SpEOMatrixF::Zero(fDS, fDS);
					area_sum[iChZ] = ImZ->get_rasterBands()[iChZ]->bandDataMatD.block(uH_sub, vH_sub, fDS, fDS);
				}
				// write area_sum to image file
				if(dSet->uLFirst==0 && dSet->uLLast==glPrms->sizeUL-1 && dSet->vLFirst==0 && dSet->vLLast==glPrms->sizeVL-1){
					for(u=0; u<fDS; u++){// all HR lines/rows in the current area
						for(v=0; v<fDS; v++){// all HR samples/columns in the current area
							// Allocate and initialize a buffer (buf_spectrum_WR) containing this->NCh
							// unsigned shorts, where the unsigned short at location iChZ will set to (unsigned int)area_sum[iChZ](u,v).
							for(iChZ=0;iChZ<this->NCh;iChZ++){
								buf_spectrum_WR[iChZ] = (unsigned short)area_sum[iChZ](u,v);
							}
							// specify the offset relative to the displacement position. In compare to displacement, the offset value is given NOT IN BYTES, but IN MULTIPLES OF etype (here etype==MPI_UNSIGNED_SHORT)
							offset = this->NCh * ( (uL*this->sizeV*fDS)+this->sizeV*u + (vL*fDS) + v);

							// Write the buffer (buf_spectrum_WR)

							result = MPI_File_write_at(fh, offset, buf_spectrum_WR, this->NCh, buftype, &status);
							if(result != MPI_SUCCESS){
								sample_error(result, (char*)"MPI_File_write_at");
							}
						}
					}
				}else{
					for(u=0; u<fDS; u++){// all HR lines/rows in the current area
						for(v=0; v<fDS; v++){// all HR samples/columns in the current area
							// Allocate and initialize a buffer (buf_spectrum_WR) containing this->NCh
							// unsigned shorts, where the unsigned short at location iChZ will set to (unsigned int)area_sum[iChZ](u,v).
							for(iChZ=0;iChZ<this->NCh;iChZ++){
								buf_spectrum_WR[iChZ] = (unsigned short)area_sum[iChZ](u,v);
							}
							// specify the offset relative to the displacement position. In compare to displacement, the offset value is given NOT IN BYTES, but IN MULTIPLES OF etype (here etype==MPI_UNSIGNED_SHORT)
							//offset = this->NCh * ( (uL*this->sizeV*fDS)+this->sizeV*u + (vL*fDS) + v);
//							offset =   this->NCh*( (uL*fDS+u)*this->sizeV + vL*fDS+v);
							//offset = this->NCh * ( (((uL-dSet->uLFirst)*fDS+u)*this->sizeV) * (dSet->vLLast-dSet->vLFirst+1)*fDS + (vL-dSet->vLFirst)*fDS+v );
							offset = this->NCh * ( (((uL-dSet->uLFirst)*fDS+u)*glPrms->sizeVH_red) + (vL-dSet->vLFirst)*fDS+v );
							// Write the buffer (buf_spectrum_WR)
							result = MPI_File_write_at(fh, offset, buf_spectrum_WR, this->NCh, buftype, &status);
							if(result != MPI_SUCCESS){
								sample_error(result, (char*)"MPI_File_write_at");
							}
						}
					}
				}
				for(iChZ=0; iChZ<this->NCh; iChZ++){
//					area_sum[iChZ] = SpEOMatrixF::Zero(fDS,fDS);
					area_sum[iChZ] = SpEOMatrixD::Zero(fDS,fDS);
				}
			}
			iL += my_processes;
		}
			MPI_Barrier(comm_write);
		result = MPI_File_close(&fh);
		if(result != MPI_SUCCESS){
			sample_error(result, (char*)"MPI_File_close");
		}
		if(my_rank==0){
			report->file.close();
		}
		delete[] area_sum;
		free(buf_spectrum_WR);
	}
}
//
//void SpEODataset::read_patch(SpEOPaths *paths, int iL, int &iP, int uP, int &pCnt, int &vP, int pszH, int uPH, int vPH, int fDS, SpEOMatrixF &myReadMat, SpEOMatrixF &patch, SpEOVectorF &pVecHR, SpEOMatrixF *area_sum, SpEODataIOSetting *dSet, int my_rank){
//	// read patch
//	char buf [paths->dir_tmp_patches.length()+42]; // length(dir) + length(/) + length(patch_u0000_v0016_iP000016.csv) + margin
//	sprintf (buf, "%s/%06d/patch_u%04d_v%04d_iP%06d.csv", paths->dir_tmp_patches.c_str(), iP, uP, vP, iP);
//	int stat_CSV_read = read_CSV(&myReadMat, buf, ',', 0);
//	for(int i=0; stat_CSV_read==-1 && i<dSet->dir_tmp_patches_additional_num; i++){
//		char buf_tmp [paths->dir_tmp_patches_additional[i].length()+42]; // length(dir) + length(/) + length(patch_u0000_v0016_iP000016.csv) + margin
//		sprintf (buf_tmp, "%s/%06d/patch_u%04d_v%04d_iP%06d.csv", paths->dir_tmp_patches_additional[i].c_str(), iP, uP, vP, iP);
//		stat_CSV_read = read_CSV(&myReadMat, buf_tmp, ',', 0);
//		cout << "try to open tmp patch file: " << buf_tmp << endl;
////		if(stat_CSV_read==1){
////			cout << "[" << my_rank << "] Problem Resolved! File was found in this directory: '"<< buf_tmp << "'!" << endl;
////		}
//	}
//	if(stat_CSV_read==-1){
//
//		//dSet->
//		//if !(uP*a<dSet->uLFirst || vP*a<dSet->vLFirst){
//		//				uP++;
//		//			}
//		//			while(vP*a<dSet->vLFirst){
//		//				vP++;
//		//			}
//
//			cout << "[" << my_rank << "] ERROR while writing pixel iL=" << iL << ": The .csv file '"<< buf << "' does not exist in any of the specified TMP patch directories!" << endl;
//		//}
//	}else{
//		// cut off negative coefficients
//		myReadMat = (myReadMat.array()>0).select(myReadMat,0);
//		pVecHR = SpEOVectorF::Zero(myReadMat.rows());
//		for(int iChZ=0; iChZ<this->NCh; iChZ++) {
//			  pVecHR = myReadMat.col(iChZ);
//			  patch = SpEOMatrixF::Map(pVecHR.data(), pszH, pszH);
//			  area_sum[iChZ] += patch.block(uPH,vPH,fDS,fDS);
//		}
//		pCnt++;
//	}
//	vP++; iP++;
//}


//void SpEODataset::dataWrite(const char * new_fname, SpEOReport *report) {
//	report->file.open(report->fileName.c_str(),
//			fstream::in | fstream::out | fstream::app);
//	report->file << "\n"
//			<< "###########################################################" << "\n"
//			<< "##                                                       ##" << "\n"
//			<< "##    Write image ImZ to file using one single process   ##" << "\n"
//			<< "##                                                       ##" << "\n"
//			<< "###########################################################"
//			<< "\n" << "\n";
//
//	cout        << "\n"
//				<< "###########################################################" << "\n"
//				<< "##                                                       ##" << "\n"
//				<< "##    Write image ImZ to file using one single process   ##" << "\n"
//				<< "##                                                       ##" << "\n"
//				<< "###########################################################"
//				<< "\n" << "\n";
//
//	this->set_fname(new_fname);
////	char **papszMetadata;
//	GDALDriver *poDriver;
//	poDriver = GetGDALDriverManager()->GetDriverByName(this->get_imFormat());
//	if (poDriver == NULL)
//		exit(1);
////	papszMetadata = poDriver->GetMetadata();
//	//papszMetadata = this->get_metadata();
//	report->file
//			<< " - Check whether data format (Driver) supports GDAL's Create() method for writing data..\n";
//	if (CSLFetchBoolean(poDriver->GetMetadata(), GDAL_DCAP_CREATE, FALSE)) {
//		report->file << "     -> Driver " << this->get_imFormat()
//				<< " supports Create() method.\n";
//	} else {
//		cout << endl << "     -> ERROR: Driver " << this->get_imFormat()
//				<< " DOES NOT SUPPORT Create() method." << endl;
//		report->file << "     -> ERROR: Driver " << this->get_imFormat()
//				<< " DOES NOT SUPPORT Create() method.\n";
//	}
//	if (CSLFetchBoolean(poDriver->GetMetadata(), GDAL_DCAP_CREATECOPY, FALSE)) {
//		report->file << "     -> Driver " << this->get_imFormat()
//				<< " supports CreateCopy() method.\n";
//	} else {
//		report->file << "     -> Driver " << this->get_imFormat()
//				<< " DOES NOT SUPPORT CreateCopy() method.\n";
//	}
//
//	report->file << " - Create GDAL Dataset..\n";
//
//	GDALDataset *poDstDS;
//	char **papszOptions = NULL;
//
////	papszOptions = CSLSetNameValue( papszOptions, "TILED", "YES" );
////	papszOptions = CSLSetNameValue( papszOptions, "COMPRESS", "PACKBITS" );
//
//	// here we assume that all bands have the same data type!
//	poDstDS =
//			poDriver->Create(this->get_fname(), (int) (this->sizeV),
//					(int) (this->sizeU), (int) (this->NCh),
//					static_cast<GDALDataType>((int) (this->rasterBands[0]->bandDataType)),
//					papszOptions);
//
//	report->file
//			<< " - Write appropriate metadata and raster data to created file...\n";
//	poDstDS->SetGeoTransform(this->get_geoTransform());
//
//	OGRSpatialReference oSRS;
//	char *pszSRS_WKT = NULL;
//	oSRS.SetUTM(11, TRUE);
//	oSRS.SetWellKnownGeogCS("NAD27");
//	oSRS.exportToWkt(&pszSRS_WKT);
//	poDstDS->SetProjection(pszSRS_WKT);
//	CPLFree(pszSRS_WKT);
//
//	for (int iCh = 1; iCh <= this->NCh; iCh++) {
//		GDALRasterBand *poBand;
//		poBand = poDstDS->GetRasterBand(iCh);
//
//		report->file << "     - Band no. " << iCh << ":  ";
//		// to be done: implement a more efficient way of converting the data to it's original data type
//		switch (this->rasterBands[iCh - 1]->bandDataType) {
//		case SpEODT_Byte: {
//			report->file << "Convert data to SpEOByte..\n";
//			SpEOByte rasterBandData[(int) (this->sizeV) * (int) (this->sizeU)];
//			for (int i = 0; i < (int) (this->sizeV) * (int) (this->sizeU); i++)
//				rasterBandData[i] =
//						this->get_rasterBands()[iCh - 1]->bandData[i];
//			poBand->RasterIO(GF_Write, 0, 0, (int) (this->sizeV),
//					(int) (this->sizeU), rasterBandData, (int) (this->sizeV),
//					(int) (this->sizeU),
//					static_cast<GDALDataType>((int) (this->rasterBands[iCh - 1]->bandDataType)),
//					0, 0);
//			break;
//		}
//		case SpEODT_UInt16: {
//			report->file << "Convert data to SpEOUInt16.." << "\n";
//			SpEOUInt16 rasterBandData[(int) (this->sizeV) * (int) (this->sizeU)];
//			for (int i = 0; i < (int) (this->sizeV) * (int) (this->sizeU); i++)
//				rasterBandData[i] =
//						this->get_rasterBands()[iCh - 1]->bandData[i];
//			poBand->RasterIO(GF_Write, 0, 0, (int) (this->sizeV),
//					(int) (this->sizeU), rasterBandData, (int) (this->sizeV),
//					(int) (this->sizeU),
//					static_cast<GDALDataType>((int) (this->rasterBands[iCh - 1]->bandDataType)),
//					0, 0);
//			break;
//		}
//		case SpEODT_Int16: {
//			report->file << "Convert data to SpEOInt16.." << "\n";
//			SpEOInt16 rasterBandData[(int) (this->sizeV) * (int) (this->sizeU)];
//			//cout << endl << "" << this->get_rasterBands()[iCh - 1]-> << endl;
//			for (int i = 0; i < (int) (this->sizeV) * (int) (this->sizeU); i++)
//				rasterBandData[i] =
//						this->get_rasterBands()[iCh - 1]->bandData[i];
//			poBand->RasterIO(GF_Write, 0, 0, (int) (this->sizeV),
//					(int) (this->sizeU), rasterBandData, (int) (this->sizeV),
//					(int) (this->sizeU),
//					static_cast<GDALDataType>((int) (this->rasterBands[iCh - 1]->bandDataType)),
//					0, 0);
//			break;
//		}
//		case SpEODT_UInt32: {
//			report->file << "Convert data to SpEOUint32.." << "\n";
//			SpEOUInt32 rasterBandData[(int) (this->sizeV) * (int) (this->sizeU)];
//			for (int i = 0; i < (int) (this->sizeV) * (int) (this->sizeU); i++)
//				rasterBandData[i] =
//						this->get_rasterBands()[iCh - 1]->bandData[i];
//			poBand->RasterIO(GF_Write, 0, 0, (int) (this->sizeV),
//					(int) (this->sizeU), rasterBandData, (int) (this->sizeV),
//					(int) (this->sizeU),
//					static_cast<GDALDataType>((int) (this->rasterBands[iCh - 1]->bandDataType)),
//					0, 0);
//			break;
//		}
//		case SpEODT_Int32: {
//			report->file << "Convert data to SpEOInt32.." << "\n";
//			SpEOInt32 rasterBandData[(int) (this->sizeV) * (int) (this->sizeU)];
//			for (int i = 0; i < (int) (this->sizeV) * (int) (this->sizeU); i++)
//				rasterBandData[i] =
//						this->get_rasterBands()[iCh - 1]->bandData[i];
//			poBand->RasterIO(GF_Write, 0, 0, (int) (this->sizeV),
//					(int) (this->sizeU), rasterBandData, (int) (this->sizeV),
//					(int) (this->sizeU),
//					static_cast<GDALDataType>((int) (this->rasterBands[iCh - 1]->bandDataType)),
//					0, 0);
//			break;
//		}
//		default: {
//			report->file << "\n\n"
//					<< "CONVERSION FAILED BEFORE WRITING DATA TO FILE !" << "\n"
//					<< "\n" << "\n";
//			SpEOByte rasterBandData[(int) (this->sizeV) * (int) (this->sizeU)];
//			for (int i = 0; i < (int) (this->sizeV) * (int) (this->sizeU); i++)
//				rasterBandData[i] =
//						this->get_rasterBands()[iCh - 1]->bandData[i];
//			poBand->RasterIO(GF_Write, 0, 0, (int) (this->sizeV),
//					(int) (this->sizeU), rasterBandData, (int) (this->sizeV),
//					(int) (this->sizeU),
//					static_cast<GDALDataType>((int) (this->rasterBands[iCh - 1]->bandDataType)),
//					0, 0);
//			break;
//		}
//		}
//	}
//	/* Once we're done, close the dataset */
//	GDALClose((GDALDatasetH) poDstDS);
//	CSLDestroy(papszOptions);
//	report->file.close();
//	chmod(this->get_fname(), S_IRWXU | S_IRWXG | S_IRWXO);
//
//}










//void SpEODataset::dataWritePar(const char * new_fname, SpEOReport *report, SpEOFusionSetting *fSet, SpEOGlobalParams *glPrms, SpEOPaths *paths) {
//
//	int my_rank; int my_processes;
//	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
//	MPI_Comm_size(MPI_COMM_WORLD, &my_processes);
//
//	cout << "my_rank==" << my_rank << " and   new_fname = " << new_fname << endl;
//
////	sleep(my_rank*30);
//
//	int iL, uL, vL, uH, vH, nL, pszL, pszH, a, uP, vP, iP, pCnt, fDS, iChY, uPH, vPH;
//	fDS = glPrms->fDS;
//	nL = glPrms->sizeUL*glPrms->sizeVL;
//	pszL = fSet->patchsize;
//	pszH = pszL*fDS;
//	a = pszL-fSet->overlap;
//
//	SpEOMatrixF myReadMat, patch;
//	//SpEOMatrixF zeroMat; // only temporarily!
//	SpEOVectorF pVecHR;
//	SpEOMatrixF *area_sum = new SpEOMatrixF[glPrms->NChY];
//
//	// read and write data in parallel:
//	iL = my_rank;
//	if(my_rank==0){
//		report->file.open(report->fileName.c_str(),
//				fstream::in | fstream::out | fstream::app);
//		report->file << "\n"
//				<< "###########################################################"
//				<< "\n"
//				<< "##                                                       ##"
//				<< "\n"
//				<< "##                     Write image: ImZ                  ##"
//				<< "\n"
//				<< "##                                                       ##"
//				<< "\n"
//				<< "###########################################################"
//				<< "\n" << "\n";
//		//	char **papszMetadata;
//	}
//	this->set_fname(new_fname);
//
//
//	if(my_rank==7){ // only one (in this case the root) process write it's part in the final image,
//		            // because GDAL does not support parallel writing.
//					// tbd: (1) write in file with BIP (Band Interleave by Pixel) format in parallel
//		            //      (2) read and evaluate the image (be careful with memory)
//
//
//				GDALDriver *poDriver;
//				poDriver = GetGDALDriverManager()->GetDriverByName(this->get_imFormat());
//				if (poDriver == NULL)
//					exit(1);
//			//	papszMetadata = poDriver->GetMetadata();
//				//papszMetadata = this->get_metadata();
//				if(my_rank==0){
//					report->file << " - Check whether data format (Driver) supports GDAL's Create() method for writing data..\n";
//				}
//				if (CSLFetchBoolean(poDriver->GetMetadata(), GDAL_DCAP_CREATE, FALSE)) {
//					if(my_rank==0){
//						report->file << "     -> Driver " << this->get_imFormat()
//								<< " supports Create() method.\n";
//					}
//				} else {
//					if(my_rank==0){
//						cout << endl << "     -> ERROR: Driver " << this->get_imFormat()
//								<< " DOES NOT SUPPORT Create() method." << endl;
//						report->file << "     -> ERROR: Driver " << this->get_imFormat()
//								<< " DOES NOT SUPPORT Create() method.\n";
//					}
//				}
//				if (CSLFetchBoolean(poDriver->GetMetadata(), GDAL_DCAP_CREATECOPY, FALSE)) {
//					if(my_rank==0){
//						report->file << "     -> Driver " << this->get_imFormat()
//									<< " supports CreateCopy() method.\n";
//					}
//				} else {
//					if(my_rank==0){
//						report->file << "     -> Driver " << this->get_imFormat()
//									<< " DOES NOT SUPPORT CreateCopy() method.\n";
//					}
//				}
//				if(my_rank==0){
//					report->file << " - Create GDAL Dataset..\n";
//				}
//				GDALDataset *poDstDS;
//				char **papszOptions = NULL;
//
//			//	papszOptions = CSLSetNameValue( papszOptions, "TILED", "YES" );
//			//	papszOptions = CSLSetNameValue( papszOptions, "COMPRESS", "PACKBITS" );
//
//				// here we assume that all bands have the same data type!
//				poDstDS =
//						poDriver->Create(this->get_fname(), (int) (this->sizeV),
//								(int) (this->sizeU), (int) (this->NCh),
//								static_cast<GDALDataType>((int) (this->rasterBands[0]->bandDataType)),
//								papszOptions);
//
//				if(my_rank==0){
//					report->file << " - Write appropriate metadata and raster data to created file...\n";
//				}
//				poDstDS->SetGeoTransform(this->get_geoTransform());
//
//				OGRSpatialReference oSRS;
//				char *pszSRS_WKT = NULL;
//				oSRS.SetUTM(11, TRUE);
//				oSRS.SetWellKnownGeogCS("NAD27");
//				oSRS.exportToWkt(&pszSRS_WKT);
//				poDstDS->SetProjection(pszSRS_WKT);
//				CPLFree(pszSRS_WKT);
//
//				while(iL<nL){//NP) { // as long as there are patches to work on
//					// calculate the coordinates of the patch
//					uL = iL / glPrms->sizeVL;
//					vL = iL % glPrms->sizeVL;
//
//
//
//					uH = uL*fDS;
//					vH = vL*fDS;
//
//					// initialization
//					for(iChY=0; iChY<glPrms->NChY; iChY++){
//						area_sum[iChY] = SpEOMatrixF::Zero(fDS, fDS);
//					}
//					pCnt = 0;
//					uP = ceil(max(0.0, ((double)(uL-pszL+1)))/a);
//					vP = ceil(max(0.0, ((double)(vL-pszL+1)))/a);
//					iP = uP*glPrms->NPV+vP;
//
//					while(uP*a <= uL  &&  uP*a+pszL-1 < glPrms->sizeUL){
//						// calc relative (local) position of area (part of patch) w.r.t. to the patch's corner coordinates (0,0)
//						uPH = uH - uP*a*fDS;
//
//						while(vP*a <= vL  &&  vP*a+pszL-1 < glPrms->sizeVL){
//							vPH = vH - vP*a*fDS; // CRITICAL: TO BE CHECKED!! PERHAPS COUNTING STARTS FROM THE BOTTOM!
//
//							// read patch
//							char buf [paths->dir_tmp_patches.length()+1+30+4]; // length(dir) + length(/) + length(patch_u0000_v0016_iP000016.csv) + margin
//							sprintf (buf, "%s/patch_u%04d_v%04d_iP%06d.csv", paths->dir_tmp_patches.c_str(), uP, vP, iP);
//							read_CSV(&myReadMat, buf, ',', 0);
//
//							// cut off negative coefficients
////							zeroMat = SpEOMatrixF::Zero(pszH*pszH,glPrms->NChY);
////							SpEOMatrixF boolMat = SpEOMatrixF::Zero(pszH*pszH,glPrms->NChY);
//
//							myReadMat = (myReadMat.array()>0).select(myReadMat,0);
//			//				boolMat = (myReadMat.array() >= 0.0 );// zeroMat.array();
//			//				myReadMat.cwiseProduct(boolMat);
//
//							pVecHR = SpEOVectorF::Zero(myReadMat.rows());
//							for(iChY=0; iChY<glPrms->NChY; iChY++) {
//								  pVecHR = myReadMat.col(iChY);
//								  patch = SpEOMatrixF::Map(pVecHR.data(), pszH, pszH);
//								  area_sum[iChY] += patch.block(uPH,vPH,fDS,fDS);
//							}
//							pCnt++; vP++; iP++;
//						}
//
//						if(vL>=glPrms->sizeVL-pszL && vP == glPrms->NPV-1){
//							vPH = vH - (glPrms->sizeVH-pszH);
//
//							// read patch
//							char buf [paths->dir_tmp_patches.length()+1+30+4]; // length(dir) + length(/) + length(patch_u0000_v0016_iP000016.csv) + margin
//							sprintf (buf, "%s/patch_u%04d_v%04d_iP%06d.csv", paths->dir_tmp_patches.c_str(), uP, vP, iP);
//							read_CSV(&myReadMat, buf, ',', 0);
//
//							// cut off negative coefficients
//							myReadMat = (myReadMat.array()>0).select(myReadMat,0);
//			//				SpEOMatrixF boolMat = myReadMat.array()>=SpEOMatrixF::Zero(pszH*pszH,glPrms->NChY);
//			//				myReadMat.cwiseProduct(boolMat);
//
//							pVecHR = SpEOVectorF::Zero(myReadMat.rows());
//							for(iChY=0; iChY<glPrms->NChY; iChY++) {
//								  pVecHR = myReadMat.col(iChY);
//								  patch = SpEOMatrixF::Map(pVecHR.data(), pszH, pszH);
//								  area_sum[iChY] += patch.block(uPH,vPH,fDS,fDS);
//							}
//							pCnt++; vP++; iP++;
//						}
//						uP++;
//					}
//					if(uL>=glPrms->sizeUL-pszL && uP == glPrms->NPU-1){
//						// calc relative (local) position of area (part of patch) w.r.t. to the patch's corner coordinates (0,0)
//						uPH = uH - (glPrms->sizeUH-pszH);
//						while(vP*a <= vL  &&  vP*a+pszL-1 < glPrms->sizeVL){
//							vPH = vH - vP*a*fDS; // CRITICAL: TO BE CHECKED!! PERHAPS COUNTING STARTS FROM THE BOTTOM!
//
//							// read patch
//							char buf [paths->dir_tmp_patches.length()+1+30+4]; // length(dir) + length(/) + length(patch_u0000_v0016_iP000016.csv) + margin
//							sprintf (buf, "%s/patch_u%04d_v%04d_iP%06d.csv", paths->dir_tmp_patches.c_str(), uP, vP, iP);
//							read_CSV(&myReadMat, buf, ',', 0);
//
//							// cut off negative coefficients
//							myReadMat = (myReadMat.array()>0).select(myReadMat,0);
//			//				SpEOMatrixF boolMat = myReadMat.array()>=SpEOMatrixF::Zero(pszH*pszH,glPrms->NChY);
//			//				myReadMat.cwiseProduct(boolMat);
//
//							pVecHR = SpEOVectorF::Zero(myReadMat.rows());
//							for(iChY=0; iChY<glPrms->NChY; iChY++) {
//								  pVecHR = myReadMat.col(iChY);
//								  patch = SpEOMatrixF::Map(pVecHR.data(), pszH, pszH);
//								  area_sum[iChY] += patch.block(uPH,vPH,fDS,fDS);
//							}
//
//							pCnt++; vP++; iP++;
//						}
//						if(vL>=glPrms->sizeVL-pszL && vP == glPrms->NPV-1){
//							vPH = vH - (glPrms->sizeVH-pszH);
//
//							// read patch
//							char buf [paths->dir_tmp_patches.length()+1+30+4]; // length(dir) + length(/) + length(patch_u0000_v0016_iP000016.csv) + margin
//							sprintf (buf, "%s/patch_u%04d_v%04d_iP%06d.csv", paths->dir_tmp_patches.c_str(), uP, vP, iP);
//							read_CSV(&myReadMat, buf, ',', 0);
//
//							// cut off negative coefficients
//							myReadMat = (myReadMat.array()>0).select(myReadMat,0);
//			//				SpEOMatrixF boolMat = myReadMat.array()>=SpEOMatrixF::Zero(pszH*pszH,glPrms->NChY);
//			//				myReadMat.cwiseProduct(boolMat);
//
//							pVecHR = SpEOVectorF::Zero(myReadMat.rows());
//							for(iChY=0; iChY<glPrms->NChY; iChY++) {
//								  pVecHR = myReadMat.col(iChY);
//								  patch = SpEOMatrixF::Map(pVecHR.data(), pszH, pszH);
//								  area_sum[iChY] += patch.block(uPH,vPH,fDS,fDS);
//							}
//
//							pCnt++; vP++; iP++;
//						}
//						uP++;
//					}
//					// write area_sum to image file
//					for(iChY = 0; iChY < glPrms->NChY; iChY++){
//						// average the area of all touching patches
//						area_sum[iChY] /= pCnt;
//
//						GDALRasterBand *poBand;
//						poBand = poDstDS->GetRasterBand(iChY+1); // GDAL counts from 1 instead of 0!
//
//						if(iL==0){
//							report->file << "     - Band no. " << iChY+1 << ":  ";
//						}
//						// to be done: implement a more efficient way of converting the data to it's original data type
//						switch (this->rasterBands[iChY]->bandDataType) {
//						case SpEODT_Byte: {
//							if(iL==0){
//								report->file << "Convert data to SpEOByte..\n";
//							}
//							SpEOByte rasterBandData[fDS*fDS];
//							for (int i = 0; i < fDS*fDS; i++){
//								rasterBandData[i] = area_sum[iChY](i/fDS,i%fDS);//this->get_rasterBands()[iChY]->bandData[i];
//							}
//							poBand->RasterIO(GF_Write, vH, uH, fDS, fDS, rasterBandData, fDS, fDS,
//									static_cast<GDALDataType>((int) (this->rasterBands[iChY]->bandDataType)),
//									0, 0);
//							break;
//						}
//						case SpEODT_UInt16: {
//							if(iL==0){
//								report->file << "Convert data to SpEOUInt16.." << "\n";
//							}
//							SpEOUInt16 rasterBandData[fDS*fDS];
//							for (int i = 0; i < fDS*fDS; i++){
//								rasterBandData[i] = area_sum[iChY](i/fDS,i%fDS);//this->get_rasterBands()[iChY]->bandData[i];
//							}
//							poBand->RasterIO(GF_Write, vH, uH, fDS, fDS, rasterBandData, fDS, fDS,
//									static_cast<GDALDataType>((int) (this->rasterBands[iChY]->bandDataType)),
//									0, 0);
//							break;
//						}
//						case SpEODT_Int16: {
//							if(iL==0){
//								report->file << "Convert data to SpEOInt16.." << "\n";
//							}
//							SpEOInt16 rasterBandData[fDS*fDS];
//							for (int i = 0; i < fDS*fDS; i++){
//								rasterBandData[i] = area_sum[iChY](i/fDS,i%fDS);//this->get_rasterBands()[iChY]->bandData[i];
//							}
//							poBand->RasterIO(GF_Write, vH, uH, fDS, fDS, rasterBandData, fDS, fDS,
//									static_cast<GDALDataType>((int) (this->rasterBands[iChY]->bandDataType)),
//									0, 0);
//							break;
//						}
//						case SpEODT_UInt32: {
//							if(iL==0){
//								report->file << "Convert data to SpEOUint32.." << "\n";
//							}
//							SpEOUInt32 rasterBandData[fDS*fDS];
//							for (int i = 0; i < fDS*fDS; i++){
//								rasterBandData[i] = area_sum[iChY](i/fDS,i%fDS);//this->get_rasterBands()[iChY]->bandData[i];
//							}
//							poBand->RasterIO(GF_Write, vH, uH, fDS, fDS, rasterBandData, fDS, fDS,
//									static_cast<GDALDataType>((int) (this->rasterBands[iChY]->bandDataType)),
//									0, 0);
//							break;
//						}
//						case SpEODT_Int32: {
//							if(iL==0){
//								report->file << "Convert data to SpEOInt32.." << "\n";
//							}
//							SpEOInt32 rasterBandData[fDS*fDS];
//							for (int i = 0; i < fDS*fDS; i++){
//								rasterBandData[i] = area_sum[iChY](i/fDS,i%fDS);//this->get_rasterBands()[iChY]->bandData[i];
//							}
//							poBand->RasterIO(GF_Write, vH, uH, fDS, fDS, rasterBandData, fDS, fDS,
//									static_cast<GDALDataType>((int) (this->rasterBands[iChY]->bandDataType)),
//									0, 0);
//							break;
//						}
//						default: {
//							if(iL==0){
//								report->file << "\n\n"
//											 << "CONVERSION FAILED BEFORE WRITING DATA TO FILE !" << "\n"
//											 << "\n" << "\n";
//							}
//							SpEOByte rasterBandData[fDS*fDS];
//							for (int i = 0; i < fDS*fDS; i++){
//								rasterBandData[i] = area_sum[iChY](i/fDS,i%fDS);//this->get_rasterBands()[iChY]->bandData[i];
//							}
//							poBand->RasterIO(GF_Write, vH, uH, fDS, fDS, rasterBandData, fDS, fDS,
//									static_cast<GDALDataType>((int) (this->rasterBands[iChY]->bandDataType)),
//									0, 0);
//							break;
//						}
//						}
//					}
//
//
//					for(iChY=0; iChY<glPrms->NChY; iChY++){
//						area_sum[iChY] = SpEOMatrixF::Zero(fDS,fDS);
//					}
////					cout << "[" << my_rank << "]: iL=" << iL << ", uL=" << uL << ", vL=" << vL << ", pCnt=" << pCnt << endl;
//					iL += my_processes;
//				}
//				delete[] area_sum;
//				/* Once we're done, close the dataset */
//				GDALClose((GDALDatasetH) poDstDS);
//				CSLDestroy(papszOptions);
//				if(iL==0){
//					report->file.close();
//				}
//
//
//	} // to be removed
//
//
//
//
//}






void SpEODataset::writeENVIHeader(SpEOReport *report, SpEODataIOSetting *dSet, SpEOFusionSetting *fSet, SpEOGlobalParams *glPrms, SpEOPaths *paths) {


	report->file.open(report->fileName.c_str(),
			fstream::in | fstream::out | fstream::app);
	report->file << "###########################################################" << "\n"
			     << "##    Write ENVI header file..                           ##" << "\n"
			     << "###########################################################"
			     << "\n" << "\n";
	cout <<    "###########################################################" << endl
			<< "##    Write ENVI header file..                           ##" << endl
			<< "###########################################################"
			<< endl << endl;


	this->set_fname(paths->fname_ImZ.c_str());

	/////////////////////////////////////////////////////////////////////////////
	////////////////////////        self made       /////////////////////////////
	/////////////////////////////////////////////////////////////////////////////

	string fileName = paths->fname_ImZ + ".hdr";
	ofstream headerfile;
	headerfile.open(fileName.c_str(),	fstream::in | fstream::out | fstream::app);

	headerfile
	<< "ENVI"
	<< "\n" << "description = {"
	<< "\n" << paths->fname_ImZ
	<< "\n" << "}"
	<< "\n" << "samples = " << (dSet->vLLast-dSet->vLFirst+1)*glPrms->fDS // this->sizeV
	<< "\n" << "lines   = " << (dSet->uLLast-dSet->uLFirst+1)*glPrms->fDS // this->sizeU
	<< "\n" << "bands   = " << this->NCh
	<< "\n" << "header offset = 0"
	<< "\n" << "file type = ENVI Standard";
	if(dSet->saveAsDouble){
		headerfile << "\n" << "data type = 5";
	}else{
		headerfile << "\n" << "data type = 12";
	}
	headerfile
	<< "\n" << "interleave = bip"
	<< "\n" << "x start = " << dSet->vLFirst*glPrms->fDS+1
	<< "\n" << "y start = " << dSet->uLFirst*glPrms->fDS+1
	<< "\n" << "byte order = 0"
//		<< "\n" << "map info = {UTM, 1, 1, 689416.5, 5337561.5, 1, 1, 11, North,North America 1927}"
//		<< "\n" << "coordinate system string = {PROJCS[\"NAD_1927_UTM_Zone_11N\",GEOGCS[\"GCS_North_American_1927\",DATUM[\"D_North_American_1927\",SPHEROID[\"Clarke_1866\",6378206.4,294.9786982]],PRIMEM[\"Greenwich\",0],UNIT[\"Degree\",0.017453292519943295]],PROJECTION[\"Transverse_Mercator\"],PARAMETER[\"latitude_of_origin\",0],PARAMETER[\"central_meridian\",-117],PARAMETER[\"scale_factor\",0.9996],PARAMETER[\"false_easting\",500000],PARAMETER[\"false_northing\",0],UNIT[\"Meter\",1]]}"
	<< "\n" << "band names = {";
	for(int iCh=0; iCh<this->NCh-1; iCh++){
		headerfile << "\n" << "Band " << iCh+1+dSet->chBundleFirst << "," ;
	}
	headerfile << "\n" << "Band " << this->NCh+dSet->chBundleFirst << "}";
	headerfile << "\n";


	report->file.close();
	headerfile.close();
	chmod(fileName.c_str(), S_IRWXU | S_IRWXG | S_IRWXO);

		/////////////////////////////////////////////////////////////////////////////


//
//
//
//	//	char **papszMetadata;
//		GDALDriver *poDriver;
//		poDriver = GetGDALDriverManager()->GetDriverByName("ENVI");
////		poDriver = GetGDALDriverManager()->GetDriverByName(this->get_imFormat());
//		if (poDriver == NULL)
//			exit(1);
//	//	papszMetadata = poDriver->GetMetadata();
//		//papszMetadata = this->get_metadata();
////		report->file << " - Check whether data format (Driver) supports GDAL's Create() method for writing data..\n";
////		if (CSLFetchBoolean(poDriver->GetMetadata(), GDAL_DCAP_CREATE, FALSE)) {
////			report->file << "     -> Driver ENVI" //report->file << "     -> Driver " << this->get_imFormat()
////					     << " supports Create() method.\n";
////		} else {
////			cout << endl << "     -> ERROR: Driver ENVI"
////					<< " DOES NOT SUPPORT Create() method." << endl;
////			report->file << "     -> ERROR: Driver ENVI"
////					<< " DOES NOT SUPPORT Create() method.\n";
////		}
////		if (CSLFetchBoolean(poDriver->GetMetadata(), GDAL_DCAP_CREATECOPY, FALSE)) {
////			report->file << "     -> Driver ENVI"
////					<< " supports CreateCopy() method.\n";
////		} else {
////			report->file << "     -> Driver ENVI"
////					<< " DOES NOT SUPPORT CreateCopy() method.\n";
////		}
//
//		report->file << " - Create GDAL Dataset..\n";
//
//		GDALDataset *poDstDS;
//		char **papszOptions = NULL;
//
//		papszOptions = CSLSetNameValue( papszOptions, "INTERLEAVE", "bip" );
////		papszOptions = CSLSetNameValue( papszOptions, "data type", "12" );//GDT_UInt16
//
//	//	papszOptions = CSLSetNameValue( papszOptions, "TILED", "YES" );
//	//	papszOptions = CSLSetNameValue( papszOptions, "COMPRESS", "PACKBITS" );
//
//		// here we assume that all bands have the same data type!
//		poDstDS =
//				poDriver->Create(this->get_fname(), (int) (this->sizeV),
//						(int) (this->sizeU), (int) (this->NCh),
//						static_cast<GDALDataType>((int) 12),
//						papszOptions);
//
//		report->file << " - Write appropriate metadata and raster data to created file...\n";
//		poDstDS->SetGeoTransform(this->get_geoTransform());
//
//		OGRSpatialReference oSRS;
//		char *pszSRS_WKT = NULL;
//		oSRS.SetUTM(11, TRUE);
//		oSRS.SetWellKnownGeogCS("NAD27");
//		oSRS.exportToWkt(&pszSRS_WKT);
//		poDstDS->SetProjection(pszSRS_WKT);
//		CPLFree(pszSRS_WKT);
//
////		for (int iCh = 1; iCh <= this->NCh; iCh++) {
////			GDALRasterBand *poBand;
////			poBand = poDstDS->GetRasterBand(iCh);
////
////			report->file << "     - Band no. " << iCh << ":  ";
////			// to be done: implement a more efficient way of converting the data to it's original data type
////			switch (this->rasterBands[iCh - 1]->bandDataType) {
////			case SpEODT_Byte: {
////				report->file << "Convert data to SpEOByte..\n";
////				SpEOByte rasterBandData[(int) (this->sizeV) * (int) (this->sizeU)];
////				for (int i = 0; i < (int) (this->sizeV) * (int) (this->sizeU); i++)
////					rasterBandData[i] =
////							this->get_rasterBands()[iCh - 1]->bandData[i];
////				poBand->RasterIO(GF_Write, 0, 0, (int) (this->sizeV),
////						(int) (this->sizeU), rasterBandData, (int) (this->sizeV),
////						(int) (this->sizeU),
////						static_cast<GDALDataType>((int) (this->rasterBands[iCh - 1]->bandDataType)),
////						0, 0);
////				break;
////			}
////			case SpEODT_UInt16: {
////				report->file << "Convert data to SpEOUInt16.." << "\n";
////				SpEOUInt16 rasterBandData[(int) (this->sizeV) * (int) (this->sizeU)];
////				for (int i = 0; i < (int) (this->sizeV) * (int) (this->sizeU); i++)
////					rasterBandData[i] =
////							this->get_rasterBands()[iCh - 1]->bandData[i];
////				poBand->RasterIO(GF_Write, 0, 0, (int) (this->sizeV),
////						(int) (this->sizeU), rasterBandData, (int) (this->sizeV),
////						(int) (this->sizeU),
////						static_cast<GDALDataType>((int) (this->rasterBands[iCh - 1]->bandDataType)),
////						0, 0);
////				break;
////			}
////			case SpEODT_Int16: {
////				report->file << "Convert data to SpEOInt16.." << "\n";
////				SpEOInt16 rasterBandData[(int) (this->sizeV) * (int) (this->sizeU)];
////				//cout << endl << "" << this->get_rasterBands()[iCh - 1]-> << endl;
////				for (int i = 0; i < (int) (this->sizeV) * (int) (this->sizeU); i++)
////					rasterBandData[i] =
////							this->get_rasterBands()[iCh - 1]->bandData[i];
////				poBand->RasterIO(GF_Write, 0, 0, (int) (this->sizeV),
////						(int) (this->sizeU), rasterBandData, (int) (this->sizeV),
////						(int) (this->sizeU),
////						static_cast<GDALDataType>((int) (this->rasterBands[iCh - 1]->bandDataType)),
////						0, 0);
////				break;
////			}
////			case SpEODT_UInt32: {
////				report->file << "Convert data to SpEOUint32.." << "\n";
////				SpEOUInt32 rasterBandData[(int) (this->sizeV) * (int) (this->sizeU)];
////				for (int i = 0; i < (int) (this->sizeV) * (int) (this->sizeU); i++)
////					rasterBandData[i] =
////							this->get_rasterBands()[iCh - 1]->bandData[i];
////				poBand->RasterIO(GF_Write, 0, 0, (int) (this->sizeV),
////						(int) (this->sizeU), rasterBandData, (int) (this->sizeV),
////						(int) (this->sizeU),
////						static_cast<GDALDataType>((int) (this->rasterBands[iCh - 1]->bandDataType)),
////						0, 0);
////				break;
////			}
////			case SpEODT_Int32: {
////				report->file << "Convert data to SpEOInt32.." << "\n";
////				SpEOInt32 rasterBandData[(int) (this->sizeV) * (int) (this->sizeU)];
////				for (int i = 0; i < (int) (this->sizeV) * (int) (this->sizeU); i++)
////					rasterBandData[i] =
////							this->get_rasterBands()[iCh - 1]->bandData[i];
////				poBand->RasterIO(GF_Write, 0, 0, (int) (this->sizeV),
////						(int) (this->sizeU), rasterBandData, (int) (this->sizeV),
////						(int) (this->sizeU),
////						static_cast<GDALDataType>((int) (this->rasterBands[iCh - 1]->bandDataType)),
////						0, 0);
////				break;
////			}
////			default: {
////				report->file << "\n\n"
////						<< "CONVERSION FAILED BEFORE WRITING DATA TO FILE !" << "\n"
////						<< "\n" << "\n";
////				SpEOByte rasterBandData[(int) (this->sizeV) * (int) (this->sizeU)];
////				for (int i = 0; i < (int) (this->sizeV) * (int) (this->sizeU); i++)
////					rasterBandData[i] =
////							this->get_rasterBands()[iCh - 1]->bandData[i];
////				poBand->RasterIO(GF_Write, 0, 0, (int) (this->sizeV),
////						(int) (this->sizeU), rasterBandData, (int) (this->sizeV),
////						(int) (this->sizeU),
////						static_cast<GDALDataType>((int) (this->rasterBands[iCh - 1]->bandDataType)),
////						0, 0);
////				break;
////			}
////			}
////		}
//		/* Once we're done, close the dataset */
//		GDALClose((GDALDatasetH) poDstDS);
//		CSLDestroy(papszOptions);
//		report->file.close();
//
//



	//////////////////////////////////////////////////////////>
	///////////////////////////////////////////////////////////
	///////      try to write data using MPI IO        ////////
	///////////////////////////////////////////////////////////

//			    /*! Sixteen bit signed integer */		SpEODT_Int16 = 3,
//			    SpEOUInt16 rasterBandData[fDS*fDS];


		    // open data file
		    /////////////////////////////////////////////>
//			    poDstDS = poDriver->Create(this->get_fname(), (int) (this->sizeV),
//								(int) (this->sizeU), (int) (this->NCh),
//								static_cast<GDALDataType>((int) (this->rasterBands[0]->bandDataType)),
//								papszOptions);
//
//			    // write meta information in header file
//				if(my_rank==0){
//					report->file << " - Write appropriate metadata and raster data to created file...\n";
//				}
//				poDstDS->SetGeoTransform(this->get_geoTransform());
//
//				OGRSpatialReference oSRS;
//				char *pszSRS_WKT = NULL;
//				oSRS.SetUTM(11, TRUE);
//				oSRS.SetWellKnownGeogCS("NAD27");
//				oSRS.exportToWkt(&pszSRS_WKT);
//				poDstDS->SetProjection(pszSRS_WKT);
//				CPLFree(pszSRS_WKT);
	///////////////////////////////////////////////////////////
	///////      try to write data using MPI IO        ////////
	///////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////



	///////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////
	///////                    part 2                  ////////
	///////////////////////////////////////////////////////////
	/*
						// writing part OLD ///////////////////////////////////////////////////////////>
	//					// write area_sum to image file
	//					for(iChY = 0; iChY < glPrms->NChY; iChY++){
	//						// average the area of all touching patches
	//						area_sum[iChY] /= pCnt;
	//
	//						GDALRasterBand *poBand;
	//						poBand = poDstDS->GetRasterBand(iChY+1); // GDAL counts from 1 instead of 0!
	//
	//						if(iL==0){
	//							report->file << "     - Band no. " << iChY+1 << ":  ";
	//						}
	//						// to be done: implement a more efficient way of converting the data to it's original data type
	//						switch (this->rasterBands[iChY]->bandDataType) {
	//						case SpEODT_UInt16: {
	//							if(iL==0){
	//								report->file << "Convert data to SpEOUInt16.." << "\n";
	//							}
	//							SpEOUInt16 rasterBandData[fDS*fDS];
	//							for (int i = 0; i < fDS*fDS; i++){
	//								rasterBandData[i] = area_sum[iChY](i/fDS,i%fDS);//this->get_rasterBands()[iChY]->bandData[i];
	//							}
	//							poBand->RasterIO(GF_Write, vH, uH, fDS, fDS, rasterBandData, fDS, fDS,
	//									static_cast<GDALDataType>((int) (this->rasterBands[iChY]->bandDataType)),
	//									0, 0);
	//							break;
	//						}
	//						default: {
	//							if(iL==0){
	//								report->file << "\n\n"
	//											 << "CONVERSION FAILED BEFORE WRITING DATA TO FILE !" << "\n"
	//											 << "\n" << "\n";
	//							}
	//							SpEOByte rasterBandData[fDS*fDS];
	//							for (int i = 0; i < fDS*fDS; i++){
	//								rasterBandData[i] = area_sum[iChY](i/fDS,i%fDS);//this->get_rasterBands()[iChY]->bandData[i];
	//							}
	//							poBand->RasterIO(GF_Write, vH, uH, fDS, fDS, rasterBandData, fDS, fDS,
	//									static_cast<GDALDataType>((int) (this->rasterBands[iChY]->bandDataType)),
	//									0, 0);
	//							break;
	//						}
	//						}
	//					}
						//</////////////////////////////////////////////////////////// writing part OLD
	///////////////////////////////////////////////////////////
	///////                    part 2                  ////////
	///////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////<
	*/




}


/*
void SpEODataset::dataWriteParMPIIO(SpEOReport *report, SpEODataIOSetting *dSet, SpEOFusionSetting *fSet, SpEOParallelSetting *pSet, SpEOGlobalParams *glPrms, SpEOPaths *paths, MPI_Comm comm_busy) {

	int my_rank; int my_processes;
	MPI_Comm_rank(comm_busy, &my_rank);
	MPI_Comm_size(comm_busy, &my_processes);

	int u,v,iL, uL, vL, uH, vH, nL, pszL, pszH, a, uP, vP, iP, pCnt, fDS, iChZ, uPH, vPH;
	fDS = glPrms->fDS;
	nL = glPrms->sizeUL*glPrms->sizeVL;
	pszL = fSet->patchsize;
	pszH = pszL*fDS;
	a = pszL-fSet->overlap;

	SpEOMatrixF myReadMat, patch;
	SpEOVectorF pVecHR;
	SpEOMatrixF *area_sum = new SpEOMatrixF[this->NCh];

	iL = my_rank;
	if(my_rank==0){
		cout << "\n"
				<< "###########################################################" << "\n"
				<< "##                                                       ##" << "\n"
				<< "##    Write image ImZ in parallel using..                ##" << "\n"
				<< "##      (1) MPI I/O                                      ##" << "\n"
				<< "##      (2) ENVI header file                             ##" << "\n"
				<< "##      (3) band interleaved by pixel (BIP)              ##" << "\n"
				<< "##      (4) 16 bit unsigned integer                      ##" << "\n"
				<< "##                                                       ##" << "\n"
				<< "###########################################################"
				<< "\n" << "\n";

		report->file.open(report->fileName.c_str(),
				fstream::in | fstream::out | fstream::app);
		report->file << "\n"
				<< "###########################################################" << "\n"
				<< "##                                                       ##" << "\n"
				<< "##    Write image ImZ in parallel using..                ##" << "\n"
				<< "##      (1) MPI I/O                                      ##" << "\n"
				<< "##      (2) ENVI header file                             ##" << "\n"
				<< "##      (3) band interleaved by pixel (BIP)              ##" << "\n"
				<< "##      (4) 16 bit unsigned integer                      ##" << "\n"
				<< "##                                                       ##" << "\n"
				<< "###########################################################"
				<< "\n" << "\n";
	}
	// MPI IO declarations and initializations
	unsigned short *buf_spectrum_WR;
	buf_spectrum_WR = (unsigned short *)malloc(this->NCh*sizeof(unsigned short));
	MPI_File fh;
	MPI_Offset disp, offset;
	MPI_Datatype etype, ftype, buftype;
	MPI_Status status;
	int result;
	etype = MPI_UNSIGNED_SHORT;
	ftype = MPI_UNSIGNED_SHORT;
	//char *datarep = (char*)"native";
	buftype = MPI_UNSIGNED_SHORT;

	//DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD
	//if(my_rank==0){ cout << "debug 1" << endl;}
	//DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD
	// communicator group comm_busy opens the data file for writing only (and creating, if necessary) 
	result = MPI_File_open(comm_busy, (char*)paths->fname_ImZ.c_str(), MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &fh);
	if(result != MPI_SUCCESS){
		sample_error(result, (char*)"MPI_File_open");
	}
	disp=0;
	result = MPI_File_set_view(fh, disp, etype, ftype, (char*)"native", MPI_INFO_NULL);
	//int MPI_File_set_view(MPI_File fh, MPI_Offset disp, MPI_Datatype etype, MPI_Datatype, filetype, char *datarep, MPI_Info info)!
	
	//DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD
	//if(my_rank==0){ cout << "debug 2" << endl;}
	//DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD

	
	if(result != MPI_SUCCESS){
		sample_error(result, (char*)"MPI_File_set_view");
	}
//	while(iL<=fSet->pLast){ // as long as there are LR pixels to work on
	iL = dSet->uLFirst*glPrms->sizeVL+dSet->vLFirst + my_rank;
	nL = dSet->uLLast*glPrms->sizeVL+dSet->vLLast +1;

	//DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD
	//if(my_rank==0){ cout << "debug 3" << endl;}
	//DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD
	while(iL<nL){ // as long as there are LR pixels to work on

		// calculate the coordinates of the patch
		uL = iL / glPrms->sizeVL;
		vL = iL % glPrms->sizeVL;
		uH = uL*fDS;
		vH = vL*fDS;

		if(uL>=dSet->uLFirst && vL>=dSet->vLFirst && uL<=dSet->uLLast && vL<=dSet->vLLast){

			//DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD
			//cout << "[" << my_rank << "] debug 4, write LR pixel: uL=" << uL << ", vL=" << vL << ", iL=" << iL << endl;
			//DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD

//			cout << "[" << my_rank << "] write LR pixel: uL=" << uL << ", vL=" << vL << ", iL=" << iL << endl;
			// initialization
			for(iChZ=0; iChZ<this->NCh; iChZ++){
				area_sum[iChZ] = SpEOMatrixF::Zero(fDS, fDS);
			}
			// first patch that influences the current part (LR pixel) of the image that is to be written
			pCnt = 0;
			uP = ceil(max(0.0, ((double)(uL-pszL+1)))/a);
			vP = ceil(max(0.0, ((double)(vL-pszL+1)))/a);
			if(uP<glPrms->uPFirst){ uP=glPrms->uPFirst; }
			if(vP<glPrms->vPFirst){	vP=glPrms->vPFirst; }
			iP = uP*glPrms->NPV+vP;

			//DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD
			//cout << "[" << my_rank << "] debug 5, iP=" << iP << ", uP=" << uP << ", vP=" << vP << endl;
			//DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD

			while(uP*a <= uL  &&  uP*a+pszL-1 < glPrms->sizeUL){
			  
				//DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD
				//cout << "[" << my_rank << "] debug 6.1" << endl;
				//DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD
				// calc relative (local) position of area (part of patch) w.r.t. to the patch's corner coordinates (0,0)
				uPH = uH - uP*a*fDS;
				vP = ceil(max(0.0, ((double)(vL-pszL+1)))/a);
				if(vP<glPrms->vPFirst){
					vP=glPrms->vPFirst;
				}
				iP = uP*glPrms->NPV+vP;
				//DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD
				//cout << "[" << my_rank << "] debug 6.2, iP=" << iP << endl;
				//DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD
				while(vP*a <= vL  &&  vP*a+pszL-1 < glPrms->sizeVL){
					vPH = vH - vP*a*fDS;
					SpEODataset::read_patch(paths, iL, iP, uP, pCnt, vP, pszH, uPH, vPH, fDS, myReadMat, patch, pVecHR, area_sum, dSet, my_rank);
				}
				if(vL>=glPrms->sizeVL-pszL && vP == glPrms->NPV-1){
					vPH = vH - (glPrms->sizeVH-pszH);
					SpEODataset::read_patch(paths, iL, iP, uP, pCnt, vP, pszH, uPH, vPH, fDS, myReadMat, patch, pVecHR, area_sum, dSet, my_rank);
				}
				uP++;
				//DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD
				//cout << "[" << my_rank << "] debug 6.3" << endl;
				//DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD
			}
			if(uL>=glPrms->sizeUL-pszL && uP == glPrms->NPU-1){
				//DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD
				//cout << "[" << my_rank << "] debug 7.1" << endl;
				//DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD
				// calc relative (local) position of area (part of patch) w.r.t. to the patch's corner coordinates (0,0)
				uPH = uH - (glPrms->sizeUH-pszH);
				vP = ceil(max(0.0, ((double)(vL-pszL+1)))/a);
				if(vP<glPrms->vPFirst){
					vP=glPrms->vPFirst;
				}
				iP = uP*glPrms->NPV+vP;
				//DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD
				//cout << "[" << my_rank << "] debug 7.2, iP=" << iP << endl;
				//DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD
				while(vP*a <= vL  &&  vP*a+pszL-1 < glPrms->sizeVL){
					vPH = vH - vP*a*fDS;
					// read patch
					SpEODataset::read_patch(paths, iL, iP, uP, pCnt, vP, pszH, uPH, vPH, fDS, myReadMat, patch, pVecHR, area_sum, dSet, my_rank);
				}
				if(vL>=glPrms->sizeVL-pszL && vP == glPrms->NPV-1){
					vPH = vH - (glPrms->sizeVH-pszH);
					// read patch
					SpEODataset::read_patch(paths, iL, iP, uP, pCnt, vP, pszH, uPH, vPH, fDS, myReadMat, patch, pVecHR, area_sum, dSet, my_rank);
				}
				uP++;
				//DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD
				//cout << "[" << my_rank << "] debug 7.3" << endl;
				//DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD
			}
			// average the area of all touching patches
			for(iChZ = 0; iChZ < this->NCh; iChZ++){
				area_sum[iChZ] /= pCnt;
			}
			//DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD
			//cout << "[" << my_rank << "] debug 8" << endl;
			//DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD
			
			// write area_sum to image file
			if(dSet->uLFirst==0 && dSet->uLLast==glPrms->sizeUL-1 && dSet->vLFirst==0 && dSet->vLLast==glPrms->sizeVL-1){
				//DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD
				//cout << "[" << my_rank << "] debug 9" << endl;
				//DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD
				for(u=0; u<fDS; u++){// all HR lines/rows in the current area
					for(v=0; v<fDS; v++){// all HR samples/columns in the current area
						// Allocate and initialize a buffer (buf_spectrum_WR) containing this->NCh
						// unsigned shorts, where the unsigned short at location iChZ will set to (unsigned int)area_sum[iChZ](u,v).
						for(iChZ=0;iChZ<this->NCh;iChZ++){
							buf_spectrum_WR[iChZ] = (unsigned short)area_sum[iChZ](u,v);
						}
						// specify the offset relative to the displacement position. In compare to displacement, the offset value is given NOT IN BYTES, but IN MULTIPLES OF etype (here etype==MPI_UNSIGNED_SHORT)
						offset = this->NCh * ( (uL*this->sizeV*fDS)+this->sizeV*u + (vL*fDS) + v);

						// Write the buffer (buf_spectrum_WR)
						result = MPI_File_write_at(fh, offset, buf_spectrum_WR, this->NCh, buftype, &status);
						if(result != MPI_SUCCESS){
							sample_error(result, (char*)"MPI_File_write_at");
						}
					}
				}
			}else{
				//DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD
				//cout << "[" << my_rank << "] debug 9" << endl;
				//DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD
				for(u=0; u<fDS; u++){// all HR lines/rows in the current area
					for(v=0; v<fDS; v++){// all HR samples/columns in the current area
						// Allocate and initialize a buffer (buf_spectrum_WR) containing this->NCh
						// unsigned shorts, where the unsigned short at location iChZ will set to (unsigned int)area_sum[iChZ](u,v).
						for(iChZ=0;iChZ<this->NCh;iChZ++){
							buf_spectrum_WR[iChZ] = (unsigned short)area_sum[iChZ](u,v);
						}
						// specify the offset relative to the displacement position. In compare to displacement, the offset value is given NOT IN BYTES, but IN MULTIPLES OF etype (here etype==MPI_UNSIGNED_SHORT)
						//offset = this->NCh * ( (uL*this->sizeV*fDS)+this->sizeV*u + (vL*fDS) + v);
//						offset =   this->NCh*( (uL*fDS+u)*this->sizeV + vL*fDS+v);
						//offset = this->NCh * ( (((uL-dSet->uLFirst)*fDS+u)*this->sizeV) * (dSet->vLLast-dSet->vLFirst+1)*fDS + (vL-dSet->vLFirst)*fDS+v );
						offset = this->NCh * ( (((uL-dSet->uLFirst)*fDS+u)*glPrms->sizeVH_red) + (vL-dSet->vLFirst)*fDS+v );
						// Write the buffer (buf_spectrum_WR)
						result = MPI_File_write_at(fh, offset, buf_spectrum_WR, this->NCh, buftype, &status);
						if(result != MPI_SUCCESS){
							sample_error(result, (char*)"MPI_File_write_at");
						}
					}
				}
			}
			//DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD
			//cout << "[" << my_rank << "] debug 10" << endl;
			//DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD
			for(iChZ=0; iChZ<this->NCh; iChZ++){
				area_sum[iChZ] = SpEOMatrixF::Zero(fDS,fDS);
			}
			//DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD
			//cout << "[" << my_rank << "] debug 11" << endl;
			//DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD
		}
		iL += my_processes;
		//iL += pSet->parWrNumProc;
		//DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD
		//cout << "[" << my_rank << "] debug 12" << endl;
		//DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD
	}
	
	//DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD
	//cout << "[" << my_rank << "] debug 13" << endl;
	//DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD
	
	MPI_Barrier(comm_busy);
	result = MPI_File_close(&fh);
	if(result != MPI_SUCCESS){
		sample_error(result, (char*)"MPI_File_close");
	}
	if(my_rank==0){
		report->file.close();
	}
	delete[] area_sum;
	free(buf_spectrum_WR);
}
*/
void SpEODataset::dataWriteParMPIIO(SpEOReport *report, SpEODataIOSetting *dSet, SpEOFusionSetting *fSet, SpEOParallelSetting *pSet, SpEOGlobalParams *glPrms, SpEOPaths *paths, MPI_Comm comm_write) {


	int my_rank; int my_processes;
	MPI_Comm_rank(comm_write, &my_rank);
	MPI_Comm_size(comm_write, &my_processes);

	int u,v,iL, uL, vL, uH, vH, nL, pszL, pszH, a, uP, vP, iP, pCnt, fDS, iChZ, uPH, vPH;
	fDS = glPrms->fDS;
	nL = glPrms->sizeUL*glPrms->sizeVL;
	pszL = fSet->patchsize;
	pszH = pszL*fDS;
	a = pszL-fSet->overlap;

	SpEOMatrixF myReadMat, patch;
	SpEOVectorF pVecHR;
	SpEOMatrixF *area_sum = new SpEOMatrixF[this->NCh];

	if(my_rank==0){
		cout << "\n"
				<< "########################################################################" << "\n"
				<< "##                                                                    ##" << "\n"
				<< "##    Write image ImZ (from tmp patches on drive) in parallel using.. ##" << "\n"
				<< "##      (1) MPI I/O                                                   ##" << "\n"
				<< "##      (2) ENVI header file                                          ##" << "\n"
				<< "##      (3) band interleaved by pixel (BIP)                           ##" << "\n"
				<< "##      (4) 16 bit unsigned integer                                   ##" << "\n"
				<< "##                                                                    ##" << "\n"
				<< "########################################################################"
				<< "\n" << "\n";

		report->file.open(report->fileName.c_str(),
				fstream::in | fstream::out | fstream::app);
		report->file << "\n"
				<< "########################################################################" << "\n"
				<< "##                                                                    ##" << "\n"
				<< "##    Write image ImZ (from tmp patches on drive) in parallel using.. ##" << "\n"
				<< "##      (1) MPI I/O                                                   ##" << "\n"
				<< "##      (2) ENVI header file                                          ##" << "\n"
				<< "##      (3) band interleaved by pixel (BIP)                           ##" << "\n"
				<< "##      (4) 16 bit unsigned integer                                   ##" << "\n"
				<< "##                                                                    ##" << "\n"
				<< "########################################################################"
				<< "\n" << "\n";
	}
	// MPI IO declarations and initializations
	unsigned short *buf_spectrum_WR;
	buf_spectrum_WR = (unsigned short *)malloc(this->NCh*sizeof(unsigned short));
	MPI_File fh;
	MPI_Offset disp, offset;
	MPI_Datatype etype, ftype, buftype;
	MPI_Status status;
	int result;
	etype = MPI_UNSIGNED_SHORT;
	ftype = MPI_UNSIGNED_SHORT;
	//char *datarep = (char*)"native";
	buftype = MPI_UNSIGNED_SHORT;

	// communicator group comm_write opens the data file for writing only (and creating, if necessary) */
	result = MPI_File_open(comm_write, (char*)paths->fname_ImZ.c_str(), MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &fh);
	if(result != MPI_SUCCESS){
		sample_error(result, (char*)"MPI_File_open");
	}
	disp=0;
	result = MPI_File_set_view(fh, disp, etype, ftype, (char*)"native", MPI_INFO_NULL);
	//int MPI_File_set_view(MPI_File fh, MPI_Offset disp, MPI_Datatype etype, MPI_Datatype, filetype, char *datarep, MPI_Info info)!
	
	
	if(result != MPI_SUCCESS){
		sample_error(result, (char*)"MPI_File_set_view");
	}
//	if(my_rank < pSet->parWrNumProc){
//		iL = my_rank;
//	while(iL<=fSet->pLast){ // as long as there are LR pixels to work on
	iL = dSet->uLFirst*glPrms->sizeVL+dSet->vLFirst + my_rank;
	nL = dSet->uLLast*glPrms->sizeVL+dSet->vLLast +1;

	while(iL<nL){ // as long as there are LR pixels to work on

		// calculate the coordinates of the patch
		uL = iL / glPrms->sizeVL;
		vL = iL % glPrms->sizeVL;
		uH = uL*fDS;
		vH = vL*fDS;

		if(uL>=dSet->uLFirst && vL>=dSet->vLFirst && uL<=dSet->uLLast && vL<=dSet->vLLast){

			//cout << "[" << my_rank << "] write LR pixel: uL=" << uL << ", vL=" << vL << ", iL=" << iL << endl;
			// initialization
			for(iChZ=0; iChZ<this->NCh; iChZ++){
				area_sum[iChZ] = SpEOMatrixF::Zero(fDS, fDS);
			}
			// first patch that influences the current part (LR pixel) of the image that is to be written
			pCnt = 0;
			uP = ceil(max(0.0, ((double)(uL-pszL+1)))/a);
			vP = ceil(max(0.0, ((double)(vL-pszL+1)))/a);
			if(uP<glPrms->uPFirst){ uP=glPrms->uPFirst; }
			if(vP<glPrms->vPFirst){	vP=glPrms->vPFirst; }
			iP = uP*glPrms->NPV+vP;

			while(uP*a <= uL  &&  uP*a+pszL-1 < glPrms->sizeUL){
				// calc relative (local) position of area (part of patch) w.r.t. to the patch's corner coordinates (0,0)
				uPH = uH - uP*a*fDS;
				vP = ceil(max(0.0, ((double)(vL-pszL+1)))/a);
				if(vP<glPrms->vPFirst){
					vP=glPrms->vPFirst;
				}
				iP = uP*glPrms->NPV+vP;
				while(vP*a <= vL  &&  vP*a+pszL-1 < glPrms->sizeVL){
					vPH = vH - vP*a*fDS;
					SpEODataset::read_patch(paths, iL, iP, uP, pCnt, vP, pszH, uPH, vPH, fDS, myReadMat, patch, pVecHR, area_sum, dSet, my_rank);
				}
				if(vL>=glPrms->sizeVL-pszL && vP == glPrms->NPV-1){
					vPH = vH - (glPrms->sizeVH-pszH);
					SpEODataset::read_patch(paths, iL, iP, uP, pCnt, vP, pszH, uPH, vPH, fDS, myReadMat, patch, pVecHR, area_sum, dSet, my_rank);
				}
				uP++;
			}
			if(uL>=glPrms->sizeUL-pszL && uP == glPrms->NPU-1){
				// calc relative (local) position of area (part of patch) w.r.t. to the patch's corner coordinates (0,0)
				uPH = uH - (glPrms->sizeUH-pszH);
				vP = ceil(max(0.0, ((double)(vL-pszL+1)))/a);
				if(vP<glPrms->vPFirst){
					vP=glPrms->vPFirst;
				}
				iP = uP*glPrms->NPV+vP;
				while(vP*a <= vL  &&  vP*a+pszL-1 < glPrms->sizeVL){
					vPH = vH - vP*a*fDS;
					// read patch
					SpEODataset::read_patch(paths, iL, iP, uP, pCnt, vP, pszH, uPH, vPH, fDS, myReadMat, patch, pVecHR, area_sum, dSet, my_rank);
				}
				if(vL>=glPrms->sizeVL-pszL && vP == glPrms->NPV-1){
					vPH = vH - (glPrms->sizeVH-pszH);
					// read patch
					SpEODataset::read_patch(paths, iL, iP, uP, pCnt, vP, pszH, uPH, vPH, fDS, myReadMat, patch, pVecHR, area_sum, dSet, my_rank);
				}
				uP++;
			}
			// average the area of all touching patches
			for(iChZ = 0; iChZ < this->NCh; iChZ++){
				area_sum[iChZ] /= pCnt;
			}
			// write area_sum to image file
			if(dSet->uLFirst==0 && dSet->uLLast==glPrms->sizeUL-1 && dSet->vLFirst==0 && dSet->vLLast==glPrms->sizeVL-1){
				for(u=0; u<fDS; u++){// all HR lines/rows in the current area
					for(v=0; v<fDS; v++){// all HR samples/columns in the current area
						// Allocate and initialize a buffer (buf_spectrum_WR) containing this->NCh
						// unsigned shorts, where the unsigned short at location iChZ will set to (unsigned int)area_sum[iChZ](u,v).
						for(iChZ=0;iChZ<this->NCh;iChZ++){
							buf_spectrum_WR[iChZ] = (unsigned short)area_sum[iChZ](u,v);
						}
						// specify the offset relative to the displacement position. In compare to displacement, the offset value is given NOT IN BYTES, but IN MULTIPLES OF etype (here etype==MPI_UNSIGNED_SHORT)
						offset = this->NCh * ( (uL*this->sizeV*fDS)+this->sizeV*u + (vL*fDS) + v);

						// Write the buffer (buf_spectrum_WR)
						result = MPI_File_write_at(fh, offset, buf_spectrum_WR, this->NCh, buftype, &status);
						if(result != MPI_SUCCESS){
							sample_error(result, (char*)"MPI_File_write_at");
						}
					}
				}
			}else{
				for(u=0; u<fDS; u++){// all HR lines/rows in the current area
					for(v=0; v<fDS; v++){// all HR samples/columns in the current area
						// Allocate and initialize a buffer (buf_spectrum_WR) containing this->NCh
						// unsigned shorts, where the unsigned short at location iChZ will set to (unsigned int)area_sum[iChZ](u,v).
						for(iChZ=0;iChZ<this->NCh;iChZ++){
							buf_spectrum_WR[iChZ] = (unsigned short)area_sum[iChZ](u,v);
						}
						// specify the offset relative to the displacement position. In compare to displacement, the offset value is given NOT IN BYTES, but IN MULTIPLES OF etype (here etype==MPI_UNSIGNED_SHORT)
						//offset = this->NCh * ( (uL*this->sizeV*fDS)+this->sizeV*u + (vL*fDS) + v);
//						offset =   this->NCh*( (uL*fDS+u)*this->sizeV + vL*fDS+v);
						//offset = this->NCh * ( (((uL-dSet->uLFirst)*fDS+u)*this->sizeV) * (dSet->vLLast-dSet->vLFirst+1)*fDS + (vL-dSet->vLFirst)*fDS+v );
						offset = this->NCh * ( (((uL-dSet->uLFirst)*fDS+u)*glPrms->sizeVH_red) + (vL-dSet->vLFirst)*fDS+v );
						// Write the buffer (buf_spectrum_WR)
						result = MPI_File_write_at(fh, offset, buf_spectrum_WR, this->NCh, buftype, &status);
						if(result != MPI_SUCCESS){
							sample_error(result, (char*)"MPI_File_write_at");
						}
					}
				}
			}
			for(iChZ=0; iChZ<this->NCh; iChZ++){
				area_sum[iChZ] = SpEOMatrixF::Zero(fDS,fDS);
			}
		}
		iL += my_processes;
	}

	MPI_Barrier(comm_write);
	result = MPI_File_close(&fh);
	if(result != MPI_SUCCESS){
		sample_error(result, (char*)"MPI_File_close");
	}
	if(my_rank==0){
		report->file.close();
	}
	delete[] area_sum;
	free(buf_spectrum_WR);
}

void SpEODataset::read_patch(SpEOPaths *paths, int iL, int &iP, int uP, int &pCnt, int &vP, int pszH, int uPH, int vPH, int fDS, SpEOMatrixF &myReadMat, SpEOMatrixF &patch, SpEOVectorF &pVecHR, SpEOMatrixF *area_sum, SpEODataIOSetting *dSet, int my_rank){
	// read patch
	char buf [paths->dir_tmp_patches.length()+42]; // length(dir) + length(/) + length(patch_u0000_v0016_iP000016.csv) + margin
	sprintf (buf, "%s/%06d/patch_u%04d_v%04d_iP%06d.csv", paths->dir_tmp_patches.c_str(), iP, uP, vP, iP);
	int stat_CSV_read = read_CSV(&myReadMat, buf, ',', 0);
	for(int i=0; stat_CSV_read==-1 && i<dSet->dir_tmp_patches_additional_num; i++){
		char buf_tmp [paths->dir_tmp_patches_additional[i].length()+50]; // length(dir) + length(/) + length(123456) + length(/) + length(patch_u1234_v1234_iP123456.csv) + 2 margin
		sprintf (buf_tmp, "%s/%06d/patch_u%04d_v%04d_iP%06d.csv", paths->dir_tmp_patches_additional[i].c_str(), iP, uP, vP, iP);
		stat_CSV_read = read_CSV(&myReadMat, buf_tmp, ',', 0);
		cout << "try to open tmp patch file: " << buf_tmp << endl;
//		if(stat_CSV_read==1){
//			cout << "[" << my_rank << "] Problem Resolved! File was found in this directory: '"<< buf_tmp << "'!" << endl;
//		}
	}
	if(stat_CSV_read==-1){

		//dSet->
		//if !(uP*a<dSet->uLFirst || vP*a<dSet->vLFirst){
		//				uP++;
		//			}
		//			while(vP*a<dSet->vLFirst){
		//				vP++;
		//			}

			cout << "[" << my_rank << "] ERROR while writing pixel iL=" << iL << ": The .csv file '"<< buf << "' does not exist in any of the specified TMP patch directories!" << endl;
		//}
	}else{
		// cut off negative coefficients
		myReadMat = (myReadMat.array()>0).select(myReadMat,0);
		pVecHR = SpEOVectorF::Zero(myReadMat.rows());
		for(int iChZ=0; iChZ<this->NCh; iChZ++) {
			  pVecHR = myReadMat.col(iChZ);
			  patch = SpEOMatrixF::Map(pVecHR.data(), pszH, pszH);
			  area_sum[iChZ] += patch.block(uPH,vPH,fDS,fDS);
		}
		pCnt++;
	}
	vP++; iP++;
}






























void SpEODataset::writePatch(SpEOVectorF vec, int psz, int posU, int posV,
		int band) {
	this->rasterBands[band]->bandDataMat.block(posU, posV, psz, psz) =
			this->rasterBands[band]->bandDataMat.block(posU, posV, psz, psz)
					+ SpEOMatrixF::Map(vec.data(), psz, psz);
}

void SpEODataset::replBandByCwiseCuotient(SpEOMatrixF * rasterMat, int band) {
	this->rasterBands[band]->bandDataMat =
			this->rasterBands[band]->bandDataMat.cwiseQuotient(*rasterMat);
}

void SpEODataset::cutNegCoeff() {
	for (int iCh = 0; iCh < this->NCh; iCh++) {
		for (int uP = 0; uP < this->sizeU; uP++) {
			for (int vP = 0; vP < this->sizeV; vP++) {
				if (this->rasterBands[iCh]->bandDataMat.coeff(uP, vP) <= 0) {
					this->rasterBands[iCh]->bandDataMat(uP, vP) = 0;
				}
			}
		}
	}
}

SpEODataset::~SpEODataset(){
//	cout << "delete dataset..." << endl;
	for (int iCh=0; iCh < this->NCh; iCh++) {
		//this->rasterBands[iCh]->bandDataMatD = SpEOMatrixD::Zero(this->sizeU, this->sizeV);
		//this->rasterBands[iCh]->bandDataMat = SpEOMatrixF::Zero(0,0);
//		cout << ".. reduce band " << iCh << "...";

//		this->rasterBands[iCh]->bandDataMat  = SpEOMatrixF::Zero(0,0); // not necessary. Implicitely done in TuxEigen destructor
//		this->rasterBands[iCh]->bandDataMatD = SpEOMatrixD::Zero(0,0); // not necessary. Implicitely done in TuxEigen destructor
//		switch(this->imFlag)
//		case imFlag_X: {// = 0,
//			this->rasterBands[iCh]->bandDataMat = SpEOMatrixF::Zero(0,0);
//			break;
//		}
//		imFlag_X_LR     = 1,
//		imFlag_Y        = 2,
//		imFlag_Z        = 3,
//		imFlag_Z_LR     = 4,
//		imFlag_Z_ref    = 5,
//		imFlag_X_sim    = 6,
//		imFlag_X_sim_LR = 7,
//		imFlag_Z_init   = 8
		delete this->rasterBands[iCh];
//		cout << "done" << endl;
	}
//	cout << " now do:  delete [] this->rasterBands ...";
	delete [] this->rasterBands;
	delete [] geoTransform;
//	cout << " done!" << endl;
}

SpEORasterBand::~SpEORasterBand(){
	this->bandDataMat  = SpEOMatrixF::Zero(0,0);
	this->bandDataMatD = SpEOMatrixD::Zero(0,0);
}

void setMetaInfo(SpEODataset *ImZ, SpEODataset *ImY, SpEODataset *ImX, SpEODataIOSetting *dSet, SpEOGlobalParams *glPrms) {
	ImZ->imFormat = ImY->imFormat;
	ImZ->imFormatLong = ImY->imFormatLong;
	ImZ->projectionRef = ImX->projectionRef;
	ImZ->sizeU = glPrms->sizeUH_red;//ImX->sizeU;
	ImZ->sizeV = glPrms->sizeVH_red;//ImX->sizeV;
	ImZ->NCh = dSet->chBundleLast-dSet->chBundleFirst+1;
	for (int i = 0; i < 6; i++) {
		ImZ->geoTransform[i] = ImX->geoTransform[i];
	}
	ImZ->rasterBands = new SpEORasterBand*[ImZ->NCh];
	for (int iChZ = 0; iChZ < ImZ->NCh; iChZ++) {
		ImZ->rasterBands[iChZ] = new SpEORasterBand();
		ImZ->rasterBands[iChZ]->bandDataMatD = SpEOMatrixD::Zero(ImZ->sizeU,	ImZ->sizeV);
		int iChY = min(dSet->chBundleFirst+iChZ,ImY->NCh-1);
		ImZ->rasterBands[iChZ]->bandDataType = ImY->rasterBands[iChY]->bandDataType;
		ImZ->rasterBands[iChZ]->nBand = ImY->rasterBands[iChY]->nBand;
		//ImZ->rasterBands[iChZ]->nBlockUSize 	= ImY->rasterBands[iChZ]->nBlockUSize;
		//ImZ->rasterBands[iChZ]->minMaxVal[0] = ImY->rasterBands[iChZ]->minMaxVal[0];
		//ImZ->rasterBands[iChZ]->minMaxVal[1] = ImY->rasterBands[iChZ]->minMaxVal[1];
		ImZ->rasterBands[iChZ]->dataType = ImY->rasterBands[iChY]->dataType;
		ImZ->rasterBands[iChZ]->colorInterp = ImY->rasterBands[iChY]->colorInterp;
	}
}
void prepGroupCalcForJSpFI( int firstBandY, int lastBandY, int panBand,
							SpEOGlobalParams &glPrms, SpEOFusionSetting &fSetting, SpEODataIOSetting &dSetting,
							SpEOMatrixD &SRF_orig, SpEOMatrixD &SRF, SpEOMatrixD &filter_coeff, SpEOMatrixD &gauss_filter,
							SpEODataset *ImX    , SpEODataset *ImX_LR,     SpEODataset *ImY,     SpEODataset *ImZ,
							SpEODataset *ImX_tmp, SpEODataset *ImX_LR_tmp, SpEODataset *ImY_tmp, SpEODataset *ImZ_tmp ){
	int my_rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

	// check input
	if (firstBandY<0 || lastBandY>=glPrms.NChY){
		if(my_rank==0){
			cerr << endl << "ERROR: The condition (firstBandY=>0 && lastBandY<glPrms.NChY) is NOT filfilled!: firstBandY=="<< firstBandY <<", lastBandY==" << lastBandY << ", and glPrms.NChY==" << glPrms.NChY << endl << endl;
		}
		MPI_Barrier(MPI_COMM_WORLD);
		exit(2);
	}
	int firstBand, lastBand;
	glPrms.NChY            = lastBandY-firstBandY+1;
	glPrms.NChZ            = glPrms.NChY;
	glPrms.Nc_vec[0]       = glPrms.NChY;
	fSetting.Nc            = glPrms.NChY;
	dSetting.chBundleFirst = 0;
	dSetting.chBundleLast  = lastBandY-firstBandY;
	SRF                    = SRF_orig.block(0,firstBandY,1,glPrms.NChY);
	/********************************
	 *  get ImX_tmp and ImX_LR_tmp  *
	 ********************************/
	// initialize ImX_tmp and ImX_LR_tmp (including meta data) from ImX and ImX_LR, respectively
	firstBand = 0;
	lastBand = 0;
	ImX_tmp->copy_from_dataset(firstBand, lastBand, ImX);
	ImX_LR_tmp->copy_from_dataset(firstBand, lastBand, ImX_LR);

	// calculate ImX_tmp and ImX_LR_tmp
	if (panBand==-1){
		// nothing to be done
	}else{
		// copy data from specified, previously reconstructed band in ImZ
		ImX_tmp->get_rasterBands()[0]->bandDataMat = ImZ->get_rasterBands()[panBand]->bandDataMatD.cast<float>();
		/******************
		 * low-pass filter and down-sample
		 ******************///
		SpEOMatrixD ImX_tmpD_F    = ImZ->get_rasterBands()[panBand]->bandDataMatD;
		SpEOMatrixD ImX_LR_tmpD_F;
		fast_filter    (ImX_LR_tmpD_F, ImX_tmpD_F, gauss_filter, filter_coeff, glPrms.fDS); // fast method (intuitive)
		ImX_LR_tmp->get_rasterBands()[0]->bandDataMat = ImX_LR_tmpD_F.cast<float>();
	}

	/********************************
	 *  get ImY_tmp                 *
	 ********************************/
	firstBand = firstBandY;
	lastBand = lastBandY;
	ImY_tmp->copy_from_dataset(firstBand, lastBand, ImY);

	/********************************
	 *  get ImZ_tmp                 *
	 ********************************/
	firstBand = firstBandY;
	lastBand = lastBandY;
	ImZ_tmp->copy_from_dataset(firstBand, lastBand, ImZ);

//	ImX; ImX_tmp
//	ImX_LR; ImX_LR_tmp
//	ImY; ImY_tmp
//	ImZ; ImZ_tmp
}
void SpEODataset::copy_from_dataset( int firstBand, int lastBand, SpEODataset *sourceIm){

	this->set_fname(sourceIm->get_fname());
	this->imFlag        = sourceIm->get_imFlag();
	this->sizeU         = sourceIm->get_sizeU();
	this->sizeV         = sourceIm->get_sizeV();
	this->NCh           = lastBand - firstBand + 1;//sourceIm->get_NCh();
	this->imFormat      = sourceIm->get_imFormat();
	this->imFormatLong  = sourceIm->get_imFormatLong();
	this->projectionRef = sourceIm->get_projectionRef();

	this->rasterBands   = new SpEORasterBand*[this->NCh];
	for (int iCh = 0; iCh < this->NCh; iCh++) {
		this->rasterBands[iCh] = new SpEORasterBand();

		this->rasterBands[iCh]->bandDataType = sourceIm->get_rasterBands()[firstBand+iCh]->get_bandDataType();
		this->rasterBands[iCh]->nBand        = iCh+1;
		this->rasterBands[iCh]->minMaxVal[0] = sourceIm->get_rasterBands()[firstBand+iCh]->get_minMaxVal()[0];
		this->rasterBands[iCh]->minMaxVal[1] = sourceIm->get_rasterBands()[firstBand+iCh]->get_minMaxVal()[1];
		this->rasterBands[iCh]->dataType     = sourceIm->get_rasterBands()[firstBand+iCh]->get_dataType();
		this->rasterBands[iCh]->colorInterp  = sourceIm->get_rasterBands()[firstBand+iCh]->get_colorInterp();
		if(sourceIm->get_imFlag()==imFlag_Z || sourceIm->get_imFlag()==imFlag_Z_LR){
			// attention: currently only the pointer is given to this->rasterBands[iCh]->bandDataMatD
			// to be done: actual data need to be copied.
			this->rasterBands[iCh]->bandDataMatD = SpEOMatrixD::Map(sourceIm->get_rasterBands()[firstBand+iCh]->get_bandDataMatD()->data(),
					                                                sourceIm->get_rasterBands()[firstBand+iCh]->get_bandDataMatD()->rows(),
					                                                sourceIm->get_rasterBands()[firstBand+iCh]->get_bandDataMatD()->cols());
		}else{
			// attention: currently only the pointer is given to this->rasterBands[iCh]->bandDataMatD
			// to be done: actual data need to be copied.
//			this->rasterBands[iCh]->bandDataMat = sourceIm->get_rasterBands()[firstBand+iCh]->get_bandDataMat();
			this->rasterBands[iCh]->bandDataMat = SpEOMatrixF::Map(sourceIm->get_rasterBands()[firstBand+iCh]->get_bandDataMat()->data(),
													               sourceIm->get_rasterBands()[firstBand+iCh]->get_bandDataMat()->rows(),
													               sourceIm->get_rasterBands()[firstBand+iCh]->get_bandDataMat()->cols());

		}
	}
}

void SpEODataset::copyMetaInfoFromDatasets(SpEODataset *sourceImSpatial, SpEODataset *sourceImSpectral, SpEODataset *sourceImFormat, SpEODataFormat dataFormat){

	this->set_fname(sourceImSpatial->get_fname());
//	this->imFlag        = SpEOImFlag;
	this->sizeU         = sourceImSpatial->get_sizeU();
	this->sizeV         = sourceImSpatial->get_sizeV();
	this->NCh           = sourceImSpectral->get_NCh();
	this->imFormat      = sourceImFormat->get_imFormat();
	this->imFormatLong  = sourceImFormat->get_imFormatLong();
	this->projectionRef = sourceImSpatial->get_projectionRef();

	this->rasterBands   = new SpEORasterBand*[this->NCh];
	for (int iCh = 0; iCh < this->NCh; iCh++) {
		int NChSoureImSpatial = sourceImSpatial->NCh;
		int iChSoureImSpatial = min(iCh,NChSoureImSpatial-1);
		this->rasterBands[iCh] = new SpEORasterBand();

		this->rasterBands[iCh]->bandDataType = sourceImSpatial->get_rasterBands()[iChSoureImSpatial]->get_bandDataType();
		this->rasterBands[iCh]->nBand        = iCh+1;
		this->rasterBands[iCh]->minMaxVal[0] = sourceImSpatial->get_rasterBands()[iChSoureImSpatial]->get_minMaxVal()[0];
		this->rasterBands[iCh]->minMaxVal[1] = sourceImSpatial->get_rasterBands()[iChSoureImSpatial]->get_minMaxVal()[1];
		this->rasterBands[iCh]->dataType     = sourceImSpatial->get_rasterBands()[iChSoureImSpatial]->get_dataType();
		this->rasterBands[iCh]->colorInterp  = sourceImSpatial->get_rasterBands()[iChSoureImSpatial]->get_colorInterp();

		if (dataFormat==SpEODouble){
		//if(sourceIm->get_imFlag()==imFlag_Z || sourceImSpatial->get_imFlag()==imFlag_Z_LR){
			// attention: currently only the pointer is given to this->rasterBands[iCh]->bandDataMatD
			// to be done: actual data need to be copied.
			this->rasterBands[iCh]->bandDataMatD = SpEOMatrixD::Zero(sourceImSpatial->get_rasterBands()[iChSoureImSpatial]->get_bandDataMatD()->rows(),
					                                                 sourceImSpatial->get_rasterBands()[iChSoureImSpatial]->get_bandDataMatD()->cols());
		}else{
			// attention: currently only the pointer is given to this->rasterBands[iCh]->bandDataMatD
			// to be done: actual data need to be copied.
//			this->rasterBands[iCh]->bandDataMat = sourceIm->get_rasterBands()[iCh]->get_bandDataMat();
//			this->rasterBands[iCh]->bandDataMat = SpEOMatrixF::Map(sourceIm->get_rasterBands()[iCh]->get_bandDataMat()->data(),
//													               sourceIm->get_rasterBands()[iCh]->get_bandDataMat()->rows(),
//													               sourceIm->get_rasterBands()[iCh]->get_bandDataMat()->cols());
			this->rasterBands[iCh]->bandDataMat = SpEOMatrixF::Zero(sourceImSpatial->get_rasterBands()[iChSoureImSpatial]->get_bandDataMatD()->rows(),
								                                    sourceImSpatial->get_rasterBands()[iChSoureImSpatial]->get_bandDataMatD()->cols());
		}
	}
}


void SpEODataset::fill_band_data(SpEODataset *sourceIm, int firstBand, int lastBand){
	int my_rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

	if(lastBand-firstBand+1 != sourceIm->NCh){
		if(my_rank==0){
			cerr << endl << "ERROR: in SpEODataset::fill_band_data(..): The number of bands of the source image has to be equal to lastBand-firstBand+1" << endl << endl;
		}
		MPI_Barrier(MPI_COMM_WORLD);
		exit(2);
	}

	for (int iCh = firstBand; iCh <= lastBand; iCh++) {
		this->rasterBands[iCh] = new SpEORasterBand();

		this->rasterBands[iCh]->bandDataType = sourceIm->get_rasterBands()[iCh-firstBand]->get_bandDataType();
		this->rasterBands[iCh]->nBand        = iCh+1;
		this->rasterBands[iCh]->minMaxVal[0] = sourceIm->get_rasterBands()[iCh-firstBand]->get_minMaxVal()[0];
		this->rasterBands[iCh]->minMaxVal[1] = sourceIm->get_rasterBands()[iCh-firstBand]->get_minMaxVal()[1];
		this->rasterBands[iCh]->dataType     = sourceIm->get_rasterBands()[iCh-firstBand]->get_dataType();
		this->rasterBands[iCh]->colorInterp  = sourceIm->get_rasterBands()[iCh-firstBand]->get_colorInterp();
		if(sourceIm->get_imFlag()==imFlag_Z || sourceIm->get_imFlag()==imFlag_Z_LR){
			// attention: currently only the pointer is given to this->rasterBands[iCh]->bandDataMatD
			// to be done: actual data need to be copied.
			this->rasterBands[iCh]->bandDataMatD = SpEOMatrixD::Map(sourceIm->get_rasterBands()[iCh-firstBand]->get_bandDataMatD()->data(),
					                                                sourceIm->get_rasterBands()[iCh-firstBand]->get_bandDataMatD()->rows(),
					                                                sourceIm->get_rasterBands()[iCh-firstBand]->get_bandDataMatD()->cols());
		}else{
			// attention: currently only the pointer is given to this->rasterBands[iCh]->bandDataMatD
			// to be done: actual data need to be copied.
//			this->rasterBands[iCh]->bandDataMat = sourceIm->get_rasterBands()[iCh-firstBand]->get_bandDataMat();
			this->rasterBands[iCh]->bandDataMat = SpEOMatrixF::Map(sourceIm->get_rasterBands()[iCh-firstBand]->get_bandDataMat()->data(),
													               sourceIm->get_rasterBands()[iCh-firstBand]->get_bandDataMat()->rows(),
													               sourceIm->get_rasterBands()[iCh-firstBand]->get_bandDataMat()->cols());



		}
	}

}

void SpEODataset::shift_image_bandwise(double *&Im_shift){
	for (int iCh = 0; iCh < this->get_NCh(); iCh++) {
		this->rasterBands[iCh]->bandDataMat.array() += (float)(Im_shift[iCh]);
	}
}






void cutRelevantInput(SpEODataset *ImX, SpEODataset *ImY, SpEODataIOSetting *dSet){

	// currently only spectral subset of ImY supported. The spatial subset is a possible future feature.

	int my_rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

	if(my_rank==0){
		cout<< "\n"
		<< "###########################################################\n"
		<< "##                                                       ##\n"
		<< "##           crop relevant spectral and spatial          ##\n"
		<< "##                  subset of input data                 ##\n"
		<< "##                                                       ##\n"
		<< "###########################################################"
		<< "\n" << "\n";
	}
	
	if(dSet->chBundleFirst < 0){
		if(my_rank==0) {
			cout << "WARNING: chBundleFirst must not be negative! It got corrected from " << dSet->chBundleFirst << " to 0." << endl;
		}
		dSet->chBundleFirst = 0;
	}
	if(dSet->chBundleFirst >= ImY->NCh){
		if(my_rank==0) {
			cout << "WARNING: chBundleFirst too great. It got corrected from " << dSet->chBundleFirst << " to " << ImY->NCh-1 << "=NChY-1," << endl;
		}
		dSet->chBundleFirst = ImY->NCh-1;
	}
	if(dSet->chBundleLast < 0){
		if(my_rank==0) {
			cout << "WARNING: chBundleLast must not be negative! It got corrected from " << dSet->chBundleLast << " to 0." << endl;
		}
		dSet->chBundleLast = 0;
	}
	if(dSet->chBundleLast >= ImY->NCh){
		if(my_rank==0) {
			cout << "WARNING: chBundleLast too great. It got corrected from " << dSet->chBundleLast << " to " << ImY->NCh-1 << "=NChY-1." << endl;
		}
		dSet->chBundleLast = ImY->NCh-1;
	}

	short NCh_buf = ImY->NCh;
	if(dSet->chBundleLast-dSet->chBundleFirst+1 < ImY->NCh){
		if(my_rank==0) {
			cout << "WARNING: NChY got corrected from " << ImY->NCh << " to " << dSet->chBundleLast-dSet->chBundleFirst+1 << "=chBundleLast-chBundleFirst+1." << endl;
		}
		ImY->NCh = dSet->chBundleLast-dSet->chBundleFirst+1;
	}
	if(!(dSet->chBundleFirst==0 && dSet->chBundleLast==NCh_buf-1)){
		for(short iChY=0; iChY<ImY->NCh; iChY++){
			ImY->rasterBands[iChY] = ImY->rasterBands[iChY+dSet->chBundleFirst];
		}
	}
}

void cutRelevantInput(SpEOResolution eRes, SpEOImFlag imFlag, SpEODataset *Im, SpEODataIOSetting *dSet, SpEOGlobalParams *glPrms, bool cutOffLRBoundaryPixels){

	short NCh_buf = Im->NCh;
	if(!(imFlag==imFlag_X || imFlag==imFlag_X_LR || imFlag==imFlag_X_sim || imFlag==imFlag_X_sim_LR)){
		if(dSet->chBundleLast-dSet->chBundleFirst+1 < Im->NCh){
			Im->NCh = dSet->chBundleLast-dSet->chBundleFirst+1;
		}
		if(!(dSet->chBundleFirst==0 && dSet->chBundleLast==NCh_buf-1)){
			for(short iCh=0; iCh<Im->NCh; iCh++){
				Im->rasterBands[iCh] = Im->rasterBands[iCh+dSet->chBundleFirst];
	//			SpEORasterBand *tmp  = Im->rasterBands[iCh+dSet->chBundleFirst];
	//			Im->rasterBands[iCh] = tmp;
			}
		}
	}

	int mySizeU_red;
	int mySizeV_red;
	int myUFirst;
	int myVFirst;

	if(eRes==LR){
		mySizeU_red = glPrms->sizeUL_red;
		mySizeV_red = glPrms->sizeVL_red;
		myUFirst = dSet->uLFirst;
		myVFirst = dSet->vLFirst;
		if(cutOffLRBoundaryPixels){
			mySizeU_red -= 2;
			mySizeV_red -= 2;
			myUFirst += 1;
			myVFirst += 1;
		}
	}else if(eRes==HR){
		mySizeU_red = glPrms->sizeUH_red;
		mySizeV_red = glPrms->sizeVH_red;
		myUFirst = dSet->uLFirst*glPrms->fDS;
		myVFirst = dSet->vLFirst*glPrms->fDS;
		if(cutOffLRBoundaryPixels){
			mySizeU_red -= 2*glPrms->fDS;
			mySizeV_red -= 2*glPrms->fDS;
			myUFirst += glPrms->fDS;
			myVFirst += glPrms->fDS;
		}
	}else{
		cout << "WARNING: First argument in the function cutRelevantInput has to specify the image resolution (set either to 'LR' of 'HR')! It got set to 'LR' by default!" << endl;
		mySizeU_red = glPrms->sizeUL_red;
		mySizeV_red = glPrms->sizeVL_red;
		myUFirst = dSet->uLFirst;
		myVFirst = dSet->vLFirst;
	}
	if(!(  dSet->uLFirst==0 && dSet->uLLast==glPrms->sizeUL-1		
		&& dSet->vLFirst==0 && dSet->vLLast==glPrms->sizeVL-1)){
		for(int iCh=0; iCh<Im->NCh; iCh++){
//				cout << "myUFirst=" << myUFirst <<", myVFirst="<< myVFirst << ", mySizeU_red="<<mySizeU_red<<", mySizeV_red="<<mySizeV_red<< endl << endl;
//				cout << "dSet->uLFirst=" << dSet->uLFirst <<", dSet->vLFirst="<< dSet->vLFirst << ", dSet->uLLast="<<dSet->uLLast<<", dSet->vLLast="<<dSet->vLLast<< endl << endl;
//				cout << "glPrms->sizeUL=" << glPrms->sizeUL << ", glPrms->sizeVL=" << glPrms->sizeVL << endl << endl;
//				cout << "Im->sizeU=" << Im->sizeU << ", Im->sizeV=" << Im->sizeV << endl << endl;

//				SpEOMatrixF tmp = Im->get_rasterBands()[iCh]->bandDataMat.block(dSet->uLFirst*glPrms->fDS,dSet->vLFirst*glPrms->fDS, glPrms->sizeUH_red, glPrms->sizeVH_red);
				SpEOMatrixF tmp = Im->get_rasterBands()[iCh]->bandDataMat.block(myUFirst, myVFirst, mySizeU_red, mySizeV_red);
				Im->get_rasterBands()[iCh]->bandDataMat = tmp;
		}
		Im->sizeU = mySizeU_red; // glPrms->sizeUH_red;
		Im->sizeV = mySizeV_red; // glPrms->sizeVH_red;
	}
}




















string get_current_time() {
	char curtime[50];
	time_t rawtime;
	struct tm * timeinfo;
	time(&rawtime);
	timeinfo = localtime(&rawtime);
	//strftime (curtime,50,"%y%m%d",timeinfo);
	strftime(curtime, 50, "%y%m%d_%H%M%S", timeinfo);
	string curtimeStr(curtime);
	return curtimeStr;
}

//void JSparseFI_reportInit(SpEOPaths *paths, SpEOFusionSetting *fSetting, SpEOSolverSetting *sSetting, SpEOGlobalParams *glPrms, ofstream& reportFile, string *reportFileName, string curTime){
void SpEOReport::initialize(SpEOPaths *paths,
							SpEODataIOSetting *dSetting,
                            SpEOFusionSetting *fSetting,
                            SpEOOutputSetting *oSetting,
                            SpEOSolverSetting *sSetting,
                            SpEOParallelSetting *pSetting,
                            int argc, char **argv) { //, ofstream& reportFile, string *reportFileName, string curTime){

	switch(fSetting->ImZ_init_type){
		  case 0:{
			  paths->fname_ImZ_init = paths->fname_ImZ_ref;
			  break;
		  }case 1:{
			  paths->fname_ImZ_init = paths->fname_ImZ_init_ImY_US;
			  break;
		  }case 2:{
			  paths->fname_ImZ_init = paths->fname_ImZ_init_rec;
			  break;
		  }
	  }


//	char *tmp_patches = new char[sizeof(paths->dir_tmp)+8];
//	tmp_patches[0] = *paths->dir_tmp;
//	tmp_patches[sizeof(paths->dir_tmp)] = "/patches";
//	mkdir(tmp_patches);
//	mkdir(paths->dir_tmp, 0777);


	int my_processes, my_rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
	MPI_Comm_size(MPI_COMM_WORLD, &my_processes);
	/*
	 * create folder and report file
	 */
	mkdir(paths->dir_out.c_str(), 0777);
	chmod(paths->dir_out.c_str(), S_IRWXU | S_IRWXG | S_IRWXO);
	this->curTimeSec = MPI_Wtime();
	this->curTime = get_current_time();
	string mydate = this->curTime.substr(0, 6);
	string mytime = this->curTime.substr(7, 6);
	paths->dir_out += "/" + mydate; // + "/" + "ws1_pGrp1";//"_ws0_NP128";
	mkdir(paths->dir_out.c_str(), 0777);
	chmod(paths->dir_out.c_str(), S_IRWXU | S_IRWXG | S_IRWXO);

	paths->dir_out += "/" + dSetting->jobName + "_" + dSetting->jobID; // + "/" + "ws1_pGrp1";//"_ws0_NP128";
	mkdir(paths->dir_out.c_str(), 0777);
	chmod(paths->dir_out.c_str(), S_IRWXU | S_IRWXG | S_IRWXO);

	stringstream dir_out_tmp;
	dir_out_tmp << paths->dir_out + "/"         // << "_NP" << glPrms->NP //			<< "_2step" << fSetting->two_step_estimation  << "_pTot" << my_processes //			<< "_ovlp"  << fSetting->overlap << "_NDP" << fSetting->NDP << "_lmd" << fSetting->lambda
			<< mytime
//			<< "_ID"    << paths->dataSetID      // << "_alg"  << fSetting->fMethod  << "_dS"   << fSetting->dictselect // << "_psz"   << fSetting->patchsize //			<< "_ovl"   << fSetting->overlap //			<< "_NDP"  << fSetting->NDP         // << "_2Step" << fSetting->two_step_estimation << "_Nc"    << fSetting->Nc << "_No" << fSetting->No
			<< "_ID"    << paths->dataSetID_str  // << "_alg"  << fSetting->fMethod  << "_dS"   << fSetting->dictselect // << "_psz"   << fSetting->patchsize //			<< "_ovl"   << fSetting->overlap //			<< "_NDP"  << fSetting->NDP         // << "_2Step" << fSetting->two_step_estimation << "_Nc"    << fSetting->Nc << "_No" << fSetting->No
//			<< "_itr"   << fSetting->iterMain // << "_NwCf"  << fSetting->useNewMethodForCalculatingZ // << "_ImXSm" << fSetting->useSimulatedImXforDictLearn
//			<< "_lX"    << fSetting->lambdaX_ABC
//			<< "_lY"    << fSetting->lambdaY_ABC
//			<< "_lZ"    << fSetting->lambdaZ_ABC //<< "_lZ0In1stIter" << fSetting->lambdaZ_ABC_in_1st_iter //		 << "_fllImOpt"<< fSetting->LQ_post_opt_im
			<< "_lXI"  << fSetting->lambdaX_im
			<< "_lYI"  << fSetting->lambdaY_im
			<< "_l"  << fSetting->lambda
			<< "_lXA"  << fSetting->lambdaX_ABC
			<< "_lYA"  << fSetting->lambdaY_ABC
			<< "_init"  << fSetting->ImZ_init_type
//			<< "_estR"   << fSetting->use_estimated_SRFs
			<< "_m"     << fSetting->ImX_sim_mode
			<< "_"      << fSetting->subspace_transform_type
			<< "_dim"   << fSetting->subspace_dim;
//			<< "_LQ"  << fSetting->LQ_post_opt
//			<< "_fA"  << sSetting->fix_Alpha
//			<< "_fm"  << sSetting->fix_delta_m
//			<< "_lY"  << fSetting->lambdaY
//			<< "_lX"  << fSetting->lambdaX; //
//			<< "_SpecNorm"  << fSetting->matrixNorm //			<< "_2Step" << fSetting->two_step_estimation << "_Nc"    << fSetting->Nc << "_No"    << fSetting->No << "_psz"   << fSetting->patchsize
//			<< "_pixwiseMean"  << fSetting->addMeanPixelwise; //			<< "_2Step" << fSetting->two_step_estimation << "_Nc"    << fSetting->Nc << "_No"    << fSetting->No << "_psz"   << fSetting->patchsize

	paths->dir_out = dir_out_tmp.str();
	mkdir(paths->dir_out.c_str(), 0777);
	chmod(paths->dir_out.c_str(), S_IRWXU | S_IRWXG | S_IRWXO);
	this->fileName = paths->dir_out + "/" + this->curTime + "_report.txt";
	this->file.open(this->fileName.c_str(),
			fstream::in | fstream::out | fstream::app);

	string datePretty = "20" + this->curTime.substr(0, 2) + "-"
			+ this->curTime.substr(2, 2) + "-" + this->curTime.substr(4, 2);
	string timePretty = this->curTime.substr(7, 2) + ":"
			+ this->curTime.substr(9, 2) + ":" + this->curTime.substr(11, 2);

	if(oSetting->saveAlphas){
		string alpha_dir = paths->dir_out + "/" + "patches";
		mkdir(alpha_dir.c_str(), 0777);
		chmod(alpha_dir.c_str(), S_IRWXU | S_IRWXG | S_IRWXO);
	}
	if(oSetting->saveDicts){
			string dict_dir = paths->dir_out + "/" + "dictCoords";
			mkdir(dict_dir.c_str(), 0777);
			chmod(dict_dir.c_str(), S_IRWXU | S_IRWXG | S_IRWXO);
		}

	/*
	 * write some information to file
	 */
	this->file
			<< "###############################################################################"
			<< "\n"
			<< "###                                                                         ###"
			<< "\n"
			<< "###                                                                         ###"
			<< "\n";
	switch (fSetting->fMethod) {
	case SparseFI: {
		this->file
				<< "###                                 SparseFI                                ###"
				<< "\n"
				<< "###             [ Sparse Fusion of Images ] for pan-sharpening              ###"
				<< "\n";
		break;
	}
	case JSparseFI: {
		this->file
				<< "###                                J-SparseFI                               ###"
				<< "\n"
				<< "###         [ Jointly Sparse Fusion of Images ] for pan-sharpening          ###"
				<< "\n";
		break;
	}
	case GroupedJSparseFI: {
		this->file
				<< "###                            GroupedJSparseFI                             ###"
				<< "\n"
				<< "###         [ Jointly Sparse Fusion of Images ] for pan-sharpening          ###"
				<< "\n";
		break;
	}
	case JSparseFIHM: {
		this->file
				<< "###                               J-SparseFI-HM                             ###"
				<< "\n"
				<< "###   [ Jointly Sparse Fusion of Hyperspectral and Multispectral Images ]   ###"
				<< "\n";
		break;
	}
	case LeastSquares: {
			this->file
				<< "###                           Least Squares                                 ###"
				<< "\n"
				<< "###                  - Dictionary contains only current patch -             ###"
				<< "\n";
			break;
		}
	}
	this->file
			<< "###                                                                         ###"
			<< "\n"
			<< "###                                                                         ###"
			<< "\n"
			<< "###############################################################################"
			<< "\n\n";

	this->file << "starting date: " << datePretty << " " << timePretty
			<< "\n";

	this->file << "\n"
			<< "###########################################################"
			<< "\n"
			<< "##                                                       ##"
			<< "\n"
			<< "##                      job details                      ##"
			<< "\n"
			<< "##                                                       ##"
			<< "\n"
			<< "###########################################################"
			<< "\n" << "\n"
			<< " - job name:  " << dSetting->jobName << "\n"
			<< " - job ID:    " << dSetting->jobID << "\n"
			<< "\n"
			<< "###########################################################"
			<< "\n"
			<< "##                                                       ##"
			<< "\n"
			<< "##                     user settings                     ##"
			<< "\n"
			<< "##                                                       ##"
			<< "\n"
			<< "###########################################################"
			<< "\n" << "\n"
			<< "*========================================*"
			<< "\n" << "* paths                                  *" << "\n"
			<< "*========================================*" << "\n"
//			<< " - dataSet ID:                " << paths->dataSetID << "\n"
			<< " - dataSet ID:                " << paths->dataSetID_str << "\n"
			<< " - directory input data:      " << paths->dir_in << "\n"
			<< " - file name ImX:             " << paths->fname_ImX << "\n"
//			<< " - file name ImX_LR:          " << paths->fname_ImX_LR << "\n"
			<< " - file name ImY:             " << paths->fname_ImY << "\n"
			<< " - file name ImZ (reference): " << paths->fname_ImZ_ref << "\n"
			<< " - file name ImZ_init:        " << paths->fname_ImZ_init << "\n"
			<< " - directory output data:     " << paths->dir_out << "\n"
			<< " - file name ImZ:             " << paths->fname_ImZ << "\n"
	        << " - ### use_estimated_SRFs:                " << fSetting->use_estimated_SRFs << "\n"
			<< " - file name SRF:                         " << paths->fname_SRF << "\n"
			<< " - file name SRF (for spectral grouping): " << paths->fname_SRF_for_Spectral_Grouping << "\n"
			<< " - ImZ_init_type = " << fSetting->ImZ_init_type << "(";
	switch(fSetting->ImZ_init_type){
		case 0:{
			this->file << "lambdaZ_ABZ=0 in 1st iter)\n";
			break;
		}
		case 1:{
			this->file << "upsampled and bilinearly interpolated low resolution image ImY)\n";
			break;
		}
		case 2:{
			this->file << "rec. result of another algorithm)\n";
			break;
		}
	}

	if(dSetting->contUnfinishedRec){
		this->file  << " - file name CSV containing iPs of incomplete set of patches from previous  unfinished reconstruction: " << "\n"
			<< "                                  " << paths->PathToIncompletePatchSetCSV << "\n"
			<< " - directories that may containe previously reconstructed TMP patches: " << "\n";
		for(int i=0; i<dSetting->dir_tmp_patches_additional_num; i++){
			this->file << "              .. directory no. " << i+1 << ": "<< paths->dir_tmp_patches_additional[i] << "\n";
		}
	}
	this->file << " - directory where tmp patches were/are stored : "<< paths->dir_tmp_patches << "\n";
	this->file 	<< "\n"
			<< "*========================================*" << "\n"
			<< "* image data I/O parameters              *" << "\n"
			<< "*========================================*" << "\n"
			<< " - first Y band to be processed (before possible correction): " << dSetting->chBundleFirst << "\n"
			<< " - last Y band to be processed (before possible correction): " << dSetting->chBundleLast << "\n"
			<< " - uLFirst: first vertical low resolution pixel to be processed (before possible correction):   " << dSetting->uLFirst << "\n"
			<< " - uLLast:  last vertical low resolution pixel to be processed (before possible correction):    " << dSetting->uLLast << "\n"
			<< " - vLFirst: first horizontal low resolution pixel to be processed (before possible correction): " << dSetting->vLFirst << "\n"
			<< " - vLLast:  last horizontal low resolution pixel to be processed (before possible correction):  " << dSetting->vLLast << "\n"
			<< " - delete_tmp_patch_folders:                                                                    " << dSetting->delete_tmp_patch_folders << "\n"
			<< " - dir_tmp_patches_additional_num: number of additional directories containing relevant tmp patches: " << dSetting->dir_tmp_patches_additional_num << "\n"
			<< " - continue unfinished reconstruction?:                                                         " << dSetting->contUnfinishedRec << "\n"
			<< "\n"
			<< "*========================================*" << "\n"
			<< "* image fusion parameters                *" << "\n"
			<< "*========================================*" << "\n"
			<< " - fusion method:                           ";
	switch (fSetting->fMethod) {
		case SparseFI: {
			this->file << "SparseFI" << "\n";
			break;
		}
		case JSparseFI: {
			this->file << "J-SparseFI" << "\n";
			break;
		}
		case JSparseFIHM: {
			this->file << "J-SparseFI-HM" << "\n";
			break;
		}
		case LeastSquares: {
			this->file << "Least Squares" << "\n";
			break;
		}
		default: {
			this->file << "unknown" << "\n";
			break;
		}
	}
	this->file << " - two-step estimation:                     " << fSetting->two_step_estimation << "\n";
	this->file << " - sparse reconstruction solver:            ";
	switch (sSetting->solver) {
	case JPFISTA: {
//		cout << "J-P-FISTA" << endl;
		this->file << "J-P-FISTA" << "\n";
		break;
	}
	default: {
//		cout << "unknown" << endl;
		this->file << "unknown" << "\n";
		break;
	}
	}
	this->file
	        << "     -> solver settings:  - maxiter_out:                         "
			<< sSetting->maxiter_out << "\n"
			<< "                          - tol:                                 "
			<< sSetting->tol << "\n"
			<< "                          - matrix norm for dict. normalization: "
			<< fSetting->matrixNorm;
	switch(fSetting->matrixNorm){
	case 0: {
			this->file << " (Frobenius norm)\n";
			break;
		}
	case 1: {
		this->file << " (spectral norm)\n";
		break;
	}
	}
 this->file << " - Least Sq. post minimiz. for Z via CGLS:  " << fSetting->LQ_post_opt << "\n"
            << "     -> Patch-wise coefficient estimation parameters:\n"
            << "        - useNewMethodForCalculatingZ:  " << fSetting->useNewMethodForCalculatingZ << "\n"
            << "        - useSimulatedImXforDictLearn:  " << fSetting->useSimulatedImXforDictLearn << "\n"
            << "                       - lambda_X_ABC:  " << fSetting->lambdaX_ABC << "\n"
            << "                       - lambda_Y_ABC:  " << fSetting->lambdaY_ABC << "\n"
            << "                       - lambda_Z_ABC:  " << fSetting->lambdaZ_ABC << "\n"
            << "        - lambdaZ_ABC_in_1st_iter:      " << fSetting->lambdaZ_ABC_in_1st_iter << "\n"
 	 	 	<< "     -> full image CGLS settings: ACTIVE=" << fSetting->LQ_post_opt_im << "\n"
            << "                       - lambda_X_im:      " << fSetting->lambdaX_im << "\n"
            << "                       - lambda_Y_im:      " << fSetting->lambdaY_im << "\n"
            << "                       - maxiter CGLS_im:  " << sSetting->maxiter_CGLS_im << "\n"
            << "                       - tol r CGLS_im:    " << sSetting->tol_r_CGLS_im << "\n"
            << "                       - subspace_transformation_type:  " << fSetting->subspace_transform_type << "\n"
            << "                       - subspace_dim:                  " << fSetting->subspace_dim << "\n"
            << "     -> Patch-wise CGLS settings (former 'Eq.3'): ACTIVE=" << fSetting->LQ_post_opt << "\n"
			<< "                               - lambda_X:        " << fSetting->lambdaX << "\n"
			<< "                               - lambda_Y:        " << fSetting->lambdaY << "\n"
			<< "                               - maxiter CGLS:    " << sSetting->maxiter_CGLS << "\n"
			<< "                               - tol r CGLS:      " << sSetting->tol_r_CGLS << "\n"
			<< "                               - Alphas fixed:    " << sSetting->fix_Alpha << "\n"
			<< "                               - delta_m fixed:   " << sSetting->fix_delta_m << "\n";

	this->file
			<< " - image normalization:                     "
			<< fSetting->nrmlIm << "\n"
			<< " - dictionary normalization:                "
			<< fSetting->nrmlDicts << "\n"
			<< " - mean subtraction:                        "
			<< fSetting->substrMean << "\n"
			<< " - patch size (LR, one side):               "
			<< fSetting->patchsize << "\n"
			<< " - overlap (LR, one side):                  "
			<< fSetting->overlap << "\n"
			<< " - lambda:                                  "
			<< fSetting->lambda << "\n"
			<< " - NDP (before possible correction):        "
			<< fSetting->NDP << "\n"
			<< " - Nc (number of bundled channels - before possible correction):      "
			<< fSetting->Nc << "\n"
			<< " - No (number of overlapping channels - before possible correction):  "
			<< fSetting->No << "\n"
			<< " - Nc_max (maximum size of spectral group above groups will" << "\n"
			<< "   be double-checked and perhaps split into subgoups (int)):          "
			<< fSetting->Nc_max << "\n"
			<< " - CC_min (minimum cross-correlation within spectral goups (double)): "
			<< fSetting->CC_min << "\n"
			<< " - tol_SRF (tolerance for decision matrix): "
			<< fSetting->tol_SRF << "\n"
			<< " - dictionary selection method ID:          "
			<< fSetting->dictselect << " (";
	switch(fSetting->dictselect){
			case 0: this->file << "Dictionary contains ONLY the current patch -> (NDP=1) & Alpha is calculated by least squares"; break;
			case 1: this->file << "Nearest Neighbors"; break;
			case 2: this->file << "PanLR norm"; break;
			case 3: this->file << "SRF approximate PanLR norm"; break;
			case 4: this->file << "PanHR norm (POSITIVE PanHR correlation)"; break;
			case 5: this->file << "PanLR-PanHR joint ranking"; break;
			case 6: this->file << "ABSOLUTE PanHR correlation"; break;
			case 7: this->file << "PanHR anti-correlation, including current patch as first atom"; break;
			case 8: this->file << "Random, including current patch as first atom"; break;
			case 9: this->file << "PanHR self uncorrelated basis approximation, including current patch as first atom"; break;
			default: this->file << "\n\n ERROR: program terminated unsuccessfully due an UNKNOWN DICTIONARY SELECTION FLAG!" << "\n";
										 cerr << "\n\n ERROR: program terminated unsuccessfully due an UNKNOWN DICTIONARY SELECTION FLAG!" << "\n";
										 exit(2);
										 break;
	}
	this->file
			<< ")" << "\n"
			<< " - evaluate:                                " << fSetting->evaluate  << "\n"
			<< " - evaluate_ImZ_init:                       " << fSetting->evaluate_ImZ_init  << "\n"
			<< "\n"
			<< "\n"
			<< "*========================================*" << "\n"
			<< "* MPI parallelization parameters         *" << "\n"
			<< "*========================================*" << "\n"
			<< " - number of processes per group:                               " << pSetting->numProcGrp << "\n"
			<< " - number of processes per patch (before possible correction):  " << pSetting->numProcPerPatch << "\n"
			<< " - total number of processes:                                   " << my_processes << "\n"
			<< " - work stealing turns (parallelization strategy):              " << pSetting->workStealingTurns << "\n"
			<< " - store rec. patches temporarily on hard drive?                " << pSetting->store_patches_tmp_on_drive << "\n"
			<< " - max number of processors used for writing final image        " << pSetting->parWrNumProc << "\n"
			<< "\n"
			<< "*========================================*" << "\n"
			<< "* assessment parameters                  *" << "\n"
			<< "*========================================*" << "\n"
			<< " - reference data available?:               " << fSetting->ImZ_ref_avlbl << "\n"
			<< " - evaluation precision (number of digits): " << oSetting->prec << "\n"
			<< "\n"
			<< "*========================================*" << "\n"
			<< "* output settings                        *" << "\n"
			<< "*========================================*" << "\n"
			<< " - write image file?:                               " << oSetting->writeImageFile << "\n"
			<< " - save coefficient vectors?:                       " << oSetting->saveAlphas << "\n"
			<< " - if so, which coefficient vectors?: first patch:  " << oSetting->pFirstAlpha << "\n"
			<< " - if so, which coefficient vectors?: last patch:   " << oSetting->pLastAlpha << "\n"
			<< " - save local dictionary coordinates?:              " << oSetting->saveDicts << "\n"
			<< " - if so, for which patches?: first patch:          " << oSetting->pFirstDict << "\n"
			<< " - if so, for which patches?: last patch:           " << oSetting->pLastDict << "\n"
			<< " - evaluation precision (number of digits):         " << oSetting->prec << "\n"
			<< "\n";

	this->file << "\n"
			<< "###########################################################"
			<< "\n"
			<< "##                                                       ##"
			<< "\n"
			<< "##                    load input data                    ##"
			<< "\n"
			<< "##                                                       ##"
			<< "\n"
			<< "###########################################################"
			<< "\n" << "\n";

//			 << "glPrms->usedPatches = " << glPrms->usedPatches << "\n"
	//			 << "glPrms->NP = " << glPrms->NP << "\n"
	//			 << "glPrms->timeMainLoop = " << glPrms->timeMainLoop << "\n\n\n\n";

	this->file.close();











	//  Write to terminal:
	//cout << "\n"<< "fSetting->fMethod = " << fSetting->fMethod << "\n";
	cout << "\n"<< "#========================================================================#\n"
				<< "#                                                                        #\n";
		switch (fSetting->fMethod) {
			case SparseFI: {
				cout    << "#                                 SparseFI                               #\n"
						<< "#             [ Sparse Fusion of Images ] for pansharpening              #\n";
				break;
			}
			case JSparseFI: {
				cout	<< "#                                J-SparseFI                              #\n"
						<< "#         [ Jointly Sparse Fusion of Images ] for pansharpening          #\n";
				break;
			}
			case GroupedJSparseFI: {
				cout	<< "#                              GroupedJSparseFI                          #\n"
						<< "#         [ Jointly Sparse Fusion of Images ] for pansharpening          #\n";
				break;
			}
			case JSparseFIHM: {
				cout    << "#                               J-SparseFI-HM                            #\n"
						<< "#   [ Jointly Sparse Fusion of Hyperspectral and Multispectral Images ]  #\n";
				break;
			}
			case LeastSquares: {
				cout    << "#                             Least Squares                              #\n"
						<< "#                - Dictionary contains only current patch -              #\n";
						break;
			}
		}
		cout    << "#                                                                        #\n"
				<< "#========================================================================#\n\n"
				<< "starting date: " << datePretty << " " << timePretty
				<< "\n\n\n"
				<< "###########################################################\n"
				<< "##                                                       ##\n"
				<< "##                      job details                      ##\n"
				<< "##                                                       ##\n"
				<< "###########################################################"
				<< "\n" << "\n"
				<< " - job name:  " << dSetting->jobName << "\n"
				<< " - job ID:    " << dSetting->jobID << "\n"
				<< "\n" << "\n"
				<< "###########################################################\n"
				<< "##                                                       ##\n"
				<< "##                     user settings                     ##\n"
				<< "##                                                       ##\n"
				<< "###########################################################"
				<< "\n" << "\n"
				<< "*========================================*" << "\n"
				<< "* paths                                  *" << "\n"
				<< "*========================================*" << "\n"
				<< " - directory input data:      " << paths->dir_in << "\n"
				<< " - file name ImX:             " << paths->fname_ImX << "\n"
//				<< " - file name ImX_LR:          " << paths->fname_ImX_LR << "\n"
				<< " - file name ImY:             " << paths->fname_ImY << "\n"
				<< " - file name ImZ (reference): " << paths->fname_ImZ_ref << "\n"
				<< " - file name ImZ_init:        " << paths->fname_ImZ_init << "\n"
				<< " - directory output data:     " << paths->dir_out << "\n"
				<< " - file name ImZ:             " << paths->fname_ImZ << "\n"
				<< " - ### use_estimated_SRFs:                " << fSetting->use_estimated_SRFs << "\n"
				<< " - file name SRF:                         " << paths->fname_SRF << "\n"
				<< " - file name SRF (for spectral grouping): " << paths->fname_SRF_for_Spectral_Grouping << "\n"
				<< " - ImZ_init_type = " << fSetting->ImZ_init_type << "(";
		switch(fSetting->ImZ_init_type){
			case 0:{
				cout << "lambdaZ_ABZ=0 in 1st iter)\n";
				break;
			}
			case 1:{
				cout << "upsampled and bilinearly interpolated low resolution image ImY)\n";
				break;
			}
			case 2:{
				cout << "rec. result of another algorithm)\n";
				break;
			}
		}
		cout	<< "\n"
				<< "*========================================*" << "\n"
				<< "* image data I/O parameters              *" << "\n"
				<< "*========================================*" << "\n"
				<< " - write image file?:                                         " << oSetting->writeImageFile << "\n"
				<< " - first Y band to be processed (before possible correction): " << dSetting->chBundleFirst << "\n"
				<< " - last Y band to be processed (before possible correction): " << dSetting->chBundleLast << "\n"
				<< " - uLFirst: first vertical low resolution pixel to be processed (before possible correction):   " << dSetting->uLFirst << "\n"
				<< " - uLLast:  last vertical low resolution pixel to be processed (before possible correction):    " << dSetting->uLLast << "\n"
				<< " - vLFirst: first horizontal low resolution pixel to be processed (before possible correction): " << dSetting->vLFirst << "\n"
				<< " - vLLast:  last horizontal low resolution pixel to be processed (before possible correction):  " << dSetting->vLLast << "\n"
				<< "\n"
				<< "*========================================*" << "\n"
				<< "* image fusion parameters                *" << "\n"
				<< "*========================================*" << "\n"
				<< " - fusion method:                           ";
		switch (fSetting->fMethod) {
		case SparseFI: {
			cout << "SparseFI" << "\n";
			break;
		}
		case JSparseFI: {
			cout << "J-SparseFI" << "\n";
			break;
		}
		case JSparseFIHM: {
			cout << "J-SparseFI-HM" << "\n";
			break;
		}
		case LeastSquares: {
			cout << "Least Squares" << "\n";
			break;
		}
		default: {
			cout << "unknown" << "\n";
			break;
		}
		}
		cout    << " - two-step estimation:                     " << fSetting->two_step_estimation << "\n"
		        << " - sparse reconstruction solver:            ";
		switch (sSetting->solver) {
		case JPFISTA: {
			cout << "J-P-FISTA" << endl;
			break;
		}
		default: {
			cout << "unknown" << endl;
			break;
		}
		}
		cout    << "     -> solver settings:  - maxiter_out:    "
				<< sSetting->maxiter_out << "\n"
				<< "                          - tol:            "
				<< sSetting->tol << "\n"
				<< " - image normalization:                     "
				<< fSetting->nrmlIm << "\n"
				<< " - dictionary normalization:                "
				<< fSetting->nrmlDicts << "\n"
				<< " - mean subtraction:                        "
				<< fSetting->substrMean << "\n"
				<< " - patch size (LR, one side):               "
				<< fSetting->patchsize << "\n"
				<< " - overlap (LR, one side):                  "
				<< fSetting->overlap << "\n"
				<< " - lambda:                                  "
				<< fSetting->lambda << "\n"
				<< " - NDP: (before possible correction)        "
				<< fSetting->NDP
				<< "\n"
				<< " - Nc (number of bundled channels - before possible correction):        "
				<< fSetting->Nc << "\n"
				<< " - No (number of overlapping channels - before possible correction):    "
				<< fSetting->No << "\n"
				<< " - Nc_max (maximum size of spectral group above groups will" << "\n"
				<< "   be double-checked and perhaps split into subgoups (int)):            "
				<< fSetting->Nc_max << "\n"
				<< " - CC_min (minimum cross-correlation within spectral goups (double)):   "
				<< fSetting->CC_min
				<< "\n"
				<< " - tol_SRF (tolerance for decision matrix): "
				<< fSetting->tol_SRF << "\n"
				<< " - evaluate:                                " << fSetting->evaluate  << "\n"
				<< " - evaluate_ImZ_init:                       " << fSetting->evaluate_ImZ_init  << "\n"
				<< " - dictionary selection method ID: "
				<< fSetting->dictselect << "\n" << "\n"
//	            << "     -> Patch-wise coefficient estimation parameters:\n"
//	            << "        - useNewMethodForCalculatingZ:      " << fSetting->useNewMethodForCalculatingZ << "\n"
//	            << "        - useSimulatedImXforDictLearn:      " << fSetting->useSimulatedImXforDictLearn << "\n"
//	            << "                       - lambda_X_ABC:      " << fSetting->lambdaX_ABC << "\n"
//	            << "                       - lambda_Y_ABC:      " << fSetting->lambdaY_ABC << "\n"
//	            << "                       - lambda_Z_ABC:      " << fSetting->lambdaZ_ABC << "\n"
//	            << "        - lambda_Z_ABC=0 in first iteration:" << fSetting->lambdaZ_ABC_in_1st_iter << "\n"
//	 	 	 	<< "     -> full image CGLS settings:\n"
//	            << "                       - lambda_X_im:       " << fSetting->lambdaX_im << "\n"
//	            << "                       - lambda_Y_im:       " << fSetting->lambdaY_im << "\n"
//	            << "                       - maxiter CGLS_im:   " << sSetting->maxiter_CGLS_im << "\n"
//	            << "                       - tol r CGLS_im:     " << sSetting->tol_r_CGLS_im << "\n"
//	            << "     -> Patch-wise CGLS settings (former 'Eq.3'):\n"
//				<< "                       - lambda_X:          " << fSetting->lambdaX << "\n"
//				<< "                       - lambda_Y:          " << fSetting->lambdaY << "\n"
//				<< "                       - maxiter CGLS:      " << sSetting->maxiter_CGLS << "\n"
//				<< "                       - tol r CGLS:        " << sSetting->tol_r_CGLS << "\n"
//				<< "                       - Alphas fixed:      " << sSetting->fix_Alpha << "\n"
//				<< "                       - delta_m fixed:     " << sSetting->fix_delta_m << "\n"
				<< "     -> Patch-wise coefficient estimation parameters:\n"
				<< "        - useNewMethodForCalculatingZ:  " << fSetting->useNewMethodForCalculatingZ << "\n"
				<< "        - useSimulatedImXforDictLearn:  " << fSetting->useSimulatedImXforDictLearn << "\n"
				<< "                       - lambda_X_ABC:  " << fSetting->lambdaX_ABC << "\n"
				<< "                       - lambda_Y_ABC:  " << fSetting->lambdaY_ABC << "\n"
				<< "                       - lambda_Z_ABC:  " << fSetting->lambdaZ_ABC << "\n"
				<< "        - lambdaZ_ABC_in_1st_iter:      " << fSetting->lambdaZ_ABC_in_1st_iter << "\n"
				<< "     -> full image CGLS settings: ACTIVE=" << fSetting->LQ_post_opt_im << "\n"
				<< "                       - lambda_X_im:      " << fSetting->lambdaX_im << "\n"
				<< "                       - lambda_Y_im:      " << fSetting->lambdaY_im << "\n"
				<< "                       - maxiter CGLS_im:  " << sSetting->maxiter_CGLS_im << "\n"
				<< "                       - tol r CGLS_im:    " << sSetting->tol_r_CGLS_im << "\n"
				<< "     -> Patch-wise CGLS settings (former 'Eq.3'): ACTIVE=" << fSetting->LQ_post_opt << "\n"
				<< "                               - lambda_X:        " << fSetting->lambdaX << "\n"
				<< "                               - lambda_Y:        " << fSetting->lambdaY << "\n"
				<< "                               - maxiter CGLS:    " << sSetting->maxiter_CGLS << "\n"
				<< "                               - tol r CGLS:      " << sSetting->tol_r_CGLS << "\n"
				<< "                               - Alphas fixed:    " << sSetting->fix_Alpha << "\n"
				<< "                               - delta_m fixed:   " << sSetting->fix_delta_m << "\n"
//
				<< "*========================================*" << "\n"
				<< "* MPI parallelization parameters         *" << "\n"
				<< "*========================================*" << "\n"
				<< " - number of processes per group:  " << pSetting->numProcGrp
				<< "\n" << " - total number of processes:      " << my_processes
				<< "\n" << " - work stealing turns (patch parallelization strategy): "
				<< pSetting->workStealingTurns << "\n" << "\n"
				<< "*========================================*" << "\n"
				<< "* assessment parameters                  *" << "\n"
				<< "*========================================*" << "\n"
				<< " - reference data available?:               " << fSetting->ImZ_ref_avlbl << "\n"
				<< " - evaluation precision (number of digits): " << oSetting->prec << "\n"
				<< "\n\n"
				<< "###########################################################\n"
				<< "##                                                       ##\n"
				<< "##                    load input data                    ##\n"
				<< "##                                                       ##\n"
				<< "###########################################################"
				<< "\n" << "\n";
}



void SpEOReport::finalize(SpEOGlobalParams *glPrms) {

	double totalTime = MPI_Wtime() - this->curTimeSec;
	glPrms->timeTotal = totalTime;
	this->curTime = get_current_time();
//	string mydate = this->curTime.substr(0, 6);
//	string mytime = this->curTime.substr(7, 6);
	this->file.open(this->fileName.c_str(),
			fstream::in | fstream::out | fstream::app);

	string datePretty = "20" + this->curTime.substr(0, 2) + "-"
			+ this->curTime.substr(2, 2) + "-" + this->curTime.substr(4, 2);
	string timePretty = this->curTime.substr(7, 2) + ":"
			+ this->curTime.substr(9, 2) + ":" + this->curTime.substr(11, 2);

	this->file << "\n"
			<< "_______________________________________________________________________________\n\n"
			<< "ending date: " << datePretty << " " << timePretty
			<< "\n\n"
			<< "  -> Method needed                           " << totalTime << " seconds in total which of \n"
			<< "         - the main loop took                " << glPrms->timeMainLoop << " seconds,\n"
			<< "         - total dictionary selection took   " << glPrms->timeDictSelect << " seconds,\n"
			<< "         - average dictionary selection took " << glPrms->timeDictSelect_avg << " seconds,\n"
			<< "         - the data-to-file writing took " << glPrms->timeFileWrite << " seconds."
			<< "\n\n";

	this->file
			<< "###############################################################################\n"
			<< "###                               __         _                              ###\n"
			<< "###                              |__  |\\ |  | \\                             ###\n"
			<< "###                              |__  | \\|  |_/                             ###\n"
			<< "###                                                                         ###\n"
			<< "###############################################################################\n";

	this->file.close();
	chmod(this->fileName.c_str(), S_IRWXU | S_IRWXG | S_IRWXO);


	cout    << "\n"
			<< "_______________________________________________________________________________\n\n"
			<< "ending date: " << datePretty << " " << timePretty
			<< "\n\n"
			<< "  -> Method needed                           " << totalTime << " seconds in total which of \n"
			<< "         - the main loop took                " << glPrms->timeMainLoop << " seconds,\n"
			<< "         - total dictionary selection took   " << glPrms->timeDictSelect << " seconds,\n"
			<< "         - average dictionary selection took " << glPrms->timeDictSelect_avg << " seconds,\n"
			<< "         - the data-to-file writing took " << glPrms->timeFileWrite << " seconds."
			<< "\n"
			<< "\n"
			<< "###############################################################################\n"
			<< "###                               __         _                              ###\n"
			<< "###                              |__  |\\ |  | \\                             ###\n"
			<< "###                              |__  | \\|  |_/                             ###\n"
			<< "###                                                                         ###\n"
			<< "###############################################################################\n\n";

}




//void JSparseFI_reportAddEval(SpEOAssessmentMetrics  *assessm_res_HR, SpEOFusionSetting *fSetting, SpEOGlobalParams *glPrms, int type, ofstream& reportFile, string reportFileName)
void SpEOReport::addEvaluation(SpEOAssessmentMetrics *assessm_res_HR,
								SpEODataIOSetting *dSetting,
								SpEOFusionSetting *fSetting,
                               SpEOOutputSetting *oSetting,
                               SpEOGlobalParams *glPrms,
                               int type) {
	SpEOAssessmentMetricsStr *assessm_res_HR_Str = new SpEOAssessmentMetricsStr;
	SpEOAssessmentMetricsStr_Sep *assessm_res_HR_StrSep = new SpEOAssessmentMetricsStr_Sep[glPrms->NChY];
//	for(int i=0; i<fSetting->maxIter; i++){

	int assessNum = NUM_OF_ASSESSMENT, temp;
	int allDigits[assessNum][glPrms->NChY + 1];
	int maxDigits[assessNum];

	assessm_res_HR_Str->RMSE_str    = double2str(assessm_res_HR->RMSE_mean,oSetting->prec);
	assessm_res_HR_Str->PSNR_str    = double2str(assessm_res_HR->PSNR_mean,oSetting->prec);
	assessm_res_HR_Str->CC_str      = double2str(assessm_res_HR->CC_mean,oSetting->prec);
	assessm_res_HR_Str->ERGAS_str   = double2str(assessm_res_HR->ERGAS_mean,oSetting->prec);
	assessm_res_HR_Str->UIQI_str    = double2str(assessm_res_HR->UIQI_mean,oSetting->prec);
	assessm_res_HR_Str->DD_str      = double2str(assessm_res_HR->DD_mean,oSetting->prec);
	assessm_res_HR_Str->SAM_str     = double2str(assessm_res_HR->SAM,oSetting->prec);
	assessm_res_HR_Str->DLambda_str = double2str(assessm_res_HR->DLambda_mean,oSetting->prec);
	assessm_res_HR_Str->AG_orig_str = double2str(assessm_res_HR->AG_orig_mean,oSetting->prec);
	assessm_res_HR_Str->AG_rec_str  = double2str(assessm_res_HR->AG_rec_mean,oSetting->prec);

	allDigits[0][0] = assessm_res_HR_Str->RMSE_str.length();
	allDigits[1][0] = assessm_res_HR_Str->CC_str.length();
	allDigits[2][0] = assessm_res_HR_Str->ERGAS_str.length();
	allDigits[3][0] = assessm_res_HR_Str->UIQI_str.length();
	allDigits[4][0] = assessm_res_HR_Str->DD_str.length();
	allDigits[5][0] = assessm_res_HR_Str->SAM_str.length();
	allDigits[6][0] = assessm_res_HR_Str->DLambda_str.length();
	allDigits[7][0] = assessm_res_HR_Str->AG_orig_str.length();
	allDigits[8][0] = assessm_res_HR_Str->AG_rec_str.length();
	allDigits[9][0] = assessm_res_HR_Str->PSNR_str.length();

	for (int i = 0; i < glPrms->NChY; i++) {

		assessm_res_HR_StrSep[i].RMSE_sep_str = double2str(
				assessm_res_HR->RMSE_sep[i], oSetting->prec);
		assessm_res_HR_StrSep[i].PSNR_sep_str = double2str(
				assessm_res_HR->PSNR_sep[i], oSetting->prec);
		assessm_res_HR_StrSep[i].CC_sep_str = double2str(
				assessm_res_HR->CC_sep[i], oSetting->prec);
		assessm_res_HR_StrSep[i].ERGAS_sep_str = double2str(
				assessm_res_HR->ERGAS_sep[i], oSetting->prec);
		assessm_res_HR_StrSep[i].UIQI_sep_str = double2str(
				assessm_res_HR->UIQI_sep[i], oSetting->prec);
		assessm_res_HR_StrSep[i].DD_sep_str = double2str(assessm_res_HR->DD_sep[i], oSetting->prec);
		assessm_res_HR_StrSep[i].AG_orig_sep_str = double2str(
				assessm_res_HR->AG_orig_sep[i], oSetting->prec);
		assessm_res_HR_StrSep[i].AG_rec_sep_str = double2str(
				assessm_res_HR->AG_rec_sep[i], oSetting->prec);

		allDigits[0][i + 1] = assessm_res_HR_StrSep[i].RMSE_sep_str.length();
		allDigits[1][i + 1] = assessm_res_HR_StrSep[i].CC_sep_str.length();
		allDigits[2][i + 1] = assessm_res_HR_StrSep[i].ERGAS_sep_str.length();
		allDigits[3][i + 1] = assessm_res_HR_StrSep[i].UIQI_sep_str.length();
		allDigits[4][i + 1] = assessm_res_HR_StrSep[i].DD_sep_str.length();
		allDigits[5][i + 1] = 0;
		allDigits[6][i + 1] = 0;
		allDigits[7][i + 1] = assessm_res_HR_StrSep[i].AG_orig_sep_str.length();
		allDigits[8][i + 1] = assessm_res_HR_StrSep[i].AG_rec_sep_str.length();
		allDigits[9][i + 1] = assessm_res_HR_StrSep[i].PSNR_sep_str.length();

	}
	temp = 0;
	for (int i = 0; i < assessNum; i++) {
		for (int j = 0; j <= glPrms->NChY; j++) {
			if (allDigits[i][j] > temp)
				temp = allDigits[i][j];
		}
		maxDigits[i] = temp;
		temp = 0;
	}
	//createTable(assessm_res_HR_Str, assessm_res_HR_StrSep, fSetting, glPrms, maxDigits, type, this->file, this->fileName);
	createTable(assessm_res_HR_Str, assessm_res_HR_StrSep, dSetting, fSetting, glPrms,
			maxDigits, type, this);
//	}
	delete assessm_res_HR_Str;
	delete[] assessm_res_HR_StrSep;

	return;
}

void SpEOReport::addGlobalParams(SpEOGlobalParams *glPrms, SpEOFusionSetting *fSetting) {

	this->file.open(this->fileName.c_str(),
			fstream::in | fstream::out | fstream::app);
	this->file << "\n"
			<< "###########################################################"
			<< "\n"
			<< "##                                                       ##"
			<< "\n"
			<< "##                    global parameters                  ##"
			<< "\n"
			<< "##                                                       ##"
			<< "\n"
			<< "###########################################################"
			<< "\n" << "\n"
			<< " - NDP after possible corrected to min(NP,NDP_init):                  "
			<< fSetting->NDP << "\n"
			<< " - Nc after possible correction to min(Nc,NChY):                      "
			<< fSetting->Nc << "\n"
			<< " - No after possible correction to either min(No,Nc) or 0):           "
			<< fSetting->No << "\n"
			<< " - Nc_max (maximum size of spectral group above groups will" << "\n"
			<< "   be double-checked and perhaps split into subgoups (int)):          "
			<< fSetting->Nc_max << "\n"
			<< " - CC_min (minimum cross-correlation within spectral goups (double)): "
			<< fSetting->CC_min << "\n"
			<< " - NChX (number of bands in image ImX):            "
			<< glPrms->NChX << "\n"
			<< " - NChY (number of bands in image ImY):            "
			<< glPrms->NChY << "\n"
			<< " - NChZ (number of bands in image ImZ):            "
			<< glPrms->NChZ << "\n"
			<< " - fDS (down-sampling factor):                     "
			<< glPrms->fDS << "\n"
			<< " - NP (total number of patches per band):                          "
			<< glPrms->NP << "\n"
			<< " - NPU (total number of patches in horizontal direction per band):  "
			<< glPrms->NPU << "\n"
			<< " - NPV (total number of patches in vertical direction per band):    "
			<< glPrms->NPV << "\n"
			<< " - NP_sub (total number of actually processed patches per band):                     "
			<< glPrms->NP_sub << "\n"
			<< " - NPU_sub (number of actually processed patches in horizontal direction per band):  "
			<< glPrms->NPU_sub << "\n"
			<< " - NPV_sub (number of actually processed patches in vertical direction per band):    "
			<< glPrms->NPV_sub << "\n"
			<< " - number of bundles of jointly processed bands:   "
			<< glPrms->Ng << "\n"
			<< " - number of problems per patch (# nontrivial entries in decision matrix):   "
			<< glPrms->numProbPerPatch << "\n"
			<< " - decision matrix C (transposed): \n"
			<< glPrms->decMat_C.transpose() << "\n" << " - N_m & k_m: \n";
	for (int iChX = 0; iChX < glPrms->NChX; iChX++) {
		this->file << "       N_m[" << iChX << "] = " << glPrms->Nm[iChX]
				<< ", \t k_m[" << iChX << "]:" << "  ";
		for (int i = 0; i < glPrms->Nm[iChX]; i++) {
			this->file << glPrms->km[iChX][i] << "\t";
		}
		this->file << "\n";
	}
	this->file << " - N_m_2 & k_m_2: \n";
	// cout Nm2 & km2
	for (int iG = 0; iG < glPrms->Ng; iG++) {
		this->file << "       N_m_2[" << iG << "] = " << glPrms->Nm2[iG]
				<< ", \t k_m_2[" << iG << "]:" << "  ";
		for (int iChX = 0; iChX < glPrms->Nm2[iG]; iChX++) {
			this->file << glPrms->km2[iG][iChX] << "\t";
		}
		this->file << "\n";
	}
	this->file << "\n";

	this->file.close();
}

string double2str(double num, int prec) {
	stringstream ss;
	ss << std::setprecision(prec) << num;
	return ss.str();
}

string int2str(int number){
   stringstream ss;
   ss << number;
   return ss.str();
}

string createStr(int number, string Str) {
	string strArr;
	for (int i = 0; i < number; i++)
		strArr += Str;
	return strArr;
}


void createTable(SpEOAssessmentMetricsStr *assessmStr,
		SpEOAssessmentMetricsStr_Sep *assessmStrSep,
		SpEODataIOSetting *dSetting, SpEOFusionSetting *fSetting, SpEOGlobalParams *glPrms, int *maxDigits,
		int type, SpEOReport *report) {

	int totalElNum = 0;
	int resNum =  2 // "# "
			    + 4 // "band"
			    + 4 + 9 * 3 // " || "  +  9* " | "
			    + 2; // " #"
	int resNumFile = 2 // "# "
	+ 4 // "band"
			+ 4 + 3 * 3 // " || "  +  3* " | "
			+ 2; // " #"
	int resNumFile2 = 2 // "# "
		+ 4 // "band"
				+ 4 + 5 * 3 // " || "  +  5* " | "
				+ 2; // " #"

	int firstColLen = 4, fileAdj = 0;

	report->file.open(report->fileName.c_str(),
			fstream::in | fstream::out | fstream::app);


	report->file << "\n"
			<< "###########################################################" << "\n"
			<< "##                                                       ##" << "\n"
			<< "##                  Evaluation results                   ##" << "\n"
			<< "##                                                       ##" << "\n"
			<< "###########################################################"
			<< "\n" << "\n";

	string temStr;
	string outLinSym = "#";
	string lineBeg = "#";
	string lineEnd = "#";
	string delimiter1 = "|";
	string delimiter2 = "||";
	string midLinSym = "-";

	// set lower bound to each entry in maxDigits
	if (maxDigits[0] < 4)
		maxDigits[0] = 4; // because length("RMSE")     = 4
	if (maxDigits[1] < 2)
		maxDigits[1] = 2; // because length("CC")       = 2
	if (maxDigits[2] < 5)
		maxDigits[2] = 5; // because length("ERGAS")    = 5
	if (maxDigits[3] < 4)
		maxDigits[3] = 4; // because length("UIQI")     = 4
	if (maxDigits[4] < 2)
		maxDigits[4] = 2; // because length("DD")       = 2
	if (maxDigits[5] < 3)
		maxDigits[5] = 3; // because length("SAM")      = 3
	if (maxDigits[6] < 8)
		maxDigits[6] = 8; // because length("D_lambda") = 8
	if (maxDigits[7] < 7)
		maxDigits[7] = 7; // because length("AG_orig")  = 7
	if (maxDigits[8] < 6)
		maxDigits[8] = 6; // because length("AG_rec")   = 6
	if (maxDigits[9] < 4)
		maxDigits[9] = 4; // because length("PSNR")     = 4

	for (int i = 0; i < NUM_OF_ASSESSMENT; i++) {
		totalElNum += maxDigits[i];
	}

	int lineLen = totalElNum + resNum;

	if (type == LR_ASSESSMENT)
		cout << "Low Resolution:" << endl;
	else if (type == HR_ASSESSMENT)
		cout << endl << "High Resolution:" << endl;

	cout << createStr(lineLen, outLinSym) << endl;

	cout << lineBeg << " " << "band" << " " << delimiter2 << " "
			<< "RMSE"     << createStr(maxDigits[0] - 4, " ") << " " << delimiter1 << " "
			<< "CC"       << createStr(maxDigits[1] - 2, " ") << " " << delimiter1 << " "
			<< "ERGAS"    << createStr(maxDigits[2] - 5, " ") << " " << delimiter1 << " "
			<< "UIQI"     << createStr(maxDigits[3] - 4, " ") << " " << delimiter1 << " "
			<< "DD"       << createStr(maxDigits[4] - 2, " ") << " " << delimiter1 << " "
			<< "SAM"      << createStr(maxDigits[5] - 3, " ") << " " << delimiter1 << " "
			<< "D_lambda" << createStr(maxDigits[6] - 8, " ") << " " << delimiter1 << " "
			<< "AG_orig"  << createStr(maxDigits[7] - 7, " ") << " " << delimiter1 << " "
			<< "AG_rec"   << createStr(maxDigits[8] - 6, " ") << " " << delimiter1 << " "
	        << "PSNR"     << createStr(maxDigits[9] - 4, " ") << " " << lineEnd << endl;

	cout << lineBeg << createStr(firstColLen + 2, midLinSym) << delimiter2
			<< createStr(totalElNum + resNum - 10, midLinSym) << lineEnd
			<< endl;

	if (type == LR_ASSESSMENT) {
		report->file << "Low Resolution:\n";
	} else if (type == HR_ASSESSMENT) {
		report->file << "High Resolution:\n";
	}

	fileAdj = maxDigits[0] + maxDigits[1] + maxDigits[2] + maxDigits[3];
	report->file << createStr(fileAdj + resNumFile, outLinSym) << "\n";

	report->file << lineBeg << " " << "band" << " " << delimiter2 << " "
			<< "RMSE"  << createStr(maxDigits[0] - 4, " ") << " " << delimiter1 << " "
			<< "CC"    << createStr(maxDigits[1] - 2, " ") << " " << delimiter1 << " "
			<< "ERGAS" << createStr(maxDigits[2] - 5, " ") << " " << delimiter1 << " "
			<< "UIQI"  << createStr(maxDigits[3] - 4, " ") << " " << lineEnd << " \n";

	report->file << lineBeg << createStr(firstColLen + 2, midLinSym)
			<< delimiter2 << createStr(fileAdj + resNumFile - 10, midLinSym)
			<< lineEnd << "\n";

	for (int j = 0; j < glPrms->NChY; j++) {

		string bandStr =
				static_cast<ostringstream*>(&(ostringstream() << j + 1))->str();
		cout << lineBeg << createStr(firstColLen - bandStr.length(), " ")
				<< j + 1 + dSetting->chBundleFirst << "  " << delimiter2 << " "
				<< assessmStrSep[j].RMSE_sep_str    << createStr( maxDigits[0] - assessmStrSep[j].RMSE_sep_str.length(),  " ") << " " << delimiter1 << " "
				<< assessmStrSep[j].CC_sep_str      << createStr( maxDigits[1] - assessmStrSep[j].CC_sep_str.length(),    " ") << " " << delimiter1 << " "
				<< assessmStrSep[j].ERGAS_sep_str   << createStr( maxDigits[2] - assessmStrSep[j].ERGAS_sep_str.length(), " ") << " " << delimiter1 << " "
				<< assessmStrSep[j].UIQI_sep_str    << createStr( maxDigits[3] - assessmStrSep[j].UIQI_sep_str.length(),  " ") << " " << delimiter1 << " "
				<< assessmStrSep[j].DD_sep_str      << createStr( maxDigits[4] - assessmStrSep[j].DD_sep_str.length(),    " ") << " " << delimiter1 << " "
				<< createStr(maxDigits[5], " ") << " " << delimiter1 << " "
				<< createStr(maxDigits[6], " ") << " " << delimiter1 << " "
				<< assessmStrSep[j].AG_orig_sep_str << createStr( maxDigits[7] - assessmStrSep[j].AG_orig_sep_str.length(), " ") << " " << delimiter1 << " "
				<< assessmStrSep[j].AG_rec_sep_str  << createStr( maxDigits[8] - assessmStrSep[j].AG_rec_sep_str.length(),  " ") << " " << delimiter1 << " "
		        << assessmStrSep[j].PSNR_sep_str    << createStr( maxDigits[9] - assessmStrSep[j].PSNR_sep_str.length(),    " ") << " " << lineEnd << endl;

		report->file << outLinSym
				<< createStr(firstColLen - bandStr.length(), " ") << j + 1 +dSetting->chBundleFirst	<< "  " << delimiter2 << " "
				<< assessmStrSep[j].RMSE_sep_str  << createStr( maxDigits[0] - assessmStrSep[j].RMSE_sep_str.length(),  " ") << " " << delimiter1 << " "
				<< assessmStrSep[j].CC_sep_str    << createStr( maxDigits[1] - assessmStrSep[j].CC_sep_str.length(),    " ") << " " << delimiter1 << " "
				<< assessmStrSep[j].ERGAS_sep_str << createStr( maxDigits[2] - assessmStrSep[j].ERGAS_sep_str.length(), " ") << " " << delimiter1 << " "
				<< assessmStrSep[j].UIQI_sep_str  << createStr( maxDigits[3] - assessmStrSep[j].UIQI_sep_str.length(),  " ") << " " << lineEnd << "\n";
	}
	cout << lineBeg << createStr(firstColLen + 2, midLinSym) << delimiter2
			<< createStr(totalElNum + resNum - 10, midLinSym) << lineEnd
			<< endl;
	report->file << lineBeg << createStr(firstColLen + 2, midLinSym)
			<< delimiter2 << createStr(fileAdj + resNumFile - 10, midLinSym)
			<< lineEnd << "\n";

	cout << lineBeg << " " << "mean" << " " << delimiter2 << " "
			<< assessmStr->RMSE_str    << createStr(maxDigits[0] - assessmStr->RMSE_str.length(), " ")    << " " << delimiter1 << " "
			<< assessmStr->CC_str      << createStr(maxDigits[1] - assessmStr->CC_str.length(), " ")      << " " << delimiter1 << " "
			<< assessmStr->ERGAS_str   << createStr(maxDigits[2] - assessmStr->ERGAS_str.length(), " ")   << " " << delimiter1 << " "
			<< assessmStr->UIQI_str    << createStr(maxDigits[3] - assessmStr->UIQI_str.length(), " ")    << " " << delimiter1 << " "
			<< assessmStr->DD_str      << createStr(maxDigits[4] - assessmStr->DD_str.length(), " ")      << " " << delimiter1 << " "
			<< assessmStr->SAM_str     << createStr(maxDigits[5] - assessmStr->SAM_str.length(), " ")     << " " << delimiter1 << " "
			<< assessmStr->DLambda_str << createStr(maxDigits[6] - assessmStr->DLambda_str.length(), " ") << " " << delimiter1 << " "
			<< assessmStr->AG_orig_str << createStr(maxDigits[7] - assessmStr->AG_orig_str.length(), " ") << " " << delimiter1 << " "
			<< assessmStr->AG_rec_str  << createStr(maxDigits[8] - assessmStr->AG_rec_str.length(), " ")  << " " << delimiter1 << " "
	        << assessmStr->PSNR_str    << createStr(maxDigits[9] - assessmStr->PSNR_str.length(), " ")    << " " << lineEnd << endl;
	cout << createStr(lineLen, outLinSym) << endl << endl;

	report->file << lineBeg << " " << "mean" << " " << delimiter2 << " "
			<< assessmStr->RMSE_str  << createStr(maxDigits[0] - assessmStr->RMSE_str.length(), " ")  << " " << delimiter1 << " "
			<< assessmStr->CC_str    << createStr(maxDigits[1] - assessmStr->CC_str.length(), " ")    << " " << delimiter1 << " "
			<< assessmStr->ERGAS_str << createStr(maxDigits[2] - assessmStr->ERGAS_str.length(), " ") << " " << delimiter1 << " "
			<< assessmStr->UIQI_str  << createStr(maxDigits[3] - assessmStr->UIQI_str.length(), " ")  << " " << lineEnd << "\n";
	report->file << createStr(fileAdj + resNumFile, outLinSym) << "\n";

	fileAdj = 0;
	fileAdj = maxDigits[4] + maxDigits[5] + maxDigits[6] + maxDigits[7] + maxDigits[8] + maxDigits[9];

	/***********First Row**************/
	report->file << createStr(fileAdj + resNumFile2, outLinSym) << "\n";
	report->file << lineBeg << " " << "band" << " " << delimiter2 << " "
			<< "DD"       << createStr(maxDigits[4] - 2, " ") << " " << delimiter1 << " "
			<< "SAM"      << createStr(maxDigits[5] - 3, " ") << " " << delimiter1 << " "
			<< "D_lambda" << createStr(maxDigits[6] - 8, " ") << " " << delimiter1 << " "
			<< "AG_orig"  << createStr(maxDigits[7] - 7, " ") << " " << delimiter1 << " "
			<< "AG_rec"   << createStr(maxDigits[8] - 6, " ") << " " << delimiter1 << " "
		    << "PSNR"     << createStr(maxDigits[9] - 4, " ") << " " << lineEnd << "\n";

	report->file << lineBeg << createStr(firstColLen + 2, midLinSym)
			<< delimiter2 << createStr(fileAdj + resNumFile2 - 10, midLinSym)
			<< lineEnd << "\n";

	/***********Band stream**************/
	for (int jj = 0; jj < glPrms->NChY; jj++) {
		string bandStr =
				static_cast<ostringstream*>(&(ostringstream() << jj + 1))->str();
		report->file << lineBeg
				<< createStr(firstColLen - bandStr.length(), " ") << jj + 1 + dSetting->chBundleFirst << "  " << delimiter2 << " "
				<< assessmStrSep[jj].DD_sep_str      << createStr( maxDigits[4] - assessmStrSep[jj].DD_sep_str.length(),      " ") << " " << delimiter1 << " "
				<< createStr(maxDigits[5], " ") << " " << delimiter1 << " "
				<< createStr(maxDigits[6], " ") << " " << delimiter1 << " "
				<< assessmStrSep[jj].AG_orig_sep_str << createStr( maxDigits[7] - assessmStrSep[jj].AG_orig_sep_str.length(), " ") << " " << delimiter1 << " "
				<< assessmStrSep[jj].AG_rec_sep_str  << createStr( maxDigits[8] - assessmStrSep[jj].AG_rec_sep_str.length(),  " ") << " " << delimiter1 << " "
				<< assessmStrSep[jj].PSNR_sep_str      << createStr( maxDigits[9] - assessmStrSep[jj].PSNR_sep_str.length(),  " ") << " " << lineEnd << "\n";
	}
	report->file << lineBeg << createStr(firstColLen + 2, midLinSym)
			<< delimiter2 << createStr(fileAdj + resNumFile2 - 10, midLinSym)
			<< lineEnd << "\n";

	/***********Mean Line**************/
	report->file << lineBeg << " " << "mean" << " " << delimiter2 << " "
			<< assessmStr->DD_str      << createStr(maxDigits[4] - assessmStr->DD_str.length(), " ")      << " " << delimiter1 << " "
			<< assessmStr->SAM_str     << createStr(maxDigits[5] - assessmStr->SAM_str.length(), " ")     << " " << delimiter1 << " "
			<< assessmStr->DLambda_str << createStr(maxDigits[6] - assessmStr->DLambda_str.length(), " ") << " " << delimiter1 << " "
			<< assessmStr->AG_orig_str << createStr(maxDigits[7] - assessmStr->AG_orig_str.length(), " ") << " " << delimiter1 << " "
			<< assessmStr->AG_rec_str  << createStr(maxDigits[8] - assessmStr->AG_rec_str.length(), " ")  << " " << delimiter1 << " "
			<< assessmStr->PSNR_str    << createStr(maxDigits[9] - assessmStr->PSNR_str.length(), " ")    << " " << lineEnd << "\n";

	report->file << createStr(fileAdj + resNumFile2, outLinSym) << "\n\n";
	report->file.close();
	return;
}

void CSVRow::readNextRow(std::istream& str, char delimiter) {
	std::string line;
	std::getline(str, line);
	this->line = line;
	std::stringstream lineStream(line);
	std::string cell;
	m_data.clear();
	while(std::getline(lineStream, cell, delimiter)){
//		if (line[0] == '#'){
//			cout << "line = " << line << ", cell = " << cell << endl;
//		}else{
			m_data.push_back(cell);
//		}
	}
}
// read and save .cvs file content to DOUBLE matrix
int read_CSV(SpEOMatrixD *Mat, const char *fname_CSV, char delimiter, int skipLns) {
	int my_rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
	std::ifstream CSV_file(fname_CSV);//std::ifstream CSV_file(fname_CSV.c_str());
	//cout << "CSV_file.is_open() = " << CSV_file.is_open() << endl;
	if (CSV_file.is_open()) {
		unsigned int CSVrow = 0;
		for(int i=0; i<skipLns; i++){
			CSV_file.ignore(500,'\n');
			//ignore the first 500 characters, or until first \n, whichever is met first
		}
		CSVRow m_row;
		m_row.readNextRow(CSV_file, delimiter);
		unsigned int CSVcol = m_row.size();
		while (CSV_file.good()) {
			if (CSVcol != m_row.size()) {
				cerr << endl << "[" << my_rank << "] ERROR while reading CSV in row number "
						<< CSVrow + 1 << ": All rows must have the same length!" << endl;
				exit(2);
			}
			CSVrow++;
			m_row.readNextRow(CSV_file, delimiter);
		}
		CSV_file.close();

		// store CSV data from .csv file to Eigen matrix :
		*Mat = SpEOMatrixD::Zero(CSVrow, CSVcol);
		CSV_file.open(fname_CSV);//CSV_file.open(fname_CSV.c_str());
		// skip as many lines as specified by skipLns
		for(int i=0; i<skipLns; i++){
			//ignore the first 500 characters, or until first \n, whichever is met first
			CSV_file.ignore(500,'\n');
		}
		for (unsigned int r=0; r<CSVrow; r++) {
			m_row.readNextRow(CSV_file, delimiter);
			for (unsigned int c=0; c<CSVcol; c++) {
				if (!(istringstream(m_row[c]) >> (*Mat)(r, c)))
					(*Mat)(r, c) = -999;
			}
		}
		return 1;
	}

        if(my_rank==0){
          cout << "["<< my_rank << "] ERROR: The file '" << fname_CSV << "' does not exist in the file system!" << endl;
          cerr << "["<< my_rank << "] ERROR: The file '" << fname_CSV << "' does not exist in the file system!" << endl;
	  //exit(2);
	}
	return -1;
}
// read and save .cvs file content to FLOAT matrix
int read_CSV(SpEOMatrixF *Mat, const char *fname_CSV, char delimiter, int skipLns) {
	int my_rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
	std::ifstream CSV_file(fname_CSV);//std::ifstream CSV_file(fname_CSV.c_str());
	//cout << "CSV_file.is_open() = " << CSV_file.is_open() << endl;
	if (CSV_file.is_open()) {
		unsigned int CSVrow = 0;
		for(int i=0; i<skipLns; i++){
			CSV_file.ignore(500,'\n');
			//ignore the first 500 characters, or until first \n, whichever is met first
		}
		CSVRow m_row;
		m_row.readNextRow(CSV_file, delimiter);
		unsigned int CSVcol = m_row.size();
		while (CSV_file.good()) {
			if (CSVcol != m_row.size()) {
				cerr << endl << "[" << my_rank << "] ERROR while reading CSV in row number "
						<< CSVrow + 1 << ": All rows must have the same length!" << endl;
				exit(2);
			}
			CSVrow++;
			m_row.readNextRow(CSV_file, delimiter);
		}
		CSV_file.close();

		// store CSV data from .csv file to Eigen matrix :
		*Mat = SpEOMatrixF::Zero(CSVrow, CSVcol);
		CSV_file.open(fname_CSV);//CSV_file.open(fname_CSV.c_str());
		// skip as many lines as specified by skipLns
		for(int i=0; i<skipLns; i++){
			//ignore the first 500 characters, or until first \n, whichever is met first
			CSV_file.ignore(500,'\n');
		}
		for (unsigned int r=0; r<CSVrow; r++) {
			m_row.readNextRow(CSV_file, delimiter);
			for (unsigned int c=0; c<CSVcol; c++) {
				if (!(istringstream(m_row[c]) >> (*Mat)(r, c)))
					(*Mat)(r, c) = -999;
			}
		}
		return 1;
	}
	return -1;
}
// read and save .cvs file content to INTEGER matrix
int read_CSV(SpEOMatrixI *Mat, const char *fname_CSV, char delimiter, int skipLns) {
	int my_rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
	std::ifstream CSV_file(fname_CSV);//std::ifstream CSV_file(fname_CSV.c_str());
	if (CSV_file.is_open()) {
		unsigned int CSVrow = 0;
		for(int i=0; i<skipLns; i++){
			CSV_file.ignore(500,'\n');
			//ignore the first 500 characters, or until first \n, whichever is met first
		}
		CSVRow m_row;
		m_row.readNextRow(CSV_file, delimiter);
		unsigned int CSVcol = m_row.size();
		while (CSV_file.good()) {
			if (CSVcol != m_row.size()) {
				cerr << endl << "[" << my_rank << "] ERROR while reading CSV in row number "
						<< CSVrow + 1 << ": All rows must have the same length!" << endl;
				exit(2);
			}
			CSVrow++;
			m_row.readNextRow(CSV_file, delimiter);
		}
		CSV_file.close();

		// store CSV data from .csv file to Eigen matrix :
		*Mat = SpEOMatrixI::Zero(CSVrow, CSVcol);
		CSV_file.open(fname_CSV);//CSV_file.open(fname_CSV.c_str());
		// skip as many lines as specified by skipLns
		for(int i=0; i<skipLns; i++){
			//ignore the first 500 characters, or until first \n, whichever is met first
			CSV_file.ignore(500,'\n');
		}
		for (unsigned int r=0; r<CSVrow; r++) {
			m_row.readNextRow(CSV_file, delimiter);
			for (unsigned int c=0; c<CSVcol; c++) {
				if (!(istringstream(m_row[c]) >> (*Mat)(r, c)))
					(*Mat)(r, c) = -999;
			}
		}
		return 1;
	}
	return -1;
}

void write_Mat_to_CSV(SpEOMatrixF *Mat, const char *fname_CSV) {
	//	std::ofstream CSV_file(fname_CSV);
	std::ofstream CSV_file;
	CSV_file.open(fname_CSV, ios::out | ios::ate | ios::app);
	if (CSV_file.is_open()) {
//		CSV_file << *Mat;
		IOFormat CSVFmt(FullPrecision, 0, ",");
		// or
//		IOFormat CSVFmt(FullPrecision, DontAlignCols, ",\t");
//		CSV_file << "# hallo welt \n" << Mat->format(CSVFmt);
		CSV_file << Mat->format(CSVFmt);
		CSV_file << "\n";
	}else{
		cout << endl << "WARNING: Matrix could not be written to .csv file!" << endl;
	}
	CSV_file.close();
}
void write_Mat_to_CSV(SpEOMatrixD *Mat, const char *fname_CSV) {
//	std::ofstream CSV_file(fname_CSV);
	std::ofstream CSV_file;
	CSV_file.open(fname_CSV, ios::out | ios::ate | ios::app);
	if (CSV_file.is_open()) {
//		CSV_file << *Mat;
		IOFormat CSVFmt(FullPrecision, 0, ",");
		// or
//		IOFormat CSVFmt(FullPrecision, DontAlignCols, ",\t");
//		CSV_file << "# hallo welt \n" << Mat->format(CSVFmt);
		CSV_file << Mat->format(CSVFmt);
		CSV_file << "\n";
	}else{
		cout << endl << "WARNING: CSV_file '" << CSV_file << "' could not be written!" << endl;
	}
	CSV_file.close();
//	system("chmod 700 fname_CSV");
	//chmod(fname_CSV, S_IRUSR | S_IWUSR | S_IRGRP | S_IWGRP | S_IROTH | S_IWOTH);
	chmod(fname_CSV, S_IRWXU | S_IRWXG | S_IRWXO);
}

void read_SRF(SpEOGlobalParams *glPrms, SpEOMatrixD *Mat, string fname_CSV, char delimiter, bool normalize_SRF) {
//void read_SRF(SpEOMatrixF *Mat, const char *fname_CSV, char delimiter) {
	int my_rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

	if(my_rank==0){
		cout << "\n"
		<< "###########################################################\n"
		<< "##                                                       ##\n"
		<< "##            read spectral response functions           ##\n"
		<< "##                                                       ##\n"
		<< "###########################################################"
		<< "\n";
	}

//	char *CSVHead;
	int stat_CSV_read = read_CSV(Mat, fname_CSV.c_str(),delimiter,0);
	if(stat_CSV_read==-1){
		if(my_rank==0){
			cout << "["<< my_rank << "] ERROR: The SRF file '" << fname_CSV.c_str() << "' does not exist in the file system!" << endl;
		}
		MPI_Barrier(MPI_COMM_WORLD);
		exit(2);
	}
	// check if all values are positive and transpose matrix.
	if ((*Mat).minCoeff() < 0.0) {
		if (my_rank == 0) {
			cerr << endl << "ERROR: SRF is not allowed to contain negative values!"	<< endl << endl ;
		}
		MPI_Barrier(MPI_COMM_WORLD);
		exit(2);
	}else{
		// check dimensions
		if(Mat->rows()==glPrms->NChX_orig && Mat->cols()==glPrms->NChY_orig){
			// That is how the SRFs are supposed to be stored.
			// -> nothing to be done!
		}else{
			if (my_rank == 0) {
				cerr << endl << "ERROR: The SRF matrix stored in the file " << endl << endl
				     << "   " << fname_CSV << endl << endl
				     << "has the wrong dimensions! It must have as many rows as number of MS channels and as many columns as number of HS channels!" << endl << endl;
			}
			MPI_Barrier(MPI_COMM_WORLD);
			exit(2);
		}
		//Mat->transposeInPlace();
	}
	// check sum to one constraint
	if(my_rank == 0){
		cout << "Check sum-to-one constraint. I.e. check if the SRFs are individually normalized..." << endl;
		cout << "The sum of the SRF entries in channel iChX is differs from 1.0 (ideal case) as follows (should be equal to 0.0)!" << endl;
	}
//	for(int iChX=0; iChX<glPrms->NChX; iChX++){
	for(int iChX=0; iChX<Mat->rows(); iChX++){
		double sum = 0.0;
		//double sum = Mat->row(iChX).sum();
//		for(int iChY=0; iChY<glPrms->NChY; iChY++){
		for(int iChY=0; iChY<Mat->cols(); iChY++){
			sum += Mat->coeff(iChX,iChY);
		}
		if(my_rank == 0){
			cout << "iChX=" << iChX << "  ->  sum(SRF("<< iChX <<"))-1 = " << sum-1 << endl;
		}
//		if(abs(sum-1>tol){
//			if(my_rank == 0){
//
//			}
//			MPI_Barrier(MPI_COMM_WORLD);
//			exit(2);
//		}
//		Mat->row(iChX) /= sum;
	}

	//cout << "SFR = " << endl << Mat->transpose() << endl << endl;
	// normalize each SRF (i.e. each row) so that the entries in each row sum up to one
	if (normalize_SRF) {
		if(my_rank == 0){
			cout << "normalize SRF, so that the sum in every row (i.e. for every channel of the image X) is equal to one..." << endl;
		}
		for(int iChX=0; iChX<Mat->rows(); iChX++){
			double sum = 0;
//			for(int iChY=0; iChY<glPrms->NChY; iChY++){
			for(int iChY=0; iChY<Mat->cols(); iChY++){
				sum += Mat->coeff(iChX,iChY);
			}
			if(sum<1e-6){
				if(my_rank == 0){
					cout << "ERROR: the sum of the SRF entries in channel " << iChX << " is close to or equal to zero. Row-wise normalization would create inf values in the SRF!" << endl;
				}
				MPI_Barrier(MPI_COMM_WORLD);
				exit(2);
			}
			Mat->row(iChX) /= sum;
		}
	}
	if(my_rank==0){
		cout << ".. spectral response functions successfully read!\n"
				<< "\n" ;
	}
}

//void createImFromTMPFiles(SpEOPaths *paths, SpEOFusionSetting *fSet, SpEOGlobalParams *glPrms, SpEODataset *ImZ){
//	// to be deleted >/////////////
//	//paths->dir_tmp_patches="/gpfs/work/pr45ne/ga39yoz2/tmp/patches/lambda1_patches0to15500/";
//        //paths->dir_tmp_patches="tmp/patches/allePatches";
//	// <////////// to be deleted
//	int uP, vP, iChY;
//	  int pszH=fSet->patchsize*glPrms->fDS;
//	  int pszH2=fSet->patchsize*fSet->patchsize*glPrms->fDS*glPrms->fDS;
//	  struct dirent *pDirent;
//	  DIR *pDir;
//	  pDir = opendir(paths->dir_tmp_patches.c_str());
//	  if (pDir == NULL) {
//		  cout << "Cannot open directory: " << paths->dir_tmp_patches << endl;
//		  exit(2);
//	  }
//	  SpEOMatrixF rmforMat = SpEOMatrixF::Zero(ImZ->get_sizeU(), ImZ->get_sizeV());
//	  for (uP=0; uP<glPrms->NPU; uP++) {
//		  for (vP=0; vP<glPrms->NPV; vP++) {
//			  rmforMat.block(glPrms->idxPUH.coeff(uP), glPrms->idxPVH.coeff(vP), pszH, pszH) =
//					  rmforMat.block(glPrms->idxPUH.coeff(uP), glPrms->idxPVH.coeff(vP), pszH, pszH)
//					  + SpEOMatrixF::Ones(pszH, pszH);
//		  }
//	  }
//	  while ((pDirent = readdir(pDir)) != NULL){
//		  if (pDirent->d_name[0] != '.'){
//			  char buf [paths->dir_tmp_patches.length()+strlen(pDirent->d_name)+4];
//			  sprintf (buf, "%s/%s",paths->dir_tmp_patches.c_str(), pDirent->d_name);
//			  SpEOMatrixF myReadMat;
//			  read_CSV(&myReadMat, buf, ',', 0);
//			  std::string uP_str(pDirent->d_name+7, pDirent->d_name+11);
//			  std::string vP_str(pDirent->d_name+13, pDirent->d_name+17);
//			  uP= atoi(uP_str.c_str());
//			  vP= atoi(vP_str.c_str());
//			  //cout << "read file: " << pDirent->d_name << ". Read " << myReadMat.cols() << " columns and " << myReadMat.rows() << "rows. uP=" << uP << " and vP=" << vP << endl;//", myReadMat=" << endl << myReadMat << endl;
//			  if(myReadMat.rows()!=pszH2){
//				  cout << "WARNING: size of temporarily stored patch does not match the actual patch size!" << endl;
//			  }
//			  SpEOVectorF patchVecHR = SpEOVectorF::Zero(myReadMat.rows());
//			  for(iChY=0; iChY<glPrms->NChY; iChY++) {
//				  SpEOVectorF patchVecHR = myReadMat.col(iChY);
//				  ImZ->writePatch(patchVecHR, pszH, glPrms->idxPUH.coeff(uP), glPrms->idxPVH.coeff(vP), iChY);
//			  }
//		  }
//	  }
//	  closedir (pDir);
//	  // divide each band of ImZ coefficient-wise by rmforMat
//	  for (iChY=0; iChY<glPrms->NChY; iChY++) {
//		  ImZ->replBandByCwiseCuotient(&rmforMat, iChY);
////			  ImZ.cwiseQuotient(rmforMat); to be done
//	  }
//	  //ImZ->cutNegCoeff();
//	  ImZ->fillArrayFromMat();
//}


void calcDecisionMat(SpEOMatrixD *SRF, SpEODataIOSetting *dSet, SpEOFusionSetting *fSet, SpEOParallelSetting *pSet,
		SpEOGlobalParams *glPrms) {
    int my_rank; int my_processes;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &my_processes);

	if(my_rank==0){
		cout<< "\n"
		<< "###########################################################\n"
		<< "##                                                       ##\n"
		<< "##                calculate decision matrix              ##\n"
		<< "##                                                       ##\n"
		<< "###########################################################"
		<< "\n";
	}

	if(fSet->Nc > glPrms->NChZ){
		if(my_rank==0){
			cout << "WARNING: Nc got corrected from Nc=" << fSet->Nc << " to Nc=" << glPrms->NChZ << "==NChZ" << endl;
		}
		fSet->Nc = glPrms->NChZ;
	}
	if(fSet->Nc == glPrms->NChZ && fSet->No !=0){
		if(my_rank==0){
			cout << "NOTICE: Since Nc==NChZ, no spectral overlap will be accounted for. Therefore No does not affect the calculation. It got set from No=" << fSet->No << " to No=0." << endl;
		}
		fSet->No = 0;
	}
	else if(fSet->No>=fSet->Nc){
		if(my_rank==0){
			cout << "WARNING: No got corrected from No=" << fSet->No << " to the maximal overlap of No=" << fSet->Nc-1 << "=Nc-1" << endl;
		}
		fSet->No = fSet->Nc-1;
	}
	int a = fSet->Nc - fSet->No;
	glPrms->Ng = ceil(((double) (glPrms->NChZ - fSet->Nc)) / ((double) (a))) + 1;
	glPrms->decMat_C = SpEOMatrixD::Zero(glPrms->NChX, glPrms->Ng);
	for (int iChX = 0; iChX < glPrms->NChX; iChX++) {
		double sum = 0;
		for (int iG = 0; iG < glPrms->Ng - 1; iG++) {
			for (int iC = 0; iC < fSet->Nc; iC++) {
				sum = sum + (double)(*SRF)(iChX, dSet->chBundleFirst + a*iG+iC); // tbd: Check offset dSet->chBundleFirst
			}
			glPrms->decMat_C(iChX, iG) = sum;
			sum = 0;
		}
		for (int iC = 0; iC < fSet->Nc; iC++) {
			sum = sum + (double)(*SRF)(iChX, dSet->chBundleFirst + glPrms->NChZ - fSet->Nc+iC); // HERE, an OFFSET of dSetting->chBundleFirst before glPrms->NChZ may have to be added!
		}
		glPrms->decMat_C(iChX, glPrms->Ng - 1) = sum;
	}
	// L1-normalize C columnwise
	for (int iG = 0; iG < glPrms->Ng; iG++) {
		double manhNorm = glPrms->decMat_C.col(iG).sum();
		if (manhNorm == 0.0) {
			glPrms->decMat_C.col(iG) = SpEOVectorD::Constant(glPrms->NChX,	1. / glPrms->NChX);
		} else {
			glPrms->decMat_C.col(iG) /= manhNorm;
		}
	}
	// Thresholding below fSet->tol_SRF
	// &  calculate the number of non-zero coefficients in each row (i.e. along each MS band)
	glPrms->Nm = new int[glPrms->NChX];
	for (int iChX = 0; iChX < glPrms->NChX; iChX++) {
		int nonZeros = 0;
		for (int iG = 0; iG < glPrms->Ng; iG++) {
			if (glPrms->decMat_C.col(iG).maxCoeff() < fSet->tol_SRF) {
				SpEOVectorD tmpp = glPrms->decMat_C.col(iG);
				SpEOMatrixD::Index maxCoeffIdx;
				double maxValue = tmpp.maxCoeff(&maxCoeffIdx);
				glPrms->decMat_C(maxCoeffIdx, iG) = (fSet->tol_SRF+1)/2;
			}
			if (glPrms->decMat_C(iChX, iG) < fSet->tol_SRF) {
				glPrms->decMat_C(iChX, iG) = 0.0;
			} else {
				nonZeros++;
			}
		}
		glPrms->Nm[iChX] = nonZeros;
	}
	// save the coordinates of all non-zero coefficients in each row (i.e. along the groups in one MS band)
	glPrms->km = new int*[glPrms->NChX];
	
	glPrms->numProbPerPatch = 0; // renamed from myChN
	for (int iChX = 0; iChX < glPrms->NChX; iChX++) {
		glPrms->km[iChX] = new int[glPrms->Nm[iChX]];
		int cnt = 0;
		for (int iG = 0; iG < glPrms->Ng; iG++) {
			if (glPrms->decMat_C(iChX, iG) > 0.0) {
				glPrms->km[iChX][cnt] = iG;
				cnt++;
				glPrms->numProbPerPatch++;
			}
		}
	}
	if(glPrms->numProbPerPatch < pSet->numProcPerPatch){
		if(my_rank == 0){
			cout << endl << "> WARNING: numProcPerPatch was set too great! It got corrected from " << pSet->numProcPerPatch << " to " << glPrms->numProbPerPatch << ", which is the total number of problems per patch <" << endl << endl;
//			cerr << endl << "> ERROR: the number of processes per patch (" << pSet->numProcPerPatch << ") must be less or equal to (and ideally a divisor of) the number of problems per patch (" << glPrms->numProbPerPatch << ") <" << endl << endl;
//			exit(2);
		}
		pSet->numProcPerPatch = glPrms->numProbPerPatch;
	}
	if(glPrms->numProbPerPatch % pSet->numProcPerPatch!=0){
		if(my_rank == 0){
			cout << endl << "> WARNING: the number of processes per patch (" << pSet->numProcPerPatch << ") should ideally be a divisor of the number of problems per patch (" << glPrms->numProbPerPatch << ") <" << endl << endl;
		}
	}
	glPrms->numPatchGroups = my_processes / pSet->numProcPerPatch;
	if(glPrms->numPatchGroups == 0){
		if(my_rank == 0){
			cerr << endl << ">>> ERROR: too few CPUs! Number of processes must be at least equal to (and ideally a multiple of) " << pSet->numProcPerPatch << " <<<" << endl << endl;
		}
		MPI_Barrier(MPI_COMM_WORLD);
		exit(2);
	}

	if(my_processes % pSet->numProcPerPatch != 0){
		if(my_rank == 0){
		 cout << endl << "WARNING: Total number of processes should be a multiple of " << pSet->numProcPerPatch << " for optimal performance! Currently, there are " << my_processes << " processes in total which of " << my_processes % pSet->numProcPerPatch << " are idling!" << endl << endl;
		}
	}
	// L1-normalize C columnwise
	for (int iG = 0; iG < glPrms->Ng; iG++) {
		double manhNorm = glPrms->decMat_C.col(iG).sum();
		glPrms->decMat_C.col(iG) /= manhNorm;
	}
	int *startProc = new int[glPrms->NChX];
        int sum=0;
	for(int iChX=0; iChX < glPrms->NChX; iChX++){
		startProc[iChX] = sum;
		sum += glPrms->Nm[iChX];
	}
	

	/*int myChX = 0;
	// while(startProc[myChX]<=my_rank % pSet->numProcPerPatch) { // semantic error + memory error
	while(startProc[myChX]+glPrms->Nm[myChX] <= my_rank % pSet->numProcPerPatch){
		myChX++;
	}

	glPrms->myChX = myChX;
	glPrms->myBundle = glPrms->km[glPrms->myChX][my_rank % pSet->numProcPerPatch-startProc[glPrms->myChX]];*/
	
	// new: we have to store more than 1 channel now
	if(my_rank % pSet->numProcPerPatch < glPrms->numProbPerPatch % pSet->numProcPerPatch){
		glPrms->myNumProbPerPatch = glPrms->numProbPerPatch / pSet->numProcPerPatch + 1;
	}
	else {
		glPrms->myNumProbPerPatch = glPrms->numProbPerPatch / pSet->numProcPerPatch;
	}

	//if(my_rank==0){cout << "glPrms->myNumProbPerPatch = " << glPrms->myNumProbPerPatch << endl;}
	//if(my_rank==0){cout << "pSet->numProcPerPatch = " << pSet->numProcPerPatch << endl;}

	glPrms->myChX = new int[glPrms->myNumProbPerPatch];
	glPrms->myBundle = new int[glPrms->myNumProbPerPatch];
	for (int ipp=0; ipp < glPrms->myNumProbPerPatch; ipp++) {
		glPrms->myChX[ipp] = 0;
		while(startProc[glPrms->myChX[ipp]]+glPrms->Nm[glPrms->myChX[ipp]] <= my_rank % pSet->numProcPerPatch + ipp*pSet->numProcPerPatch){
			glPrms->myChX[ipp]++;
		}
		glPrms->myBundle[ipp] = glPrms->km[glPrms->myChX[ipp]][my_rank % pSet->numProcPerPatch + ipp*pSet->numProcPerPatch - startProc[glPrms->myChX[ipp]]];

	}
	// calculate Nm2 and km2
	glPrms->Nm2 = new int[glPrms->Ng];
	glPrms->km2 = new int*[glPrms->Ng];
	for (int iG = 0; iG < glPrms->Ng; iG++) {
		glPrms->Nm2[iG] = 0;
		for (int iChX = 0; iChX < glPrms->NChX; iChX++) {
			if (glPrms->decMat_C(iChX, iG) > 0) {
				glPrms->Nm2[iG]++;
			}
		}
		// km2: save the coordinates of all non-zero coefficients in each column (i.e. along the MS bands in one group)
		glPrms->km2[iG] = new int[glPrms->Nm2[iG]];
		int cnt = 0;
		for (int iChX = 0; iChX < glPrms->NChX; iChX++) {
			if (glPrms->decMat_C(iChX, iG) > 0) {
				glPrms->km2[iG][cnt] = iChX;
				cnt++;
			}
		}
	}
	// What is that good for? Do we need to adapt numProbPerPatch again?? Isn't that doubled here? [No, that's true. It looks redundant. I'll delete it after double checking]
	//glPrms->numProbPerPatch = 0;
	//for(int iChX=0; iChX<glPrms->NChX; iChX++){
	//	glPrms->numProbPerPatch += glPrms->Nm[iChX];
	//}
	
	if (my_rank == 0) {
		cout << ".. decision matrix successfully calculated!" << endl << "glPrms->decMat_C.transpose() = " << endl << glPrms->decMat_C.transpose() << endl;
	}

	// ###################################
	// #  Calculate Nc_vec, P_lmd, etc.  #
	// ###################################
	int iG, iChZ, iC, ipp;
//	if(my_rank==0) cout << endl << "bp1:" << endl << "  glPrms->numProbPerPatch   = " << glPrms->numProbPerPatch << endl
//			                                      << "  glPrms->myNumProbPerPatch = " << glPrms->myNumProbPerPatch << endl
//                                                  << "  glPrms->Ng = "                << glPrms->Ng << endl;
	// Nc_vec: Numbers of HS channels in channel groups
	// (in case of Spectral Grouping, all values in Nc_vec are identical and equal to the Nc value specified by the user in fSettings)
	glPrms->Nc_vec = new int[glPrms->numProbPerPatch];
//	glPrms->Nc_vec = new int[glPrms->Ng];
	for(ipp=0; ipp<glPrms->numProbPerPatch; ipp++){
//	for(iG=0; iG<glPrms->Ng; iG++){
		// in case of Spectral Grouping (will need to be adapted to varying numbers later, when the correlation based channel grouping is implemented)
//		glPrms->Nc_vec[iG] = fSet->Nc;
		glPrms->Nc_vec[ipp] = fSet->Nc;
	}
	//if(my_rank==0){cout << "bp1"<< endl;}
//	if(my_rank==0) cout << endl << "bp2:" << endl; for(ipp=0; ipp<glPrms->numProbPerPatch; ipp++) { cout << "   glPrms->Nc_vec[ipp="<< ipp <<"] = " << glPrms->Nc_vec[ipp] << endl;}
	// initialize the vectors (containing the non-trivial entries along the diagonal in the corresponding blocks) in P_lmd
//	glPrms->P_lmd_vecs = new SpEOVectorD[glPrms->Ng];
	glPrms->P_lmd_vecs = new SpEOVectorD[glPrms->numProbPerPatch];
//	for(iG=0; iG<glPrms->Ng; iG++){
	for(ipp=0; ipp<glPrms->numProbPerPatch; ipp++){
//		glPrms->P_lmd_vecs[iG] = SpEOVectorD::Constant(glPrms->Nc_vec[iG],1);
		glPrms->P_lmd_vecs[ipp] = SpEOVectorD::Ones(glPrms->Nc_vec[ipp])*glPrms->decMat_C(glPrms->myChX[ipp], glPrms->myBundle[ipp]);
	}
	//if(my_rank==0){cout << "bp2"<< endl;}
//	if(my_rank==0) cout << endl << "bp3:" << endl; for(ipp=0; ipp<glPrms->numProbPerPatch; ipp++) { cout << "   glPrms->P_lmd_vecs[ipp="<< ipp <<"].transpose() = " << glPrms->P_lmd_vecs[ipp].transpose() << endl;}
	// For every block (1,...,Ng) calculate the absolute (total) row- and column index in the actual sparse matrix P_lmd. Row indexes are stored in the first row and the corresponding column indexes are stored in the second row.
//	glPrms->P_lmd_idx_bl = new int*[glPrms->Ng];
	glPrms->P_lmd_idx_bl = new int*[glPrms->numProbPerPatch];
//	array = new int[arg1][arg2]
//	int P_lmd_idx_bl[glPrms->Ng][2];

	int col_idx=0, row_idx;

//	for(iG=0; iG<glPrms->Ng-1; iG++){
	for(ipp=0; ipp<glPrms->numProbPerPatch-glPrms->Nm2[glPrms->Ng-1]; ipp++){
		iG = glPrms->myBundle[ipp];
		glPrms->P_lmd_idx_bl[ipp] = new int[2];

		// =====> row_idx = index of the first HS channel in the current  group (in the simple case of Spectral Grouping, all groups (possibly except for the last group) are equally spaced, i.e. row_idx is linearly increasing)
		row_idx = iG*(fSet->Nc - fSet->No);
		// <=====
		glPrms->P_lmd_idx_bl[ipp][0] = row_idx;
		glPrms->P_lmd_idx_bl[ipp][1] = col_idx;
		col_idx += glPrms->Nc_vec[ipp];
	}
	//if(my_rank==0){cout << "bp3"<< endl;}
	// treat last group separately:
//	iG = glPrms->Ng-1;
//	glPrms->P_lmd_idx_bl[iG] = new int[2];
//	row_idx = glPrms->NChZ - glPrms->Nc_vec[iG];
//	glPrms->P_lmd_idx_bl[iG][0] = row_idx;
//	glPrms->P_lmd_idx_bl[iG][1] = col_idx;
//	if(my_rank==0) cout << endl << "bp4:" << endl;
	for(ipp=glPrms->numProbPerPatch-glPrms->Nm2[glPrms->Ng-1]; ipp<glPrms->numProbPerPatch; ipp++){
//		iG = glPrms->Ng-1;
		glPrms->P_lmd_idx_bl[ipp] = new int[2];
		row_idx = glPrms->NChZ - glPrms->Nc_vec[ipp];
		glPrms->P_lmd_idx_bl[ipp][0] = row_idx;
		glPrms->P_lmd_idx_bl[ipp][1] = col_idx;
		col_idx += glPrms->Nc_vec[ipp];
	}
	//if(my_rank==0){cout << "bp4"<< endl;}

//	for(int ipp=0; ipp<glPrms->numProbPerPatch-glPrms->Nm2[glPrms->Ng-1]; ipp++){
//		cout << "bp1.. ipp=" << ipp << endl;
//		cout << "glPrms->P_lmd_idx_bl[ipp][0] = " << glPrms->P_lmd_idx_bl[ipp][0] << endl;
//		cout << "glPrms->P_lmd_idx_bl[ipp][1] = " << glPrms->P_lmd_idx_bl[ipp][1] << endl;
//	}

//	if(my_rank==0) cout << endl << "bp5:" << endl; for(ipp=0; ipp<glPrms->numProbPerPatch; ipp++) { cout << "   glPrms->P_lmd_idx_bl[ipp="<< ipp <<"] = [" << glPrms->P_lmd_idx_bl[ipp][0] << "," << glPrms->P_lmd_idx_bl[ipp][1] << "]" << endl;}

	//col_idx += glPrms->Nc_vec[glPrms->Ng-1];

	// for every row (each corresponding to one HS channel iChY) calculate the relevant corresponding block indexes (first row) and the corresponding entry indexes relative to the corresponding block's origin each starting at 0.
	glPrms->P_lmd_idx_row = new SpEOMatrixI[glPrms->NChZ];
	// first the number of non-trivial entries in each row (iChZ=0,...,NChZ-1) needs to be counted
	SpEOVectorI entries_cnt     = SpEOVectorI::Zero(glPrms->NChZ);
	SpEOVectorI entries_divider = SpEOVectorI::Zero(glPrms->NChZ);
//	for(iG=0; iG<glPrms->Ng; iG++){
	bool incr_entries_divider = false;
	bool iG_done[glPrms->Ng];
	for(iG=0; iG<glPrms->Ng; iG++){
		iG_done[iG]=false;
	}



	/*
	if(my_rank==0){
		cout << endl
			 << "glPrms->Ng = " <<glPrms->Ng << endl
			 << "glPrms->numProbPerPatch = " << glPrms->numProbPerPatch << endl;
		for(ipp=0; ipp<glPrms->numProbPerPatch; ipp++){
			cout << "  glPrms->Nc_vec[ipp"<<ipp<<"] = " << glPrms->Nc_vec[ipp] << endl
			     << "  iG = glPrms->myBundle[ipp="<<ipp<<"] = " << glPrms->myBundle[ipp] << endl;
		}

		cout << "glPrms->Ng = " <<glPrms->Ng << endl
			 << "###########################" << endl
			 << "# print P_lmd block-wise  #" << endl
			 << "###########################" << endl;
		for(iG=0; iG<glPrms->Ng; iG++){
			cout << "iG=" << iG << " and P_lmd_idx_bl[iG][0]=row_idx="<< glPrms->P_lmd_idx_bl[iG][0]
								<< " and P_lmd_idx_bl[iG][1]=col_idx="<< glPrms->P_lmd_idx_bl[iG][1] << endl;
	//			 << " P_lmd in this block = ";
	//		for (iC=0; iC<glPrms->Nc_vec[iG];iC++){
	//			cout << " " << glPrms->P_lmd_vecs[iG](iC);
	//		}
	//		cout << endl;
		}
		// save the coordinates of all non-zero coefficients in each row (i.e. along the groups in one MS band)
		for (int iChX = 0; iChX < glPrms->NChX; iChX++) {
			int cnt = 0;
			for (int iG = 0; iG < glPrms->Ng; iG++) {
				if (glPrms->decMat_C(iChX, iG) > 0.0) {
					cout << "glPrms->km[iChX="<<iChX<<"][cnt="<<cnt<<"] = iG = " << glPrms->km[iChX][cnt] << endl;
					cnt++;
				}
			}
		}
	}
	*/





	//if(my_rank==0){cout << "bp5"<< ", for glPrms->myBundle: glPrms->myNumProbPerPatch=" << glPrms->myNumProbPerPatch << ", for Nc_vec: glPrms->numProbPerPatch=" << glPrms->numProbPerPatch << ", for entries_divider: glPrms->NChZ=" << glPrms->NChZ << ", glPrms->NChY=" << glPrms->NChY << endl;}
	for(ipp=0; ipp<glPrms->numProbPerPatch; ipp++){
		//if(my_rank==0){cout << "bp5.ipp=" << ipp <<  ", glPrms->myBundle[ipp]=" << glPrms->myBundle[ipp] << endl;}
		iG = glPrms->myBundle[ipp];
		if(!iG_done[iG]){
			incr_entries_divider = true;
			iG_done[iG] = true;
		}else{
			incr_entries_divider = false;
		}
//		if(my_rank==0){cout << "bp6: ipp="<< ipp << " and iG="<<iG << endl;}
		for(iC=0; iC<glPrms->Nc_vec[ipp]; iC++){

			//if(my_rank==0){cout << "bp5.ipp=" << ipp <<  ", glPrms->myBundle[ipp]=" << glPrms->myBundle[ipp] << ", glPrms->P_lmd_idx_bl[ipp][0]+iC=" << glPrms->P_lmd_idx_bl[ipp][0]+iC << "=iChZ"<< endl;}

			iChZ = glPrms->P_lmd_idx_bl[ipp][0]+iC;
			if(incr_entries_divider){//if(glPrms->Nm2[iG]>1 && newgroup){
				//if(my_rank==0){cout << "inside.." << endl;}
				entries_divider(iChZ)++;
			}
			//if(my_rank==0){cout << "outside.." << endl;}
			entries_cnt(iChZ)++;
			//if(my_rank==0){cout << "done!" << endl;}
		}
	}
	//if(my_rank==0){cout << "bp6"<< endl;}
//	if(my_rank==0) cout << endl << "bp7:" << endl << "   entries_cnt = " << endl << entries_cnt.transpose() << endl << "entries_divider = " << endl << entries_divider.transpose() << endl << endl;
	for(iChZ=0; iChZ<glPrms->NChZ; iChZ++){
		glPrms->P_lmd_idx_row[iChZ] = SpEOMatrixI::Zero(2,entries_cnt.coeff(iChZ));
	}
	SpEOVectorI entries_cnt_tmp = SpEOVectorI::Zero(glPrms->NChZ);
//	for(iG=0; iG<glPrms->Ng; iG++){
	for(ipp=0; ipp<glPrms->numProbPerPatch; ipp++){
		for(iC=0; iC<glPrms->Nc_vec[ipp]; iC++){
			iChZ = glPrms->P_lmd_idx_bl[ipp][0]+iC;//glPrms->P_lmd_idx_bl[iG][0]+iC;
//			for(int ii=0; ii<glPrms->P_lmd_idx_row[iChZ].cols(); ii++){
				// store block index
				glPrms->P_lmd_idx_row[iChZ](0,entries_cnt_tmp(iChZ)) = ipp;//iG;
				// store coefficient index relative to corresponding block origin
				glPrms->P_lmd_idx_row[iChZ](1,entries_cnt_tmp(iChZ)) = iC;
//			}
			// calculate value in P_lmd (here implemented is only simple averaging!! If desired, that can be changed later)
//			glPrms->P_lmd_vecs[iG](iC) /= entries_divider(iChZ);//entries_cnt(iChZ);
			glPrms->P_lmd_vecs[ipp](iC) /= entries_divider(iChZ);//entries_cnt(iChZ);
//			// store block index
//			glPrms->P_lmd_idx_row[iChZ](0,entries_cnt_tmp(iChZ)) = iG;
//			// store coefficient index relative to corresponding block origin
//			glPrms->P_lmd_idx_row[iChZ](1,entries_cnt_tmp(iChZ)) = iC;
//			// calculate value in P_lmd (here implemented is only simple averaging!! If desired, that can be changed later)
//			glPrms->P_lmd_vecs[iG](iC) /= entries_cnt(iChZ);
			entries_cnt_tmp(iChZ)++;
		}
	}
	//if(my_rank==0){cout << "bp7"<< endl;}
	//if(my_rank==0) cout << endl << "bp8:" << endl; for(iChZ=0; iChZ<glPrms->NChZ; iChZ++)        { cout << "   glPrms->P_lmd_idx_row[iChZ"<<iChZ<<"] = (row0=):" << glPrms->P_lmd_idx_row[iChZ].row(0) << " (row1=): "<<glPrms->P_lmd_idx_row[iChZ].row(1) << endl;}
	//if(my_rank==0) cout << endl << "bp9:" << endl; for(ipp=0; ipp<glPrms->numProbPerPatch; ipp++){ cout << "   glPrms->P_lmd_vecs[ipp"<<ipp<<"] = " << glPrms->P_lmd_vecs[ipp].transpose() << endl;}


	///*
	if(my_rank==0){
		//
		// print P_lmd to check for correctness
		/*
		cout << endl
			 << "###########################" << endl
			 << "# print P_lmd block-wise  #" << endl
			 << "###########################" << endl;
		for(iG=0; iG<glPrms->Ng; iG++){
			cout << "iG=" << iG << " and P_lmd_idx_bl[iG][0]=row_idx="<< glPrms->P_lmd_idx_bl[iG][0]
								<< " and P_lmd_idx_bl[iG][1]=col_idx="<< glPrms->P_lmd_idx_bl[iG][1]
				 << " P_lmd in this block = ";
			for (iC=0; iC<glPrms->Nc_vec[iG];iC++){
				cout << " " << glPrms->P_lmd_vecs[iG](iC);
			}
			cout << endl;
		}
		*/
		/*
		cout << endl
			 << "###########################" << endl
			 << "# print P_lmd row-wise #" << endl
			 << "###########################" << endl;
		for(iChZ=0; iChZ<glPrms->NChZ; iChZ++){
			cout << "iChZ=" << iChZ << "   ->   P_lmd_idx_row[iChZ]= ";
			for(int ii=0; ii<glPrms->P_lmd_idx_row[iChZ].cols(); ii++){
				int block  = glPrms->P_lmd_idx_row[iChZ].coeff(0,ii);
				int relidx = glPrms->P_lmd_idx_row[iChZ].coeff(1,ii);
				cout << "    (block=" << block
				   << ", relidx=" << relidx
				   << ", val=" << glPrms->P_lmd_vecs[block](relidx)
				   << ") ";
			}
			cout << endl;
		}
		*/
		/*cout << endl
			 << "###########################" << endl
			 << "# print P_lmd block-wise  #" << endl
			 << "###########################" << endl;
		for(ipp=0; ipp<glPrms->numProbPerPatch; ipp++){
			cout << "ipp=" << ipp << " and P_lmd_idx_bl[ipp][0]=row_idx="<< glPrms->P_lmd_idx_bl[ipp][0]
								<< " and P_lmd_idx_bl[ipp][1]=col_idx="<< glPrms->P_lmd_idx_bl[ipp][1]
				 << " P_lmd in this block = ";
			for (iC=0; iC<glPrms->Nc_vec[ipp];iC++){
				cout << " " << glPrms->P_lmd_vecs[ipp](iC);
			}
			cout << endl;
		}
		cout << endl
			 << "###########################" << endl
			 << "# print P_lmd row-wise #" << endl
			 << "###########################" << endl;
		for(iChZ=0; iChZ<glPrms->NChZ; iChZ++){
			cout << "iChZ=" << iChZ << " -> P_lmd_idx_row[iChZ="<< iChZ <<"]= ";
			for(int ii=0; ii<glPrms->P_lmd_idx_row[iChZ].cols(); ii++){
				int block  = glPrms->P_lmd_idx_row[iChZ].coeff(0,ii);
				int relidx = glPrms->P_lmd_idx_row[iChZ].coeff(1,ii);
				cout << "(block=" << block
				   << ", relidx=" << relidx
				   << ", val=" << glPrms->P_lmd_vecs[block](relidx)
				   << ") ";
			}
			cout << endl;
		}
		*/
	}




	/******************************************************************************
	 *  Naive calculation of P_lambda as large sparse matrix including all zeros  *
	 ******************************************************************************
	// Nc_start_pos: vector containing the index of the first HS channel in each group
	// (in case of Spectral Grouping, all groups (possibly except for the last group) are equally spaced, i.e. Nc_start_pos is linearly increasing)
	glPrms->Nc_start_pos = new int[glPrms->Ng];
	// in case of Spectral Grouping (will need to be adapted to varying numbers later, when the correlation based channel grouping is implemented)
	//glPrms->Nc_start_pos[]
	a = fSet->Nc - fSet->No;
	for (int iG = 0; iG < glPrms->Ng - 1; iG++) {
		glPrms->Nc_start_pos[iG] = a*iG;
	}
	glPrms->Nc_start_pos[glPrms->Ng-1] = glPrms->NChZ - fSet->Nc;

	glPrms->sum_Nc_vec = 0;
	for(int iG=0; iG<glPrms->Ng; iG++){
		glPrms->sum_Nc_vec += glPrms->Nc_vec[iG];
		cout << "iG=" << iG << ", glPrms->Nc_vec[iG="<<iG<<"]=" << glPrms->Nc_vec[iG] << ", and Nc_start_pos[iG="<<iG<<"]=" << glPrms->Nc_start_pos[iG] << endl;
	}
	cout << endl << "glPrms->Ng="<< glPrms->Ng <<", glPrms->sum_Nc_vec=" << glPrms->sum_Nc_vec << endl<< endl;
	glPrms->P_lmd = SpEOMatrixF::Zero(glPrms->NChY,glPrms->sum_Nc_vec);
	SpEOVectorI entry_cnt = SpEOVectorI::Zero(glPrms->NChZ);
	int startpos_grp = 0;
	for(int iG=0; iG<glPrms->Ng; iG++){
		for(int iC=0; iC<glPrms->Nc_vec[iG]; iC++){
			// only for the simple spectral grouping case:
			glPrms->P_lmd(glPrms->Nc_start_pos[iG]+iC, startpos_grp+iC) = 1;
			entry_cnt(glPrms->Nc_start_pos[iG]+iC)++;
		}
		startpos_grp += glPrms->Nc_vec[iG];
	}
	// find indexes of all non-trivial entries for efficient implementation
	startpos_grp = 0;

	glPrms->P_lmd_idx = new SpEOVectorI[glPrms->NChZ];
	for (int iChZ=0; iChZ<glPrms->NChZ; iChZ++){
		glPrms->P_lmd_idx[iChZ] = SpEOVectorI::Zero(entry_cnt.coeff(iChZ));
	}
	entry_cnt = SpEOVectorI::Zero(glPrms->NChZ);
	for(int iG=0; iG<glPrms->Ng; iG++){
		for(int iC=0; iC<glPrms->Nc_vec[iG]; iC++){
			// only for the simple spectral grouping case:
			//glPrms->P_lmd(glPrms->Nc_start_pos[iG]+iC, startpos_grp+iC) = 1;
			glPrms->P_lmd_idx[glPrms->Nc_start_pos[iG]+iC](entry_cnt(glPrms->Nc_start_pos[iG]+iC)) = startpos_grp+iC;
			entry_cnt(glPrms->Nc_start_pos[iG]+iC)++;
		}
		startpos_grp += glPrms->Nc_vec[iG];
	}
	// only for the simple spectral grouping case: simply average all overlapping channels:
	int sum_vec[glPrms->NChZ];
	for(int iChZ=0; iChZ<glPrms->NChZ; iChZ++){
		sum_vec[iChZ] = 0;
		for(int ii=0; ii<glPrms->P_lmd_idx[iChZ].rows(); ii++){
			sum_vec[iChZ] += glPrms->P_lmd(iChZ,glPrms->P_lmd_idx[iChZ](ii));
		}
		if(sum_vec[iChZ]==0){
			if(my_rank==0){
				cerr << "ERROR: Devision by zero during the calculation of P_lmd!" << endl << endl;
			}
			MPI_Barrier(MPI_COMM_WORLD);
			exit(2);
		}
		for(int ii=0; ii<glPrms->P_lmd_idx[iChZ].rows(); ii++){
			glPrms->P_lmd(iChZ,glPrms->P_lmd_idx[iChZ](ii)) /= sum_vec[iChZ];
		}
	}
	//cout << "glPrms->P_lmd.block(0,0,10,glPrms->sum_Nc_vec) = " << endl << glPrms->P_lmd.block(0,0,10,glPrms->sum_Nc_vec) << endl << endl;
	//cout << "glPrms->P_lmd.block(150,0,10,glPrms->sum_Nc_vec) = " << endl << glPrms->P_lmd.block(150,0,10,glPrms->sum_Nc_vec) << endl << endl;

	*/

	delete[] startProc;
}

int remove_dir(const char *path){
   DIR *d = opendir(path);
   size_t path_len = strlen(path);
   int r = -1;
   if(d){
      struct dirent *p;
      r = 0;
      while(!r && (p=readdir(d))){
          int r2 = -1;
          char *buf;
          size_t len;
          /* Skip the names "." and ".." as we don't want to recurse on them. */
          if (!strcmp(p->d_name, ".") || !strcmp(p->d_name, "..")){
             continue;
          }
          len = path_len + strlen(p->d_name) + 2;
          buf = (char*)malloc(len);
          if(buf){
             struct stat statbuf;
             snprintf(buf, len, "%s/%s", path, p->d_name);
             if (!stat(buf, &statbuf)){
                if (S_ISDIR(statbuf.st_mode)){
                   r2 = remove_dir(buf);
                }else{
                   r2 = unlink(buf);
                }
             }
             free(buf);
          }
          r = r2;
      }
      closedir(d);
   }
   if(!r){
      r = rmdir(path);
   }

   return r;
}

void save_dSetting(SpEOPaths *paths, SpEODataIOSetting *dSetting){
	string dir_dSetting = paths->dir_out + "/" + "dSetting";
	cout << "write data I/O settings to files in directory: " << endl << "     " << dir_dSetting << " .. ";
	mkdir(dir_dSetting.c_str(), 0777);
	chmod(dir_dSetting.c_str(), S_IRWXU | S_IRWXG | S_IRWXO);
	SpEOMatrixD tmp_mat = SpEOMatrixD::Constant(1,1,-99999);
	string fname_tmp="";

	tmp_mat(0,0) = dSetting->chBundleFirst;
	fname_tmp = dir_dSetting + "/" + "chBundleFirst.csv";
	write_Mat_to_CSV(&tmp_mat, fname_tmp.c_str());

	tmp_mat(0,0) = dSetting->chBundleLast;
	fname_tmp = dir_dSetting + "/" + "chBundleLast.csv";
	write_Mat_to_CSV(&tmp_mat, fname_tmp.c_str());

	tmp_mat(0,0) = dSetting->uLFirst;
	fname_tmp = dir_dSetting + "/" + "uLFirst.csv";
	write_Mat_to_CSV(&tmp_mat, fname_tmp.c_str());

	tmp_mat(0,0) = dSetting->uLLast;
	fname_tmp = dir_dSetting + "/" + "uLLast.csv";
	write_Mat_to_CSV(&tmp_mat, fname_tmp.c_str());

	tmp_mat(0,0) = dSetting->vLFirst;
	fname_tmp = dir_dSetting + "/" + "vLFirst.csv";
	write_Mat_to_CSV(&tmp_mat, fname_tmp.c_str());

	tmp_mat(0,0) = dSetting->vLLast;
	fname_tmp = dir_dSetting + "/" + "vLLast.csv";
	write_Mat_to_CSV(&tmp_mat, fname_tmp.c_str());

	cout << "done!" << endl;
}

void save_fSetting(SpEOPaths *paths, SpEOFusionSetting *fSetting){
	string dir_fSetting = paths->dir_out + "/" + "fSetting";
	cout << "write fusion settings to files in directory: " << endl << "     " << dir_fSetting << " .. ";
	mkdir(dir_fSetting.c_str(), 0777);
	chmod(dir_fSetting.c_str(), S_IRWXU | S_IRWXG | S_IRWXO);
	//SpEOMatrixD tmp_mat = SpEOMatrixD::Constant(1,1,-99999);
	//string fname_tmp="";

	//tmp_mat(0,0) = fSetting->nrmlIm;
	//fname_tmp = dir_fSetting + "/" + "nrmlIm.csv";
	//write_Mat_to_CSV(&tmp_mat, fname_tmp.c_str());

	//tmp_mat(0,0) = fSetting->nrmlDicts;
	//fname_tmp = dir_fSetting + "/" + "nrmlDicts.csv";
	//write_Mat_to_CSV(&tmp_mat, fname_tmp.c_str());

	//tmp_mat(0,0) = fSetting->substrMean;
	//fname_tmp = dir_fSetting + "/" + "substrMean.csv";
	//write_Mat_to_CSV(&tmp_mat, fname_tmp.c_str());

	//tmp_mat(0,0) = fSetting->ImZ_ref_avlbl;
	//fname_tmp = dir_fSetting + "/" + "ImZ_ref_avlbl.csv";
	//write_Mat_to_CSV(&tmp_mat, fname_tmp.c_str());

	//tmp_mat(0,0) = fSetting->Nc;
	//fname_tmp = dir_fSetting + "/" + "Nc.csv";
	//write_Mat_to_CSV(&tmp_mat, fname_tmp.c_str());

	//tmp_mat(0,0) = fSetting->No;
	//fname_tmp = dir_fSetting + "/" + "No.csv";
	//write_Mat_to_CSV(&tmp_mat, fname_tmp.c_str());

	//tmp_mat(0,0) = fSetting->tol_SRF;
	//fname_tmp = dir_fSetting + "/" + "tol_SRF.csv";
	//write_Mat_to_CSV(&tmp_mat, fname_tmp.c_str());

	//tmp_mat(0,0) = fSetting->Nc_max;
	//fname_tmp = dir_fSetting + "/" + "Nc_max.csv";
	//write_Mat_to_CSV(&tmp_mat, fname_tmp.c_str());

	//tmp_mat(0,0) = fSetting->CC_min;
	//fname_tmp = dir_fSetting + "/" + "CC_min.csv";
	//write_Mat_to_CSV(&tmp_mat, fname_tmp.c_str());

	//tmp_mat(0,0) = fSetting->patchsize;
	//fname_tmp = dir_fSetting + "/" + "patchsize.csv";
	//write_Mat_to_CSV(&tmp_mat, fname_tmp.c_str());

	//tmp_mat(0,0) = fSetting->overlap;
	//fname_tmp = dir_fSetting + "/" + "overlap.csv";
	//write_Mat_to_CSV(&tmp_mat, fname_tmp.c_str());

	//tmp_mat(0,0) = fSetting->lambda;
	//fname_tmp = dir_fSetting + "/" + "lambda.csv";
	//write_Mat_to_CSV(&tmp_mat, fname_tmp.c_str());

	//tmp_mat(0,0) = fSetting->NDP;
	//fname_tmp = dir_fSetting + "/" + "NDP.csv";
	//write_Mat_to_CSV(&tmp_mat, fname_tmp.c_str());

	//tmp_mat(0,0) = fSetting->two_step_estimation;
	//fname_tmp = dir_fSetting + "/" + "2step.csv";
	//write_Mat_to_CSV(&tmp_mat, fname_tmp.c_str());

	//tmp_mat(0,0) = fSetting->dictselect;
	//fname_tmp = dir_fSetting + "/" + "dictselect.csv";
	//write_Mat_to_CSV(&tmp_mat, fname_tmp.c_str());

	//tmp_mat(0,0) = fSetting->matrixNorm;
	//fname_tmp = dir_fSetting + "/" + "matrixNorm.csv";
	//write_Mat_to_CSV(&tmp_mat, fname_tmp.c_str());

	//tmp_mat(0,0) = fSetting->addMeanPixelwise;
	//fname_tmp = dir_fSetting + "/" + "addMeanPixelwise.csv";
	//write_Mat_to_CSV(&tmp_mat, fname_tmp.c_str());

	//tmp_mat(0,0) = fSetting->LQ_post_opt;
	//fname_tmp = dir_fSetting + "/" + "LQ_post_opt.csv";
	//write_Mat_to_CSV(&tmp_mat, fname_tmp.c_str());

	//tmp_mat(0,0) = fSetting->lambdaX;
	//fname_tmp = dir_fSetting + "/" + "lambdaX.csv";
	//write_Mat_to_CSV(&tmp_mat, fname_tmp.c_str());

	//tmp_mat(0,0) = fSetting->lambdaY;
	//fname_tmp = dir_fSetting + "/" + "lambdaY.csv";
	//write_Mat_to_CSV(&tmp_mat, fname_tmp.c_str());

	//tmp_mat(0,0) = fSetting->LQ_post_opt_im;
	//fname_tmp = dir_fSetting + "/" + "LQ_post_opt_im.csv";
	//write_Mat_to_CSV(&tmp_mat, fname_tmp.c_str());

	//tmp_mat(0,0) = fSetting->lambdaX_im;
	//fname_tmp = dir_fSetting + "/" + "lambdaX_im.csv";
	//write_Mat_to_CSV(&tmp_mat, fname_tmp.c_str());

	//tmp_mat(0,0) = fSetting->lambdaY_im;
	//fname_tmp = dir_fSetting + "/" + "lambdaY_im.csv";
	//write_Mat_to_CSV(&tmp_mat, fname_tmp.c_str());

	//tmp_mat(0,0) = fSetting->iterMain;
	//fname_tmp = dir_fSetting + "/" + "iterMain.csv";
	//write_Mat_to_CSV(&tmp_mat, fname_tmp.c_str());

	//tmp_mat(0,0) = fSetting->ImZ_init_type;
	//fname_tmp = dir_fSetting + "/" + "ImZ_init_type.csv";
	//write_Mat_to_CSV(&tmp_mat, fname_tmp.c_str());

	//tmp_mat(0,0) = fSetting->doFullImOptWithoutPatRec;
	//fname_tmp = dir_fSetting + "/" + "doFullImOptWithoutPatRec.csv";
	//write_Mat_to_CSV(&tmp_mat, fname_tmp.c_str());

	//tmp_mat(0,0) = fSetting->set_neg_to_0;
	//fname_tmp = dir_fSetting + "/" + "set_neg_to_0.csv";
	//write_Mat_to_CSV(&tmp_mat, fname_tmp.c_str());

	//tmp_mat(0,0) = fSetting->use_estimated_SRFs;
	//fname_tmp = dir_fSetting + "/" + "use_estimated_SRFs.csv";
	//write_Mat_to_CSV(&tmp_mat, fname_tmp.c_str());

	//cout << "done!" << endl;


	//##################################################################
	// write to .m file
	//string fname_matlab;
	string fname_fSetting = dir_fSetting + "/" + "fusionSettings.m";
	std::ofstream fSetting_file_matlab(fname_fSetting.c_str());
	if (fSetting_file_matlab.is_open()) {
		fSetting_file_matlab << setiosflags(ios::fixed);
		//fSetting_file_matlab << setprecision(10);
//		fSetting_file_matlab << "% Cumulative dictionary learning time: "           << endl << "dict_learn_time_cumulative=" << glPrms->timeDictSelect      << ";"<< endl;
//		fSetting_file_matlab << "% Average dictionary learning time per problem: "  << endl << "dict_learn_time_avg="        << glPrms->timeDictSelect_avg  << ";"<< endl;
		fSetting_file_matlab << "  fSet.nrmlIm="<< fSetting->nrmlIm << ";" << endl;
		fSetting_file_matlab << "  fSet.nrmlDicts="<< fSetting->nrmlDicts << ";" << endl;
		fSetting_file_matlab << "  fSet.subtrMean="<< fSetting->substrMean << ";" << endl;
		fSetting_file_matlab << "  fSet.ImZ_ref_avlbl="<< fSetting->ImZ_ref_avlbl << ";" << endl;
		fSetting_file_matlab << "  fSet.Nc="<< fSetting->Nc << ";" << endl;
		fSetting_file_matlab << "  fSet.No="<< fSetting->No << ";" << endl;
		fSetting_file_matlab << "  fSet.tol_SRF="<< fSetting->tol_SRF << ";" << endl;
		fSetting_file_matlab << "  fSet.Nc_max="<< fSetting->Nc_max << ";" << endl;
		fSetting_file_matlab << "  fSet.CC_min="<< fSetting->CC_min << ";" << endl;
		fSetting_file_matlab << "  fSet.patchsize="<< fSetting->patchsize << ";" << endl;
		fSetting_file_matlab << "  fSet.overlap="<< fSetting->overlap << ";" << endl;
		fSetting_file_matlab << "  fSet.lambda="<< fSetting->lambda << ";" << endl;
		fSetting_file_matlab << "  fSet.NDP="<< fSetting->NDP << ";" << endl;
		fSetting_file_matlab << "  fSet.two_step_estimation="<< fSetting->two_step_estimation << ";" << endl;
		fSetting_file_matlab << "  fSet.dictselect="<< fSetting->dictselect << ";" << endl;
		fSetting_file_matlab << "  fSet.matrixNorm="<< fSetting->matrixNorm << ";" << endl;
		fSetting_file_matlab << "  fSet.addMeanPixelwise="<< fSetting->addMeanPixelwise << ";" << endl;
		fSetting_file_matlab << "  fSet.LQ_post_opt="<< fSetting->LQ_post_opt << ";" << endl;
		fSetting_file_matlab << "  fSet.lambdaX="<< fSetting->lambdaX << ";" << endl;
		fSetting_file_matlab << "  fSet.lambdaY="<< fSetting->lambdaY << ";" << endl;
		fSetting_file_matlab << "  fSet.LQ_post_opt_im="<< fSetting->LQ_post_opt_im << ";" << endl;
		fSetting_file_matlab << "  fSet.lambdaX_im="<< fSetting->lambdaX_im << ";" << endl;
		fSetting_file_matlab << "  fSet.lambdaY_im="<< fSetting->lambdaY_im << ";" << endl;
		fSetting_file_matlab << "  fSet.iterMain="<< fSetting->iterMain << ";" << endl;
		fSetting_file_matlab << "  fSet.ImZ_init_type="<< fSetting->ImZ_init_type << ";" << endl;
		fSetting_file_matlab << "  fSet.doFullImOptWithoutPatRec="<< fSetting->doFullImOptWithoutPatRec << ";" << endl;
		fSetting_file_matlab << "  fSet.set_neg_to_0="<< fSetting->set_neg_to_0 << ";" << endl;
		fSetting_file_matlab << "  fSet.use_estimated_SRFs="<< fSetting->use_estimated_SRFs << ";" << endl;
		fSetting_file_matlab << "  fSet.ImX_sim_mode="<< fSetting->ImX_sim_mode << ";" << endl;
		fSetting_file_matlab << "  fSet.use_LRnorm_for_dic_normalization="<< fSetting->use_LRnorm_for_dic_normalization << ";" << endl;

		//fSetting_file_matlab << setprecision(2);

		cout << endl << "fusionSettings.m written." << endl;
	}else{
		cout << endl << "WARNING: statistics_matlab.m could not be written!" << endl;
	}
	fSetting_file_matlab.close();
	chmod(fname_fSetting.c_str(), S_IRWXU | S_IRWXG | S_IRWXO);
	

	//#################################################################
}

void save_oSetting(SpEOPaths *paths, SpEOOutputSetting *oSetting){
	string dir_oSetting = paths->dir_out + "/" + "oSetting";
	cout << "write output settings to files in directory: " << endl << "     " << dir_oSetting << " .. ";
	mkdir(dir_oSetting.c_str(), 0777);
	chmod(dir_oSetting.c_str(), S_IRWXU | S_IRWXG | S_IRWXO);
	SpEOMatrixD tmp_mat = SpEOMatrixD::Constant(1,1,-99999);
	string fname_tmp="";

	tmp_mat(0,0) = oSetting->prec;
	fname_tmp = dir_oSetting + "/" + "precision.csv";
	write_Mat_to_CSV(&tmp_mat, fname_tmp.c_str());

	tmp_mat(0,0) = oSetting->saveAlphas;
	fname_tmp = dir_oSetting + "/" + "saveAlphas.csv";
	write_Mat_to_CSV(&tmp_mat, fname_tmp.c_str());

	tmp_mat(0,0) = oSetting->pFirstAlpha;
	fname_tmp = dir_oSetting + "/" + ".pFirstAlpha.csv";
	write_Mat_to_CSV(&tmp_mat, fname_tmp.c_str());

	tmp_mat(0,0) = oSetting->pLastAlpha;
	fname_tmp = dir_oSetting + "/" + "pLastAlpha.csv";
	write_Mat_to_CSV(&tmp_mat, fname_tmp.c_str());

	tmp_mat(0,0) = oSetting->saveDicts;
	fname_tmp = dir_oSetting + "/" + "saveDicts.csv";
	write_Mat_to_CSV(&tmp_mat, fname_tmp.c_str());

	tmp_mat(0,0) = oSetting->pFirstDict;
	fname_tmp = dir_oSetting + "/" + ".pFirstDict.csv";
	write_Mat_to_CSV(&tmp_mat, fname_tmp.c_str());

	tmp_mat(0,0) = oSetting->pLastDict;
	fname_tmp = dir_oSetting + "/" + "pLastDict.csv";
	write_Mat_to_CSV(&tmp_mat, fname_tmp.c_str());

	cout << "done!" << endl;
}


void save_sSetting(SpEOPaths *paths, SpEOSolverSetting *sSetting){
	string dir_sSetting = paths->dir_out + "/" + "sSetting";
	cout << "write solver settings to files in directory: " << endl << "     " << dir_sSetting << " .. ";
	mkdir(dir_sSetting.c_str(), 0777);
	chmod(dir_sSetting.c_str(), S_IRWXU | S_IRWXG | S_IRWXO);
	SpEOMatrixD tmp_mat = SpEOMatrixD::Constant(1,1,-99999);
	string fname_tmp="";

	tmp_mat(0,0) = sSetting->solver;
	fname_tmp = dir_sSetting + "/" + "solver.csv";
	write_Mat_to_CSV(&tmp_mat, fname_tmp.c_str());

//	tmp_mat(0,0) = sSetting->maxiter_in;
//	fname_tmp = dir_sSetting + "/" + "maxiter_in.csv";
//	write_Mat_to_CSV(&tmp_mat, fname_tmp.c_str());

	tmp_mat(0,0) = sSetting->maxiter_out;
	fname_tmp = dir_sSetting + "/" + "maxiter_out.csv";
	write_Mat_to_CSV(&tmp_mat, fname_tmp.c_str());

	tmp_mat(0,0) = sSetting->tol;
	fname_tmp = dir_sSetting + "/" + "tol.csv";
	write_Mat_to_CSV(&tmp_mat, fname_tmp.c_str());

	tmp_mat(0,0) = sSetting->maxiter_CGLS;
	fname_tmp = dir_sSetting + "/" + "maxiter_CGLS.csv";
	write_Mat_to_CSV(&tmp_mat, fname_tmp.c_str());

	tmp_mat(0,0) = sSetting->tol_r_CGLS;
	fname_tmp = dir_sSetting + "/" + "tol_r_CGLS.csv";
	write_Mat_to_CSV(&tmp_mat, fname_tmp.c_str());

	tmp_mat(0,0) = sSetting->fix_Alpha;
	fname_tmp = dir_sSetting + "/" + "fix_Alpha.csv";
	write_Mat_to_CSV(&tmp_mat, fname_tmp.c_str());

	tmp_mat(0,0) = sSetting->fix_delta_m;
	fname_tmp = dir_sSetting + "/" + "fix_delta_m.csv";
	write_Mat_to_CSV(&tmp_mat, fname_tmp.c_str());

	tmp_mat(0,0) = sSetting->maxiter_CGLS_im;
	fname_tmp = dir_sSetting + "/" + "maxiter_CGLS_im.csv";
	write_Mat_to_CSV(&tmp_mat, fname_tmp.c_str());

	tmp_mat(0,0) = sSetting->tol_r_CGLS_im;
	fname_tmp = dir_sSetting + "/" + "tol_r_CGLS_im.csv";
	write_Mat_to_CSV(&tmp_mat, fname_tmp.c_str());

//	tmp_mat(0,0) = sSetting->weight;
//	fname_tmp = dir_sSetting + "/" + "weight.csv";
//	write_Mat_to_CSV(&tmp_mat, fname_tmp.c_str());

	cout << "done!" << endl;
}


void save_pSetting(SpEOPaths *paths, SpEOParallelSetting *pSetting){
	string dir_pSetting = paths->dir_out + "/" + "pSetting";
	cout << "write parallelization settings to files in directory: " << endl << "     " << dir_pSetting << " .. ";
	mkdir(dir_pSetting.c_str(), 0777);
	chmod(dir_pSetting.c_str(), S_IRWXU | S_IRWXG | S_IRWXO);
	SpEOMatrixD tmp_mat = SpEOMatrixD::Constant(1,1,-99999);
	string fname_tmp="";

	tmp_mat(0,0) = pSetting->numProcGrp;
	fname_tmp = dir_pSetting + "/" + "numProcGrp.csv";
	write_Mat_to_CSV(&tmp_mat, fname_tmp.c_str());

	tmp_mat(0,0) = pSetting->numProcTot;
	fname_tmp = dir_pSetting + "/" + "numProcTot.csv";
	write_Mat_to_CSV(&tmp_mat, fname_tmp.c_str());

	tmp_mat(0,0) = pSetting->store_patches_tmp_on_drive;
	fname_tmp = dir_pSetting + "/" + "store_patches_tmp_on_drive.csv";
	write_Mat_to_CSV(&tmp_mat, fname_tmp.c_str());

	tmp_mat(0,0) = pSetting->workStealingTurns;
	fname_tmp = dir_pSetting + "/" + "workStealingTurns.csv";
	write_Mat_to_CSV(&tmp_mat, fname_tmp.c_str());

	tmp_mat(0,0) = pSetting->numProcPerPatch;
	fname_tmp = dir_pSetting + "/" + "numProcPerPatch.csv";
	write_Mat_to_CSV(&tmp_mat, fname_tmp.c_str());

	cout << "done!" << endl;
}





void save_glPrms(SpEOPaths *paths, SpEOGlobalParams *glPrms){
	string dir_glPrms = paths->dir_out + "/" + "glPrms";
	cout << "write global parameters to files in directory: " << endl << "     " << dir_glPrms << " .. ";
	mkdir(dir_glPrms.c_str(), 0777);
	chmod(dir_glPrms.c_str(), S_IRWXU | S_IRWXG | S_IRWXO);

	SpEOMatrixD tmp_mat = SpEOMatrixD::Constant(1,1,-99999);
	string fname_tmp="";

	tmp_mat(0,0) = glPrms->NChX;
	fname_tmp = dir_glPrms + "/" + "NChX.csv";
	write_Mat_to_CSV(&tmp_mat, fname_tmp.c_str());

	tmp_mat(0,0) = glPrms->NChY;
	fname_tmp = dir_glPrms + "/" + "NChY.csv";
	write_Mat_to_CSV(&tmp_mat, fname_tmp.c_str());

	tmp_mat(0,0) = glPrms->NPU;
	fname_tmp = dir_glPrms + "/" + "NPU.csv";
	write_Mat_to_CSV(&tmp_mat, fname_tmp.c_str());

	tmp_mat(0,0) = glPrms->NPV;
	fname_tmp = dir_glPrms + "/" + "NPV.csv";
	write_Mat_to_CSV(&tmp_mat, fname_tmp.c_str());

	tmp_mat(0,0) = glPrms->NP;
	fname_tmp = dir_glPrms + "/" + "NP.csv";
	write_Mat_to_CSV(&tmp_mat, fname_tmp.c_str());

	tmp_mat(0,0) = glPrms->NPU_sub;
	fname_tmp = dir_glPrms + "/" + "NPU_sub.csv";
	write_Mat_to_CSV(&tmp_mat, fname_tmp.c_str());

	tmp_mat(0,0) = glPrms->NPV_sub;
	fname_tmp = dir_glPrms + "/" + "NPV_sub.csv";
	write_Mat_to_CSV(&tmp_mat, fname_tmp.c_str());

	tmp_mat(0,0) = glPrms->NP_sub;
	fname_tmp = dir_glPrms + "/" + "NP_sub.csv";
	write_Mat_to_CSV(&tmp_mat, fname_tmp.c_str());

	tmp_mat(0,0) = glPrms->uPFirst;
	fname_tmp = dir_glPrms + "/" + "uPFirst.csv";
	write_Mat_to_CSV(&tmp_mat, fname_tmp.c_str());

	tmp_mat(0,0) = glPrms->uPLast;
	fname_tmp = dir_glPrms + "/" + "uPLast.csv";
	write_Mat_to_CSV(&tmp_mat, fname_tmp.c_str());

	tmp_mat(0,0) = glPrms->vPFirst;
	fname_tmp = dir_glPrms + "/" + "vPFirst.csv";
	write_Mat_to_CSV(&tmp_mat, fname_tmp.c_str());

	tmp_mat(0,0) = glPrms->vPLast;
	fname_tmp = dir_glPrms + "/" + "vPLast.csv";
	write_Mat_to_CSV(&tmp_mat, fname_tmp.c_str());

	tmp_mat(0,0) = glPrms->Ng;
	fname_tmp = dir_glPrms + "/" + "Ng.csv";
	write_Mat_to_CSV(&tmp_mat, fname_tmp.c_str());

	tmp_mat(0,0) = glPrms->fDS;
	fname_tmp = dir_glPrms + "/" + "fDS.csv";
	write_Mat_to_CSV(&tmp_mat, fname_tmp.c_str());

//	tmp_mat(0,0) = glPrms->usedPatches;
//	fname_tmp = dir_glPrms + "/" + "usedPatches.csv";
//	write_Mat_to_CSV(&tmp_mat, fname_tmp.c_str());

	tmp_mat(0,0) = glPrms->timeMainLoop;
	fname_tmp = dir_glPrms + "/" + "timeMainLoop.csv";
	write_Mat_to_CSV(&tmp_mat, fname_tmp.c_str());

	tmp_mat(0,0) = glPrms->timeFileWrite;
	fname_tmp = dir_glPrms + "/" + "timeFileWrite.csv";
	write_Mat_to_CSV(&tmp_mat, fname_tmp.c_str());

	tmp_mat(0,0) = glPrms->timeDictSelect;
	fname_tmp = dir_glPrms + "/" + "timeDictSelect.csv";
	write_Mat_to_CSV(&tmp_mat, fname_tmp.c_str());

	tmp_mat(0,0) = glPrms->timeDictSelect_avg;
	fname_tmp = dir_glPrms + "/" + "timeDictSelect_avg.csv";
	write_Mat_to_CSV(&tmp_mat, fname_tmp.c_str());


	//		  tmp_mat(0,0) = glPrms->myChX;
	//		  fname_tmp = dir_glPrms + "/" + "myChX.csv";
	//		  write_Mat_to_CSV(&tmp_mat, fname_tmp.c_str());

	//		  tmp_mat(0,0) = glPrms->myBundle;
	//		  fname_tmp = dir_glPrms + "/" + "myBundle.csv";
	//		  write_Mat_to_CSV(&tmp_mat, fname_tmp.c_str());

	tmp_mat(0,0) = glPrms->numPatchGroups;
	fname_tmp = dir_glPrms + "/" + "numPatchGroups.csv";
	write_Mat_to_CSV(&tmp_mat, fname_tmp.c_str());

	tmp_mat(0,0) = glPrms->timeTotal;
	fname_tmp = dir_glPrms + "/" + "timeTotal.csv";
	write_Mat_to_CSV(&tmp_mat, fname_tmp.c_str());

	tmp_mat(0,0) = glPrms->numProbPerPatch;
	fname_tmp = dir_glPrms + "/" + "numProbPerPatch.csv";
	write_Mat_to_CSV(&tmp_mat, fname_tmp.c_str());

	SpEOMatrixD idxPUH_tmp = glPrms->idxPUH.cast<double>();
	fname_tmp = dir_glPrms + "/" + "idxPUH.csv";
	write_Mat_to_CSV(&idxPUH_tmp, fname_tmp.c_str());

	SpEOMatrixD idxPVH_tmp = glPrms->idxPVH.cast<double>();
	fname_tmp = dir_glPrms + "/" + "idxPVH.csv";
	write_Mat_to_CSV(&idxPUH_tmp, fname_tmp.c_str());

	//		  SpEOMatrixD idxPUL_tmp = glPrms->idxPUL.cast<double>();
	//		  fname_tmp = dir_glPrms + "/" + "idxPUL.csv";
	//		  write_Mat_to_CSV(&(idxPUL_tmp), fname_tmp.c_str());
	//
	//		  SpEOMatrixD idxPVL_tmp = glPrms->idxPVL.cast<double>();
	//		  fname_tmp = dir_glPrms + "/" + "idxPVL.csv";
	//		  write_Mat_to_CSV(&(idxPUL_tmp), fname_tmp.c_str());

	//		  SpEOMatrixD idxPVH_tmp = glPrms->idxPVH.cast<double>();
	fname_tmp = dir_glPrms + "/" + "decMat_C.csv";
	write_Mat_to_CSV(&glPrms->decMat_C, fname_tmp.c_str());

	cout << "done!" << endl;
}




void save_fusion_setup(SpEOPaths *paths, SpEODataIOSetting *dSetting, SpEOFusionSetting *fSetting, SpEOParallelSetting *pSetting, SpEOSolverSetting *sSetting, SpEOOutputSetting *oSetting, SpEOGlobalParams *glPrms){

	//################################
	//# save fusion setup
	//################################
	string fname_fusionSetup = paths->dir_out + "/" + "fusionSetup.m";
	cout << "write fusion setup to file: " << endl << "     " << fname_fusionSetup << " .. ";
	std::ofstream fusionSetup_file_matlab(fname_fusionSetup.c_str());
	if (fusionSetup_file_matlab.is_open()) {
		fusionSetup_file_matlab << setiosflags(ios::fixed);
		//fusionSetup_file_matlab << setprecision(10);
	//	fusionSetup_file_matlab << "  dSet.uLast="<< dSetting->uLLast << ";" << endl;
	//	fusionSetup_file_matlab << "  fSet.use_estimated_SRFs="<< fSetting->use_estimated_SRFs << ";" << endl;
	//	fusionSetup_file_matlab << "  sSet.tol_r_CGLS_im="<< sSetting->tol_r_CGLS_im << ";" << endl;
	//	fusionSetup_file_matlab << "  oSet.writeImageFileAfterEveryIter="<< oSetting->writeImageFileAfterEveryIter << ";" << endl;
	//	fusionSetup_file_matlab << "  pSet.numProcTotal="<< pSetting->numProcTot << ";" << endl;
	//	fusionSetup_file_matlab << "  glPrms.fDS="<< glPrms->fDS << ";" << endl;

		 fusionSetup_file_matlab << "function [fSet,dSet,sSet,oSet,glPrms,paths]=fusionSetup(~)" << endl;
		 fusionSetup_file_matlab << "%% Fusion Settings" << endl 
					 << "fSet.fMethod="<<fSetting->fMethod << ";" << endl
					 << "fSet.two_step_estimation="<<fSetting->two_step_estimation << ";" << endl
					 << "fSet.nrmlIm="<<fSetting->nrmlIm << ";" << endl
					 << "fSet.nrmlDicts="<<fSetting->nrmlDicts << ";" << endl
					 << "fSet.substrMean="<<fSetting->substrMean << ";" << endl
					 << "fSet.ImZ_ref_avlbl="<<fSetting->ImZ_ref_avlbl << ";" << endl
					 << "fSet.Nc="<<fSetting->Nc << ";" << endl
					 << "fSet.No="<<fSetting->No << ";" << endl
					 << "fSet.tol_SRF="<<fSetting->tol_SRF << ";" << endl
					 << "fSet.patchsize="<<fSetting->patchsize << ";" << endl
					 << "fSet.winSize="<<fSetting->winSize << ";" << endl
					 << "fSet.overlap="<<fSetting->overlap << ";" << endl
					 << "fSet.lambda="<<fSetting->lambda << ";" << endl
					 << "fSet.lambdaX="<<fSetting->lambdaX << ";" << endl
					 << "fSet.lambdaY="<<fSetting->lambdaY << ";" << endl
					 << "fSet.lambdaX_im="<<fSetting->lambdaX_im << ";" << endl
					 << "fSet.lambdaY_im="<<fSetting->lambdaY_im << ";" << endl
					 << "fSet.NDP="<<fSetting->NDP << ";" << endl
					 << "fSet.evaluate="<<fSetting->evaluate << ";" << endl
					 << "fSet.evaluate_ImZ_init="<<fSetting->evaluate_ImZ_init << ";" << endl
					 << "fSet.dictselect="<<fSetting->dictselect << ";" << endl
					 << "fSet.matrixNorm="<<fSetting->matrixNorm << ";" << endl
					 << "fSet.addMeanPixelwise="<<fSetting->addMeanPixelwise << ";" << endl
					 << "fSet.LQ_post_opt="<<fSetting->LQ_post_opt << ";" << endl
					 << "fSet.LQ_post_opt_im="<<fSetting->LQ_post_opt_im << ";" << endl
					 << "fSet.useSimulatedImXforDictLearn="<<fSetting->useSimulatedImXforDictLearn << ";" << endl
					 << "fSet.useNewMethodForCalculatingZ="<<fSetting->useNewMethodForCalculatingZ << ";" << endl
					 << "fSet.lambdaX_ABC="<<fSetting->lambdaX_ABC << ";" << endl
					 << "fSet.lambdaY_ABC="<<fSetting->lambdaY_ABC << ";" << endl
					 << "fSet.lambdaZ_ABC="<<fSetting->lambdaZ_ABC << ";" << endl
					 << "fSet.lambdaZ_ABC_in_1st_iter="<<fSetting->lambdaZ_ABC_in_1st_iter << ";" << endl
					 << "fSet.ImZ_init_type="<<fSetting->ImZ_init_type << ";" << endl
					 << "fSet.iterMain="<<fSetting->iterMain << ";" << endl
					 << "fSet.doFullImOptWithoutPatRec="<<fSetting->doFullImOptWithoutPatRec << ";" << endl
					 << "fSet.Nc_max="<<fSetting->Nc_max << ";" << endl
					 << "fSet.CC_min="<<fSetting->CC_min << ";" << endl
					 << "fSet.set_neg_to_0="<<fSetting->set_neg_to_0 << ";" << endl
					 << "fSet.use_estimated_SRFs="<<fSetting->use_estimated_SRFs << ";" << endl
					 << "fSet.fullImOptOnSubspace="<<fSetting->fullImOptOnSubspace << ";" << endl
					 << "fSet.use_init_value_Eq1Unmixing="<<fSetting->use_init_value_Eq1Unmixing << ";" << endl
					 << "fSet.subspace_transform_type=\'"<<fSetting->subspace_transform_type << "\';" << endl
					 << "fSet.subspace_dim="<<fSetting->subspace_dim << ";" << endl
					 << "fSet.ImX_sim_mode="<<fSetting->ImX_sim_mode << ";" << endl
					 << "fSet.SNR_normalization="<<fSetting->SNR_normalization << ";" << endl
					 << "fSet.balance_ImX_term_coef="<<fSetting->balance_ImX_term_coef << ";" << endl
					 << endl;

		 fusionSetup_file_matlab << "%% dataIO Settings" << endl
					 << "dSet.jobName=\'"<<dSetting->jobName << "\';" << endl
					 << "dSet.jobID=\'"<<dSetting->jobID << "\';" << endl
					 << "dSet.chBundleFirst="<<dSetting->chBundleFirst << ";" << endl
					 << "dSet.chBundleLast="<<dSetting->chBundleLast << ";" << endl
					 << "dSet.uLFirst="<<dSetting->uLFirst << ";" << endl
					 << "dSet.uLLast="<<dSetting->uLLast << ";" << endl
					 << "dSet.vLFirst="<<dSetting->vLFirst << ";" << endl
					 << "dSet.vLLast="<<dSetting->vLLast << ";" << endl
					 << "dSet.delete_tmp_patch_folders="<<dSetting->delete_tmp_patch_folders << ";" << endl
					 << "dSet.dir_tmp_patches_additional_num="<<dSetting->dir_tmp_patches_additional_num << ";" << endl
					 << "dSet.imageConstructionOnly="<<dSetting->imageConstructionOnly << ";" << endl
					 << "dSet.contUnfinishedRec="<<dSetting->contUnfinishedRec << ";" << endl
					 << "dSet.platformID="<<dSetting->platformID << ";" << endl
					 << "dSet.saveAsDouble="<<dSetting->saveAsDouble << ";" << endl
					 << endl;
		 fusionSetup_file_matlab << "%% Parrallilization Settings" << endl 
					 << "pSet.numProcTot="<<pSetting->numProcTot << ";" << endl
					 << "pSet.numProcGrp="<<pSetting->numProcGrp << ";" << endl
					 << "pSet.numProcPerPatch="<<pSetting->numProcPerPatch << ";" << endl
					 << "pSet.workStealingTurns="<<pSetting->workStealingTurns << ";" << endl
					 << "pSet.store_patches_tmp_on_drive="<<pSetting->store_patches_tmp_on_drive << ";" << endl
					 << "pSet.parWrNumProc="<<pSetting->parWrNumProc << ";" << endl
					 << endl;
		 fusionSetup_file_matlab << "%% Solver Settings" << endl 
					 << "sSet.solver="<<sSetting->solver << ";" << endl
					 << "sSet.maxiter_out="<<sSetting->maxiter_out << ";" << endl
					 << "sSet.tol="<<sSetting->tol << ";" << endl
					 << "sSet.maxiter_CGLS="<<sSetting->maxiter_CGLS << ";" << endl
					 << "sSet.tol_r_CGLS="<<sSetting->tol_r_CGLS << ";" << endl
					 << "sSet.fix_Alpha="<<sSetting->fix_Alpha << ";" << endl
					 << "sSet.fix_delta_m="<<sSetting->fix_delta_m << ";" << endl
					 << "sSet.maxiter_CGLS_im="<<sSetting->maxiter_CGLS_im << ";" << endl
					 << "sSet.tol_r_CGLS_im="<<sSetting->tol_r_CGLS_im << ";" << endl
					 << endl;
		 fusionSetup_file_matlab << "%% Output Settings" << endl 
					 << "oSet.prec="<<oSetting->prec << ";" << endl
					 << "oSet.saveAlphas="<<oSetting->saveAlphas << ";" << endl
					 << "oSet.pFirstAlpha="<<oSetting->pFirstAlpha << ";" << endl
					 << "oSet.pLastAlpha="<<oSetting->pLastAlpha << ";" << endl
					 << "oSet.saveDicts="<<oSetting->saveDicts << ";" << endl
					 << "oSet.pFirstDict="<<oSetting->pFirstDict << ";" << endl
					 << "oSet.pLastDict="<<oSetting->pLastDict << ";" << endl
					 << "oSet.writeImageFile="<<oSetting->writeImageFile << ";" << endl
					 << "oSet.writeImageFileAfterEveryIter="<<oSetting->writeImageFileAfterEveryIter << ";" << endl
					 << endl;
		 fusionSetup_file_matlab << "%% Global Parameters" << endl 
					 << "glPrms.NChX="<<glPrms->NChX << ";" << endl
					 << "glPrms.NChY="<<glPrms->NChY << ";" << endl
					 << "glPrms.NChY_subspace="<<glPrms->NChY_subspace << ";" << endl
					 << "glPrms.NChX_orig="<<glPrms->NChX_orig << ";" << endl
					 << "glPrms.NChY_orig="<<glPrms->NChY_orig << ";" << endl
					 << "glPrms.NChZ="<<glPrms->NChZ << ";" << endl
					 << "glPrms.NPU="<<glPrms->NPU << ";" << endl
					 << "glPrms.NPV="<<glPrms->NPV << ";" << endl
					 << "glPrms.NPU_sub="<<glPrms->NPU_sub << ";" << endl
					 << "glPrms.NPV_sub="<<glPrms->NPV_sub << ";" << endl
					 << "glPrms.uPFirst="<<glPrms->uPFirst << ";" << endl
					 << "glPrms.uPLast="<<glPrms->uPLast << ";" << endl
					 << "glPrms.vPFirst="<<glPrms->vPFirst << ";" << endl
					 << "glPrms.vPLast="<<glPrms->vPLast << ";" << endl
					 << "glPrms.fDS="<<glPrms->fDS << ";" << endl
					 << "glPrms.NP="<<glPrms->NP << ";" << endl
					 << "glPrms.Ng="<<glPrms->Ng << ";" << endl
					 << "glPrms.timeMainLoop="<<glPrms->timeMainLoop << ";" << endl
					 << "glPrms.timeFileWrite="<<glPrms->timeFileWrite << ";" << endl
					 << "glPrms.timeTotal="<<glPrms->timeTotal << ";" << endl
					 << "glPrms.timeDictSelect="<<glPrms->timeDictSelect << ";" << endl
					 << "glPrms.timeDictSelect_avg="<<glPrms->timeDictSelect_avg << ";" << endl
					 <<  endl;
		 fusionSetup_file_matlab << "%% paths" << endl 
					 << "paths.dir_in=\'"<<paths->dir_in << "\';" << endl
					 << "paths.dataSetID_str=\'"<< paths->dataSetID_str << "\';" <<  endl
					 << endl;

		cout << endl << "fusionSetup.m written." << endl;
	}else{
		cout << endl << "WARNING: fusionSetup.m could not be written!" << endl;
	}
	fusionSetup_file_matlab.close();
	chmod(fname_fusionSetup.c_str(), S_IRWXU | S_IRWXG | S_IRWXO);
	//#################################################################

//	//################################
//	//# save SpEOFusionSetting
//	//################################
//	string dir_fSetting = paths->dir_out + "/" + "fSetting";
//	cout << "write fusion settings to files in directory: " << endl << "     " << dir_fSetting << " .. ";
//	mkdir(dir_fSetting.c_str(), 0777);
//	chmod(dir_fSetting.c_str(), S_IRWXU | S_IRWXG | S_IRWXO);
//	string fname_fSetting = dir_fSetting + "/" + "fusionSettings.m";
//	std::ofstream fSetting_file_matlab(fname_fSetting.c_str());
//	if (fSetting_file_matlab.is_open()) {
//		fSetting_file_matlab << setiosflags(ios::fixed);
//		//fSetting_file_matlab << setprecision(10);
////		fSetting_file_matlab << "% Cumulative dictionary learning time: "           << endl << "dict_learn_time_cumulative=" << glPrms->timeDictSelect      << ";"<< endl;
////		fSetting_file_matlab << "% Average dictionary learning time per problem: "  << endl << "dict_learn_time_avg="        << glPrms->timeDictSelect_avg  << ";"<< endl;
//		fSetting_file_matlab << "  fSet.nrmlIm="<< fSetting->nrmlIm << ";" << endl;
//		fSetting_file_matlab << "  fSet.nrmlDicts="<< fSetting->nrmlDicts << ";" << endl;
//		fSetting_file_matlab << "  fSet.subtrMean="<< fSetting->substrMean << ";" << endl;
//		fSetting_file_matlab << "  fSet.ImZ_ref_avlbl="<< fSetting->ImZ_ref_avlbl << ";" << endl;
//		fSetting_file_matlab << "  fSet.Nc="<< fSetting->Nc << ";" << endl;
//		fSetting_file_matlab << "  fSet.No="<< fSetting->No << ";" << endl;
//		fSetting_file_matlab << "  fSet.tol_SRF="<< fSetting->tol_SRF << ";" << endl;
//		fSetting_file_matlab << "  fSet.Nc_max="<< fSetting->Nc_max << ";" << endl;
//		fSetting_file_matlab << "  fSet.CC_min="<< fSetting->CC_min << ";" << endl;
//		fSetting_file_matlab << "  fSet.patchsize="<< fSetting->patchsize << ";" << endl;
//		fSetting_file_matlab << "  fSet.overlap="<< fSetting->overlap << ";" << endl;
//		fSetting_file_matlab << "  fSet.lambda="<< fSetting->lambda << ";" << endl;
//		fSetting_file_matlab << "  fSet.NDP="<< fSetting->NDP << ";" << endl;
//		fSetting_file_matlab << "  fSet.two_step_estimation="<< fSetting->two_step_estimation << ";" << endl;
//		fSetting_file_matlab << "  fSet.dictselect="<< fSetting->dictselect << ";" << endl;
//		fSetting_file_matlab << "  fSet.matrixNorm="<< fSetting->matrixNorm << ";" << endl;
//		fSetting_file_matlab << "  fSet.addMeanPixelwise="<< fSetting->addMeanPixelwise << ";" << endl;
//		fSetting_file_matlab << "  fSet.LQ_post_opt="<< fSetting->LQ_post_opt << ";" << endl;
//		fSetting_file_matlab << "  fSet.lambdaX="<< fSetting->lambdaX << ";" << endl;
//		fSetting_file_matlab << "  fSet.lambdaY="<< fSetting->lambdaY << ";" << endl;
//		fSetting_file_matlab << "  fSet.LQ_post_opt_im="<< fSetting->LQ_post_opt_im << ";" << endl;
//		fSetting_file_matlab << "  fSet.lambdaX_im="<< fSetting->lambdaX_im << ";" << endl;
//		fSetting_file_matlab << "  fSet.lambdaY_im="<< fSetting->lambdaY_im << ";" << endl;
//		fSetting_file_matlab << "  fSet.iterMain="<< fSetting->iterMain << ";" << endl;
//		fSetting_file_matlab << "  fSet.ImZ_init_type="<< fSetting->ImZ_init_type << ";" << endl;
//		fSetting_file_matlab << "  fSet.doFullImOptWithoutPatRec="<< fSetting->doFullImOptWithoutPatRec << ";" << endl;
//		fSetting_file_matlab << "  fSet.set_neg_to_0="<< fSetting->set_neg_to_0 << ";" << endl;
//		fSetting_file_matlab << "  fSet.use_estimated_SRFs="<< fSetting->use_estimated_SRFs << ";" << endl;
//
//		//fSetting_file_matlab << setprecision(2);
//
//		cout << endl << "fusionSettings.m written." << endl;
//	}else{
//		cout << endl << "WARNING: statistics_matlab.m could not be written!" << endl;
//	}
//	fSetting_file_matlab.close();
//	chmod(fname_fSetting.c_str(), S_IRWXU | S_IRWXG | S_IRWXO);
	

//	//################################
//	//# save SpEODataIOSetting
//	//################################
//	string dir_dSetting = paths->dir_out + "/" + "dSetting";
//	cout << "write dataIO settings to files in directory: " << endl << "     " << dir_dSetting << " .. ";
//	mkdir(dir_dSetting.c_str(), 0777);
//	chmod(dir_dSetting.c_str(), S_IRWXU | S_IRWXG | S_IRWXO);
//	string fname_dSetting = dir_fSetting + "/" + "dataIOSettings.m";
//	std::ofstream dSetting_file_matlab(fname_fSetting.c_str());
//	if (dSetting_file_matlab.is_open()) {
//		dSetting_file_matlab << setiosflags(ios::fixed);
//		//dSetting_file_matlab << setprecision(10);
//		dSetting_file_matlab << "  fSet.use_estimated_SRFs="<< fSetting->use_estimated_SRFs << ";" << endl;
//
//		cout << endl << "dataIOSettings.m written." << endl;
//	}else{
//		cout << endl << "WARNING: statistics_matlab.m could not be written!" << endl;
//	}
//	dSetting_file_matlab.close();
//	chmod(fname_dSetting.c_str(), S_IRWXU | S_IRWXG | S_IRWXO);
//	//#################################################################
//
//	//################################
//	//# save SpEOFusionSetting
//	//################################
//	string dir_fSetting = paths->dir_out + "/" + "fSetting";
//	cout << "write fusion settings to files in directory: " << endl << "     " << dir_fSetting << " .. ";
//	mkdir(dir_fSetting.c_str(), 0777);
//	chmod(dir_fSetting.c_str(), S_IRWXU | S_IRWXG | S_IRWXO);
//	string fname_fSetting = dir_fSetting + "/" + "fusionSettings.m";
//	std::ofstream fSetting_file_matlab(fname_fSetting.c_str());
//	if (fSetting_file_matlab.is_open()) {
//		fSetting_file_matlab << setiosflags(ios::fixed);
//		//fSetting_file_matlab << setprecision(10);
//		fSetting_file_matlab << "  fSet.use_estimated_SRFs="<< fSetting->use_estimated_SRFs << ";" << endl;
//
//		cout << endl << "fusionSettings.m written." << endl;
//	}else{
//		cout << endl << "WARNING: statistics_matlab.m could not be written!" << endl;
//	}
//	fSetting_file_matlab.close();
//	chmod(fname_fSetting.c_str(), S_IRWXU | S_IRWXG | S_IRWXO);
//	//#################################################################
//
//	//################################
//	//# save SpEOFusionSetting
//	//################################
//	string dir_fSetting = paths->dir_out + "/" + "fSetting";
//	cout << "write fusion settings to files in directory: " << endl << "     " << dir_fSetting << " .. ";
//	mkdir(dir_fSetting.c_str(), 0777);
//	chmod(dir_fSetting.c_str(), S_IRWXU | S_IRWXG | S_IRWXO);
//	string fname_fSetting = dir_fSetting + "/" + "fusionSettings.m";
//	std::ofstream fSetting_file_matlab(fname_fSetting.c_str());
//	if (fSetting_file_matlab.is_open()) {
//		fSetting_file_matlab << setiosflags(ios::fixed);
//		//fSetting_file_matlab << setprecision(10);
//		fSetting_file_matlab << "  fSet.use_estimated_SRFs="<< fSetting->use_estimated_SRFs << ";" << endl;
//
//		cout << endl << "fusionSettings.m written." << endl;
//	}else{
//		cout << endl << "WARNING: statistics_matlab.m could not be written!" << endl;
//	}
//	fSetting_file_matlab.close();
//	chmod(fname_fSetting.c_str(), S_IRWXU | S_IRWXG | S_IRWXO);
//	//#################################################################
//
//	//################################
//	//# save SpEOFusionSetting
//	//################################
//	string dir_fSetting = paths->dir_out + "/" + "fSetting";
//	cout << "write fusion settings to files in directory: " << endl << "     " << dir_fSetting << " .. ";
//	mkdir(dir_fSetting.c_str(), 0777);
//	chmod(dir_fSetting.c_str(), S_IRWXU | S_IRWXG | S_IRWXO);
//	string fname_fSetting = dir_fSetting + "/" + "fusionSettings.m";
//	std::ofstream fSetting_file_matlab(fname_fSetting.c_str());
//	if (fSetting_file_matlab.is_open()) {
//		fSetting_file_matlab << setiosflags(ios::fixed);
//		//fSetting_file_matlab << setprecision(10);
//		fSetting_file_matlab << "  fSet.use_estimated_SRFs="<< fSetting->use_estimated_SRFs << ";" << endl;
//
//		cout << endl << "fusionSettings.m written." << endl;
//	}else{
//		cout << endl << "WARNING: statistics_matlab.m could not be written!" << endl;
//	}
//	fSetting_file_matlab.close();
//	chmod(fname_fSetting.c_str(), S_IRWXU | S_IRWXG | S_IRWXO);
//	//#################################################################
//
//	//################################
//	//# save SpEOFusionSetting
//	//################################
//	string dir_fSetting = paths->dir_out + "/" + "fSetting";
//	cout << "write fusion settings to files in directory: " << endl << "     " << dir_fSetting << " .. ";
//	mkdir(dir_fSetting.c_str(), 0777);
//	chmod(dir_fSetting.c_str(), S_IRWXU | S_IRWXG | S_IRWXO);
//	string fname_fSetting = dir_fSetting + "/" + "fusionSettings.m";
//	std::ofstream fSetting_file_matlab(fname_fSetting.c_str());
//	if (fSetting_file_matlab.is_open()) {
//		fSetting_file_matlab << setiosflags(ios::fixed);
//		//fSetting_file_matlab << setprecision(10);
//		fSetting_file_matlab << "  fSet.use_estimated_SRFs="<< fSetting->use_estimated_SRFs << ";" << endl;
//
//		cout << endl << "fusionSettings.m written." << endl;
//	}else{
//		cout << endl << "WARNING: statistics_matlab.m could not be written!" << endl;
//	}
//	fSetting_file_matlab.close();
//	chmod(fname_fSetting.c_str(), S_IRWXU | S_IRWXG | S_IRWXO);
//	//#################################################################
}







void save_evalResults(SpEOPaths *paths, SpEOGlobalParams *glPrms, SpEOFusionSetting *fSetting, SpEOAssessmentMetrics *assMetrics_HR, int iterMain, int numIterMain, bool beforeFullImOpt, bool doFullImOptWithoutPatRec, bool init_image_eval, bool final_evaluation){
//  	  string numStr = std::to_string(iterMain);
//  	  string numStr = lexical_cast<string>(iterMain);
	string dir_eval = paths->dir_out + "/" + "eval";
//	if(init_image_eval){//(iterMain==0 && (beforeFullImOpt || (!beforeFullImOpt && doFullImOptWithoutPatRec))){
//		mkdir(dir_eval.c_str(), 0777);
//		chmod(dir_eval.c_str(), S_IRWXU | S_IRWXG | S_IRWXO);
//	}
	//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
	//if(!final_evaluation){
	//	stringstream numStrSS;
	//	numStrSS << iterMain;
	//	string numStr = numStrSS.str();
	//	if(init_image_eval){
	//		dir_eval=dir_eval + "/" + "iter" + numStr + "_0_init_image";
	//	}else{
	//		if(beforeFullImOpt){
	//			dir_eval=dir_eval + "/" + "iter" + numStr + "_1_after_patch_rec";
	//		}else{
	//			dir_eval=dir_eval + "/" + "iter" + numStr + "_2_after_full_im_opt";
	//		}
	//	}
	//	//mkdir(dir_eval.c_str(), 0777);
	//	//chmod(dir_eval.c_str(), S_IRWXU | S_IRWXG | S_IRWXO);
	//}else{
	//	stringstream numStrSS;
        //        numStrSS << iterMain;
        //        string numStr = numStrSS.str();
	//	if(beforeFullImOpt){
        //                dir_eval=dir_eval + "/" + "iter" + numStr + "_1_after_patch_rec";
        //        }else{
        //                dir_eval=dir_eval + "/" + "iter" + numStr + "_2_after_full_im_opt";
        //        }
	//}
	//=======================================
        stringstream numStrSS;
        numStrSS << iterMain;
        string numStr = numStrSS.str();
        //dir_eval=dir_eval + "/" + "iter" + numStr + ;

	string endOfName;
	if(!final_evaluation){
                if(init_image_eval){
                         endOfName = "iter" + numStr + "_0_init_image";
                }else{
                        if(beforeFullImOpt){
                                 endOfName = "iter" + numStr + "_1_after_patch_rec";
                        }else{
                                 endOfName = "iter" + numStr + "_2_after_full_im_opt";
                        } 
                }
                //mkdir(dir_eval.c_str(), 0777);
                //chmod(dir_eval.c_str(), S_IRWXU | S_IRWXG | S_IRWXO);
        }else{
                if(beforeFullImOpt){
                         endOfName = "iter" + numStr + "_1_after_patch_rec";
                }else{
                        endOfName = "iter" + numStr + "_2_after_full_im_opt";
                }
        }
	dir_eval=dir_eval + "/" + endOfName;
	//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
	cout << "write assessment results to files in directory: " << endl << "     " << dir_eval << " .. ";






	//################################
	//# save SpEOFusionSetting
	//################################
	cout << "write assessment results to files to file" << endl << "     " << dir_eval << " .. ";
	//mkdir(dir_fSetting.c_str(), 0777);
	//chmod(dir_fSetting.c_str(), S_IRWXU | S_IRWXG | S_IRWXO);
	string fname_eval = dir_eval + ".m";
	std::ofstream eval_file_matlab(fname_eval.c_str());
	if (eval_file_matlab.is_open()) {
		eval_file_matlab << setiosflags(ios::fixed);
		eval_file_matlab << setprecision(6);
		eval_file_matlab << "function [e]=" << endOfName << "(~)" << endl;
		eval_file_matlab << "  e.RMSE_mean="<< assMetrics_HR->RMSE_mean << ";" << endl
				 << "  e.PSNR_mean="<< assMetrics_HR->PSNR_mean << ";" << endl
				 << "  e.CC_mean="<< assMetrics_HR->CC_mean << ";" << endl
				 << "  e.ERGAS_mean="<< assMetrics_HR->ERGAS_mean << ";" << endl
				 << "  e.UIQI_mean="<< assMetrics_HR->UIQI_mean << ";" << endl
				 << "  e.DD_mean="<< assMetrics_HR->DD_mean << ";" << endl
				 << "  e.SAM="<< assMetrics_HR->SAM << ";" << endl
				 << "  e.DLambda_mean="<< assMetrics_HR->DLambda_mean << ";" << endl
				 << "  e.AG_orig_mean="<< assMetrics_HR->AG_orig_mean << ";" << endl
				 << "  e.AG_rec_mean="<< assMetrics_HR->AG_rec_mean << ";" << endl;

	        int iChY, iChYU, iChYV;
	        int NChY = glPrms->NChY;
	      
	        eval_file_matlab << "  e.RMSE_sep=[";
	        for(iChY=0; iChY<glPrms->NChY-1; iChY++){
	      	  eval_file_matlab << assMetrics_HR->RMSE_sep[iChY] << ",";
	        }
	        eval_file_matlab << assMetrics_HR->RMSE_sep[NChY-1] << "];" << endl;

	        eval_file_matlab << "  e.PSNR_sep=[";
	        for(iChY=0; iChY<glPrms->NChY-1; iChY++){
	      	  eval_file_matlab << assMetrics_HR->PSNR_sep[iChY] << ",";
	        }
	        eval_file_matlab << assMetrics_HR->PSNR_sep[NChY-1] << "];" << endl;

	        eval_file_matlab << "  e.CC_sep=[";
	        for(iChY=0; iChY<glPrms->NChY-1; iChY++){
	      	  eval_file_matlab << assMetrics_HR->CC_sep[iChY] << ",";
	        }
	        eval_file_matlab << assMetrics_HR->CC_sep[NChY-1] << "];" << endl;

	        eval_file_matlab << "  e.ERGAS_sep=[";
	        for(iChY=0; iChY<glPrms->NChY-1; iChY++){
	      	  eval_file_matlab << assMetrics_HR->ERGAS_sep[iChY] << ",";
	        }
	        eval_file_matlab << assMetrics_HR->ERGAS_sep[NChY-1] << "];" << endl;

	        eval_file_matlab << "  e.UIQI_sep=[";
	        for(iChY=0; iChY<glPrms->NChY-1; iChY++){
	      	  eval_file_matlab << assMetrics_HR->UIQI_sep[iChY] << ",";
	        }
	        eval_file_matlab << assMetrics_HR->UIQI_sep[NChY-1] << "];" << endl;

	        eval_file_matlab << "  e.DD_sep=[";
	        for(iChY=0; iChY<glPrms->NChY-1; iChY++){
	      	  eval_file_matlab << assMetrics_HR->DD_sep[iChY] << ",";
	        }
	        eval_file_matlab << assMetrics_HR->DD_sep[NChY-1] << "];" << endl;

		/*
	        eval_file_matlab << "  e.AG_orig_sep=[";
	        for(iChY=0; iChY<glPrms->NChY-1; iChY++){
	      	  eval_file_matlab << assMetrics_HR->AG_orig_sep[iChY] << ",";
	        }
	        eval_file_matlab << assMetrics_HR->AG_orig_sep[NChY-1] << "];" << endl;

	        eval_file_matlab << "  e.AG_rec_sep=[";
	        for(iChY=0; iChY<glPrms->NChY-1; iChY++){
	      	  eval_file_matlab << assMetrics_HR->AG_rec_sep[iChY] << ",";
	        }
	        eval_file_matlab << assMetrics_HR->AG_rec_sep[NChY-1] << "];" << endl;

		eval_file_matlab << setprecision(3);
	        eval_file_matlab << "  e.DLambda_mat=[" << endl;
	        for(iChYU=0; iChYU<glPrms->NChY-1; iChYU++){
	       	 for(iChYV=0; iChYV<glPrms->NChY-1; iChYV++){
	       	         eval_file_matlab << assMetrics_HR->DLambda_mat[iChYU][iChYV] << ",";
	       	 }
	       	 eval_file_matlab << assMetrics_HR->DLambda_mat[iChYU][NChY-1] << ";" << endl;
	        }
	        for(iChYV=0; iChYV<glPrms->NChY-1; iChYV++){
                        eval_file_matlab << assMetrics_HR->DLambda_mat[NChY-1][iChYV] << ",";
                }
                eval_file_matlab << assMetrics_HR->DLambda_mat[NChY-1][NChY-1] << "];" << endl;
		*/

		//eval_file_matlab << endl << "end";

	      	cout << endl << "evaluation file written." << endl;
	}else{
		cout << endl << "WARNING: statistics_matlab.m could not be written!" << endl;
	}
	eval_file_matlab.close();
	chmod(fname_eval.c_str(), S_IRWXU | S_IRWXG | S_IRWXO);
	//#################################################################







//	  SpEOMatrixD tmp_mat = SpEOMatrixD::Constant(1,1,-99999);
//	  string fname_tmp="";
//
//	  // mean values:
//	  tmp_mat(0,0) = assMetrics_HR->RMSE_mean;
//	  fname_tmp = dir_eval + "/" + "RMSE_mean.csv";
//	  write_Mat_to_CSV(&tmp_mat, fname_tmp.c_str());
//
//	  tmp_mat(0,0) = assMetrics_HR->PSNR_mean;
//	  fname_tmp = dir_eval + "/" + "PSNR_mean.csv";
//	  write_Mat_to_CSV(&tmp_mat, fname_tmp.c_str());
//
//	  tmp_mat(0,0) = assMetrics_HR->CC_mean;
//	  fname_tmp = dir_eval + "/" + "CC_mean.csv";
//	  write_Mat_to_CSV(&tmp_mat, fname_tmp.c_str());
//
//	  tmp_mat(0,0) = assMetrics_HR->ERGAS_mean;
//	  fname_tmp = dir_eval + "/" + "ERGAS_mean.csv";
//	  write_Mat_to_CSV(&tmp_mat, fname_tmp.c_str());
//
//	  tmp_mat(0,0) = assMetrics_HR->UIQI_mean;
//	  fname_tmp = dir_eval + "/" + "UIQI_mean.csv";
//	  write_Mat_to_CSV(&tmp_mat, fname_tmp.c_str());
//
//	  tmp_mat(0,0) = assMetrics_HR->DD_mean;
//	  fname_tmp = dir_eval + "/" + "DD_mean.csv";
//	  write_Mat_to_CSV(&tmp_mat, fname_tmp.c_str());
//
//	  tmp_mat(0,0) = assMetrics_HR->SAM;
//	  fname_tmp = dir_eval + "/" + "SAM.csv";
//	  write_Mat_to_CSV(&tmp_mat, fname_tmp.c_str());
//
//	  tmp_mat(0,0) = assMetrics_HR->DLambda_mean;
//	  fname_tmp = dir_eval + "/" + "DLambda_mean.csv";
//	  write_Mat_to_CSV(&tmp_mat, fname_tmp.c_str());
//
//	  tmp_mat(0,0) = assMetrics_HR->AG_orig_mean;
//	  fname_tmp = dir_eval + "/" + "AG_orig_mean.csv";
//	  write_Mat_to_CSV(&tmp_mat, fname_tmp.c_str());
//
//	  tmp_mat(0,0) = assMetrics_HR->AG_rec_mean;
//	  fname_tmp = dir_eval + "/" + "AG_rec_mean.csv";
//	  write_Mat_to_CSV(&tmp_mat, fname_tmp.c_str());
//
//	  // band-wise values:
//
//	  SpEOMatrixD tmp_eval_sep = SpEOMatrixD::Constant(glPrms->NChY,1,-9999);
//	  SpEOMatrixD tmp_eval_mat = SpEOMatrixD::Constant(glPrms->NChY,glPrms->NChY,-9999);
//
//	  int iChY, jChY;
//
//	  for(iChY=0; iChY<glPrms->NChY; iChY++){
//		  tmp_eval_sep(iChY,0) = assMetrics_HR->RMSE_sep[iChY];
//	  }
//	  fname_tmp = dir_eval + "/" + "RMSE_sep.csv";
//	  write_Mat_to_CSV(&tmp_eval_sep, fname_tmp.c_str());
//
//	  for(iChY=0; iChY<glPrms->NChY; iChY++){
//		  tmp_eval_sep(iChY,0) = assMetrics_HR->PSNR_sep[iChY];
//	  }
//	  fname_tmp = dir_eval + "/" + "PSNR_sep.csv";
//	  write_Mat_to_CSV(&tmp_eval_sep, fname_tmp.c_str());
//
//	  for(iChY=0; iChY<glPrms->NChY; iChY++){
//		  tmp_eval_sep(iChY,0) = assMetrics_HR->CC_sep[iChY];
//	  }
//	  fname_tmp = dir_eval + "/" + "CC_sep.csv";
//	  write_Mat_to_CSV(&tmp_eval_sep, fname_tmp.c_str());
//
//	  for(iChY=0; iChY<glPrms->NChY; iChY++){
//		  tmp_eval_sep(iChY,0) = assMetrics_HR->ERGAS_sep[iChY];
//	  }
//	  fname_tmp = dir_eval + "/" + "ERGAS_sep.csv";
//	  write_Mat_to_CSV(&tmp_eval_sep, fname_tmp.c_str());
//
//	  for(iChY=0; iChY<glPrms->NChY; iChY++){
//		  tmp_eval_sep(iChY,0) = assMetrics_HR->UIQI_sep[iChY];
//	  }
//	  fname_tmp = dir_eval + "/" + "UIQI_sep.csv";
//	  write_Mat_to_CSV(&tmp_eval_sep, fname_tmp.c_str());
//
//	  for(iChY=0; iChY<glPrms->NChY; iChY++){
//		  for(jChY=0; jChY<glPrms->NChY; jChY++){
//			  tmp_eval_mat(iChY,jChY) = assMetrics_HR->DLambda_mat[iChY][jChY];
//		  }
//	  }
//	  fname_tmp = dir_eval + "/" + "DLambda_mat.csv";
//	  write_Mat_to_CSV(&tmp_eval_mat, fname_tmp.c_str());
//
//	  for(iChY=0; iChY<glPrms->NChY; iChY++){
//		  tmp_eval_sep(iChY,0) = assMetrics_HR->DD_sep[iChY];
//	  }
//	  fname_tmp = dir_eval + "/" + "DD_sep.csv";
//	  write_Mat_to_CSV(&tmp_eval_sep, fname_tmp.c_str());
//
//	  for(iChY=0; iChY<glPrms->NChY; iChY++){
//		  tmp_eval_sep(iChY,0) = assMetrics_HR->AG_orig_sep[iChY];
//	  }
//	  fname_tmp = dir_eval + "/" + "AG_orig_sep.csv";
//	  write_Mat_to_CSV(&tmp_eval_sep, fname_tmp.c_str());
//
//	  for(iChY=0; iChY<glPrms->NChY; iChY++){
//		  tmp_eval_sep(iChY,0) = assMetrics_HR->AG_rec_sep[iChY];
//	  }
//	  fname_tmp = dir_eval + "/" + "AG_rec_sep.csv";
//	  write_Mat_to_CSV(&tmp_eval_sep, fname_tmp.c_str());
//
////	  if(iterMain ==  numIterMain-1 && !beforeFullImOpt){
//////		  string dir_eval_main = paths->dir_out + "/" + "eval";
////		  ifstream in(dir_eval);
////		  ofstream out(paths->dir_out);
////		  out << in.rdbuf();
////	  }
//
	  cout << "done!" << endl;
}

// for controlled MPI I/O error handling.
void sample_error(int error, char *string)
{
	fprintf(stderr, "Error %d in %s\n", error, string);
	MPI_Finalize();
	exit(-1);
}




// ====================================
// Dictionary Selection Function
// ====================================
/* Select dictionary atoms on the basis of Euclidean Norm comparison for ImX_LR, ImX, and ImY patches
 */

void dictSelectFunc(SpEOMatrixF *patchCompSubsetfloat, SpEODataset *ImX_LR, SpEODataset *ImX, SpEODataset *ImY, SpEOFusionSetting *fSet, SpEOGlobalParams *glPrms, int NP, int iP, int uP, int vP, SpEOVectorI *idxPUH, SpEOVectorI *idxPVH, SpEOVectorI *idxPUL, SpEOVectorI *idxPVL, SpEOMatrixD *SRF, SpEOMatrixF *patchComp, int my_rank, int ipp, int &NDP){	// Current Patch and Comparison Patch
	int pszH = fSet->patchsize*glPrms->fDS;
	int pszL = fSet->patchsize;
	int pszL2 = pszL*pszL;
	int pszH2 = pszH*pszH;
	int NPV = glPrms->NPV;
	bool sortON;

//	// FLOAT
	SpEOMatrixF currenttmp = SpEOMatrixF::Zero(pszL,pszL);
	SpEOMatrixF currP = SpEOMatrixF::Zero(pszL,pszL);
	SpEOMatrixF compP = SpEOMatrixF::Zero(pszL,pszL);
	// Current Patch and Comparison Patch in vector format
	SpEOVectorF currV, compV;
	// Patch indices and norms: SpEOMatrixF patchComp(NP,3);
	// Ordered patch indices and norms
	SpEOMatrixF patchComp2; // Necessary for dictselect=4 Do Not Delete!
	SpEOMatrixF patchCompOrdered(NP,3);
	SpEOMatrixF patchCompSubset(NDP,3);

    // DOUBLE
//	SpEOMatrixD currenttmp = SpEOMatrixD::Zero(pszL,pszL);
//	SpEOMatrixD currP = SpEOMatrixD::Zero(pszL,pszL);
//	SpEOMatrixD compP = SpEOMatrixD::Zero(pszL,pszL);
//	// Current Patch and Comparison Patch in vector format
//	SpEOVectorF currV, compV;
//	// Patch indices and norms: SpEOMatrixD patchComp(NP,3);
//	// Ordered patch indices and norms
//	SpEOMatrixD patchComp2; // Necessary for dictselect=4 Do Not Delete!
//	SpEOMatrixD patchCompOrdered(NP,3);
//	SpEOMatrixD patchCompSubset(NDP,3);


	int iG;
	if(fSet->useNewMethodForCalculatingZ){
		iG= glPrms->myChX[ipp];
	}else{
		iG= ipp;
	}

	switch(fSet->dictselect){

	case 0:{// Case (0): dictionary contains only current patch under reconstruction
			sortON = false;
			(*patchCompSubsetfloat)(0,0) = uP;
			(*patchCompSubsetfloat)(0,1) = vP;
			break;
	}
	case 1:{// Case (1): Nearest Neighbor patches according to the maximum norm of the patch coordinates
		sortON = false;
		(*patchCompSubsetfloat)(0,0) = uP;
		(*patchCompSubsetfloat)(0,1) = vP;
		int pntCnt = 0;
		// (du, dv) is a vector - direction in which we move right now
		int du = 1;
		int dv = 0;
		int segment_length = 1;
		// current position (u, v) and how much of current segment we passed
		int u = uP;
		int v = vP;
		int segment_passed = 0;
		int k = 1;
		int tot_patch_cnt = 0;
		while (k < NDP && tot_patch_cnt<glPrms->NP_sub) {
			// make a step, add 'direction' vector (du, dv) to current position (u, v)
			u += du;
			v += dv;
			++segment_passed;
			if (segment_passed == segment_length) {
				// done with current segment
				segment_passed = 0;
				// 'rotate' directions
				int buffer = du;
				du = -dv;
				dv = buffer;
				// increase segment length if necessary
				if (dv == 0) {
					++segment_length;
				}
			}

			if (   u >= 0 && u <= glPrms->NPU-1
				&& v >= 0 && v <= glPrms->NPV-1) {
				tot_patch_cnt ++;
				bool patch_of_zeros_found = false;
				for(int ig=0; ig<ImX->get_NCh(); ig++){
//					cout << "jhfidsa 1 ig=" << ig << endl;
					SpEOMatrixD patchHR_tmp = (ImX->get_rasterBands()[ig]->bandDataMatD.block(idxPUH->coeff(u), idxPVH->coeff(v), pszH, pszH));
//					cout << "jhfidsa 2 ig=" << ig << endl;
					SpEOMatrixD patchLR_tmp = (ImX_LR->get_rasterBands()[ig]->bandDataMatD.block(idxPUL->coeff(u), idxPVL->coeff(v), pszL, pszL));
					if(patchHR_tmp.norm()<1e-8 || patchLR_tmp.norm()<1e-8){
						patch_of_zeros_found = true;
						cout << "In either ImX_sim or in ImX_sim_LR: patch (uP,vP)=("<<uP<<","<<vP<<") inf band/group ig="<<ig<<" contains only zeros and was therefore excluded from the dictionary!" << endl;
					}
				}
				if(!patch_of_zeros_found){
					(*patchCompSubsetfloat)(k,0) = u;
					(*patchCompSubsetfloat)(k,1) = v;
					pntCnt++;
					k++;
				}
			}
		}
		if(k < NDP-1){
			NDP = k+1;
		}

		break;

	}

	case 2: {// Case (2): Compute norms based on |ImX_LRcurrent-ImX_LR|
		// Set the currP as the current patch in the ImX_LR image, for norm comparison.
		sortON = true;
//		int bp_cnt=0;
		//cout << "bp" << ++bp_cnt << endl;

//		if(fSet->useNewMethodForCalculatingZ){
		currP = ImX_LR->get_rasterBands()[iG]->get_bandDataMat()->block((*idxPUL)(uP), (*idxPVL)(vP), pszL, pszL);
//		currP = ImX_LR->get_rasterBands()[iG]->bandDataMatD.block(idxPUL->coeff(uP), idxPVL->coeff(vP), pszL, pszL);



		// Arrange current patch representation, currP, into a vector, currV.
		//cout << "bp" << ++bp_cnt << endl;
		currV = SpEOVectorF::Map(currP.data(), pszL2);
		//cout << "bp" << ++bp_cnt << endl;
		// Subtract mean and normalise
		currV.array() -= currV.mean();
		//cout << "bp" << ++bp_cnt << endl;
		float norm_tmp = currV.norm();
		//cout << "bp" << ++bp_cnt << endl;
		currV /= norm_tmp;
		//cout << "bp" << ++bp_cnt << endl;
		// Metric comparison
		for(int iN = 0; iN<NP; iN++){
			// Case of LR norm comparison
			// Get Indices for all NP patches
			int uN = iN / NPV;
			int vN = iN % NPV;
			//cout << "bp" << ++bp_cnt << "iN=" << iN << "_1" << endl;
			compP = ImX_LR->get_rasterBands()[iG]->get_bandDataMat()->block((*idxPUL)(uN), (*idxPVL)(vN), pszL, pszL);
//			compP = ImX_LR->get_rasterBands()[iG]->bandDataMatD.block(idxPUL->coeff(uN), idxPVL->coeff(vN), pszL, pszL);
			//cout << "bp" << ++bp_cnt << "iN=" << iN << "_2" << endl;
			compV = SpEOVectorF::Map(compP.data(), pszL2);
			//cout << "bp" << ++bp_cnt << "iN=" << iN << "_3" << endl;
			// Subtract mean and normalise
			compV.array() -= compV.mean();
			float norm_tmp = compV.norm();
			compV /= norm_tmp;
			//cout << "bp" << ++bp_cnt << "iN=" << iN << "_4" << endl;
			(*patchComp)(iN,2) = (currV-compV).norm();
			//cout << "bp" << ++bp_cnt << "iN=" << iN << "_5" << endl;
		}

		break;
	}

	case 3: {// Case (3): Compute norms based on |ImX_approx-ImX|, ImX_approx = sum_Nc(ImY_channel*SRF_channel)
		// Use sum of ImY_channel*SRF_channel to form LR reconstruction as the currP
		sortON = true;
		for(int channel = 0; channel<fSet->Nc; channel++){
			currenttmp = (ImY->get_rasterBands()[channel]->get_bandDataMat()->block((*idxPUL)(uP), (*idxPVL)(vP), pszL, pszL) )*(*SRF)(0,channel);
			currP = currP + currenttmp;
		}
		// Arrange current patch representation, currP, into a vector.
		currV = SpEOVectorF::Map(currP.data(), pszL2);
		// Subtract mean and normalise
		currV.array() -= currV.mean();
		float norm_tmp = currV.norm();
		currV /= norm_tmp;


		// Metric comparison
		for(int iN = 0; iN<NP; iN++){
			// Case of LR norm comparison
			// Get Indices for all NP patches
			int uN = iN / NPV;
			int vN = iN % NPV;
			compP = ImX_LR->get_rasterBands()[iG]->get_bandDataMat()->block((*idxPUL)(uN), (*idxPVL)(vN), pszL, pszL);
			compV = SpEOVectorF::Map(compP.data(), pszL2);
			// Subtract mean and normalise
			compV.array() -= compV.mean();
			float norm_tmp = compV.norm();
			compV /= norm_tmp;
			(*patchComp)(iN,2) = (currV-compV).norm();
		}

		break;
	}

	case 4: {// Case(4): HR norm comparison
		// Set the currP as the current patch in the ImX_HR image, for norm comparison.
		sortON = true;
		currP.resize(pszH,pszH);
		compP.resize(pszH,pszH);
		currP = ImX->get_rasterBands()[iG]->get_bandDataMat()->block((*idxPUH)(uP), (*idxPVH)(vP), pszH, pszH);
		// Arrange current patch representation, currP, into a vector.
		currV = SpEOVectorF::Map(currP.data(), pszH2);
		// Subtract mean and normalise
		currV.array() -= currV.mean();
		float norm_tmp = currV.norm();
		currV /= norm_tmp;

		for(int iN = 0; iN<NP; iN++){
			int uN = iN / NPV;
			int vN = iN % NPV;
			// Case of HR Correlation comparison (less efficient than HR norm Comparison)
			// Get Indices for all NP patches
			compP = ImX->get_rasterBands()[iG]->get_bandDataMat()->block((*idxPUH)(uN), (*idxPVH)(vN), pszH, pszH);
			compV = SpEOVectorF::Map(compP.data(), pszH2);
			// Subtract mean and normalise
			compV.array() -= compV.mean();
			float norm_tmp = compV.norm();
			compV /= norm_tmp;
			(*patchComp)(iN,2) = (currV-compV).norm();
		}

		break;
	}
	case 5:{// Case(5): HR & LR Norm combined ranking !!! Inefficient !!!
		sortON = true;


		//################################   LR   ##################################
		// Set the currP as the current patch in the ImX_LR image, for norm comparison.
//		cout<<"iP??? "<<iP<<endl;
//		cout<<"uP??? "<<uP<<endl;
//		cout<<"vP??? "<<vP<<endl;
//		cout<<"patchComp.row(0) "<<patchComp->row(0)<<endl;




		currP = ImX_LR->get_rasterBands()[iG]->get_bandDataMat()->block((*idxPUL)(uP), (*idxPVL)(vP), pszL, pszL);
		// Arrange current patch representation, currP, into a vector.
		currV = SpEOVectorF::Map(currP.data(), pszL2);
		// Subtract mean and normalise
		currV.array() -= currV.mean();
		float norm_tmp = currV.norm();
		currV /= norm_tmp;
		// Case of LR norm comparison
		// Get Indices for all NP patches
		SpEOMatrixF LRComp=(*patchComp);
		SpEOMatrixF HRComp=(*patchComp);
		SpEOMatrixF LRCompOrdered(NP,3);
		for(int iN = 0; iN<NP; iN++){
			int uN = iN / NPV;
			int vN = iN % NPV;
			compP = ImX_LR->get_rasterBands()[iG]->get_bandDataMat()->block((*idxPUL)(uN), (*idxPVL)(vN), pszL, pszL);
			compV = SpEOVectorF::Map(compP.data(), pszL2);
			compV.array() -= compV.mean();
			float norm_tmp = compV.norm();
			compV /= norm_tmp;

			LRComp(iN,2) = (currV-compV).norm();
		}

		// Order indices according to norm metric
		VectorXi sort_orderLR(NP);
		std::vector<argsort_pair> dataLR(NP);
		for(int i=0;i<NP;i++) {
			dataLR[i].first = i;
			dataLR[i].second = LRComp.col(2)(i);
		}
		// Sort data matrix according to the norms
		std::sort(dataLR.begin(), dataLR.end(), argsort_comp);
		// Get permutations ---> sort_order vector
		for(int i=0;i< NP;i++) {
			sort_orderLR(dataLR[i].first) = i;
		}

		// Generate matrix of indices and norms, in ascending order according to the Norm column
		LRCompOrdered = sort_orderLR.asPermutation()*LRComp;
		// In case the user wants to check the norms and ordering of patches are appropriate print the lines below:
//				cout<<"LRCompOrdered"<<endl<<LRCompOrdered.block(0,0,5,3)<<endl;
//				cout<<"---------"<<endl;

		//################################   HR   ##################################
		// Set the currP as the current patch in the ImX_HR image, for norm comparison.
		currP.resize(pszH,pszH);
		compP.resize(pszH,pszH);
		currP = ImX->get_rasterBands()[iG]->get_bandDataMat()->block((*idxPUH)(uP), (*idxPVH)(vP), pszH, pszH);
		// Arrange current patch representation, currP, into a vector.
		currV = SpEOVectorF::Map(currP.data(), pszH2);
		// Subtract mean and normalise
		currV.array() -= currV.mean();
		norm_tmp = currV.norm();
		currV /= norm_tmp;
		// Case of HR norm comparison
		// Get Indices for all NP patches
		SpEOMatrixF HRCompOrdered(NP,3);
		for(int iN = 0; iN<NP; iN++){
			int uN = iN / NPV;
			int vN = iN % NPV;
			compP = ImX->get_rasterBands()[iG]->get_bandDataMat()->block((*idxPUH)(uN), (*idxPVH)(vN), pszH, pszH);
			compV = SpEOVectorF::Map(compP.data(), pszH2);
			compV.array() -= compV.mean();
			float norm_tmp = compV.norm();
			compV /= norm_tmp;
			HRComp(iN,2) = (currV-compV).norm();
		}

		// Order indices according to norm metric
		VectorXi sort_orderHR(NP);
		std::vector<argsort_pair> dataHR(NP);
		for(int i=0;i<NP;i++) {
			dataHR[i].first = i;
			dataHR[i].second = HRComp.col(2)(i);
		}
		// Sort data matrix according to the norms
		std::sort(dataHR.begin(), dataHR.end(), argsort_comp);
		// Get permutations ---> sort_order vector
		for(int i=0;i< NP;i++) {
			sort_orderHR(dataHR[i].first) = i;
		}

		// Generate matrix of indices and norms, in ascending order according to the Norm column
		HRCompOrdered = sort_orderHR.asPermutation()*HRComp;
		// In case the user wants to check the norms and ordering of patches are appropriate print the lines below:
//				cout<<"HRCompOrdered"<<endl<<HRCompOrdered.block(0,0,5,3)<<endl;
//				cout<<"---------"<<endl;


		// To avoid cross referenceing use patchComp2 (IMPORTANT - Don't delete)
		patchComp2 = LRCompOrdered;

		for(int i=0;i<NP;i++){
			// Compute combined ranking = sum of the ordered HR and LR row positions for each patch.
			for(int j=0;j<NP;j++){
				if(LRCompOrdered.row(i).col(0) == HRCompOrdered.row(j).col(0) && LRCompOrdered.row(i).col(1) == HRCompOrdered.row(j).col(1) ){
					// Combine HR an LR rankings
					patchComp2(i,2) =  i+j;
				}

			}
		}

		break;
	}
	case 6:{ // Case (6): Compute metric based on ABSOLUTE correlation between ImX_HRcurrent and ImX_HR
		// Set the currP as the current patch in the ImX_HR image, for norm comparison.
		sortON = true;
		currP.resize(pszH,pszH);
		compP.resize(pszH,pszH);
		currP = ImX->get_rasterBands()[iG]->get_bandDataMat()->block((*idxPUH)(uP), (*idxPVH)(vP), pszH, pszH);
		// Arrange current patch representation, currP, into a vector.
		currV = SpEOVectorF::Map(currP.data(), pszH2);
		// Subtract mean and normalise
		currV.array() -= currV.mean();
		float norm_tmp = currV.norm();
		currV /= norm_tmp;

		for(int iN = 0; iN<NP; iN++){
			int uN = iN / NPV;
			int vN = iN % NPV;
			// Case of HR Correlation comparison
			// Get Indices for all NP patches
			compP = ImX->get_rasterBands()[iG]->get_bandDataMat()->block((*idxPUH)(uN), (*idxPVH)(vN), pszH, pszH);
			compV = SpEOVectorF::Map(compP.data(), pszH2);
			// Subtract mean and normalise
			compV.array() -= compV.mean();
			float norm_tmp = compV.norm();
			compV /= norm_tmp;
			(*patchComp)(iN,2) = -abs(currV.dot(compV));
		}

		break;
	}
	case 7: {// Case (7): Compute metric based on **MINIMUM** ABSOLUTE correlation between ImX_HRcurrent and ImX_HR
		// Set the currP as the current patch in the ImX_HR image, for norm comparison.
		sortON = true;
		currP.resize(pszH,pszH);
		compP.resize(pszH,pszH);
		currP = ImX->get_rasterBands()[iG]->get_bandDataMat()->block((*idxPUH)(uP), (*idxPVH)(vP), pszH, pszH);
		// Arrange current patch representation, currP, into a vector.
		currV = SpEOVectorF::Map(currP.data(), pszH2);
		// Subtract mean and normalise
		currV.array() -= currV.mean();
		float norm_tmp = currV.norm();
		currV /= norm_tmp;

		for(int iN = 0; iN<NP; iN++){
			int uN = iN / NPV;
			int vN = iN % NPV;
			// Case of HR Correlation comparison
			// Get Indices for all NP patches
			compP = ImX->get_rasterBands()[iG]->get_bandDataMat()->block((*idxPUH)(uN), (*idxPVH)(vN), pszH, pszH);
			compV = SpEOVectorF::Map(compP.data(), pszH2);
			// Subtract mean and normalise
			compV.array() -= compV.mean();
			float norm_tmp = compV.norm();
			compV /= norm_tmp;

			// Insure current patch, iP, is the first element in the dictionary
			if(iN==iP){
				(*patchComp)(iN,2) = 0;
			}
			else{
				(*patchComp)(iN,2) = abs(currV.dot(compV));
			}
		}
		break;
	}
	case 8: {// Case (8): Compute metric based on random patch subset, though retaining the current patch as first dictionary atom
		// Set current patch as the first dictionary patch
		sortON = false;
		(*patchCompSubsetfloat)(0,0) = uP;
		(*patchCompSubsetfloat)(0,1) = vP;
		(*patchCompSubsetfloat)(0,2) = 0;

		// Loop through the remaining dictionary atoms and select random patches
		for(int iD = 1; iD<NDP; iD++){


			// Compute pseudo hashing function for seed (this could be improved by including processor rank).
			int time = MPI_Wtime();
			//int mod = 104729+iP;
			//int seed = abs(((time*181)*((iD-83)*iP*359))) % mod; // TO BE MODIFIED
			int seed = abs(time*121 - my_rank*12345 + iP*iD - uP +vP + iG*79) % 104729;
			// Seed random number
			srand(seed);
			// Compute random index
			int iR = rand() % NP;


			/* requires
			 * 1) to include -std=c++11 in the compiler options: e.g. CFLAGS   = -O3 -w -std=c++11
			 * 2) to uncomment #include <random> in the includes.h file
			int min = 0;
			int max = NP-1;
			random_device rd;     // only used once to initialise (seed) engine
			mt19937 rng(rd());    // random-number engine used (Mersenne-Twister in this case)
			uniform_int_distribution<int> uni(min,max); // guaranteed unbiased
			auto random_integer = uni(rng);
			int iR = random_integer;
			*/


			int uR = iR / NPV;
			int vR = iR % NPV;

			(*patchCompSubsetfloat)(iD,0) = uR;
			(*patchCompSubsetfloat)(iD,1) = vR;
			(*patchCompSubsetfloat)(iD,2) = iR;
		}

		break;
	}
	case 9: {// Case (9): Uncorrelated dictionary - i.e dictionary  attempts not to self correlate, to approximate a basis
		// WARNING: Case 9 involves an extremely expensive combinatorial search - NOT RECOMMENDED
		sortON = false;
		currP.resize(pszH,pszH);
		compP.resize(pszH,pszH);

		// First patch in the dictionary is the current patch
		patchCompOrdered.col(0).row(0) << uP;
		patchCompOrdered.col(1).row(0) << vP;
		patchCompOrdered.col(2).row(0) << 0;
		VectorXi selectedPatches(NDP);
		selectedPatches.row(0) << iP;

		// For the remainder atoms of the dictionary
		for(int iD = 1; iD<NDP; iD++){
			// sum: net summation of correlation coefficients between patch iN and current patches in the dictionary
			float sum = 0;
			// cnt: counter for the selfcorrelation matrix
			int cnt = 0;
			//sizeSC: number of rows of the selfcorrelation matrix
			int sizeSC = NP-iD;
			SpEOMatrixF selfcorrelation(sizeSC,4);
			SpEOMatrixF selfcorrelationOrdered(sizeSC,4);

			// Loop Through all patches - find the net correlation of each patch to the current dictionary
			for(int iN = 0; iN<NP; iN++){
				int uN = iN / NPV;
				int vN = iN % NPV;
				sum = 0;
				// anyof: switch deciding whether the patch iN is any of the patches in the dictionary
				int anyof=0;
				// If patch is not in selected patches (sP)
				// anyof switch
				//cout<<"selectedPatches first"<<selectedPatches<<endl;
				for(int sP=0;sP<iD;sP++){
					//cout<<"iN: "<<iN<<endl;
					//cout<<"selectedPatches.rows() : "<<selectedPatches.rows()<<endl;
					//cout<<"sP: "<<sP<<endl;
					// cout<<"selectedPatches(sP): "<<selectedPatches(sP)<<endl;

					if(iN == selectedPatches(sP)){
						anyof=anyof+1;
					}
				}

				//cout<<"anyof: "<<anyof<<endl;
				//cout<<"--------- "<<endl;
				if(anyof==0){
					// This is the comparison patch
					compP = ImX->get_rasterBands()[iG]->get_bandDataMat()->block((*idxPUH)(uN), (*idxPVH)(vN), pszH, pszH);
					compV = SpEOVectorF::Map(compP.data(), pszH2);
					// Subtract mean and normalise
					compV.array() -= compV.mean();
					float norm_tmp = compV.norm();
					compV /= norm_tmp;

					for(int sP=0;sP<iD;sP++){

						// Find correlation between compare patch and patches that are ALREADY SELECTED IN THE DICTIONARY
						// Sum these correlations
						// iiN: index of patch inside the vector selectedPatches
						int iiN = selectedPatches(sP);
						int uiN = iiN / NPV;
						int viN = iiN % NPV;
						SpEOMatrixF alreadySelectedP = ImX->get_rasterBands()[iG]->get_bandDataMat()->block((*idxPUH)(uiN), (*idxPVH)(viN), pszH, pszH);
						SpEOVectorF alreadySelectedV = SpEOVectorF::Map(alreadySelectedP.data(), pszH2);
						alreadySelectedV.array() -= compV.mean();
						norm_tmp = alreadySelectedV.norm();
						alreadySelectedV /= norm_tmp;

						//find the absolute correlation between the patch in question and the current patch ALREADY SELECTED IN THE DICTIONARY
						float cc = abs(alreadySelectedV.dot(compV));
						sum = sum + cc;

					}
					//cout<<"sum: "<<sum<<endl;
					//cout<<"selectedPatches"<<selectedPatches<<endl;
					//cout<<"size selfcorrelation: "<<sizeSC<<endl;
					//cout<<"size selected patches: "<<selectedPatches.rows()<<endl;
					//					cout<<"iN: "<<iN<<endl;
					//					cout<<"iD: "<<iD<<endl;
					//					cout<<"cnt: "<<cnt<<endl;
					selfcorrelation.col(0).row(cnt) <<  uN;
					selfcorrelation.col(1).row(cnt) <<  vN;
					selfcorrelation.col(2).row(cnt) <<  sum;
					selfcorrelation.col(3).row(cnt) <<  iN;
					cnt = cnt+1;


				}//<--- END if(iN != selectedPatches.row(sP))
				//cout<<"--------- "<<endl;


			} //<-- for looping through iN patchs to form selfcorrelation matrix
			// Order selfcorrelation matrix

			VectorXi sort_order(sizeSC); // VectorXi sort_order(patchComp.rows());
			std::vector<argsort_pair> data(sizeSC); //std::vector<argsort_pair> data(patchComp.col(2).size());
			for(int i=0;i<sizeSC;i++) {//for(int i=0;i<patchComp.col(2).size();i++) {
				data[i].first = i;
				data[i].second = selfcorrelation.col(2)(i);//data[i].second = patchComp(2,i);//
			}

			// Sort data matrix according to the norms
			std::sort(data.begin(), data.end(), argsort_comp);
			// Get permutations ---> sort_order vector
			for(int i=0;i<sizeSC;i++) { //for(int i=0;i< data.size();i++) {
				sort_order(data[i].first) = i;
			}


			selfcorrelationOrdered = sort_order.asPermutation()*selfcorrelation;
			// the most uncorrelated patch is placed in patch comp ordered
			patchCompOrdered.col(0).row(iD) << selfcorrelationOrdered.col(0).row(0);
			patchCompOrdered.col(1).row(iD) << selfcorrelationOrdered.col(1).row(0);
			patchCompOrdered.col(2).row(iD) << selfcorrelationOrdered.col(2).row(0);
			// Update Selected patches list of indices

			int indexadd = selfcorrelationOrdered(0,3);
			//cout<<"indexadd"<<indexadd<<endl;
			selectedPatches.row(iD)<< indexadd;

		}
		*patchCompSubsetfloat = patchCompOrdered.block(0,0,NDP,3);

		break;
	}

	default :{
		sortON = false;
		if(my_rank == 0 && iP == 0){
			cout << "##### WARNING: Unknown dictionary selection method: " << fSet->dictselect << "! Choose a method between 0 and 8! The nearest neighbour dictionary selection method is used as the default!" << endl;
		}

		// Nearest Neighbor Default
		(*patchCompSubsetfloat)(0,0) = uP;
		(*patchCompSubsetfloat)(0,1) = vP;
		int pntCnt = 0;
		// (du, dv) is a vector - direction in which we move right now
		int du = 1;
		int dv = 0;
		int segment_length = 1;
		// current position (u, v) and how much of current segment we passed
		int u = uP;
		int v = vP;
		int segment_passed = 0;
		int k = 1;
		while (k < NDP) {
			// make a step, add 'direction' vector (du, dv) to current position (u, v)
			u += du;
			v += dv;
			++segment_passed;
			if (segment_passed == segment_length) {
				// done with current segment
				segment_passed = 0;
				// 'rotate' directions
				int buffer = du;
				du = -dv;
				dv = buffer;
				// increase segment length if necessary
				if (dv == 0) {
					++segment_length;
				}
			}
			if (u >= 0 && u <= glPrms->NPU-1 && v >= 0 && v <= glPrms->NPV-1) {
				(*patchCompSubsetfloat)(k,0) = u;
				(*patchCompSubsetfloat)(k,1) = v;
				pntCnt++;
				k++;
			}
		}
		break;
	}

	}// <----------------- END OF SWITCH STRUCTURE

	// ##### Sorting #####
	// Order the patchComp matrix of the form |NPIndices.coeff(0)|NPIndices.coeff(1)|Norm|
	// such that the metric column is arranged in increasing order
	if(sortON==true){

		// Dependencies:
		// typedef std::pair<int, float> argsort_pair ---> Defined in dataIO.h
		// argsort_comp ---> Defined in dataIO.cpp, dataIO.h

		// Note: patchcomp2 must be sorted for dictselect=4

		if(fSet->dictselect==4){
			VectorXi sort_order(NP); // VectorXi sort_order(patchComp.rows());
			std::vector<argsort_pair> data(NP); //std::vector<argsort_pair> data(patchComp.col(2).size());
			for(int i=0;i<NP;i++) {
				data[i].first = i;
				data[i].second = patchComp2.col(2)(i);
			}
			// Sort data matrix according to the norms
			std::sort(data.begin(), data.end(), argsort_comp);
			// Get permutations ---> sort_order vector
			for(int i=0;i<NP;i++) {
				sort_order(data[i].first) = i;
			}

			// Generate matrix of indices and norms, in ascending order according to the Norm column
			patchCompOrdered = sort_order.asPermutation()*patchComp2;

			// Take first NDP rows of Nonlocal patchCompOrdered table to form dictionary of size NDP
			*patchCompSubsetfloat = patchCompOrdered.block(0,0,NDP,3);
		}
		else{
			VectorXi sort_order(NP); // VectorXi sort_order(patchComp.rows());
			std::vector<argsort_pair> data(NP); //std::vector<argsort_pair> data(patchComp.col(2).size());
			for(int i=0;i<NP;i++) {//for(int i=0;i<patchComp.col(2).size();i++) {
				data[i].first = i;
				data[i].second = patchComp->col(2)(i);//data[i].second = patchComp(2,i);//
			}
			// Sort data matrix according to the norms
			std::sort(data.begin(), data.end(), argsort_comp);
			// Get permutations ---> sort_order vector
			for(int i=0;i<NP;i++) {
				sort_order(data[i].first) = i;
			}

			// Generate matrix of indices and norms, in ascending order according to the Norm column
			patchCompOrdered = sort_order.asPermutation()*(*patchComp);

			// Take first NDP rows of Nonlocal patchCompOrdered table to form dictionary of size NDP
			*patchCompSubsetfloat = patchCompOrdered.block(0,0,NDP,3);

		}


	}


	// In case the user wants to check the norms and ordering of patches are appropriate print the lines below:
	//		VectorXf normsOrdered = patchCompOrdered.col(2);
	// [
	//cout<<"patchCompSubset"<<endl<<patchCompSubsetfloat->block(0,0,15,3)<<endl;
	//	 cout<<"normsOrdered"<<endl<<normsOrdered.block(0,0,30,1)<<endl;//(NP-30,0,29,1)<<endl;
	//	 cout<<"---------"<<endl;
	// ]

	//return SpEOMatrixF(patchCompSubset);
	currP.resize(pszL,pszL);
	compP.resize(pszL,pszL);
	currP = SpEOMatrixF::Zero(pszL,pszL);
	compP = SpEOMatrixF::Zero(pszL,pszL);
	currenttmp = SpEOMatrixF::Zero(pszL,pszL);

}




























void dictSelectFunc(SpEOMatrixD *patchCompSubsetDouble, SpEODataset *ImX_LR, SpEODataset *ImX, SpEODataset *ImY, SpEOFusionSetting *fSet, SpEOGlobalParams *glPrms, int NP, int iP, int uP, int vP, SpEOVectorI *idxPUH, SpEOVectorI *idxPVH, SpEOVectorI *idxPUL, SpEOVectorI *idxPVL, SpEOMatrixD *SRF, SpEOMatrixD *patchComp, int my_rank, int ipp, int &NDP){	// Current Patch and Comparison Patch
	int pszH = fSet->patchsize*glPrms->fDS;
	int pszL = fSet->patchsize;
	int pszL2 = pszL*pszL;
	int pszH2 = pszH*pszH;
	int NPV = glPrms->NPV;
	bool sortON;

//	// FLOAT
//	SpEOMatrixF currenttmp = SpEOMatrixF::Zero(pszL,pszL);
//	SpEOMatrixF currP = SpEOMatrixF::Zero(pszL,pszL);
//	SpEOMatrixF compP = SpEOMatrixF::Zero(pszL,pszL);
//	// Current Patch and Comparison Patch in vector format
//	SpEOVectorF currV, compV;
//	// Patch indices and norms: SpEOMatrixF patchComp(NP,3);
//	// Ordered patch indices and norms
//	SpEOMatrixF patchComp2; // Necessary for dictselect=4 Do Not Delete!
//	SpEOMatrixF patchCompOrdered(NP,3);
//	SpEOMatrixF patchCompSubset(NDP,3);

    // DOUBLE
	SpEOMatrixD currenttmp = SpEOMatrixD::Zero(pszL,pszL);
	SpEOMatrixD currP = SpEOMatrixD::Zero(pszL,pszL);
	SpEOMatrixD compP = SpEOMatrixD::Zero(pszL,pszL);
	// Current Patch and Comparison Patch in vector format
	SpEOVectorD currV, compV;
	// Patch indices and norms: SpEOMatrixD patchComp(NP,3);
	// Ordered patch indices and norms
	SpEOMatrixD patchComp2; // Necessary for dictselect=4 Do Not Delete!
	SpEOMatrixD patchCompOrdered(NP,3);
	SpEOMatrixD patchCompSubset(NDP,3);


	int iG;
	if(fSet->useNewMethodForCalculatingZ){
		iG= ipp;
	}else{
		iG= glPrms->myChX[ipp];
	}

	switch(fSet->dictselect){

	case 0:{// Case (1): Nearest Neighbor patches according to the maximum norm of the patch coordinates
			sortON = false;
			(*patchCompSubsetDouble)(0,0) = uP;
			(*patchCompSubsetDouble)(0,1) = vP;
			break;
	}
	case 1:{// Case (1): Nearest Neighbor patches according to the maximum norm of the patch coordinates
		sortON = false;
		(*patchCompSubsetDouble)(0,0) = uP;
		(*patchCompSubsetDouble)(0,1) = vP;
		int pntCnt = 0;
		// (du, dv) is a vector - direction in which we move right now
		int du = 1;
		int dv = 0;
		int segment_length = 1;
		// current position (u, v) and how much of current segment we passed
		int u = uP;
		int v = vP;
		int segment_passed = 0;
		int k = 1;
		int tot_patch_cnt = 0;
		while (k < NDP && tot_patch_cnt<glPrms->NP_sub) {
			// make a step, add 'direction' vector (du, dv) to current position (u, v)
			u += du;
			v += dv;
			++segment_passed;
			if (segment_passed == segment_length) {
				// done with current segment
				segment_passed = 0;
				// 'rotate' directions
				int buffer = du;
				du = -dv;
				dv = buffer;
				// increase segment length if necessary
				if (dv == 0) {
					++segment_length;
				}
			}

			if (   u >= 0 && u <= glPrms->NPU-1
				&& v >= 0 && v <= glPrms->NPV-1) {
				tot_patch_cnt ++;
				bool patch_of_zeros_found = false;
				for(int ig=0; ig<ImX->get_NCh(); ig++){
//					cout << "jhfidsa 1 ig=" << ig << endl;
					SpEOMatrixD patchHR_tmp = (ImX->get_rasterBands()[ig]->bandDataMatD.block(idxPUH->coeff(u), idxPVH->coeff(v), pszH, pszH));
//					cout << "jhfidsa 2 ig=" << ig << endl;
					SpEOMatrixD patchLR_tmp = (ImX_LR->get_rasterBands()[ig]->bandDataMatD.block(idxPUL->coeff(u), idxPVL->coeff(v), pszL, pszL));
					if(patchHR_tmp.norm()<1e-8 || patchLR_tmp.norm()<1e-8){
						patch_of_zeros_found = true;
						cout << "In either ImX_sim or in ImX_sim_LR: patch (uP,vP)=("<<uP<<","<<vP<<") inf band/group ig="<<ig<<" contains only zeros and was therefore excluded from the dictionary!" << endl;
					}
				}
				if(!patch_of_zeros_found){
					(*patchCompSubsetDouble)(k,0) = u;
					(*patchCompSubsetDouble)(k,1) = v;
					pntCnt++;
					k++;
				}
			}
		}
		if(k < NDP-1){
			NDP = k+1;
		}

		break;

	}

	case 2: {// Case (2): Compute norms based on |ImX_LRcurrent-ImX_LR|
		// Set the currP as the current patch in the ImX_LR image, for norm comparison.
		sortON = true;
//		int bp_cnt=0;
		//cout << "bp" << ++bp_cnt << endl;

//		if(fSet->useNewMethodForCalculatingZ){
//		currP = ImX_LR->get_rasterBands()[iG]->get_bandDataMat()->block((*idxPUL)(uP), (*idxPVL)(vP), pszL, pszL);
		currP = ImX_LR->get_rasterBands()[iG]->bandDataMatD.block(idxPUL->coeff(uP), idxPVL->coeff(vP), pszL, pszL);



		// Arrange current patch representation, currP, into a vector, currV.
		//cout << "bp" << ++bp_cnt << endl;
		currV = SpEOVectorD::Map(currP.data(), pszL2);
		//cout << "bp" << ++bp_cnt << endl;
		// Subtract mean and normalise
		currV.array() -= currV.mean();
		//cout << "bp" << ++bp_cnt << endl;
		double norm_tmp = currV.norm();
		//cout << "bp" << ++bp_cnt << endl;
		currV /= norm_tmp;
		//cout << "bp" << ++bp_cnt << endl;
		// Metric comparison
		for(int iN = 0; iN<NP; iN++){
			// Case of LR norm comparison
			// Get Indices for all NP patches
			int uN = iN / NPV;
			int vN = iN % NPV;
			//cout << "bp" << ++bp_cnt << "iN=" << iN << "_1" << endl;
//			compP = ImX_LR->get_rasterBands()[iG]->get_bandDataMat()->block((*idxPUL)(uN), (*idxPVL)(vN), pszL, pszL);
			compP = ImX_LR->get_rasterBands()[iG]->bandDataMatD.block(idxPUL->coeff(uN), idxPVL->coeff(vN), pszL, pszL);
			//cout << "bp" << ++bp_cnt << "iN=" << iN << "_2" << endl;
			compV = SpEOVectorD::Map(compP.data(), pszL2);
			//cout << "bp" << ++bp_cnt << "iN=" << iN << "_3" << endl;
			// Subtract mean and normalise
			compV.array() -= compV.mean();
			double norm_tmp = compV.norm();
			compV /= norm_tmp;
			//cout << "bp" << ++bp_cnt << "iN=" << iN << "_4" << endl;
			(*patchComp)(iN,0) = uN;
			(*patchComp)(iN,1) = vN;
			(*patchComp)(iN,2) = (currV-compV).norm();
			//cout << "bp" << ++bp_cnt << "iN=" << iN << "_5" << endl;
		}

		break;
	}

	case 3: {// Case (3): Compute norms based on |ImX_approx-ImX|, ImX_approx = sum_Nc(ImY_channel*SRF_channel)
		// Use sum of ImY_channel*SRF_channel to form LR reconstruction as the currP
		sortON = true;
		for(int channel = 0; channel<fSet->Nc; channel++){
			SpEOMatrixF tttmp = (ImY->get_rasterBands()[channel]->get_bandDataMat()->block((*idxPUL)(uP), (*idxPVL)(vP), pszL, pszL) );
			currenttmp = (tttmp.cast<double>()*(*SRF)(0,channel));
			currP = currP + currenttmp;
		}
		// Arrange current patch representation, currP, into a vector.
		currV = SpEOVectorD::Map(currP.data(), pszL2);
		// Subtract mean and normalise
		currV.array() -= currV.mean();
		double norm_tmp = currV.norm();
		currV /= norm_tmp;


		// Metric comparison
		for(int iN = 0; iN<NP; iN++){
			// Case of LR norm comparison
			// Get Indices for all NP patches
			int uN = iN / NPV;
			int vN = iN % NPV;
			compP = ImX_LR->get_rasterBands()[iG]->bandDataMatD.block((*idxPUL)(uN), (*idxPVL)(vN), pszL, pszL);
			compV = SpEOVectorD::Map(compP.data(), pszL2);
			// Subtract mean and normalise
			compV.array() -= compV.mean();
			double norm_tmp = compV.norm();
			compV /= norm_tmp;
			(*patchComp)(iN,0) = uN;
			(*patchComp)(iN,1) = vN;
			(*patchComp)(iN,2) = (currV-compV).norm();
		}

		break;
	}

	case 4: {// Case(4): HR norm comparison
		// Set the currP as the current patch in the ImX_HR image, for norm comparison.
		sortON = true;
		currP.resize(pszH,pszH);
		compP.resize(pszH,pszH);
		currP = ImX->get_rasterBands()[iG]->bandDataMatD.block((*idxPUH)(uP), (*idxPVH)(vP), pszH, pszH);
		// Arrange current patch representation, currP, into a vector.
		currV = SpEOVectorD::Map(currP.data(), pszH2);
		// Subtract mean and normalise
		currV.array() -= currV.mean();
		double norm_tmp = currV.norm();
		currV /= norm_tmp;

		for(int iN = 0; iN<NP; iN++){
			int uN = iN / NPV;
			int vN = iN % NPV;
			// Case of HR Correlation comparison (less efficient than HR norm Comparison)
			// Get Indices for all NP patches
			compP = ImX->get_rasterBands()[iG]->bandDataMatD.block((*idxPUH)(uN), (*idxPVH)(vN), pszH, pszH);
			compV = SpEOVectorD::Map(compP.data(), pszH2);
			// Subtract mean and normalise
			compV.array() -= compV.mean();
			double norm_tmp = compV.norm();
			compV /= norm_tmp;
			(*patchComp)(iN,0) = uN;
			(*patchComp)(iN,1) = vN;
			(*patchComp)(iN,2) = (currV-compV).norm();
		}

		break;
	}
	case 5:{// Case(5): HR & LR Norm combined ranking !!! Inefficient !!!
		sortON = true;


		//################################   LR   ##################################
		// Set the currP as the current patch in the ImX_LR image, for norm comparison.
//		cout<<"iP??? "<<iP<<endl;
//		cout<<"uP??? "<<uP<<endl;
//		cout<<"vP??? "<<vP<<endl;
//		cout<<"patchComp.row(0) "<<patchComp->row(0)<<endl;

		currP = ImX_LR->get_rasterBands()[iG]->bandDataMatD.block((*idxPUL)(uP), (*idxPVL)(vP), pszL, pszL);
		// Arrange current patch representation, currP, into a vector.
		currV = SpEOVectorD::Map(currP.data(), pszL2);
		// Subtract mean and normalise
		currV.array() -= currV.mean();
		double norm_tmp = currV.norm();
		currV /= norm_tmp;
		// Case of LR norm comparison
		// Get Indices for all NP patches
		SpEOMatrixD LRComp=(*patchComp);
		SpEOMatrixD HRComp=(*patchComp);
		SpEOMatrixD LRCompOrdered(NP,3);
		for(int iN = 0; iN<NP; iN++){
			int uN = iN / NPV;
			int vN = iN % NPV;
			compP = ImX_LR->get_rasterBands()[iG]->bandDataMatD.block((*idxPUL)(uN), (*idxPVL)(vN), pszL, pszL);
			compV = SpEOVectorD::Map(compP.data(), pszL2);
			compV.array() -= compV.mean();
			double norm_tmp = compV.norm();
			compV /= norm_tmp;

			LRComp(iN,0) = uN;
			LRComp(iN,1) = vN;
			LRComp(iN,2) = (currV-compV).norm();
		}

		// Order indices according to norm metric
		VectorXi sort_orderLR(NP);
		std::vector<argsort_pair> dataLR(NP);
		for(int i=0;i<NP;i++) {
			dataLR[i].first = i;
			dataLR[i].second = LRComp.col(2)(i);
		}
		// Sort data matrix according to the norms
		std::sort(dataLR.begin(), dataLR.end(), argsort_comp);
		// Get permutations ---> sort_order vector
		for(int i=0;i< NP;i++) {
			sort_orderLR(dataLR[i].first) = i;
		}

		// Generate matrix of indices and norms, in ascending order according to the Norm column
		LRCompOrdered = sort_orderLR.asPermutation()*LRComp;
		// In case the user wants to check the norms and ordering of patches are appropriate print the lines below:
//				cout<<"LRCompOrdered"<<endl<<LRCompOrdered.block(0,0,5,3)<<endl;
//				cout<<"---------"<<endl;

		//################################   HR   ##################################
		// Set the currP as the current patch in the ImX_HR image, for norm comparison.
		currP.resize(pszH,pszH);
		compP.resize(pszH,pszH);
		currP = ImX->get_rasterBands()[iG]->bandDataMatD.block((*idxPUH)(uP), (*idxPVH)(vP), pszH, pszH);
		// Arrange current patch representation, currP, into a vector.
		currV = SpEOVectorD::Map(currP.data(), pszH2);
		// Subtract mean and normalise
		currV.array() -= currV.mean();
		norm_tmp = currV.norm();
		currV /= norm_tmp;
		// Case of HR norm comparison
		// Get Indices for all NP patches
		SpEOMatrixD HRCompOrdered(NP,3);
		for(int iN = 0; iN<NP; iN++){
			int uN = iN / NPV;
			int vN = iN % NPV;
			compP = ImX->get_rasterBands()[iG]->bandDataMatD.block((*idxPUH)(uN), (*idxPVH)(vN), pszH, pszH);
			compV = SpEOVectorD::Map(compP.data(), pszH2);
			compV.array() -= compV.mean();
			double norm_tmp = compV.norm();
			compV /= norm_tmp;
			HRComp(iN,2) = (currV-compV).norm();
		}

		// Order indices according to norm metric
		VectorXi sort_orderHR(NP);
		std::vector<argsort_pair> dataHR(NP);
		for(int i=0;i<NP;i++) {
			dataHR[i].first = i;
			dataHR[i].second = HRComp.col(2)(i);
		}
		// Sort data matrix according to the norms
		std::sort(dataHR.begin(), dataHR.end(), argsort_comp);
		// Get permutations ---> sort_order vector
		for(int i=0;i< NP;i++) {
			sort_orderHR(dataHR[i].first) = i;
		}

		// Generate matrix of indices and norms, in ascending order according to the Norm column
		HRCompOrdered = sort_orderHR.asPermutation()*HRComp;
		// In case the user wants to check the norms and ordering of patches are appropriate print the lines below:
//				cout<<"HRCompOrdered"<<endl<<HRCompOrdered.block(0,0,5,3)<<endl;
//				cout<<"---------"<<endl;


		// To avoid cross referenceing use patchComp2 (IMPORTANT - Don't delete)
		patchComp2 = LRCompOrdered;

		for(int i=0;i<NP;i++){
			// Compute combined ranking = sum of the ordered HR and LR row positions for each patch.
			for(int j=0;j<NP;j++){
				if(LRCompOrdered.row(i).col(0) == HRCompOrdered.row(j).col(0) && LRCompOrdered.row(i).col(1) == HRCompOrdered.row(j).col(1) ){
					// Combine HR an LR rankings
					patchComp2(i,2) =  i+j;
				}

			}
		}

		break;
	}
	case 6:{ // Case (6): Compute metric based on ABSOLUTE correlation between ImX_HRcurrent and ImX_HR
		// Set the currP as the current patch in the ImX_HR image, for norm comparison.
		sortON = true;
		currP.resize(pszH,pszH);
		compP.resize(pszH,pszH);
		currP = ImX->get_rasterBands()[iG]->bandDataMatD.block((*idxPUH)(uP), (*idxPVH)(vP), pszH, pszH);
		// Arrange current patch representation, currP, into a vector.
		currV = SpEOVectorD::Map(currP.data(), pszH2);
		// Subtract mean and normalise
		currV.array() -= currV.mean();
		double norm_tmp = currV.norm();
		currV /= norm_tmp;

		for(int iN = 0; iN<NP; iN++){
			int uN = iN / NPV;
			int vN = iN % NPV;
			// Case of HR Correlation comparison
			// Get Indices for all NP patches
			compP = ImX->get_rasterBands()[iG]->bandDataMatD.block((*idxPUH)(uN), (*idxPVH)(vN), pszH, pszH);
			compV = SpEOVectorD::Map(compP.data(), pszH2);
			// Subtract mean and normalise
			compV.array() -= compV.mean();
			double norm_tmp = compV.norm();
			compV /= norm_tmp;

			(*patchComp)(iN,0) = uN;
			(*patchComp)(iN,1) = vN;
			(*patchComp)(iN,2) = -abs(currV.dot(compV));
		}

		break;
	}
	case 7: {// Case (7): Compute metric based on **MINIMUM** ABSOLUTE correlation between ImX_HRcurrent and ImX_HR
		// Set the currP as the current patch in the ImX_HR image, for norm comparison.
		sortON = true;
		currP.resize(pszH,pszH);
		compP.resize(pszH,pszH);
		currP = ImX->get_rasterBands()[iG]->bandDataMatD.block((*idxPUH)(uP), (*idxPVH)(vP), pszH, pszH);
		// Arrange current patch representation, currP, into a vector.
		currV = SpEOVectorD::Map(currP.data(), pszH2);
		// Subtract mean and normalise
		currV.array() -= currV.mean();
		double norm_tmp = currV.norm();
		currV /= norm_tmp;

		for(int iN = 0; iN<NP; iN++){
			int uN = iN / NPV;
			int vN = iN % NPV;
			// Case of HR Correlation comparison
			// Get Indices for all NP patches
			compP = ImX->get_rasterBands()[iG]->bandDataMatD.block((*idxPUH)(uN), (*idxPVH)(vN), pszH, pszH);
			compV = SpEOVectorD::Map(compP.data(), pszH2);
			// Subtract mean and normalise
			compV.array() -= compV.mean();
			double norm_tmp = compV.norm();
			compV /= norm_tmp;

			(*patchComp)(iN,0) = uN;
			(*patchComp)(iN,1) = vN;
			// Insure current patch, iP, is the first element in the dictionary
			if(iN==iP){
				(*patchComp)(iN,2) = 0;
			}
			else{
				(*patchComp)(iN,2) = abs(currV.dot(compV));
			}
		}
		break;
	}
	case 8: {// Case (8): Compute metric based on random patch subset, though retaining the current patch as first dictionary atom
		// Set current patch as the first dictionary patch
		sortON = false;
		(*patchCompSubsetDouble)(0,0) = uP;
		(*patchCompSubsetDouble)(0,1) = vP;
		(*patchCompSubsetDouble)(0,2) = 0;

		// Loop through the remaining dictionary atoms and select random patches
		for(int iD = 1; iD<NDP; iD++){

			// Compute pseudo hashing function for seed (this could be improved by including processor rank).
			int time = MPI_Wtime();
			//int mod = 104729+iP;
			//int seed = abs(((time*181)*((iD-83)*iP*359))) % mod; // TO BE MODIFIED
			int seed = abs(time*121 - (my_rank+1)*12345 + (iP+1)*iD - uP +vP + (iG+1)*79) % 104729;
			// Seed random number
			srand(seed);
			// Compute random index
			int iR = rand() % NP;


			/* requires
			 * 1) to include -std=c++11 in the compiler options: e.g. CFLAGS   = -O3 -w -std=c++11
			 * 2) to uncomment #include <random> in the includes.h file
			int min = 0;
			int max = NP-1;
			random_device rd;     // only used once to initialise (seed) engine
			mt19937 rng(rd());    // random-number engine used (Mersenne-Twister in this case)
			uniform_int_distribution<int> uni(min,max); // guaranteed unbiased
			auto random_integer = uni(rng);
			int iR = random_integer;
			*/


			int uR = iR / NPV;
			int vR = iR % NPV;

			(*patchCompSubsetDouble)(iD,0) = uR;
			(*patchCompSubsetDouble)(iD,1) = vR;
			(*patchCompSubsetDouble)(iD,2) = iR;
		}

		break;
	}
	case 9: {// Case (9): Uncorrelated dictionary - i.e dictionary  attempts not to self correlate, to approximate a basis
		// WARNING: Case 9 involves an extremely expensive combinatorial search - NOT RECOMMENDED
		sortON = false;
		currP.resize(pszH,pszH);
		compP.resize(pszH,pszH);

		// First patch in the dictionary is the current patch
		patchCompOrdered.col(0).row(0) << uP;
		patchCompOrdered.col(1).row(0) << vP;
		patchCompOrdered.col(2).row(0) << 0;
		VectorXi selectedPatches(NDP);
		selectedPatches.row(0) << iP;

		// For the remainder atoms of the dictionary
		for(int iD = 1; iD<NDP; iD++){
			// sum: net summation of correlation coefficients between patch iN and current patches in the dictionary
			double sum = 0;
			// cnt: counter for the selfcorrelation matrix
			int cnt = 0;
			//sizeSC: number of rows of the selfcorrelation matrix
			int sizeSC = NP-iD;
			SpEOMatrixD selfcorrelation(sizeSC,4);
			SpEOMatrixD selfcorrelationOrdered(sizeSC,4);

			// Loop Through all patches - find the net correlation of each patch to the current dictionary
			for(int iN = 0; iN<NP; iN++){
				int uN = iN / NPV;
				int vN = iN % NPV;
				sum = 0;
				// anyof: switch deciding whether the patch iN is any of the patches in the dictionary
				int anyof=0;
				// If patch is not in selected patches (sP)
				// anyof switch
				//cout<<"selectedPatches first"<<selectedPatches<<endl;
				for(int sP=0;sP<iD;sP++){
					//cout<<"iN: "<<iN<<endl;
					//cout<<"selectedPatches.rows() : "<<selectedPatches.rows()<<endl;
					//cout<<"sP: "<<sP<<endl;
					// cout<<"selectedPatches(sP): "<<selectedPatches(sP)<<endl;

					if(iN == selectedPatches(sP)){
						anyof=anyof+1;
					}
				}

				//cout<<"anyof: "<<anyof<<endl;
				//cout<<"--------- "<<endl;
				if(anyof==0){
					// This is the comparison patch
					compP = ImX->get_rasterBands()[iG]->bandDataMatD.block((*idxPUH)(uN), (*idxPVH)(vN), pszH, pszH);
					compV = SpEOVectorD::Map(compP.data(), pszH2);
					// Subtract mean and normalise
					compV.array() -= compV.mean();
					double norm_tmp = compV.norm();
					compV /= norm_tmp;

					for(int sP=0;sP<iD;sP++){

						// Find correlation between compare patch and patches that are ALREADY SELECTED IN THE DICTIONARY
						// Sum these correlations
						// iiN: index of patch inside the vector selectedPatches
						int iiN = selectedPatches(sP);
						int uiN = iiN / NPV;
						int viN = iiN % NPV;
						SpEOMatrixD alreadySelectedP = ImX->get_rasterBands()[iG]->bandDataMatD.block((*idxPUH)(uiN), (*idxPVH)(viN), pszH, pszH);
						SpEOVectorD alreadySelectedV = SpEOVectorD::Map(alreadySelectedP.data(), pszH2);
						alreadySelectedV.array() -= compV.mean();
						norm_tmp = alreadySelectedV.norm();
						alreadySelectedV /= norm_tmp;

						//find the absolute correlation between the patch in question and the current patch ALREADY SELECTED IN THE DICTIONARY
						double cc = abs(alreadySelectedV.dot(compV));
						sum = sum + cc;

					}
					//cout<<"sum: "<<sum<<endl;
					//cout<<"selectedPatches"<<selectedPatches<<endl;
					//cout<<"size selfcorrelation: "<<sizeSC<<endl;
					//cout<<"size selected patches: "<<selectedPatches.rows()<<endl;
					//					cout<<"iN: "<<iN<<endl;
					//					cout<<"iD: "<<iD<<endl;
					//					cout<<"cnt: "<<cnt<<endl;
					selfcorrelation.col(0).row(cnt) <<  uN;
					selfcorrelation.col(1).row(cnt) <<  vN;
					selfcorrelation.col(2).row(cnt) <<  sum;
					selfcorrelation.col(3).row(cnt) <<  iN;
					cnt = cnt+1;


				}//<--- END if(iN != selectedPatches.row(sP))
				//cout<<"--------- "<<endl;


			} //<-- for looping through iN patchs to form selfcorrelation matrix
			// Order selfcorrelation matrix

			VectorXi sort_order(sizeSC); // VectorXi sort_order(patchComp.rows());
			std::vector<argsort_pair> data(sizeSC); //std::vector<argsort_pair> data(patchComp.col(2).size());
			for(int i=0;i<sizeSC;i++) {//for(int i=0;i<patchComp.col(2).size();i++) {
				data[i].first = i;
				data[i].second = selfcorrelation.col(2)(i);//data[i].second = patchComp(2,i);//
			}

			// Sort data matrix according to the norms
			std::sort(data.begin(), data.end(), argsort_comp);
			// Get permutations ---> sort_order vector
			for(int i=0;i<sizeSC;i++) { //for(int i=0;i< data.size();i++) {
				sort_order(data[i].first) = i;
			}


			selfcorrelationOrdered = sort_order.asPermutation()*selfcorrelation;
			// the most uncorrelated patch is placed in patch comp ordered
			patchCompOrdered.col(0).row(iD) << selfcorrelationOrdered.col(0).row(0);
			patchCompOrdered.col(1).row(iD) << selfcorrelationOrdered.col(1).row(0);
			patchCompOrdered.col(2).row(iD) << selfcorrelationOrdered.col(2).row(0);
			// Update Selected patches list of indices

			int indexadd = selfcorrelationOrdered(0,3);
			//cout<<"indexadd"<<indexadd<<endl;
			selectedPatches.row(iD)<< indexadd;

		}
		*patchCompSubsetDouble = patchCompOrdered.block(0,0,NDP,3);

		break;
	}

	default :{
		sortON = false;
		if(my_rank == 0 && iP == 0){
			cout << "##### WARNING: Unknown dictionary selection method: " << fSet->dictselect << "! Choose a method between 0 and 8! The nearest neighbour dictionary selection method is used as the default!" << endl;
		}

		// Nearest Neighbor Default
		(*patchCompSubsetDouble)(0,0) = uP;
		(*patchCompSubsetDouble)(0,1) = vP;
		int pntCnt = 0;
		// (du, dv) is a vector - direction in which we move right now
		int du = 1;
		int dv = 0;
		int segment_length = 1;
		// current position (u, v) and how much of current segment we passed
		int u = uP;
		int v = vP;
		int segment_passed = 0;
		int k = 1;
		while (k < NDP) {
			// make a step, add 'direction' vector (du, dv) to current position (u, v)
			u += du;
			v += dv;
			++segment_passed;
			if (segment_passed == segment_length) {
				// done with current segment
				segment_passed = 0;
				// 'rotate' directions
				int buffer = du;
				du = -dv;
				dv = buffer;
				// increase segment length if necessary
				if (dv == 0) {
					++segment_length;
				}
			}
			if (u >= 0 && u <= glPrms->NPU-1 && v >= 0 && v <= glPrms->NPV-1) {
				(*patchCompSubsetDouble)(k,0) = u;
				(*patchCompSubsetDouble)(k,1) = v;
				pntCnt++;
				k++;
			}
		}
		break;
	}

	}// <----------------- END OF SWITCH STRUCTURE

	// ##### Sorting #####
	// Order the patchComp matrix of the form |NPIndices.coeff(0)|NPIndices.coeff(1)|Norm|
	// such that the metric column is arranged in increasing order
	if(sortON==true){

		// Dependencies:
		// typedef std::pair<int, double> argsort_pair ---> Defined in dataIO.h
		// argsort_comp ---> Defined in dataIO.cpp, dataIO.h

		// Note: patchcomp2 must be sorted for dictselect=4

		if(fSet->dictselect==4){
			VectorXi sort_order(NP); // VectorXi sort_order(patchComp.rows());
			std::vector<argsort_pair> data(NP); //std::vector<argsort_pair> data(patchComp.col(2).size());
			for(int i=0;i<NP;i++) {
				data[i].first = i;
				data[i].second = patchComp2.col(2)(i);
			}
			// Sort data matrix according to the norms
			std::sort(data.begin(), data.end(), argsort_comp);
			// Get permutations ---> sort_order vector
			for(int i=0;i<NP;i++) {
				sort_order(data[i].first) = i;
			}

			// Generate matrix of indices and norms, in ascending order according to the Norm column
			patchCompOrdered = sort_order.asPermutation()*patchComp2;

			// Take first NDP rows of Nonlocal patchCompOrdered table to form dictionary of size NDP
			*patchCompSubsetDouble = patchCompOrdered.block(0,0,NDP,3);
		}
		else{
			VectorXi sort_order(NP); // VectorXi sort_order(patchComp.rows());
			std::vector<argsort_pair> data(NP); //std::vector<argsort_pair> data(patchComp.col(2).size());
			for(int i=0;i<NP;i++) {//for(int i=0;i<patchComp.col(2).size();i++) {
				data[i].first = i;
				data[i].second = patchComp->col(2)(i);//data[i].second = patchComp(2,i);//
			}
			// Sort data matrix according to the norms
			std::sort(data.begin(), data.end(), argsort_comp);
			// Get permutations ---> sort_order vector
			for(int i=0;i<NP;i++) {
				sort_order(data[i].first) = i;
			}

			// Generate matrix of indices and norms, in ascending order according to the Norm column
			patchCompOrdered = sort_order.asPermutation()*(*patchComp);

			// Take first NDP rows of Nonlocal patchCompOrdered table to form dictionary of size NDP
			*patchCompSubsetDouble = patchCompOrdered.block(0,0,NDP,3);

		}


	}


	// In case the user wants to check the norms and ordering of patches are appropriate print the lines below:
	//		VectorXf normsOrdered = patchCompOrdered.col(2);
	// [
	//cout<<"patchCompSubset"<<endl<<patchCompSubsetDouble->block(0,0,15,3)<<endl;
	//	 cout<<"normsOrdered"<<endl<<normsOrdered.block(0,0,30,1)<<endl;//(NP-30,0,29,1)<<endl;
	//	 cout<<"---------"<<endl;
	// ]

	//return SpEOMatrixD(patchCompSubset);
	currP.resize(pszL,pszL);
	compP.resize(pszL,pszL);
	currP = SpEOMatrixD::Zero(pszL,pszL);
	compP = SpEOMatrixD::Zero(pszL,pszL);
	currenttmp = SpEOMatrixD::Zero(pszL,pszL);

}





























// Use this function for easy sorting of tables using the std::sort() function
bool argsort_comp(const argsort_pair& left, const argsort_pair& right) {
    return left.second < right.second;
}

/////////////////////  END DICTIONARY SELECTION FUNCTIONS /////////////////////////////



void checkInputForJSpFI(SpEOGlobalParams &glPrms, SpEOFusionSetting &fSetting, SpEOParallelSetting &pSetting){
	int my_rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
//
	  if(pSetting.store_patches_tmp_on_drive){
		  if (my_rank==0){
			  cerr << endl << "ERROR: For the J-SparseFI regime, the option 'store_patches_tmp_on_drive' has to be set to 0!" << endl << endl;
		  }
		  MPI_Barrier(MPI_COMM_WORLD);
		  exit(2);
	  }
	  if (glPrms.NChY<glPrms.NChY_orig) {
		  if (my_rank==0){
			  cerr << endl << "ERROR: For the J-SparseFI regime with spectral groups there is currently no option to calculate only a spectral subset!" << endl << endl;
		  }
		  MPI_Barrier(MPI_COMM_WORLD);
		  exit(2);
	  }
	  if (glPrms.NChY!=8) {
		  if (my_rank==0){
			  cerr << endl << "ERROR: For the moment, the J-SparseFI regime is only possible for 8-band WorldView-2 data!" << endl << endl;
		  }
		  MPI_Barrier(MPI_COMM_WORLD);
		  exit(2);
	  }
	  if (glPrms.NChX!=1) {
		  if (my_rank==0){
			  cerr << endl << "ERROR: The J-SparseFI regime requires the high resolution image ImX to be a a single channel panchromatic image!" << endl << endl;
		  }
		  MPI_Barrier(MPI_COMM_WORLD);
		  exit(2);
	  }
	  if (glPrms.Nc_vec[0]!=glPrms.NChY || fSetting.Nc!=glPrms.NChY) {
		  if (my_rank==0){
			  cerr << endl << "ERROR: The J-SparseFI regime it is required to have glPrms.Nc_vec[0]==glPrms.NChY==fSetting.Nc! However, glPrms.NChY=" << glPrms.NChY << ", " << "fSetting.Nc=" << fSetting.Nc << ", and glPrms.Nc_vec[0]=" << glPrms.Nc_vec[0] << endl << endl;
		  }
		  MPI_Barrier(MPI_COMM_WORLD);
		  exit(2);
	  }
	  if (fSetting.No!=0) {
		  if (my_rank==0){
			  cerr << endl << "ERROR: The J-SparseFI regime does not allow spectral overlap! The argument 'No' has to be set to 0!" << endl << endl;
		  }
		  MPI_Barrier(MPI_COMM_WORLD);
		  exit(2);
	  }
	  if (glPrms.decMat_C.rows()!=1 || glPrms.decMat_C.cols()!=1) {
		  if (my_rank==0){
			  cerr << endl << "ERROR: The decision matrix C is supposed to contain only a single entry, which is equal to one, but is has the dimension " << glPrms.decMat_C.rows() << "x" << glPrms.decMat_C.cols() << endl << endl;
		  }
		  MPI_Barrier(MPI_COMM_WORLD);
		  exit(2);
	  }
	  if (glPrms.Ng!=1) {
		  if (my_rank==0){
			  cerr << endl << "ERROR: The J-SparseFI regime requires glPrms.Ng to be equal to 1!" << endl << endl;
		  }
		  MPI_Barrier(MPI_COMM_WORLD);
		  exit(2);
	  }
	  if(glPrms.myNumProbPerPatch != 1){
		  if (my_rank==0){
			  cerr << endl << "ERROR: glPrms->myNumProbPerPatch is supposed to be equal to 1! However, it is equal to " << glPrms.myNumProbPerPatch << endl << endl;
		  }
		  MPI_Barrier(MPI_COMM_WORLD);
		  exit(2);
	  }
	  if(glPrms.myBundle[0] != 0){
		  if (my_rank==0){
			  cerr << endl << "ERROR: For any processor, the quantity glPrms.myBundle[0] is supposed to be equal to 0! However, for process ["<< my_rank << "] myBundle = " << glPrms.myBundle[0] << endl << endl;
		  }
		  MPI_Barrier(MPI_COMM_WORLD);
		  exit(2);
	  }
	  if(glPrms.myChX[0] != 0){
		  if (my_rank==0){
			  cerr << endl << "ERROR: For any processor, the quantity glPrms.myChX[0] is supposed to be equal to 0! However, for process ["<< my_rank << "] myChX = " << glPrms.myChX[0] << endl << endl;
		  }
		  MPI_Barrier(MPI_COMM_WORLD);
		  exit(2);
	  }
	  if(pSetting.numProcPerPatch != 1){
		  if (my_rank==0){
			  cerr << endl << "ERROR: pSetting.numProcPerPatch has to be set to 1! However, numProcPerPatch==" << pSetting.numProcPerPatch << endl << endl;
		  }
		  MPI_Barrier(MPI_COMM_WORLD);
		  exit(2);
	  }

}


//template<typename Derived>
//inline bool is_finite(const Eigen::MatrixBase<Derived>& x)
bool is_inf_or_nan(double x){
   return !( x-x == x-x);
}
bool is_inf_or_nan(float x){
   return !( x-x == x-x);
}
bool is_inf_or_nan(int x){
   return !( x-x == x-x);
}
bool contains_inf_or_nan(SpEOMatrixD x){
   return !(( (x - x).array() == (x - x).array()).all());
}
bool contains_inf_or_nan(SpEOMatrixF x){
   return !(( (x - x).array() == (x - x).array()).all());
}
bool contains_inf_or_nan(SpEOMatrixI x){
   return !(( (x - x).array() == (x - x).array()).all());
}
bool contains_inf_or_nan(SpEOVectorD x){
   return !(( (x - x).array() == (x - x).array()).all());
}
bool contains_inf_or_nan(SpEOVectorF x){
   return !(( (x - x).array() == (x - x).array()).all());
}
bool contains_inf_or_nan(SpEOVectorI x){
   return !(( (x - x).array() == (x - x).array()).all());
}
void check_for_inf_or_nan(int my_rank,double x, const char *numDescr, int number, const char *descr){
	if(is_inf_or_nan(x)){
		cerr << endl
			 << "EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE" << endl
			 << "E   [my_rank=" << my_rank << "]:   ";
		if(number != -123){
			cerr << numDescr << number;
		}
		cerr << descr << "  is  NaN or INF!" << endl
			 << "EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE" << endl
			 << endl;
		//sleep(0.1);
		//exit(2);
	}
}
void check_for_inf_or_nan(int my_rank,float x, const char *numDescr, int number, const char *descr){
	if(is_inf_or_nan(x)){
		cerr << endl
			 << "EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE" << endl
			 << "E   [my_rank=" << my_rank << "]:   ";
		if(number != -123){
			cerr << numDescr << number;
		}
		cerr << descr << "  is  NaN or INF!" << endl
			 << "EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE" << endl
			 << endl;
		//sleep(0.1);
		//exit(2);
	}
}
void check_for_inf_or_nan(int my_rank,int x, const char *numDescr, int number, const char *descr){
	if(is_inf_or_nan(x)){
		cerr << endl
			 << "EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE" << endl
			 << "E   [my_rank=" << my_rank << "]:   ";
		if(number != -123){
			cerr << numDescr << number;
		}
		cerr << descr << "  is  NaN or INF!" << endl
			 << "EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE" << endl
			 << endl;
		//sleep(0.1);
		//exit(2);
	}
}
void check_for_inf_or_nan(int my_rank, SpEOMatrixD x, const char *numDescr, int number, const char *descr){
	if(contains_inf_or_nan(x)){
		cerr << endl
			 << "EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE" << endl
			 << "E   [my_rank=" << my_rank << "]:   ";
		if(number != -123){
			cerr << numDescr << number;
		}
		cerr << descr << "  contains  NaN or INF!" << endl
			 << "EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE" << endl
			 << endl;
		//sleep(0.1);
		//exit(2);
	}
}
void check_for_inf_or_nan(int my_rank, SpEOMatrixF x, const char *numDescr, int number, const char *descr){
	if(contains_inf_or_nan(x)){
		cerr << endl
			 << "EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE" << endl
			 << "E   [my_rank=" << my_rank << "]:   ";
		if(number != -123){
			cerr << numDescr << number;
		}
		cerr << descr << "  contains  NaN or INF!" << endl
			 << "EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE" << endl
			 << endl;
		//sleep(0.1);
		//exit(2);
	}
}
void check_for_inf_or_nan(int my_rank, SpEOMatrixI x, const char *numDescr, int number, const char *descr){
	if(contains_inf_or_nan(x)){
		cerr << endl
			 << "EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE" << endl
			 << "E   [my_rank=" << my_rank << "]:   ";
		if(number != -123){
			cerr << numDescr << number;
		}
		cerr << descr << "  contains  NaN or INF!" << endl
			 << "EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE" << endl
			 << endl;
		//sleep(0.1);
		//exit(2);
	}
}
void check_for_inf_or_nan(int my_rank, SpEOVectorD x, const char *numDescr, int number, const char *descr){
	if(contains_inf_or_nan(x)){
		cerr << endl
			 << "EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE" << endl
			 << "E   [my_rank=" << my_rank << "]:   ";
		if(number != -123){
			cerr << numDescr << number;
		}
		cerr << descr << "  contains  NaN or INF!" << endl
			 << "EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE" << endl
			 << endl;
		//sleep(0.1);
		//exit(2);
	}
}
void check_for_inf_or_nan(int my_rank, SpEOVectorF x, const char *numDescr, int number, const char *descr){
	if(contains_inf_or_nan(x)){
		cerr << endl
			 << "EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE" << endl
			 << "E   [my_rank=" << my_rank << "]:   ";
		if(number != -123){
			cerr << numDescr << number;
		}
		cerr << descr << "  contains  NaN or INF!" << endl
			 << "EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE" << endl
			 << endl;
		//sleep(0.1);
		//exit(2);
	}
}
void check_for_inf_or_nan(int my_rank, SpEOVectorI x, const char *numDescr, int number, const char *descr){
	if(contains_inf_or_nan(x)){
		cerr << endl
			 << "EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE" << endl
			 << "E   [my_rank=" << my_rank << "]:   ";
		if(number != -123){
			cerr << numDescr << number;
		}
		cerr << descr << "  contains  NaN or INF!" << endl
			 << "EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE" << endl
			 << endl;
		//sleep(0.1);
		//exit(2);
	}
}

//template<typename Derived2>
//inline bool is_nan(const Eigen::MatrixBase<Derived>& x)
bool contains_nan(SpEOMatrixD x){
 return !(( x.array() == x.array() ).all());
//   return ( (x - x).array() == (x - x).array()).all();
}
bool contains_nan(SpEOVectorD x){
   return !((x.array() == x.array()).all());
}

double calc_Corr_Coeff(SpEOMatrixD mat1, SpEOMatrixD mat2){
	int numEl = mat1.rows()*mat1.cols();
  	double mean1 = mat1.mean();
  	double mean2 = mat2.mean();
	SpEOMatrixD mat1_zrMn = mat1.array()-mean1;
	SpEOMatrixD mat2_zrMn = mat2.array()-mean2;
	mat1_zrMn.resize(numEl,1);
	mat2_zrMn.resize(numEl,1);
	double CC = (mat1_zrMn.transpose()*mat2_zrMn)(0,0)/(sqrt((mat1_zrMn.transpose()*mat1_zrMn)(0,0))*sqrt((mat2_zrMn.transpose()*mat2_zrMn)(0,0)));
	return CC;
}

void lowPassFilter_and_downSample(SpEODataset *Im,SpEODataset *Im_LR, SpEODataFormat input_format, SpEODataFormat output_format, SpEOGlobalParams &glPrms){
	int NCh = Im->get_NCh();
	SpEOMatrixD gauss_filter;
	int filter_size = 2*glPrms.fDS - glPrms.fDS%2;
	create_Gaussian_filter(gauss_filter, filter_size);
	SpEOMatrixD filter_coeff;
	calc_filter_boundary_coeff(filter_coeff, gauss_filter, glPrms.fDS);
	for(int iCh=0; iCh<NCh; iCh++){
		SpEOMatrixD mat_LR_tmp, mat_HR_tmp;
		mat_LR_tmp = SpEOMatrixD::Zero(glPrms.sizeUL,glPrms.sizeVL);
		if(input_format==SpEODouble){
			mat_HR_tmp = Im->get_rasterBands()[iCh]->bandDataMatD;
		}else{
			mat_HR_tmp = Im->get_rasterBands()[iCh]->get_bandDataMat()->cast<double>();
		}
		fast_filter(mat_LR_tmp, mat_HR_tmp, gauss_filter, filter_coeff, glPrms.fDS); // fast method (intuitive)
		if(output_format==SpEODouble){
			Im_LR->get_rasterBands()[iCh]->bandDataMatD = mat_LR_tmp;
		}else{
			Im_LR->get_rasterBands()[iCh]->bandDataMat = mat_LR_tmp.cast<float>();
		}
	}

}

void calc_ImX_sim(SpEODataset *ImX_sim, SpEODataset *ImX, SpEODataset *ImX_LR, SpEODataset *ImY,
		SpEOFusionSetting *fSetting, SpEODataIOSetting *dSetting, SpEOGlobalParams *glPrms,
		SpEOMatrixD *patX_LR, SpEOMatrixD *patY, SpEOVectorI *idxPUL, SpEOVectorI *idxPVL, int uP, int vP,
		bool localCalculation, int Ng, int *Nc_vec, int *idxChY, int sim_mode, MPI_Comm comm_busy){

	  int my_rank,iChX,iChY;
	  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
	  int chBundleFirst_tmp = dSetting->chBundleFirst; dSetting->chBundleFirst = 0;
	  int chBundleLast_tmp;
	  if(fSetting->useNewMethodForCalculatingZ){
		  chBundleLast_tmp  = dSetting->chBundleLast;  dSetting->chBundleLast  = Ng-1;
	  }else{
		  chBundleLast_tmp  = dSetting->chBundleLast;  dSetting->chBundleLast  = glPrms->numProbPerPatch-1;
	  }
	  int sizeUH_red_tmp = glPrms->sizeUH_red; glPrms->sizeUH_red = glPrms->sizeUH;
	  int sizeVH_red_tmp = glPrms->sizeVH_red; glPrms->sizeVH_red = glPrms->sizeVH;
	  setMetaInfo(ImX_sim, ImY, ImX, dSetting, glPrms);
//	  int sizeUH_red_tmp = glPrms->sizeUH_red; glPrms->sizeUH_red = glPrms->sizeUL_red;
//	  int sizeVH_red_tmp = glPrms->sizeVH_red; glPrms->sizeVH_red = glPrms->sizeVL_red;
//	  glPrms->sizeUH_red = glPrms->sizeUL;
//	  glPrms->sizeVH_red = glPrms->sizeVL;
	  //setMetaInfo(ImX_sim_LR, ImY, ImX_LR, &dSetting, &glPrms);
	  glPrms->sizeUH_red = sizeUH_red_tmp;
	  glPrms->sizeVH_red = sizeVH_red_tmp;
	  dSetting->chBundleFirst = chBundleFirst_tmp;
	  dSetting->chBundleLast  = chBundleLast_tmp;


	  //fSetting->useSimulatedImXforDictLearn = false;
	  //fSetting->useSimulatedImXforDictLearn = false; // TO BE OUTSOURCED
	  if(fSetting->useSimulatedImXforDictLearn){
//		  if(my_rank==0){
//			  cout << "   use simulated ImX for dictionary learning.." << endl;
//		  }
		  //>>>>>>>>>>>>>>> TO BE DONE >>>>>>>>>>>>>>>>>>

		  //int *Nc_vec = new int[Ng];

//		  SpEOVectorI idxChY = SpEOVectorI::Zero(Ng,1);
//		  int aY = fSetting->Nc-fSetting->No;
//		  int iG;
//		  for(iG=0; iG<Ng-1; iG++){
//			  idxChY(iG) = aY*iG;
//		  }
//		  idxChY(Ng-1) = glPrms->NChY-glPrms->Nc_vec[glPrms->numProbPerPatch-1];
		  //<<<<<<<<<<<<<<< TO BE DONE <<<<<<<<<<<<<<<<<<

//		  if(my_rank==0){cout << "glPrms->numProbPerPatch = " << glPrms->numProbPerPatch << endl;}
//		  if(my_rank==0){cout << "Ng = " << Ng << endl;}

//			if(my_rank==0){
//				for(int ig=0; ig<Ng; ig++){
//					cout << "ig="<<ig<<", ( idxChY[ig=" << ig << "], Nc_vec[ig=" << ig << "] ) = (" << idxChY[ig] << ", " << Nc_vec[ig] << " )" << endl;
////						cout << "Nc_vec[ig=" << ig << "]=" << Nc_vec[ig] << endl;
//				}
//			}
			//if(my_rank==0){ cout << "bp1" << endl << endl;}



		  // get window correct indeces
		  int winSizeL = fSetting->winSize;
                  //int winSizeH = winSizeL*glPrms->fDS;
                  int idxWUL = max(0, idxPUL->coeff(uP)-(int)(0.5*(double)(winSizeL-fSetting->patchsize)) );
                  int idxWVL = max(0, idxPVL->coeff(vP)-(int)(0.5*(double)(winSizeL-fSetting->patchsize)) );
                  int winSizeUL = min((int)ImY->get_sizeU()-1,  idxPUL->coeff(uP)-(int)(0.5*(double)(winSizeL-fSetting->patchsize))  + winSizeL-1   )-idxWUL+1;
                  int winSizeVL = min((int)ImY->get_sizeV()-1,  idxPVL->coeff(vP)-(int)(0.5*(double)(winSizeL-fSetting->patchsize))  + winSizeL-1   )-idxWVL+1;
                  //int idxWUH = glPrms->fDS*idxWUL;
                  //int idxWVH = glPrms->fDS*idxWVL;
                  //int winSizeUH = glPrms->fDS*winSizeUL;
                  //int winSizeVH = glPrms->fDS*winSizeVL;
		  //cout << "["<<my_rank<<"] uP=" << uP << " vP=" << vP <<" winSizeL=" << winSizeL << " idxWUL=" << idxWUL << " idxWVL=" << idxWVL << " winSizeUL=" << winSizeUL << " winSizeVL=" << winSizeUL << endl;
		  //MPI_Barrier(comm_busy);

		  for(int iG=0; iG<Ng; iG++){
//			  cout << "**************************" << endl
//				   << "*    iG=" << iG << " (of Ng=" << Ng << endl;
//		  for(int ipp=0; ipp<glPrms->numProbPerPatch; ipp++){
//			  int iG  = glPrms->myBundle[ipp];

//			  if(my_rank==0){cout << "ipp = " << ipp << ", iG = " << iG << endl;}
//			  if(my_rank==0){
//				  cout << "------------------------" << endl
//					   << "ipp=" << ipp << ", iG=" << iG << endl
//					   << "------------------------" << endl;
//			  }
			  SpEOMatrixD patX_LR_mod = *patX_LR;

			  // calculate average of all HS bands in current group
			  SpEOVectorD ImY_g_avg_vec;// = SpEOVectorD::Zero(fSetting->patchsize*fSetting->patchsize);
			  SpEOVectorD CCs = SpEOVectorD::Zero(glPrms->NChX); //double CCs[glPrms->NChX];

			  //cout << "bp1" << "iG="<<iG << endl;
			  double c_hat[glPrms->NChX]; // this is the vector of weights for different ImX bands estimated via one of the 3 following modes
			  int b_tilde[glPrms->NChX]; // ordering of ImX bands

			  switch(sim_mode){
				case 0: {
			  		if(localCalculation){

			  		        //if(my_rank==0){ cout << "bp2" << endl << endl;}
//			  		        for(iChX=0; iChX<patX_LR->rows(); iChX++){
//			  		      	  patX_LR->row(iChX).array()-=patX_LR->row(iChX).mean();
//			  		        }
//			  		        for(iChY=0; iChY<patY->rows(); iChY++){
//			  		      	  patY->row(iChY).array() -=patY->row(iChY).mean();
//			  		        }
					//cout << "[" << my_rank << "] " << "bp1" << endl;
					//MPI_Barrier(comm_busy);
			  		        ImY_g_avg_vec = SpEOVectorD::Zero(patY->cols());
			  		        int iC;
			  		        for(iC=0; iC<Nc_vec[iG]; iC++){
			  		      	  iChY = idxChY[iG]+iC;
			  		      	  ImY_g_avg_vec += patY->row(iChY);
			  		        }
			  		        ImY_g_avg_vec/=(double)(Nc_vec[iG]);

					//cout << "[" << my_rank << "] " << "bp2" << endl;
					//MPI_Barrier(comm_busy);
					
			  		        //if(my_rank==0){ cout << "bp3" << endl << endl;}
			  		        for(iChX=0; iChX<glPrms->NChX; iChX++){
			  		      	  CCs(iChX) = ImY_g_avg_vec.dot(patX_LR_mod.row(iChX));// CCs[iChX] = ImY_g_avg_vec.dot(patX_LR_mod.row(iChX));//calc_Corr_Coeff(ImY_g_avg, ImX_LR->get_rasterBands()[iChX]->get_bandDataMat()->cast<double>());
//			  		      	  double CC = (mat1_zrMn.transpose()*mat2_zrMn)(0,0)/(sqrt((mat1_zrMn.transpose()*mat1_zrMn)(0,0))*sqrt((mat2_zrMn.transpose()*mat2_zrMn)(0,0)));
			  		        }
			  		        // if all CCs are close to zero, then the corresponding band in ImX_sim will contain only
			  		        // entries that are close to zero which in turn results in dictionary atoms which have norm
			  		        // close to zero. This would cause NaN numbers as the dictionary atoms will be (and have to be)
			  		        //normalized. Hence, we enlarge the window size for this group until this problem is resolved
					//cout << "[" << my_rank << "] " << "bp3" << endl;
					//MPI_Barrier(comm_busy);
			  		        if (CCs.maxCoeff()<(1e-8) && CCs.minCoeff()>-(1e-8)){
			  		      	  winSizeL = fSetting->winSize;

			  		      	  int NChY = patY->rows();
			  		      	  int NChX = patX_LR_mod.rows();
			  		      	  bool extraRoundCompleted = false;
					//cout << "[" << my_rank << "] " << "bp4" << endl;
					//MPI_Barrier(comm_busy);
			  		      	  cout << "         ["<<my_rank<<"] WWWW Warning:"
			  		      	  								   << " iG="<<iG << " [ImY bands "<< idxChY[iG] << " to " << idxChY[iG]+Nc_vec[iG]-1 << "] -> all CCs (and probably the entire ImY patch in this group) are close to zero!"// << endl
//			  		      	  								   << " CCs = " << CCs.transpose() << endl
//			  		      	  								   << " ImY_g_avg_vec = " << ImY_g_avg_vec.transpose() << endl
			  		      	  								   << " => Enlarge winSizeL from " << winSizeL;
			  		      	  do{
			  		      		  if ( !(CCs.maxCoeff()<(1e-8) && CCs.minCoeff()>-(1e-8))){
			  		      			  extraRoundCompleted = true;
			  		      		  }
			  		      		  cout << " to " << winSizeL+4;

			  		      		  winSizeL += 4;
			  		      		  //int winSizeH = winSizeL*glPrms->fDS;
			  		      		  idxWUL = max(0, idxPUL->coeff(uP)-(int)(0.5*(double)(winSizeL-fSetting->patchsize)) );
			  		      		  idxWVL = max(0, idxPVL->coeff(vP)-(int)(0.5*(double)(winSizeL-fSetting->patchsize)) );
			  		      		  winSizeUL = min((int)ImY->get_sizeU()-1,  idxPUL->coeff(uP)-(int)(0.5*(double)(winSizeL-fSetting->patchsize))  + winSizeL-1   )-idxWUL+1;
			  		      		  winSizeVL = min((int)ImY->get_sizeV()-1,  idxPVL->coeff(vP)-(int)(0.5*(double)(winSizeL-fSetting->patchsize))  + winSizeL-1   )-idxWVL+1;
			  		      		  int idxWUH = glPrms->fDS*idxWUL;
			  		      		  int idxWVH = glPrms->fDS*idxWVL;
			  		      		  int winSizeUH = glPrms->fDS*winSizeUL;
			  		      		  int winSizeVH = glPrms->fDS*winSizeVL;
			  		      		  // extract window winY from ImY in each channel individually
			  		      		  SpEOMatrixD winY = SpEOMatrixD::Zero(NChY,winSizeUL*winSizeVL);
			  		      		  for(iChY=0; iChY<NChY; iChY++){
			  		      			  SpEOMatrixD windowBand = (ImY->get_rasterBands()[iChY]->get_bandDataMat()->block(idxWUL, idxWVL, winSizeUL, winSizeVL)).cast<double>();
			  		      			  SpEOVectorD winY_row = SpEOVectorD::Map(windowBand.data(), winSizeUL*winSizeVL);
			  		      			  double winY_row_mean = winY_row.mean();
			  		      			  winY_row.array() -= winY_row_mean;
			  		      			  double winY_row_norm = winY_row.norm();
			  		      			  if(winY_row_norm>1e-8){
			  		      				  winY_row /= winY_row_norm;
			  		      			  }
			  		      			  winY.row(iChY) = winY_row;
			  		      		  }
			  		      		  check_for_inf_or_nan(my_rank,winY, " ", -123, "winY");
//			  		      		  SpEOMatrixD winX_LR = SpEOMatrixD::Zero(NChX,winSizeUL*winSizeVL);
//			  		      		  for(iChX=0; iChX<NChX; iChX++){
//			  		      			  SpEOMatrixD windowBand = (ImX_LR->get_rasterBands()[iChX]->get_bandDataMat()->block(idxWUL, idxWVL, winSizeUL, winSizeVL)).cast<double>();
//			  		      			  winX_LR.row(iChX) = SpEOVectorD::Map(windowBand.data(), winSizeUL*winSizeVL);
//			  		      			  winX_LR.row(iChX).array() -= winX_LR.row(iChX).mean();
//			  		      			  double winX_LR_row_norm = winX_LR.row(iChX).norm();
//			  		      			  if(winX_LR_row_norm>1e-8){
//			  		      				  winX_LR.row(iChX) /= winX_LR_row_norm;
//			  		      			  }
//			  		      		  }
			  		      		  patX_LR_mod = SpEOMatrixD::Zero(NChX,winSizeUL*winSizeVL);
			  		      		  for(iChX=0; iChX<NChX; iChX++){
			  		      			  SpEOMatrixD windowBand = (ImX_LR->get_rasterBands()[iChX]->get_bandDataMat()->block(idxWUL, idxWVL, winSizeUL, winSizeVL)).cast<double>();
			  		      			  patX_LR_mod.row(iChX) = SpEOVectorD::Map(windowBand.data(), winSizeUL*winSizeVL);
			  		      			  patX_LR_mod.row(iChX).array() -= patX_LR_mod.row(iChX).mean();
			  		      			  double patX_LR_mod_row_norm = patX_LR_mod.row(iChX).norm();
			  		      			  if(patX_LR_mod_row_norm>1e-8){
			  		      				  patX_LR_mod.row(iChX) /= patX_LR_mod_row_norm;
			  		      			  }
			  		      		  }

			  		      		  check_for_inf_or_nan(my_rank,patX_LR_mod, " ", -123, "patX_LR_mod");

//			  		      		  for(iChX=0; iChX<NChX; iChX++){
//			  		      			  patX_LR_mod.row(iChX).array()-=patX_LR_mod.row(iChX).mean();
//			  		      		  }
		//	  		      		  for(iChY=0; iChY<winY->rows(); iChY++){
		//	  		      			  winY->row(iChY).array() -=winY->row(iChY).mean();
		//	  		      		  }
			  		      		  ImY_g_avg_vec = SpEOVectorD::Zero(winY.cols());
			  		      		  for(iC=0; iC<Nc_vec[iG]; iC++){
			  		      			  iChY = idxChY[iG]+iC;
			  		      			  ImY_g_avg_vec += winY.row(iChY);
//			  		      			  cout << "iC="<<iC<<", iChY="<<iChY<<", winY.row(iChY).mean() = " << winY.row(iChY).mean() << endl;
			  		      		  }
			  		      		  ImY_g_avg_vec/=(double)(Nc_vec[iG]);

			  		      		  //if(my_rank==0){ cout << "bp3" << endl << endl;}
					//cout << "[" << my_rank << "] " << "bp5" << endl;
					//MPI_Barrier(comm_busy);
			  		      		  for(iChX=0; iChX<glPrms->NChX; iChX++){
			  		      			  CCs(iChX) = ImY_g_avg_vec.dot(patX_LR_mod.row(iChX));// CCs[iChX] = ImY_g_avg_vec.dot(patX_LR_mod->row(iChX));//calc_Corr_Coeff(ImY_g_avg, ImX_LR->get_rasterBands()[iChX]->get_bandDataMat()->cast<double>());
		//	  		      			  double CC = (mat1_zrMn.transpose()*mat2_zrMn)(0,0)/(sqrt((mat1_zrMn.transpose()*mat1_zrMn)(0,0))*sqrt((mat2_zrMn.transpose()*mat2_zrMn)(0,0)));
			  		      		  }
//			  		      		  cout << " CCs = " << CCs.transpose() << endl
//			  		      			   << " ImY_g_avg_vec = " << ImY_g_avg_vec.transpose() << endl
//			  		      			   << " ImY_g_avg_vec.mean() = " << ImY_g_avg_vec.mean() << endl
//			  		      			   << " patX_LR_mod = " << endl
//			  		      			   <<   patX_LR_mod << endl;

			  		      		  if(winSizeUL>ImY->get_sizeU() && winSizeVL>ImY->get_sizeV()){
			  		      			  cerr << "         EEEE Error: Error occurred while enlarging the window for generating ImX_sim.." << endl
			  		      				   << "                     window exceeds the size of the image, which can only mean that " << endl
			  		      				   << "                     ImY contains only zeros in the current group (iG=" << iG << ") of bands.." << endl
			  		      				   << "                     In fact, ImY_g_avg_vec.norm() = " << ImY_g_avg_vec.norm() << endl
			  		      				   << endl;
			  		      			  break;
			  		      		  }
			  		      	  }while ( (CCs.maxCoeff()<(1e-8) && CCs.minCoeff()>-(1e-8)) || !extraRoundCompleted );
			  		      	  cout << endl;// << "<=== WWWWWW Warning" << endl;
			  		        }
			  		}else{
			  		        //cout << "bp6" << "iG="<<iG << endl;
					//MPI_Barrier(comm_busy);
			  		        //ImY_g_avg_vec = SpEOVectorD::Zero(glPrms->sizeUL*glPrms->sizeUL);
			  		        SpEOMatrixD ImY_g_avg = SpEOMatrixD::Zero(glPrms->sizeUL,glPrms->sizeVL);
			  		        int iC, iChY;
			  		        for(iC=0; iC<Nc_vec[iG]; iC++){
			  		      	  iChY = idxChY[iG]+iC;
			  		      	  SpEOMatrixD ttmp = ImY->get_rasterBands()[iChY]->get_bandDataMat()->cast<double>();
			  		      	  ttmp.array() -= ttmp.mean();
			  		      	  ttmp.array() /= ttmp.norm();
			  		      	  ImY_g_avg += ttmp;
			  		        }
			  		        //cout << "bp7" << "iG="<<iG << endl;
					//MPI_Barrier(comm_busy);
			  		        ImY_g_avg/=(double)(Nc_vec[iG]);
			  		        // calculate correlation between ImY_g_avg and all bands in ImX_LR
			  		        //int iChX;
			  		        for(iChX=0; iChX<glPrms->NChX; iChX++){
			  		      	  SpEOMatrixD ttmp = ImX_LR->get_rasterBands()[iChX]->get_bandDataMat()->cast<double>();
			  		      	  ttmp.array() -= ttmp.mean();
			  		      	  ttmp.array() /= ttmp.norm();
			  		      	  CCs(iChX) = calc_Corr_Coeff(ImY_g_avg, ttmp); //CCs[iChX] = calc_Corr_Coeff(ImY_g_avg, ttmp);
			  		        }
			  		        ImY_g_avg_vec = SpEOVectorD::Map(ImY_g_avg.data(),ImY_g_avg.rows()*ImY_g_avg.cols());
			  		        //cout << "bp8" << "iG="<<iG << endl;
					//MPI_Barrier(comm_busy);
			  		}

			  		//if(my_rank==0){ //cout << "bp9" << endl << endl;}
					//MPI_Barrier(comm_busy);
			  		// sort the CCs in a decreasing manner to find b_tilde[0],...,b_tilde[NChX-1]
			  		//int b_tilde[glPrms->NChX];
			  		std::vector<argsort_pair> CCs_for_sorting(glPrms->NChX);
			  		for(int iChX=0; iChX<glPrms->NChX; iChX++){
			  		        CCs_for_sorting[iChX].first = iChX;
			  		        CCs_for_sorting[iChX].second = CCs(iChX); //CCs[iChX];
			  		}
			  		//if(my_rank==0){ //cout << "bp10" << endl << endl;}
					//MPI_Barrier(comm_busy);
			  		std::sort(CCs_for_sorting.begin(), CCs_for_sorting.end(), argsort_comp);
			  		for(int iChX=0; iChX< glPrms->NChX; iChX++) {
			  		        b_tilde[CCs_for_sorting[glPrms->NChX-1-iChX].first] = iChX;
			  		}
			  		//cout << "[" << my_rank << "] " << "bp11" << "iG="<<iG << endl;
					//MPI_Barrier(comm_busy);
			  		//if(my_rank==0){ cout << "bp6" << endl << endl;}
//			  		for(int iChX=0; iChX< glPrms->NChX; iChX++) {
//			  		        cout << "CCs[iChX=" << iChX << "]=" << CCs[iChX] << ", CCs[b_tilde[" << iChX << "]=" << b_tilde[iChX] << "]=" << CCs[b_tilde[iChX]] << endl;
//			  		}
			  		// find a linear combination of the bands in ImX_LR to approximate ImY_g_avg while favoring the highly correlated bands
			  		//double c_hat[glPrms->NChX];

			  		SpEOVectorD vec2;
			  		//if(my_rank==0){ cout << "bp12" << endl << endl;}
					//MPI_Barrier(comm_busy);
			  		if(localCalculation){
			  		        vec2 = patX_LR_mod.row(b_tilde[0]);
			  		        //if(my_rank==0){ cout << "bp8" << endl << endl;}
			  		}else{
			  		        //cout << "bp5" << "iG="<<iG << endl;
//			  		        vec1 = SpEOVectorD::Map(ImY_g_avg.data(),ImY_g_avg.rows()*ImY_g_avg.cols());
			  		        SpEOMatrixD mat2 = ImX_LR->get_rasterBands()[b_tilde[0]]->get_bandDataMat()->cast<double>();
			  		        vec2 = SpEOVectorD::Map(mat2.data(),mat2.rows()*mat2.cols());
			  		        vec2.array() -= vec2.mean();
			  		        vec2.array() /= vec2.norm();
			  		}

			  		//if(my_rank==0){ cout << "bp13" << endl << endl;}
					//MPI_Barrier(comm_busy);

			  		c_hat[0] = 1/(vec2.dot(vec2))*vec2.dot(ImY_g_avg_vec); // least squares solution
			  		SpEOVectorD sum_ci_ImXLR_vec = c_hat[0]*vec2;
			  		for(iChX=1; iChX<glPrms->NChX; iChX++){

			  		        //>>>>>>>>>>>>>>> TO BE DONE: make it local >>>>>>>>>>>>>>>>>>
			  		        SpEOVectorD vec3;
			  		        if(localCalculation){
			  		      	  //if(my_rank==0){ cout << "bp14.iChX="<< iChX << endl << endl;}
					//MPI_Barrier(comm_busy);
			  		      	  vec3 = patX_LR_mod.row(b_tilde[iChX]);
			  		        }else{
			  		      	  SpEOMatrixD mat3 = ImX_LR->get_rasterBands()[b_tilde[iChX]]->get_bandDataMat()->cast<double>();
			  		      	  vec3 = SpEOVectorD::Map(mat3.data(),mat3.rows()*mat3.cols());
			  		      	  vec3.array() -= vec3.mean();
			  		      	  vec3.array() /= vec3.norm();
			  		        }
			  		        //<<<<<<<<<<<<<<< TO BE DONE: make it local <<<<<<<<<<<<<<<<<<


			  		        c_hat[iChX] = 1/(vec3.dot(vec3))*vec3.dot(ImY_g_avg_vec-sum_ci_ImXLR_vec);  // least squares solution
			  		        sum_ci_ImXLR_vec += c_hat[iChX]*vec3;
			  		}
					break;
			  	}
			  	case 1: { // unconstrained least squares
					//cout << "[" << my_rank << "] " << "bp1" << endl;
					//MPI_Barrier(comm_busy);
					// extract window winY from ImY in each channel of the current group 
                                  	//SpEOMatrixD winY_g = SpEOMatrixD::Zero(iG,winSizeUL*winSizeVL);
                                  	SpEOMatrixD winY_g = SpEOMatrixD::Zero(Nc_vec[iG],winSizeUL*winSizeVL);
					int iC;
                                  	for(iC=0; iC<Nc_vec[iG]; iC++){	
						iChY = idxChY[iG]+iC;
			  		      	SpEOMatrixD windowBand = (ImY->get_rasterBands()[iChY]->get_bandDataMat()->block(idxWUL, idxWVL, winSizeUL, winSizeVL)).cast<double>();
			  		      	winY_g.row(iC) = SpEOVectorD::Map(windowBand.data(), winSizeUL*winSizeVL);
                                  	}
					//cout << "[" << my_rank << "] " << "bp2" << endl;
					//MPI_Barrier(comm_busy);
					// calculate average ImY patch for current group
					SpEOVectorD winY_g_avg_vec = SpEOVectorD::Zero(winY_g.cols());
					for(iC=0; iC<Nc_vec[iG]; iC++){
						winY_g_avg_vec += winY_g.row(iC);
					}
					//cout << "[" << my_rank << "] " << "bp3" << endl;
					//MPI_Barrier(comm_busy);
					ImY_g_avg_vec/=(double)(Nc_vec[iG]);
					// extract window winX_LR from ImX_LR
                                   	SpEOMatrixD winX_LR = SpEOMatrixD::Zero(glPrms->NChX,winSizeUL*winSizeVL);
					//cout << "[" << my_rank << "] " << "bp4" << endl;
					//MPI_Barrier(comm_busy);
                                   	for(iChX=0; iChX<glPrms->NChX; iChX++){
                                   	        SpEOMatrixD windowBand = (ImX_LR->get_rasterBands()[iChX]->get_bandDataMat()->block(idxWUL, idxWVL, winSizeUL, winSizeVL)).cast<double>();
                                   	        winX_LR.row(iChX) = SpEOVectorD::Map(windowBand.data(), winSizeUL*winSizeVL);
                                   	}
					//cout << "[" << my_rank << "] " << "bp5" << endl;
					//MPI_Barrier(comm_busy);
					// calculate weights via unconstrained least squares
					SpEOMatrixD A = winX_LR.transpose();
					SpEOVectorD b = winY_g_avg_vec;
					//cout << "[" << my_rank << "] " << "bp6" << endl;
					//MPI_Barrier(comm_busy);
					SpEOVectorD x = A.jacobiSvd(ComputeThinU | ComputeThinV).solve(b);
					//cout << "[" << my_rank << "] " << "bp7" << endl;
					//MPI_Barrier(comm_busy);
					for(iChX=0; iChX<glPrms->NChX; iChX++){
						c_hat[iChX] = x(iChX);
						b_tilde[iChX]=iChX;
					}
					//if(my_rank==0) cout << "bp8" << endl;
					//MPI_Barrier(comm_busy);
					break;
			  	}
			  	case 2: { // nonnegative least squares
					// extract window winY from ImY in each channel of the current group 
                                  	SpEOMatrixD winY_g = SpEOMatrixD::Zero(Nc_vec[iG],winSizeUL*winSizeVL);
					//cout << "[" << my_rank << "] " << "bp1 iG=" << iG << " Nc_vec[iG]=" << Nc_vec[iG] << endl;
					//MPI_Barrier(comm_busy);
					int iC;
                                  	for(iC=0; iC<Nc_vec[iG]; iC++){	
						//cout << "[" << my_rank << "] iC=" << iC << " idxChY[iG]=" << idxChY[iG] << " " << endl;
						iChY = idxChY[iG]+iC;
						//cout << "[" << my_rank << "] bp1.1=" << endl;
			  		      	SpEOMatrixD windowBand = (ImY->get_rasterBands()[iChY]->get_bandDataMat()->block(idxWUL, idxWVL, winSizeUL, winSizeVL)).cast<double>();
						//cout << "[" << my_rank << "] bp1.2=" << endl;
			  		      	winY_g.row(iC) = SpEOVectorD::Map(windowBand.data(), winSizeUL*winSizeVL);
						//cout << "[" << my_rank << "] bp1.3=" << endl;
                                  	}
					// calculate average ImY patch for current group
					SpEOVectorD winY_g_avg_vec = SpEOVectorD::Zero(winY_g.cols());
					//cout << "[" << my_rank << "] " << "bp2" << endl;
					//MPI_Barrier(comm_busy);
					for(iC=0; iC<Nc_vec[iG]; iC++){
						winY_g_avg_vec += winY_g.row(iC);
					}
					ImY_g_avg_vec/=(double)(Nc_vec[iG]);
					//cout << "[" << my_rank << "] " << "bp3" << endl;
					//MPI_Barrier(comm_busy);
					// extract window winX_LR from ImX_LR
                                   	SpEOMatrixD winX_LR = SpEOMatrixD::Zero(glPrms->NChX,winSizeUL*winSizeVL);
                                   	for(iChX=0; iChX<glPrms->NChX; iChX++){
                                   	        SpEOMatrixD windowBand = (ImX_LR->get_rasterBands()[iChX]->get_bandDataMat()->block(idxWUL, idxWVL, winSizeUL, winSizeVL)).cast<double>();
                                   	        winX_LR.row(iChX) = SpEOVectorD::Map(windowBand.data(), winSizeUL*winSizeVL);
                                   	}
					// calculate weights via unconstrained least squares
					//cout << "[" << my_rank << "] " << "bp4" << endl;
					//MPI_Barrier(comm_busy);
					SpEOMatrixD A = winX_LR.transpose();
					SpEOVectorD b = winY_g_avg_vec;
					SpEOVectorD x = A.jacobiSvd(ComputeThinU | ComputeThinV).solve(b);

					//cout << "[" << my_rank << "] " << "bp5" << endl;
					//MPI_Barrier(comm_busy);
      					//if(my_rank == 0) cout << " iG=" << iG << ": ";
					//MPI_Barrier(comm_busy);
      					bool all_good = NNLS<MatrixXd>::solve(A, b, x);
					if(!all_good){
						cout << "[" << my_rank << "] WARNING: NNLS failed at uP=" << uP << " vP=" << vP << " iG=" << iG << endl 
						     << "The problem failed to solve was Ax=b where x>=0 and " << "A = " << A << endl << endl << "b'=" << endl << b.transpose() << endl; 
					}
      					//if(my_rank == 0) cout << ".. done! all_good=" << all_good << endl;

					for(iChX=0; iChX<glPrms->NChX; iChX++){
						c_hat[iChX] = x(iChX);
						b_tilde[iChX]=iChX;
					}
					//cout << "[" << my_rank << "] " << "bp6" << endl;
					//MPI_Barrier(comm_busy);
					break;
			  	}
				default: {
					break;
				}
			  }
 

//			  for(iChX=0; iChX<glPrms->NChX; iChX++){
//				  cout << "c_hat[iChX="<< iChX <<"]="<<c_hat[iChX]<< endl;
//			  }
			  // simulate HR image band ImX_sim so that ImX_sim*B*S is (locally at current patch) highly correlated with ImY_g_avg.
			  SpEOMatrixD ImX_sim_band = SpEOMatrixD::Zero(ImX->get_rasterBands()[0]->get_bandDataMat()->rows(),ImX->get_rasterBands()[0]->get_bandDataMat()->cols());
			  for(int iChX=0; iChX<glPrms->NChX; iChX++){
//				  if(iG==56){
//					  cout << "c_hat[iChX="<<iChX<<"] = " << c_hat[iChX] << endl;
//				  }
				  ImX_sim_band += c_hat[iChX]*ImX->get_rasterBands()[b_tilde[iChX]]->get_bandDataMat()->cast<double>();
				  //if(my_rank==0){
//					  if(contains_nan(ImX_sim_band)){
//						  cout << endl;
//						  cout << "nan nan nan nan nan nan nan nan nan nan nan nan " << endl
//							   << "nan"
//							   << "nan    found in     ImX_sim_band    corresp. to group iG=" << iG << endl
//							   << "nan"
//							   << "nan nan nan nan nan nan nan nan nan nan nan nan " << endl << endl;
//					  }
				  //}
			  }
			  //cout << "bp7" << "iG="<<iG << endl;
//			  if(my_rank==0 && iG==56){
//				  cout << "WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW" << endl
//					   << "WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW" << endl
//					   << "WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW" << endl
//					   << "WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW" << endl
//					   << "WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW" << endl;
//				  cout << endl << endl << endl;
//			  }
//		  }
			  //if(my_rank==0){ cout << "bp11" << "   ->   ImX_sim_band.rows()="<< ImX_sim_band.rows() << ", ImX_sim_band.cols()="<< ImX_sim_band.cols() << endl << endl;}
//			  ImX_sim->get_rasterBands()[ipp]->bandDataMatD = ImX_sim_band;

			  ImX_sim->get_rasterBands()[iG]->bandDataMatD = ImX_sim_band;
			  //cout << "bp8" << "iG="<<iG << endl;
			  //if(my_rank==0){ cout << "bp12" << "   ->   ImX_sim->get_rasterBands()[iG]->bandDataMatD.rows()=" << ImX_sim->get_rasterBands()[iG]->bandDataMatD.rows() << ", ImX_sim->get_rasterBands()[iG]->bandDataMatD.cols()=" << ImX_sim->get_rasterBands()[iG]->bandDataMatD.cols() << endl << endl;}
		  }
	  }else{
		  if(my_rank==0){
			  cout << "use original (not simulated) ImX for dictionary learning.." << endl;
		  }
		  for(int ipp=0; ipp<glPrms->numProbPerPatch; ipp++){
			  ImX_sim->get_rasterBands()[ipp]->bandDataMatD = ImX->get_rasterBands()[glPrms->myChX[ipp]]->get_bandDataMat()->cast<double>();
		  }
	  }
	  //cout << "bp9" << endl;
}

void calc_P_matrices(SpEOVectorD* P_lmd_vecs_loc, int **P_lmd_idx_bl_loc, SpEOMatrixI* P_lmd_idx_row_loc, int Ng, int *Nc_vec, int NChZ, int *idxChY,int my_rank){
	int iG, iChZ, iC, ig;

	for(ig=0; ig<Ng; ig++){
//					P_lmd_vecs_loc[ig] = SpEOVectorD::Ones(Nc_vec[ig])*glPrms->decMat_C(glPrms->myChX[ig], glPrms->myBundle[ig]);
		P_lmd_vecs_loc[ig] = SpEOVectorD::Ones(Nc_vec[ig]);
	}
	int col_idx=0, row_idx;
	for(ig=0; ig<Ng; ig++){
//					P_lmd_idx_bl_loc[ig] = new int[2];
		// =====> row_idx = index of the first HS channel in the current  group (in the simple case of Spectral Grouping, all groups (possibly except for the last group) are equally spaced, i.e. row_idx is linearly increasing)
//					row_idx = ig*(fSet->Nc - fSet->No);
		row_idx = idxChY[ig];
		// <=====
		P_lmd_idx_bl_loc[ig][0] = row_idx;
		P_lmd_idx_bl_loc[ig][1] = col_idx;
		col_idx += Nc_vec[ig];
	}
	// for every row (each corresponding to one HS channel iChY) calculate the relevant corresponding block indexes (first row) and the corresponding entry indexes relative to the corresponding block's origin each starting at 0.
	// first the number of non-trivial entries in each row (iChZ=0,...,NChZ-1) needs to be counted
	SpEOVectorI entries_cnt     = SpEOVectorI::Zero(NChZ);
	SpEOVectorI entries_divider = SpEOVectorI::Zero(NChZ);
	bool incr_entries_divider = false;
	bool iG_done[Ng];
	for(iG=0; iG<Ng; iG++){
		iG_done[iG]=false;
	}
	for(ig=0; ig<Ng; ig++){
		if(!iG_done[ig]){
			incr_entries_divider = true;
			iG_done[ig] = true;
		}else{
			incr_entries_divider = false;
		}
		for(iC=0; iC<Nc_vec[ig]; iC++){
			iChZ = P_lmd_idx_bl_loc[ig][0]+iC;
			if(incr_entries_divider){
				entries_divider(iChZ)++;
			}
			entries_cnt(iChZ)++;
		}
	}
	for(iChZ=0; iChZ<NChZ; iChZ++){
		P_lmd_idx_row_loc[iChZ] = SpEOMatrixI::Zero(2,entries_cnt.coeff(iChZ));
	}
	SpEOVectorI entries_cnt_tmp = SpEOVectorI::Zero(NChZ);
	for(ig=0; ig<Ng; ig++){
		for(iC=0; iC<Nc_vec[ig]; iC++){
			iChZ = P_lmd_idx_bl_loc[ig][0]+iC;
				// store block index
				P_lmd_idx_row_loc[iChZ](0,entries_cnt_tmp(iChZ)) = ig;//iG;
				// store coefficient index relative to corresponding block origin
				P_lmd_idx_row_loc[iChZ](1,entries_cnt_tmp(iChZ)) = iC;
			P_lmd_vecs_loc[ig](iC) /= entries_divider(iChZ);//entries_cnt(iChZ);
			entries_cnt_tmp(iChZ)++;
		}
	}




//	if(my_rank==0){
//		//
//		// print P_lmd to check for correctness
//		///*
//		cout << endl
//			 << "###########################" << endl
//			 << "# print P_lmd block-wise  #" << endl
//			 << "###########################" << endl;
//		for(iG=0; iG<Ng; iG++){
//			cout << "iG=" << iG << " and P_lmd_idx_bl_loc[iG][0]=row_idx="<< P_lmd_idx_bl_loc[iG][0]
//								<< " and P_lmd_idx_bl_loc[iG][1]=col_idx="<< P_lmd_idx_bl_loc[iG][1]
//				 << " P_lmd in this block = ";
//			for (iC=0; iC<Nc_vec[iG];iC++){
//				cout << " " << P_lmd_vecs_loc[iG](iC);
//			}
//			cout << endl;
//		}
//		//*/
//		///*
//		cout << endl
//			 << "###########################" << endl
//			 << "# print P_lmd row-wise #" << endl
//			 << "###########################" << endl;
//		for(iChZ=0; iChZ<NChZ; iChZ++){
//			cout << "iChZ=" << iChZ << "   ->   P_lmd_idx_row_loc[iChZ]= ";
//			for(int ii=0; ii<P_lmd_idx_row_loc[iChZ].cols(); ii++){
//				int block  = P_lmd_idx_row_loc[iChZ].coeff(0,ii);
//				int relidx = P_lmd_idx_row_loc[iChZ].coeff(1,ii);
//				cout << "    (block=" << block
//				   << ", relidx=" << relidx
//				   << ", val=" << P_lmd_vecs_loc[block](relidx)
//				   << ") ";
//			}
//			cout << endl;
//		}
//		//*/
//		///*
//		cout << endl
//			 << "###########################" << endl
//			 << "# print P_lmd block-wise  #" << endl
//			 << "###########################" << endl;
//		for(ig=0; ig<Ng; ig++){
//			cout << "ig=" << ig << " and P_lmd_idx_bl_loc[ig][0]=row_idx="<< P_lmd_idx_bl_loc[ig][0]
//								<< " and P_lmd_idx_bl_loc[ig][1]=col_idx="<< P_lmd_idx_bl_loc[ig][1]
//				 << " P_lmd in this block = ";
//			for (iC=0; iC<Nc_vec[ig];iC++){
//				cout << " " << P_lmd_vecs_loc[ig](iC);
//			}
//			cout << endl;
//		}
//		cout << endl
//			 << "###########################" << endl
//			 << "# print P_lmd row-wise #" << endl
//			 << "###########################" << endl;
//		for(iChZ=0; iChZ<NChZ; iChZ++){
//			cout << "iChZ=" << iChZ << " -> P_lmd_idx_row_loc[iChZ="<< iChZ <<"]= ";
//			for(int ii=0; ii<P_lmd_idx_row_loc[iChZ].cols(); ii++){
//				int block  = P_lmd_idx_row_loc[iChZ].coeff(0,ii);
//				int relidx = P_lmd_idx_row_loc[iChZ].coeff(1,ii);
//				cout << "(block=" << block
//				   << ", relidx=" << relidx
//				   << ", val=" << P_lmd_vecs_loc[block](relidx)
//				   << ") ";
//			}
//			cout << endl;
//		}
//		//*/
//		}

}


void CSG_corr_based_spectral_grouping(int &Ng_final, int *idxChY_final, int *Nc_vec_final, SpEOMatrixD &winY, double CCmin, int Nc_max, int my_rank){

	if(contains_nan(winY)){
		cout << endl << endl
			 << " EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE" << endl
			 << " E" << endl
			 << " E    NaN    found in     winY" << endl
			 << " E" << endl
			 << " EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE"
			 << endl << endl;
	}

	SpEOMatrixD winY_TMP = winY; //winY.block(0,0,24,winY.cols()); //
	bool printOutAll = false;

	int iChYU, iChYV;
	double q = -8;
	int NChY    = winY_TMP.rows();
	int winY_TMP_sz = winY_TMP.cols();
	//##################################################################################
	//#     1) Calculate correlation matrix of current patch across all channels       #
	//##################################################################################
	SpEOMatrixD CCmat = SpEOMatrixD::Ones(NChY,NChY);
	for(iChYU=0; iChYU<NChY; iChYU++){
		for(iChYV=iChYU+1; iChYV<NChY; iChYV++){
			CCmat(iChYU,iChYV) = winY_TMP.row(iChYU).dot(winY_TMP.row(iChYV));
			// make it symmetric
			CCmat(iChYV,iChYU) = CCmat(iChYU,iChYV);
		}
	}
	if(my_rank==0 && printOutAll){cout << "CCmat1 = " << endl << CCmat << endl << endl;}
	//##################################################################################
	//#     2) Threshold CC matrix by setting all entries to -1 that are smaller       #
	//#        than a pre-defined threshold CCmin                                      #
	//##################################################################################
	threshold_CCmat(CCmat, CCmin, q);
//	for(iChYU=0; iChYU<NChY; iChYU++){
//		for(iChYV=iChYU+1; iChYV<NChY; iChYV++){
//			if(CCmat(iChYU,iChYV) < CCmin){
//				CCmat(iChYU,iChYV) = q;
//			}
//			// make it symmetric
//			CCmat(iChYV,iChYU) = CCmat(iChYU,iChYV);
//		}
//	}
	if(my_rank==0 && printOutAll){cout << "CCmat2 = " << endl << CCmat << endl << endl;}
	//##################################################################################
	//#     3) “Clean up” CC_mat matrix by setting all entries (i,j) to -2 for which   #
	//#        there is an entry CC_mat(k,l) = -2  s.t.                                #
	//#        EITHER      l=j && j<k<i                                                #
	//#        OR          k=i && i<l<j                                                #
	//#        This step limits the relevant entries in the matrix to possibly         #
	//#        overlapping blocks along the diagonal.                                  #
	//##################################################################################
	cleanUp_CCmat(CCmat,q);
	if(my_rank==0 && printOutAll){cout << "CCmat3 = " << endl << CCmat << endl << endl;}
	//##################################################################################
	//#     4) Identify all blocks along the diagonal and associate every block with a #
	//#        group of adjacent bands. Moreover, distinguish two types of groups:     #
	//#        - Groups that are composed of fewer than or equal to NC,max bands       #
	//#          (colored in green) and                                                #
	//#        - Groups that are composed of more than NC,max                          #
	//#          bands (“large” groups)                                                #
	//##################################################################################
	//--------------------
	//   calculate Ng
	//--------------------
	int Ng = calc_Ng(CCmat, q);
	if(my_rank==0 && printOutAll){cout << "Ng = " << Ng << endl << endl;}
	//--------------------
	//   calculate idxChY
	//   and Nc_vec
	//--------------------
	int idxChY[Ng], Nc_vec[Ng];
	calc_idxChY_and_Nc_vec(idxChY, Nc_vec, Ng, CCmat, q, 0);
	for(int iG=0; iG<Ng; iG++){
		if(my_rank==0 && printOutAll){cout << "###################" << endl
		     << "# group iG=" << iG   << endl
		     << "###################" << endl
		     << "idxChY[iG] = " << idxChY[iG] << endl
		     << "Nc_vec[iG] = " << Nc_vec[iG] << endl
		     << "CCmat.block(idxChY[iG],idxChY[iG],Nc_vec[iG],Nc_vec[iG]) = " << endl
		     << CCmat.block(idxChY[iG],idxChY[iG],Nc_vec[iG],Nc_vec[iG])
		     << endl << endl;}
	}

	//##################################################################################
	//#     5) Split each large group into possibly overlapping smaller sub-groups as  #
	//#        follows:                                                                #
	//#        a) Refine the pre-defined threshold CCmin (for each group individually) #
	//#           by calculating the following three potential thresholds of which the #
	//#           smallest one will be selected (as a minimum requirement so to say):  #
	//#           Let A be the number of all entries in the upper half (incl. the      #
	//#           diagonal) of the symmetric current large block in CC_mat.            #
	//#           i)  Calculating CCmin,1: Let B<A be the number of those entries      #
	//#               (i,j) in the upper half (incl. the diagonal) of the symmetric    #
	//#               current large block in CC_mat that are less than or equal to     #
	//#               NC,max entries away from the diagonal, i.e. i <= j < i+NC,max or #
	//#               j <= i < j+NC,max. Next, put this upper part of the current      #
	//#               block in a vector VEC (of length A) and sort it according to its #
	//#               entries (correlation values) in a descending manner.             #
	//#               Set CCmin,1 = VEC(B). In other words, CCmin,1 is the Bth highest #
	//#               correlation value in the upper half of the current correlation   #
	//#               sub-matrix.                                                      #
	//#           ii) Calculating CCmin,2: CCmin,2 = VEC(k*), where k* is calculated   #
	//#               via                                                              #
	//#               sum_k=1^(k*) VEC(k) approx.= (1/10)*sum_k=1^A VEC(k)             #
	//#           iii) Calculating CCmin,3: CCmin,3 = CCmin + (2/3)*(1-CCmin)          #
	//#        b) - Select the smallest of the three potential thresholds. This way all#
	//#             of the conditions 5)a)i)-iii) are sufficient, but only the weakest #
	//#             one is necessary (which one can vary between datasets, patches and #
	//#             even blocks).                                                      #
	//#           - Threshold the sub-matrix like in step 2)                           #
	//#           - Clean up the resulting matrix like in step 3)                      #
	//#           - Eliminate single-band “groups” by grouping every possibly          #
	//#             occurring single band 1<=i<=NChY with the neighboring band j=i-1   #
	//#             or j=i+1 that has the higher correlation with channel i.           #
	//#                                                                                #
	//##################################################################################
	int iG;
	int num_large_grps = 0;
	for(iG=0; iG<Ng; iG++){
		if (Nc_vec[iG]>Nc_max){
			num_large_grps++;
		}
	}


	int num_sml_grps = Ng-num_large_grps;
	int iG_lrg_all[num_large_grps];
	int iG_sml_all[num_sml_grps];
	int Ng_grp[num_large_grps];
	int **idxChY_grp = new int*[num_large_grps];
	int **Nc_vec_grp = new int*[num_large_grps];

	int iG_lrg =- 1;
	int iG_sml =- 1;


	if(my_rank==0 && printOutAll){cout << "*******************************************************" << endl
		 << "*                                                     *" << endl
		 << "*                                                     *" << endl
		 << "*                                                     *" << endl;}
	for(iG=0; iG<Ng; iG++){
		if (Nc_vec[iG]>Nc_max){
			if(my_rank==0 && printOutAll){cout << "*********************" << endl
				 << "*    iG="<<iG<<"            *" << endl
				 << "*                   *" << endl
			     << "*                   *" << endl;}
			iG_lrg ++;
			iG_lrg_all[iG_lrg] = iG;
			SpEOMatrixD CCmat_grp = CCmat.block(idxChY[iG],idxChY[iG],Nc_vec[iG],Nc_vec[iG]);
			if(my_rank==0 && printOutAll){cout << "CCmat_grp = " << endl
				 << CCmat_grp << endl << endl;}
			double CCmin1, CCmin2, CCmin3;
			calc_refined_CCmins(CCmin1, CCmin2, CCmin3, CCmat_grp, CCmin, Nc_max);
			if(my_rank==0 && printOutAll){cout << "CCmin = " << CCmin << endl
			     << "CCmin1 = " << CCmin1 << endl
			     << "CCmin2 = " << CCmin2 << endl
			     << "CCmin3 = " << CCmin3 << endl;}
			double CCmin_grp = min(CCmin1, CCmin2);
			CCmin_grp = min(CCmin_grp, CCmin3);
			if(my_rank==0 && printOutAll){cout << "CCmin_grp = " << CCmin_grp << endl << endl;}
			threshold_CCmat(CCmat_grp, CCmin_grp, q);
			if(my_rank==0 && printOutAll){cout << "thresholded CCmat_grp = " << endl
							 << CCmat_grp << endl << endl;}
			cleanUp_CCmat(CCmat_grp,q);
			if(my_rank==0 && printOutAll){cout << "thresholded cleaned CCmat_grp = " << endl
										 << CCmat_grp << endl << endl;}
			Ng_grp[iG_lrg] = calc_Ng(CCmat_grp, q);
			idxChY_grp[iG_lrg] = new int[Ng_grp[iG_lrg]];
			Nc_vec_grp[iG_lrg] = new int[Ng_grp[iG_lrg]];
			calc_idxChY_and_Nc_vec(idxChY_grp[iG_lrg], Nc_vec_grp[iG_lrg], Ng_grp[iG_lrg], CCmat_grp, q, idxChY[iG]);

			if(my_rank==0 && printOutAll){cout << "Ng_grp[iG_lrg="<<iG_lrg<<"] = " << Ng_grp[iG_lrg] << endl;}
			for(int iCh=0; iCh<Ng_grp[iG_lrg]; iCh++){
				if(my_rank==0 && printOutAll){cout << "  iCh=" << iCh << "  -  "
					 << "idxChY_grp[iG_lrg][iCh] = " << idxChY_grp[iG_lrg][iCh]
					 << "  -  "
					 << "Nc_vec_grp[iG_lrg][iCh] = " << Nc_vec_grp[iG_lrg][iCh]
					 << endl;}
			}
			grp_sngl_bnd_grps_with_high_corr_nbrs(Nc_vec_grp, idxChY_grp, Ng_grp, iG_lrg, CCmat);

			if(my_rank==0 && printOutAll){cout << endl << "After single band groups were groupt with high correlated neighbors:" << endl;}
			for(int iCh=0; iCh<Ng_grp[iG_lrg]; iCh++){
				if(my_rank==0 && printOutAll){cout << "  iCh=" << iCh << "  -  "
					 << "idxChY_grp[iG_lrg][iCh] = " << idxChY_grp[iG_lrg][iCh]
					 << "  -  "
					 << "Nc_vec_grp[iG_lrg][iCh] = " << Nc_vec_grp[iG_lrg][iCh]
					 << endl;}
			}

			if(my_rank==0 && printOutAll){cout << "*                   *" << endl
				 << "*                   *" << endl
				 << "*********************" << endl;}

		}else{
			iG_sml++;
			iG_sml_all[iG_sml] = iG;
		}
	}
	if(my_rank==0 && printOutAll){cout << "*                                                     *" << endl
		 << "*                                                     *" << endl
		 << "*                                                     *" << endl
		 << "*******************************************************" << endl;}

//	int Ng_grp[num_large_grps];
//	int **idxChY_grp = new int*[num_large_grps];
//	int **Nc_vec_grp = new int*[num_large_grps];
	if(my_rank==0 && printOutAll){cout << "*******************************************************" << endl
		 << "*                                                     *" << endl
		 << "*                                                     *" << endl
		 << "*                                                     *" << endl;}
	for(iG_lrg=0; iG_lrg<num_large_grps; iG_lrg++){
		if(my_rank==0 && printOutAll){cout << "iG_lrg=" << iG_lrg << ", iG_lrg_all[iG_lrg]=" << iG_lrg_all[iG_lrg] << endl;}
		for(int iG_sub=0; iG_sub<Ng_grp[iG_lrg]; iG_sub++){
			if(my_rank==0 && printOutAll){cout << "iG_sub=" << iG_sub << "  -  "
				 << "idxChY_grp[iG_lrg][iG_sub] = " << idxChY_grp[iG_lrg][iG_sub]
				 << "  -  "
				 << "Nc_vec_grp[iG_lrg][iG_sub] = " << Nc_vec_grp[iG_lrg][iG_sub]
				 << endl;}
		}

	}
	if(my_rank==0 && printOutAll){cout << endl;}
	for(iG_sml=0; iG_sml<num_sml_grps; iG_sml++){
		if(my_rank==0 && printOutAll){cout << "iG_sml=" << iG_sml << ", iG_sml_all[iG_sml]=" << iG_sml_all[iG_sml] << endl;}
	}

	if(my_rank==0 && printOutAll){cout << "*                                                     *" << endl
		 << "*                                                     *" << endl
		 << "*                                                     *" << endl
		 << "*******************************************************" << endl;}

	//##################################################################################
	//#     6) Remove redundant groups that might have been produced by individually   #
	//#        splitting possibly overlapping large groups in step. More precisely,    #
	//#        remove every group that is entirely contained in another group.         #
	//##################################################################################

	int Ng_all_grps = num_sml_grps;
	for(iG_lrg=0; iG_lrg<num_large_grps; iG_lrg++){
		Ng_all_grps += Ng_grp[iG_lrg];
	}

	int all_grps_Nc_vec[Ng_all_grps];
	std::vector<argsort_pair> all_grps_idxChY(Ng_all_grps);
//	std::vector<argsort_pair> all_grps_Nc_vec(Ng_all_grps);

	if(my_rank==0 && printOutAll){cout << endl
		     << "&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&" << endl
		     << "&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&" << endl;}
	int iG_all_grps = -1;
	for(iG_sml=0; iG_sml<num_sml_grps; iG_sml++){
		iG_all_grps++;
		all_grps_idxChY[iG_all_grps].first  = iG_all_grps;
		all_grps_idxChY[iG_all_grps].second = idxChY[iG_sml_all[iG_sml]];
		all_grps_Nc_vec[iG_all_grps] = Nc_vec[iG_sml_all[iG_sml]];
		if(my_rank==0 && printOutAll){cout << "&  iG_all_grps=" << iG_all_grps << ": idxChY[iG_sml_all[iG_sml=" << iG_sml << "]=" << iG_sml_all[iG_sml] << "] = " << idxChY[iG_sml_all[iG_sml]]
		                                                                                  << "  #  " << Nc_vec[iG_sml_all[iG_sml]] << " = Nc_vec" << endl;}
	}
	if(my_rank==0 && printOutAll){cout << endl << "-------------------------------------" << endl;}

	for(iG_lrg=0; iG_lrg<num_large_grps; iG_lrg++){
		//if(my_rank==0 && printOutAll){cout << "iG_lrg=" << iG_lrg << ", iG_lrg_all[iG_lrg]=" << iG_lrg_all[iG_lrg] << endl;}
		for(int iG_sub=0; iG_sub<Ng_grp[iG_lrg]; iG_sub++){
//			if(my_rank==0 && printOutAll){cout << "iG_sub=" << iG_sub << "  -  "
//				 << "idxChY_grp[iG_lrg][iG_sub] = " << idxChY_grp[iG_lrg][iG_sub]
//				 << "  -  "
//				 << "Nc_vec_grp[iG_lrg][iG_sub] = " << Nc_vec_grp[iG_lrg][iG_sub]
//				 << endl;}
			iG_all_grps++;
			all_grps_idxChY[iG_all_grps].first  = iG_all_grps;
			all_grps_idxChY[iG_all_grps].second = idxChY_grp[iG_lrg][iG_sub];
			all_grps_Nc_vec[iG_all_grps] = Nc_vec_grp[iG_lrg][iG_sub];
			if(my_rank==0 && printOutAll){cout << "&  iG_all_grps=" << iG_all_grps << ": idxChY_grp[iG_lrg=" << iG_lrg << "][iG_sub=" << iG_sub << "] = " << idxChY_grp[iG_lrg][iG_sub]
			                                                                      << "  #  " << Nc_vec_grp[iG_lrg][iG_sub] << " = Nc_vec" << endl;}
		}
		if(my_rank==0 && printOutAll){cout << endl;}
	}
	// Ng_all_grps==iG_all_grps+1 must be fulfilled after loop.
	if(Ng_all_grps!=iG_all_grps+1){
		cerr << "Ng_all_grps=="<< Ng_all_grps << " != iG_all_grps+1==" << iG_all_grps+1 << " !!! -> Fatal Error!" << endl << endl << endl;
		exit(2);
	}

	if(my_rank==0 && printOutAll){cout << endl
	     << "&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&" << endl
		 << "&   before sorting " << endl
		 << "&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&" << endl;}
	for(iG_all_grps=0; iG_all_grps<Ng_all_grps; iG_all_grps++){
		if(my_rank==0 && printOutAll){cout << "&  all_grps_idxChY[iG_all_grps=="<<iG_all_grps<<"] = ("   << all_grps_idxChY[iG_all_grps].first
															   << " ," << all_grps_idxChY[iG_all_grps].second
															   << ")  #  " << all_grps_Nc_vec[iG_all_grps] << " = Nc_vec"
															   << endl;}
	}
	std::sort(all_grps_idxChY.begin(), all_grps_idxChY.end(), argsort_comp);


	int grps_to_keep[Ng_all_grps];
	int num_grps_to_keep = Ng_all_grps;
	int iG_all_grps_run;
	int curr_frst_bnd;
	int iG_all_grps_biggest;
	iG_all_grps=0;
	while (iG_all_grps<Ng_all_grps){
		// index if largest of all groups starting at all_grps_idxChY[iG_all_grps].second
		iG_all_grps_biggest = iG_all_grps;
		curr_frst_bnd = all_grps_idxChY[iG_all_grps].second;
		iG_all_grps_run         = iG_all_grps+1;

		while (iG_all_grps_run < Ng_all_grps &&  (int)(all_grps_idxChY[iG_all_grps_run].second) == curr_frst_bnd){
			if(  all_grps_Nc_vec[all_grps_idxChY[iG_all_grps_run].first]
			   > all_grps_Nc_vec[all_grps_idxChY[iG_all_grps].first]){
				iG_all_grps_biggest = iG_all_grps_run;
			}
			iG_all_grps_run ++;
		}
		for(int iG_all_grps_run_tmp=iG_all_grps; iG_all_grps_run_tmp < iG_all_grps_run; iG_all_grps_run_tmp++){
			if(iG_all_grps_run_tmp == iG_all_grps_biggest){
				grps_to_keep[iG_all_grps_run_tmp] = true;
			}else{
				grps_to_keep[iG_all_grps_run_tmp] = false;
				num_grps_to_keep -= 1;
			}
		}
	    iG_all_grps = iG_all_grps_run;
	}
	if(my_rank==0 && printOutAll){cout << endl
		 << "&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&" << endl
		 << "&   after sorting and identifying redundant groups that START with the same band (keep only largest)" << endl
		 << "&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&" << endl;}
	for(iG_all_grps=0; iG_all_grps<Ng_all_grps; iG_all_grps++){
		if(my_rank==0 && printOutAll){cout << "&  all_grps_idxChY[iG_all_grps=="<<iG_all_grps<<"] = ("   << all_grps_idxChY[iG_all_grps].first
															   << " ," << all_grps_idxChY[iG_all_grps].second
															   << ")  #  " << all_grps_Nc_vec[all_grps_idxChY[iG_all_grps].first] << " = Nc_vec"
															   << " # " << grps_to_keep[iG_all_grps]
															   << endl;}
	}

	if(my_rank==0 && printOutAll){cout << endl
		 << "&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&" << endl
		 << "&   after identifying remaining redundant groups " << endl
		 << "&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&" << endl;}
	for(iG_all_grps=0; iG_all_grps<Ng_all_grps; iG_all_grps++){
		if(grps_to_keep[iG_all_grps]){
			// check if group contains other groups. If so, mark as "to be removed" / "not to keep"
			iG_all_grps_run = iG_all_grps+1;
			while(   iG_all_grps_run < Ng_all_grps
				  &&   all_grps_idxChY[iG_all_grps_run].second   <  all_grps_idxChY[iG_all_grps].second + (float)(all_grps_Nc_vec[all_grps_idxChY[iG_all_grps].first]) ){
				if(   grps_to_keep[iG_all_grps_run]
				   &&    all_grps_idxChY[iG_all_grps_run].second + (float)(all_grps_Nc_vec[all_grps_idxChY[iG_all_grps_run].first])
				      <= all_grps_idxChY[iG_all_grps].second     + (float)(all_grps_Nc_vec[all_grps_idxChY[iG_all_grps].first])      ){
					grps_to_keep[iG_all_grps_run] = false;
					num_grps_to_keep -=1;
				}
				iG_all_grps_run++;
			}
		}
	}
	for(iG_all_grps=0; iG_all_grps<Ng_all_grps; iG_all_grps++){
		if(my_rank==0 && printOutAll){cout << "&  all_grps_idxChY[iG_all_grps=="<<iG_all_grps<<"] = ("   << all_grps_idxChY[iG_all_grps].first
															   << " ," << all_grps_idxChY[iG_all_grps].second
															   << ")  #  " << all_grps_Nc_vec[all_grps_idxChY[iG_all_grps].first] << " = Nc_vec"
															   << " # " << grps_to_keep[iG_all_grps]
															   << endl;}
	}
	if(my_rank==0 && printOutAll){cout << endl
		 << "&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&" << endl
		 << "&   after removing redundant groups " << endl
		 << "&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&" << endl;}
	int Ng_new = num_grps_to_keep;
	int idxChY_new[Ng_new];
	int Nc_vec_new[Ng_new];

	int iG_new = -1;
	for(iG_all_grps=0; iG_all_grps<Ng_all_grps; iG_all_grps++){
		if(grps_to_keep[iG_all_grps]){
			iG_new ++;
			idxChY_new[iG_new] = (int)(all_grps_idxChY[iG_all_grps].second);
			Nc_vec_new[iG_new] = all_grps_Nc_vec[all_grps_idxChY[iG_all_grps].first];
		}
	}
	for(iG_new=0; iG_new<Ng_new; iG_new++){
		if(my_rank==0 && printOutAll){cout << "&  idxChY_new[iG_new=="<<iG_new<<"] = ("   << idxChY_new[iG_new]
															   << " ," << idxChY_new[iG_new]
															   << ")  #  " << Nc_vec_new[iG_new] << " = Nc_vec"
															   << endl;}
	}




	//##################################################################################
	//#   7.1) remove redundant center groups. I.e. groups that have both a right and  #
	//#        a left neighbor group which sufficiently overlap without the center     #
	//#        group                                                                   #
	//##################################################################################
	if(my_rank==0 && printOutAll){cout << endl << "################## 7.1 stated .." << endl;}
	bool rm_cntr_grps=true;
	if(rm_cntr_grps){
		double myEps=1e-6;
		SpEOVectorD o_fw = SpEOVectorD::Zero(Ng_new);
		SpEOVectorD o_bw = SpEOVectorD::Zero(Ng_new);
		SpEOVectorD Nc_half_ceil = SpEOVectorD::Zero(Ng_new);
		for(iG_new=0; iG_new<Ng_new-2; iG_new++){
			o_fw(iG_new) = idxChY_new[iG_new]+Nc_vec_new[iG_new] - idxChY_new[iG_new+2];
			o_bw(iG_new+2) = o_fw(iG_new);
		}
		for(iG_new=0; iG_new<Ng_new; iG_new++){
			Nc_half_ceil(iG_new) = ceil(0.5*(double)(Nc_vec_new[iG_new]));
		}
		SpEOVectorD o_fw_min_Nc_half_ceil = o_fw-Nc_half_ceil;
		SpEOVectorD o_bw_min_Nc_half_ceil = o_bw-Nc_half_ceil;
		while(o_fw_min_Nc_half_ceil.maxCoeff() > 0. || o_bw_min_Nc_half_ceil.maxCoeff() > 0.){
			// forward overlapping
			double fw_max = -1;
			for(iG_new=0; iG_new<Ng_new; iG_new++){
				if(    (o_fw(iG_new) - Nc_half_ceil(iG_new) > 0)
					&& (o_fw(iG_new)/Nc_vec_new[iG_new] > fw_max) ){
					fw_max = o_fw(iG_new)/Nc_vec_new[iG_new];
				}
			}
			int fw_max_cnt = 0;
			SpEOVectorI I_fw = SpEOVectorI::Zero(Ng_new); // that's a conservative memory allocation. It will contain only one element in most cases
			for(iG_new=0; iG_new<Ng_new; iG_new++){
				if ( ( o_fw(iG_new)-Nc_half_ceil(iG_new) > 0 )  &&  ( o_fw(iG_new)/Nc_vec_new[iG_new] > fw_max-myEps ) ){
					I_fw(fw_max_cnt) = iG_new;
					fw_max_cnt ++;
				}
			}
			// backward overlapping
			double bw_max = -1;
			for(iG_new=0; iG_new<Ng_new; iG_new++){
				if(    (o_bw(iG_new) - Nc_half_ceil(iG_new) > 0)
					&& (o_bw(iG_new)/Nc_vec_new[iG_new] > bw_max) ){
					bw_max = o_bw(iG_new)/Nc_vec_new[iG_new];
				}
			}
			int bw_max_cnt = 0;
			SpEOVectorI I_bw = SpEOVectorI::Constant(Ng_new,Ng_new); // that's a conservative memory allocation. It will contain only one element in most cases
			for(iG_new=0; iG_new<Ng_new; iG_new++){
				if ( ( o_bw(iG_new)-Nc_half_ceil(iG_new) > 0 )  &&  ( o_bw(iG_new)/Nc_vec_new[iG_new] > bw_max-myEps ) ){
					I_bw(bw_max_cnt) = iG_new;
					bw_max_cnt ++;
				}
			}
			if(my_rank==0 && printOutAll){
				cout << " ___________________________________________" << endl;
				cout << " |    Ng_new = " << Ng_new <<	endl
				     << " |    idxChY_new = ";
				for(int iG_tmp=0; iG_tmp<Ng_new; iG_tmp++){
					cout << idxChY_new[iG_tmp] << ", ";
				}
				cout << endl;
				cout << " |    Nc_vec_new = ";
				for(int iG_tmp=0; iG_tmp<Ng_new; iG_tmp++){
					cout << Nc_vec_new[iG_tmp] << ", ";
				}
				cout << endl
					 << " |    o_fw = " << o_fw.transpose() << endl
					 << " |    o_fw_min_Nc_half_ceil = " << o_fw_min_Nc_half_ceil.transpose() << endl
					 << " |   bw_max = " << bw_max << "   ->   I_bw = ";
				for(int ibla=0; ibla<bw_max_cnt; ibla++){
					cout << I_bw(ibla) << ", ";
				}
				cout << endl;
				cout << " |   fw_max = " << fw_max << "   ->   I_fw = ";
				for(int ibla=0; ibla<fw_max_cnt; ibla++){
					cout << I_fw(ibla) << ", ";
				}
				cout << endl;
			}
			int irun;
			if(bw_max_cnt==0 && fw_max_cnt==0){
				if(my_rank==0 && printOutAll){cout << "overlapping reduction completed successfully.. (bw_max_cnt==0 && fw_max_cnt==0)" << endl;}
			}else{
				if(my_rank==0 && printOutAll){cout << "bw_max_cnt = " << bw_max_cnt << ", fw_max_cnt=" << fw_max_cnt << endl;}
				if(fw_max >= bw_max){
					int iii = I_fw.maxCoeff();
					if(my_rank==0 && printOutAll){cout << "   ACTION: fw_max >= bw_max   &&   remove group iii+1=" << iii+1;}
					if(my_rank==0 && printOutAll){cout << "           (idxChY_new[iii+1]=" << idxChY_new[iii+1] << ", and  Nc_vec_new[iii+1]=" << Nc_vec_new[iii+1] << ")" << endl;}
					// remove group iii+1
					for(irun=iii+1; irun<Ng_new-1; irun++){
						idxChY_new[irun] = idxChY_new[irun+1];
						Nc_vec_new[irun] = Nc_vec_new[irun+1];
					}
					Ng_new -= 1;
				}else{
					int iii = I_bw.minCoeff();
					if(my_rank==0 && printOutAll){cout << "   ACTION: fw_max < bw_max   &&   remove group iii-1=" << iii-1;}
					if(my_rank==0 && printOutAll){cout << "           (idxChY_new[iii-1]=" << idxChY_new[iii-1] << ", and  Nc_vec_new[iii-1]=" << Nc_vec_new[iii-1] << ")" << endl;}
					// remove group iii-1
					for(irun=iii-1; irun<Ng_new-1; irun++){
						idxChY_new[irun] = idxChY_new[irun+1];
						Nc_vec_new[irun] = Nc_vec_new[irun+1];
					}
					Ng_new -= 1;
				}
				o_fw = SpEOVectorD::Zero(Ng_new);
				o_bw = SpEOVectorD::Zero(Ng_new);
				I_fw = SpEOVectorI::Zero(Ng_new);
				//I_bw = SpEOVectorI::Zero(Ng_new);
				I_bw = SpEOVectorI::Constant(Ng_new,Ng_new);
				for(iG_new=0; iG_new<Ng_new-2; iG_new++){
					o_fw(iG_new) = idxChY_new[iG_new]+Nc_vec_new[iG_new] - idxChY_new[iG_new+2];
					o_bw(iG_new+2) = o_fw(iG_new);
				}
	
				Nc_half_ceil = SpEOVectorD::Zero(Ng_new);
				for(iG_new=0; iG_new<Ng_new; iG_new++){
					Nc_half_ceil(iG_new) = ceil(0.5*(double)(Nc_vec_new[iG_new]));
				}
				o_fw_min_Nc_half_ceil = o_fw-Nc_half_ceil;
				o_bw_min_Nc_half_ceil = o_bw-Nc_half_ceil;
				if(my_rank==0 && printOutAll){
					cout << " |    ->  Ng_new = " << Ng_new <<	endl
						 << " |    ->  idxChY_new = ";
					for(int iG_tmp=0; iG_tmp<Ng_new; iG_tmp++){
						cout << idxChY_new[iG_tmp] << ", ";
					}
					cout << endl;
					cout << " |    ->  Nc_vec_new = ";
					for(int iG_tmp=0; iG_tmp<Ng_new; iG_tmp++){
						cout << Nc_vec_new[iG_tmp] << ", ";
					}
					cout << endl
						 << " |    ->  o_fw = " << o_fw.transpose() << endl
						 << " |    ->  o_fw_min_Nc_half_ceil = " << o_fw_min_Nc_half_ceil.transpose() << endl << endl;
				}
			}
		}
		if(my_rank==0 && printOutAll){
			for(iG_new=0; iG_new<Ng_new; iG_new++){
				cout << "&  idxChY_new[iG_new=="<<iG_new<<"] = (" << idxChY_new[iG_new]
										  << " ," << idxChY_new[iG_new]
										  << ")  #  " << Nc_vec_new[iG_new] << " = Nc_vec"
										  << endl;
			}
			cout << endl
				 << "&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&" << endl
				 << "&   after removing center groups &" << endl
				 << "&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&" << endl;
		}
	}	

	//##################################################################################
	//#   7.2) Reduce possible high overlapping and redundancy of neighboring groups   #
	//#        by an iterative reduction of groups in terms of both the size reduction #
	//#        and elimination, following Algorithm1.                                  #
	//##################################################################################
	if(my_rank==0 && printOutAll){cout << "################## 7.2 stated .." << endl;}
	bool reduce_larger_grps = true;
	//bool remove_redundant_grps_apriori = true;
	double myEps=1e-6;
	SpEOVectorD o_fw = SpEOVectorD::Zero(Ng_new);
	SpEOVectorD o_bw = SpEOVectorD::Zero(Ng_new);
	SpEOVectorD Nc_half_ceil = SpEOVectorD::Zero(Ng_new);
	for(iG_new=0; iG_new<Ng_new-1; iG_new++){
		o_fw(iG_new) = idxChY_new[iG_new]+Nc_vec_new[iG_new] - idxChY_new[iG_new+1];
		o_bw(iG_new+1) = o_fw(iG_new);
	}
	for(iG_new=0; iG_new<Ng_new; iG_new++){
		Nc_half_ceil(iG_new) = ceil(0.5*(double)(Nc_vec_new[iG_new]));
	}
	SpEOVectorD o_fw_min_Nc_half_ceil = o_fw-Nc_half_ceil;
	SpEOVectorD o_bw_min_Nc_half_ceil = o_bw-Nc_half_ceil;
	while(o_fw_min_Nc_half_ceil.maxCoeff() > 0. || o_bw_min_Nc_half_ceil.maxCoeff() > 0.){
		// forward overlapping
		double fw_max = -1;
		for(iG_new=0; iG_new<Ng_new; iG_new++){
			if(    (o_fw(iG_new) - Nc_half_ceil(iG_new) > 0)
				&& (o_fw(iG_new)/Nc_vec_new[iG_new] > fw_max) ){
				fw_max = o_fw(iG_new)/Nc_vec_new[iG_new];
			}
		}
		int fw_max_cnt = 0;
		SpEOVectorI I_fw = SpEOVectorI::Zero(Ng_new); // that's a conservative memory allocation. It will contain only one element in most cases
		for(iG_new=0; iG_new<Ng_new; iG_new++){
			if ( ( o_fw(iG_new)-Nc_half_ceil(iG_new) > 0 )  &&  ( o_fw(iG_new)/Nc_vec_new[iG_new] > fw_max-myEps ) ){
				I_fw(fw_max_cnt) = iG_new;
				fw_max_cnt ++;
			}
		}
		// backward overlapping
		double bw_max = -1;
		for(iG_new=0; iG_new<Ng_new; iG_new++){
			if(    (o_bw(iG_new) - Nc_half_ceil(iG_new) > 0)
				&& (o_bw(iG_new)/Nc_vec_new[iG_new] > bw_max) ){
				bw_max = o_bw(iG_new)/Nc_vec_new[iG_new];
			}
		}
		int bw_max_cnt = 0;
		SpEOVectorI I_bw = SpEOVectorI::Constant(Ng_new,Ng_new); // that's a conservative memory allocation. It will contain only one element in most cases
		for(iG_new=0; iG_new<Ng_new; iG_new++){
			if ( ( o_bw(iG_new)-Nc_half_ceil(iG_new) > 0 )  &&  ( o_bw(iG_new)/Nc_vec_new[iG_new] > bw_max-myEps ) ){
				I_bw(bw_max_cnt) = iG_new;
				bw_max_cnt ++;
			}
		}


//		for(iG_new=0; iG_new<Ng_new; iG_new++){
//			if(my_rank==0 && printOutAll){cout << "HALLO &  idxChY_new[iG_new=="<<iG_new<<"] = ("   << idxChY_new[iG_new]
//																   << " ," << idxChY_new[iG_new]
//																   << ")  #  " << Nc_vec_new[iG_new] << " = Nc_vec"
//																   << endl;}
//		}

		if(my_rank==0 && printOutAll){
			cout << " ___________________________________________" << endl;
			cout << " |    Ng_new = " << Ng_new <<	endl
			     << " |    idxChY_new = ";
			for(int iG_tmp=0; iG_tmp<Ng_new; iG_tmp++){
				cout << idxChY_new[iG_tmp] << ", ";
			}
			cout << endl;
			cout << " |    Nc_vec_new = ";
			for(int iG_tmp=0; iG_tmp<Ng_new; iG_tmp++){
				cout << Nc_vec_new[iG_tmp] << ", ";
			}
			cout << endl
				 << " |    o_fw = " << o_fw.transpose() << endl
				 << " |    o_fw_min_Nc_half_ceil = " << o_fw_min_Nc_half_ceil.transpose() << endl
				 << " |   bw_max = " << bw_max << "   ->   I_bw = ";
			for(int ibla=0; ibla<bw_max_cnt; ibla++){
				cout << I_bw(ibla) << ", ";
			}
			cout << endl;
			cout << " |   fw_max = " << fw_max << "   ->   I_fw = ";
			for(int ibla=0; ibla<fw_max_cnt; ibla++){
				cout << I_fw(ibla) << ", ";
			}
			cout << endl;
		}
		int irun;
		if(bw_max_cnt==0 && fw_max_cnt==0){
			if(my_rank==0 && printOutAll){cout << "overlapping reduction completed successfully.. (bw_max_cnt==0 && fw_max_cnt==0)" << endl;}
		}else{
			if(my_rank==0 && printOutAll){cout << "bw_max_cnt = " << bw_max_cnt << ", fw_max_cnt=" << fw_max_cnt << endl;}
			if(fw_max >= bw_max){
				int iii = I_fw.maxCoeff();
				if( (iii>0) && ((idxChY_new[iii-1]+Nc_vec_new[iii-1]-idxChY_new[iii+1])==(o_fw(iii)-1)) ){
					if(my_rank==0 && printOutAll){cout << "   ACTION: fw_max >= bw_max   &&   remove group iii=" << iii;}
					if(my_rank==0 && printOutAll){cout << "           (idxChY_new[iii]=" << idxChY_new[iii] << ", and  Nc_vec_new[iii]=" << Nc_vec_new[iii] << ")" << endl;}
					// remove group iii
					for(irun=iii; irun<Ng_new-1; irun++){
						idxChY_new[irun] = idxChY_new[irun+1];
						Nc_vec_new[irun] = Nc_vec_new[irun+1];
					}
					Ng_new -= 1;
				}else{
					if(reduce_larger_grps && iii<Ng_new-1){
						// shrink group iii+1
						if(my_rank==0 && printOutAll){cout << "   ACTION: fw_max >= bw_max   &&   shrink group iii+1=" << iii << "+1 by setting  Nc_vec_new[iii+1] -= 1 AND idxChY_new[iii+1] += 1";}
						if(my_rank==0 && printOutAll){cout << "           (idxChY_new[iii]=" << idxChY_new[iii] << ", and  Nc_vec_new[iii]=" << Nc_vec_new[iii] << ")" << endl;}
						Nc_vec_new[iii+1] -= 1;
						idxChY_new[iii+1] += 1;
					}else{
						// shrink group iii
						if(my_rank==0 && printOutAll){cout << "   ACTION: fw_max >= bw_max   &&   shrink group iii=" << iii << " by setting  Nc_vec_new[iii] -= 1";}
						if(my_rank==0 && printOutAll){cout << "           (idxChY_new[iii]=" << idxChY_new[iii] << ", and  Nc_vec_new[iii]=" << Nc_vec_new[iii] << ")" << endl;}
						Nc_vec_new[iii] -= 1;
					}
				}
			}else{
				int iii = I_bw.minCoeff();
				if( (iii<Ng_new-1) && ((idxChY_new[iii-1]+Nc_vec_new[iii-1]-idxChY_new[iii+1])==(o_bw(iii)-1)) ){
					if(my_rank==0 && printOutAll){cout << "   ACTION: fw_max < bw_max   &&   remove group iii=" << iii;}
					if(my_rank==0 && printOutAll){cout << "           (idxChY_new[iii]=" << idxChY_new[iii] << ", and  Nc_vec_new[iii]=" << Nc_vec_new[iii] << ")" << endl;}
					// remove group iii
					for(irun=iii; irun<Ng_new-1; irun++){
						idxChY_new[irun] = idxChY_new[irun+1];
						Nc_vec_new[irun] = Nc_vec_new[irun+1];
					}
					Ng_new -= 1;
				}else{
					if(reduce_larger_grps && iii>0){
						// shrink group iii-1
						if(my_rank==0 && printOutAll){cout << "   ACTION: fw_max < bw_max   &&   shrink group iii-1=" << iii << "-1 by setting  Nc_vec_new[iii-1] -= 1";}
						if(my_rank==0 && printOutAll){cout << "           (idxChY_new[iii]=" << idxChY_new[iii] << ", and  Nc_vec_new[iii]=" << Nc_vec_new[iii] << ")" << endl;}
						Nc_vec_new[iii-1] -= 1;
					}else{
						// shrink group iii
						if(my_rank==0 && printOutAll){cout << "   ACTION: fw_max < bw_max   &&   shrink group iii=" << iii << " by setting  Nc_vec_new[iii] -= 1;  AND  idxChY_new[iii] += 1;";}
						if(my_rank==0 && printOutAll){cout << "           (idxChY_new[iii]=" << idxChY_new[iii] << ", and  Nc_vec_new[iii]=" << Nc_vec_new[iii] << ")" << endl;}
						Nc_vec_new[iii] -= 1;
						idxChY_new[iii] += 1;
					}
				}
			}
			o_fw = SpEOVectorD::Zero(Ng_new);
			o_bw = SpEOVectorD::Zero(Ng_new);
			I_fw = SpEOVectorI::Zero(Ng_new);
			//I_bw = SpEOVectorI::Zero(Ng_new);
			I_bw = SpEOVectorI::Constant(Ng_new,Ng_new);
			for(iG_new=0; iG_new<Ng_new-1; iG_new++){
				o_fw(iG_new) = idxChY_new[iG_new]+Nc_vec_new[iG_new] - idxChY_new[iG_new+1];
				o_bw(iG_new+1) = o_fw(iG_new);
			}

			Nc_half_ceil = SpEOVectorD::Zero(Ng_new);
			for(iG_new=0; iG_new<Ng_new; iG_new++){
				Nc_half_ceil(iG_new) = ceil(0.5*(double)(Nc_vec_new[iG_new]));
			}
			o_fw_min_Nc_half_ceil = o_fw-Nc_half_ceil;
			o_bw_min_Nc_half_ceil = o_bw-Nc_half_ceil;
			if(my_rank==0 && printOutAll){
				cout << " |    ->  Ng_new = " << Ng_new <<	endl
					 << " |    ->  idxChY_new = ";
				for(int iG_tmp=0; iG_tmp<Ng_new; iG_tmp++){
					cout << idxChY_new[iG_tmp] << ", ";
				}
				cout << endl;
				cout << " |    ->  Nc_vec_new = ";
				for(int iG_tmp=0; iG_tmp<Ng_new; iG_tmp++){
					cout << Nc_vec_new[iG_tmp] << ", ";
				}
				cout << endl
					 << " |    ->  o_fw = " << o_fw.transpose() << endl
					 << " |    ->  o_fw_min_Nc_half_ceil = " << o_fw_min_Nc_half_ceil.transpose() << endl << endl;
			}
		}

//		if(~isempty(I_bw) || ~isempty(I_fw))
//			 if(fw_max >= bw_max)
//				iii = max(I_fw);
//				if (iii>1)  && ((first_band_in_group(iii-1)+Nc(iii-1)-first_band_in_group(iii+1))==(o_fw(iii)-1));
//					% remove group iii
//					Ng = Ng-1;
//					first_band_in_group = first_band_in_group([1:iii-1,iii+1:end]);
//					Nc = Nc([1:iii-1,iii+1:end]);
//				else
//					% shrink group iii
//					Nc(iii) = Nc(iii)-1;
//				end
//			 else
//				iii = min(I_bw);
//				if (iii<Ng) && ((first_band_in_group(iii-1)+Nc(iii-1)-first_band_in_group(iii+1))==(o_bw(iii)-1));
//					% remove group iii
//					Ng = Ng-1;
//					first_band_in_group = first_band_in_group([1:iii-1,iii+1:end]);
//					Nc = Nc([1:iii-1,iii+1:end]);
//				else
//					% shrink group iii
//					Nc(iii) = Nc(iii)-1;
//					first_band_in_group(iii) = first_band_in_group(iii)+1;
//				end
//			 end
//			 o_fw = zeros(1,Ng);
//			 o_bw = zeros(1,Ng);
//			 I_fw = [];
//			 I_bw = [];
//			 for iG=1:Ng-1
//				 o_fw(iG) = first_band_in_group(iG)+Nc(iG)-first_band_in_group(iG+1);
//			 end
//			 o_bw(2:Ng) = o_fw(1:Ng-1);
//			 Nc_half_ceil = ceil(Nc/2);
//		 else
//			 disp('overlapping reduction completed successfully!')
//			 break
//		 end


	}
	if(my_rank==0 && printOutAll){
		for(iG_new=0; iG_new<Ng_new; iG_new++){
			cout << "&  idxChY_new[iG_new=="<<iG_new<<"] = (" << idxChY_new[iG_new]
									  << " ," << idxChY_new[iG_new]
									  << ")  #  " << Nc_vec_new[iG_new] << " = Nc_vec"
									  << endl;
		}
		cout << endl
			 << "&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&" << endl
			 << "&   after shrinking and possibly removing highly overlapping groups " << endl
			 << "&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&" << endl;
		cout << endl
			 << "&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&" << endl
			 << "&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&" << endl;
	}




	// only temporarily:
//	Ng_final     = Ng_new;
//	for(iG_new=0; iG_new<Ng_new; iG_new++){
//		idxChY_final[iG_new] = idxChY_new[iG_new];
//		Nc_vec_final[iG_new] = Nc_vec_new[iG_new];
//	}




//	std::vector<argsort_pair> CCmat_grp_vec_for_sorting(A);
//	for(int ii=0; ii<A; ii++){
//		CCmat_grp_vec_for_sorting[ii].first = ii;
//		CCmat_grp_vec_for_sorting[ii].second = CCmat_grp_vec[ii];
//	}
//	std::sort(CCmat_grp_vec_for_sorting.begin(), CCmat_grp_vec_for_sorting.end(), argsort_comp);
//	for(int ii=0; ii<A; ii++){
//		cout << "CCmat_grp_vec_for_sorting[ii=="<<ii<<"] = ("   << CCmat_grp_vec_for_sorting[ii].first
//				                                       << " , " << CCmat_grp_vec_for_sorting[ii].second
//				                                       << ")"   << endl;
//	}



	for(iG_lrg=0; iG_lrg<num_large_grps; iG_lrg++){
		delete[] idxChY_grp[iG_lrg];
		delete[] Nc_vec_grp[iG_lrg];
	}
	delete[] idxChY_grp;
	delete[] Nc_vec_grp;


	Ng_final = Ng_new;
	int iG_fin;
	for(iG_fin=0; iG_fin<Ng_final; iG_fin++){
		idxChY_final[iG_fin] = idxChY_new[iG_fin];
		Nc_vec_final[iG_fin] = Nc_vec_new[iG_fin];
	}

//	Ng_final = 4;
//	idxChY_final[0] = 0;
//	idxChY_final[1] = 40;
//	idxChY_final[2] = 80;
//	idxChY_final[3] = 120;
//	Nc_vec_final[0] = 41;
//	Nc_vec_final[1] = 41;
//	Nc_vec_final[2] = 41;
//	Nc_vec_final[3] = 40;

}

void threshold_CCmat(SpEOMatrixD &CCmat, double CCmin, double q){
	int NChY = CCmat.rows();
	for(int iChYU=0; iChYU<NChY; iChYU++){
		for(int iChYV=iChYU+1; iChYV<NChY; iChYV++){
			if(CCmat(iChYU,iChYV) < CCmin){
				CCmat(iChYU,iChYV) = q;
			}
			// make it symmetric
			CCmat(iChYV,iChYU) = CCmat(iChYU,iChYV);
		}
	}
}

void cleanUp_CCmat(SpEOMatrixD &CCmat, double q){
	double myeps = 0.0001;
	int iChYV, iChYU2, iChYV2;
	bool jumpToNextU;
	int NChY = CCmat.rows();
	for(int iChYU=0; iChYU<NChY-1; iChYU++){
		jumpToNextU = false;
		iChYV=iChYU+1;
		while(iChYV<NChY && !jumpToNextU){ // while(iChYV<NChY-1 && !jumpToNextU){
			if (CCmat(iChYU, iChYV) < q+myeps){
				SpEOMatrixD q_mat_tmp = SpEOMatrixD::Constant(iChYU+1, NChY-iChYV,q);
				CCmat.block(0,iChYV, iChYU+1, NChY-iChYV) = q_mat_tmp;
				CCmat.block(iChYV,0, NChY-iChYV, iChYU+1) = q_mat_tmp.transpose();
				jumpToNextU = true;
			}
			iChYV++;
		}
	}
}
int calc_Ng(SpEOMatrixD &CCmat, double q){
	double myeps = 0.0001;
	int iChYU, iChYV;
	int NChY = CCmat.rows();
	int Ng = 0;
	int iChYV_last_grp = -1;
	for(iChYU=0; iChYU<NChY; iChYU++){
		iChYV = iChYU;
		while(iChYV<NChY && CCmat(iChYU,iChYV)>q+myeps){
			iChYV++;
		}
		iChYV -= 1;
		if(iChYV != iChYV_last_grp){
			Ng++;
			iChYV_last_grp = iChYV;
		}
	}
	return Ng;
}
void calc_idxChY_and_Nc_vec(int *idxChY, int *Nc_vec, int Ng, SpEOMatrixD &CCmat, double q, int channelOffset){
	double myeps = 0.0001;
	int NChY = CCmat.rows();
	int iChYU, iChYV;
	int iChYV_last_grp = -1;
	int iG=-1;
	for(iChYU=0; iChYU<NChY; iChYU++){
		iChYV = iChYU;
		while(iChYV<NChY && CCmat(iChYU,iChYV)>q+myeps){
			iChYV++;
		}
		iChYV -= 1;
		if(iChYV != iChYV_last_grp){
			iG++;
			idxChY[iG] = iChYU + channelOffset;
			Nc_vec[iG] = iChYV-iChYU+1;
			iChYV_last_grp = iChYV;
		}
	}
}
void calc_refined_CCmins(double &CCmin1, double &CCmin2, double &CCmin3, SpEOMatrixD &CCmat_grp, double CCmin_orig, int Nc_max){
	int N = CCmat_grp.rows();
	int A = (N*N + N)/2; // number of entries in upper half incl. diagonal


	//###########################
	// CCmin1
	//###########################
	int N_sub = N-Nc_max;
	int B = (N_sub*N_sub + N_sub)/2; // number of entries in upper half of upper sub-matrix incl. "diagonal"
	int desired_num_entries = A-B;  // index corresponding to CCmin1 in sorted vec which contains the upper half of CCmin_orig

	double CCmat_grp_sum = 0;
	double CCmat_grp_vec[A];
	int a=-1;
	for(int u=0; u<N; u++){
		for(int v=u; v<N; v++){
			a++;
			CCmat_grp_vec[a] = CCmat_grp(u,v);
			CCmat_grp_sum += CCmat_grp_vec[a];
		}
	}
	// A==a-1 must be fulfilled after loop.
	if(A!=a+1){
		cerr << "A!=a+1 !!! -> Fatal Error!" << endl << endl << endl;
		exit(2);
	}
	// sort CCmat_grp_vec in a decreasing manner
//	int b_tilde[A];
	std::vector<argsort_pair> CCmat_grp_vec_for_sorting(A);
	for(int ii=0; ii<A; ii++){
		CCmat_grp_vec_for_sorting[ii].first = ii;
		CCmat_grp_vec_for_sorting[ii].second = CCmat_grp_vec[ii];
	}

	std::sort(CCmat_grp_vec_for_sorting.begin(), CCmat_grp_vec_for_sorting.end(), argsort_comp);
//	for(int ii=0; ii< A; ii++) {
//		b_tilde[CCmat_grp_vec_for_sorting[A-1-ii].first] = ii;
//	}
	// find a linear combination of the bands in ImX_LR to approximate ImY_g_avg while favoring the highly correlated bands
//	double c_hat[A];

//	for(int ii=0; ii<A; ii++){
//		cout << "CCmat_grp_vec_for_sorting[ii=="<<ii<<"] = ("
//				<<  CCmat_grp_vec_for_sorting[ii].first
//				<< ", "
//				<<  CCmat_grp_vec_for_sorting[ii].second
//				<< ")"
//				<< endl;
//	}
	CCmin1 = CCmat_grp_vec_for_sorting[A-desired_num_entries].second;
	//###########################
	// CCmin2
	//###########################
	int kk=0;
	double sub_sum = 0;
	while( sub_sum < (1./10)*CCmat_grp_sum && kk<A){
		sub_sum += CCmat_grp_vec_for_sorting[kk].second;
		kk++;
	}
	CCmin2 = CCmat_grp_vec_for_sorting[kk].second;

//	cout << "CCmat_grp_sum = " << CCmat_grp_sum << endl
//		 << "sub_sum = " << sub_sum << endl
//		 << "kk=" << kk << ", A=" << A << endl
//		 << "CCmin2=" << CCmin2 << endl;
	//###########################
	// CCmin3
	//###########################
	CCmin3 = CCmin_orig + (2./3)*(1.-CCmin_orig);
}

void grp_sngl_bnd_grps_with_high_corr_nbrs(int **Nc_vec_grp, int **idxChY_grp, int *Ng_grp, int iG_lrg, SpEOMatrixD &CCmat){
	// group each single channel group with its higher correlated neighbor channel
	// case 0: single channel group has no left neighbor (first group):
	if(Nc_vec_grp[iG_lrg][0]==1){
		Nc_vec_grp[iG_lrg][0] = 2;
	}
	// case 1: single channel group has no right neighbor (last group):
	if(Nc_vec_grp[iG_lrg][Ng_grp[iG_lrg]-1]==1){
		idxChY_grp[iG_lrg][Ng_grp[iG_lrg]-1] -= 1;
		Nc_vec_grp[iG_lrg][Ng_grp[iG_lrg]-1] =2;
	}
	// case 2: single channel groups has both left and right neighbor:
	for(int iG_sub=1; iG_sub<Ng_grp[iG_lrg]-1; iG_sub++){
		if(Nc_vec_grp[iG_lrg][iG_sub]==1){
//			cout << "  idxChY_grp[iG_lrg=" << iG_lrg << "][iG_sub=" << iG_sub << "] = " << idxChY_grp[iG_lrg][iG_sub] << endl
//			     << "    ## CCmat(idxChY_grp[iG_lrg][iG_sub],idxChY_grp[iG_lrg][iG_sub]-1) = " << CCmat(idxChY_grp[iG_lrg][iG_sub],idxChY_grp[iG_lrg][iG_sub]-1) << endl
//			     << "    ## CCmat(idxChY_grp[iG_lrg][iG_sub],idxChY_grp[iG_lrg][iG_sub]+1) = " << CCmat(idxChY_grp[iG_lrg][iG_sub],idxChY_grp[iG_lrg][iG_sub]+1) << endl;
			if(  CCmat(idxChY_grp[iG_lrg][iG_sub],idxChY_grp[iG_lrg][iG_sub]-1)
			   > CCmat(idxChY_grp[iG_lrg][iG_sub],idxChY_grp[iG_lrg][iG_sub]+1) ){
				idxChY_grp[iG_lrg][iG_sub] -= 1;
			}
			Nc_vec_grp[iG_lrg][iG_sub] = 2;
		}
	}
}

//##########################################
//#                NNLS
//##########################################
template <typename MatrixType>
bool testNNLS(const MatrixType &A,
              const Matrix<typename MatrixType::Scalar, MatrixType::RowsAtCompileTime, 1> &b,
              const Matrix<typename MatrixType::Scalar, MatrixType::ColsAtCompileTime, 1> &x,
              typename MatrixType::Scalar eps=1e-10)
//              double eps)
//              typename MatrixType::Scalar eps=1e-10)
{
  NNLS<MatrixType> nnls(A, 30, eps);
  if (! nnls.solve(b)) {
    std::cerr << __FILE__ << ": Convergence failed!" << std::endl; return false;
  }
  Array<typename MatrixType::Scalar, MatrixType::ColsAtCompileTime, 1> err
      = (x-nnls.x()).array().abs();
  if (err.maxCoeff() > 1e-6){
    std::cerr << __FILE__ << ": Precision error: Expected (" << x.transpose()
              << ") got: (" << nnls.x().transpose()
              << "), err: " << err.maxCoeff() << std::endl;
    return false;
  }

  if (! nnls.check(b)) {
    std::cerr << __FILE__ << ": check() KKT returned false!" << std::endl; return false;
  }

  return true;
}

// 4x2 problem, unconstrained solution positive
bool case_1 () {
  Matrix<double, 4, 2> A(4,2);
  Matrix<double, 4, 1> b(4);
  Matrix<double, 2, 1> x(2);
  A << 1, 1,  2, 4,  3, 9,  4, 16;
  b << 0.6, 2.2, 4.8, 8.4;
  x << 0.1, 0.5;

  return testNNLS(A, b, x);
}

// 4x3 problem, unconstrained solution positive
bool case_2 () {
  Matrix<double, 4, 3> A(4,3);
  Matrix<double, 4, 1> b(4);
  Matrix<double, 3, 1> x(3);

  A << 1,  1,  1,
       2,  4,  8,
       3,  9, 27,
       4, 16, 64;
  b << 0.73, 3.24, 8.31, 16.72;
  x << 0.1, 0.5, 0.13;

  return testNNLS(A, b, x);
}

// Simple 4x4 problem, unconstrained solution non-negative
bool case_3 () {
  Matrix<double, 4, 4> A(4, 4);
  Matrix<double, 4, 1> b(4);
  Matrix<double, 4, 1> x(4);

  A <<1,  1,  1,   1,
      2,  4,  8,  16,
      3,  9, 27,  81,
      4, 16, 64, 256;
  b << 0.73, 3.24, 8.31, 16.72;
  x << 0.1, 0.5, 0.13, 0;

  return testNNLS(A, b, x);
}

// Simple 4x3 problem, unconstrained solution non-negative
bool case_4 () {
  Matrix<double, 4, 3> A(4, 3);
  Matrix<double, 4, 1> b(4);
  Matrix<double, 3, 1> x(3);

  A <<1,  1,  1,
      2,  4,  8,
      3,  9, 27,
      4, 16, 64;
  b << 0.23, 1.24, 3.81, 8.72;
  x << 0.1, 0, 0.13;

  return testNNLS(A, b, x);
}

// Simple 4x3 problem, unconstrained solution indefinite
bool case_5 () {
  Matrix<double, 4, 3> A(4, 3);
  Matrix<double, 4, 1> b(4);
  Matrix<double, 3, 1> x(3);

  A <<1,  1,  1,
      2,  4,  8,
      3,  9, 27,
      4, 16, 64;
  b << 0.13, 0.84, 2.91, 7.12;
  // Solution obtained by original nnls() implementation
  // in Fortran
  x << 0.0, 0.0, 0.1106544;

  return testNNLS(A, b, x);
}

// 200x100 random problem
bool case_6 () {
 // MatrixXd A(200,100); A.setRandom();
 // VectorXd b(200); b.setRandom();
 // VectorXd x(100);

  SpEOMatrixD A(200,100); A.setRandom();
  SpEOVectorD b(200); b.setRandom();
  SpEOVectorD x(100);

  return NNLS<MatrixXd>::solve(A, b, x);
}


//int main(int argc, char *argv[])
int NNLS_testing_CG(){
  // Run test cases...
  bool ok = true;
  cout << "case_1.." << endl;
  ok &= case_1();
  cout << "case_2.." << endl;
  ok &= case_2();
  cout << "case_3.." << endl;
  ok &= case_3();
  cout << "case_4.." << endl;
  ok &= case_4();
  cout << "case_5.." << endl;
  ok &= case_5();
  cout << "case_6.." << endl;
  ok &= case_6();

  if (ok) return 0;
  else return -1;

}
void transform_SpEODataset_to_2D(SpEODataset *Im ,SpEOMatrixD &Im_2D){
  int iCh,u,v,uv;
  for(iCh=0; iCh<Im->get_NCh(); iCh++){
    for(u=0; u<Im->get_sizeU(); u++){
      for(v=0; v<Im->get_sizeV(); v++){
        uv = u*Im->get_sizeV() + v;
        Im_2D(iCh,uv) = (double)(Im->get_rasterBands()[iCh]->get_bandDataMat()->coeff(u,v));
      }
    }
  }
}
bool estimate_SRFs_and_ImX_shift(SpEOMatrixD &SRF, double *&ImX_shift, SpEOMatrixD &ImY_2D, SpEOMatrixD &ImX_LR_2D,int my_rank){
    
    bool ok = true;
    double my_scaling_factor = ImX_LR_2D.mean();
    ImY_2D /= my_scaling_factor;
    ImX_LR_2D /= my_scaling_factor;
 
    // tmp:
//    SpEOMatrixD SRFtmp(SRF.rows(),SRF.cols());
//    SpEOMatrixD SRFtmp_noUniqueness(SRF.rows(),SRF.cols());
//    SpEOMatrixD SRF_est_incl_shift(SRF.rows(),SRF.cols());

    int NChY = ImY_2D.rows();
    int NChX = ImX_LR_2D.rows();
    int NPix = ImY_2D.cols();
    // generate A
//    SpEOMatrixD A_plain = ImY_2D.transpose();
//    cout << "A_plain.bottomRightCorner(5,5) = " << endl << A_plain.bottomRightCorner(5,5) << endl << endl;
    SpEOMatrixD A(NPix,NChY+2);
    A.block(0,0,NPix,NChY) = ImY_2D.transpose();
    A.rightCols(2) = -SpEOMatrixD::Ones(NPix,2);
    A.rightCols(1) *= -1;
//    cout << "A.bottomRightCorner(5,5) = " << endl << A.bottomRightCorner(5,5) << endl << endl;
  
    // generate A0 which is equal to A but with an extra row to ensure uniqueness of the solution    
//    SpEOMatrixD A0(NPix+1,NChY+2);
//    A0.topRows(NPix) = A;
//    A0.bottomRows(1) = SpEOMatrixD::Zero(1,NChY+2);
//    A0(NPix,NChY+1) = 1.0;
//    cout << "A0.bottomRightCorner(5,5) = " << endl << A0.bottomRightCorner(5,5) << endl << endl;

    int iChX;
    for(iChX=0; iChX<NChX; iChX++){
      // generate b
//      SpEOVectorD b_plain = ImX_LR_2D.row(iChX);
      SpEOVectorD b = ImX_LR_2D.row(iChX);
//  //  generate b0 which is equal to A but with an extra row to ensure uniqueness of the solution
//      SpEOVectorD b0(NPix+1);
//      b0.head(NPix) = b;
//      b0(NPix) = 0.0;

      // generate x
      SpEOVectorD x(NChY+2);
//      SpEOVectorD x0(NChY+2);
//      SpEOVectorD x_plain(NChY);

      // least squares without ImX shift
//      cout << "dim(A_plain)= ["<< A_plain.rows() <<" x "<< A_plain.cols() <<"]" << endl;
//      cout << "dim(b_plain)= ["<< b_plain.rows() <<" x "<< b_plain.cols() <<"] and b_plain.tail(5) = " << b_plain.tail(5).transpose() << endl;
//      cout << "dim(x_plain)= ["<< x_plain.rows() <<" x "<< x_plain.cols() <<"]" << endl;
//      bool all_good0 = NNLS<MatrixXd>::solve(A_plain, b_plain, x_plain);
//      cout << "x_plain.tail(5) = " << x_plain.tail(5).transpose() << endl << endl << endl;
//      SRF.row(iChX) = x_plain;

      // with last row in A and b (uniqueness guaranteed)
//      cout << "dim(A0)= ["<< A0.rows() <<" x "<< A0.cols() <<"]" << endl;
//      cout << "dim(b0)= ["<< b0.rows() <<" x "<< b0.cols() <<"] and b0.tail(5) = " << b0.tail(5).transpose() << endl;
//      cout << "dim(x0)= ["<< x0.rows() <<" x "<< x0.cols() <<"]" << endl;
//      bool all_good1 = NNLS<MatrixXd>::solve(A0, b0, x0);
//      SRFtmp.row(iChX) = x0.head(NChY);
//      ImX_shift[iChX] = x0.coeff(NChY)-x0.coeff(NChY+1);
//      cout << "x0.tail(5) = " << x0.tail(5).transpose() << endl << endl << endl;

      // without last row in A and b (no uniqueness guaranteed)
      //SpEOMatrixD A1 = A0.topRows(NPix);
      //SpEOVectorD b1 = b;
     // cout << "dim(A)= ["<< A.rows() <<" x "<< A.cols() <<"]" << endl;
     //cout << "dim(b)= ["<< b.rows() <<" x "<< b.cols() <<"] and b.tail(5) = " << b.tail(5).transpose() << endl;
     // cout << "dim(x)= ["<< x.rows() <<" x "<< x.cols() <<"]" << endl;
      if(my_rank == 0) cout << " iChX=" << iChX << ": ";
      bool all_good = NNLS<MatrixXd>::solve(A, b, x);
      if(my_rank == 0) cout << ".. done! ";
      ok &= all_good;

      SRF.row(iChX) = x.head(NChY);
      ImX_shift[iChX] = x.coeff(NChY)-x.coeff(NChY+1);
      if(my_rank == 0) cout << " ImX_shift = " << ImX_shift[iChX]*my_scaling_factor << endl;
      //cout << "x.tail(5) = " << x.tail(5).transpose() << endl << endl << endl;


      //bool all_good = NNLS<MatrixXd>::solve(A, b, x);
   //   cout << "dim(A)= ["<< A.rows() <<" x "<< A.cols() <<"]" << endl;
   //   cout << "dim(b)= ["<< b.rows() <<" x "<< b.cols() <<"]" << endl;
   //   cout << "dim(x)= ["<< x.rows() <<" x "<< x.cols() <<"]" << endl;
 //     cout << "dim(AAA0)= ["<< AAA0.rows() <<" x "<< AAA0.cols() <<"]" << endl;
 //     cout << "dim(bbb0)= ["<< bbb0.rows() <<" x "<< bbb0.cols() <<"]" << endl;
 //     cout << "dim(x)= ["<< x.rows() <<" x "<< x.cols() <<"]" << endl;
 //     bool all_good = NNLS<MatrixXd>::solve(AAA0, bbb0, x);
 //   cout << endl << "bp6.1 all_good = " << all_good << endl << "x = " << x.transpose() << endl << endl;
 //   cout << endl << "x.head(NChY) = " << x.head(NChY) << endl << endl;

 //     cout << "dim(AAA)= ["<< AAA.rows() <<" x "<< AAA.cols() <<"]" << endl;
 //     cout << "dim(bbb)= ["<< bbb.rows() <<" x "<< bbb.cols() <<"]" << endl;
 //     cout << "dim(x)= ["<< x.rows() <<" x "<< x.cols() <<"]" << endl;
 //     bool all_good2 = NNLS<MatrixXd>::solve(AAA, bbb, x);

 //   cout << endl << "bp6.2 all_good2 = " << all_good2 << endl << "x = " << x.transpose() << endl << endl;
 //   cout << endl << "x.head(NChY) = " << x.head(NChY) << endl << endl;
 //     SRF.row(iChX) = x.head(NChY);
 //  // cout << endl << "bp7 " << endl << endl;
 //  // cout << endl << "bp8 " << endl << endl;
    }
    
    //SpEOVectroD SRF_band(ImY_2D)
   

//    cout << endl << "SRF (estimated while ImX shift was accounted for) = " << endl << SRF << endl << endl << endl;
//    cout << endl << "SRFtmp (estimated) = " << endl << SRFtmp << endl << endl << endl << endl;
//    cout << endl << "SRFtmp_noUniqueness (estimated) = " << endl << SRFtmp_noUniqueness << endl << endl << endl << endl;
//    cout << endl << "SRF_est_incl_shift (estimated) = " << endl << SRF_est_incl_shift << endl << endl << endl << endl;

    ImY_2D *= my_scaling_factor;
    ImX_LR_2D *= my_scaling_factor;

    for(iChX=0; iChX<NChX; iChX++){
      ImX_shift[iChX] *= my_scaling_factor;
   //   cout << "ImX_shift[iChX="<< iChX << "] = " << ImX_shift[iChX]<< ", while mean in this band is = " << ImX_LR_2D.row(iChX).mean() << endl;
    }

//    SpEOMatrixD AA(200,100); AA.setRandom();
//    SpEOVectorD bb(200); bb.setRandom();
//    SpEOVectorD xx(100);   
//    return NNLS<MatrixXd>::solve(AA, bb, xx);

    return ok;
}


