/*
 * dataIO.cpp
 *
 *  Created on: Apr 25, 2013
 *      Author: claas Grohnfeldt
 */

#include "dataIO.h"
#include "filter.h"

using namespace std;
using namespace Eigen;

SpEODictionary::SpEODictionary(SpEODataset *corrDS,
		SpEOFusionSetting *fSetting) {
	cout << "Create dictionary corresponding to dataset:" << corrDS->get_fname()
			<< std::endl;
	this->correspDataset = corrDS;
}

void printOneAsTwoDimArray(int nVSz, int nUSz, float *oneDimArray, int numRows,
		int numCols, int *rows, int *cols) {
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


void SpEODataset::dataRead(string fnm, SpEOReport *report) {
	int my_rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
	this->set_fname(fnm);
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
	}

	GDALDataset *poDataset;
	GDALAllRegister();
	poDataset = (GDALDataset *) GDALOpen(this->get_fname(), GA_ReadOnly);

	this->sizeU = poDataset->GetRasterYSize();
	this->sizeV = poDataset->GetRasterXSize();
	this->NCh = poDataset->GetRasterCount();
	this->imFormat = poDataset->GetDriver()->GetDescription();
	this->imFormatLong = poDataset->GetDriver()->GetMetadataItem(
			GDAL_DMD_LONGNAME);
	this->projectionRef = poDataset->GetProjectionRef();

	if (my_rank == 0) {
		report->file << " - Information about the (GDAL-)dataset:\n"
				<< "     - Driver:\t" << this->imFormat << "/"
				<< this->imFormatLong << "\n" << "     - Image Size:\t"
				<< this->sizeV << "(VSize) x " << this->sizeU << "(USize) x "
				<< this->NCh << "(RasterCount)\n";
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
	this->rasterBands = new SpEORasterBand*[this->NCh];
	int bGotMin, bGotMax;

	for (int iCh = 1; iCh < this->NCh + 1; iCh++) {
		this->rasterBands[iCh - 1] = new SpEORasterBand();

		// ====================================
		// fetching a raster band
		// ====================================
		if (my_rank == 0) {
			report->file << "     - Raster band no. " << iCh << ":  ";
		}

		GDALRasterBand *poBand;
		poBand = poDataset->GetRasterBand(iCh);
		

		this->rasterBands[iCh - 1]->bandDataType =
				static_cast<SpEODataType>((int) (poBand->GetRasterDataType()));

		this->rasterBands[iCh - 1]->nBand = iCh;
		this->rasterBands[iCh - 1]->minMaxVal[0] = poBand->GetMinimum(&bGotMin);
		this->rasterBands[iCh - 1]->minMaxVal[1] = poBand->GetMaximum(&bGotMax);
		if (!(bGotMin && bGotMax)) {
			GDALComputeRasterMinMax(poBand, TRUE,	this->rasterBands[iCh - 1]->minMaxVal);
		}
		this->rasterBands[iCh - 1]->dataType = GDALGetDataTypeName(
				poBand->GetRasterDataType());
		this->rasterBands[iCh - 1]->colorInterp =
				GDALGetColorInterpretationName(
						poBand->GetColorInterpretation());
		if (my_rank == 0) {
			report->file << "Type = " << this->rasterBands[iCh - 1]->dataType
					<< ",  ";
			report->file << "ColorInterp = "
					<< this->rasterBands[iCh - 1]->colorInterp << ",  ";
			report->file << "(Min, Max) = ("
					<< this->rasterBands[iCh - 1]->minMaxVal[0] << ", "
					<< this->rasterBands[iCh - 1]->minMaxVal[1] << ")\n";
			if (poBand->GetOverviewCount() > 0)
				report->file << "Band has " << poBand->GetOverviewCount()
						<< " overviews.\n";
			if (poBand->GetColorTable() != NULL)
				report->file << "Band has a color table with "
						<< poBand->GetColorTable()->GetColorEntryCount()
						<< "entries.\n";
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
		if(this->get_imFlag()==imFlag_Z || this->get_imFlag()==imFlag_Z_init){
			double *bandData = (double *) CPLMalloc( sizeof(double) * this->sizeV * this->sizeU);
			poBand->RasterIO(GF_Read, 0, 0, this->sizeV, this->sizeU, bandData, this->sizeV, this->sizeU, GDT_Float64, 0, 0);
			this->rasterBands[iCh - 1]->bandDataMatD = SpEOMatrixD::Zero(this->sizeU, this->sizeV);
			this->rasterBands[iCh - 1]->bandDataMatD = SpEOMatrixD::Map(bandData, this->sizeU, this->sizeV);
			CPLFree(bandData);
			bandData = NULL;
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
		}else{
			float *bandData = (float *) CPLMalloc( sizeof(float) * this->sizeV * this->sizeU);
			poBand->RasterIO(GF_Read, 0, 0, this->sizeV, this->sizeU, bandData, this->sizeV, this->sizeU, GDT_Float32, 0, 0);
			this->rasterBands[iCh - 1]->bandDataMat = SpEOMatrixF::Zero(this->sizeU, this->sizeV);
			this->rasterBands[iCh - 1]->bandDataMat = SpEOMatrixF::Map(bandData, this->sizeU, this->sizeV);
			CPLFree(bandData);
			bandData = NULL;


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
	}
	GDALClose(poDataset);
	if (my_rank == 0) {
		report->file << "\n";
		report->file.close();
		cout << "done!" << endl;
	}
}


void SpEODataset::dataWrite(SpEOReport *report, SpEODataIOSetting *dSet, SpEOFusionSetting *fSet, SpEOParallelSetting *pSet, SpEOGlobalParams *glPrms, SpEOPaths *paths, MPI_Comm comm_write, SpEODataset *ImZ) {

	int my_rank; int my_processes;
	MPI_Comm_rank(comm_write, &my_rank);
	MPI_Comm_size(comm_write, &my_processes);

	int u,v,iL, uL, vL, uH, vH, nL, fDS, iChZ;
	fDS = glPrms->fDS;
	nL = glPrms->sizeUL*glPrms->sizeVL;
	SpEOMatrixD *area_sum = new SpEOMatrixD[this->NCh];

	if(my_rank==0){
		cout << "\n"
				<< "###########################################################" << "\n"
				<< "##    Write image ImZ (from memory) in parallel using..  ##" << "\n"
				<< "##        - MPI I/O                                      ##" << "\n"
				<< "##        - ENVI header file                             ##" << "\n"
				<< "##        - band interleaved by pixel (BIP)              ##" << "\n"
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
				<< "###########################################################"
				<< "\n" << "\n";
	}
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
		buftype = MPI_DOUBLE;
		// communicator group comm_write opens the data file for writing only (and creating, if necessary) */
		result = MPI_File_open(comm_write, (char*)paths->fname_ImZ_out.c_str(), MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &fh);
		if(result != MPI_SUCCESS){
			sample_error(result, (char*)"MPI_File_open");
		}
		disp=0;
		result = MPI_File_set_view(fh, disp, etype, ftype, (char*)"native", MPI_INFO_NULL);
		if(result != MPI_SUCCESS){
			sample_error(result, (char*)"MPI_File_set_view");
		}
		iL = dSet->uLFirst*glPrms->sizeVL+dSet->vLFirst + my_rank;
		nL = dSet->uLLast*glPrms->sizeVL+dSet->vLLast +1;
		while(iL<nL){ // as long as there are LR pixels to work on
			// calculate the coordinates of the patch
			uL = iL / glPrms->sizeVL;
			vL = iL % glPrms->sizeVL;
			uH = uL*fDS;
			vH = vL*fDS;
			int uH_sub = uH - dSet->uLFirst*fDS;
			int vH_sub = vH - dSet->vLFirst*fDS;
			if(uL>=dSet->uLFirst && vL>=dSet->vLFirst && uL<=dSet->uLLast && vL<=dSet->vLLast){
				// initialization
				for(iChZ=0; iChZ<this->NCh; iChZ++){
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
		buftype = MPI_UNSIGNED_SHORT;

		// communicator group comm_write opens the data file for writing only (and creating, if necessary) */
		result = MPI_File_open(comm_write, (char*)paths->fname_ImZ_out.c_str(), MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &fh);
		if(result != MPI_SUCCESS){
			sample_error(result, (char*)"MPI_File_open");
		}
		disp=0;
		result = MPI_File_set_view(fh, disp, etype, ftype, (char*)"native", MPI_INFO_NULL);
		if(result != MPI_SUCCESS){
			sample_error(result, (char*)"MPI_File_set_view");
		}
		iL = dSet->uLFirst*glPrms->sizeVL+dSet->vLFirst + my_rank;
		nL = dSet->uLLast*glPrms->sizeVL+dSet->vLLast +1;
		while(iL<nL){ // as long as there are LR pixels to work on
			// calculate the coordinates of the patch
			uL = iL / glPrms->sizeVL;
			vL = iL % glPrms->sizeVL;
			uH = uL*fDS;
			vH = vL*fDS;
			int uH_sub = uH - dSet->uLFirst*fDS;
			int vH_sub = vH - dSet->vLFirst*fDS;
			if(uL>=dSet->uLFirst && vL>=dSet->vLFirst && uL<=dSet->uLLast && vL<=dSet->vLLast){
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


	this->set_fname(paths->fname_ImZ_out.c_str());

	string fileName = paths->fname_ImZ_out + ".hdr";
	ofstream headerfile;
	headerfile.open(fileName.c_str(),	fstream::in | fstream::out | fstream::app);

	headerfile
	<< "ENVI"
	<< "\n" << "description = {"
	<< "\n" << paths->fname_ImZ_out
	<< "\n" << "}"
	<< "\n" << "samples = " << (dSet->vLLast-dSet->vLFirst+1)*glPrms->fDS
	<< "\n" << "lines   = " << (dSet->uLLast-dSet->uLFirst+1)*glPrms->fDS
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
	<< "\n" << "band names = {";
	for(int iCh=0; iCh<this->NCh-1; iCh++){
		headerfile << "\n" << "Band " << iCh+1+dSet->chBundleFirst << "," ;
	}
	headerfile << "\n" << "Band " << this->NCh+dSet->chBundleFirst << "}";
	headerfile << "\n";


	report->file.close();
	headerfile.close();
	chmod(fileName.c_str(), S_IRWXU | S_IRWXG | S_IRWXO);
}


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
	result = MPI_File_open(comm_write, (char*)paths->fname_ImZ_out.c_str(), MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &fh);
	if(result != MPI_SUCCESS){
		sample_error(result, (char*)"MPI_File_open");
	}
	disp=0;
	result = MPI_File_set_view(fh, disp, etype, ftype, (char*)"native", MPI_INFO_NULL);
	if(result != MPI_SUCCESS){
		sample_error(result, (char*)"MPI_File_set_view");
	}
	iL = dSet->uLFirst*glPrms->sizeVL+dSet->vLFirst + my_rank;
	nL = dSet->uLLast*glPrms->sizeVL+dSet->vLLast +1;

	while(iL<nL){ // as long as there are LR pixels to work on

		// calculate the coordinates of the patch
		uL = iL / glPrms->sizeVL;
		vL = iL % glPrms->sizeVL;
		uH = uL*fDS;
		vH = vL*fDS;

		if(uL>=dSet->uLFirst && vL>=dSet->vLFirst && uL<=dSet->uLLast && vL<=dSet->vLLast){

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
	char buf [paths->dir_tmp_patches.length()+42];
	sprintf (buf, "%s/%06d/patch_u%04d_v%04d_iP%06d.csv", paths->dir_tmp_patches.c_str(), iP, uP, vP, iP);
	int stat_CSV_read = read_CSV(&myReadMat, buf, ',', 0);
	/*
	for(int i=0; stat_CSV_read==-1 && i<dSet->dir_tmp_patches_additional_num; i++){
		char buf_tmp [paths->dir_tmp_patches_additional[i].length()+50];
		sprintf (buf_tmp, "%s/%06d/patch_u%04d_v%04d_iP%06d.csv", paths->dir_tmp_patches_additional[i].c_str(), iP, uP, vP, iP);
		stat_CSV_read = read_CSV(&myReadMat, buf_tmp, ',', 0);
		cout << "try to open tmp patch file: " << buf_tmp << endl;
	}
	*/
	if(stat_CSV_read==-1){
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
	for (int iCh=0; iCh < this->NCh; iCh++) {
		delete this->rasterBands[iCh];
	}
	delete [] this->rasterBands;
	delete [] geoTransform;
}

SpEORasterBand::~SpEORasterBand(){
	this->bandDataMat  = SpEOMatrixF::Zero(0,0);
	this->bandDataMatD = SpEOMatrixD::Zero(0,0);
}

void setMetaInfo(SpEODataset *ImZ, SpEODataset *ImY, SpEODataset *ImX, SpEODataIOSetting *dSet, SpEOGlobalParams *glPrms) {
	ImZ->imFormat = ImY->imFormat;
	ImZ->imFormatLong = ImY->imFormatLong;
	ImZ->projectionRef = ImX->projectionRef;
	ImZ->sizeU = glPrms->sizeUH_red;
	ImZ->sizeV = glPrms->sizeVH_red;
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
	//glPrms.Nc_vec[0]       = glPrms.NChY;
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
}


void SpEODataset::copy_from_dataset( int firstBand, int lastBand, SpEODataset *sourceIm){

	this->set_fname(sourceIm->get_fname());
	this->imFlag        = sourceIm->get_imFlag();
	this->sizeU         = sourceIm->get_sizeU();
	this->sizeV         = sourceIm->get_sizeV();
	this->NCh           = lastBand - firstBand + 1;
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
			this->rasterBands[iCh]->bandDataMatD = SpEOMatrixD::Map(sourceIm->get_rasterBands()[firstBand+iCh]->get_bandDataMatD()->data(),
					                                                sourceIm->get_rasterBands()[firstBand+iCh]->get_bandDataMatD()->rows(),
					                                                sourceIm->get_rasterBands()[firstBand+iCh]->get_bandDataMatD()->cols());
		}else{
			this->rasterBands[iCh]->bandDataMat = SpEOMatrixF::Map(sourceIm->get_rasterBands()[firstBand+iCh]->get_bandDataMat()->data(),
													               sourceIm->get_rasterBands()[firstBand+iCh]->get_bandDataMat()->rows(),
													               sourceIm->get_rasterBands()[firstBand+iCh]->get_bandDataMat()->cols());

		}
	}
}

void SpEODataset::copyMetaInfoFromDatasets(SpEODataset *sourceImSpatial, SpEODataset *sourceImSpectral, SpEODataset *sourceImFormat, SpEODataFormat dataFormat){

	this->set_fname(sourceImSpatial->get_fname());
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
			this->rasterBands[iCh]->bandDataMatD = SpEOMatrixD::Zero(sourceImSpatial->get_rasterBands()[iChSoureImSpatial]->get_bandDataMatD()->rows(),
					                                                 sourceImSpatial->get_rasterBands()[iChSoureImSpatial]->get_bandDataMatD()->cols());
		}else{
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
			this->rasterBands[iCh]->bandDataMatD = SpEOMatrixD::Map(sourceIm->get_rasterBands()[iCh-firstBand]->get_bandDataMatD()->data(),
					                                                sourceIm->get_rasterBands()[iCh-firstBand]->get_bandDataMatD()->rows(),
					                                                sourceIm->get_rasterBands()[iCh-firstBand]->get_bandDataMatD()->cols());
		}else{
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
			cout << "chBundleFirst has been corrected from " << dSet->chBundleFirst << " to " << ImY->NCh-1 << "=NChY-1," << endl;
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
			cout << "chBundleLast has been corrected from " << dSet->chBundleLast << " to " << ImY->NCh-1 << "=NChY-1." << endl;
		}
		dSet->chBundleLast = ImY->NCh-1;
	}

	short NCh_buf = ImY->NCh;
	if(dSet->chBundleLast-dSet->chBundleFirst+1 < ImY->NCh){
		if(my_rank==0) {
			cout << "NChY has been corrected from " << ImY->NCh << " to " << dSet->chBundleLast-dSet->chBundleFirst+1 << "=chBundleLast-chBundleFirst+1." << endl;
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
				SpEOMatrixF tmp = Im->get_rasterBands()[iCh]->bandDataMat.block(myUFirst, myVFirst, mySizeU_red, mySizeV_red);
				Im->get_rasterBands()[iCh]->bandDataMat = tmp;
		}
		Im->sizeU = mySizeU_red;
		Im->sizeV = mySizeV_red;
	}
}


string get_current_time() {
	char curtime[50];
	time_t rawtime;
	struct tm * timeinfo;
	time(&rawtime);
	timeinfo = localtime(&rawtime);
	strftime(curtime, 50, "%y%m%d_%H%M%S", timeinfo);
	string curtimeStr(curtime);
	return curtimeStr;
}

void SpEOReport::initialize(SpEOPaths *paths,
			    SpEODataIOSetting *dSetting,
                            SpEOFusionSetting *fSetting,
                            SpEOOutputSetting *oSetting,
                            SpEOSolverSetting *sSetting,
                            SpEOParallelSetting *pSetting,
                            int argc, char **argv) {

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
	paths->dir_out += "/" + mydate + "_" + mytime + "_" + dSetting->jobName;
	mkdir(paths->dir_out.c_str(), 0777);
	chmod(paths->dir_out.c_str(), S_IRWXU | S_IRWXG | S_IRWXO);

	// paths->dir_out += "/" + dSetting->jobName;// + "_" + dSetting->jobID;
	// mkdir(paths->dir_out.c_str(), 0777);
	// chmod(paths->dir_out.c_str(), S_IRWXU | S_IRWXG | S_IRWXO);

	stringstream dir_out_tmp;
	dir_out_tmp << paths->dir_out; //+ "/"         
/*			<< mytime
//                      << "_NP" << glPrms->NP 
//                      << "_pTot" << my_processes 
//			<< "_ID"    << paths->dataSetID
//                      << "_ID"    << paths->dataSetID_str  
//                      << "_psz"   << fSetting->patchsize 
//			<< "_ovlp"  << fSetting->overlap 
//                      << "_NDP" << fSetting->NDP 
//                      << "_lmd" << fSetting->lambda
//                      << "_alg"  << fSetting->fMethod  
//                      << "_dS"   << fSetting->dictselect 
//                      << "_2Step" << fSetting->two_step_estimation 
//                      << "_Nc"    << fSetting->Nc 
//                      << "_No" << fSetting->No
//                      << "_alg"  << fSetting->fMethod  
//                      << "_dS"   << fSetting->dictselect 
//			<< "_itr"   << fSetting->iterMain 
//                      << "_NwCf"  << fSetting->useNewMethodForCalculatingZ 
//                      << "_ImXSm" << fSetting->useSimulatedImXforDictLearn
//			<< "_lX"    << fSetting->lambdaX_ABC
//			<< "_lY"    << fSetting->lambdaY_ABC
//			<< "_lZ"    << fSetting->lambdaZ_ABC 
//                      << "_lZ0In1stIter" << fSetting->lambdaZ_ABC_in_1st_iter 
//                      << "_fllImOpt"<< fSetting->LQ_post_opt_im
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
//			<< "_SpecNorm"  << fSetting->matrixNorm 
//			<< "_pixwiseMean"  << fSetting->addMeanPixelwise; 
*/
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
			<< " - file name ImX:             " << paths->fname_ImX << "\n"
			<< " - file name ImY:             " << paths->fname_ImY << "\n"
			<< " - file name ImZ (reference): " << paths->fname_ImZ_ref << "\n"
			<< " - file name ImZ_init:        " << paths->fname_ImZ_init << "\n"
			<< " - directory output data:     " << paths->dir_out << "\n"
			<< " - file name ImZ:             " << paths->fname_ImZ_out << "\n\n"
	        << " - use_estimated_SRFs:        " << fSetting->use_estimated_SRFs << "\n"
			<< " - file name SRF:             " << paths->fname_SRF << "\n";
/*			
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
*/
	/*
	if(dSetting->contUnfinishedRec){
		this->file  << " - file name CSV containing iPs of incomplete set of patches from previous  unfinished reconstruction: " << "\n"
			<< "                                  " << paths->PathToIncompletePatchSetCSV << "\n"
			<< " - directories that may containe previously reconstructed TMP patches: " << "\n";
		for(int i=0; i<dSetting->dir_tmp_patches_additional_num; i++){
			this->file << "              .. directory no. " << i+1 << ": "<< paths->dir_tmp_patches_additional[i] << "\n";
		}
	}*/
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
			<< " - vLLast:  last horizontal low resolution pixel to be processed (before possible correction):  " << dSetting->vLLast << "\n";
			//<< " - delete_tmp_patch_folders:                                                                    " << dSetting->delete_tmp_patch_folders << "\n"
			//<< " - dir_tmp_patches_additional_num: number of additional directories containing relevant tmp patches: " << dSetting->dir_tmp_patches_additional_num << "\n"
			//<< " - continue unfinished reconstruction?:                                                         " << dSetting->contUnfinishedRec << "\n"
	this->file << "\n"
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
		this->file << "J-P-FISTA" << "\n";
		break;
	}
	default: {
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
 //this->file << " - Least Sq. post minimiz. for Z via CGLS:  " << fSetting->LQ_post_opt << "\n"
 this->file << "     -> Patch-wise coefficient estimation parameters:\n"
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
            << "                       - subspace_dim:                  " << fSetting->subspace_dim << "\n";
   //         << "     -> Patch-wise CGLS settings (former 'Eq.3'): ACTIVE=" << fSetting->LQ_post_opt << "\n"
			//<< "                               - lambda_X:        " << fSetting->lambdaX << "\n"
			//<< "                               - lambda_Y:        " << fSetting->lambdaY << "\n"
			//<< "                               - maxiter CGLS:    " << sSetting->maxiter_CGLS << "\n"
			//<< "                               - tol r CGLS:      " << sSetting->tol_r_CGLS << "\n"
			//<< "                               - Alphas fixed:    " << sSetting->fix_Alpha << "\n"
			//<< "                               - delta_m fixed:   " << sSetting->fix_delta_m << "\n";

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
			<< " - theta (minimum cross-correlation within spectral goups (double)): "
			<< fSetting->theta << "\n"
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

	this->file.close();


	//  Write to terminal:
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
				<< " - file name ImX:             " << paths->fname_ImX << "\n"
				<< " - file name ImY:             " << paths->fname_ImY << "\n"
				<< " - file name ImZ (reference): " << paths->fname_ImZ_ref << "\n"
				<< " - file name ImZ_init:        " << paths->fname_ImZ_init << "\n"
				<< " - directory output data:     " << paths->dir_out << "\n"
				<< " - file name ImZ:             " << paths->fname_ImZ_out << "\n"
				<< " - use_estimated_SRFs:                    " << fSetting->use_estimated_SRFs << "\n"
				<< " - file name SRF:                         " << paths->fname_SRF << "\n";
/*				<< " - ImZ_init_type = " << fSetting->ImZ_init_type << "(";
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
*/
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
				<< " - theta (minimum cross-correlation within spectral goups (double)):   "
				<< fSetting->theta
				<< "\n"
				<< " - tol_SRF (tolerance for decision matrix): "
				<< fSetting->tol_SRF << "\n"
				<< " - evaluate:                                " << fSetting->evaluate  << "\n"
				<< " - evaluate_ImZ_init:                       " << fSetting->evaluate_ImZ_init  << "\n"
				<< " - dictionary selection method ID: "
				<< fSetting->dictselect << "\n" << "\n"
				<< "     -> Patch-wise coefficient estimation parameters:\n"
				<< "                       - lambda_X_ABC:  " << fSetting->lambdaX_ABC << "\n"
				<< "                       - lambda_Y_ABC:  " << fSetting->lambdaY_ABC << "\n"
				<< "                       - lambda_Z_ABC:  " << fSetting->lambdaZ_ABC << "\n"
				<< "        - lambdaZ_ABC_in_1st_iter:      " << fSetting->lambdaZ_ABC_in_1st_iter << "\n"
				<< "     -> full image CGLS settings: ACTIVE=" << fSetting->LQ_post_opt_im << "\n"
				<< "                       - lambda_X_im:      " << fSetting->lambdaX_im << "\n"
				<< "                       - lambda_Y_im:      " << fSetting->lambdaY_im << "\n"
				<< "                       - maxiter CGLS_im:  " << sSetting->maxiter_CGLS_im << "\n"
				<< "                       - tol r CGLS_im:    " << sSetting->tol_r_CGLS_im << "\n";
				//<< "     -> Patch-wise CGLS settings (former 'Eq.3'): ACTIVE=" << fSetting->LQ_post_opt << "\n"
				//<< "                               - lambda_X:        " << fSetting->lambdaX << "\n"
				//<< "                               - lambda_Y:        " << fSetting->lambdaY << "\n"
				//<< "                               - maxiter CGLS:    " << sSetting->maxiter_CGLS << "\n"
				//<< "                               - tol r CGLS:      " << sSetting->tol_r_CGLS << "\n"
				//<< "                               - Alphas fixed:    " << sSetting->fix_Alpha << "\n"
				//<< "                               - delta_m fixed:   " << sSetting->fix_delta_m << "\n"
		cout	<< "*========================================*" << "\n"
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

	// double totalTime = MPI_Wtime() - this->curTimeSec;
	// //glPrms->timeTotal = totalTime;
	// this->curTime = get_current_time();
	this->file.open(this->fileName.c_str(),
			fstream::in | fstream::out | fstream::app);

	// string datePretty = "20" + this->curTime.substr(0, 2) + "-"
	// 		+ this->curTime.substr(2, 2) + "-" + this->curTime.substr(4, 2);
	// string timePretty = this->curTime.substr(7, 2) + ":"
	// 		+ this->curTime.substr(9, 2) + ":" + this->curTime.substr(11, 2);

	// this->file << "\n"
			// << "_______________________________________________________________________________\n\n"
			// << "ending date: " << datePretty << " " << timePretty
			// << "\n\n"
			// << "  -> Method needed                           " << totalTime << " seconds in total which of \n"
			// << "         - the main loop took                " << glPrms->timeMainLoop << " seconds,\n"
			// << "         - total dictionary selection took   " << glPrms->timeDictSelect << " seconds,\n"
			// << "         - average dictionary selection took " << glPrms->timeDictSelect_avg << " seconds,\n"
			// << "         - the data-to-file writing took " << glPrms->timeFileWrite << " seconds."
			// << "\n\n";

	this->file
			<< "###############################################################################\n"
			<< "###                               __         _                              ###\n"
			<< "###                              |__  |\\ |  | \\                             ###\n"
			<< "###                              |__  | \\|  |_/                             ###\n"
			<< "###                                                                         ###\n"
			<< "###############################################################################\n";

	this->file.close();
	chmod(this->fileName.c_str(), S_IRWXU | S_IRWXG | S_IRWXO);


	// cout    << "\n"
			// << "_______________________________________________________________________________\n\n"
			// << "ending date: " << datePretty << " " << timePretty
			// << "\n\n"
			// << "  -> Method needed                           " << totalTime << " seconds in total which of \n"
			// << "         - the main loop took                " << glPrms->timeMainLoop << " seconds,\n"
			// << "         - total dictionary selection took   " << glPrms->timeDictSelect << " seconds,\n"
			// << "         - average dictionary selection took " << glPrms->timeDictSelect_avg << " seconds,\n"
			// << "         - the data-to-file writing took " << glPrms->timeFileWrite << " seconds."
	cout	<< "\n"
			<< "\n"
			<< "###############################################################################\n"
			<< "###                               __         _                              ###\n"
			<< "###                              |__  |\\ |  | \\                             ###\n"
			<< "###                              |__  | \\|  |_/                             ###\n"
			<< "###                                                                         ###\n"
			<< "###############################################################################\n\n";

}

void SpEOReport::addEvaluation(SpEOAssessmentMetrics *assessm_res_HR,
								SpEODataIOSetting *dSetting,
								SpEOFusionSetting *fSetting,
                               SpEOOutputSetting *oSetting,
                               SpEOGlobalParams *glPrms,
                               int type) {
	SpEOAssessmentMetricsStr *assessm_res_HR_Str = new SpEOAssessmentMetricsStr;
	SpEOAssessmentMetricsStr_Sep *assessm_res_HR_StrSep = new SpEOAssessmentMetricsStr_Sep[glPrms->NChY];

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
	createTable(assessm_res_HR_Str, assessm_res_HR_StrSep, dSetting, fSetting, glPrms,
			maxDigits, type, this);
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
			<< " - theta (minimum cross-correlation within spectral goups (double)): "
			<< fSetting->theta << "\n"
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
			<< glPrms->NPV_sub << "\n";
	/*	
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
	*/
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
		m_data.push_back(cell);
	}
}
// read and save .cvs file content to DOUBLE matrix
int read_CSV(SpEOMatrixD *Mat, const char *fname_CSV, char delimiter, int skipLns) {
	int my_rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
	std::ifstream CSV_file(fname_CSV);
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
	}
	return -1;
}
// read and save .cvs file content to FLOAT matrix
int read_CSV(SpEOMatrixF *Mat, const char *fname_CSV, char delimiter, int skipLns) {
	int my_rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
	std::ifstream CSV_file(fname_CSV);
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
	std::ofstream CSV_file;
	CSV_file.open(fname_CSV, ios::out | ios::ate | ios::app);
	if (CSV_file.is_open()) {
		IOFormat CSVFmt(FullPrecision, 0, ",");
		CSV_file << Mat->format(CSVFmt);
		CSV_file << "\n";
	}else{
		cout << endl << "WARNING: Matrix could not be written to .csv file!" << endl;
	}
	CSV_file.close();
}


void write_Mat_to_CSV(SpEOMatrixD *Mat, const char *fname_CSV) {
	std::ofstream CSV_file;
	CSV_file.open(fname_CSV, ios::out | ios::ate | ios::app);
	if (CSV_file.is_open()) {
		IOFormat CSVFmt(FullPrecision, 0, ",");
		CSV_file << Mat->format(CSVFmt);
		CSV_file << "\n";
	}else{
		cout << endl << "WARNING: CSV_file '" << fname_CSV << "' could not be written!" << endl;
	}
	CSV_file.close();
	chmod(fname_CSV, S_IRWXU | S_IRWXG | S_IRWXO);
}


void read_SRF(SpEOGlobalParams *glPrms, SpEOMatrixD *Mat, string fname_CSV, char delimiter, bool normalize_SRF) {
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
	}
	// check sum to one constraint
	if(my_rank == 0){
		cout << "Check sum-to-one constraint. I.e. check if the SRFs are individually normalized..." << endl;
		cout << "The sum of the SRF entries in channel iChX is differs from 1.0 (ideal case) as follows (should be equal to 0.0)!" << endl;
	}
	for(int iChX=0; iChX<Mat->rows(); iChX++){
		double sum = 0.0;
		for(int iChY=0; iChY<Mat->cols(); iChY++){
			sum += Mat->coeff(iChX,iChY);
		}
		if(my_rank == 0){
			cout << "iChX=" << iChX << "  ->  sum(SRF("<< iChX <<"))-1 = " << sum-1 << endl;
		}
	}
	// normalize each SRF (i.e. each row) so that the entries in each row sum up to one
	if (normalize_SRF) {
		if(my_rank == 0){
			cout << "normalize SRF, so that the sum in every row (i.e. for every channel of the image X) is equal to one..." << endl;
		}
		for(int iChX=0; iChX<Mat->rows(); iChX++){
			double sum = 0;
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

/*
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
				sum = sum + (double)(*SRF)(iChX, dSet->chBundleFirst + a*iG+iC);
			}
			glPrms->decMat_C(iChX, iG) = sum;
			sum = 0;
		}
		for (int iC = 0; iC < fSet->Nc; iC++) {
			sum = sum + (double)(*SRF)(iChX, dSet->chBundleFirst + glPrms->NChZ - fSet->Nc+iC);
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
	
	glPrms->numProbPerPatch = 0;
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
	}glPrms->numPatchGroups = my_processes / pSet->numProcPerPatch;

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
	
	if(my_rank % pSet->numProcPerPatch < glPrms->numProbPerPatch % pSet->numProcPerPatch){
		glPrms->myNumProbPerPatch = glPrms->numProbPerPatch / pSet->numProcPerPatch + 1;
	}
	else {
		glPrms->myNumProbPerPatch = glPrms->numProbPerPatch / pSet->numProcPerPatch;
	}


	//glPrms->numPatchGroups = my_processes / pSet->numProcPerPatch;
	glPrms->numPatchGroups = my_processes;// / pSet->numProcPerPatch;
	
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
	
	if (my_rank == 0) {
		cout << ".. decision matrix successfully calculated!" << endl; 
	}


	// ###################################
	// #  Calculate Nc_vec, P_lmd, etc.  #
	// ###################################
	int iG, iChZ, iC, ipp;
	// (in case of Spectral Grouping, all values in Nc_vec are identical and equal to the Nc value specified by the user in fSettings)
	glPrms->Nc_vec = new int[glPrms->numProbPerPatch];
	for(ipp=0; ipp<glPrms->numProbPerPatch; ipp++){
		glPrms->Nc_vec[ipp] = fSet->Nc;
	}
	glPrms->P_lmd_vecs = new SpEOVectorD[glPrms->numProbPerPatch];
	for(ipp=0; ipp<glPrms->numProbPerPatch; ipp++){
		glPrms->P_lmd_vecs[ipp] = SpEOVectorD::Ones(glPrms->Nc_vec[ipp])*glPrms->decMat_C(glPrms->myChX[ipp], glPrms->myBundle[ipp]);
	}
	glPrms->P_lmd_idx_bl = new int*[glPrms->numProbPerPatch];
	int col_idx=0, row_idx;
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
	for(ipp=glPrms->numProbPerPatch-glPrms->Nm2[glPrms->Ng-1]; ipp<glPrms->numProbPerPatch; ipp++){
		glPrms->P_lmd_idx_bl[ipp] = new int[2];
		row_idx = glPrms->NChZ - glPrms->Nc_vec[ipp];
		glPrms->P_lmd_idx_bl[ipp][0] = row_idx;
		glPrms->P_lmd_idx_bl[ipp][1] = col_idx;
		col_idx += glPrms->Nc_vec[ipp];
	}
	glPrms->P_lmd_idx_row = new SpEOMatrixI[glPrms->NChZ];
	// first the number of non-trivial entries in each row (iChZ=0,...,NChZ-1) needs to be counted
	SpEOVectorI entries_cnt     = SpEOVectorI::Zero(glPrms->NChZ);
	SpEOVectorI entries_divider = SpEOVectorI::Zero(glPrms->NChZ);
	bool incr_entries_divider = false;
	bool iG_done[glPrms->Ng];
	for(iG=0; iG<glPrms->Ng; iG++){
		iG_done[iG]=false;
	}
	for(ipp=0; ipp<glPrms->numProbPerPatch; ipp++){
		iG = glPrms->myBundle[ipp];
		if(!iG_done[iG]){
			incr_entries_divider = true;
			iG_done[iG] = true;
		}else{
			incr_entries_divider = false;
		}
		for(iC=0; iC<glPrms->Nc_vec[ipp]; iC++){
			iChZ = glPrms->P_lmd_idx_bl[ipp][0]+iC;
			if(incr_entries_divider){
				entries_divider(iChZ)++;
			}
			entries_cnt(iChZ)++;
		}
	}
	for(iChZ=0; iChZ<glPrms->NChZ; iChZ++){
		glPrms->P_lmd_idx_row[iChZ] = SpEOMatrixI::Zero(2,entries_cnt.coeff(iChZ));
	}
	SpEOVectorI entries_cnt_tmp = SpEOVectorI::Zero(glPrms->NChZ);
	for(ipp=0; ipp<glPrms->numProbPerPatch; ipp++){
		for(iC=0; iC<glPrms->Nc_vec[ipp]; iC++){
			iChZ = glPrms->P_lmd_idx_bl[ipp][0]+iC;
				// store block index
				glPrms->P_lmd_idx_row[iChZ](0,entries_cnt_tmp(iChZ)) = ipp;
				// store coefficient index relative to corresponding block origin
				glPrms->P_lmd_idx_row[iChZ](1,entries_cnt_tmp(iChZ)) = iC;
			// calculate value in P_lmd (here implemented is only simple averaging!! If desired, that can be changed later)
			glPrms->P_lmd_vecs[ipp](iC) /= entries_divider(iChZ);
			entries_cnt_tmp(iChZ)++;
		}
	}
	delete[] startProc;
}
*/

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

	//##################################################################
	// write to .m file
	string fname_fSetting = dir_fSetting + "/" + "fusionSettings.m";
	std::ofstream fSetting_file_matlab(fname_fSetting.c_str());
	if (fSetting_file_matlab.is_open()) {
		fSetting_file_matlab << setiosflags(ios::fixed);
		fSetting_file_matlab << "  fSet.nrmlIm="<< fSetting->nrmlIm << ";" << endl;
		fSetting_file_matlab << "  fSet.nrmlDicts="<< fSetting->nrmlDicts << ";" << endl;
		fSetting_file_matlab << "  fSet.subtrMean="<< fSetting->substrMean << ";" << endl;
		fSetting_file_matlab << "  fSet.ImZ_ref_avlbl="<< fSetting->ImZ_ref_avlbl << ";" << endl;
		fSetting_file_matlab << "  fSet.Nc="<< fSetting->Nc << ";" << endl;
		fSetting_file_matlab << "  fSet.No="<< fSetting->No << ";" << endl;
		fSetting_file_matlab << "  fSet.tol_SRF="<< fSetting->tol_SRF << ";" << endl;
		fSetting_file_matlab << "  fSet.Nc_max="<< fSetting->Nc_max << ";" << endl;
		fSetting_file_matlab << "  fSet.theta="<< fSetting->theta << ";" << endl;
		fSetting_file_matlab << "  fSet.patchsize="<< fSetting->patchsize << ";" << endl;
		fSetting_file_matlab << "  fSet.overlap="<< fSetting->overlap << ";" << endl;
		fSetting_file_matlab << "  fSet.lambda="<< fSetting->lambda << ";" << endl;
		fSetting_file_matlab << "  fSet.NDP="<< fSetting->NDP << ";" << endl;
		fSetting_file_matlab << "  fSet.two_step_estimation="<< fSetting->two_step_estimation << ";" << endl;
		fSetting_file_matlab << "  fSet.dictselect="<< fSetting->dictselect << ";" << endl;
		fSetting_file_matlab << "  fSet.matrixNorm="<< fSetting->matrixNorm << ";" << endl;
		fSetting_file_matlab << "  fSet.addMeanPixelwise="<< fSetting->addMeanPixelwise << ";" << endl;
		//fSetting_file_matlab << "  fSet.LQ_post_opt="<< fSetting->LQ_post_opt << ";" << endl;
		fSetting_file_matlab << "  fSet.lambdaX="<< fSetting->lambdaX << ";" << endl;
		fSetting_file_matlab << "  fSet.lambdaY="<< fSetting->lambdaY << ";" << endl;
		fSetting_file_matlab << "  fSet.LQ_post_opt_im="<< fSetting->LQ_post_opt_im << ";" << endl;
		fSetting_file_matlab << "  fSet.lambdaX_im="<< fSetting->lambdaX_im << ";" << endl;
		fSetting_file_matlab << "  fSet.lambdaY_im="<< fSetting->lambdaY_im << ";" << endl;
		fSetting_file_matlab << "  fSet.iterMain="<< fSetting->iterMain << ";" << endl;
		//fSetting_file_matlab << "  fSet.ImZ_init_type="<< fSetting->ImZ_init_type << ";" << endl;
		fSetting_file_matlab << "  fSet.doFullImOptWithoutPatRec="<< fSetting->doFullImOptWithoutPatRec << ";" << endl;
		fSetting_file_matlab << "  fSet.set_neg_to_0="<< fSetting->set_neg_to_0 << ";" << endl;
		fSetting_file_matlab << "  fSet.use_estimated_SRFs="<< fSetting->use_estimated_SRFs << ";" << endl;
		fSetting_file_matlab << "  fSet.ImX_sim_mode="<< fSetting->ImX_sim_mode << ";" << endl;
		fSetting_file_matlab << "  fSet.use_LRnorm_for_dic_normalization="<< fSetting->use_LRnorm_for_dic_normalization << ";" << endl;
		cout << endl << "fusionSettings.m written." << endl;
	}else{
		cout << endl << "WARNING: statistics_matlab.m could not be written!" << endl;
	}
	fSetting_file_matlab.close();
	chmod(fname_fSetting.c_str(), S_IRWXU | S_IRWXG | S_IRWXO);
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

	// tmp_mat(0,0) = pSetting->numProcPerPatch;
	// fname_tmp = dir_pSetting + "/" + "numProcPerPatch.csv";
	// write_Mat_to_CSV(&tmp_mat, fname_tmp.c_str());

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

	// tmp_mat(0,0) = glPrms->Ng;
	// fname_tmp = dir_glPrms + "/" + "Ng.csv";
	// write_Mat_to_CSV(&tmp_mat, fname_tmp.c_str());

	tmp_mat(0,0) = glPrms->fDS;
	fname_tmp = dir_glPrms + "/" + "fDS.csv";
	write_Mat_to_CSV(&tmp_mat, fname_tmp.c_str());

	// tmp_mat(0,0) = glPrms->timeMainLoop;
	// fname_tmp = dir_glPrms + "/" + "timeMainLoop.csv";
	// write_Mat_to_CSV(&tmp_mat, fname_tmp.c_str());

	// tmp_mat(0,0) = glPrms->timeFileWrite;
	// fname_tmp = dir_glPrms + "/" + "timeFileWrite.csv";
	// write_Mat_to_CSV(&tmp_mat, fname_tmp.c_str());

	// tmp_mat(0,0) = glPrms->timeDictSelect;
	// fname_tmp = dir_glPrms + "/" + "timeDictSelect.csv";
	// write_Mat_to_CSV(&tmp_mat, fname_tmp.c_str());

	// tmp_mat(0,0) = glPrms->timeDictSelect_avg;
	// fname_tmp = dir_glPrms + "/" + "timeDictSelect_avg.csv";
	// write_Mat_to_CSV(&tmp_mat, fname_tmp.c_str());

	tmp_mat(0,0) = glPrms->numPatchGroups;
	fname_tmp = dir_glPrms + "/" + "numPatchGroups.csv";
	write_Mat_to_CSV(&tmp_mat, fname_tmp.c_str());

	// tmp_mat(0,0) = glPrms->timeTotal;
	// fname_tmp = dir_glPrms + "/" + "timeTotal.csv";
	// write_Mat_to_CSV(&tmp_mat, fname_tmp.c_str());

	// tmp_mat(0,0) = glPrms->numProbPerPatch;
	// fname_tmp = dir_glPrms + "/" + "numProbPerPatch.csv";
	// write_Mat_to_CSV(&tmp_mat, fname_tmp.c_str());

	SpEOMatrixD idxPUH_tmp = glPrms->idxPUH.cast<double>();
	fname_tmp = dir_glPrms + "/" + "idxPUH.csv";
	write_Mat_to_CSV(&idxPUH_tmp, fname_tmp.c_str());

	SpEOMatrixD idxPVH_tmp = glPrms->idxPVH.cast<double>();
	fname_tmp = dir_glPrms + "/" + "idxPVH.csv";
	write_Mat_to_CSV(&idxPUH_tmp, fname_tmp.c_str());

	// fname_tmp = dir_glPrms + "/" + "decMat_C.csv";
	// write_Mat_to_CSV(&glPrms->decMat_C, fname_tmp.c_str());

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
					 //<< "fSet.LQ_post_opt="<<fSetting->LQ_post_opt << ";" << endl
					 << "fSet.LQ_post_opt_im="<<fSetting->LQ_post_opt_im << ";" << endl
					 << "fSet.lambdaX_ABC="<<fSetting->lambdaX_ABC << ";" << endl
					 << "fSet.lambdaY_ABC="<<fSetting->lambdaY_ABC << ";" << endl
					 << "fSet.lambdaZ_ABC="<<fSetting->lambdaZ_ABC << ";" << endl
					 << "fSet.lambdaZ_ABC_in_1st_iter="<<fSetting->lambdaZ_ABC_in_1st_iter << ";" << endl;
				fusionSetup_file_matlab << //<< "fSet.ImZ_init_type="<<fSetting->ImZ_init_type << ";" << endl
					    "fSet.iterMain="<<fSetting->iterMain << ";" << endl
					 << "fSet.doFullImOptWithoutPatRec="<<fSetting->doFullImOptWithoutPatRec << ";" << endl
					 << "fSet.Nc_max="<<fSetting->Nc_max << ";" << endl
					 << "fSet.theta="<<fSetting->theta << ";" << endl
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
					 << "dSet.chBundleFirst="<<dSetting->chBundleFirst << ";" << endl
					 << "dSet.chBundleLast="<<dSetting->chBundleLast << ";" << endl
					 << "dSet.uLFirst="<<dSetting->uLFirst << ";" << endl
					 << "dSet.uLLast="<<dSetting->uLLast << ";" << endl
					 << "dSet.vLFirst="<<dSetting->vLFirst << ";" << endl
					 << "dSet.vLLast="<<dSetting->vLLast << ";" << endl
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
					 <<  endl;
					 // << "glPrms.Ng="<<glPrms->Ng << ";" << endl
					 // << "glPrms.timeMainLoop="<<glPrms->timeMainLoop << ";" << endl
					 // << "glPrms.timeFileWrite="<<glPrms->timeFileWrite << ";" << endl
					 // << "glPrms.timeTotal="<<glPrms->timeTotal << ";" << endl
					 // << "glPrms.timeDictSelect="<<glPrms->timeDictSelect << ";" << endl
					 // << "glPrms.timeDictSelect_avg="<<glPrms->timeDictSelect_avg << ";" << endl
					 // <<  endl;
		 fusionSetup_file_matlab << "%% paths" << endl 
					 << endl;

		cout << endl << "fusionSetup.m written." << endl;
	}else{
		cout << endl << "WARNING: fusionSetup.m could not be written!" << endl;
	}
	fusionSetup_file_matlab.close();
	chmod(fname_fusionSetup.c_str(), S_IRWXU | S_IRWXG | S_IRWXO);
}


void save_evalResults(SpEOPaths *paths, SpEOGlobalParams *glPrms, SpEOFusionSetting *fSetting, SpEOAssessmentMetrics *assMetrics_HR, int iterMain, int numIterMain, bool beforeFullImOpt, bool doFullImOptWithoutPatRec, bool init_image_eval, bool final_evaluation){
	string dir_eval = paths->dir_out + "/" + "eval";
        stringstream numStrSS;
        numStrSS << iterMain;
        string numStr = numStrSS.str();

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
        }else{
                if(beforeFullImOpt){
                         endOfName = "iter" + numStr + "_1_after_patch_rec";
                }else{
                        endOfName = "iter" + numStr + "_2_after_full_im_opt";
                }
        }
	dir_eval=dir_eval + "/" + endOfName;
	cout << "write assessment results to files in directory: " << endl << "     " << dir_eval << " .. ";

	//################################
	//# save SpEOFusionSetting
	//################################
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
	      	
		cout << endl << "evaluation file written." << endl;
	}else{
		cout << endl << "WARNING: statistics_matlab.m could not be written!" << endl;
	}
	eval_file_matlab.close();
	chmod(fname_eval.c_str(), S_IRWXU | S_IRWXG | S_IRWXO);
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

	int iG;
	//if(fSet->useNewMethodForCalculatingZ){
	//	iG= glPrms->myChX[ipp];
	//}else{
		iG= ipp;
	//}

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
					SpEOMatrixD patchHR_tmp = (ImX->get_rasterBands()[ig]->bandDataMatD.block(idxPUH->coeff(u), idxPVH->coeff(v), pszH, pszH));
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
		currP = ImX_LR->get_rasterBands()[iG]->get_bandDataMat()->block((*idxPUL)(uP), (*idxPVL)(vP), pszL, pszL);

		// Arrange current patch representation, currP, into a vector, currV.
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
			int seed = abs(time*121 - my_rank*12345 + iP*iD - uP +vP + iG*79) % 104729;
			// Seed random number
			srand(seed);
			// Compute random index
			int iR = rand() % NP;
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
					if(iN == selectedPatches(sP)){
						anyof=anyof+1;
					}
				}

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
					selfcorrelation.col(0).row(cnt) <<  uN;
					selfcorrelation.col(1).row(cnt) <<  vN;
					selfcorrelation.col(2).row(cnt) <<  sum;
					selfcorrelation.col(3).row(cnt) <<  iN;
					cnt = cnt+1;
				}
			} //<-- for looping through iN patchs to form selfcorrelation matrix
			// Order selfcorrelation matrix
			VectorXi sort_order(sizeSC);
			std::vector<argsort_pair> data(sizeSC);
			for(int i=0;i<sizeSC;i++) {
				data[i].first = i;
				data[i].second = selfcorrelation.col(2)(i);
			}
			// Sort data matrix according to the norms
			std::sort(data.begin(), data.end(), argsort_comp);
			// Get permutations ---> sort_order vector
			for(int i=0;i<sizeSC;i++) {
				sort_order(data[i].first) = i;
			}
			selfcorrelationOrdered = sort_order.asPermutation()*selfcorrelation;
			// the most uncorrelated patch is placed in patch comp ordered
			patchCompOrdered.col(0).row(iD) << selfcorrelationOrdered.col(0).row(0);
			patchCompOrdered.col(1).row(iD) << selfcorrelationOrdered.col(1).row(0);
			patchCompOrdered.col(2).row(iD) << selfcorrelationOrdered.col(2).row(0);
			// Update Selected patches list of indices
			int indexadd = selfcorrelationOrdered(0,3);
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
			VectorXi sort_order(NP);
			std::vector<argsort_pair> data(NP);
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
			VectorXi sort_order(NP);
			std::vector<argsort_pair> data(NP);
			for(int i=0;i<NP;i++) {
				data[i].first = i;
				data[i].second = patchComp->col(2)(i);
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
	currP.resize(pszL,pszL);
	compP.resize(pszL,pszL);
	currP = SpEOMatrixF::Zero(pszL,pszL);
	compP = SpEOMatrixF::Zero(pszL,pszL);
	currenttmp = SpEOMatrixF::Zero(pszL,pszL);
}


void dictSelectFunc(SpEOMatrixD *patchCompSubsetDouble, SpEODataset *ImX_LR, SpEODataset *ImX, SpEODataset *ImY, SpEOFusionSetting *fSet, SpEOGlobalParams *glPrms, int NP, int iP, int uP, int vP, SpEOVectorI *idxPUH, SpEOVectorI *idxPVH, SpEOVectorI *idxPUL, SpEOVectorI *idxPVL, SpEOMatrixD *SRF, SpEOMatrixD *patchComp, int my_rank, int ipp, int &NDP){
	int pszH = fSet->patchsize*glPrms->fDS;
	int pszL = fSet->patchsize;
	int pszL2 = pszL*pszL;
	int pszH2 = pszH*pszH;
	int NPV = glPrms->NPV;
	bool sortON;

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
	//if(fSet->useNewMethodForCalculatingZ){
		iG= ipp;
	//}else{
	//	iG= glPrms->myChX[ipp];
	//}
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
					SpEOMatrixD patchHR_tmp = (ImX->get_rasterBands()[ig]->bandDataMatD.block(idxPUH->coeff(u), idxPVH->coeff(v), pszH, pszH));
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
		currP = ImX_LR->get_rasterBands()[iG]->bandDataMatD.block(idxPUL->coeff(uP), idxPVL->coeff(vP), pszL, pszL);
		// Arrange current patch representation, currP, into a vector, currV.
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
			compP = ImX_LR->get_rasterBands()[iG]->bandDataMatD.block(idxPUL->coeff(uN), idxPVL->coeff(vN), pszL, pszL);
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
			int seed = abs(time*121 - (my_rank+1)*12345 + (iP+1)*iD - uP +vP + (iG+1)*79) % 104729;
			// Seed random number
			srand(seed);
			// Compute random index
			int iR = rand() % NP;
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
				for(int sP=0;sP<iD;sP++){
					if(iN == selectedPatches(sP)){
						anyof=anyof+1;
					}
				}
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
					selfcorrelation.col(0).row(cnt) <<  uN;
					selfcorrelation.col(1).row(cnt) <<  vN;
					selfcorrelation.col(2).row(cnt) <<  sum;
					selfcorrelation.col(3).row(cnt) <<  iN;
					cnt = cnt+1;
				}
			}
			// Order selfcorrelation matrix
			VectorXi sort_order(sizeSC);
			std::vector<argsort_pair> data(sizeSC);
			for(int i=0;i<sizeSC;i++) {
				data[i].first = i;
				data[i].second = selfcorrelation.col(2)(i);
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
			VectorXi sort_order(NP);
			std::vector<argsort_pair> data(NP);
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
			VectorXi sort_order(NP);
			std::vector<argsort_pair> data(NP);
			for(int i=0;i<NP;i++) {
				data[i].first = i;
				data[i].second = patchComp->col(2)(i);
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
/*
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
*/

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
	}
}


bool contains_nan(SpEOMatrixD x){
 return !(( x.array() == x.array() ).all());
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
	  //if(fSetting->useNewMethodForCalculatingZ){
		  chBundleLast_tmp  = dSetting->chBundleLast;  
		  dSetting->chBundleLast  = Ng-1;
	  //}else{
		 // chBundleLast_tmp  = dSetting->chBundleLast;  
		 // dSetting->chBundleLast  = glPrms->numProbPerPatch-1;
	  //}
	  int sizeUH_red_tmp = glPrms->sizeUH_red; glPrms->sizeUH_red = glPrms->sizeUH;
	  int sizeVH_red_tmp = glPrms->sizeVH_red; glPrms->sizeVH_red = glPrms->sizeVH;
	  setMetaInfo(ImX_sim, ImY, ImX, dSetting, glPrms);
	  glPrms->sizeUH_red = sizeUH_red_tmp;
	  glPrms->sizeVH_red = sizeVH_red_tmp;
	  dSetting->chBundleFirst = chBundleFirst_tmp;
	  dSetting->chBundleLast  = chBundleLast_tmp;

	  //if(fSetting->useSimulatedImXforDictLearn){
		  // get window correct indeces
		  int winSizeL = fSetting->winSize;
                  int idxWUL = max(0, idxPUL->coeff(uP)-(int)(0.5*(double)(winSizeL-fSetting->patchsize)) );
                  int idxWVL = max(0, idxPVL->coeff(vP)-(int)(0.5*(double)(winSizeL-fSetting->patchsize)) );
                  int winSizeUL = min((int)ImY->get_sizeU()-1,  idxPUL->coeff(uP)-(int)(0.5*(double)(winSizeL-fSetting->patchsize))  + winSizeL-1   )-idxWUL+1;
                  int winSizeVL = min((int)ImY->get_sizeV()-1,  idxPVL->coeff(vP)-(int)(0.5*(double)(winSizeL-fSetting->patchsize))  + winSizeL-1   )-idxWVL+1;
		  for(int iG=0; iG<Ng; iG++){
			  SpEOMatrixD patX_LR_mod = *patX_LR;
			  // calculate average of all HS bands in current group
			  SpEOVectorD ImY_g_avg_vec;
			  SpEOVectorD CCs = SpEOVectorD::Zero(glPrms->NChX);
			  double c_hat[glPrms->NChX]; // this is the vector of weights for different ImX bands estimated via one of the 3 following modes
			  int b_tilde[glPrms->NChX]; // ordering of ImX bands
			  switch(sim_mode){
				case 0: {
			  		if(localCalculation){
			  		        ImY_g_avg_vec = SpEOVectorD::Zero(patY->cols());
			  		        int iC;
			  		        for(iC=0; iC<Nc_vec[iG]; iC++){
			  		      	  iChY = idxChY[iG]+iC;
			  		      	  ImY_g_avg_vec += patY->row(iChY);
			  		        }
			  		        ImY_g_avg_vec/=(double)(Nc_vec[iG]);
			  		        for(iChX=0; iChX<glPrms->NChX; iChX++){
			  		      	  CCs(iChX) = ImY_g_avg_vec.dot(patX_LR_mod.row(iChX));
			  		        }
			  		        // if all CCs are close to zero, then the corresponding band in ImX_sim will contain only
			  		        // entries that are close to zero which in turn results in dictionary atoms which have norm
			  		        // close to zero. This would cause NaN numbers as the dictionary atoms will be (and have to be)
			  		        // normalized. Hence, we enlarge the window size for this group until this problem is resolved
			  		        if (CCs.maxCoeff()<(1e-8) && CCs.minCoeff()>-(1e-8)){
			  		      	  winSizeL = fSetting->winSize;
			  		      	  int NChY = patY->rows();
			  		      	  int NChX = patX_LR_mod.rows();
			  		      	  bool extraRoundCompleted = false;
			  		      	  cout << "         ["<<my_rank<<"] WWWW Warning:"
			  		      	  								   << " iG="<<iG << " [ImY bands "<< idxChY[iG] << " to " << idxChY[iG]+Nc_vec[iG]-1 << "] -> all CCs (and probably the entire ImY patch in this group) are close to zero!"
			  		      	  								   << " => Enlarge winSizeL from " << winSizeL;
			  		      	  do{
			  		      		  if ( !(CCs.maxCoeff()<(1e-8) && CCs.minCoeff()>-(1e-8))){
			  		      			  extraRoundCompleted = true;
			  		      		  }
			  		      		  cout << " to " << winSizeL+4;
			  		      		  winSizeL += 4;
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
			  		      		  ImY_g_avg_vec = SpEOVectorD::Zero(winY.cols());
			  		      		  for(iC=0; iC<Nc_vec[iG]; iC++){
			  		      			  iChY = idxChY[iG]+iC;
			  		      			  ImY_g_avg_vec += winY.row(iChY);
			  		      		  }
			  		      		  ImY_g_avg_vec/=(double)(Nc_vec[iG]);
			  		      		  for(iChX=0; iChX<glPrms->NChX; iChX++){
			  		      			  CCs(iChX) = ImY_g_avg_vec.dot(patX_LR_mod.row(iChX));
			  		      		  }

			  		      		  if(winSizeUL>ImY->get_sizeU() && winSizeVL>ImY->get_sizeV()){
			  		      			  cerr << "         EEEE Error: Error occurred while enlarging the window for generating ImX_sim.." << endl
			  		      				   << "                     window exceeds the size of the image, which can only mean that " << endl
			  		      				   << "                     ImY contains only zeros in the current group (iG=" << iG << ") of bands.." << endl
			  		      				   << "                     In fact, ImY_g_avg_vec.norm() = " << ImY_g_avg_vec.norm() << endl
			  		      				   << endl;
			  		      			  break;
			  		      		  }
			  		      	  }while ( (CCs.maxCoeff()<(1e-8) && CCs.minCoeff()>-(1e-8)) || !extraRoundCompleted );
			  		      	  cout << endl;
			  		        }
			  		}else{
			  		        SpEOMatrixD ImY_g_avg = SpEOMatrixD::Zero(glPrms->sizeUL,glPrms->sizeVL);
			  		        int iC, iChY;
			  		        for(iC=0; iC<Nc_vec[iG]; iC++){
			  		      	  iChY = idxChY[iG]+iC;
			  		      	  SpEOMatrixD ttmp = ImY->get_rasterBands()[iChY]->get_bandDataMat()->cast<double>();
			  		      	  ttmp.array() -= ttmp.mean();
			  		      	  ttmp.array() /= ttmp.norm();
			  		      	  ImY_g_avg += ttmp;
			  		        }
			  		        ImY_g_avg/=(double)(Nc_vec[iG]);
			  		        // calculate correlation between ImY_g_avg and all bands in ImX_LR
			  		        for(iChX=0; iChX<glPrms->NChX; iChX++){
			  		      	  SpEOMatrixD ttmp = ImX_LR->get_rasterBands()[iChX]->get_bandDataMat()->cast<double>();
			  		      	  ttmp.array() -= ttmp.mean();
			  		      	  ttmp.array() /= ttmp.norm();
			  		      	  CCs(iChX) = calc_Corr_Coeff(ImY_g_avg, ttmp);
			  		        }
			  		        ImY_g_avg_vec = SpEOVectorD::Map(ImY_g_avg.data(),ImY_g_avg.rows()*ImY_g_avg.cols());
			  		}
			  		// sort the CCs in a decreasing manner to find b_tilde[0],...,b_tilde[NChX-1]
			  		std::vector<argsort_pair> CCs_for_sorting(glPrms->NChX);
			  		for(int iChX=0; iChX<glPrms->NChX; iChX++){
			  		        CCs_for_sorting[iChX].first = iChX;
			  		        CCs_for_sorting[iChX].second = CCs(iChX);
			  		}
			  		std::sort(CCs_for_sorting.begin(), CCs_for_sorting.end(), argsort_comp);
			  		for(int iChX=0; iChX< glPrms->NChX; iChX++) {
			  		        b_tilde[CCs_for_sorting[glPrms->NChX-1-iChX].first] = iChX;
			  		}
			  		SpEOVectorD vec2;
			  		if(localCalculation){
			  		        vec2 = patX_LR_mod.row(b_tilde[0]);
			  		}else{
			  		        SpEOMatrixD mat2 = ImX_LR->get_rasterBands()[b_tilde[0]]->get_bandDataMat()->cast<double>();
			  		        vec2 = SpEOVectorD::Map(mat2.data(),mat2.rows()*mat2.cols());
			  		        vec2.array() -= vec2.mean();
			  		        vec2.array() /= vec2.norm();
			  		}

			  		c_hat[0] = 1/(vec2.dot(vec2))*vec2.dot(ImY_g_avg_vec); // least squares solution
			  		SpEOVectorD sum_ci_ImXLR_vec = c_hat[0]*vec2;
			  		for(iChX=1; iChX<glPrms->NChX; iChX++){
			  		        SpEOVectorD vec3;
			  		        if(localCalculation){
			  		      	  vec3 = patX_LR_mod.row(b_tilde[iChX]);
			  		        }else{
			  		      	  SpEOMatrixD mat3 = ImX_LR->get_rasterBands()[b_tilde[iChX]]->get_bandDataMat()->cast<double>();
			  		      	  vec3 = SpEOVectorD::Map(mat3.data(),mat3.rows()*mat3.cols());
			  		      	  vec3.array() -= vec3.mean();
			  		      	  vec3.array() /= vec3.norm();
			  		        }
			  		        c_hat[iChX] = 1/(vec3.dot(vec3))*vec3.dot(ImY_g_avg_vec-sum_ci_ImXLR_vec);  // least squares solution
			  		        sum_ci_ImXLR_vec += c_hat[iChX]*vec3;
			  		}
					break;
			  	}

			  	case 1: { // unconstrained least squares
					// extract window winY from ImY in each channel of the current group 
                                  	SpEOMatrixD winY_g = SpEOMatrixD::Zero(Nc_vec[iG],winSizeUL*winSizeVL);
					int iC;
                                  	for(iC=0; iC<Nc_vec[iG]; iC++){	
						iChY = idxChY[iG]+iC;
			  		      	SpEOMatrixD windowBand = (ImY->get_rasterBands()[iChY]->get_bandDataMat()->block(idxWUL, idxWVL, winSizeUL, winSizeVL)).cast<double>();
			  		      	winY_g.row(iC) = SpEOVectorD::Map(windowBand.data(), winSizeUL*winSizeVL);
                                  	}
					// calculate average ImY patch for current group
					SpEOVectorD winY_g_avg_vec = SpEOVectorD::Zero(winY_g.cols());
					for(iC=0; iC<Nc_vec[iG]; iC++){
						winY_g_avg_vec += winY_g.row(iC);
					}
					ImY_g_avg_vec/=(double)(Nc_vec[iG]);
					// extract window winX_LR from ImX_LR
                                   	SpEOMatrixD winX_LR = SpEOMatrixD::Zero(glPrms->NChX,winSizeUL*winSizeVL);
                                   	for(iChX=0; iChX<glPrms->NChX; iChX++){
                                   	        SpEOMatrixD windowBand = (ImX_LR->get_rasterBands()[iChX]->get_bandDataMat()->block(idxWUL, idxWVL, winSizeUL, winSizeVL)).cast<double>();
                                   	        winX_LR.row(iChX) = SpEOVectorD::Map(windowBand.data(), winSizeUL*winSizeVL);
                                   	}
					// calculate weights via unconstrained least squares
					SpEOMatrixD A = winX_LR.transpose();
					SpEOVectorD b = winY_g_avg_vec;
					SpEOVectorD x = A.jacobiSvd(ComputeThinU | ComputeThinV).solve(b);
					for(iChX=0; iChX<glPrms->NChX; iChX++){
						c_hat[iChX] = x(iChX);
						b_tilde[iChX]=iChX;
					}
					break;
			  	}
			  	case 2: { // nonnegative least squares
					// extract window winY from ImY in each channel of the current group 
                                  	SpEOMatrixD winY_g = SpEOMatrixD::Zero(Nc_vec[iG],winSizeUL*winSizeVL);
					int iC;
                                  	for(iC=0; iC<Nc_vec[iG]; iC++){	
						iChY = idxChY[iG]+iC;
			  		      	SpEOMatrixD windowBand = (ImY->get_rasterBands()[iChY]->get_bandDataMat()->block(idxWUL, idxWVL, winSizeUL, winSizeVL)).cast<double>();
			  		      	winY_g.row(iC) = SpEOVectorD::Map(windowBand.data(), winSizeUL*winSizeVL);
                                  	}
					// calculate average ImY patch for current group
					SpEOVectorD winY_g_avg_vec = SpEOVectorD::Zero(winY_g.cols());
					for(iC=0; iC<Nc_vec[iG]; iC++){
						winY_g_avg_vec += winY_g.row(iC);
					}
					ImY_g_avg_vec/=(double)(Nc_vec[iG]);
					// extract window winX_LR from ImX_LR
                                   	SpEOMatrixD winX_LR = SpEOMatrixD::Zero(glPrms->NChX,winSizeUL*winSizeVL);
                                   	for(iChX=0; iChX<glPrms->NChX; iChX++){
                                   	        SpEOMatrixD windowBand = (ImX_LR->get_rasterBands()[iChX]->get_bandDataMat()->block(idxWUL, idxWVL, winSizeUL, winSizeVL)).cast<double>();
                                   	        winX_LR.row(iChX) = SpEOVectorD::Map(windowBand.data(), winSizeUL*winSizeVL);
                                   	}
					// calculate weights via unconstrained least squares
					SpEOMatrixD A = winX_LR.transpose();
					SpEOVectorD b = winY_g_avg_vec;
					SpEOVectorD x = A.jacobiSvd(ComputeThinU | ComputeThinV).solve(b);
      					bool all_good = NNLS<MatrixXd>::solve(A, b, x);
					if(!all_good){
						cout << "[" << my_rank << "] WARNING: NNLS failed at uP=" << uP << " vP=" << vP << " iG=" << iG << endl 
						     << "The problem failed to solve was Ax=b where x>=0 and " << "A = " << A << endl << endl << "b'=" << endl << b.transpose() << endl; 
					}
					for(iChX=0; iChX<glPrms->NChX; iChX++){
						c_hat[iChX] = x(iChX);
						b_tilde[iChX]=iChX;
					}
					break;
			  	}
				default: {
					break;
				}
			  }
			  // simulate HR image band ImX_sim so that ImX_sim*B*S is (locally at current patch) highly correlated with ImY_g_avg.
			  SpEOMatrixD ImX_sim_band = SpEOMatrixD::Zero(ImX->get_rasterBands()[0]->get_bandDataMat()->rows(),ImX->get_rasterBands()[0]->get_bandDataMat()->cols());
			  for(int iChX=0; iChX<glPrms->NChX; iChX++){
				  ImX_sim_band += c_hat[iChX]*ImX->get_rasterBands()[b_tilde[iChX]]->get_bandDataMat()->cast<double>();
			  }
			  ImX_sim->get_rasterBands()[iG]->bandDataMatD = ImX_sim_band;
		  }
	  // }else{
		 //  if(my_rank==0){
			//   cout << "use original (not simulated) ImX for dictionary learning.." << endl;
		 //  }
		 //  for(int ipp=0; ipp<glPrms->numProbPerPatch; ipp++){
			//   ImX_sim->get_rasterBands()[ipp]->bandDataMatD = ImX->get_rasterBands()[glPrms->myChX[ipp]]->get_bandDataMat()->cast<double>();
		 //  }
	  // }
}

void calc_P_matrices(SpEOVectorD* P_lmd_vecs_loc, int **P_lmd_idx_bl_loc, SpEOMatrixI* P_lmd_idx_row_loc, int Ng, int *Nc_vec, int NChZ, int *idxChY,int my_rank){
	int iG, iChZ, iC, ig;

	for(ig=0; ig<Ng; ig++){
		P_lmd_vecs_loc[ig] = SpEOVectorD::Ones(Nc_vec[ig]);
	}
	int col_idx=0, row_idx;
	for(ig=0; ig<Ng; ig++){
		row_idx = idxChY[ig];
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
	SpEOMatrixD winY_TMP = winY;
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
		for(int iG_sub=0; iG_sub<Ng_grp[iG_lrg]; iG_sub++){
			iG_all_grps++;
			all_grps_idxChY[iG_all_grps].first  = iG_all_grps;
			all_grps_idxChY[iG_all_grps].second = idxChY_grp[iG_lrg][iG_sub];
			all_grps_Nc_vec[iG_all_grps] = Nc_vec_grp[iG_lrg][iG_sub];
			if(my_rank==0 && printOutAll){cout << "&  iG_all_grps=" << iG_all_grps << ": idxChY_grp[iG_lrg=" << iG_lrg << "][iG_sub=" << iG_sub << "] = " << idxChY_grp[iG_lrg][iG_sub]
			                                                                      << "  #  " << Nc_vec_grp[iG_lrg][iG_sub] << " = Nc_vec" << endl;}
		}
		if(my_rank==0 && printOutAll){cout << endl;}
	}
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
		while(iChYV<NChY && !jumpToNextU){
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
	std::vector<argsort_pair> CCmat_grp_vec_for_sorting(A);
	for(int ii=0; ii<A; ii++){
		CCmat_grp_vec_for_sorting[ii].first = ii;
		CCmat_grp_vec_for_sorting[ii].second = CCmat_grp_vec[ii];
	}

	std::sort(CCmat_grp_vec_for_sorting.begin(), CCmat_grp_vec_for_sorting.end(), argsort_comp);
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
  SpEOMatrixD A(200,100); A.setRandom();
  SpEOVectorD b(200); b.setRandom();
  SpEOVectorD x(100);

  return NNLS<MatrixXd>::solve(A, b, x);
}

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

    int NChY = ImY_2D.rows();
    int NChX = ImX_LR_2D.rows();
    int NPix = ImY_2D.cols();
    // generate A
    SpEOMatrixD A(NPix,NChY+2);
    A.block(0,0,NPix,NChY) = ImY_2D.transpose();
    A.rightCols(2) = -SpEOMatrixD::Ones(NPix,2);
    A.rightCols(1) *= -1;
    int iChX;
    for(iChX=0; iChX<NChX; iChX++){
      // generate b
      SpEOVectorD b = ImX_LR_2D.row(iChX);
      // generate x
      SpEOVectorD x(NChY+2);
      if(my_rank == 0) cout << " iChX=" << iChX << ": ";
      bool all_good = NNLS<MatrixXd>::solve(A, b, x);
      if(my_rank == 0) cout << ".. done! ";
      ok &= all_good;

      SRF.row(iChX) = x.head(NChY);
      ImX_shift[iChX] = x.coeff(NChY)-x.coeff(NChY+1);
      if(my_rank == 0) cout << " ImX_shift = " << ImX_shift[iChX]*my_scaling_factor << endl;
    }

    ImY_2D *= my_scaling_factor;
    ImX_LR_2D *= my_scaling_factor;

    for(iChX=0; iChX<NChX; iChX++){
      ImX_shift[iChX] *= my_scaling_factor;
    }
    return ok;
}
