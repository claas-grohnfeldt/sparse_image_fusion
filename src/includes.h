/*
 * includes.h
 *
 *  Created on: Feb 19, 2014
 *      Author: Claas Grohnfeldt
 */

#ifndef INCLUDES_H_
#define INCLUDES_H_

// standard headers
#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <time.h>
#include <unistd.h>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>

#include <Eigen/SVD>

// only needed for solvertest:
#include <algorithm>

//// GDAL headers
#include "cpl_conv.h"    // for CPLMalloc()
#include "cpl_string.h"  // for creating & copying files (writing output images)
#include "gdal_priv.h"
// .. and:
#include <ogr_spatialref.h>  // for creating GDALDataset -> geotiff image we need 'OGRSpatialReferenc'

// Eigen headers
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Geometry>

// for CSV parser:
#include <iterator>
#include <vector>

// for removing tmp directories
#include <dirent.h>

// this comment can be removed

#endif /* INCLUDES_H_ */
