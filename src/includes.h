/*
 * includes.h
 *
 *  Created on: Feb 19, 2014
 *      Author: Claas Grohnfeldt
 */

#ifndef INCLUDES_H_
#define INCLUDES_H_

// standard headers
#include <mpi.h>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <iomanip>
#include <time.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <string>
#include <fstream>
#include <sstream>
#include <math.h>

#include <Eigen/SVD>

// only needed for solvertest:
#include <algorithm>

//// GDAL headers
#include "gdal_priv.h"
#include "cpl_conv.h"   // for CPLMalloc()
#include "cpl_string.h" // for creating & copying files (writing output images)
// .. and:
#include <ogr_spatialref.h> // for creating GDALDataset -> geotiff image we need 'OGRSpatialReferenc'

// Eigen headers
#include <Eigen/Dense>
#include <Eigen/Core>
#include <Eigen/Geometry>


// for CSV parser:
#include <iterator>
#include <vector>

// for removing tmp directories
#include <dirent.h>

// this comment can be removed

#endif /* INCLUDES_H_ */
