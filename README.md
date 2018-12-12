# Sparse image fusion

This document is structured as follows:

- [Sparse image fusion](#sparse-image-fusion)
  - [Description and literature](#description-and-literature)
  - [Getting started](#getting-started)
    - [Create directories at an accessible location with sufficient storage space](#create-directories-at-an-accessible-location-with-sufficient-storage-space)
    - [Install third party libraries](#install-third-party-libraries)
      - [Eigen: C++ template library for linear algebra](#eigen-c-template-library-for-linear-algebra)
      - [GDAL: Geospatial Data Abstraction Library](#gdal-geospatial-data-abstraction-library)
    - [Link paths to repository's main directory](#link-paths-to-repositorys-main-directory)
  - [How to add a new data set](#how-to-add-a-new-data-set)
  - [platform-dependency](#platform-dependency)
  - [SuperMUC specifics](#supermuc-specifics)
  
## Description and literature

This software suite comprises the \[ **SparseFI** - **J-SparseFI** - **J-SparseFI-HM** \] family of multi-sensor image fusion algorithms for multi- and hyperspectral image super-resolution based on sparse representations.

<!--- ## Literature with detailed description of the algorithms --->

A detailed description of the algotihms implemented in this software suite is provided in my dissertation:
> C. Grohnfeldt, "Multi-sensor Data Fusion for Multi- and Hyperspectral Resolution Enhancement Based on Sparse Representations ", Ph.D. Dissertation, Technical University of Munich, 2017; doi:10.14459/2017md1366768

## Getting started

<!---### Setup directories, links and external libraries--->

### Create directories at an accessible location with sufficient storage space

The following directories should be created at a location, say `<path-to-storage-dir>`, with sufficient storage space, which is accessible from the repository's main directory.

```bash
mkdir <path-to-storage-dir>/lib      # (third party libraries)
mkdir <path-to-storage-dir>/data     # (input data)
mkdir <path-to-storage-dir>/results  # (output data)
mkdir <path-to-storage-dir>/tmp      # (temporary data)
```

Principally, those folders can be placed in different parent directories if preferred.

### Install third party libraries

This software suite depends on two external libraries, which need to be installed and linked to the repository's main directory:

#### Eigen: C++ template library for linear algebra

The project is hosted on [http://eigen.tuxfamily.org](http://eigen.tuxfamily.org). A git mirrow is available on [GitHub](https://github.com/eigenteam/eigen-git-mirrow). We'll clone that into the above-created lib directory as follows:

```bash
git clone https://github.com/eigenteam/eigen-git-mirror.git <path-to-storage-dir>/lib/eigen
```

#### GDAL: Geospatial Data Abstraction Library

download and install from www.gdal.org)

### Link paths to repository's main directory

In order to preserve the directories structure, which is partially hard-coded in `src/paths.cpp`, it is recommended to link those directories as follows:

```bash
ln -s <PATH_TO_your_input_data_dir/HS_MS> <PATH_TO_sparse_image_fusion/data/HS_MS>
ln -s <PATH_TO_your_input_data_dir/MS_PAN> <PATH_TO_sparse_image_fusion/data/MS_PAN>
ln -s <PATH_TO_your_output_data_dir> <PATH_TO_sparse_image_fusion/results>
ln -s <PATH_TO_your_temporary_data_dir> <PATH_TO_sparse_image_fusion/tmp>
```

links to external libraries should be set to the following locations:

```bash
ln -s <PATH_TO_gdal_include_dir> <PATH_TO_sparse_image_fusion/lib/gdal/inc>
ln -s <PATH_TO_gdal_lib_dir> <PATH_TO_sparse_image_fusion/lib/gdal/lib>
ln -s <PATH_TO_eigen_lib_dir> <PATH_TO_sparse_image_fusion/lib/eigen/lib>
```


## How to add a new data set

1. give the data set a unique ID number (can be any string unique to this data set. One option is to use the 12-digit encription that is descriped in the file ```src/paths.cpp```, which corresponds to the sensors used, resolution ratio, SNR, etc.)
2. The directory of the data set should be structured like the ones used in the demo. Pay attention to the symbolic links in the directory "links" with generic names, which should be created for any new data set.
3. modify paths.cpp by adding a line such as

   ```cpp
   }else if(paths->dataSetID_str == "665211108350"){   paths->dir_in = maindir_path + "/" + "HS_MS"  + "/" + "665211108350_ROSIS_Pavia_Univeristy"             + "/" + "InputData" + "/" + "links";
   ```

4. re-compile
5. use the data set ID as a program argument while calling the binary (running the program)

## platform-dependency

1. Makefile: \
   compiler: on the supermuc, the compiler should be set to ```CXX = mpiCC```, while on other platforms you might want to set ```CXX = mpic++```.

## SuperMUC specifics

- Preferred machine: haswell (hw.supermuc.lrz.de)
- LoadLeveler Scripts: modify the following lines:

  ```bash
  #@ initialdir = <your absolut path to the main directory sparse_image_fusion>
  #@ notify_user = <your e-mail address> 
  ```

- compilation: Before compiling the code, make sure to load Intel's impementation of the MPI compiler instead of the standard IBM one. To do so, enter the following two command lines:

  ```bash
  module unload mpi.ibm
  module load mpi.intel/2017
  ```

  Note that those lines are also included in the demo LoadLeveler scripts and should be used in your LoadLeveler scripts as well.
