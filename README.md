# Sparse Image Fusion

### Description
This software suite comprises the \[SparseFI - J-SparseFI - J-SparseFI-HM\]-family of multi-sensor image super-resolution algorithms for multi- and hyperspectral image super-resolution based on sparse representations

### Literature with detailed description
A detailed description of the algotihms implemented in this software suite is provided in my dissertation: \
>   C. Grohnfeldt, "Multi-sensor Data Fusion for Multi- and Hyperspectral Resolution Enhancement Based on Sparse Representations ", Ph.D. Dissertation, Technical University of Munich, 2017; http://doi.org/10.14459/2017md1366768

#### Link paths to directories with sufficient storage space
The following 3 sub-directories should be placed at a location with sufficient storage space that is accessible from the repository's main directory:
- `sparse_image_fusion/data` (input data)
- `sparse_image_fusion/results` (output data)
- `sparse_image_fusion/tmp` (temporary data)

In order to preserve the directories structure, which is partially hard-coded in the file `src/paths.cpp`, it is recommended to link those directories as follows:
```bash
ln -s <PATH_TO_your_input_data_dir/HS_MS> <PATH_TO_sparse_image_fusion/data/HS_MS>
ln -s <PATH_TO_your_input_data_dir/MS_PAN> <PATH_TO_sparse_image_fusion/data/MS_PAN>
ln -s <PATH_TO_your_output_data_dir> <PATH_TO_sparse_image_fusion/results>
ln -s <PATH_TO_your_temporary_data_dir> <PATH_TO_sparse_image_fusion/tmp>
```
#### libraries
This software suite depends on external libraries, which need to be linked to the repository's main directory:
1. GDAL - Geospatial Data Abstraction Library (download and install from www.gdal.org)
2. eigen - c++ template library for linear algebra (download and install from http://eigen.tuxfamily.org)

links should be set to the following locations:
```bash
ln -s <PATH_TO_>
```


#### How to add a new data set
1. give the data set a unique ID number (can be any string unique to this data set. One option is to use the 12-digit encription that is descriped in the file ```src/paths.cpp```, which corresponds to the sensors used, resolution ratio, SNR, etc.)

2. The directory of the data set should be structured like the ones used in the demo.
Pay attention to the symbolic links in the directory "links" with generic names, which should be created for any new data set.

3. modify paths.cpp by adding a line such as
```cpp
}else if(paths->dataSetID_str == "665211108350"){   paths->dir_in = maindir_path + "/" + "HS_MS"  + "/" + "665211108350_ROSIS_Pavia_Univeristy"             + "/" + "InputData" + "/" + "links";
```

4. re-compile

5. use the data set ID as a program argument while calling the binary (running the program)

#### platform-dependency
1. Makefile: \ 
   compiler: on the supermuc, the compiler should be set to ```CXX = mpiCC```, while on other platforms you might want to set ```CXX = mpic++```.
   
#### SuperMUC specifics
- Preferred machine: haswell (hw.supermuc.lrz.de)
- LoadLeveler Scripts: modify the following lines: \ 
``` #@ initialdir = <your absolut path to the main directory sparse_image_fusion> ```
``` #@ notify_user = <your e-mail address> ```
- compilation: Before compiling the code, make sure to load Intel's impementation of the MPI compiler instead of the standard IBM one. To do so, enter the following two command lines:
```bash
module unload mpi.ibm
module load mpi.intel/2017
```
Note that those lines are also included in the demo LoadLeveler scripts and should be used in your LoadLeveler scripts as well.
