# Sparse Image Fusion

## The \[SparseFI / J-SparseFI / J-SparseFI-HM\] family of multi-sensor image super-resolution algorithms

### Multi-sensor data fusion for multi- and hyperspectral image super-resolution based on sparse representations


#### Set data paths
```bash
ln -s <PATH_TO_data_JSparseFIHM/HS_MS> <PATH_TO_sparse_image_fusion/data/HS_MS>
ln -s <PATH_TO_data_JSparseFIHM/MS_PAN> <PATH_TO_sparse_image_fusion/data/MS_PAN>
```

#### How to add new data set
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
