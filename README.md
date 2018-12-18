# Sparse image fusion

This document is structured as follows:

- [Sparse image fusion](#sparse-image-fusion)
  - [Description and literature](#description-and-literature)
  - [Getting started](#getting-started)
    - [Prerequisites](#prerequisites)
    - [Installing third party libraries](#installing-third-party-libraries)
      - [Eigen: C++ template library for linear algebra](#eigen-c-template-library-for-linear-algebra)
      - [GDAL: Geospatial Data Abstraction Library](#gdal-geospatial-data-abstraction-library)
      - [Link libraries to repository's main directory](#link-libraries-to-repositorys-main-directory)
    - [Compile the program](#compile-the-program)
  - [Run the program](#run-the-program)
  - [Flowchart](#flowchart)
  
## Description and literature

This software suite comprises the \[ **SparseFI** - **J-SparseFI** - **J-SparseFI-HM** \] family of multi-sensor image fusion algorithms for multi- and hyperspectral image super-resolution based on sparse representations.

<!--- ## Literature with detailed description of the algorithms --->

A detailed description of the algotihms implemented in this software suite is provided in my dissertation:
> C. Grohnfeldt, "Multi-sensor Data Fusion for Multi- and Hyperspectral Resolution Enhancement Based on Sparse Representations ", Ph.D. Dissertation, Technical University of Munich, 2017; doi:10.14459/2017md1366768

## Getting started

### Prerequisites

This code has been tested on UNIX-based machines only. Its dependencies are as follows: `gcc/g++`, `OpenMPI`, `curl`, `wget`, `grep`.

### Installing third party libraries

This software suite depends on two external libraries, which need to be installed and linked to the repository's main directory.

It is recommended to install those libraries outside of the repo and link them via symbolic links, as instructed further below.

Create a directory, somewhere outside of the repository, into which third party libraries will be installed.

```bash
mkdir <chosen-path-to-thirdparty-library-dir>
```

#### Eigen: C++ template library for linear algebra

The *Eigen* project is hosted on [http://eigen.tuxfamily.org](http://eigen.tuxfamily.org). A git mirrow is available on [GitHub](https://github.com/eigenteam/eigen-git-mirrow). We'll clone that into the above-created lib directory as follows:

```bash
git clone https://github.com/eigenteam/eigen-git-mirror.git <chosen-path-to-thirdparty-library-dir>/eigen
```

That's it for Eigen. [Further below](#Link-paths-to-repository's-main-directory), we will link this library to our repository's directory.

#### GDAL: Geospatial Data Abstraction Library

GDAL can be installed either on a system level or locally following the descriptions on [www.gdal.org](www.gdal.org). Here, we'll build it from source code and install it locally to avoid dependency on sudo permissions.

- Download the source code of the latest stable release:

  ```bash
  cd <chosen-path-to-thirdparty-library-dir>
  mkdir downloads
  cd downloads
  tmp=$(curl http://download.osgeo.org/gdal/CURRENT/ | grep -o "gdal-[2-9].[0-9].[0-9].tar.gz")
  filename=${tmp/.tar.gz*/.tar.gz}
  wget "http://download.osgeo.org/gdal/CURRENT/$filename"
  tar zxf $filename
  ```

- Install locally as follows:

  ```bash
  cd ${filename/.tar.gz/}
  gdal_prefix="<chosen-path-to-thirdparty-library-dir>/gdal"
  ./configure --prefix=$gdal_prefix
  # for compilation on the SuperMUC, you need to load a gcc
  # compiler that is more recent that the default one:
  # module load gcc/8
  make
  make install
  ```

- Test if local installation process was successfull:

  ```bash
  export PATH=${gdal_prefix}/bin:$PATH
  export LD_LIBRARY_PATH=${gdal_prefix}/lib:$LD_LIBRARY_PATH
  export GDAL_DATA=${gdal_prefix}/share/gdal
  gdalinfo --version
  ```

  The last command should output something like the following:\
  `GDAL 2.3.2, released 2018/09/21`.

- Clean up: Navigate to repository and remove downloaded files.

  ```bash
  cd <path-to-repo(sparse_image_fusion)>
  rm -r <chosen-path-to-thirdparty-library-dir>/downloads
  ```

#### Link libraries to repository's main directory

```bash
ln -s <chosen-path-to-thirdparty-library-dir> <path-to-repo(sparse_image_fusion)/lib>
```

### Compile the program

Depending on the underlying maching this program is compiled and run on, you need to modify the `Makefile` by setting the compiler variable `CXX` to either `mpic++` (PC) or `mpiCC` (server).

```bash
module unload mpi.ibm           # <- on SuperMUC only
module load mpi.intel/2018_gcc  # <- on SuperMUC only
module load gcc/8               # <- on SuperMUC only
make
```

## Run the program

Run the script `run.sh`, in which all program arguments are set including specification of paths to the input data set.

```bash
bash run.sh
```

On the SuperMUC, corresponding LoadLeveler scripts need to be created and submitted.

## Flowchart

A technical flowchart describing the J-SpareFI-HM algorithm is given below.

<!---
<object data="http://yoursite.com/the.pdf" type="application/pdf" width="700px" height="700px">
    <embed src="http://yoursite.com/the.pdf">
        <p>This browser does not support PDFs. Please download the PDF to view it: <a href="http://yoursite.com/the.pdf">Download PDF</a>.</p>
    </embed>
</object>
--->
<object data="Mitgliedsbescheinigung.pdf" type="application/pdf" width="700px" height="700px">
    <embed src="Mitgliedsbescheinigung.pdf">
        <p>This browser does not support PDFs. Please download the PDF to view it: <a href="Mitgliedsbescheinigung.pdf">Download PDF</a>.</p>
    </embed>
</object>

other method:

![Alt](Mitgliedsbescheinigung.pdf)