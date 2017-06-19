hapbinX
=======

`hapbinX` is a modification of the hapbin suite of tools (https://github.com/evotools/hapbin) for calculating an extended version of the [Integrated Haplotype Score (iHS)](http://dx.doi.org/10.1371/journal.pbio.0040072) evolutionary statistic. This new statistic (mmiHS) extends iHS to multiple markers enabling the signatures of selection across pairs of variants to be studied. 

## Tools ##

To run `hapbinX` across a set of variant pairs:

   * `ihs2binsub --hap [.hap/.hapbin file] --map [.map file] --list [list of pairs] --minmaf [0-1] --out [output prefix]` - calculate the multi-marker iHS across pairs of loci specified in a list file
  
For additional options, see `ihs2binsub --help`.

## Copyright and License ##

This code is licensed under the GPL v3. Copyright is retained by the original authors, Colin Maclean and the University of Edinburgh.

## Building from source code ##

### Dependencies ###

   * A C++11 capable compiler (GCC >= 4.7 for required features). OpenMP support required for threaded execution.

   * Optional dependency: MPI for execution on distributed memory systems (clusters/supercomputers).

If any of these are not already installed on your system then for the main Linux distributions they can simply be added via their package managers.

For example to install on **Ubuntu** (tested on 14.04 LTS):

     sudo apt-get update
     sudo apt-get install git cmake libcr-dev mpich libmpich-dev

On **openSUSE** (tested on Enterprise Server 12):

     sudo zypper install cmake git-core gcc-c++ openmpi openmpi-devel
     export PATH=$PATH:/usr/lib64/mpi/gcc/openmpi/lib64:/usr/lib64/mpi/gcc/openmpi/bin

On **Red Hat** (tested on Enterprise Linux 7.1):

     sudo yum install cmake git gcc-c++

For those servers where CMake is not installed, and you do not have the necessary permissions to add it as above, precompiled binaries can be downloaded from [here](http://www.cmake.org/download/).

### Building the source code ###

An out of source build is suggested in order to keep the source directory clean. To do this, check out the hapbin source, move to the build directory, then run `cmake [path to src directory]`. Once CMake has finished generating the necessary files, simply run `make`.

For example:

     git clone https://github.com/evotools/hapbinX.git
     cd hapbinX/build/
     cmake ../src/
     make

### Installing with an alternate toolchain ###

First, check out the hapbin source:
    
    git clone https://github.com/evotools/hapbinX.git
    cd hapbinX/build/

Next, create a `toolchain.cmake` file with the necessary overrides:

    ...
    SET(CMAKE_C_COMPILER "/path/to/c/compiler")
    SET(CMAKE_CXX_COMPILER "/path/to/cxx/compiler")
    ...

`MPI_C_LIBRARIES`, `MPI_CXX_LIBRARIES`, `MPI_C_INCLUDE_PATH`, and `MPI_CXX_INCLUDE_PATH` can be set in this file, too, if necessary.

Then, tell cmake to use this toolchain and build:

    cmake ../src/ -DCMAKE_TOOLCHAIN_FILE=toolchain.cmake
    make

### Input file formats ###

The hap files (`--hap`), containing phased haplotypes, should be in IMPUTE [hap format](https://mathgen.stats.ox.ac.uk/impute/impute_v2.html#-h). IMPUTE provides phased haplotypes in this format for several publically available human cohorts [here](https://mathgen.stats.ox.ac.uk/impute/impute_v2.html#reference). If your data is in VCF format it can be converted to IMPUTE format using [vcftools](https://vcftools.github.io).

The map files (`--map`) should be in the same format as used by [Selscan](https://github.com/szpiech/selscan) with one row per variant and four space-separated columns specifiying chromosome, locus ID, genetic position and physical position.

The list file comprises three tab-separated columns: Chromosome, variant 1 location, variant 2 location.

### Output file formats ###

- Each row of the output file contains information on the variant pairs, the frequencies of each haplotype spanning the variants (AF) and the corresponding IHH and (unstandardised) iHS values. Haplotype alleles are coded using the same notation as in the hap file i.e. the 00 haplotype corresponds to the alleles coded as 0 at the two variants in the hap input file.


