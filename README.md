# AncInf

Human ancestry inference using RFmix

## Requirements

### RFMix
Download/install RFMix:

    git clone https://github.com/slowkoni/rfmix
    cd rfmix
    autoreconf --force --install # creates the configure script and all its dependencies
    ./configure                  # generates the Makefile
    make
    
 Add to PATH:
 
    export PATH=<path_to_rfmix_executable>:$PATH
    
### Plink2
Can be loaded as a module on MSI via:

    module load plink/2.00-alpha-091019
    
### Shapeit4
This program can be very difficult to compile correctly.  The following steps eventually worked for me, but the program will ONLY run on Mesabi nodes.

Clone GitHub repo:

    git clone https://github.com/odelaneau/shapeit4.git
    cd shapeit4
    
Open makefile and change line 2 to:

    CXX=/panfs/roc/msisoft/gcc/7.2.0/bin/g++ -std=c++11

Change lines 5-6 to:

    HTSLIB_INC=/panfs/roc/msisoft/htslib/1.9_gcc-7.2.0_haswell/include/htslib/
    HTSLIB_LIB=/panfs/roc/msisoft/htslib/1.9_gcc-7.2.0_haswell/lib/libhts.a

Change lines 9-11 to:

    BOOST_INC=/panfs/roc/msisoft/boost/1.65.1/gnu-7.2.0/include/
    BOOST_LIB_IO=/panfs/roc/msisoft/boost/1.65.1/gnu-7.2.0/lib/libboost_iostreams.a
    BOOST_LIB_PO=/panfs/roc/msisoft/boost/1.65.1/gnu-7.2.0/lib/libboost_program_options.a
    
To line 36 add:

    -lcurl -lcrypto
    
After making this changes, compile via:

    make
From within the shapeit4 directory.

### Snakemake
The workflow described below uses Snakemake to coordinate the parallel execution of the component tasks.  Snakemake can be installed via pip by entering:

    pip3 install --user snakemake pyaml
at the command line on MSI.

## Running the workflow

