# AncInf

Human local ancestry inference using RFmix (Maples 2013).  

The pipeline is coordinated and run on an HPC using Snakemake.  On UMN's HPC, snakemake can be installed by:

    module load python3/3.6.3_anaconda5.0.1
    conda install -c conda-forge -c bioconda snakemake python=3.6

The 'module load' command will likely need to be run each time prior to use of Snakemake.


Submit jobs to the HPC via:

    snakemake -s workflow/Snakefile -j <XXXX> --cluster 'qsub -l {cluster.l} -M {cluster.M} -A {cluster.A} -m {cluster.m} -o {cluster.o} -e {cluster.e} -r {cluster.r}' --cluster-config workflow/cluster.yaml

where -j specifies the maximum number of jobs to submit at a time.  Importantly, you will need to modify several lines in the *workflow/config.yml* file and the *workflow/cluster.yaml* file.  I have indicated the necessary lines in each file with #MODIFY.

The installation of required programs can be completely avoid by utilizing a Singularity image.  This image is too large to be hosted on Github, although you can find the definitions file used to create the image at *./singularity/Singularity_defs.def*.  Building of images is still not currently supported at MSI, so I used a Vagrant virtual machine, which comes with Singularity pre-configured/installed (https://app.vagrantup.com/singularityware/boxes/singularity-2.4/versions/2.4).  I can also share the img file directly upon request.

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

## Data preparation

We assume that the reference VCF file has been filtered, phased, and indexed.  

The query data for ancestry inference should be provided in PLINK format.  Additionally, SNPs with substantial missing data (e.g. >10%) should be removed from the data.  

## Running the workflow
Clone this repository to the location where you want to store the output of the pipeline.

    git clone https://github.com/pmonnahan/AncInf.git rfmix_test
    cd rfmix_test
    
The critical files responsible for executing the pipeline are contained in the *./workflow* subdirectory contained within the cloned repo.  They are: 

* Snakefile
* config.yml
* cluster.yaml  

The **Snakefile** is the primary workhouse of _snakemake_, which specifies the dependencies of various parts of the pipeline and coordinates their submission as jobs to the MSI cluster.  No modifications to the **Snakefile** are necessary.  However, in order for the **Snakefile** to locate all of the necessary input and correctly submit jobs to the cluster, **both** the config.yaml and cluster.yaml need to be modified.  Open these files and change the entries that are indicated with 'MODIFY'.  

Once these files have been modified, the entire pipeline can be run from within the cloned folder via:

    snakemake --cluster "qsub -l {cluster.l} -M {cluster.M} -A {cluster.A} -m {cluster.m} -o {cluster.o} -e {cluster.e} -r {cluster.r}" --cluster-config workflow/cluster.yaml -j 32
where -j specifies the number of jobs that can be submitted at once.  

The pipeline is currently set up to run on the _small_ queue on the _mesabi_ cluster, which has a per-user submission limit of 500 jobs.  This is more than enough for the entire pipeline, so running with -n 500 will submit all necessary jobs as soon as possible.  If -j is small (e.g. 32), snakemake will submit the first 32 jobs and then submit subsequent jobs as these first jobs finish.

The attractive feature of _snakemake_ is its ability to keep track of the progress and dependencies of the different stages of the pipeline.  Specifically, if an error is encountered or the pipeline otherwise stops before the final step, _snakemake_ can resume the pipeline where it left off, avoiding redundant computation for previously completed tasks.  

To run a specific part of the pipeline, do:

    snakemake -R <rule_name> --cluster "qsub -l {cluster.l} -M {cluster.M} -A {cluster.A} -m {cluster.m} -o {cluster.o} -e {cluster.e} -r {cluster.r}" --cluster-config workflow/cluster.yaml -j 20 --rerun-incomplete
where _rule\_name_ indicates the 'rule' (i.e. job) in the Snakefile that you wish to run.

## Citations
Maples, Brian K et al. “RFMix: a discriminative modeling approach for rapid and robust local-ancestry inference.” American journal of human genetics vol. 93,2 (2013): 278-88. doi:10.1016/j.ajhg.2013.06.020
