#  Ancestry Inference

Human local ancestry inference using [RFmix](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3738819/).  The basic workflow is to parse the input PLINK file by chromosome, perform reference-based haplotype phasing on the data using [ShapeIt4](https://odelaneau.github.io/shapeit4/), and, finally, perform local ancestry inference with RFMix.  More information is provided in the _Pipeline Overview_ below.  With the RFMix output, admixture mapping (i.e. associating local ancestry with phenotype) can be accomplished via a separate pipeline found [here](https://github.com/pmonnahan/admixMap).

<!-- START doctoc generated TOC please keep comment here to allow auto update -->
<!-- DON`T EDIT THIS SECTION, INSTEAD RE-RUN doctoc TO UPDATE -->
**Table of Contents**

 - [Requirements](#requirements)
   - [Snakemake](#snakemake)
   - [Singularity](#singularity)
 - [Running the workflow](#running-the-workflow)
   - [Other Notes](#other-notes)
    - [Debugging and error reports](#debugging-and-error-reports)
 - [Pipeline Overview](#pipeline-overview)
   - [Input Data](#input-data)
   - [Output](#output)
   - [Reference population](#reference-population)
   - [Phasing](#phasing)
   - [Local Ancestry Inference](#local-ancestry-inference)

<!-- END doctoc generated TOC please keep comment here to allow auto update -->

![Pipeline DAG](https://github.com/pmonnahan/AncInf/blob/master/Pipeline_DAG.png)

## Requirements

### Snakemake
The pipeline is coordinated and run on an HPC (or locally) using _Snakemake_.  To install snakemake, first create a virtual environment via:
  
    module load python3/3.6.3_anaconda5.0.1
    conda install -c conda-forge mamba
    mamba create -c conda-forge -c bioconda -n <your_environment_name> snakemake
  
This will create a new virtual environment and install `snakemake`.  Then, activate this environment and perform following installations:

    conda activate <your_environment_name>
    conda install numpy yaml pandas

Anytime you need to run the pipeline, activate this environment beforehand via:

    conda activate <environment_name>

If you choose not to create an environment, you must ensure that these packages are installed and available for your python installation.

### Singularity

The installation of the individual programs used throughout this pipeline can be completely avoid by utilizing a Singularity image.  This image is too large to be hosted on Github, although you can find the definitions file used to create the image [here](https://github.com/pmonnahan/AncInf/blob/master/singularity/Singularity_defs.def).  Building of images is still not currently supported at MSI, so I used a Vagrant virtual machine, which comes with Singularity pre-configured/installed (https://app.vagrantup.com/singularityware/boxes/singularity-2.4/versions/2.4).  I can also share the img file directly upon request.

However, in order to utilize the singularity image, _singularity_ must be installed on the HPC.  Currently, the pipeline assumes that _singularity_ will be available as a module and can be loaded into the environment via the command specified in the config.yml file, in the `module` entry under the  `singularity` section.  The default setting will work for MSI at UMN.

Singularity settings in config.yml

    singularity:
      use_singularity: 'true'
      image: '/home/pmonnaha/pmonnaha/singularity/AncestryInference.sif
      module: 'module load singularity'

## Running the workflow

Clone this repository to the location where you want to store the output of the pipeline.

    git clone https://github.com/pmonnahan/AncInf.git rfmix_test
    cd rfmix_test
    
The critical files responsible for executing the pipeline are contained in the _./workflow_ subdirectory contained within the cloned repo.  They are: 

* Snakefile
* config.yml
* cluster.yaml  

The _Snakefile_ is the primary workhouse of snakemake, which specifies the dependencies of various parts of the pipeline and coordinates execution.  No modifications to the _Snakefile_ are necessary.  

In order for the _Snakefile_ to locate all of the necessary input and correctly submit jobs to the cluster, **both** the _config.yaml_ and _cluster.yaml_ need to be modified. Open these files and change the required entries that are indicated with 'MODIFY'.  Other fields do not require modification, although this may be desired given the particulars of the run you wish to implement.  Details on each entry in the config file (e.g. what the program expects in each entry as well as the purpose of the entry) are provided in the _Pipeline Overview_ at the bottom.

The entire pipeline can be executed on a local machine (not recommended) or on an HPC, and the _cluster.yaml_ file is required only for the latter.  For a local run, change the `local_run` entry to `true` under the `run_settings` section of the config file, and launch snakemake from within the parent directory by the simple command:

    snakemake

However, multiple steps in the pipeline have high resource demands, and so are unlikely to be able to be run locally.  This option exists primarily for testing and troubleshooting, so the remainder of the  documentation assumes that the pipeline will be executed on an HPC.  In order to coordinate the use of the HPC, the following modifications to the snakemake command are required:

    snakemake --cluster "qsub -l {cluster.l} -M {cluster.M} -A {cluster.A} -m {cluster.m} -o {cluster.o} -e {cluster.e} -r {cluster.r}" --cluster-config workflow/cluster.yaml -j 32

where -j specifies the number of jobs that can be submitted at once.  Note that the 'qsub' command is specific to the commonly-used PBS scheduler.  To run on a different HPC scheduler, the command would need to be modified accordingly.  For example, to coordinate submission to a slurm scheduler, the following command would be used:

    snakemake --cluster "sbatch --no-requeue --partition={cluster.p} --time={cluster.time} --mem={cluster.mem} --ntasks={cluster.ntasks} --nodes={cluster.nodes} --mail-user={cluster.mail-user} --mail-type={cluster.mail-type} -o {cluster.o} -e {cluster.e} -A {cluster.A}" --cluster-config workflow/cluster_yale.yaml -j 32

Note also that a different _cluster.yaml_ file is required for the different scheduler.  If you open and inspect the _cluster.yaml_ file vs the _cluster_yale.yaml_ file, you will see syntax that is specific to PBS and slurm schedulers, respectively.  

### Other notes

It is recommended that _snakemake_ is run as an interactive session on an HPC.  _Snakemake_ will launch the specified number (via the -j flag) of jobs, and then will hang and wait for them to finish.  As jobs finish (and assuming no errors), _snakemake_ will launch additional jobs keeping the total running jobs at whatever -j is set for.  Although _snakemake_ should not use a lot of memory, it could have long run times, which is generally not advisable on login nodes.  

One attractive feature of _snakemake_ is its ability to keep track of the progress and dependencies of the different stages of the pipeline.  Specifically, if an error is encountered or the pipeline otherwise stops before the final step, _snakemake_ can resume the pipeline where it left off, avoiding redundant computation for previously completed tasks.  To do so, simply resubmit the original _snakemake_ command.

To run a specific part of the pipeline, do:

    snakemake -R <rule_name> --cluster "qsub -l {cluster.l} -M {cluster.M} -A {cluster.A} -m {cluster.m} -o {cluster.o} -e {cluster.e} -r {cluster.r}" --cluster-config workflow/cluster.yaml -j 20 --rerun-incomplete

where _rule\_name_ indicates the 'rule' (i.e. job) in the Snakefile that you wish to run.

Also, it is often very helpful to do a 'dry-run' of the pipeline in which the different steps and dependencies are printed to screen, but no actual jobs are executed.  This can be helpful to ensure that config entries are correct, etc.  To perform a dry-run, do:

    snakemake -nrp

#### Unlocking the working directory

When _snakemake_ is launched it will place a lock on the working directory, such that other _snakemake_ runs are prohibited from starting.  When _snakemake_ finishes or errors out, it will remove this lock.  However, sometimes this lock is not correctly removed.  This can occur, for example, if the VPN drops connection while _snakemake_ is running.  If you receive a "Directory cannot be locked..." error message from _snakemake_ and you are sure that no other _snakemake_ processes are currently running, you can unlock the directory by:

    snakemake --unlock
    
Then, you can run the usual _snakemake_ command to restart the pipeline.
  
#### Debugging and error reports

Should an error be encountered in a job, snakemake will halt the pipeline and indicate in the terminal that an error has occurred.  The offending job will also be printed in red in the terminal window.  More information on why the job failed can be found in the 'stdout' and 'stderr' files that are output to the _'OandE'_ directory and will be labelled with the jobname.

## Pipeline Overview

### Input Data

The pipeline expects as input a single set of PLINK files (.bed, .fam, .bim) that has gone through basic QC steps (missingness, hwe, maf, etc).  I have written QC pipelines for non-imputed and imputed data, which are available [here](https://github.com/pmonnahan/DataPrep) and [here](https://github.com/pmonnahan/DataPrep/tree/master/postImpute), respectively.  It is technically possible to use imputed data in ancestry inference, although this is not widely seen throughout the literature.  

The input PLINK files are specified in the `query` entry within the config file.

    query: "PATH_TO_PLINK_PREFIX" 
    samples: "all"  

The user can also provide a path to a file in the `samples` entry, in which case the program will subset the `query` dataset to include only the samples in the file (one sample per line).  

It is assumed that the query coordinates and chromosome names are consistent with those used in the reference VCF (see below). 

### Output

All output is labelled using the prefix specified in the `outname` entry in the config file.

    outname: "AncInf"

The RFMix results will be output to the _rfmix_ directory that is automatically created. RFMix outputs a number of files, but the most relevant files are those ending in _.Q_ (which contain the global ancestry percentage estimates for each individual) and the files ending in _.msp.tsv_ (which contain the maximum-likelihood ancestry state in each window analyzed; i.e. local ancestry).  The _.Q_ files can be easily filtered to isolate individuals of a given ethnicity, based on user-provided thresholds.

A set of phased BCF files (separated by chromosome) are generated as an intermediate step and are saved to the _phased_ directory.  This directory will also contain the phased BCF of the individuals from the reference population.  

A good initial check that the results make sense is to simply look at the average local ancestry along a chromosome.  A full collection of these images (one for each chromosome) will be created and output into the _chrom_plots_ folder within the master run directory.  These averages should remain fairly stable across the chromosome.  Any large, sudden changes in the dominant ancestral component are indicative of issues in phasing or ancestry inference.  Furthermore, these chromosome plots should be inspected to identify areas of suspect inference.  For example, drastic changes in average ancestry is often observed near centromeres or telomeres.  These can also likely be flagged by low SNP counts in the inferred windows (which is reported in the _.msp.tsv_ files).   

### Reference population 

The reference VCF to be used for phasing as well as for ancestry inference is provided under the `reference` section of the config file.  The pipeline is currently set up to use the 1000Genomes VCF (available [here](https://www.internationalgenome.org/) or by request) for the reference population.  However, any VCF should work in theory as long as the necessary accessory files are provided.

    reference:
      vcf: "PATH_TO_REFERENCE_VCF"
      subpops: "accessory/1000G_PopLabels.txt"
      genmap: "PATH_TO_DATA_SUBDIRECTORY/genetic_map_hg19.txt"
      phased_bcf: 'none`  

There are two required files that need to accompany the reference VCF, and these are provided at the `subpops` and `genmap` entries.  The `subpops` file should be a text file with two columns: sample ID as it appears in the VCF in the first column and the subpopulation label for that sample in the second column.  If using the 1000Genomes VCF, then the `subpop` file was automatically downloaded to the _accessory_ subdirectory, and the `subpop` entry does not need to be changed.  The `genmap` file specifies the genetic map for the reference genome and is too large to be hosted on GitHub.  However, the hg19 genetic map is available [here](https://mathgen.stats.ox.ac.uk/impute/1000GP_Phase3.html) or by request.  The file contains 3 space-delimited columns: chromosome, base position, genetic position.  

It is assumed that the reference VCF file has been filtered, phased, and indexed. The VCF does NOT need to be subsetted to include only the individuals from the desired reference subpopulations.  This is accomplished by the initial steps of the pipeline, using the `subpops` file described above along with the comma-separated lists (no spaces!) in the `ref_pops` and `pop_names` entries under the `rfmix` section of the config file.  

    rfmix:
      ref_pops: "YRI,GWD,ESN,CEU,IBS,TSI" # No spaces!!
      pop_names: "AFR,AFR,AFR,EUR,EUR,EUR" 
  
 Based on the information contained in the `subpops` file described above, individuals corresponding to the subpopulation names provided in `ref_pops` entry are extracted from the reference VCF.  In addition, a new file is created at:
 
     accessory/Population_Map_File.txt
  
  , which re-labels the subsetted individuals with the corresponding value in the `pop_names`.  
  
  There is expected to be a 1:1 ordered correspondence between the subpopulation labels `ref_pops` and the superpopulation names in `pop_names`.  In this example where we are interesting in inferring 2-way admixture between AFR and EUR populations, all YRI, GWD, and ESN individuals would be extracted and re-labelled as AFR individuals, while the CEU, IBS, and TSI individuals would be labelled as EUR individuals.  This scheme was developed to allow for flexibility in the inclusion/exclusion of particular subpopulations.  
   
  RFMix will sample randomly from within these superpopulations to generate the training/test sets needed for the machine learning algorithm.  It is best if the reference individuals from a superpopulation are evenly distributed across subpopulations, so that a single subpopulation does not dominate during the resampling.  

### Phasing
   
The config file has the following options for modifying the behavior of haplotype phasing in ShapeIt4:
   
    phase:
      threads: "12"
      pbwt_depth: "4"
      sequence_data: 'true'

Increasing the `pbwt_depth` may increase the phasing accuracy, but comes at a substantial computational cost.  The `sequence_data` entry should be set to false if the data comes from an array.

### Local Ancestry Inference

In addition to the `ref_pops` and `pop_names`, the `rfmix` section of the config file provides a number of options for modifying the behavior of RFMix.  

    rfmix:
      ref_pops: "YRI,GWD,ESN,CEU,IBS,TSI" # No spaces!!
      pop_names: "AFR,AFR,AFR,EUR,EUR,EUR" 
      generations: "8"
      reanalyze_reference: "true" 
      window_size: "0.02" 
      threads: "12"

 The `generations` entry specifies the number of generations in the past when admixture between the superpopulations is assumed to have begun.  Values used in the literature are typically approximations based off of historical events or genomic dating methods.  [Bryc et al 2015](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4289685/#:~:text=Patterns%20of%20Genetic%20Ancestry%20of%20Self%2DReported%20Latinos&text=On%20average%2C%20we%20estimate%20that,%2C%20and%206.2%25%20African%20ancestry) provide a good reference for African American and Latinx ancestry inference.  For both scenarios, they modelled admixture between Europeans and Native Americans at 11-12 generations ago and subsequent admixture with Africans 6-8 generations ago.  Unfortunately, RFMix only allows the user to specify a single value, so I have used '8' for African Americans (modelling 2-way admixture between AFR and EUR) and '12' for Latinx individuals (modelling 3-way admixture between AFR, EUR, and AMR)
 
 In the case that a set of reference haplotypes may not be of "pure" ancestry and may themselves be somewhat admixed, the option --reanalyze-reference will cause the program to iteratively analyze the reference haplotypes as if they were query haplotypes, in addition to analyzing the query input (see the [RFmix](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3738819/) paper for a more thorough explanation of this procedure).  This is often advised for inferring local ancestry in Latinx populations, where a 3-way AFR, EUR, and AMR admixture is modelled.  However, it is likely not necessary for inferring ancestry in African American populations, where the ancestral populations likely do not contain any admixed individuals.  

The last relevant option is the window size in which ancestry is to be inferred.  This value is specified in centiMorgans (cM).  Default is 0.2 cM, which corresponds to ~100 - 150 kb windows.  For a given window, there is a minimum requirement on the number of SNPs, and windows will be expanded to meet this requirement regardless of the specified window size.   

 
