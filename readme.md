# Germline Enrichment Snakemake

## Introduction

An rewrite of the GermlineEnrichment pipeline incorporating several new features.

* Snakemake workflow manager
* GATK4
* Parallelised variant calling

This should give improved speed as well as giving the ability to restart pipelines from arbitrary points.

## Requirements

Runs on Centos 6/7

## Install

The pipeline uses Conda environments in order to deploy all dependencies. In order to install the base requirements Miniconda should be installed first [1]. Once this has been installed follow the below steps:

```
git clone {repository}

```
```
cd {repository}

```

```
conda env create -f envs/main.yaml

```
```
source activate germline_enrichment_main

```

The pipeline also uses GATK3 for some steps. Installing this requires that you register a jar file to the conda enviroment in question. First download the GATK3 (3.8-1-0) jar from [2]. Once you have done this follow the below steps once you have activated the germline\_enrichment\_main environment.

```

gatk3-register /path/to/GenomeAnalysisTK[-$PKG_VERSION.tar.bz2|.jar]

```

## Run

### Input Requirements

The pipeline requires gzipped fastq files in a directory structure as follows:

```
input
│
└───sample1
│   │   
│   │   sample1_S1_L001_R1_001.fastq.gz
│   │   sample1_S1_L001_R2_001.fastq.gz
│   │   sample1_S1_L002_R1_001.fastq.gz
│   │   sample1_S1_L002_R2_001.fastq.gz
│   │   sample1.variables 
│   │   
└───sample2
│   │   
│   │   sample2_S2_L001_R1_001.fastq.gz
│   │   sample2_S2_L001_R2_001.fastq.gz
│   │   sample2_S2_L002_R1_001.fastq.gz
│   │   sample2_S2_L002_R2_001.fastq.gz
│   │   sample2.variables


```

There can be a variable amount of lanes although all samples must have the same number of lanes. The *.variables file is the output from the AWGL tool https://github.com/mcgml/MakeVariableFiles

An example variables file is included in the repository examples/ directory.

### Configuration

The pipeline requires a YAML configuration file to run. This file specifies the location of resources that the pipeline requires. For example the location of the reference genome and the worksheet id. Examples of configuration files can be found in the configs/ directory.

TO DO: Add description of each variable in config file

To prepare a configuration file for a pipeline run the following steps should be followed.

* Create or get a template config file for a particular assay (e.g. TSO). Example templates are shown in the config directory. This template config contains information about the location of resources required by the pipeline, but does not contain any run specific information such as the worksheet id.

* Run the scripts/make_config.py script to add the run specific information to the template config file. This creates a config file at the location specified by the --output argument.

```
python scripts/make_config.py --config {template_config_location} --input {sample_input_directory} --output {run_id}_config.yaml


```

### Cluster Configuration

The configs/cluster/cluster_config.json file contains information specific to running the pipeline on a cluster. For example how many cores to allocate to a specific rule/job.

### Starting Pipeline

The pipeline can be run locally or on a compute cluster with a job scheduler installed. Currently we use Torque/PBS but the commands shown can be configured for most job scheduling programs [3].

#### Local

To do a dry run of the pipeline first enter the directory containing the pipeline and input and activate the germline\_enrichment\_main environment.

```
snakemake -np

```

To created an image of the Directed Acyclic Graph of the execution steps in the pipeline enter:

```
snakemake --dag | dot -Tsvg > dag.svg

```

To run the whole pipeline enter: 

```
snakemake -p -j {cores} --use-conda

```

Where {cores} is equal to the number of cores you mish to use on the local machine.

To run the pipeline up to a specific output file enter:

```
snakemake -p -j {cores} --use-conda {output_file} 

```
Where {output_file} is the name of the ouput file required.

For example to run just the SNP/Indel calling part of the pipeline with 2 cores type:

```
snakemake -p -j 2 --use-conda output/jointvcf_all_variants_filtered_genotype_vep_roi_meta_nomt/181212_D00501_0264_BHMHMLBCX2_all_variants_filtered_genotype_vep_roi_meta_nomt.vcf

```
Where 181212\_D00501\_0264\_BHMHMLBCX2 is replaced with your run id.

#### Cluster

To run the whole pipeline enter:

```
export PATH="$PATH:/home/transfer/miniconda3/bin/"; snakemake -p --jobs {n_jobs} \
--cluster "qsub -V -o {log_dir} -e {log_dir} -l ncpus={cluster.ncpus} -l walltime={cluster.walltime} -d {run_folder}" \
--directory {run_folder} -s {snakefile} --cluster-config {cluster_config} --latency-wait {latency} --use-conda --restart-times {restart_times}

```

Where:

{jobs} = number of jobs to submit to cluster at same time.

{log_dir} = absolute path of a directory to store log files in.

{run_folder} = The absolute path of the directory with the snakefile and data in.

{snakefile} = The absolute path to the Snakefile being executed.

{cluster_config} = The absoloute path to the cluster configuration file.

{latency} = How long to wait for files to appear due to latency on the shared file system.

{restart_times} = How many times to retry a failed job.

The path to the conda executables should also be edited if necessary.

For example

```
export PATH="$PATH:/home/transfer/miniconda3/bin/"; snakemake -p --jobs 20 \
--cluster "qsub -V -o /share/data/results/181214_D00501_0264_BHMHMLBCX2/IlluminaTruSightOne/gatk_snakemake/logs -e /share/data/results/181214_D00501_0264_BHMHMLBCX2/IlluminaTruSightOne/gatk_snakemake/logs -l ncpus={cluster.ncpus} -l walltime={cluster.walltime} -d /share/data/results/181214_D00501_0264_BHMHMLBCX2/IlluminaTruSightOne/gatk_snakemake/" \
--directory /share/data/results/181214_D00501_0264_BHMHMLBCX2/IlluminaTruSightOne/gatk_snakemake/ \
-s /share/data/results/181214_D00501_0264_BHMHMLBCX2/IlluminaTruSightOne/gatk_snakemake/Snakefile \
--cluster-config configs/cluster/cluster_config.json \
--latency-wait 160 \
--use-conda \
--restart-times 0


```

To avoid redownloading the conda environment for manta on each pipeline run the --conda-prefix flag can also be set with a absolute path to a location in which to store the snakemake created conda environment.

##### Cluster Troubleshooting

* The paths to resources shoule generally be set to the path to the file from the node. For example instead of /data/db/human/gatk/2.8/b37/human\_g1k\_v37.fasta use /state/partition1/db/human/gatk/2.8/b37/human\_g1k\_v37.fasta 

* The path to the input data should be set with the prefix /share/data instead of /data

* Adjust the file latency if jobs fail - i.e. make it larger.


### Restarting Pipeline

If a pipeline stops for whatever reason it can be restarted from the last succesfully completed rule. To run part of a pipeline again delete the files up to the point you wish to rerun and enter the commands specified in the Cluster part of the readme.

## Output

## Testing and Validation

## References

[1] https://conda.io/miniconda.html

[2] https://software.broadinstitute.org/gatk/download/archive

[3] https://snakemake.readthedocs.io/en/stable/executable.html#cluster-execution




