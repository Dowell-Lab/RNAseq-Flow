# Nextflow Implementation of the Dowell Lab Steady State (RNA-seq) Pipeline

### Usage

#### Download and Installation

Clone this repository in your home directory:

    $ git clone git@github.com:Dowell-Lab/RNAseq-Flow.git

Install Nextflow:

    $ curl -s https://get.nextflow.io | bash
    
#### Slurm-Specific Usage Requirements
##### Primary Run Settings

If you are using Linux, this will install nextflow to your home directory. As such, to run Nextflow, you will need add your user home directory to your PATH. Use the following command to set your home directory to your PATH as a variable so you can still access other paths on your cluster without conflict:

    $export PATH=~:$PATH

Secondly, edit the appropriate config file, e.g. `conf/slurm_grch38.config`, to ensure the proper paths are set for genome reference files and other executables (look for all mentions of `COMPLETE_*`). Variable names should hopefully be self-explanatory. An example run with the required arguments is as follows:

```
    $ nextflow run main.nf -profile slurm_grch38 --workdir '</nextflow/work/temp/>'  --outdir '</my/project/>' --email <john.doe@themailplace.com> --sras '</dir/to/sras/*>'
    
```

Directory paths for sras/fastqs must be enclosed in quotes. Notice the name of the configuration file. It's generally a good idea to keep separate configuration files for samples using different reference genomes, and different organisms. The pipeline runs ***paired-end by default***. The --singleEnd flag must be added for all single-end data. While most nascent data is single-end, Groovy configurations make paired-end processing an easier default.

If anything went wrong, you don't need to restart the pipeline from scratch. Instead...

    $ nextflow run main.nf -profile slurm_grch38 -resume
    
To see a full list of options and pipeline version, enter:
    
    $ nextflow run main.nf -profile slurm_grch38 --help

##### Python Package Requirements

***IMPORTANT: For individual users, we highly recommend installing all python packages in a virtual environment***

This pipeline requires a number of optional python packages for qc and analysis. To install RSeQC and MultiQC, you can run the following:

```
$ pip3 install MultiQC --user
$ pip3 install RSeQC --user
```

Note that all packages are Python3.

##### Running Nextflow Using an sbatch script

The best way to run Nextflow is using an sbatch script using the same command specified above. It's advisable to execute the workflow at least in a `screen` session, so you can log out of your cluster and check the progress and any errors in standard output more easily. Nextflow does a great job at keeping logs of every transaction, anyway, should you lose access to the console. The memory requirements do not exceed 8GB, so you do not need to request more RAM than this. SRAs must be downloaded prior to running the pipeline.

## Arguments

**Required Arguments**

| Arugment  | Usage                            | Description                                                          |
|-----------|----------------------------------|----------------------------------------------------------------------|
| -profile  | \<base,slurm\>                    | Configuration profile to use.                                       |
| --fastqs  | \</project/\*\_{R1,R2}\*.fastq.gz\> | Directory pattern for fastq files (gzipped).                      |
| --sras    | \</project/\*.sra\>              | Directory pattern for sra files.                                     |
| --workdir | \</project/tmp/\>                | Nextflow working directory where all intermediate files are saved.   |
| --email   | \<EMAIL\>                        | Where to send workflow report email.                                 |

**Save Options**

| Arguments  | Usage         | Description                                               |
|------------|---------------|-----------------------------------------------------------|
| --outdir   | \</project/\> | Specifies where to save the output from the nextflow run. |
| --savefq   |               | Compresses and saves raw fastq reads.                     |
| --saveTrim |               | Compresses and saves trimmed fastq reads.                 |
| --saveAll  |               | Compresses and saves all fastq reads.                     |
| --skipBAM  |               | Skips saving BAM files (only save CRAM). Default=False    |

**Input File Options**

| Arguments    | Usage       | Description                                                                  |
|--------------|-------------|------------------------------------------------------------------------------|
| --singleEnd  |             | Specifies that the input files are not paired reads (default is paired-end). |
| --flip       |             | Reverse complements each strand. Necessary for some library preps.           |
| --flipR2     |             | Reverse complements R2 only (will not work in singleEnd mode).               |

**Performance Options**

| Arguments       | Usage       | Description                                             |
|-----------------|-------------|---------------------------------------------------------|
| --threadfqdump  |             | Runs multi-threading for fastq-dump for sra processing. |

**QC Options**

| Arguments       | Usage       | Description                                             |
|-----------------|-------------|---------------------------------------------------------|
| --skipMultiQC   |             | Skip running MultiQC.                                   |
| --skipRSeQC     |             | Skip running RSeQC.                                     |

**Analysis Options**

| Arguments       | Usage       | Description                                                                         |
|-----------------|-------------|-------------------------------------------------------------------------------------|
| --count       |               | Count reads (FPKM normalized) over RefSeq gene file. ***Should not be used as stand-alone analysis! Only to be used as a quick first pass.*** |

### Credits

* Margaret Gruca ([@magruca](https://github.com/magruca)): Nextflow pipeline optimization, original pipeline design and optimization

