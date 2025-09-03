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

**WARNING**: For fastqs and sras, make sure you have the path in quotations: e.g. `--fastqs "/path/ot/fastqs/*.fastq.gz"`. Right now, it assumes R1/R2 and 1/2 notation at the END of the file name (e.g. _R1.fastq.gz and _R2.fastq.gz).

## Arguments

**Required Arguments**

| Arugment  | Usage                            | Description                                                          |
|-----------|----------------------------------|----------------------------------------------------------------------|
| -profile  | \<base,slurm\>                    | Configuration profile to use.                                       |
| --fastqs  | \</project/\*\_{R1,R2 or 1,2}.fastq.gz\> | Directory pattern for fastq files (gzipped).                      |
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

**Strandness Options**

| Arguments             | Usage       | Description                                                                  |
|-----------------------|-------------|------------------------------------------------------------------------------|
| --unStranded          |             | Input data will be procssed in HISAT2 as unstranded (default).               |
| --forwardStranded     |             | Indicates data is forward first-stranded.                                    |
| --reverseStranded     |             | Indicates data is reverse first-stranded.                                    |

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

### How do I figure out if mine is forwardStranded, reverseStranded, or unStranded
Run the pipeline on just 1 of your samples (two fastqs if paired). reverseStranded is the most common so try that first. THen do both of the following:
1. Download the TDF file for the sample and visualize on IGV. Visualize a few genes and double check that the reads are following the correct strand of the gene.
2. Check out the file in outdir/qc/rseqc/infer_experiment/ and see where most of the reads fall. If most (>90%) reads are found in one of the following cases, it's likely that stranded type of library. If there is a large mix, you probably have multiple types.

| RSeQC pattern (dominant)        | Library type name                                                                 |
|---------------------------------|-----------------------------------------------------------------------------------|
| 1++, 1--, 2+-, 2-+              | **fr-unstranded** (should use --unStranded)                                                   |
| 1--, 1-+, 2++, 2--              | **fr-firststrand** (also called *reverse*, dUTP-based, should use --reverseStranded)                            |
| 1++, 1--, 2+-, 2-+  | **fr-secondstrand** (forward, rare these days, should use --forwardStranded)                             |

### Example Script
```
#!/bin/bash
#SBATCH --job-name=RNAnextflow # Job name
#SBATCH -p long
#SBATCH --mail-type=ALL # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=hope.townsend@colorado.edu # Where to send mail
#SBATCH --nodes=1 # Run on a single node
#SBATCH --ntasks=1     # Number of CPU (processer cores i.e. tasks) In this example I use 1. I only need one, since none of the commands I run are parallelized.
#SBATCH --mem=8gb # Memory limit
#SBATCH --time=200:00:00 # Time limit hrs:min:sec
#SBATCH --output=/scratch/Shares/clauset/Clauset_ABNexus/e_and_o/EGA_8.27.25_RNAflow.%j.out # Standard output
#SBATCH --error=/scratch/Shares/clauset/Clauset_ABNexus/e_and_o/EGA_8.27.25_RNAflow.%j.err # Standard error log

#activate virtual environment
source /scratch/Shares/clauset/Clauset_ABNexus/venv/bin/activate
echo "Activated environment"
#load modules 
module load sra/2.8.0
module load bbmap/38.05
module load fastqc/0.11.8
module load hisat2/2.1.0
module load samtools/1.8
module load preseq
module load igvtools/2.3.75
module load mpich/3.2.1
module load bedtools/2.28.0
module load openmpi/1.6.4
module load gcc/7.1.0
module load python/3.6.3/rseqc/3.0.0


#Get the Nextflow paths (for script and Nextflow Executive)
SRC=/Users/hoto7260/Flows/RNAseq-Flow
NF_EXE='/scratch/Shares/dowell/dbnascent/pipeline_assets/nextflow'


${NF_EXE} ${SRC}/main.nf -profile slurm --workdir '/scratch/Shares/clauset/Clauset_ABNexus/RNAseq_flow_out/tmp/' --genome_id 'hg38' --outdir '/scratch/Shares/clauset/Clauset_ABNexus/RNAseq_flow_out/EGA_8.27/' --email hope.townsend@colorado.edu --fastqs '/scratch/Shares/clauset/Clauset_ABNexus/EGA/fastqs/*.fastq.gz' --reverseStranded
```


### Credits

* Margaret Gruca ([@magruca](https://github.com/magruca)): Nextflow pipeline optimization, original pipeline design and optimization
* Hope Townsend: Debugging and Nextflow pipeline optimization

