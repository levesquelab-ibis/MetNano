```
-----------------------------------------
 __  __      _   _   _
|  \/  | ___| |_| \ | | __ _ _ __   ___
| |\/| |/ _ \ __|  \| |/ _` |  _ \ / _ \
| |  | |  __/ |_| |\  | (_| | | | | (_) |
|_|  |_|\___|\__|_| \_|\__,_|_| |_|\___/

-----------------------------------------
```

# Description:
`MetNano` is a metagenomic analysis pipeline to analyze FASTQ files from long read Oxford Nanopore Technologies.

# Pipeline
`MetNano` consists of a pipeline made of multiple bioinformatics software from read quality control to taxonomic abundance estimation and community analyses.

<img width="2299" height="2164" alt="Flowchart_MetNano_taxonomie_EN_portrait_v1 4" src="https://github.com/user-attachments/assets/fc432567-14ee-431d-bc48-078e0193cf35" />

# Installation:
`MetNano` and all of its dependencies will be available on bioconda. In the meantime, it is possible to build and use a `MetNano` conda environment with the following commands:
```
# Create dedicated MetNano environment
conda create -n metnano_test
conda activate metnano_test

# Install MetNano from conda recipe
conda build -c conda-forge -c bioconda conda_recipe_MetNano/
conda install --use-local metnano_test -c conda-forge -c bioconda
```
The conda recipe to build `MetNano` will be distributed with the software until available on bioconda.

# Usage:
```
./MetNano -f [FASTQ_dir] -o [output_dir] -d [kraken2_database]

REQUIRED ARGUMENTS:
-f [FASTQ_dir]          : Path to directory containing Nanopore FASTQ files
-o [output_dir]         : Path to the output directory where results will be stored
-d [kraken2_database]   : Path to directory containing database used by Kraken2

OPTIONAL ARGUMENTS:
-u                      : Unzip raw FASTQs from MinKNOW folder output
-r                      : Rename FASTQs using association file (-a)
-a [association_file]   : Path to the association file containing the barcode IDs and associated sample names (required if -r is provided)
-q [min_read_quality]   : minimum Phred average quality score with chopper (default  = 15)
-l [min_read_length]    : minimum read length with chopper (default  = 300)
-H [headcrop]           : Trim N nucleotides from the start of a read with chopper (default: 0)
-T [tailcrop]           : Trim N nucleotides from the end of a read with chopper (default: 0)
-m [reference genome]   : Path to FASTA file of host reference genome (default = none)
-c [confidence kraken2] : Confidence score used by Kraken2 (default = 0.1)
-n [path_to_copy_db]    : Copy kraken2 database into temporay in-memory filesystem for use with --memore-mapping option in Kraken2. Typically use /dev/shm/ if available
-R [min_read_threshold] : number of reads required PRIOR to abundance estimation with Bracken (default: 10)
-s [step]               : Set the step after which MetNano will stop. Possible options are: QC (stops after nanoQC), filtering (stops after chopper), and complete (default = complete)
-t [INT]                : Number of threads (default = 2)

ADDITIONAL OPTIONS:
-h                	     : Display this help message and exit
```
- If the `-u` argument is provided, `MetNano` expects as input FASTQ directory (`-f`) the default output directory from MinKNOW for barcoded samples. See the "**Input FASTQs**" section for more information
- If the `-r` argument is provided, `MetNano` will rename the FASTQ files provided with `-f` using the association file provided with `-a`. The association file must be a TAB-delimited file with the first column being the FASTQ file name without extensions and the second column being the new names to give to the FASTQs. Using the example association file below, `barcode01.fastq` would be renamed `sample14.fastq`.
```
barcode01 sample14
barcode02 sample52
```
- For more information about the minimum read quality (`-q`), minimum read length (`-l`), headcrop (`-H`), and tailcrop (`-T`) parameters please refer to the `chopper` manual.
- If the `-m` argument is provided, the filtered reads from `chopper` will be mapped onto the provided FASTA reference genome using `minimap2` with `-x lr:hp` preset for accurate long reads (error rate <1%). Unmapped reads from the resulting alignment will then be extracted with `samtools` and used for the subsequent steps. This option is useful do deplete the reads from a host genome in a microbiome metagenomic experiment. If this option is not provided, `SeqKit` will directly use the filtered reads from `chopper`.
- `SeqKit` segments the reads with the minimum read length specified with `-l`
- For more information about `Kraken2` databases, confidence score (`-c`) and memory mapping (`-n`) parameters please refer to the `Kraken2` manual.
- `Bracken` estimates taxa relative abundance for the following taxonomic levels: Phylum, Class, Order, Family, and Species. If no reads where directly attributed to any taxon of the specificied taxonomic levels, `Bracken` output files for this specific sample and taxonomic level combination will not be created, but `MetNano` will continue.
- For more information about the minum read threshold used by Bracken (`-T`) parameter please refer to the `Bracken` manual.
- A one-sample BIOM table will be created for every sample and taxonomic level combination.
- The step parameter (`-s`) allows to stop `MetNano` at three different poinst: 1) following the Quality Control steps just after `nanoQC`, 2) following the read filtration step just after `chopper`, and 3) at the end of the pipeline just after `kraken-biom`.
- When analysing multiple samples, each sample will be analysed in parallel from `NanoPlot` to `SeqKit` by dividing the number of provided threads (`-t`) equally between samples. For all other steps every samples are analysed together.
- BIOM tables can be imported in the `phyloseq` `R` package for statistics analyses and visualization. An example `Rscript` is provided to display the kind of statistics and visualization offered by the `phyloseq` package.

# Input FASTQs
`MetNano` can either take as input:
- unzipped FASTQs file with one FASTQ per sample
- raw sequencing output folder from MinKNOW from barcoded samples that have been demultiplexed. The folder must at least have the following structure. WARNING: Provide the entire folder as input (`run_name/` from the example below):
```
run_name/
└── no_sample_id/
    └── run_name_XYZ/
        └── fastq_pass/
            ├── barcode01/
            └── barcode02/
```
All FASTQs/samples must have different names.

# Output directory structure
```
MetNano_results/
├── 00-demultiplexed_fastqs
├── 01-assembly-stats
├── 02-NanoPlot
├── 03-nanoQC
├── 04-chopper
├── 04b-minimap2
├── 05-SeqKit
├── 06-Kraken2
├── 07-Bracken
└── 08-kraken-biom
```
# How to cite:
Coming soon...

# Authors:
Anthony Piot, Jeff Gauthier, and Guillaume Nguyen
