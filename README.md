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
MetNano runs a metagenomic analysis pipeline using provided FASTQ files from Nanopore sequencing.
It is the first step to process the raw basecall data compiled in one file fastQ folder

# Pipeline
`MetNano` consists of a pipeline made of multiple bioinformatics software from read quality control to taxonomic abundance visualization.

<img width="2244" height="1535" alt="Figure 1 - Flowchart" src="https://github.com/user-attachments/assets/07a08090-15c2-4a17-827a-6ca62278f2df" />


# Installation:
MetNano and all of its dependencies are available on bioconda. You can create and use a MetNano conda environment with the following commands:
```
conda create -n metnano -c conda-forge -c bioconda metnano
conda activate metnano
```

# Usage:
```
./MetNano -f [FASTQ_dir] -o [output_dir] -d [kraken2_database]

REQUIRED ARGUMENTS:
-f [FASTQ_directory]    : Path to directory containing Nanopore FASTQ files
-o [output_directory]   : Path to the output directory where results will be stored
-d [kraken2_database]   : Path to directory containing database used by Kraken2

OPTIONAL ARGUMENTS:
-u                      : Unzip raw FASTQs from MinKNOW folder output
-r                      : Rename FASTQs using association file (-a)
-a [association_file]   : Path to the association file containing the barcode IDs and associated sample names (required if -r is provided)
-q [min_read_quality]   : Nanopore quality score filtering parameter used for chopper (default  = 15)
-l [min_read_length]    : Nanopore sequence length filtering parameter used for chopper (default  = 300)
-m [reference genome]   : Path to FASTA file of host reference genome (default = none)
-c [confidence kraken2] : Confidence score used for Kraken2 (default = 0.1)
-n [RAM_disk]           : Path where kraken2 database will be copied into temporay in-memory filesystem for use with --memore-mapping option in Kraken2. Typically use /dev/shm/ if available
-t [INT]                : Number of threads (default = 2)
-h                      : Display this help message and exit
```
- If the `-u` argument is provided, MetNano expects as input FASTQ directory (`-f`) the default output directory from MinKNOW with at least the following structure present. WARNING: Provide the entire folder as input (`run_name/` from the example below):
```
run_name/
└── no_sample_id
    └── run_name_XYZ
        └── fastq_pass
            ├── barcode01
            └── barcode02
```
- If the `-r` argument is provided, MetNano will rename the FASTQ files provided with `-f` using the association file provided with `-a`. The association file must be a TAB-delimited file with the first column being the FASTQ file name without extensions and the second column being the new names you want to give to the FASTQs. Using the example association file below, `barcode01.fastq` would be renamed `sample14.fastq`.
```
barcode01 sample14
barcode02 sample52
```
- For more information about the minimum read quality (`-q`) and minimum read length (`-l`) parameters please refer to the `chopper` manual.
- If the `-m` argument is provided, the filtered reads from `chopper` will be mapped onto the provided FASTA reference genome using `minimap2` with `-x lr:hp` preset for accurate long reads (error rate <1%). Unmapped reads from the resulting alignment will then be extracted with `samtools` and used for the subsequent steps. This option is useful do deplete the reads from a host genome in a microbiome metagenomic experiment. If this option is not provided, `SeqKit` will directly use the filtered reads from `chopper`.
- `SeqKit` segments the reads with the minimum read length specified with `-l`
- For more information about the confidence score (`-c`) and memory mapping (`-n`) parameters please refer to the `Kraken2` manual.
- `Bracken` estimates taxa relative abundance for the following taxonomic levels: Phylum, Class, Order, Family, and Species. If no reads where directly attributed to any taxon of the specificied taxonomic levels, `Bracken` output files for this specific sample and taxonomic level combination will not be created, but MetNano will continue.
- A one-sample BIOM table will be created for every sample and taxonomic level combination.
- BIOM tables can be imported in the `phyloseq` `R` package for statistics analyses and visualization. An example `Rscript` is provided to display the kind of statistics and visualization offered by the `phyloseq` package.


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
Guillaume Nguyen, Jeff Gauthier and Anthony Piot (Levesque Lab)
