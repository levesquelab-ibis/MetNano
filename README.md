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

# Installation:


# Usage:
```
./MetNano -L [association_file] -o [output_dir] -n [dir_with_nanopore_fastqs] -d [kraken2_database]

MANDATORY ARGUMENTS:

-L [association_file]  		: Path to the association file containing the barcode IDs and associated sample names

-o [output_dir]			: Path to the output directory where results will be stored

-n [input_folder_nanopore] 	: Directory containing Nanopore FASTQ files

-d [kraken2_database]   	: Database used by Kraken2

OPTIONAL ARGUMENTS:

-f [T or F]    		: Use previous folder regardless of the date (default = F)

-Q [T or F]    		: Perform the quality check on FASTQs (default = T)

-q [minQ]      		: Nanopore quality score filtering parameter used for chopper (default  = 15)

-l [minlength] 		: Nanopore sequence length filtering parameter used for chopper (default  = 300)

-c [confidence kraken2] : Confidence score used for Kraken2 (default = 0.1)

-m [reference genome]	: FASTA file of host reference genome (default = none)

-u [T or F]		: Unzip raw FASTQs from MinKNOW folder output

ADDITIONAL OPTIONS:

-h                	: Display this help message and exit

EXAMPLE:

./MetNano -L metadata.txt -o output_dir -n /path/to/nanopore/ -d /path/to/kraken2/database/
```
# How to cite:


# Authors:
Guillaume Nguyen, Jeff Gauthier and Anthony Piot (Levesque Lab)
