#!/bin/bash

FASTQ_PATH=$(realpath $1 2>/dev/null)
CPUS=$2

### Greeter #####################################################

function usage {

	echo ""
	echo "Unzips barcoded Oxford Nanopore reads produced by MinKNOW"
	echo "(assuming they were already basecalled on the GridION or the P2solo)"
	echo "_____________________________________________________________________________"
	echo "Usage: ./unzip_barcoded_reads.sh /mnt/ibis/rclevesq/raw_data/[RAW_ONT_FOLDER] [N_CPUS_FOR_UNZIPPING]"
	echo ""
	echo "*** NOTE: Use the whole folder as input, not one of the subfolders! ***"
	echo "*** NOTE: Specify 1 to 40 CPUs for the script to run"
	echo ""

}

# check for missing args and enforce CPU value between 1 and 40
[[ -z $FASTQ_PATH || -z $CPUS ]] && { usage; exit 1; }
[[ 1 -ge $CPUS || 40 -lt $CPUS ]] && { usage; exit 1; }
[[ -z "$(ls $FASTQ_PATH/*/*/fastq_pass)" ]] && { usage; exit 1; }

### Step 1 - unzip and concatenate ##############################

# CREATE OUTPUT DIR
mkdir $(basename $FASTQ_PATH)_unzipped

# UNZIP FILES
for BARCODE in $FASTQ_PATH/*/*/fastq_pass/barcode*; do
	NAME=$(basename $BARCODE)
	srun --time 4:00:00 -c $CPUS pigz -d -p $CPUS -cv $BARCODE/*.fastq.gz >> $(basename $FASTQ_PATH)_unzipped/$NAME.fastq
done
