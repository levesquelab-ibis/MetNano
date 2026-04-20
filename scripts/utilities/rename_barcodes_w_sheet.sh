#!/bin/bash

if [ -z "$1" ]; then
    echo "USAGE: rename_barcodes_w_sheet.sh [two-column tab-separated sheet with barcodes and names]"
    echo "EXAMPLE: rename_barcodes_w_sheet.sh sample_sheet.txt in a directory"
    exit 1
fi

DIR=$2

while IFS=$'\t' read -r BC ID; do
    # Match files ending with $BC.fastq
    for f in $DIR/*"$BC".fastq; do
        if [ -e "$f" ]; then
            mv "$f" $DIR/"${ID}.fastq"
        fi
    done

done < "$1"
