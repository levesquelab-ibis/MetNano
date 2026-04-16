#!/bin/bash

if [ -z "$1" ]; then
    echo "USAGE: rename_barcodes_w_sheet.sh [two-column tab-separated sheet with barcodes and names]"
    echo "EXAMPLE: rename_barcodes_w_sheet.sh sample_sheet.txt in a directory"
    exit 1
fi

while IFS=$'\t' read -r BC ID; do
    # Match files ending with $BC.fastq
    for f in *"$BC".fastq; do
        if [ -e "$f" ]; then
            mv "$f" "${ID}.fastq"
        fi
    done

done < "$1"

mkdir -p removed
shopt -s nullglob
files=( *barcode* )
if [ ${#files[@]} -gt 0 ]; then
    mv "${files[@]}" removed/
fi
