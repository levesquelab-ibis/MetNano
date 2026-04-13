import os
import sys

## Utilities folder localisation
dir = sys.argv[0].replace("/run_container_MetNano.py", "")
utilities_f = f"{dir}/utilities/"

## Unzip and demultiplex
unzip_barcode_script= os.path.join(utilities_f,  "unzip_barcoded_reads.sh")
rename_barcodes_w_sheet =  os.path.join(utilities_f, "rename_barcodes_w_sheet.sh")
assembly_stats =  os.path.join(utilities_f, "assembly-stats*")

## Default values
# chopper
head_trim=20
tail_trim=20

# minimap2 mapping preset
mapping_preset="lr:hq"
