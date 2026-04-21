"""
Store some fixed parameters used by the MetNano pipeline
"""

import os
import sys

## Utilities folder localisation
directory = sys.argv[0].replace("/run_container_MetNano.py", "")
utilities_f = f"{directory}/utilities/"

## Unzip and demultiplex
unzip_barcode_script= os.path.join(utilities_f,  "unzip_barcoded_reads.sh")
rename_barcodes_w_sheet =  os.path.join(utilities_f, "rename_barcodes_w_sheet.sh")

## Default values
# chopper
head_trim=20
tail_trim=20

# minimap2 mapping preset
mapping_preset="lr:hq"
