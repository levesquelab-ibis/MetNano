import os

## Conda settings
CONDA_SOURCE_private = "source /project/rclevesq/software/miniconda3/gqngu_manitou/etc/profile.d/conda.sh"
CONDA_SOURCE_general = "source /project/rclevesq/software/miniconda3/manitou/etc/profile.d/conda.sh" ## RGI
CONDA_SOURCE_other = "source /mnt/ibis/rclevesq/software/miniconda3/jegau138/etc/profile.d/conda.sh" #Mobsuite / GTDBTK

## Conda environments
nanoQC_env="/project/rclevesq/software/miniconda3/manitou/envs/nanoqc_0.10.0"
chopper_env="/project/rclevesq/software/miniconda3/manitou/envs/chopper_env"
nanoplot_env="/project/rclevesq/software/miniconda3/manitou/envs/nanoplot_1.46.2"
minimap2_env="/project/rclevesq/software/miniconda3/manitou/envs/minimap_2.28"
kraken_env='/project/rclevesq/software/miniconda3/gqngu_manitou/envs/kraken2_2_17_1_env'
braken_env='/project/rclevesq/software/miniconda3/manitou/envs/bracken_3.1_ap'
kraken_biom_env='/project/rclevesq/software/miniconda3/manitou/envs/kraken-biom_1.2.0'

## Utilities folder localisation
utilities_f='/mnt/ibis/rclevesq/software/metnano/current/src/utilities/'

## Unzip and demultiplex
unzip_barcode_script= os.path.join(utilities_f,  "unzip_barcoded_reads.sh")
rename_barcodes_w_sheet =  os.path.join(utilities_f, "rename_barcodes_w_sheet.sh")
assembly_stats =  os.path.join(utilities_f, "assembly-stats*")

## Default values
# Chopper
head_trim=20
tail_trim=20

# Mapping to reference
mapping_preset="lr:hq"
