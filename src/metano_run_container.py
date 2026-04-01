import sys
import os
import argparse



print("Running metano route",flush=True)


# Create the parser
parser = argparse.ArgumentParser(description="Process some parameters.")


# Positional arguments (not required)
parser.add_argument("-n", "--folder_nanopore", help="Path to the nanopore folder.")
parser.add_argument("-f", "--force_set", help="Force set.")
parser.add_argument("-L", "--list_association", help="association file.")

# Optional arguments with default values
parser.add_argument("-l", "--minlength", help="minlength seq.")
parser.add_argument("-q", "--qualityscore", help="quality seq.")
parser.add_argument("-c", "--confidence", help="quality seq.")
parser.add_argument("-d", "--database", help="quality seq.")
parser.add_argument("-m", "--readmap", help="mapping to ref.")


# Parse the arguments
args = parser.parse_args()

# Access the arguments
folder_nanopore = args.folder_nanopore
list_association=args.list_association
minlength=args.minlength
qualityscore=args.qualityscore
force_set = args.force_set
run_quality = args.list_association
confidence= args.confidence
database= args.database
readmap= args.readmap

database_name= os.path.basename(os.path.normpath(database))

# For debugging, you can print all the variables to verify
print(f"Folder Nanopore: {folder_nanopore}")
print(f"list_association: {list_association}")
print(f"minlength: {minlength}")
print(f"qualityscore: {qualityscore}")
print(f"confidence: {confidence}")
print(f"database: {database}")
print(f"readmap: {readmap}")


# Define your conda environments

# Add the directory containing your scripts to the Python path
sys.path.append('/mnt/ibis/rclevesq/software/nanobatch_wgs/lspq_validation/src')


from class_module_container_metano import *

from parameter_metano import *

##### demultiplexing and QC

unzup = Demultiplexing_qc( folder_nanopore,6,list_association,readmap,qualityscore,minlength)
unzup.run()


sliding_folder = unzup.sliding_dir
trimmed_folder = unzup.trimmed_dir
outputunzfile = unzup.outputunzfile

print(f"Trimmed FASTQ files are in: {trimmed_folder}")
print(f"Sliding FASTQ files are in: {sliding_folder}")
print(f"Output folder for kraken run FASTQ files are in: {outputunzfile}")



kraken_input_dir = os.path.join(outputunzfile)

runner = Kraken2Runner(
    db_path=database,
    db_name=database_name,
    fasta_dir=sliding_folder,
    out_dir=kraken_input_dir,
    threads=6,
    confidence=confidence,
    kraken_env = kraken_env,
    braken_env = braken_env
)

runner.run_all()
