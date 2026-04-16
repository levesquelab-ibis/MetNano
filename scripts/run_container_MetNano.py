import os
import sys
import argparse
from datetime import datetime
from parameter_MetNano import *
from class_module_container_MetNano import *

# Create the parser
parser = argparse.ArgumentParser(description="Process some parameters.")

# Positional arguments (required)
parser.add_argument("-n", "--folder_nanopore", help="Path to the nanopore folder.")
parser.add_argument("-L", "--association_file", help="association file.")
parser.add_argument("-d", "--database", help="quality seq.")

# Optional arguments with default values
parser.add_argument("-l", "--minlength", help="minlength seq.")
parser.add_argument("-q", "--qualityscore", help="quality seq.")
parser.add_argument("-c", "--confidence", help="kraken2 confidence score.")
parser.add_argument("-m", "--mapping", help="mapping to ref.")
parser.add_argument("-f", "--force_set", help="Force set.")
parser.add_argument("-u", "--unzip", help="unzip fastq.")
parser.add_argument("-t", "--threads", help="number of threads.")

# Parse the arguments
args = parser.parse_args()

# Access the arguments
folder_nanopore = args.folder_nanopore
association_file = args.association_file
minlength = args.minlength
qualityscore = args.qualityscore
force_set = args.force_set
run_quality = args.association_file
confidence = args.confidence
database = args.database
mapping = args.mapping
unzip = args.unzip
threads = args.threads

database_name= os.path.basename(os.path.normpath(database))

print("# Starting MetNano with the following parameters:", flush=True)

# For debugging, you can print all the variables to verify
print(f"{datetime.now().isoformat(sep=' ', timespec='seconds')} [INFO] ONT folder: {folder_nanopore}")
print(f"{datetime.now().isoformat(sep=' ', timespec='seconds')} [INFO] Association file: {association_file}")
print(f"{datetime.now().isoformat(sep=' ', timespec='seconds')} [INFO] Read minimum length: {minlength}")
print(f"{datetime.now().isoformat(sep=' ', timespec='seconds')} [INFO] Read minimum quality: {qualityscore}")
print(f"{datetime.now().isoformat(sep=' ', timespec='seconds')} [INFO] Reference genome: {mapping}")
print(f"{datetime.now().isoformat(sep=' ', timespec='seconds')} [INFO] Kraken2 confidence score: {confidence}")
print(f"{datetime.now().isoformat(sep=' ', timespec='seconds')} [INFO] Kraken2 database: {database}")
print(f"{datetime.now().isoformat(sep=' ', timespec='seconds')} [INFO] Unzip FASTQs: {unzip}")
print(f"{datetime.now().isoformat(sep=' ', timespec='seconds')} [INFO] Number of threads: {threads}")

total_cpus = os.process_cpu_count()
print(f"{datetime.now().isoformat(sep=' ', timespec='seconds')} [INFO] Number of CPUs: {total_cpus}")

# Run Demultiplexing and QC pipeline
preprocessed = Preprocessing(
    folder_nanopore,
    threads,
    association_file,
    mapping,
    qualityscore,
    minlength,
    unzip
)

preprocessed.run_all()

# Run read taxonomic classification and taxa abundance estimation pipeline
classified = TaxonomicClassification(
    database,
    database_name,
    preprocessed.seqkit_dir,
    threads,
    confidence,
    minlength
)

classified.run_all()
