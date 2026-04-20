import os
import sys
import argparse
from datetime import datetime
from parameter_MetNano import *
from class_module_container_MetNano import *

# Create the parser
parser = argparse.ArgumentParser(description="Process some parameters.")

# Positional arguments (required)
parser.add_argument("-f", "--FASTQ_dir", help="Path to the ONT folder.")
parser.add_argument("-d", "--database", help="quality seq.")

# Optional arguments with default values
parser.add_argument("-l", "--minlength", help="minlength seq.")
parser.add_argument("-q", "--qualityscore", help="quality seq.")
parser.add_argument("-c", "--confidence", help="kraken2 confidence score.")
parser.add_argument("-m", "--mapping", help="mapping to ref.")
parser.add_argument("-r", "--rename", help="rename FASTQs.")
parser.add_argument("-a", "--association_file", help="association file.")
parser.add_argument("-u", "--unzip", help="unzip FASTQs.")
parser.add_argument("-t", "--threads", help="number of threads.")
parser.add_argument("-n", "--memory_mapping", help="Folder to copy Kraken2 database for --memory_mapping.")

# Parse the arguments
args = parser.parse_args()

# Access the arguments
FASTQ_dir = args.FASTQ_dir
association_file = args.association_file
minlength = args.minlength
qualityscore = args.qualityscore
run_quality = args.association_file
confidence = args.confidence
database = args.database
mapping = args.mapping
rename = args.rename
unzip = args.unzip
threads = args.threads
memory_mapping = args.memory_mapping

print("# Starting MetNano with the following parameters:", flush=True)

# For debugging, you can print all the variables to verify
print(f"{datetime.now().isoformat(sep=' ', timespec='seconds')} [INFO] FASTQ folder: {FASTQ_dir}")
print(f"{datetime.now().isoformat(sep=' ', timespec='seconds')} [INFO] Unzip FASTQs: {unzip}")
print(f"{datetime.now().isoformat(sep=' ', timespec='seconds')} [INFO] Rename FASTQs: {rename}")
print(f"{datetime.now().isoformat(sep=' ', timespec='seconds')} [INFO] Association file: {association_file}")
print(f"{datetime.now().isoformat(sep=' ', timespec='seconds')} [INFO] Read minimum length: {minlength}")
print(f"{datetime.now().isoformat(sep=' ', timespec='seconds')} [INFO] Read minimum quality: {qualityscore}")
print(f"{datetime.now().isoformat(sep=' ', timespec='seconds')} [INFO] Reference genome: {mapping}")
print(f"{datetime.now().isoformat(sep=' ', timespec='seconds')} [INFO] Kraken2 confidence score: {confidence}")
print(f"{datetime.now().isoformat(sep=' ', timespec='seconds')} [INFO] Kraken2 database: {database}")
print(f"{datetime.now().isoformat(sep=' ', timespec='seconds')} [INFO] Kraken2 memory mapping: {memory_mapping}")
total_cpus = os.process_cpu_count()
print(f"{datetime.now().isoformat(sep=' ', timespec='seconds')} [INFO] Number of CPUs: {total_cpus}")
print(f"{datetime.now().isoformat(sep=' ', timespec='seconds')} [INFO] Number of threads: {threads}")

# Run Demultiplexing and QC pipeline
preprocessed = Preprocessing(
    FASTQ_dir,
    threads,
    association_file,
    mapping,
    qualityscore,
    minlength,
    unzip,
    rename
)

preprocessed.run_all()

# Run read taxonomic classification and taxa abundance estimation pipeline
classified = TaxonomicClassification(
    database,
    preprocessed.seqkit_dir,
    threads,
    confidence,
    minlength,
    memory_mapping
)

classified.run_all()
