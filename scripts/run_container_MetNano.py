import sys
import os
import argparse
from class_module_container_MetNano import *
from parameter_MetNano import *
from datetime import datetime

print("# Starting MetNano with the following parameters:", flush=True)

# Create the parser
parser = argparse.ArgumentParser(description="Process some parameters.")

# Positional arguments (required)
parser.add_argument("-n", "--folder_nanopore", help="Path to the nanopore folder.")
parser.add_argument("-L", "--list_association", help="association file.")
parser.add_argument("-d", "--database", help="quality seq.")

# Optional arguments with default values
parser.add_argument("-l", "--minlength", help="minlength seq.")
parser.add_argument("-q", "--qualityscore", help="quality seq.")
parser.add_argument("-c", "--confidence", help="quality seq.")
parser.add_argument("-m", "--readmap", help="mapping to ref.")
parser.add_argument("-f", "--force_set", help="Force set.")

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
print(f"{datetime.now().isoformat(sep=' ', timespec='seconds')} [INFO] ONT folder: {folder_nanopore}")
print(f"{datetime.now().isoformat(sep=' ', timespec='seconds')} [INFO] association file: {list_association}")
print(f"{datetime.now().isoformat(sep=' ', timespec='seconds')} [INFO] minlength: {minlength}")
print(f"{datetime.now().isoformat(sep=' ', timespec='seconds')} [INFO] minqualityscore: {qualityscore}")
print(f"{datetime.now().isoformat(sep=' ', timespec='seconds')} [INFO] read mapping on host genome: {readmap}")
print(f"{datetime.now().isoformat(sep=' ', timespec='seconds')} [INFO] kraken2 confidence score: {confidence}")
print(f"{datetime.now().isoformat(sep=' ', timespec='seconds')} [INFO] kraken2 database: {database}")
print()

# Add the directory containing your scripts to the Python path
sys.path.append('/mnt/ibis/rclevesq/software/nanobatch_wgs/lspq_validation/src')

# Demultiplexing and QC
preprocessed = Preprocessing(
    folder_nanopore,
    6,
    list_association,
    readmap,
    qualityscore,
    minlength)

preprocessed.run_all()

sliding_folder = preprocessed.sliding_dir
trimmed_folder = preprocessed.trimmed_dir
outputunzfile = preprocessed.outputunzfile

kraken_input_dir = os.path.join(outputunzfile)

classified = TaxonomicClassification(
    db_path=database,
    db_name=database_name,
    fasta_dir=sliding_folder,
    out_dir=kraken_input_dir,
    threads=6,
    confidence=confidence,
    kraken_env = kraken_env,
    bracken_env = bracken_env
)

classified.run_all()
