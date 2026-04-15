import os
import sys
import subprocess
import shutil
import glob
import pandas as pd
import tempfile
import math
import re
import pandas as pd
from datetime import datetime
from concurrent.futures import ThreadPoolExecutor, as_completed
from parameter_MetNano import *

############### Class for preprocessing raw ONT reads ###############
class Preprocessing:
    """
    Python class allowing the computation of sequencing statistics, quality
    control, and filtering from from raw reads, followed by read
    segmentation to homogenous length. Optionnaly, filtered reads can also
    be mapped onto a reference genome to deplete host reads.
    """
    def __init__(self, folder_no_sample_id, thread, metadata, readmap, qualityscore=15, minlength=300, unzip = False):
        self.folder_no_sample_id = folder_no_sample_id.rstrip("/")
        self.last_folder = os.path.basename(self.folder_no_sample_id.rstrip("/"))
        self.qualityscore=qualityscore
        self.minlength=minlength
        self.thread = thread
        self.mapping = readmap
        self.metadata = metadata
        self.unzip=unzip
        self.current_dir = os.path.abspath(".")
        self.date = datetime.now().strftime('%Y-%m-%d')
        self.outputunzfile = f"{self.current_dir}"
        self.outputunzfilefastq= os.path.join(self.outputunzfile, "00_demultiplexed_fastqs")

    def log(self, type, message):
            print(f"{datetime.now().isoformat(sep=' ', timespec='seconds')} [{type}] {message}", flush=True)

    def parallel_execution(self, function, fastqs):
        # ---------- Parallel execution ----------
        with ThreadPoolExecutor(max_workers=self.thread) as exe:
                      future_map = {
                           exe.submit(function, fq): fq
                           for fq in fastqs
                          }
                      for future in as_completed(future_map):
                           fq = future_map[future]
                           future.result()  # crash fast, loud, and clear

    def run_unzip(self):
        # ---------- Unzip and demultiplex barcoded samples ----------
        if not os.path.exists(self.outputunzfilefastq):
            self.log("INFO", "Folder unzip {self.outputunzfilefastq} is missing, running {unzip_barcode_script}")
            self.log("INFO", "Running unzip_barcode_script ...")
            unz_cmd = f"{unzip_barcode_script} {self.folder_no_sample_id} {self.thread}"
            self.log("RUN", f"{unz_cmd}")
            subprocess.run(unz_cmd, shell=True, check=True)
            self.log("INFO", "Unzip and demultiplexing done")
            os.makedirs(self.outputunzfilefastq, exist_ok=True)
            unz_cmd = f"mv {self.outputunzfilefastq}/../*_unzipped/*.fastq {self.outputunzfilefastq}"
            subprocess.run(unz_cmd, shell=True, check=True)
            rm_cmd = f"rmdir {self.outputunzfilefastq}/../*_unzipped/"
            subprocess.run(rm_cmd, shell=True, check=True)
            self.log("INFO", f"Moved file to {self.outputunzfilefastq}\n")
        else:
            self.log("INFO", "Folder exists, skipping step")

    def create_ln_file(self) :
        # ---------- Skip demultiplexing ----------
        if not os.path.exists(self.outputunzfilefastq):
            self.log("INFO", "Demultiplexing step not set with --unzip parameter.")
            self.log("INFO", "Direct copy link from original file to working directory")
            os.makedirs(self.outputunzfilefastq, exist_ok=True)
            ln_cmd= f"ln -s {self.folder_no_sample_id}/*.fastq* {self.outputunzfilefastq}/ "
            self.log("RUN", f"{ln_cmd} ")
            subprocess.run(ln_cmd, shell=True, check=True)
            self.log("INFO", "Symbolic link created into working directory\n")
        elif os.path.exists(self.outputunzfilefastq):
            self.log("INFO", "Symbolic link already exists\n")

    def run_rename(self):
        # ---------- Rename file with metadata provided ----------
        barcode_fastqs = glob.glob(os.path.join(self.outputunzfilefastq, "*barcode*.fastq*"))
        fastqs = glob.glob(os.path.join(self.outputunzfilefastq, "*.fastq"))
        if barcode_fastqs:
            self.log("INFO", "Renaming samples with metadata...")
            rename_cmd = f"{rename_barcodes_w_sheet} {self.metadata}"
            self.log("RUN", f"{rename_cmd}")
            subprocess.run(rename_cmd, shell=True, check=True, cwd=self.outputunzfilefastq)
            self.log("INFO", "Renaming samples with metadata done\n")
        elif fastqs and not barcode_fastqs:
            self.log("INFO", "FASTQ files already renamed")
        else:
            self.log("ERROR", "No barcode file exists")
            self.log("ERROR", "Check your FASTQ directory")
            self.log("END", "Exiting MetNano")
            sys.exit(1)

    def assembly_stat(self):
        # ---------- Generate assembly stats files for all samples in ONT folder ----------
        out_dir = os.path.join(self.outputunzfile,"01-assembly-stats" )
        out_file = os.path.join(out_dir ,"assembly-stats.tsv")
        os.makedirs(out_dir, exist_ok=True)
        if os.path.isfile(out_file) and os.path.getsize(out_file) > 0:
            self.log("INFO", "File assembly-stats.tsv already exists")
            return
        self.log("INFO", f"Computing sequencing reads statistics...")
        cmd = f"assembly-stats -t {self.outputunzfilefastq}/*.fastq > {out_file}  "
        self.log("RUN", f" {cmd} in {self.outputunzfile}")
        subprocess.run(cmd, shell=True, check=True, cwd=self.outputunzfile)
        self.log("INFO", f"Sequencing statistics written to {out_file}")

    def run_NanoPlot_parallel(self, fastqs):
        # ---------- Run individual NanoPlot for each sample ----------
        if not fastqs:
            self.log("WARN", "No FASTQ files found, skipping NanoPlot\n")
            return

        def run_NanoPlot_single(fq):
            sample = os.path.basename(fq).replace(".fastq", "")
            nanoplot_out = os.path.join(self.outputunzfile,"02-NanoPlot" ,f"NanoPlot_{sample}")

            # ---------- NanoPlot ----------
            if os.path.isdir(nanoplot_out):
                self.log("INFO", f"NanoPlot already exists for sample {sample}")
            else:
                self.log("INFO", f"Analyzing sample {sample}...")
                nanoplot_cmd = (
                    f"NanoPlot --fastq {fq} "
                    f"--outdir {nanoplot_out} "
                    f"--prefix {sample}_ "
                    f"-t 1"
                )
                self.log("RUN", f"{nanoplot_cmd}")
                subprocess.run(nanoplot_cmd, shell=True, check=True)
                self.log("INFO", f"NanoPlot completed for sample {sample}")

        # ---------- Run in parallel ----------
        self.parallel_execution(run_NanoPlot_single, fastqs)

    def run_nanoQC_parallel(self, fastqs):
        # ---------- Run individual nanoQC for each sample ----------
        fastqs = glob.glob(os.path.join(self.outputunzfilefastq, "*.fastq"))
        if not fastqs:
            self.log("WARN", "No FASTQ files found, skipping nanoQC\n")
            return

        def run_nanoQC_single(fq):
            sample = os.path.basename(fq).replace(".fastq", "")
            nanoQC_out = os.path.join(self.outputunzfile,"03-nanoQC" ,f"nanoQC_{sample}")

            if os.path.isdir(nanoQC_out):
                self.log("INFO", f"nanoQC already exists for sample {sample}")
            else:
                self.log("INFO", f"Analyzing sample {sample}...")
                nanoQC_cmd = (
                    f"nanoQC {fq} --outdir {nanoQC_out}"
                )
                self.log("RUN", f"{nanoQC_cmd}")
                subprocess.run(nanoQC_cmd, shell=True, check=True)
                self.log("INFO", f"nanoQC completed for sample {sample}")

        # ---------- Run in parallel ----------
        self.parallel_execution(run_nanoQC_single, fastqs)

    def run_chopper_parallel(self, fastqs):
        # ---------- Run chopper with selected parameters (see parameter_MetNano.py file) ----------
        # Output directory
        self.chopper_dir = os.path.join(self.outputunzfile, "04-chopper")
        os.makedirs(self.chopper_dir, exist_ok=True)

        def run_chopper_single(fq):
            base_name = os.path.basename(fq)
            if fq.endswith(".fastq.gz"):
                out_fq = os.path.join(
                    self.chopper_dir,
                    base_name.replace(".fastq.gz", "_trimmed.fastq")
                )
            else:
                out_fq = os.path.join(
                    self.chopper_dir,
                    base_name.replace(".fastq", "_trimmed.fastq")
                )

            if not os.path.exists(out_fq):
                self.log("INFO", f"Trimming file {base_name}...")
                cmd = (
                    f"cat {fq} | "
                    f"chopper --headcrop {head_trim} --tailcrop {tail_trim} -t 1 | "
                    f"chopper -q {self.qualityscore} -l {self.minlength} -t 1 "
                    f"> {out_fq}"
                )
                self.log("RUN", f"{cmd}")
                subprocess.run(cmd, shell=True, check=True)
                self.log("INFO", f"Trimming completed for {base_name}")
            else:
                self.log("INFO", f"Trimmed file already exists for {base_name}")

        # ---------- Run in parallel ----------
        self.parallel_execution(run_chopper_single, fastqs)

    def run_minimap2_parallel(self, fastqs):
        # ---------- Run minimap2 with selected parameters (see parameter_MetNano.py file) ----------
        # Output directory
        self.minimap2_dir = os.path.join(self.outputunzfile, "04b-minimap2")
        os.makedirs(self.minimap2_dir, exist_ok=True)

        def run_minimap2_single(fq):
            base_name = os.path.basename(fq)
            if fq.endswith(".fastq.gz"):
                out_fq = os.path.join(
                    self.chopper_dir,
                    base_name.replace(".fastq.gz", "_trimmed.fastq")
                )
            else:
                out_fq = os.path.join(
                    self.chopper_dir,
                    base_name.replace(".fastq", "_trimmed.fastq")
                )
            unmapped_fq = os.path.join(
                self.minimap2_dir,
                base_name.replace(".fastq.gz", "")
                         .replace(".fastq", "")
                         + "_unmapped.fastq"
            )
            if not os.path.exists(unmapped_fq):
                self.log("INFO", f"Mapping reads from {base_name}...")
                cmd_map = (
                    f"minimap2 -a -x {mapping_preset} -t 1 "
                    f"{self.mapping} {out_fq} "
                    f"> {unmapped_fq.replace('_unmapped.fastq', '_alignment.sam')}"
                )
                self.log("RUN", f"{cmd_map}")
                subprocess.run(cmd_map, shell=True, check=True)
                self.log("INFO", f"Mapping completed for {out_fq}")

        # ---------- Run in parallel ----------
        self.parallel_execution(run_minimap2_single, fastqs)

    def run_samtools_parallel(self, fastqs):
        # ---------- Extract unmapped reads with samtools ----------
        # Output directory
        self.minimap2_dir = os.path.join(self.outputunzfile, "04b-minimap2")

        def run_samtools_single(fq):
            base_name = os.path.basename(fq)
            unmapped_fq = os.path.join(
                self.minimap2_dir,
                base_name.replace(".fastq.gz", "")
                         .replace(".fastq", "")
                         + "_unmapped.fastq"
            )
            if not os.path.exists(unmapped_fq):
                self.log("INFO", f"Extracting reads from {os.path.basename(unmapped_fq)}...")
                cmd_samtools = (
                    f"samtools view -bf 4 "
                    f"{unmapped_fq.replace('_unmapped.fastq', '_alignment.sam')} "
                    f"> {unmapped_fq.replace('_unmapped.fastq', '_unmapped.bam')} && "
                    f"samtools fastq "
                    f"{unmapped_fq.replace('_unmapped.fastq', '_unmapped.bam')} "
                    f"> {unmapped_fq}"
                )
                self.log("RUN", f"{cmd_samtools}")
                subprocess.run(cmd_samtools, shell=True, check=True)
                self.log("INFO", f"Reads extracted from {unmapped_fq.replace('_unmapped.fastq', '_unmapped.bam')}")
            else:
                self.log("INFO", f"samtools output already exists: {unmapped_fq}")

        # ---------- Run in parallel ----------
        self.parallel_execution(run_samtools_single, fastqs)

    def run_SeqKit_parallel(self, fastqs):
        # ---------- Run  seqkit with selected parameters (see parameter_MetNano.py file) ----------
        # Output directory
        self.seqkit_dir = os.path.join(self.outputunzfile, "05-SeqKit")
        os.makedirs(self.seqkit_dir, exist_ok=True)

        def run_SeqKit_single(fq):
            base_name = os.path.basename(fq)
            if self.mapping in [None, "None"]:
                if fq.endswith(".fastq.gz"):
                    input_fq = os.path.join(
                        self.chopper_dir,
                        base_name.replace(".fastq.gz", "_trimmed.fastq")
                    )
                else:
                    input_fq = os.path.join(
                        self.chopper_dir,
                        base_name.replace(".fastq", "_trimmed.fastq")
                    )
            else:
                if fq.endswith(".fastq.gz"):
                    input_fq = os.path.join(
                        self.minimap2_dir,
                        base_name.replace(".fastq.gz", "_unmapped.fastq")
                    )
                else:
                    input_fq = os.path.join(
                        self.minimap2_dir,
                        base_name.replace(".fastq", "_unmapped.fastq")
                    )
            output_fq = os.path.join(
                self.seqkit_dir,
                base_name.replace(".fastq.gz", "")
                         .replace(".fastq", "")
                         + "_segmented.fastq"
            )
            if os.path.exists(output_fq):
                self.log("INFO", f"Segmented file already exists for {input_fq}")
                return
            else:
                self.log("INFO", f"Segmenting reads for {input_fq}...")
                cmd = (
                    f"seqkit sliding "
                    f"-W {self.minlength} "
                    f"-s {self.minlength} "
                    f"{input_fq} "
                    f"> {output_fq}"
                )
                self.log("RUN", f"{cmd}")
                subprocess.run(cmd, shell=True, check=True)
                self.log("INFO", f"Segmenting completed for {input_fq}")

        # ---------- Run in parallel ----------
        self.parallel_execution(run_SeqKit_single, fastqs)

    def run_all(self):

        # Unzip and concatenate FASTQ in raw ONT folder if --unzip True
        if self.unzip == "True":
            print("# Unzip FASTQ files", flush = True)
            self.run_unzip()
        else:
            print("# Demultiplexing", flush=True)
            self.create_ln_file()

        # Rename FASTQ files
        print(f"# Renaming barcodes", flush=True)
        self.run_rename()

        # Get FASTQ file list for task parallelization
        fastqs = glob.glob(os.path.join(self.outputunzfilefastq, "*.fastq*"))
        self.log("INFO", f"FASTQ file found: {len(fastqs)}\n")

        # Run assembly-stats
        print(f"# Compute sequencing statistics with assembly-stats", flush=True)
        self.assembly_stat()

        # Run NanoPlot
        print("\n# Quality control with NanoPlot", flush = True)
        self.run_NanoPlot_parallel(fastqs)

        # Run nanoQC
        print("\n# Quality control with nanoQC", flush = True)
        self.run_nanoQC_parallel(fastqs)

        # Run chopper
        print("\n# Read filtering and trimming with chopper", flush=True)
        self.run_chopper_parallel(fastqs)

        # If a reference genome is specified
        if hasattr(self, 'mapping') and self.mapping not in [None, "None"]:
            # Run minimap2
            print(f"\n# Mapping reads to reference genome with minimap2", flush = True)
            self.run_minimap2_parallel(fastqs)

            # Run samtools
            print(f"# Extraction of unmapped reads with samtools", flush = True)
            self.run_samtools_parallel(fastqs)

        # Run SeqKit
        print("\n# Read segmentation with SeqKit", flush=True)
        self.run_SeqKit_parallel(fastqs)

########## Class for read taxonomic classification and abundance estimation ##########
class TaxonomicClassification:
        """
        Python class allowing the assignment of taxonomic classification to reads with
        Kraken2 using a temporary DB copy in RAM (/dev/shm) for speed, the estimation
        of taxa abundance with Bracken and the conversion of Kraken2-style reports into
        BIOM tables with kraken-biom.
        """
        def __init__(self, db_path, db_name, fastq_dir, out_dir, threads=6, confidence=0.1, minlength=300):
            self.db_path = db_path
            self.db_name = db_name
            self.fastq_dir = fastq_dir
            self.out_dir = out_dir
            self.threads = threads
            self.confidence = confidence
            self.minlength = minlength
            self.temp_db_path = os.path.join("/dev/shm", self.db_name)
            self.base_kraken = os.path.join(self.out_dir, "06-Kraken2", self.db_name)
            self.base_bracken = os.path.join(self.out_dir, "07-Bracken", self.db_name)
            self.base_biom = os.path.join(self.out_dir, "08-kraken-biom", self.db_name)

            #os.makedirs(self.out_dir, exist_ok=True)

        def log(self, type, message):
            print(f"{datetime.now().isoformat(sep=' ', timespec='seconds')} [{type}] {message}", flush=True)

        def copy_db_to_ramdisk(self):
            # ---------- Load Kraken2 database in ram disk (/dev/shm) ----------
            print(f"# Load Kraken2 database into RAM disk", flush = True)
            os.makedirs(self.temp_db_path, exist_ok=True)
            self.log("INFO", "Copying database to /dev/shm/...")
            cmd= (
                 f"rsync -vP "
                 f"{self.db_path}/*  "
                 f"{self.temp_db_path}"
            )
            required_files = ["hash.k2d", "opts.k2d", "taxo.k2d"]
            subprocess.run(
                cmd, shell=True, check=True
            )
            for f in required_files:
               if not os.path.exists(os.path.join(self.temp_db_path, f)):
                   raise RuntimeError(f"Database copy incomplete: missing {f}")
            self.log("INFO", "Database copy completed\n")

        def run_kraken2(self, fastq_file):
            # --------- Kraken2 command ----------
            """Run Kraken2 on a single FASTQ file"""
            print("# Read taxonomic classification with Kraken2", flush = True)
            os.makedirs(self.base_kraken, exist_ok=True)
            sample_name = os.path.splitext(os.path.basename(fastq_file))[0]
            sample_name = re.sub(r"_sliding(_[0-9.]+)?", "", sample_name, flags=re.IGNORECASE)
            prefix = f"{sample_name}_{self.confidence}"
            kraken_output = os.path.join(self.base_kraken, f"{prefix}.all")
            kraken_report = os.path.join(self.base_kraken, f"{prefix}.report")
            bracken_output = os.path.join(self.base_bracken, f"{prefix}.bracken")
            bracken_report = os.path.join(self.base_bracken, f"{prefix}.bracken.report")

            def file_exists_and_not_empty(path):
                return os.path.exists(path) and os.path.getsize(path) > 0

            # Skip logic
            if file_exists_and_not_empty(kraken_output) and file_exists_and_not_empty(kraken_report):
                self.log("INFO", f"Kraken2 output already exists for {prefix}\n")
            else:

                # Check that required files exist
                required_files = ["hash.k2d", "opts.k2d", "taxo.k2d"]
                missing_files = [f for f in required_files if not os.path.exists(os.path.join(self.temp_db_path, f))]
                if missing_files:
                    raise RuntimeError(f"[ERROR] Missing required DB files in RAM: {missing_files}")

                # ----- Kraken2 -----
                kraken_cmd = (
                    f"kraken2 "
                    f"--db {self.temp_db_path} "
                    f"--threads 1 "
                    f"--output {kraken_output} "
                    f"--report {kraken_report} "
                    f"--confidence {self.confidence} "
                    f"--memory-mapping "
                    f"{fastq_file}"
                )

                self.log("INFO", f"Running Kraken2 for {prefix}")
                subprocess.run(kraken_cmd, shell=True, check=True)
                self.log("END", f"Kraken2 ended\n")

                # Post-run check
                if not file_exists_and_not_empty(kraken_output) or not file_exists_and_not_empty(kraken_report):
                    raise RuntimeError(f"[ERROR] Kraken2 run did not produce expected output files for {prefix}")

        def run_bracken(self, fastq_file):
            # ----- Bracken (all taxonomic levels) -----
            print("# Taxon abundance estimation with Bracken", flush = True)
            os.makedirs(self.base_bracken, exist_ok=True)
            sample_name = os.path.splitext(os.path.basename(fastq_file))[0]
            sample_name = re.sub(r"_sliding(_[0-9.]+)?", "", sample_name, flags=re.IGNORECASE)
            prefix = f"{sample_name}_{self.confidence}"
            kraken_output = os.path.join(self.base_kraken, f"{prefix}.all")
            kraken_report = os.path.join(self.base_kraken, f"{prefix}.report")
            bracken_output = os.path.join(self.base_bracken, f"{prefix}.bracken")
            bracken_report = os.path.join(self.base_bracken, f"{prefix}.bracken.report")

            BRACKEN_LEVELS = {
                "species": "S",
                "genus": "G",
                "family": "F",
                "order": "O",
                "class": "C",
                "phylum": "P"
            }
            exceptions = []
            for level_name, level_code in BRACKEN_LEVELS.items():
                bracken_output_s = os.path.join(
                    self.base_bracken,
                    f"{prefix}.bracken.{level_name}.tsv"
                )
                bracken_report_s = os.path.join(
                    self.base_bracken,
                    f"{prefix}.bracken.{level_name}.report"
                )

                # Skip if already exists
                if os.path.exists(bracken_output_s) and os.path.getsize(bracken_output_s) > 0:
                    self.log("INFO", f"Bracken {level_name} file already exists for sample {sample_name}")
                    continue

                bracken_cmd = (
                    f"bracken "
                    f"-d {self.db_path} "
                    f"-i {kraken_report} "
                    f"-o {bracken_output_s} "
                    f"-w {bracken_report_s} "
                    f"-r {self.minlength} "
                    f"-l {level_code}"
                )
                self.log("INFO", f"Running Bracken at {level_name} level for sample {sample_name}...")
                self.log("RUN", bracken_cmd)
                try:
                    subprocess.run(bracken_cmd, shell=True, check=True)
                except subprocess.CalledProcessError as e:
                    self.log("ERROR", f"Bracken failed at {level_name} level for sample {sample_name}")
                    exceptions.append(level_name)
                print()
            if not exceptions:
                self.log("INFO", "Bracken completed for all taxonomic levels\n")
            else:
                self.log("INFO", f"Bracken completed for all taxonomic levels, except for {', '.join(exceptions)} levels\n")

        def run_krakenbiom(self):
            # ---------- Generate BIOM tables, one per report ----------
            print("# Converting Kraken2-style reports to BIOM format with kraken-biom", flush = True)
            os.makedirs(self.base_biom, exist_ok=True)
            report_files = glob.glob(os.path.join(self.base_bracken, "*.report"))
            bracken_reports = report_files
            if not bracken_reports:
                self.log("ERROR", "No Bracken reports found. Cannot generate BIOM tables")
                sys.exit(1)
            OUTPUT_DIR_biom = os.path.join(self.out_dir, "08-kraken-biom", f"{self.db_name}")
            os.makedirs(OUTPUT_DIR_biom, exist_ok=True)
            for report in bracken_reports:
                sample_name = os.path.splitext(os.path.basename(report))[0]
                output_file = os.path.join(OUTPUT_DIR_biom, f"{sample_name}.biom")

                # Skip if BIOM already exists
                if os.path.exists(output_file) and os.path.getsize(output_file) > 0:
                    self.log("INFO", f"BIOM for {sample_name} already exists")
                    continue

                # Build command for this single file
                cmd_biom = (
                    f"kraken-biom {report} "
                    f"--fmt json "
                    f"--output_fp {output_file} "
                    f"--max D --min S"
                )
                self.log("INFO", f"Generating BIOM from Kraken2-style reports for {sample_name}...")
                subprocess.run(cmd_biom, shell=True, check=True)

            self.log("INFO", f"BIOM tables generated\n")

        def run_all(self, max_workers=1):
                fastq_files = glob.glob(os.path.join(self.fastq_dir, "*.fastq*"))
                if not fastq_files:
                    self.log("ERROR", "No FASTQ files found.")
                    return
                files_to_run = []
                for f in fastq_files:
                    sample_name = os.path.splitext(os.path.basename(f))[0]
                    sample_name = re.sub(r"_sliding(_[0-9.]+)?", "", sample_name, flags=re.IGNORECASE)
                    expected_report = os.path.join(self.base_kraken, f"{sample_name}_{self.confidence}.report")
                    if not os.path.exists(expected_report) or os.path.getsize(expected_report) == 0:
                        files_to_run.append(f)

                # Only load RAM if needed
                if files_to_run:
                    self.log("INFO", f"{len(files_to_run)} samples need processing\n")
                    self.copy_db_to_ramdisk()
                    # Run Kraken2
                    for f in files_to_run:
                        try:
                            self.run_kraken2(f)
                        except subprocess.CalledProcessError as e:
                            self.log("ERROR", f"ERROR processing {f}")
                            self.log("ERROR", f"{e}")
                            sys.exit(1)
                else:
                    self.log("INFO", "Kraken2 output files already exist")
                    self.log("INFO", "Not loading DB into memory")
                    self.log("INFO", "Skipping Kraken2\n")
                # Optional DB cleanup here

                # Run Braken
                for f in fastq_files:
                    try:
                        self.run_bracken(f)
                    except subprocess.CalledProcessError as e:
                        self.log("ERROR", f"ERROR processing {f}")
                        self.log("ERROR", f"{e}")
                        sys.exit(1)

                # Run kraken-biom
                try:
                    self.run_krakenbiom()
                except subprocess.CalledProcessError as e:
                    self.log("ERROR", f"ERROR running kraken-biom")
                    self.log("ERROR", f"{e}")
                    sys.exit(1)

                self.log("END", "MetNano is now complete. Bye!")
