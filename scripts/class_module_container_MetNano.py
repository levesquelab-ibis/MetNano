import os
import sys
import subprocess
import shutil
import glob
import pandas as pd
import tempfile
import math
import re
from datetime import datetime
from concurrent.futures import ThreadPoolExecutor, as_completed
from parameter_MetNano import *

############### Class for preprocessing of raw ONT reads ###############
class Demultiplexing_qc:

    def __init__(self, folder_no_sample_id, tread, metadata, readmap, qualityscore=10, minlength=1000, unzip=True):
        self.folder_no_sample_id = folder_no_sample_id.rstrip("/")
        self.last_folder = os.path.basename(self.folder_no_sample_id.rstrip("/"))
        self.qualityscore=qualityscore
        self.minlength=minlength
        self.tread = tread
        self.readmap = readmap
        self.metadata = metadata
        self.nanoQC_env = nanoQC_env
        self.unzip=unzip
        self.current_dir = os.path.abspath(".")
        self.date = datetime.now().strftime('%Y-%m-%d')
        self.activate_commandnanoQC_env= f"{CONDA_SOURCE_private} && conda activate {nanoQC_env}"
        self.activate_commandchopper_env= f"{CONDA_SOURCE_private} && conda activate {chopper_env}"
        self.outputunzfile = f"{self.current_dir}"
        self.outputunzfilefastq= os.path.join(self.outputunzfile, "00_fastq_demultiplexed")
        self.activate_commandnanoplot_env = f"{CONDA_SOURCE_private} && conda activate {nanoplot_env}"
        self.activate_commandminimap_env = f"{CONDA_SOURCE_private} && conda activate {minimap2_env}"

    def log(self, type, message):
            print(f"{datetime.now().isoformat(sep=' ', timespec='seconds')} [{type}] {message}", flush=True)

    def run_unzip(self):
        # ---------- Unzip and demultiplex barcoded samples ----------
        if not os.path.exists(self.outputunzfilefastq):
            self.log("INFO", "Folder unzip {self.outputunzfilefastq} is missing, running {unzip_barcode_script}")
            self.log("INFO", "Running unzip_barcode_script ...")
            unz_cmd = f"{unzip_barcode_script} {self.folder_no_sample_id} {self.tread}"
            self.log("RUN", "{unz_cmd}")
            subprocess.run(unz_cmd, shell=True, check=True)
            self.log("INFO", "Unzip and demultiplexing done")
            os.makedirs(self.outputunzfilefasta, exist_ok=True)
            unz_cmd = f"mv *.fastq  {self.outputunzfilefastq}"
            subprocess.run(unz_cmd, shell=True, check=True)
            self.log("INFO", "move file to {self.outputunzfilefastq}")
        else:
            self.log("INFO", "Folder exists, skipping step")

    def create_ln_file(self) :
        # ---------- Skip demultiplexing ----------
        print("# Demultiplexing", flush=True)
        if not os.path.exists(self.outputunzfilefastq):
            self.log("INFO", "Demultiplexing step not set with --unzip parameter.")
            self.log("INFO Direct copy link from original file to working directory")
            os.makedirs(self.outputunzfilefastq, exist_ok=True)
            ln_cmd= f"ln -s {self.folder_no_sample_id}/*.fastq* {self.outputunzfilefastq}/ "
            self.log("RUN", f"{ln_cmd} ")
            subprocess.run(ln_cmd, shell=True, check=True)
            self.log("INFO", "Symbolic link created  into working directory\n")
        elif os.path.exists(self.outputunzfilefastq):
            self.log("INFO", "Symbolic link already exists\n")

    def run_rename(self):
        # ---------- Rename file with metadata provided ----------
        barcode_fastqs = glob.glob(os.path.join(self.outputunzfilefastq, "*barcode*.fastq*"))
        print(f"# Renaming barcodes", flush=True)
        if barcode_fastqs:
            rename_cmd = f"{rename_barcodes_w_sheet} {self.metadata}"
            self.log("RUN", f"{rename_cmd}")
            subprocess.run(rename_cmd, shell=True, check=True, cwd=self.outputunzfilefastq)
            self.log("INFO", "Renaming samples with metadata done\n")
        else:
            self.log("INFO", "No barcode file exists")

    def run_nanoplot_parallel(self):
        # ---------- Run individual NanoPlot and nanoQC for each samples ----------
        fastqs = glob.glob(os.path.join(self.outputunzfilefastq, "*.fastq"))
        self.log("INFO", f"FASTQ file found: {len(fastqs)}\n")

        if not fastqs:
            self.log("WARN", "No FASTQ files found, skipping NanoPlot/nanoQC\n")
            return

        def run_one(fq):
            sample = os.path.basename(fq).replace(".fastq", "")
            nanoplot_out = os.path.join(self.outputunzfile,"01_nanoplot" ,f"nanoplot_{sample}")
            nanoQC_out = os.path.join(self.outputunzfile,"02_nanoQC" ,f"nanoQC_{sample}")

            # ---------- NanoPlot ----------
            print("# Quality control with NanoPlot", flush = True)
            self.log("INFO", f"Analyzing sample {sample}...")
            if os.path.isdir(nanoplot_out):
                self.log("INFO", f"NanoPlot already exists for sample {sample}\n")
            else:
                nanoplot_cmd = (
                    f"{self.activate_commandnanoplot_env}   && NanoPlot --fastq {fq} "
                    f"--outdir {nanoplot_out} "
                    f"--prefix {sample}_ "
                    f"-t 1"
                )
                self.log("RUN", "{nanoplot_cmd}")
                subprocess.run(nanoplot_cmd, shell=True, check=True)
                self.log("INFO", f"NanoPlot completed for sample {sample}")

            # ---------- nanoQC ----------
            print("# Quality control with nanoQC", flush = True)
            self.log("INFO", f"Analyzing sample {sample}...")
            if os.path.isdir(nanoQC_out):
                self.log("INFO", f"nanoQC already exists for sample {sample}\n")
            else:
                nanoQC_cmd = (
                    f"{self.activate_commandnanoQC_env} && "
                    f"nanoQC {fq} --outdir {nanoQC_out}"
                )
                self.log("RUN", f"{nanoQC_cmd}")
                subprocess.run(nanoQC_cmd, shell=True, check=True)
                self.log("INFO", f"nanoQC completed for sample {sample}\n")

        # ---------- Parallel execution ----------
        with ThreadPoolExecutor(max_workers=self.tread) as exe:
            futures = [exe.submit(run_one, fq) for fq in fastqs]
            for f in as_completed(futures):
                f.result()  # surface exceptions

    def assembly_stat(self):
        # ---------- Generate assembly stats files for all samples in ONT folder ----------
        print(f"# Sequencing statistics with assembly-stats", flush=True)
        self.log("INFO", f"Computing sequencing reads statistics...")
        out_dir = os.path.join(self.outputunzfile,"03_statistics" )
        out_file = os.path.join(out_dir ,"fastq_metrics.txt")
        os.makedirs(out_dir, exist_ok=True)
        if os.path.isfile(out_file):
            self.log("INFO", "assembly_stats.txt already exists\n")
            return
        cmd = f"{assembly_stats} -t {self.outputunzfilefastq}/*.fastq > {out_file}  "
        self.log("RUN", f" {cmd} in {self.outputunzfile}")
        subprocess.run(cmd, shell=True, check=True, cwd=self.outputunzfile)
        self.log("INFO", f"Sequencing stats written to {out_file}\n")

    def run_chopper_parallel(self):
        # ---------- Run chopper and seqkit  with selected parameters (see parameter_MetNano.py file) ----------
        print("# Read filtering and trimming with chopper", flush=True)
        filelist_fastQ = glob.glob(os.path.join(self.outputunzfilefastq, "*.fastq*"))

        # Output directories
        self.trimmed_dir = os.path.join(self.outputunzfile, "04_trimmed")
        self.unmapped_dir = os.path.join(self.outputunzfile, "04a_unmapped")
        self.sliding_dir = os.path.join(self.outputunzfile, "05_sliding")
        os.makedirs(self.trimmed_dir, exist_ok=True)
        os.makedirs(self.sliding_dir, exist_ok=True)

        def process_fastq(fq):
            base_name = os.path.basename(fq)

            # Trimmed output
            if fq.endswith(".fastq.gz"):
                out_fq = os.path.join(
                    self.trimmed_dir,
                    base_name.replace(".fastq.gz", "_trimmed.fastq")
                )
            else:
                out_fq = os.path.join(
                    self.trimmed_dir,
                    base_name.replace(".fastq", "_trimmed.fastq")
                )

            # Segmented output
            out_fq_sliding = os.path.join(
                self.sliding_dir,
                base_name.replace(".fastq.gz", "")
                         .replace(".fastq", "")
                         + "_sliding.fastq"
            )

            # ---- STEP 1: chopper ----
            self.log("INFO", f"Trimming file {base_name}...")
            if not os.path.exists(out_fq):
                cmd = (
                    f"{self.activate_commandchopper_env} && "
                    f"cat {fq} | "
                    f"chopper --headcrop {head_trim} --tailcrop {tail_trim} -t 1 | "
                    f"chopper -q {self.qualityscore} -l {self.minlength} -t 1 "
                    f"> {out_fq}"
                )
                self.log("RUN", f"{cmd}")
                subprocess.run(cmd, shell=True, check=True)
                self.log("INFO", f"Trimming completed for {base_name}\n")
            else:
                self.log("INFO", f"Trimmed file already exists for {base_name}\n")

            # ----- STEP 1b: optional minimap2 + samtools if readmap parameter is provided -----
            if hasattr(self, 'readmap') and self.readmap not in [None, "None"]:
                unmapped_fq = os.path.join(
                    self.unmapped_dir,
                    base_name.replace(".fastq.gz", "")
                             .replace(".fastq", "")
                             + "_unmapped.fastq"
                )

                if not os.path.exists(unmapped_fq):
                    # ----- Run minimap2 alignment -----
                    print(f"# Mapping reads to reference genome with minimap2", flush = True)
                    cmd_map = (
                        f" {self.activate_commandminimap_env} &&  minimap2 -a -x {mapping_preset} -t 1 "
                        f"{self.readmap} {out_fq} "
                        f"> {unmapped_fq.replace('_unmapped.fastq', '_alignment.sam')}"
                    )
                    self.log("RUN", f"{cmd_map}")
                    subprocess.run(cmd_map, shell=True, check=True)
                    self.log("INFO", f"Mapping completed for {out_fq}\n")

                    # ----- Extract unmapped reads -----
                    print(f"# Extraction of unmapped reads with samtools", flush = True)
                    cmd_samtools = (
                        f"conda deactivate && conda activate samtools && "
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
                    self.log("INFO", f"Readmap output already exists: {unmapped_fq}\n")

                # Use unmapped reads for segmentation step
                fq_for_sliding = unmapped_fq
            else:
                fq_for_sliding = out_fq

            # ----- STEP 2: SeqKit -----
            print("# Read segmentation with SeqKit", flush=True)
            self.log("INFO", f"Segmenting reads for {fq_for_sliding}...")
            if os.path.exists(out_fq_sliding):
                self.log("INFO", f"Segmented file already exists for {fq_for_sliding}\n")
                return

            cmd = (
                f"module load seqkit && "
                f"seqkit sliding "
                f"-W {self.minlength} "
                f"-s {self.minlength} "
                f"{fq_for_sliding} "
                f"> {out_fq_sliding}"
            )
            self.log("RUN", f"{cmd}")
            subprocess.run(cmd, shell=True, check=True)
            self.log("INFO", f"Segmenting completed for {fq_for_sliding}\n")

        with ThreadPoolExecutor(max_workers=self.tread) as exe:
                      future_map = {
                           exe.submit(process_fastq, fq): fq
                           for fq in filelist_fastQ
                          }
                      for future in as_completed(future_map):
                           fq = future_map[future]
                           future.result()  # crash fast, loud, and clear

    def run(self):
        if not self.unzip:
            self.run_unzip()
            self.log("INFO", "Skipping unzip")  # runs if False
        else:
            self.create_ln_file()
        self.run_rename()
        self.run_nanoplot_parallel()
        self.assembly_stat()
        self.run_chopper_parallel()

########## Class for read taxonomic classification and abundance estimation ##########
class Kraken2Runner:
        """
        Python class to run Kraken2 on multiple FASTA/FASTQ files,
        using a temporary DB copy in RAM (/dev/shm) for speed.
        """
        def __init__(self, db_path, db_name, fasta_dir, out_dir,kraken_env, bracken_env,
                     threads=6, confidence=0.01 ,minlength=300):
            self.db_path = db_path
            self.db_name = db_name
            self.fasta_dir = fasta_dir
            self.out_dir = out_dir
            self.threads = threads
            self.confidence = confidence
            self.minlength=minlength
            self.temp_db_path = os.path.join("/dev/shm", self.db_name)
            self.base_kraken = os.path.join(self.out_dir, "06_kraken2", self.db_name)
            self.base_bracken = os.path.join(self.out_dir, "07_bracken", self.db_name)
            self.base_biom = os.path.join(self.out_dir, "08_biom", self.db_name)

            os.makedirs(self.base_kraken, exist_ok=True)
            os.makedirs(self.base_bracken, exist_ok=True)
            os.makedirs(self.base_biom, exist_ok=True)

            self.activate_commandkraken_env= f"{CONDA_SOURCE_private} && conda activate {kraken_env}"
            self.activate_commandbracken_env= f"{CONDA_SOURCE_private} && conda activate {bracken_env}"
            self.activate_commandbrabiom_env = f"{CONDA_SOURCE_private} && conda activate {kraken_biom_env}"

            os.makedirs(self.out_dir, exist_ok=True)

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

        def concatenate_bracken(self, bracken_files, fasta_file):
            # ---------- Set environment variables for kraken2 and bracken ----------
            expected_cols = [
                "name",
                "taxonomy_id",
                "taxonomy_lvl",
                "kraken_assigned_reads",
                "added_reads",
                "new_est_reads",
                "fraction_total_reads"]
            dfs = []

            for f in bracken_files:
                level = os.path.basename(f).split(".bracken.")[1].replace(".tsv", "")
                df = pd.read_csv(f, sep="\t")
                # Ensure consistent columns
                for col in expected_cols:
                    if col not in df.columns:
                        df[col] = 0
                df = df[expected_cols]
                df["sample"] = os.path.basename(fasta_file)
                df["level"] = level
                dfs.append(df)
            final_df = pd.concat(dfs, ignore_index=True)
            out_concat = os.path.join(
                self.out_dir,"07_bracken" ,f"{self.db_name}",
                f"{os.path.basename(fasta_file)}_{self.confidence}.bracken.ALL.tsv"
            )
            final_df.to_csv(out_concat, sep="\t", index=False)
            return out_concat

        def run_kraken2(self, fasta_file):
            import pandas as pd
            # --------- Kraken2 command ----------
            """Run Kraken2 on a single FASTA/FASTQ file"""
            print("# Read taxonomic classification with Kraken2", flush = True)
            sample_name = os.path.splitext(os.path.basename(fasta_file))[0]
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
                    #f"{self.activate_commandkraken_env} && "
                    f"module load kraken && "
                    f"kraken2 "
                    f"--db {self.temp_db_path} "
                    f"--threads 1 "
                    f"--output {kraken_output} "
                    f"--report {kraken_report} "
                    f"--confidence {self.confidence} "
                    f"--memory-mapping "
                    f"{fasta_file}"
                )

                self.log("INFO", f"Running Kraken2 for {prefix}")
                subprocess.run(kraken_cmd, shell=True, check=True)
                self.log("END", f"Kraken2 ended\n")

                # Post-run check
                if not file_exists_and_not_empty(kraken_output) or not file_exists_and_not_empty(kraken_report):
                    raise RuntimeError(f"[ERROR] Kraken2 run did not produce expected output files for {prefix}")

            # ----- Bracken (all taxonomic levels) -----
            print("# Taxon abundance estimation with Bracken", flush = True)
            BRACKEN_LEVELS = {
                "species": "S",
                "genus": "G",
                "family": "F",
                "order": "O",
                "class": "C",
                "phylum": "P"
            }

            bracken_outputs = []
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
                    self.log("INFO", f"Bracken {level_name} file already exists")
                    bracken_outputs.append(bracken_output)
                    continue

                bracken_cmd = (
                    f"{self.activate_commandbracken_env} && bracken "
                    f"-d {self.db_path} "
                    f"-i {kraken_report} "
                    f"-o {bracken_output_s} "
                    f"-w {bracken_report_s} "
                    f"-r {self.minlength} "
                    f"-l {level_code}"
                )
                self.log("INFO", f"Running Bracken ({level_name})...")
                self.log("RUN", bracken_cmd)
                subprocess.run(bracken_cmd, shell=True, check=True)
                bracken_outputs.append(bracken_output_s)
                print()
            self.log("INFO", "Bracken completed for all taxonomic levels")
            concat_file = self.concatenate_bracken(bracken_outputs, fasta_file)
            self.log("INFO", f"Final Bracken table written to {concat_file}\n")

        def run_all(self, max_workers=1):

                fasta_files = glob.glob(os.path.join(self.fasta_dir, "*.fastq*"))

                if not fasta_files:
                    self.log("WARN", "No FASTQ files found.")
                    return

                files_to_run = []

                for f in fasta_files:
                    sample_name = os.path.splitext(os.path.basename(f))[0]
                    sample_name = re.sub(r"_sliding(_[0-9.]+)?", "", sample_name, flags=re.IGNORECASE)
                    expected_report = os.path.join(self.base_kraken, f"{sample_name}.report")

                    if not os.path.exists(expected_report) or os.path.getsize(expected_report) == 0:
                        files_to_run.append(f)

                # Only load RAM + run Kraken2 if needed
                if files_to_run:
                    self.log("INFO", f"{len(files_to_run)} samples need processing\n")
                    self.copy_db_to_ramdisk()

                    for f in files_to_run:
                        try:
                            self.run_kraken2(f)
                        except subprocess.CalledProcessError as e:
                            self.log("ERROR", f"ERROR processing {f}: {e}")

                    self.log("INFO", "Kraken2 processing finished.\n")

                    # Optional cleanup here

                else:
                    self.log("INFO", "All samples already processed. Skipping Kraken2.\n")

                # ---------- Generate BIOM tables: one per report ----------
                print("# Converting Kraken2-style reports to BIOM format with kraken-biom", flush = True)
                report_files = glob.glob(os.path.join(self.base_bracken, "*.report"))
                bracken_reports = report_files
                if not bracken_reports:
                    self.log("WARN", "No Bracken reports found. Cannot generate BIOM tables")
                    return

                OUTPUT_DIR_biom = os.path.join(self.out_dir, "08_biom","bracken",f"{self.db_name}")
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
                        f"{self.activate_commandbrabiom_env} && "
                        f"kraken-biom {report} "
                        f"--fmt json "
                        f"--output_fp {output_file} "
                        f"--max D --min S"
                    )

                    self.log("INFO", f"Generating BIOM from Kraken2-style reports for {sample_name}...")
                    subprocess.run(cmd_biom, shell=True, check=True)

                self.log("INFO", "All BIOM tables are generated\n")
                self.log("END", "MetNano is now completed. Bye!")
