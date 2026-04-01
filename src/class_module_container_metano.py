import os
import sys
import subprocess
import shutil
from datetime import datetime
import glob
import pandas as pd
import tempfile
import math
import re
from concurrent.futures import ThreadPoolExecutor, as_completed

from parameter_metano import *
#from parameter import *
#


#class for demulti and QC
class Demultiplexing_qc:

    def __init__(self, folder_no_sample_id, tread, metadata,readmap,qualityscore=10,minlength=1000,unzip=True   ):
        self.folder_no_sample_id = folder_no_sample_id.rstrip("/")
        
        self.last_folder = os.path.basename(self.folder_no_sample_id.rstrip("/"))
        print(f"Debug  self.last_folder = { self.last_folder}", flush=True)
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
       
 
    def run_unzip(self):
        # ---------- unzip and demultiplexe barcoded sample ----------

        if not os.path.exists(self.outputunzfilefastq):
            print(f"Folder unzip {self.outputunzfilefastq} is missing, running {unzip_barcode_script}", flush=True)
            print("Running unzip_barcode_script ...", flush=True)
            unz_cmd = f"{unzip_barcode_script} {self.folder_no_sample_id} {self.tread}"
            print(f"Command executed: {unz_cmd}", flush=True)
            subprocess.run(unz_cmd, shell=True, check=True)
            print("Unzip and demultiplexing done", flush=True)
            os.makedirs(self.outputunzfilefasta, exist_ok=True)
            unz_cmd = f"mv *.fastq  {self.outputunzfilefastq}"
            subprocess.run(unz_cmd, shell=True, check=True)
            print("move file to {self.outputunzfilefastq}", flush=True)
           
        else:
            print("Folder exists, skipping step", flush=True)

    def create_ln_file(self) :
       ### skip demultiplexing
        if not os.path.exists(self.outputunzfilefastq):
               print("Demultiplexing step not set with --unzip parameter. Direct copy ln link from original file to working directory", flush=True)
               os.makedirs(self.outputunzfilefastq, exist_ok=True)
               ln_cmd= f" ln -s {self.folder_no_sample_id }/*.fastq* {self.outputunzfilefastq}/ "
               print(f"Command executed : {ln_cmd} ", flush=True)
               subprocess.run(ln_cmd, shell=True, check=True)
               print("Symbolic link created  into workink directory ", flush=True)
            

        

    def run_rename(self):
        # ---------- Rename file with metadata provided ----------
        barcode_fastqs = glob.glob(os.path.join(self.outputunzfilefastq, "*barcode*.fastq*"))
        print(f"Check barcode file fastq in {self.outputunzfilefastq}", flush=True)
        if barcode_fastqs:
            print(f"Next step running {rename_barcodes_w_sheet} using metadata file provided {self.metadata}", flush=True)
            rename_cmd = f"{rename_barcodes_w_sheet} {self.metadata}"
            print(f"Command executed: {rename_cmd} in folder {self.outputunzfilefastq}", flush=True)
            subprocess.run(rename_cmd, shell=True, check=True, cwd=self.outputunzfilefastq)
            print("Rename with metadata done", flush=True)
        else:
            print("No barcodefile exists, skipping step run_rename", flush=True)
            

    def run_nanoplot_parallel(self):
         # ---------- Run individual nanoplot and mano QC for each conserved samples ----------
        print("Next step running nanoplot + nanoQC", flush=True)

        fastqs = glob.glob(os.path.join(self.outputunzfilefastq, "*.fastq"))
        print(f"[INFO] FASTQ files found: {len(fastqs)}", flush=True)

        if not fastqs:
            print("[WARN] No FASTQ files found, skipping NanoPlot/nanoQC", flush=True)
            return

        def run_one(fq):
            sample = os.path.basename(fq).replace(".fastq", "")
            nanoplot_out = os.path.join(self.outputunzfile,"01_nanoplot" ,f"nanoplot_{sample}")
            nanoQC_out   = os.path.join(self.outputunzfile,"02_nanoQC" ,f"nanoQC_{sample}")

            # ---------- NanoPlot ----------
            if os.path.isdir(nanoplot_out):
                print(f"[SKIP] NanoPlot already exists for {sample}", flush=True)
            else:
                nanoplot_cmd = (
                    f"{self.activate_commandnanoplot_env}   && NanoPlot --fastq {fq} "
                    f"--outdir {nanoplot_out} "
                    f"--prefix {sample} "
                    f"-t 1"
                )
                print(f"[RUN] {nanoplot_cmd}", flush=True)
                subprocess.run(nanoplot_cmd, shell=True, check=True)
                print(f"[END] NanoPlot {sample}", flush=True)

            # ---------- nanoQC ----------
            if os.path.isdir(nanoQC_out):
                print(f"[SKIP] nanoQC already exists for {sample}", flush=True)
            else:
                nanoQC_cmd = (
                    f"{self.activate_commandnanoQC_env} && "
                    f"nanoQC {fq} --outdir {nanoQC_out}"
                )
                print(f"[RUN] {nanoQC_cmd}", flush=True)
                subprocess.run(nanoQC_cmd, shell=True, check=True)
                print(f"[END] nanoQC {sample}", flush=True)

        # ---------- Parallel execution ----------
        with ThreadPoolExecutor(max_workers=self.tread) as exe:
            futures = [exe.submit(run_one, fq) for fq in fastqs]
            for f in as_completed(futures):
                f.result()  # surface exceptions



  

    def assembly_stat(self):
         # ---------- Generate assembly stats file for all sample in main folder ----------
        print(f"Next step running assembly stats", flush=True)
        out_dir = os.path.join(self.outputunzfile,"03_statistics" )
        out_file = os.path.join(out_dir ,"fastq_metrics.txt")
        os.makedirs(out_dir, exist_ok=True)
        if os.path.isfile(out_file):
            print("[SKIP] assembly_stats.txt already exists", flush=True)
            return
        cmd = f"{assembly_stats} -t  {self.outputunzfilefastq}/*.fastq > {out_file}  "
        print(f"Command executed: {cmd} in {self.outputunzfile} ", flush=True)
        subprocess.run(cmd, shell=True, check=True, cwd=self.outputunzfile)
      
    def run_chopper_parallel(self):
        # ---------- Run chopper and seqkit  with selected parameters (see parameter_metnano.txt file.) ----------        
        print("Next step running chopper", flush=True)
    
        filelist_fastQ = glob.glob(os.path.join(self.outputunzfilefastq, "*.fastq*"))
    
        # output directories
        self.trimmed_dir = os.path.join(self.outputunzfile, "04_trimmed")
        self.unmapped_dir = os.path.join(self.outputunzfile, "04a_unmapped")
        self.sliding_dir = os.path.join(self.outputunzfile, "05_sliding")
    
        os.makedirs(self.trimmed_dir, exist_ok=True)
        os.makedirs(self.sliding_dir, exist_ok=True)
    
        def process_fastq(fq):
            base_name = os.path.basename(fq)
    
            # ---- trimmed output ----
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
    
            # ---- sliding output ----
            out_fq_sliding = os.path.join(
                self.sliding_dir,
                base_name.replace(".fastq.gz", "")
                         .replace(".fastq", "")
                         + "_sliding.fastq"
            )
    
            # ---- STEP 1: chopper ----
            if not os.path.exists(out_fq):
                cmd = (
                    f"{self.activate_commandchopper_env} && "
                    f"cat {fq} | "
                    f"chopper --headcrop {head_trim} --tailcrop {tail_trim} -t 1 | "
                    f"chopper -q {self.qualityscore} -l {self.minlength} -t 1 "
                    f"> {out_fq}"
                )
                print(f"Running chopper: {cmd}", flush=True)
                subprocess.run(cmd, shell=True, check=True)
            else:
                print(f"Skipping chopper (exists): {out_fq}", flush=True)

                     # ---- STEP 1b: optional minimap2 + samtools if readmap is provided ----
            if hasattr(self, 'readmap') and self.readmap not in [None, "None"]:
                print(f"preparinf running mapping to reference", flush=True)
                unmapped_fq = os.path.join(
                    self.unmapped_dir,
                    base_name.replace(".fastq.gz", "")
                             .replace(".fastq", "")
                             + "_unmapped.fastq"
                )

                if not os.path.exists(unmapped_fq):
                    # Run minimap2 alignment
                    cmd_map = (
                        f" {self.activate_commandminimap_env} &&  minimap2 -a -x {mapping_preset} -t 1 "
                        f"{self.readmap} {out_fq} "
                        f"> {unmapped_fq.replace('_unmapped.fastq', '_alignment.sam')}"
                    )
                    print(f"Running minimap2: {cmd_map}", flush=True)
                    subprocess.run(cmd_map, shell=True, check=True)

                    # Extract unmapped reads
                    cmd_samtools = (
                        f"conda deactivate && conda activate samtools && "
                        f"samtools view -bf 4 "
                        f"{unmapped_fq.replace('_unmapped.fastq', '_alignment.sam')} "
                        f"> {unmapped_fq.replace('_unmapped.fastq', '_unmapped.bam')} && "
                        f"samtools fastq "
                        f"{unmapped_fq.replace('_unmapped.fastq', '_unmapped.bam')} "
                        f"> {unmapped_fq}"
                    )
                    print(f"Extracting unmapped reads: {cmd_samtools}", flush=True)
                    subprocess.run(cmd_samtools, shell=True, check=True)
                else:
                    print(f"Skipping readmap filtering (exists): {unmapped_fq}", flush=True)

                # Use unmapped reads for sliding step
                fq_for_sliding = unmapped_fq
            else:
                fq_for_sliding = out_fq
    
            # ---- STEP 2: seqkit sliding ----
            if os.path.exists(out_fq_sliding):
                print(f"Skipping sliding (exists): {out_fq_sliding}", flush=True)
                return
    
            print("Chopper done, starting seqkit sliding", flush=True)
    
            cmd = (
                f"module load seqkit && "
                f"seqkit sliding "
                f"-W {self.minlength} "
                f"-s {self.minlength} "
                f"{fq_for_sliding} "
                f"> {out_fq_sliding}"
            )
    
            print(f"Running seqkit: {cmd}", flush=True)
            subprocess.run(cmd, shell=True, check=True)


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
            print("Skipping unzip")  # runs if False
        else:
            self.create_ln_file()
        self.run_rename()
        self.run_nanoplot_parallel()
        self.assembly_stat()
        self.run_chopper_parallel()
        print(f"Demultiplexing_qc is Done", flush=True)

#####################################################################################################

class Kraken2Runner:
        """
        Python class to run Kraken2 on multiple FASTA/FASTQ files,
        using a temporary DB copy in RAM (/dev/shm) for speed.
        """
    
        def __init__(self, db_path, db_name, fasta_dir, out_dir,kraken_env, braken_env,
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
            self.base_braken = os.path.join(self.out_dir, "07_braken", self.db_name)
            self.base_biom = os.path.join(self.out_dir, "08_biom", self.db_name)
          
            os.makedirs(self.base_kraken, exist_ok=True)
            os.makedirs(self.base_braken, exist_ok=True)
            os.makedirs(self.base_biom, exist_ok=True)
           
            self.activate_commandkraken_env= f"{CONDA_SOURCE_private} && conda activate {kraken_env}"
            self.activate_commandbraken_env= f"{CONDA_SOURCE_private} && conda activate {braken_env}"
            self.activate_commandbrabiom_env = f"{CONDA_SOURCE_private} && conda activate {kraken_biom_env}"
    
            os.makedirs(self.out_dir, exist_ok=True)
    
        def log(self, message):
            print(f"{datetime.now()} - INFO - {message}", flush=True)
    
        def copy_db_to_ramdisk(self):
         # ---------- load db in ram disk ----------        
            """Copy the Kraken2 database to /dev/shm"""
            self.log(f"Copying database {self.db_name} to RAM (/dev/shm)...")
            os.makedirs(self.temp_db_path, exist_ok=True)
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

            self.log("Database copy completed.")

        def concatenate_bracken(self, bracken_files, fasta_file):
             # ----------set environnement variables for kraken and braken ----------        
                expected_cols = [
                    "name",
                    "taxonomy_id",
                    "taxonomy_lvl",
                    "kraken_assigned_reads",
                    "added_reads",
                    "new_est_reads",
                    "fraction_total_reads"
                ]
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
                    self.out_dir,"07_braken" ,f"{self.db_name}",
                    f"{os.path.basename(fasta_file)}_{self.confidence}.bracken.ALL.tsv"
                ) 
                final_df.to_csv(out_concat, sep="\t", index=False)
                return out_concat    
    
        def run_kraken2(self, fasta_file):
        # ----------import panda command ----------     
            import pandas as pd
        # ----------Kraken2 command ----------        
           # """Run Kraken2 on a single FASTA/FASTQ file"""
            sample_name = os.path.splitext(os.path.basename(fasta_file))[0]
            print(f"BEFORE: {sample_name}")
            sample_name = re.sub(r"_sliding(_[0-9.]+)?", "", sample_name, flags=re.IGNORECASE)
            print(f"AFTER: {sample_name}")
            
            prefix = f"{sample_name}_{self.confidence}"
     
           
            
            kraken_output = os.path.join(self.base_kraken, f"{prefix}.all")
            kraken_report = os.path.join(self.base_kraken, f"{prefix}.report")
            
            bracken_output = os.path.join(self.base_braken, f"{prefix}.bracken")
            bracken_report = os.path.join(self.base_braken, f"{prefix}.bracken.report")

    
            
            def file_exists_and_not_empty(path):
                return os.path.exists(path) and os.path.getsize(path) > 0  
            
           # ---- Skip logic ----
            if file_exists_and_not_empty(kraken_output) and file_exists_and_not_empty(kraken_report):
                print(f"[INFO] Kraken2 output already exists for {prefix}. Skipping.")
            else:
          
                # Vérifie que les fichiers essentiels existent
                required_files = ["hash.k2d", "opts.k2d", "taxo.k2d"]
                missing_files = [f for f in required_files if not os.path.exists(os.path.join(self.temp_db_path, f))]
                if missing_files:
                    raise RuntimeError(f"[ERROR] Missing required DB files in RAM: {missing_files}")
            
                # ---- Kraken2 ----
                kraken_cmd = (
                    
                    #f"{self.activate_commandkraken_env} && "
                    f" module load kraken && "
                    f"kraken2 "
                    f"--db {self.temp_db_path} "
                    f"--threads 1 "
                    f"--output {kraken_output} "
                    f"--report {kraken_report} "
                    f"--confidence {self.confidence} "
                    f"--memory-mapping "
                    f"{fasta_file} "
                )
            
                print(f"[INFO] Running Kraken2 for {prefix}")
                subprocess.run(kraken_cmd, shell=True, check=True)
                print(f"[INFO] Kraken2 ended, start braken")

                # ---- Post-run check ----
                if not file_exists_and_not_empty(kraken_output) or not file_exists_and_not_empty(kraken_report):
                    raise RuntimeError(f"[ERROR] Kraken2 run did not produce expected output files for {prefix}")

            # ---- Bracken (all taxonomic levels) ----
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
                    self.base_braken,
                    f"{prefix}.bracken.{level_name}.tsv"
                )

                bracken_report_s = os.path.join(
                    self.base_braken,
                    f"{prefix}.bracken.{level_name}.report"
                )
            
                # Skip if already exists
                if os.path.exists(bracken_output_s) and os.path.getsize(bracken_output_s) > 0:
                    self.log(f"looking in {bracken_output} if file exist")
                    self.log(f"Skipping Bracken {level_name} — file exists.")
                    bracken_outputs.append(bracken_output)
                    continue
            
                bracken_cmd = (
                    f"{self.activate_commandbraken_env} && bracken "
                    f"-d {self.db_path} "
                    f"-i {kraken_report} "
                    f"-o {bracken_output_s} "
                    f"-w {bracken_report_s} "  
                    f"-r {self.minlength} "
                    f"-l {level_code}"
                )
            
                self.log(f"Running Bracken ({level_name})...")
                self.log(bracken_cmd)
                subprocess.run(bracken_cmd, shell=True, check=True)          
                bracken_outputs.append(bracken_output_s)          
            self.log("Bracken completed for all taxonomic levels.")
            concat_file = self.concatenate_bracken(bracken_outputs, fasta_file)
            self.log(f"Final Bracken table written to {concat_file}")
        
    
        def run_all(self, max_workers=1):
            
                fasta_files = glob.glob(os.path.join(self.fasta_dir, "*.fastq*"))
            
                if not fasta_files:
                    self.log("No FASTQ files found.")
                    return
            
                files_to_run = []
            
                for f in fasta_files:
                    sample_name = os.path.splitext(os.path.basename(f))[0]
                    print(f"BEFORE: {sample_name}")
                    sample_name = re.sub(r"_sliding(_[0-9.]+)?", "", sample_name, flags=re.IGNORECASE)
                    print(f"AFTER: {sample_name}")
                    expected_report = os.path.join(self.base_kraken, f"{sample_name}.report")
            
                    if not os.path.exists(expected_report) or os.path.getsize(expected_report) == 0:
                        files_to_run.append(f)
            
                # Only load RAM + run Kraken if needed
                if files_to_run:
                    self.log(f"{len(files_to_run)} samples need processing. Loading DB into RAM...")
                    self.copy_db_to_ramdisk()
            
                    for f in files_to_run:
                        try:
                            self.run_kraken2(f)
                        except subprocess.CalledProcessError as e:
                            self.log(f"ERROR processing {f}: {e}")
            
                    self.log("Kraken processing finished.")
                    # optional cleanup here
            
                else:
                    self.log("All samples already processed. Skipping Kraken2.")

                # ---------- Generate BIOM tables: one per report ----------
                report_files = glob.glob(os.path.join(self.base_braken, "*.report"))
                braken_reports = report_files
                if not braken_reports:
                    self.log("No Braken reports found. Cannot generate BIOM tables.")
                    return

                OUTPUT_DIR_biom = os.path.join(self.out_dir, "08_biom","braken",f"{self.db_name}")
                os.makedirs(OUTPUT_DIR_biom, exist_ok=True)

                for report in braken_reports:
                    sample_name = os.path.splitext(os.path.basename(report))[0]
                    output_file = os.path.join(OUTPUT_DIR_biom, f"{sample_name}.biom")

                    # Skip if BIOM already exists
                    if os.path.exists(output_file) and os.path.getsize(output_file) > 0:
                        self.log(f"BIOM for {sample_name} already exists. Skipping.")
                        continue

                    # Build command for this single file
                    cmd_biom = (
                        f"{self.activate_commandbrabiom_env} && "
                        f"kraken-biom {report} "
                        f"--fmt json "
                        f"--output_fp {output_file}"
                        f" --max D  --min S"
                    )

                    self.log(f"Generating BIOM Braken for {sample_name}...")
                    subprocess.run(cmd_biom, shell=True, check=True)

                self.log("All BIOM tables generated. Bye!")

             
                    





       
