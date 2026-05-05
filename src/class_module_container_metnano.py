import os
import sys
import glob
import subprocess
from datetime import datetime
from concurrent.futures import ThreadPoolExecutor, as_completed
from parameter_metnano import *

class Preprocessing:
    """
    Python class allowing the computation of sequencing statistics, quality
    control, and filtering from from raw reads, followed by read
    segmentation to homogenous length. Optionnaly, filtered reads can also
    be mapped onto a reference genome to deplete host reads.
    """
    def __init__(self, ont_folder_path, threads, metadata, readmap, qualityscore, minlength, unzip, rename):
        self.ont_folder_path = ont_folder_path.rstrip("/")
        self.ont_folder = os.path.basename(self.ont_folder_path.rstrip("/"))
        self.qualityscore = qualityscore
        self.minlength = minlength
        self.parameters = f"min_{minlength}bp_{qualityscore}Q"
        self.threads = int(threads)
        self.ref = readmap
        self.metadata = metadata
        self.unzip  = unzip
        self.rename = rename
        self.output_folder = os.path.abspath(".")
        self.output_folder_fastq= os.path.join(self.output_folder, "00-raw_FASTQs")
        self.chopper_dir = os.path.join(self.output_folder, "04-chopper", self.parameters)
        self.minimap2_dir = os.path.join(self.output_folder, "04b-minimap2", self.parameters)
        self.seqkit_dir = os.path.join(self.output_folder, "05-SeqKit", self.parameters)

    def log(self, kind, message):
        """Print formatted log message"""
        print(f"{datetime.now().isoformat(sep=' ', timespec='seconds')} [{kind}] {message}", flush=True)

    def parallel_execution(self, function, fastqs):
        """Parallel execution"""
        with ThreadPoolExecutor(max_workers = self.max_workers) as exe:
            future_map = {
               exe.submit(function, fq): fq
               for fq in fastqs
              }
            for future in as_completed(future_map):
                fq = future_map[future]
                future.result()  # crash fast, loud, and clear

    def unzip_fastq(self):
        """Unzip and concatenate archived FASTQ files"""
        if not os.path.exists(self.output_folder_fastq) or not os.listdir(self.output_folder_fastq):
            os.makedirs(self.output_folder_fastq, exist_ok=True)
            self.log("INFO", "Unzipping FASTQ files...")
            unzip_cmd = f"{unzip_barcode_script} {self.ont_folder_path} {self.threads}"
            self.log("RUN", f"{unzip_cmd}")
            subprocess.run(unzip_cmd, shell=True, check=True)
            self.log("INFO", "Unzip and demultiplexing done")
            mv_cmd = f"mv {self.output_folder_fastq}/../*_unzipped/*.fastq {self.output_folder_fastq}"
            subprocess.run(mv_cmd, shell=True, check=True)
            rm_cmd = f"rmdir {self.output_folder_fastq}/../*_unzipped/"
            subprocess.run(rm_cmd, shell=True, check=True)
            self.log("INFO", f"FASTQ files are in {self.output_folder_fastq}")
        else:
            self.log("INFO", "FASTQ folder exists, skipping unzipping")

    def create_ln_file(self) :
        """Create symbolic link toward raw FASTQs"""
        if not os.path.exists(self.output_folder_fastq) or not os.listdir(self.output_folder_fastq):
            self.log("INFO", "Linking original files to working directory...")
            os.makedirs(self.output_folder_fastq, exist_ok=True)
            ln_cmd= f"ln -s {self.ont_folder_path}/*.fastq {self.output_folder_fastq}/ "
            self.log("RUN", f"{ln_cmd} ")
            subprocess.run(ln_cmd, shell=True, check=True)
            self.log("INFO", "Symbolic link created into working directory")
        elif os.path.exists(self.output_folder_fastq):
            self.log("INFO", "Symbolic links already exist")

    def rename_fastq(self):
        """Rename file with provided association file"""
        barcode_fastqs = glob.glob(os.path.join(self.output_folder_fastq, "*barcode*.fastq*"))
        fastqs = glob.glob(os.path.join(self.output_folder_fastq, "*.fastq"))
        if barcode_fastqs:
            self.log("INFO", "Renaming samples with metadata...")
            rename_cmd = f"{rename_barcodes_w_sheet} {self.metadata} {self.output_folder_fastq}"
            self.log("RUN", f"{rename_cmd}")
            subprocess.run(rename_cmd, shell=True, check=True)
            self.log("INFO", "Renaming samples with metadata done")
        elif fastqs and not barcode_fastqs:
            self.log("INFO", "FASTQ files already renamed")
        else:
            self.log("ERROR", "No barcode file exists")
            self.log("ERROR", "Check your FASTQ directory")
            self.log("END", "Exiting MetNano")
            sys.exit(1)

    def assembly_stats(self):
        """Generate assembly stats files for all samples in ONT folder"""
        assemblystats_dir = os.path.join(self.output_folder,"01-assembly-stats" )
        out_file = os.path.join(assemblystats_dir ,"assembly-stats.tsv")
        os.makedirs(assemblystats_dir, exist_ok=True)
        if not os.path.isfile(out_file) or os.path.getsize(out_file) == 0:
            self.log("INFO", "Computing sequencing reads statistics...")
            assemblystats_cmd = f"assembly-stats -t {self.output_folder_fastq}/*.fastq > {out_file}  "
            self.log("RUN", f"{assemblystats_cmd}")
            subprocess.run(assemblystats_cmd, shell=True, check=True)
            rename_cmd = f"sed -i 's|{self.output_folder_fastq}/||' {out_file}"
            subprocess.run(rename_cmd, shell=True, check=True)
            self.log("INFO", "Sequencing statistics computed for all FASTQs")
        else:
            self.log("INFO", "File assembly-stats.tsv already exists")

    def NanoPlot_parallel(self, fastqs):
        """Quality Control with NanoPlot"""
        nanoplot_dir = os.path.join(self.output_folder, "02-NanoPlot")
        os.makedirs(nanoplot_dir, exist_ok=True)

        def NanoPlot_single(fq):
            sample = os.path.basename(fq).replace(".fastq", "")
            nanoplot_out = os.path.join(nanoplot_dir, f"{sample}")
            nanoplot_html = os.path.join(nanoplot_out, f"{sample}_NanoPlot-report.html")
            nanoplot_log = os.path.join(nanoplot_out, f"{sample}.log")
            os.makedirs(nanoplot_out, exist_ok=True)
            if not os.path.isfile(nanoplot_html) or os.path.getsize(nanoplot_html) == 0:
                self.log("INFO", f"Analyzing sample {sample}...")
                nanoplot_cmd = (
                    f"NanoPlot --fastq {fq} "
                    f"--outdir {nanoplot_out} "
                    f"--prefix {sample}_ "
                    f"-t {self.threads_per_task}"
                )
                self.log("RUN", f"{nanoplot_cmd}")
                with open(nanoplot_log, "w") as log:
                    subprocess.run(
                        nanoplot_cmd,
                        shell=True,
                        check=True,
                        stdout=log,
                        stderr=log
                    )
                    self.log("INFO", f"NanoPlot completed for sample {sample}")
            else:
                self.log("INFO", f"NanoPlot output already exists for sample {sample}")

        # ---------- Run in parallel ----------
        self.parallel_execution(NanoPlot_single, fastqs)

    def nanoQC_parallel(self, fastqs):
        """Quality Control with nanoQC"""
        nanoQC_dir = os.path.join(self.output_folder, "03-nanoQC")
        os.makedirs(nanoQC_dir, exist_ok=True)

        def nanoQC_single(fq):
            sample = os.path.basename(fq).replace(".fastq", "")
            nanoQC_out = os.path.join(nanoQC_dir, f"{sample}")
            nanoQC_html = os.path.join(nanoQC_out, "nanoQC.html")
            os.makedirs(nanoQC_out, exist_ok=True)
            if not os.path.isfile(nanoQC_html) or os.path.getsize(nanoQC_html) == 0:
                self.log("INFO", f"Analyzing sample {sample}...")
                nanoQC_cmd = (
                    f"nanoQC {fq} --outdir {nanoQC_out}"
                )
                self.log("RUN", f"{nanoQC_cmd}")
                subprocess.run(nanoQC_cmd, shell=True, check=True)
                self.log("INFO", f"nanoQC completed for sample {sample}")
            else:
                self.log("INFO", f"nanoQC output already exists for sample {sample}")

        # ---------- Run in parallel ----------
        self.parallel_execution(nanoQC_single, fastqs)

    def chopper_parallel(self, fastqs):
        """Run chopper with selected parameters (see parameter_MetNano.py file)"""
        os.makedirs(self.chopper_dir, exist_ok=True)

        def chopper_single(fq):
            sample = os.path.basename(fq).replace(".fastq.gz", "").replace(".fastq", "")
            trimmed_fastq = os.path.join(self.chopper_dir, sample + "_trimmed.fastq")
            chopper_log = os.path.join(self.chopper_dir, f"{sample}.log")
            if not os.path.isfile(trimmed_fastq) or os.path.getsize(trimmed_fastq) == 0:
                self.log("INFO", f"Trimming file {sample}...")
                chopper_cmd = (
                    f"cat {fq} | "
                    f"chopper --headcrop {head_trim} --tailcrop {tail_trim} -t {self.threads_per_task} | "
                    f"chopper -q {self.qualityscore} -l {self.minlength} -t {self.threads_per_task} "
                    f"> {trimmed_fastq}"
                )
                self.log("RUN", f"{chopper_cmd}")
                with open(chopper_log, "w") as log:
                    subprocess.run(
                        chopper_cmd,
                        shell=True,
                        check=True,
                        stdout=log,
                        stderr=log
                    )
                self.log("INFO", f"Trimming completed for {sample}")
            else:
                self.log("INFO", f"Trimmed file already exists for {sample}")

        # ---------- Run in parallel ----------
        self.parallel_execution(chopper_single, fastqs)

    def minimap2_parallel(self, fastqs):
        """Run minimap2 with selected parameters (see parameter_MetNano.py file)"""
        os.makedirs(self.minimap2_dir, exist_ok=True)

        def minimap2_single(fq):
            sample = os.path.basename(fq).replace(".fastq.gz", "").replace(".fastq", "")
            trimmed_fastq = os.path.join(self.chopper_dir, sample + "_trimmed.fastq")
            sam_alignment = os.path.join(self.minimap2_dir, sample + "_alignment.sam")
            minimap2_log = os.path.join(self.minimap2_dir, f"{sample}.log")
            if not os.path.isfile(sam_alignment) or os.path.getsize(sam_alignment) == 0:
                self.log("INFO", f"Mapping reads from {sample}...")
                minimap2_cmd = (
                    f"minimap2 -a -x {mapping_preset} -t {self.threads_per_task} "
                    f"{self.ref} {trimmed_fastq} "
                    f"> {sam_alignment}"
                )
                self.log("RUN", f"{minimap2_cmd}")
                with open(minimap2_log, "w") as log:
                    subprocess.run(
                        minimap2_cmd,
                        shell=True,
                        check=True,
                        stdout=log,
                        stderr=log
                    )
                self.log("INFO", f"Mapping completed for {sample}")
            else:
                self.log("INFO", f"Alignment file already exists for {sample}")

        # ---------- Run in parallel ----------
        self.parallel_execution(minimap2_single, fastqs)

    def samtools_parallel(self, fastqs):
        """Extract unmapped reads with samtools"""
        def samtools_single(fq):
            sample = os.path.basename(fq).replace(".fastq.gz", "").replace(".fastq", "")
            sam_alignment = os.path.join(self.minimap2_dir, sample + "_alignment.sam")
            unmapped_bam = sam_alignment.replace("_alignment.sam", "_unmapped.bam")
            unmapped_fq = unmapped_bam.replace("_unmapped.bam", "_unmapped.fastq")
            samtools_log = os.path.join(self.minimap2_dir, sample + ".log")
            if not os.path.isfile(unmapped_fq) or os.path.getsize(unmapped_fq) == 0:
                self.log("INFO", f"Extracting reads from {os.path.basename(unmapped_fq)}...")
                samtools_cmd = (
                    f"samtools view -bf 4 "
                    f"-T {self.ref} "
                    f"-@ {self.threads_per_task} "
                    f"{sam_alignment} "
                    f"> {unmapped_bam} && "
                    f"samtools fastq "
                    f"{unmapped_bam} "
                    f"> {unmapped_fq}"
                )
                self.log("RUN", f"{samtools_cmd}")
                with open(samtools_log, "w") as log:
                    subprocess.run(
                        samtools_cmd,
                        shell=True,
                        check=True,
                        stdout=log,
                        stderr=log
                    )
                self.log("INFO", f"Unmapped reads extracted for {sample}")
            else:
                self.log("INFO", f"Unmapped reads already exist for {sample}")

        # ---------- Run in parallel ----------
        self.parallel_execution(samtools_single, fastqs)

    def SeqKit_parallel(self, fastqs):
        """Segment reads at minimum read length with seqkit"""
        os.makedirs(self.seqkit_dir, exist_ok=True)

        def SeqKit_single(fq):
            sample = os.path.basename(fq).replace(".fastq.gz", "").replace(".fastq", "")
            if self.ref in [None, "None"]:
                input_fq = os.path.join(self.chopper_dir,sample + "_trimmed.fastq")
            else:
                input_fq = os.path.join(self.minimap2_dir, sample + "_unmapped.fastq")
            segmented_fastq = os.path.join(self.seqkit_dir, sample + "_segmented.fastq")
            if not os.path.isfile(segmented_fastq) or os.path.getsize(segmented_fastq) == 0:
                self.log("INFO", f"Segmenting reads for {sample}...")
                seqkit_cmd = (
                    f"seqkit sliding "
                    f"-W {self.minlength} "
                    f"-s {self.minlength} "
                    f"{input_fq} "
                    f"> {segmented_fastq}"
                )
                self.log("RUN", f"{seqkit_cmd}")
                subprocess.run(seqkit_cmd, shell=True, check=True)
                self.log("INFO", f"Segmenting completed for {sample}")
            else:
                self.log("INFO", f"Segmented file already exists for {sample}")

        # ---------- Run in parallel ----------
        self.parallel_execution(SeqKit_single, fastqs)

    def run_all(self):
        # Unzip and concatenate FASTQ in raw ONT folder
        if self.unzip == "True":
            print("\n# Unzip FASTQ files", flush = True)
            self.unzip_fastq()
        else:
            print("\n# Creating link to FASTQs", flush=True)
            self.create_ln_file()

        # Rename FASTQ files
        if self.rename == "True":
            print("\n# Renaming barcodes", flush=True)
            self.rename_fastq()

        # Change .fq extension to .fastq
        for fq in glob.glob(os.path.join(self.output_folder_fastq, "*.fq")):
            base = os.path.splitext(fq)[0]
            os.rename(os.path.join(self.output_folder_fastq, fq), os.path.join(self.output_folder_fastq, base + ".fastq"))

        # Get FASTQ file list for task parallelization
        print("\n# Adding parameters", flush=True)
        fastqs = glob.glob(os.path.join(self.output_folder_fastq, "*.fastq"))
        self.log("INFO", f"FASTQ file found: {len(fastqs)}")
        if not fastqs:
            self.log("ERROR", f"No FASTQ files found in {self.output_folder_fastq}")
            self.log("END", "Exiting MetNano")
            sys.exit(1)

        # Calculate optimal number of threads per task
        self.max_workers = min(len(fastqs), self.threads)
        self.threads_per_task = max(1, self.threads // self.max_workers)
        self.log("INFO", f"Threads per parallel task = {self.threads_per_task}\n")

        # Run assembly-stats
        print("# Compute sequencing statistics with assembly-stats", flush=True)
        self.assembly_stats()

        # Run NanoPlot
        print("\n# Quality control with NanoPlot", flush = True)
        self.NanoPlot_parallel(fastqs)

        # Run nanoQC
        print("\n# Quality control with nanoQC", flush = True)
        self.nanoQC_parallel(fastqs)

        # Run chopper
        print("\n# Read filtering and trimming with chopper", flush=True)
        self.chopper_parallel(fastqs)

        # If a reference genome is specified
        if hasattr(self, 'ref') and self.ref not in [None, "None"]:
            # Run minimap2
            print("\n# Mapping reads to reference genome with minimap2", flush = True)
            self.minimap2_parallel(fastqs)

            # Run samtools
            print("\n# Extraction of unmapped reads with samtools", flush = True)
            self.samtools_parallel(fastqs)

        # Run SeqKit
        print("\n# Read segmentation with SeqKit", flush=True)
        self.SeqKit_parallel(fastqs)
        print()

class TaxonomicClassification:
    """
    Python class allowing the assignment of taxonomic classification to reads with
    Kraken2 using a temporary DB copy in RAM (/dev/shm) for speed, the estimation
    of taxa abundance with Bracken and the conversion of Kraken2-style reports into
    BIOM tables with kraken-biom.
    """
    def __init__(self, db_path, fastq_dir, parameters, threads, confidence, minlength, memory_mapping):
        self.db_path = db_path
        self.db_name = os.path.basename(os.path.normpath(self.db_path))
        self.fastq_dir = fastq_dir
        self.parameters = parameters
        self.out_dir = os.path.abspath(".")
        self.threads = threads
        self.confidence = confidence
        self.minlength = minlength
        self.memory_mapping = memory_mapping
        self.temp_db_path = os.path.join(memory_mapping, self.db_name)
        self.kraken2_dir = os.path.join(self.out_dir, "06-Kraken2", self.parameters, self.db_name, self.confidence)
        self.required_files = ["hash.k2d", "opts.k2d", "taxo.k2d"]
        self.bracken_dir = os.path.join(self.out_dir, "07-Bracken", self.parameters, self.db_name, self.confidence)
        self.krakenbiom_dir = os.path.join(self.out_dir, "08-kraken-biom", self.parameters, self.db_name, self.confidence)

    def log(self, kind, message):
        """Print formatted log message"""
        print(f"{datetime.now().isoformat(sep=' ', timespec='seconds')} [{kind}] {message}", flush=True)

    def copy_db_to_ramdisk(self):
        """Load Kraken2 database in ram disk"""
        print("# Load Kraken2 database into RAM disk", flush = True)
        os.makedirs(self.temp_db_path, exist_ok=True)
        self.log("INFO", "Copying database to /dev/shm/...")
        cmd= (
            f"rsync "
            f"{self.db_path}/*  "
            f"{self.temp_db_path}"
        )
        subprocess.run(
            cmd, shell=True, check=True
        )
        for f in self.required_files:
            if not os.path.exists(os.path.join(self.temp_db_path, f)):
                self.log("ERROR", f"Database copy incomplete: missing {f}")
                self.log("END", "Exiting MetNano")
                sys.exit(1)
        self.log("INFO", "Database copy completed\n")

    def kraken2(self, fastq_file):
        """Assign taxonomic classification to reads with Kraken2"""
        os.makedirs(self.kraken2_dir, exist_ok=True)
        sample_name = os.path.basename(fastq_file).replace("_segmented.fastq", "")
        kraken2_output = os.path.join(self.kraken2_dir, f"{sample_name}.all")
        kraken2_report = os.path.join(self.kraken2_dir, f"{sample_name}.report")
        kraken2_log = os.path.join(self.kraken2_dir, f"{sample_name}.log")
        if self.memory_mapping != "False":
            missing_files = [f for f in self.required_files if not os.path.exists(os.path.join(self.temp_db_path, f))]
            if missing_files:
                self.log("ERROR", f"Missing required DB files in RAM: {missing_files}")
                self.log("END", "Exiting MetNano")
                sys.exit(1)
            kraken2_cmd = (
                f"kraken2 "
                f"--db {self.temp_db_path} "
                f"--threads {self.threads} "
                f"--output {kraken2_output} "
                f"--report {kraken2_report} "
                f"--confidence {self.confidence} "
                f"--memory-mapping "
                f"{fastq_file}"
            )
        else:
            missing_files = [f for f in self.required_files if not os.path.exists(os.path.join(self.db_path, f))]
            if missing_files:
                self.log("ERROR", f"Missing required DB files in RAM: {missing_files}")
                self.log("END", "Exiting MetNano")
                sys.exit(1)
            kraken2_cmd = (
                f"kraken2 "
                f"--db {self.db_path} "
                f"--threads {self.threads} "
                f"--output {kraken2_output} "
                f"--report {kraken2_report} "
                f"--confidence {self.confidence} "
                f"{fastq_file}"
            )
        self.log("INFO", f"Running Kraken2 for {sample_name}")
        with open(kraken2_log, "a") as log:
            subprocess.run(
                kraken2_cmd,
                shell=True,
                check=True,
                stdout=log,
                stderr=log
            )

    def bracken(self, fastq_file):
        """Estimate relative abundance of all taxonomic levels with Bracken"""
        os.makedirs(self.bracken_dir, exist_ok=True)
        sample_name = os.path.basename(fastq_file).replace("_segmented.fastq", "")
        kraken2_report = os.path.join(self.kraken2_dir, f"{sample_name}.report")
        bracken_log = os.path.join(self.bracken_dir, f"{sample_name}.log")
        self.BRACKEN_LEVELS = {
            "Species": "S",
            "Genus": "G",
            "Family": "F",
            "Order": "O",
            "Class": "C",
            "Phylum": "P"
        }
        exceptions = []
        for level_name, level_code in self.BRACKEN_LEVELS.items():
            bracken_dir_taxrank = os.path.join(self.bracken_dir, level_name)
            os.makedirs(bracken_dir_taxrank, exist_ok = True)
            bracken_output = os.path.join(
                bracken_dir_taxrank,
                f"{sample_name}.bracken"
            )
            bracken_report = os.path.join(
                bracken_dir_taxrank,
                f"{sample_name}.breport"
            )
            if not os.path.exists(bracken_output) or os.path.getsize(bracken_output) == 0:
                self.log("INFO", f"Running Bracken at {level_name} level for sample {sample_name}...")
                bracken_cmd = (
                    f"bracken "
                    f"-d {self.db_path} "
                    f"-i {kraken2_report} "
                    f"-o {bracken_output} "
                    f"-w {bracken_report} "
                    f"-r {self.minlength} "
                    f"-l {level_code}"
                )
                self.log("RUN", bracken_cmd)
                try:
                    with open(bracken_log, "a") as log:
                        subprocess.run(
                            bracken_cmd,
                            shell=True,
                            check=True,
                            stdout=log,
                            stderr=log
                        )
                except subprocess.CalledProcessError:
                    self.log("ERROR", f"Bracken failed at {level_name} level for sample {sample_name}")
                    exceptions.append(level_name)
            else:
                self.log("INFO", f"Bracken {level_name} file already exists for sample {sample_name}")
        if not exceptions:
            self.log("INFO", f"Bracken completed for sample {sample_name} for all taxonomic levels\n")
        else:
            self.log("INFO", f"Bracken completed for sample {sample_name} for all taxonomic levels, except for {', '.join(exceptions)} levels\n")

    def krakenbiom(self, fastq_file):
        """Generate BIOM tables, one per report"""
        os.makedirs(self.krakenbiom_dir, exist_ok=True)
        sample_name = os.path.basename(fastq_file).replace("_segmented.fastq", "")
        self.log("INFO", f"Generating BIOM from Kraken2-style reports for {sample_name}...")
        exceptions = []
        for level_name, level_code in self.BRACKEN_LEVELS.items():
            bracken_dir_taxrank = os.path.join(self.bracken_dir, level_name)
            krakenbiom_dir_taxrank = os.path.join(self.krakenbiom_dir, level_name)
            os.makedirs(krakenbiom_dir_taxrank, exist_ok = True)
            bracken_report = os.path.join(
                bracken_dir_taxrank,
                f"{sample_name}.breport"
            )
            krakenbiom_table = os.path.join(
                krakenbiom_dir_taxrank,
                f"{sample_name}.biom"
            )
            if not os.path.exists(krakenbiom_table) or os.path.getsize(krakenbiom_table) == 0:
                self.log("INFO", f"Running Bracken at {level_name} level for sample {sample_name}...")
                krakenbiom_cmd = (
                    f"kraken-biom {bracken_report} "
                    f"--fmt json "
                    f"--output_fp {krakenbiom_table} "
                    f"--max D --min S"
                )
                self.log("RUN", krakenbiom_cmd)
                try:
                    subprocess.run(krakenbiom_cmd, shell=True, check=True)
                except subprocess.CalledProcessError:
                    self.log("ERROR", f"kraken-biom failed at {level_name} level for sample {sample_name}")
                    exceptions.append(level_name)
            else:
                self.log("INFO", f"BIOM for {sample_name} already exists")
        if not exceptions:
            self.log("INFO", f"kraken-biom completed for sample {sample_name} for all taxonomic levels\n")
        else:
            self.log("INFO", f"kraken-biom completed for sample {sample_name} for all taxonomic levels, except for {', '.join(exceptions)} levels\n")

    def run_all(self):
        fastq_files = glob.glob(os.path.join(self.fastq_dir, "*.fastq*"))
        if not fastq_files:
            self.log("ERROR", "No FASTQ files found.")
            self.log("END", "Exiting MetNano")
            sys.exit(1)
        files_to_run = []
        for f in fastq_files:
            sample_name = os.path.basename(f).replace("_segmented.fastq", "")
            expected_report = os.path.join(self.kraken2_dir, f"{sample_name}.report")
            if not os.path.exists(expected_report) or os.path.getsize(expected_report) == 0:
                files_to_run.append(f)

        # Only load RAM if needed
        if files_to_run:
            if self.memory_mapping != "False":
                self.copy_db_to_ramdisk()

            # Run Kraken2
            print("# Read taxonomic classification with Kraken2", flush = True)
            for f in files_to_run:
                try:
                    self.kraken2(f)
                except subprocess.CalledProcessError as e:
                    self.log("ERROR", f"ERROR processing {f}")
                    self.log("ERROR", f"{e}")
                    self.log("END", "Exiting MetNano")
                    sys.exit(1)
            self.log("END", "Kraken2 ended\n")
        else:
            print("# Read taxonomic classification with Kraken2", flush = True)
            self.log("INFO", "Kraken2 output files already exist")
            self.log("INFO", "Not loading DB into memory")
            self.log("INFO", "Skipping Kraken2\n")
        # Optional DB cleanup here

        # Run Bracken
        print("# Taxon abundance estimation with Bracken", flush = True)
        for f in fastq_files:
            try:
                self.bracken(f)
            except subprocess.CalledProcessError as e:
                self.log("ERROR", f"ERROR processing {f}")
                self.log("ERROR", f"{e}")
                self.log("END", "Exiting MetNano")
                sys.exit(1)

        # Run kraken-biom
        print("# Converting Kraken2-style reports to BIOM format with kraken-biom", flush = True)
        for f in fastq_files:
            try:
                self.krakenbiom(f)
            except subprocess.CalledProcessError as e:
                self.log("ERROR", f"ERROR processing {f}")
                self.log("ERROR", f"{e}")
                self.log("END", "Exiting MetNano")
                sys.exit(1)

        self.log("END", "MetNano is now complete. Bye!")
