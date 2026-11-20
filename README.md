<<<<<<< HEAD
# APSCALE-Nanopore
=======
# APSCALE
Advanced Pipeline for Simple yet Comprehensive AnaLysEs of DNA metabarcoding data
>>>>>>> 97f05d31f8015ae398ed876484a6c08d18cc529f

[![Downloads](https://static.pepy.tech/badge/apscale-nanopore)](https://pepy.tech/project/apscale-nanopore)

Metabarcoding pipeline for **Oxford Nanopore** data, based on [APSCALE](https://github.com/DominikBuchner/apscale). APSCALE-Nanopore takes you from the raw Nanopore reads to the read table and the taxonomic assignment: index demultiplexing, primer trimming, quality filtering, clustering/denoising, read table and taxonomic assignment with [APSCALE-blast](https://github.com/TillMacher/apscale_blast). It can be used with a graphical user interface (runs in your browser) or from the command line, and also supports live processing during sequencing.

<<<<<<< HEAD
* **Input:** Nanopore sequencing data (`.fastq` or `.fastq.gz`), demultiplexed or not
* **Output:** read table (`.xlsx`/`.parquet`), OTU/ESV sequences (`.fasta`), taxonomy table (`.xlsx`) and quality/summary reports (`.html`) in your project folder, ready for [TaxonTableTools2](https://github.com/TillMacher/TaxonTableTools2)

## Installation

The easiest way is the conda environment of the [APSCALE installer](https://github.com/TillMacher/apscale_installer), which installs APSCALE-Nanopore and all its dependencies (APSCALE-blast, cutadapt, VSEARCH, SWARM, BLAST+).

Alternatively, install via pip (Python ≥ 3.10):

```
pip install apscale_nanopore
```

With pip, [VSEARCH](https://github.com/torognes/vsearch), [SWARM](https://github.com/torognes/swarm) and [BLAST+](https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/) must be installed separately and be available in your PATH. APSCALE-Nanopore checks whether they are installed and tells you what is missing.
=======
# APSCALE-Nanopore

## Introduction
Apscale-nanopore is a modified version of the metabarcoding pipeline [apscale](https://github.com/DominikBuchner/apscale/tree/main) and is used
for the processing of Oxford Nanopore data.

Programs used:
* [cutadapt](https://github.com/marcelm/cutadapt) 
* [vsearch](https://github.com/torognes/vsearch)
* [swarm](https://github.com/torognes/swarm)
* [blast+](https://blast.ncbi.nlm.nih.gov/doc/blast-help/downloadblastdata.html) (blastn module)

Input:
* Non-demultiplexed Nanopore sequence data in .fastq format.
* Demultiplexed Nanopore sequence data in .fastq format.

Output:
* read table, taxonomy table, log files, report

## Installation

Apscale-nanopore can be installed on all common operating systems (Windows, Linux, MacOS).
Apscale-nanopore requires Python 3.10 or higher.
>>>>>>> 97f05d31f8015ae398ed876484a6c08d18cc529f

To download all databases at once from the GUI, the browser used for the download must be installed once:

```
playwright install chromium
```

<<<<<<< HEAD
To update:
=======
`conda activate apscale4`

Alternatively, you can install Apscale-nanopore via pip:

`pip install apscale_nanopore`

To update apscale-nanopore run:

`pip install --upgrade apscale_nanopore`

## Graphical-User-Infercace

The GUI can be started with the following command:

`apscale_nanopore gui`
>>>>>>> 97f05d31f8015ae398ed876484a6c08d18cc529f

```
pip install --upgrade apscale_nanopore
```

## Graphical user interface

```
conda activate apscale4
apscale_nanopore gui
```

The graphical interface opens in your browser. Enter the path to your projects folder (e.g. an empty folder named `APSCALE_projects`) in the sidebar, create the database folder, and create or select a project. Click **"🔄 Refresh files and folders"** in the sidebar when you have added new files.

### Pages

**Start**
* **Home:** overview of the active project (raw files, demultiplexed files, FASTA files, databases), installed tools and package versions, links and citations

**Raw Data Processing**
* **APSCALE-Nanopore:** download databases, adjust the settings (index demultiplexing, primer trimming, quality filtering, clustering/denoising, read table, taxonomic assignment) and run the selected steps. The settings are saved to the project's settings file when the analysis starts. For demultiplexed data, the raw files can be imported into the Demultiplexing sheet with one click

**Reports**
* **Quality Report:** create and view quality plots (mean PHRED score, read length, quality vs. length) for the raw, demultiplexed or primer-trimmed data
* **Summary Report:** reads and OTUs/ESVs per sample, created automatically after each run

Live processing during sequencing is only available from the command line (see below).

## Command line

### Create a project

```
apscale_nanopore create -p PATH/TO/My_project
```

A new project folder (`My_project_apscale_nanopore`) is created. You are asked whether your data is already demultiplexed, and the matching settings file is created. Existing projects are never overwritten.

### Run APSCALE-Nanopore

APSCALE-Nanopore can be run in four different ways:

**1) Non-demultiplexed data**

Copy your non-demultiplexed `.fastq(.gz)` files to `1_raw_data/data`. All files are demultiplexed according to the Demultiplexing sheet.

```
apscale_nanopore run -p PATH/TO/My_project_apscale_nanopore
```

**2) Demultiplexed data**

Copy your demultiplexed `.fastq(.gz)` files to `1_raw_data/data` (one file per sample) and add them to the Demultiplexing sheet:

```
apscale_nanopore import_raw_data -p PATH/TO/My_project_apscale_nanopore
```

This writes one row per file (the file name is the sample ID) and copies the primers and lengths of the first row to all samples. Check the sheet, then run with `-sd` (skip demultiplexing):

```
apscale_nanopore run -p PATH/TO/My_project_apscale_nanopore -sd
```

**3) Live processing of non-demultiplexed data**

Let the sequencer write your `.fastq(.gz)` files to `1_raw_data/data`. APSCALE-Nanopore scans the folder for new files and processes them automatically. Press `Ctrl+C` to stop.

```
apscale_nanopore run -p PATH/TO/My_project_apscale_nanopore -live
```

**4) Live processing of demultiplexed data**

```
apscale_nanopore run -p PATH/TO/My_project_apscale_nanopore -live -sd
```

### Run individual steps

Run a single step (`-step X`) or all steps from a specific step onwards (`-steps X`):

* 1 = Index demultiplexing
* 2 = Primer trimming
* 3 = Quality filtering
* 4 = Clustering/denoising
* 5 = Read table
* 6 = Taxonomic assignment

Example: run only the clustering/denoising

```
apscale_nanopore run -p PATH/TO/My_project_apscale_nanopore -step 4
```

Example: run all steps from the quality filtering onwards

```
apscale_nanopore run -p PATH/TO/My_project_apscale_nanopore -steps 3
```

Re-running the index demultiplexing replaces the previous results (reads are not added twice).

### Change settings on the fly

Settings can be changed for a single run without editing the settings file:

| Option | Setting |
|---|---|
| `-e1` | allowed errors index |
| `-e2` | allowed errors primer |
| `-minlen` / `-maxlen` | minimum / maximum length (for all samples) |
| `-minq` | minimum quality |
| `-mode` | clustering/denoising mode |
| `-percid` | vsearch clustering percentage identity |
| `-alpha` | vsearch denoising alpha |
| `-d` | swarm's d value |
| `-rep` | representative sequence (centroid or consensus) |
| `-minreads` | read filter threshold |

Example:

```
apscale_nanopore run -p PATH/TO/My_project_apscale_nanopore -steps 4 -mode "Swarm OTUs" -percid 0.98
```

### Quality control

Create quality plots for the raw, demultiplexed or primer-trimmed data (you are asked which folder):

```
apscale_nanopore qc -p PATH/TO/My_project_apscale_nanopore
```

## Project folder

```
APSCALE_projects/
├── APSCALE_databases/
└── My_project_apscale_nanopore/
    ├── 1_raw_data/data                <- raw data (non-demultiplexed or demultiplexed)
    ├── 2_index_demultiplexing/data
    ├── 3_primer_trimming/data
    ├── 4_quality_filtering/data
    ├── 5_clustering_denoising/data
    ├── 6_read_table                   <- read table (.xlsx/.parquet)
    │   └── data                       <- OTU/ESV sequences (.fasta)
    ├── 7_taxonomic_assignment         <- taxonomy table
    ├── 8_nanopore_report              <- quality reports (data) and summary report (main)
    └── My_project_settings.xlsx
```

## Settings file

The settings file has two sheets.

### Demultiplexing sheet

One row per sample. The index columns are only required for non-demultiplexed data and can be left empty otherwise. The minimum and maximum length are set per sample, so different markers can be processed in the same run.

| Forward index 5-3 | Forward primer 5-3 | Reverse index 5-3 | Reverse primer 5-3 | ID | Minimum length | Maximum length |
|---|---|---|---|---|---|---|
| AGAACGACTTCCATACTCGTGTGA | AAACTCGTGCCAGCCACC | AGAACGACTTCCATACTCGTGTGA | GGGTATCTAATCCCAGTTTG | Sample_1 | 157 | 177 |
| AACGAGTCTCTTGGGACCCATAGA | AAACTCGTGCCAGCCACC | AACGAGTCTCTTGGGACCCATAGA | GGGTATCTAATCCCAGTTTG | Sample_2 | 157 | 177 |

### Settings sheet (defaults)

| Category | Default | Comment |
|---|---|---|
| cpu count | auto | Number of cores to use (auto = all cores minus two) |
| allowed errors index | 3 | Allowed errors during index demultiplexing |
| allowed errors primer | 3 | Allowed errors during primer trimming |
| minimum quality | 20 | Reads below this average PHRED quality score are discarded |
| mode | denoised OTUs | ESVs, Swarms, denoised OTUs or Swarm OTUs |
| percid | 0.97 | vsearch clustering percentage identity |
| alpha | 1 | vsearch denoising alpha value |
| d | 1 | swarm's d value |
| representative | centroid | Representative sequence of OTUs: centroid or consensus |
| minimum reads | 0.0001 | Read filter per sample: < 1 = relative (0.0001 = 0.01% of the sample's reads), ≥ 1 = absolute number of reads |
| apscale blast | Yes | Run the taxonomic assignment (Yes or No) |
| task | megablast | BLAST task (megablast, blastn, dc-megablast) |
| apscale db | | Path to the APSCALE database (selected automatically in the GUI) |

## Bioinformatic workflow

**1) Index demultiplexing** (`cutadapt`) – sorts the reads into samples by their index combination. Reads are searched in both orientations (5'-3' and 3'-5').

**2) Primer trimming** (`cutadapt`) – removes the primers (both orientations); reads without primers are discarded.

**3) Quality filtering** (`python`, `vsearch`) – removes reads below the minimum mean PHRED score and outside the length range of the sample, then dereplicates the reads.

**4) Clustering/denoising** (`vsearch`, `swarm`) – per sample, followed by de novo chimera removal (`vsearch --uchime3_denovo`):
* **ESVs:** denoising (UNOISE3) to exact sequence variants
* **Swarms:** SWARM clustering; the seed (= centroid) is always used as representative
* **denoised OTUs:** denoising followed by similarity clustering
* **Swarm OTUs:** SWARM clustering followed by similarity clustering

The read count of each OTU/ESV includes the reads of all its cluster members.

**5) Read table** (`python`) – builds the read table (OTUs/ESVs × samples) and removes low-abundance OTUs/ESVs per sample (minimum reads).

**6) Taxonomic assignment** (`BLASTn` via [APSCALE-blast](https://github.com/TillMacher/apscale_blast)) – assigns taxonomy using a local reference database.

**7) Reports** (`python`, `plotly`) – summary report after each run and optional quality reports.

## Databases

Pre-compiled databases for APSCALE-blast are [available here](https://seafile.rlp.net/d/c172d076de1e4c45b594/) and are updated regularly (the same databases as in [APSCALE-GUI](https://github.com/TillMacher/apscale_gui)). They are stored in `APSCALE_projects/APSCALE_databases` and can be downloaded directly in the GUI. Each database folder contains a `citation.txt`. The database path must not contain spaces.

## Citation

APSCALE-Nanopore: coming soon. Until then, please cite APSCALE:

* **APSCALE:** Buchner, Macher & Leese (2022) Bioinformatics 38: 4817–4819. https://doi.org/10.1093/bioinformatics/btac588

Please also cite the programs and data you used:

* **cutadapt:** Martin (2011) EMBnet.journal 17: 10–12. https://doi.org/10.14806/ej.17.1.200
* **VSEARCH:** Rognes et al. (2016) PeerJ 4: e2584. https://doi.org/10.7717/peerj.2584
* **SWARM:** Mahé et al. (2022) Bioinformatics 38: 267–269. https://doi.org/10.1093/bioinformatics/btab493
* **BLAST+:** Camacho et al. (2009) BMC Bioinformatics 10: 421. https://doi.org/10.1186/1471-2105-10-421
* **Databases:** see the `citation.txt` in the database folder and the [APSCALE-blast citations](https://github.com/TillMacher/apscale_blast#citation)

## Reporting issues

Please open an [issue](https://github.com/TillMacher/apscale_nanopore/issues) or contact me via email (macher@uni-trier.de).
