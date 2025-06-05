# apscale
Advanced Pipeline for Simple yet Comprehensive AnaLysEs of DNA metabarcoding data

[![Downloads](https://static.pepy.tech/badge/apscale)](https://pepy.tech/project/apscale)  - apscale

[![Downloads](https://static.pepy.tech/badge/apscale-nanopore)](https://pepy.tech/project/apscale-nanopore)  - apscale-nanopore

# apscale-nanopore

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
Apscale-nanopore requires Python 3.10 or higher and can be easily installed via pip in any command line:

`pip install apscale_nanopore`

To update apscale-blast run:

`pip install --upgrade apscale_nanopore`

The easiest installation option is the [Conda apscale environment](https://github.com/TillMacher/apscale_installer). This way, all dependencies will automatically be installed.

Then activate the conda environment.

`conda activate apscale`

## Create project

First, create a new project:

`apscale_nanopore create -p PATH/TO/PROJECT`

A new project will be created. Follow the instructions and fill out the settings file accordingly.

<pre>
/YOUR_PROJECT_PATH/My_new_project/
├───1_raw_data
│   └───data
├───2_index_demultiplexing
│   └───data
├───3_primer_trimming
│   └───data
├───4_quality_filtering
│   └───data
├───5_clustering_denoising
│   └───data
├───6_read_table
│   └───data
├───7_taxonomic_assignment
│   └───data
├───8_nanopore_report
My_new_project_settings.xlsx
</pre>

## Settings file

### Sample index and primer combinations (Example)

| Forward index 5'-3'         | Forward primer 5'-3'                  | Reverse index 5'-3'         | Reverse primer 5'-3'                 | ID        |
|-----------------------------|----------------------------------------|-----------------------------|---------------------------------------|-----------|
| AGAACGACTTCCATACTCGTGTGA    | RGCHTTYCCHCGWATAAAYAAYATAAG            | AGAACGACTTCCATACTCGTGTGA    | GRGGRTAWACWGTTCAWCCWGTNCC              | Sample_1  |
| AACGAGTCTCTTGGGACCCATAGA    | RGCHTTYCCHCGWATAAAYAAYATAAG            | AACGAGTCTCTTGGGACCCATAGA    | GRGGRTAWACWGTTCAWCCWGTNCC              | Sample_2  |
| AGGTCTACCTCGCTAACACCACTG    | RGCHTTYCCHCGWATAAAYAAYATAAG            | AGGTCTACCTCGCTAACACCACTG    | GRGGRTAWACWGTTCAWCCWGTNCC              | Sample_3  |
| CGTCAACTGACAGTGGTTCGTACT    | RGCHTTYCCHCGWATAAAYAAYATAAG            | CGTCAACTGACAGTGGTTCGTACT    | GRGGRTAWACWGTTCAWCCWGTNCC              | Sample_4  |
| ACCCTCCAGGAAAGTACCTCTGAT    | RGCHTTYCCHCGWATAAAYAAYATAAG            | ACCCTCCAGGAAAGTACCTCTGAT    | GRGGRTAWACWGTTCAWCCWGTNCC              | Sample_5  |
| CCAAACCCAACAACCTAGATAGGC    | RGCHTTYCCHCGWATAAAYAAYATAAG            | CCAAACCCAACAACCTAGATAGGC    | GRGGRTAWACWGTTCAWCCWGTNCC              | Sample_6  |

### Apscale-nanopore settings  (Example)

| Step                     | Category            | Variable                   | Comment                                                                 |
|--------------------------|---------------------|----------------------------|-------------------------------------------------------------------------|
| General                  | cpu count           | 7                          | Number of cores to use                                                 |
| demultiplexing (index)   | allowed errors index| 3                          | Allowed errors during index demultiplexing                             |
| primer trimming          | allowed errors primer| 4                         | Allowed errors during primer trimming                                  |
| quality filtering        | minimum length      | 54                         | Reads below this length will be discarded                              |
| quality filtering        | maximum length      | 74                         | Reads above this length will be discarded                              |
| quality filtering        | minimum quality     | 20                         | Reads below this average PHRED quality score will be discarded         |
| clustering/denoising     | mode                | denoised OTUs              | Choose clustering/denoising algorithm                                  |
| clustering/denoising     | percid              | 0.97                       | Vsearch clustering percentage identity                                 |
| clustering/denoising     | alpha               | 1                          | Vsearch denoising alpha value                                          |
| clustering/denoising     | d                   | 1                          | Swarm's d value                                                        |
| read table               | minimum reads       | 10                         | Discard reads below this threshold                                     |
| taxonomic assignment     | apscale blast       | Yes                        | Run APSCALE megablast (yes or no)                                      |
| taxonomic assignment     | apscale db          | ...                        | Path to local database                                                 |

## Run apscale-nanopore

Apscale-nanopore operates in four different ways:

### 1) Raw-data processing of non-demultiplexed data

* Copy your non-demultiplexed .fastq(.gz) files to the "1_raw_data/data" folder.
   
`apscale_nanopore run -p PATH/TO/PROJECT`

* Apscale-nanopore will demultiplex all your files according to the demultiplexing sheet.

### 2) Raw-data processing of demultiplexed data
   
* Copy your demultiplexed .fastq(.gz) files to the "1_raw_data/data" folder.
   
`apscale_nanopore run -p PATH/TO/PROJECT -sd`

* Apscale-nanopore will skip the demultiplexing and immediately start with the raw-data processing.
* Important: Enter the primer sequences (5'-3') in the first row of the demultiplexing sheet. The index columns can be left blank.

### 3) Live raw-data processing of non-demultiplexed data
   
* Output your non-demultiplexed .fastq(.gz) files to the "1_raw_data/data" folder during sequencing.
* Apscale-nanopore will automatically scan the folder for incoming files and automatically process them.
* Press Ctrl+C to interupt the live-calling.
   
`apscale_nanopore run -p PATH/TO/PROJECT -l`

* Apscale-nanopore will demultiplex all your files according to the demultiplexing sheet.

### 4) Live raw-data processing of demultiplexed data
   
* Output your demultiplexed .fastq(.gz) files to the "1_raw_data/data" folder during sequencing.
* Apscale-nanopore will automatically scan the folder for incoming files and automatically process them.
* Press Ctrl+C to interupt the live-calling.
   
`apscale_nanopore run -p PATH/TO/PROJECT -l -sd`

* Apscale-nanopore will skip the demultiplexing and immediately start with the raw-data processing.
* Important: Enter the primer sequences (5'-3') in the first row of the demultiplexing sheet. The index columns can be left blank.

### Run individual steps

* Apscale can run individual steps (-step X) or all steps after a specific module (-steps X).

Step indices:
* 1 = Index demultiplexing
* 2 = Primer trimming
* 3 = Quality filtering
* 4 = Clustering/denoising
* 5 = Read table
* 6 = Taxonomic assignment

Example: Run "clustering/denoising"

`apscale_nanopore run -p PATH/TO/PROJECT -step 4`

Example: Run all steps after the "quality filtering":

`apscale_nanopore run -p PATH/TO/PROJECT -steps 3`

## Quality control

A quality control can be conducted for all fastq files. Simply run:

`apscale_nanopore qc -p PATH/TO/PROJECT`

## Bioinformatics Workflow Overview

### 1) Demultiplexing
**Tool:** `cutadapt`
**Settings:** `Allowed errors (default=3)`

Demultiplex raw sequencing reads based on barcode sequences to generate sample-specific FASTQ files.

---

### 2) Primer Trimming  
**Tool:** `cutadapt`
**Settings:** `Allowed errors (default=4)`

Remove primer sequences from demultiplexed reads to retain only target regions.

---

### 3) Quality Filtering  
**Tools:** `python`, `vsearch`
**Settings:** `Min. mean Q-Score (default=20), Min. and max. length (fragment-specific)`

Filter reads based on:
- Mean PHRED quality score  
- Minimum and maximum fragment length  

This step ensures only high-quality reads are retained for downstream processing.

---

### 4) Clustering / Denoising  

**Tool:** `vsearch`  
**Settings:** `d (default=1), percentage identity (default=0.97), alpha (default=1)`

Choose from the following processing strategies:

- **Swarm denoising**: Local clustering using the Swarm algorithm for fine-scale resolution.
- **Swarm OTUs**: Swarm denoising followed by similarity clustering.  
- **ESV denoising**: Error-correction to obtain Exact Sequence Variants.
- **Denoised OTUs**: Denoising followed by similarity clustering.  


---

### 5) Read Table Construction and Filtering  

**Tool:** `python`  
**Settings:** `minimum reads (default=10)`

Construct an abundance table (ESVs/OTUs × samples).  

Apply a minimum read threshold to remove low-abundance features.

---

### 6) Taxonomic Assignment

**Tool:** `BLASTn` via [`apscale-blast`](https://github.com/TillMacher/apscale_blast)

**Settings:** `Apscale-blast database`

Assign taxonomy to representative sequences using a local reference database.

---

### 7) Quality Control and Reporting

**Tool:** `python`

Generate summary statistics and visual diagnostics.
