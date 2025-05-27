# apscale
Advanced Pipeline for Simple yet Comprehensive AnaLysEs of DNA metabarcoding data

[![Downloads](https://static.pepy.tech/badge/apscale)](https://pepy.tech/project/apscale)  - apscale

[![Downloads](https://static.pepy.tech/badge/apscale-blast)](https://pepy.tech/project/apscale-blast)  - apscale_blast

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

Output:
* read table, taxonomy table, log files

## Installation

Apscale-nanopore can be installed on all common operating systems (Windows, Linux, MacOS).
Apscale-nanopore requires Python 3.10 or higher and can be easily installed via pip in any command line:

`pip install apscale_nanopore`

To update apscale-blast run:

`pip install --upgrade apscale_nanopore`

The easiest installation option is the [Conda apscale environment](https://github.com/TillMacher/apscale_installer). This way, all dependencies will automatically be installed.
