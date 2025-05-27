import argparse
import glob
import multiprocessing
import os
import sys
from datetime import datetime
from pathlib import Path
import time
import shutil
import gzip
import pandas as pd
from tqdm import tqdm
import subprocess
from joblib import Parallel, delayed
from Bio import SeqIO
import hashlib
from collections import defaultdict
import multiprocessing
import os
import platform
import subprocess
import time
from Bio.SeqIO import SeqRecord
from Bio.Seq import Seq
import itertools
import time

# All settings
index_error_values = [i for i in range(1,4)]
primer_error_values = [i for i in range(1,4)]
minq = [i for i in range(10,41,10)]
target_len = 64
plus_minus = [i for i in range(10,21,10)]
maxmin_values = [[target_len+i, target_len-i] for i in plus_minus]
d_values = [1,2,3]
min_read_values = [2,10,20]

# Generate all combinations
combinations = list(itertools.product(
    index_error_values,
    primer_error_values,
    minq,
    maxmin_values,
    d_values,
    min_read_values
))

n_combinations = len(combinations) # 1350 combinations
print(f'Number of test combinations: {n_combinations}')
res = []

combinations_file = Path('/Users/tillmacher/Desktop/APSCALE_projects/test_dataset_apscale_nanopore/10_benchmark/benchmark_results.xlsx')
if combinations_file.exists():
    combinations_df = pd.read_excel(combinations_file).fillna('')
    existing_ids = combinations_df['id'].values.tolist()
else:
    existing_ids = []

for combination in combinations:

    # record start time
    start_time = time.time()

    # Collect settings
    e_index_value = combination[0]
    e_primer_value = combination[1]
    e_tag_value = combination[2]
    t_value = combination[3]
    max_value = combination[4][0]
    min_value = combination[4][1]
    maxee_value = combination[5]
    d_value = combination[6]
    minreads_value = combination[7]

    # Save settings
    id ='-'.join([str(i) for i in [e_index_value, e_primer_value, e_tag_value, t_value, max_value, min_value, maxee_value, d_value, minreads_value]])
    run = [id]

    if id in existing_ids:
        print(f'Combination "{id}" already exists.')
    else:
        run.extend([e_index_value, e_primer_value, e_tag_value, t_value, max_value, min_value, maxee_value, d_value, minreads_value])

        # Run apscale nanopore
        command = f"python3 /Users/tillmacher/Documents/GitHub/apscale_nanopore/apscale_nanopore/__main__.py run -p /Users/tillmacher/Desktop/APSCALE_projects/test_dataset_apscale_nanopore -e1 {e_index_value} -e2 {e_primer_value} -e3 {e_tag_value} -minlen {min_value} -maxlen {max_value} -maxee {maxee_value} -d {d_value} -t {t_value} -minreads {minreads_value}"
        process = subprocess.Popen(command, shell=True, text=True)
        process.wait()

        # Compare results and create report
        # Reference
        reference_df = pd.read_excel('/Users/tillmacher/Desktop/APSCALE_projects/test_dataset_apscale_nanopore/10_benchmark/merged_main_new/merged_main_new_taxonomy.xlsx').fillna('')
        reference_ESVs = reference_df['unique ID'].values.tolist()
        reference_ESVs_n = len(reference_ESVs)
        reference_species = [i for i in reference_df['Species'].drop_duplicates().values.tolist() if i != '']
        reference_genera = [i for i in reference_df['Genus'].drop_duplicates().values.tolist() if i != '']
        reference_families = [i for i in reference_df['Family'].drop_duplicates().values.tolist() if i != '']

        # Nanopore
        nanopore_df = pd.read_excel('/Users/tillmacher/Desktop/APSCALE_projects/test_dataset_apscale_nanopore/8_taxonomic_assignment/test_dataset/test_dataset_taxonomy.xlsx').fillna('')
        nanopore_ESVs = nanopore_df['unique ID'].values.tolist()
        nanopore_ESVs_n = len(nanopore_ESVs)
        nanopore_species = [i for i in nanopore_df['Species'].drop_duplicates().values.tolist() if i != '']
        nanopore_genera = [i for i in nanopore_df['Genus'].drop_duplicates().values.tolist() if i != '']
        nanopore_families = [i for i in nanopore_df['Family'].drop_duplicates().values.tolist() if i != '']

        # Comparison
        # ESVs
        reference_only = len(set(reference_ESVs) - set(nanopore_ESVs))
        shared = len(set(reference_ESVs) & set(nanopore_ESVs))
        nanopore_only = len(set(nanopore_ESVs) - set(reference_ESVs))
        # Save settings
        run.extend([reference_only, shared, nanopore_only])

        # Species
        reference_only = len(set(reference_species) - set(nanopore_species))
        shared = len(set(reference_species) & set(nanopore_species))
        nanopore_only = len(set(nanopore_species) - set(reference_species))
        # Save settings
        run.extend([reference_only, shared, nanopore_only])

        # Genera
        reference_only = len(set(reference_genera) - set(nanopore_genera))
        shared = len(set(reference_genera) & set(nanopore_genera))
        nanopore_only = len(set(nanopore_genera) - set(reference_genera))
        # Save settings
        run.extend([reference_only, shared, nanopore_only])

        # Families
        reference_only = len(set(reference_families) - set(nanopore_families))
        shared = len(set(reference_families) & set(nanopore_families))
        nanopore_only = len(set(nanopore_families) - set(reference_families))
        # Save settings
        run.extend([reference_only, shared, nanopore_only])

        # record end time
        end_time = time.time()
        # calculate elapsed time
        elapsed_time = end_time - start_time
        run.append(elapsed_time)

        res.append(run)

        # Clean-up
        demultiplexing_folder = Path('/Users/tillmacher/Desktop/APSCALE_projects/test_dataset_apscale_nanopore/2_index_demultiplexing/data')
        if demultiplexing_folder.exists():
            shutil.rmtree(demultiplexing_folder)
            os.makedirs(demultiplexing_folder, exist_ok=True)

        demultiplexing_folder = Path('/Users/tillmacher/Desktop/APSCALE_projects/test_dataset_apscale_nanopore/4_tag_demultiplexing/data')
        if demultiplexing_folder.exists():
            shutil.rmtree(demultiplexing_folder)
            os.makedirs(demultiplexing_folder, exist_ok=True)

        main_file = Path('/Users/tillmacher/Desktop/APSCALE_projects/test_dataset_apscale_nanopore/1_raw_data/tmp/merged_nanopore_data.fastq.gz')
        moved_main_file = Path('/Users/tillmacher/Desktop/APSCALE_projects/test_dataset_apscale_nanopore/1_raw_data/data/merged_nanopore_data.fastq.gz')
        shutil.move(main_file, moved_main_file)

        # Write Excel
        df = pd.DataFrame(res)
        df.columns = ['id', 'index_error', 'primer_error', 'tag_error' ,'truncation', 'maxlen', 'minlen', 'maxee', 'd', 'minreads', 'ESVs reference', 'ESVs shared', 'ESVs nanopore', 'Species reference', 'Species shared', 'Species nanopore', 'Genus reference', 'Genus shared', 'Genus nanopore', 'Family reference', 'Family shared', 'Family nanopore', 'elapsed_time']
        df2 = pd.concat([df, combinations_df], ignore_index=True)

        res = []
        for _, row in df2.iterrows():
            res_row = []
            for test in ['ESVs', 'Species', 'Genus', 'Family']:
                # calculate sensitivity and precision
                # -> s = TP/(TP+FN)
                # -> p = TP/(TP+FP)
                # ESV sensitivity and precision
                TP = row[f'{test} reference'] + row[f'{test} shared']
                FN = row[f'{test} reference']
                FP = row[f'{test} nanopore']
                sensitivity = round(TP/(TP+FN),2)
                precision = round(TP/(TP+FP), 2)
                res_row.extend([sensitivity,precision])
            res.append(res_row)

        df3 = pd.DataFrame(res, columns=['ESVs sensitivity', 'ESVs precision', 'Species sensitivity', 'Species precision', 'Genus sensitivity', 'Genus precision', 'Family sensitivity', 'Family precision'])
        df4 = pd.concat([df2, df3], axis=1)
        df4.to_excel('/Users/tillmacher/Desktop/APSCALE_projects/test_dataset_apscale_nanopore/10_benchmark/benchmark_results.xlsx', index=False)

