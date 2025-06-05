# Imports
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
from Bio.SeqIO import SeqRecord
from Bio.Seq import Seq
import itertools
import random

# =====================
# === USER SETTINGS ===
# =====================

# Base Project Path
BASE_PROJECT_PATH = Path('/Users/tillmacher/Desktop/APSCALE_projects/test_dataset_apscale_nanopore')

# Input/Output Files
REFERENCE_FILE = BASE_PROJECT_PATH / '9_benchmark/merged_main_new/merged_main_new_taxonomy.xlsx'
NANOPORE_FILE = BASE_PROJECT_PATH / '7_taxonomic_assignment/test_dataset/test_dataset_taxonomy.xlsx'
RESULTS_FILE = BASE_PROJECT_PATH / '9_benchmark/benchmark_results.xlsx'

# Intermediate Files & Folders
DEMUX_INDEX_PATH = BASE_PROJECT_PATH / '2_index_demultiplexing/data'
MERGED_FASTQ_SRC = BASE_PROJECT_PATH / '1_raw_data/tmp/merged_nanopore_data.fastq.gz'
MERGED_FASTQ_DST = BASE_PROJECT_PATH / '1_raw_data/data/merged_nanopore_data.fastq.gz'

# =====================
# === SETTINGS LOGIC ===
# =====================

index_error_values = [i for i in range(1, 4)]
primer_error_values = [i for i in range(1, 4)]
minq = [i for i in range(10, 41, 10)]
target_len = 64
plus_minus = [i for i in range(10, 21, 10)]
maxmin_values = [[target_len + i, target_len - i] for i in plus_minus]
mode = ['ESVs', 'Swarms', 'Denoised OTUs', 'Swarm OTUs']
percid_values = [0.97]
alpha_values = [1, 2, 3]
d_values = [1, 2, 3]
min_read_values = [2, 10]

combinations = []

# Generate valid combinations
for idx_err in index_error_values:
    for prim_err in primer_error_values:
        for q in minq:
            for maxmin in maxmin_values:
                for m in mode:
                    for min_reads in min_read_values:
                        if m == 'Swarm OTUs':
                            for percid in percid_values:
                                for d in d_values:
                                    combinations.append((idx_err, prim_err, q, maxmin[1], maxmin[0], m, percid, 1, d, min_reads))
                        elif m == 'ESVs':
                            for alpha in alpha_values:
                                combinations.append((idx_err, prim_err, q, maxmin[1], maxmin[0], m, 1, alpha, 1, min_reads))
                        elif m == 'Swarms':
                            for d in d_values:
                                combinations.append((idx_err, prim_err, q, maxmin[1], maxmin[0], m, 1, 1, d, min_reads))
                        elif m == 'Denoised OTUs':
                            for percid in percid_values:
                                for alpha in alpha_values:
                                    combinations.append((idx_err, prim_err, q, maxmin[1], maxmin[0], m, percid, alpha, 1, min_reads))

print(f'Number of test combinations: {len(combinations)}')

res = []
res_df = pd.read_excel(RESULTS_FILE).fillna('').drop_duplicates()
existing_ids = res_df['id'].tolist() if RESULTS_FILE.exists() else []

random.shuffle(combinations)

for combination in combinations:
    start_time = time.time()

    e_index_value, e_primer_value, minq_value, min_value, max_value, mode_value, percid_value, alpha_value, d_value, minreads_value = combination

    id = '-'.join([str(i) for i in combination])
    run = [id]

    if id in existing_ids:
        print(f'Combination "{id}" already exists.')
        continue

    run += list(combination)

    command = (
        f"apscale_nanopore run -p {BASE_PROJECT_PATH} "
        f"-e1 {e_index_value} -e2 {e_primer_value} -minq {minq_value} "
        f"-minlen {min_value} -maxlen {max_value} -mode \"{mode_value}\" "
        f"-percid {percid_value} -alpha {alpha_value} -d {d_value} -minreads {minreads_value}"
    )
    subprocess.run(command, shell=True, text=True)

    # Load and compare
    ref_df = pd.read_excel(REFERENCE_FILE).fillna('')
    test_df = pd.read_excel(NANOPORE_FILE).fillna('')

    def get_unique(df, col):
        return [i for i in df[col].drop_duplicates().values.tolist() if i != '']

    comparisons = []
    for col in ['unique ID', 'Species', 'Genus', 'Family']:
        ref_set = set(get_unique(ref_df, col))
        test_set = set(get_unique(test_df, col))
        comparisons += [len(ref_set - test_set), len(ref_set & test_set), len(test_set - ref_set)]

    run += comparisons
    run.append(round(time.time() - start_time, 2))
    res.append(run)

    # Cleanup
    if DEMUX_INDEX_PATH.exists():
        shutil.rmtree(DEMUX_INDEX_PATH)
        os.makedirs(DEMUX_INDEX_PATH, exist_ok=True)

    # shutil.move(str(MERGED_FASTQ_SRC), str(MERGED_FASTQ_DST))

    # Write progress
    df = pd.DataFrame(res)
    df.columns = ['id', 'index_error', 'primer_error', 'minq', 'minlen', 'maxlen', 'mode', 'percid', 'alpha', 'd', 'minreads',
                  'ESVs_ref', 'ESVs_shared', 'ESVs_only', 'Species_ref', 'Species_shared', 'Species_only',
                  'Genus_ref', 'Genus_shared', 'Genus_only', 'Family_ref', 'Family_shared', 'Family_only', 'elapsed_time']
    df_comb = pd.concat([df, res_df], ignore_index=True)

    df_comb.to_excel(RESULTS_FILE, index=False)





