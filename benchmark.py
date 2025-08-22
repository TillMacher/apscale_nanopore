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

# BASE_PROJECT_PATH = Path('/Users/tillmacher/Desktop/APSCALE_projects/test_dataset_apscale_nanopore')
# target_len = 64

print('Enter the FULL PATH to apscale-nanopore project:')
BASE_PROJECT_PATH = input()
if not Path(BASE_PROJECT_PATH).exists():
    print('Warning: This folder does not exists!')
    sys.exit()

print('Enter the TARGET LENGTH of your fragment:')
target_len = input()
try:
    target_len = int(target_len)
except:
    print('Warning: This is not a proper target length value!')
    sys.exit()

# =====================
# === FILES ===
# =====================

# Input/Output Files
name = BASE_PROJECT_PATH.name.replace('_apscale_nanopore', '')
DEMUX_INDEX_PATH = BASE_PROJECT_PATH / '2_index_demultiplexing/data'
NANOPORE_FILE = BASE_PROJECT_PATH / f'7_taxonomic_assignment/test_dataset/{name}_taxonomy.xlsx'
READ_TABLE_FILE = BASE_PROJECT_PATH / f'6_read_table/{name}_read_table.xlsx'
BENCHMARK_FOLDER = BASE_PROJECT_PATH / '9_benchmark'
os.makedirs(BENCHMARK_FOLDER, exist_ok=True)

# =====================
# === SETTINGS LOGIC ===
# =====================

index_error_values = [i for i in range(1, 4)]
primer_error_values = [i for i in range(1, 4)]
minq = [i for i in range(10, 41, 10)]
maxmin_values = [[target_len + 10, target_len - 10]]
mode = ['ESVs', 'Swarms', 'Denoised OTUs', 'Swarm OTUs']
percid_values = [0.97, 0.98]
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
existing_ids = [Path(i).name for i in glob.glob(str(BENCHMARK_FOLDER / '*'))]

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

    # save runtime
    end_time = round(time.time() - start_time, 2)
    time_df = pd.DataFrame([end_time], columns=['time'])

    # save files
    res_folder =  BENCHMARK_FOLDER / id
    os.makedirs(res_folder, exist_ok=True)

    taxonomy_table = BASE_PROJECT_PATH / f'7_taxonomic_assignment/{name}/{name}_taxonomy.xlsx'
    dest = res_folder / f'{id}_taxonomy.xlsx'
    shutil.move(taxonomy_table, dest)

    read_table = BASE_PROJECT_PATH / f'6_read_table/{name}_read_table.xlsx'
    dest = res_folder / f'{id}_read_table.xlsx'
    shutil.move(read_table, dest)

    seq_fasta = BASE_PROJECT_PATH / f'6_read_table/data/{name}_centroids.fasta'
    dest = res_folder / f'{id}_centroids.fasta'
    shutil.move(seq_fasta, dest)

    time_file = res_folder / f'{id}_time.xlsx'
    time_df.to_excel(time_file, index=False)








