import argparse
import glob
import math
import multiprocessing
import os
import statistics
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
from Bio.Seq import reverse_complement
from Bio.Seq import Seq
import statistics
import numpy as np
import math
import plotly.graph_objects as go

# project_folder = Path('/Users/tillmacher/Desktop/APSCALE_projects/test_dataset_apscale_nanopore/')
# settings_df = pd.read_excel('/Users/tillmacher/Desktop/APSCALE_projects/test_dataset_apscale_nanopore/test_dataset_settings.xlsx', sheet_name='Settings')
# demultiplexing_df = pd.read_excel('/Users/tillmacher/Desktop/APSCALE_projects/test_dataset_apscale_nanopore/test_dataset_settings.xlsx', sheet_name='Demultiplexing').fillna('')
# cpu_count = 7

def quality_report(project_folder, sub_folder, compressed):
    # Set parameters
    sub_folder = '2_index_demultiplexing'
    folder = project_folder / sub_folder / 'data'

    # Get all .fastq or .fastq.gz files
    suffix = '*.fastq.gz' if compressed else '*.fastq'
    fastq_files = folder.glob(suffix)

    def analyse_read(filename: Path, compressed: bool) -> dict:
        """Analyze read qualities and lengths and create summary plots."""

        compressed = True
        filename = Path('/Volumes/Coruscant/APSCALE_projects/naturalis_test_apscale_nanopore/3_primer_trimming/data/naturalis_sample_1_trimmed.fastq.gz')
        outdir = Path('/Volumes/Coruscant/APSCALE_projects/naturalis_test_apscale_nanopore/9_nanopore_report/')

        # Choose appropriate open method
        open_func = gzip.open if compressed else open

        mean_qualities = []
        read_lengths = []

        with open_func(filename, "rt") as handle:
            for record in SeqIO.parse(handle, "fastq"):
                phred_scores = record.letter_annotations["phred_quality"]
                mean_qualities.append(math.ceil(np.mean(phred_scores)))
                read_lengths.append(len(record.seq))

        # -----------------------------
        # Plot 1: Mean Phred Score Distribution
        # -----------------------------
        fig1 = go.Figure()
        x_vals = list(range(1, 61))
        y_vals = [mean_qualities.count(i) for i in x_vals]
        fig1.add_trace(go.Scatter(x=x_vals, y=y_vals, mode='markers', name='Mean Quality'))
        fig1.update_layout(title='Distribution of Mean Phred Scores',
                           xaxis_title='Mean Phred Score',
                           yaxis_title='Read Count')
        fig1.write_image(outdir / 'mean_phred_distribution.pdf')

        # -----------------------------
        # Plot 2: Read Length Distribution
        # -----------------------------
        fig2 = go.Figure()
        max_len = max(read_lengths)
        x_vals = list(range(0, max_len + 1))
        y_vals = [read_lengths.count(i) for i in x_vals]
        fig2.add_trace(go.Scatter(x=x_vals, y=y_vals, mode='markers', name='Read Length'))
        fig2.update_layout(title='Read Length Distribution',
                           xaxis_title='Read Length (bp)',
                           yaxis_title='Read Count')
        fig2.write_image(outdir / 'read_length_distribution.pdf')

        # -----------------------------
        # Plot 3: Mean Quality vs. Read Length
        # -----------------------------
        fig3 = go.Figure()
        fig3.add_trace(go.Scatter(
            x=read_lengths,
            y=mean_qualities,
            mode='markers',
            name='Length vs Quality',
            marker=dict(size=4, opacity=0.5)
        ))
        fig3.update_layout(title='Mean Quality vs. Read Length',
                           xaxis_title='Read Length (bp)',
                           yaxis_title='Mean Phred Quality')
        fig3.write_image(outdir / 'quality_vs_length.pdf')

    # Run in parallel
    res = Parallel(n_jobs=cpu_count, backend='loky')(delayed(calculate_phred_score)(fastq_file, compressed) for fastq_file in fastq_files)

def is_file_still_writing(filepath, wait_time=1.0):
    initial_size = os.path.getsize(filepath)
    time.sleep(wait_time)
    current_size = os.path.getsize(filepath)
    return initial_size != current_size

def open_file(filepath):
    system = platform.system()
    if system == 'Windows':
        os.startfile(filepath)
    elif system == 'Darwin':  # macOS
        subprocess.run(['open', filepath])
    else:  # Linux and others
        subprocess.run(['xdg-open', filepath])

def create_project(project_folder, project_name):

    # Create subfolders
    sub_folders = ['1_raw_data', '2_index_demultiplexing', '3_tag_demultiplexing', '4_primer_trimming', '5_quality_filtering', '6_denoising', '7_ESV_table', '8_taxonomic_assignment', '9_nanopore_report']
    for folder in sub_folders:
        folder_path = project_folder.joinpath(folder)
        os.makedirs(folder_path, exist_ok=True)
        data_folder_path = folder_path.joinpath('data')
        os.makedirs(data_folder_path, exist_ok=True)
        print(f'{datetime.now().strftime("%H:%M:%S")} - Created "{folder}" folder.')

    # Define path to settings file
    settings_file = project_folder.joinpath(project_name + '_settings.xlsx')

    # Create demultiplexing sheet
    cols = ['Forward index 5-3', 'Forward tag 5-3', 'Forward primer 5-3', 'Forward index 5-3', 'Forward tag 5-3', 'Forward primer 5-3', 'ID']
    rows = [['CTGT', '', 'AAACTCGTGCCAGCCACC', 'GTCCTA', '', 'GGGTATCTAATCCCAGTTTG', 'example_1_only_barcodes']]
    demultipexing_df_empty = pd.DataFrame(rows, columns=cols)

    # Create settings sheet
    cols = ['Step', 'Category', 'Variable', 'Comment']
    rows = [['General', 'cpu count', multiprocessing.cpu_count()-1, 'Number of cores to use'],
             ['demultiplexing (index)',
              'allowed errors index',
              3,
              'Allowed errors during index demultiplexing'],
            ['primer trimming',
             'allowed errors primer',
             4,
             'Allowed errors during primer trimming'],
            ['demultiplexing (tag)',
             'allowed errors tag',
             1,
             'Allowed errors during tag demultiplexing'],
             ['quality filtering',
              'minimum length',
              '',
              'Reads below this length will be discarded'],
             ['quality filtering',
              'maximum length',
              '',
              'Reads above this length will be discarded'],
            ['quality filtering',
             'minimum quality',
             20,
             'Reads below this average PHRED quality score will be discarded'],
             ['swarm denoising', 'd', 1, 'Stringency of denoising'],
            ['read table', 'minimum reads', 10, 'Discard reads below this threshold'],
            ['taxonomic assignment', 'apscale blast', 'yes', 'Run apscale megablast (yes or no)'],
            ['taxonomic assignment', 'apscale db', '', 'Path to local database'],
            ]
    settings_df_empty = pd.DataFrame(rows, columns=cols)

    # Write to multiple sheets
    with pd.ExcelWriter(settings_file, engine='openpyxl') as writer:
        demultipexing_df_empty.to_excel(writer, sheet_name='Demultiplexing', index=False)
        settings_df_empty.to_excel(writer, sheet_name='Settings', index=False)

    print(f'{datetime.now().strftime("%H:%M:%S")} - Created settings file.')
    print('')
    print(f'{datetime.now().strftime("%H:%M:%S")} - Copy your data into the "1_raw_data/data" folder.')
    print(f'{datetime.now().strftime("%H:%M:%S")} - Adjust the settings file.')
    res = input('Open in settings file Excel? (y/n): ')
    if res.upper() == 'Y':
        open_file(settings_file)
    print(f'{datetime.now().strftime("%H:%M:%S")} - Then run:')
    print(f'          $ apscale_nanopore run -p {project_name}_apscale_nanopore')
    print('')

def watch_folder(project_folder, settings_df, demultiplexing_df, live_calling, steps):

    try:
        while True:
            # Define folders
            raw_data_folder = project_folder.joinpath('1_raw_data', 'data')
            raw_tmp_folder = project_folder.joinpath('1_raw_data', 'tmp')
            os.makedirs(raw_tmp_folder, exist_ok=True)

            # Scan for files
            print(f'{datetime.now().strftime("%H:%M:%S")} - Scanning for files...')
            main_files = [i for i in glob.glob(str(raw_data_folder.joinpath('*.fastq*')))]
            main_files = {Path(file).name:Path(file) for file in main_files}
            batch = 0

            # Collect number of available CPUs
            cpu_count = settings_df[settings_df['Category'] == 'cpu count']['Variable'].values.tolist()[0]

            # Gzip files if required
            for name, file in main_files.items():
                suffix = file.suffix
                if suffix == ".fastq":
                    print(f'{datetime.now().strftime("%H:%M:%S")} - Zipping {name}...')
                    file_gz = Path(str(file) + '.gz')
                    with open(file, 'rb') as f_in:
                        with gzip.open(file_gz, 'wb') as f_out:
                            shutil.copyfileobj(f_in, f_out)
                    time.sleep(0.1)
                    os.remove(file)
                    main_files[name] = file_gz

            # Sleep if no files are present
            if len(main_files) == 0:
                print(f'{datetime.now().strftime("%H:%M:%S")} - Could not find any files to process! Waiting for new files...')
                time.sleep(2)

            # Analyse files if present
            else:
                print(f'{datetime.now().strftime("%H:%M:%S")} - Found {len(main_files)} file(s) to process!\n')
                batch += 1

                # Analyse the files
                i = 0
                for name, main_file in main_files.items():
                    # Check if file is still being written
                    while is_file_still_writing(main_file):
                        print("Waiting for file to finish writing...")
                        time.sleep(1)

                    # Start processing of the file
                    name = name.replace('.fastq.gz', '')
                    main_file = Path(main_file)
                    print(f'{datetime.now().strftime("%H:%M:%S")} - Starting analysis for: {name} ({i+1}/{len(main_files)})')

                    #=======# Index demultiplexing #=======#
                    if "Index demultiplexing" in steps:
                        print(f'{datetime.now().strftime("%H:%M:%S")} - Starting cutadapt index demultiplexing...')
                        cutadapt_index_demultiplexing(project_folder, main_file, settings_df, demultiplexing_df)
                        print(f'{datetime.now().strftime("%H:%M:%S")} - Finished cutadapt index demultiplexing!')
                        print('')

                    #=======# Tag demultiplexing #=======#
                    if "Tag demultiplexing" in steps:
                        print(f'{datetime.now().strftime("%H:%M:%S")} - Starting cutadapt tag demultiplexing...')
                        cutadapt_tag_demultiplexing(project_folder, settings_df, demultiplexing_df)
                        print(f'{datetime.now().strftime("%H:%M:%S")} - Finished cutadapt tag demultiplexing!')
                        print('')

                    #=======# Primer trimming #=======#
                    if "Primer trimming" in steps:
                        print(f'{datetime.now().strftime("%H:%M:%S")} - Starting cutadapt primer trimming...')
                        fastq_files = glob.glob(str(project_folder.joinpath('3_tag_demultiplexing', 'data', '*.fastq')))
                        Parallel(n_jobs=cpu_count, backend='threading')(delayed(cutadapt_primer_trimming)(project_folder, fastq_file, settings_df, demultiplexing_df) for fastq_file in fastq_files)
                        print(f'{datetime.now().strftime("%H:%M:%S")} - Finished cutadapt primer trimming!')
                        print('')

                    #=======# Quality filtering #=======#
                    if "Quality filtering" in steps:
                        print(f'{datetime.now().strftime("%H:%M:%S")} - Starting quality filtering...')
                        fastq_files = glob.glob(str(project_folder.joinpath('4_primer_trimming', 'data', '*.fastq.gz')))
                        Parallel(n_jobs=cpu_count, backend='loky')(delayed(python_quality_filtering)(project_folder, fastq_file, settings_df) for fastq_file in fastq_files)
                        print(f'{datetime.now().strftime("%H:%M:%S")} - Finished vsearch quality filtering!')
                        print('')

                    #=======# Denoising #=======#
                    if "Denoising" in steps:
                        print(f'{datetime.now().strftime("%H:%M:%S")} - Starting swarm denoising...')
                        fasta_files = glob.glob(str(project_folder.joinpath('5_quality_filtering', 'data', '*.fasta')))
                        Parallel(n_jobs=cpu_count, backend='loky')(delayed(swarm_denoising)(project_folder, fasta_file, settings_df) for fasta_file in fasta_files)
                        print(f'{datetime.now().strftime("%H:%M:%S")} - Finished swarm denoising...')
                        print('')

                    #=======# Read table #=======#
                    if "ESV table" in steps:
                        print(f'{datetime.now().strftime("%H:%M:%S")} - Starting to build read table...')
                        create_read_table(project_folder, settings_df)
                        print(f'{datetime.now().strftime("%H:%M:%S")} - Finished building read table!')
                        print('')

                    #=======# Taxonomic assignment #=======#
                    if "Tax. assignment" in steps:
                        print(f'{datetime.now().strftime("%H:%M:%S")} - Starting taxonomic assignment...')
                        apscale_taxonomic_assignment(project_folder, settings_df)
                        print(f'{datetime.now().strftime("%H:%M:%S")} - Finished building read table!')
                        print('')

                    if len(steps) == 7:
                        # Move file to finish analysis for the file
                        new_file = Path(str(raw_tmp_folder.joinpath(name)) + '.fastq.gz')
                        shutil.move(main_file, new_file)
                        print(f'{datetime.now().strftime("%H:%M:%S")} - Moved {name}...')

                    # Finish file
                    print(f'{datetime.now().strftime("%H:%M:%S")} - Finished analysis for: {name}\n')
                    time.sleep(1)
                    i += 1

                # =======# Create report #=======#
                # create_report(project_folder)

                print(f'{datetime.now().strftime("%H:%M:%S")} - Finished analysis for: {name} ({i}/{len(main_files)})')
                print('')

            if live_calling == False:
                break

    except KeyboardInterrupt:
        print('Stopping apscale nanopore live processing.')

def cutadapt_index_demultiplexing(project_folder, main_file, settings_df, demultiplexing_df):

    # Preprare output files
    # main_file = Path("/Users/tillmacher/Desktop/APSCALE_projects/test_dataset_apscale_nanopore/1_raw_data/data/merged_nanopore_data.fastq.gz")
    input_file = main_file
    name = input_file.name.replace('.fastq.gz', '')
    output_folder_tmp = project_folder.joinpath('2_index_demultiplexing', 'tmp')
    output_folder_data = project_folder.joinpath('2_index_demultiplexing', 'data')
    output_file = output_folder_tmp.joinpath("{name}_fwd.fastq")

    # Create tmp folder
    tmp_folder = project_folder.joinpath('2_index_demultiplexing', 'tmp')
    os.makedirs(tmp_folder, exist_ok=True)

    # Create reverse complement of untrimmed
    output_file_rc = output_folder_tmp.joinpath("{name}_rc.fastq")
    untrimmed_folder = project_folder.joinpath('2_index_demultiplexing', 'untrimmed')
    os.makedirs(untrimmed_folder, exist_ok=True)
    untrimmed_fastq = untrimmed_folder.joinpath('untrimmed.fastq')
    untrimmed_rc_fastq = untrimmed_folder.joinpath('untrimmed_rc.fastq')

    # Collect required settings
    number_of_errors = settings_df[settings_df['Category'] == 'allowed errors index']['Variable'].values.tolist()[0]
    cpu_count = settings_df[settings_df['Category'] == 'cpu count']['Variable'].values.tolist()[0]

    ##======## Demultuplexing forward reads ##======##
    # Run cutadapt demultiplexing
    g_args = []
    for _, row in demultiplexing_df.iterrows():
        # Create forward sequence
        fwd_seq = row['Forward index 5-3']
        # Create reverse sequence
        rvs_seq = reverse_complement(row['Forward index 5-3'])
        # Combine to search sequence
        search_seq = f'{fwd_seq}...{rvs_seq}'
        g_args.extend(['-g', search_seq])

    # Run cutadapt demultiplexing
    command = f"cutadapt -e {number_of_errors} {' '.join(g_args)} --cores {cpu_count} -o {output_file} --untrimmed-output {untrimmed_fastq} --report=minimal {input_file}"
    process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    stdout, stderr = process.communicate()
    in_reads1 = int(stdout.split()[11])
    out_reads1 = int(stdout.split()[-3])

    print(f'{datetime.now().strftime("%H:%M:%S")} - Finished demultiplexing in 5\'-3\' orientation!')

    ##======## Demultuplexing RC reads ##======##
    # Vsearch reverse complement
    command = f"vsearch --fastx_revcomp {untrimmed_fastq} --fastqout {untrimmed_rc_fastq}"
    process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    stdout, stderr = process.communicate()

    if untrimmed_fastq.exists():
        os.remove(untrimmed_fastq)

    # Run cutadapt again
    command = f"cutadapt -e {number_of_errors} {' '.join(g_args)} --cores {cpu_count} -o {output_file_rc} --discard-untrimmed --report=minimal {untrimmed_rc_fastq}"
    process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    stdout, stderr = process.communicate()
    out_reads2 = int(stdout.split()[-3])

    print(f'{datetime.now().strftime("%H:%M:%S")} - Finished demultiplexing in 3\'-5\' orientation!')

    try:
        reads_perc = round((out_reads1 + out_reads2) / in_reads1 * 100, 2)
    except ZeroDivisionError:
        reads_perc = 0
    print(f'{datetime.now().strftime("%H:%M:%S")} - Finished demultiplexing of {name}: {in_reads1:,} -> {out_reads1 + out_reads2:,} reads ({reads_perc}% passed).')

    if untrimmed_folder.exists():
        shutil.rmtree(untrimmed_folder)

    ##======## Merge files ##======##
    # Collect all .fastq files (uncompressed)
    demultiplexed_fwd_files = sorted(glob.glob(str(output_folder_tmp.joinpath('*_fwd.fastq'))))
    demultiplexed_rc_files = sorted(glob.glob(str(output_folder_tmp.joinpath('*_rc.fastq'))))

    # Only continue if files were demultiplexed
    if len(demultiplexed_fwd_files) == 0:
        print(f'{datetime.now().strftime("%H:%M:%S")} - Error: Could not find any demultiplexed files!')
        return

    # Define merging function (uncompressed)
    def merge_fastq_files(tmp_file_fwd, tmp_file_rc, output_folder_data):
        index = int(Path(tmp_file_fwd).name.replace('_fwd.fastq', '')) - 1
        sample_id = demultiplexing_df['ID'][index]
        output_file = output_folder_data.joinpath(f'{sample_id}.fastq')

        with open(output_file, "w") as out_handle:
            for in_file in [tmp_file_fwd, tmp_file_rc]:
                with open(in_file, "r") as in_handle:
                    shutil.copyfileobj(in_handle, out_handle, length=1024 * 1024)  # 1MB buffer

        os.remove(tmp_file_fwd)
        os.remove(tmp_file_rc)
        print(f'{datetime.now().strftime("%H:%M:%S")} - Merged and saved {sample_id}!')

    Parallel(n_jobs=-1, backend='loky')(delayed(merge_fastq_files)(tmp_file_fwd, tmp_file_rc, output_folder_data) for tmp_file_fwd, tmp_file_rc in zip(demultiplexed_fwd_files, demultiplexed_rc_files))

def cutadapt_tag_demultiplexing(project_folder, settings_df, demultiplexing_df):

    # Check if tags are used
    print(f'{datetime.now().strftime("%H:%M:%S")} - Samples will be copied to the "3_tag_demultiplexing" folder.')
    # Still copy files to respective folder
    files = glob.glob(str(project_folder.joinpath('2_index_demultiplexing', 'data', '*.fastq')))
    for file in files:
        new_file = file.replace('2_index_demultiplexing', '3_tag_demultiplexing')
        shutil.move(file, new_file)

    if '' in demultiplexing_df['Forward tag 5-3'].values.tolist():
        print(f'{datetime.now().strftime("%H:%M:%S")} - No tagging information was found - skipping tag demultiplexing.')
        return

    # Preprare output files
    output_folder_tmp = project_folder.joinpath('3_tag_demultiplexing', 'tmp')
    output_folder_data = project_folder.joinpath('3_tag_demultiplexing', 'data')
    output_file = output_folder_tmp.joinpath("{name}.fastq")

    # Create tmp folder
    tmp_folder = project_folder.joinpath('3_tag_demultiplexing', 'tmp')
    os.makedirs(tmp_folder, exist_ok=True)

    # Collect required settings
    number_of_errors = settings_df[settings_df['Category'] == 'allowed errors tag']['Variable'].values.tolist()[0]
    cpu_count = settings_df[settings_df['Category'] == 'cpu count']['Variable'].values.tolist()[0]

    # Collect all files
    files = glob.glob(str(project_folder.joinpath('4_primer_trimming', 'data', '*.fastq.gz')))
    merged_file = project_folder.joinpath('3_tag_demultiplexing', 'tmp', 'merged_samples.fastq.gz')
    with gzip.open(merged_file, "wb") as outfile:
        for f in files:
            with gzip.open(f, "rb") as infile:
                shutil.copyfileobj(infile, outfile)

    # Run cutadapt demultiplexing
    g_args = []
    for _, row in demultiplexing_df.iterrows():
        # Create forward sequence
        fwd_seq = row['Forward tag 5-3']
        # Create reverse sequence
        rvs_seq = row['Forward tag 5-3']
        # Combine to search sequence
        search_seq = f'{fwd_seq}...{rvs_seq}'
        g_args.extend(['-g', search_seq])

    # Run cutadapt demultiplexing and primer trimming
    command = f"cutadapt -e {number_of_errors} {' '.join(g_args)} --cores 1 -o {output_file} --discard-untrimmed --report=minimal {merged_file}"
    process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    stdout, stderr = process.communicate()

    # You can now use `stdout` and `stderr` as variables
    in_reads = int(stdout.split()[11])
    out_reads = int(stdout.split()[-3])
    out_reads_perc = round(out_reads / in_reads * 100, 2)
    print(f'{datetime.now().strftime("%H:%M:%S")} - Finished tag demultiplexing of: {out_reads_perc}% of reads kept.')

    # Collect all .fastq files (uncompressed)
    demultiplexed_files = glob.glob(str(output_folder_tmp.joinpath('*.fastq')))

    if len(demultiplexed_files) == 0:
        print(f'{datetime.now().strftime("%H:%M:%S")} - Error: Could not find any demultiplexed files!')
        return False

    for tmp_file in tqdm(demultiplexed_files, desc='Writing sample files'):
        tmp_file = Path(tmp_file)
        index = int(tmp_file.name.replace('.fastq', '')) - 1
        sample_id = demultiplexing_df['ID'][index]
        sample_file = output_folder_data.joinpath(f'{sample_id}.fastq.gz')

        # Append compressed data to the .fastq.gz file
        with open(tmp_file, "rb") as in_handle, gzip.open(sample_file, "ab") as out_handle:
            shutil.copyfileobj(in_handle, out_handle)

    # Remove tmp files
    shutil.rmtree(output_folder_tmp)

def cutadapt_primer_trimming(project_folder, file, settings_df, demultiplexing_df):

    # Preprare output files
    # file = '/Users/tillmacher/Desktop/APSCALE_projects/test_dataset_apscale_nanopore/2_index_demultiplexing/data/Sample_3.fastq'
    input_file = Path(file)
    name = input_file.name.replace('.fastq', '')
    output_folder_data = project_folder.joinpath('4_primer_trimming', 'data')
    output_file = output_folder_data.joinpath(f"{name}_trimmed.fastq.gz")

    # Also save the untrimmed reads
    # Create reverse complement of untrimmed
    untrimmed_folder = project_folder.joinpath('4_primer_trimming', 'untrimmed')
    os.makedirs(untrimmed_folder, exist_ok=True)
    output_file_untrimmed = untrimmed_folder.joinpath(f"{name}_untrimmed.fastq")
    output_file_untrimmed_rc = untrimmed_folder.joinpath(f"{name}_untrimmed_rc.fastq")
    output_file_rc = output_folder_data.joinpath(f"{name}_trimmed_rc.fastq.gz")
    tmp_file = untrimmed_folder.joinpath(f"{name}_tmp.fastq.gz")

    # Collect required settings
    number_of_errors = settings_df[settings_df['Category'] == 'allowed errors primer']['Variable'].values.tolist()[0]

    # Run cutadapt demultiplexing
    # Create forward sequence
    sub_df = demultiplexing_df[demultiplexing_df['ID'] == name]
    fwd_seq = sub_df['Forward primer 5-3'].values.tolist()[0]
    # Create reverse sequence
    rvs_seq_rc = reverse_complement(sub_df['Reverse primer 5-3'].values.tolist()[0])
    adapter = f'{fwd_seq}...{rvs_seq_rc}'

    ##======## Trimming of reads in 5'-3' orientation ##======##
    # Run cutadapt demultiplexing and primer trimming
    command = f"cutadapt -e {number_of_errors} -g {adapter} --cores 1 -o {output_file} --untrimmed-output {output_file_untrimmed} --report=minimal {input_file}"
    process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    stdout, stderr = process.communicate()
    # You can now use `stdout` and `stderr` as variables
    in_reads1 = int(stdout.split()[11])
    out_reads1 = int(stdout.split()[-3])

    ##======## Trimming of reads in 3'-5' orientation ##======##
    # Vsearch reverse complement
    command = f"vsearch --fastx_revcomp {output_file_untrimmed} --fastqout {output_file_untrimmed_rc}"
    process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    stdout, stderr = process.communicate()

    if output_file_untrimmed.exists():
        os.remove(output_file_untrimmed)

    # Run cutadapt primer trimming
    command = f"cutadapt -e {number_of_errors} -g {adapter} --cores 1 -o {output_file_rc} --discard-untrimmed --report=minimal {output_file_untrimmed_rc}"
    process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    stdout, stderr = process.communicate()
    # You can now use `stdout` and `stderr` as variables
    in_reads2 = int(stdout.split()[11])
    out_reads2 = int(stdout.split()[-3])

    # Combine and overwrite
    with gzip.open(tmp_file, 'wb') as out_f:
        with gzip.open(output_file, 'rb') as f1:
            shutil.copyfileobj(f1, out_f)
        with gzip.open(output_file_rc, 'rb') as f2:
            shutil.copyfileobj(f2, out_f)

    # Replace original file with combined content
    os.replace(tmp_file, output_file)
    # Remove tmp files
    if output_file_untrimmed_rc.exists():
        os.remove(output_file_untrimmed_rc)
    if output_file_rc.exists():
        os.remove(output_file_rc)

    try:
        reads_perc = round((out_reads1 + out_reads2) / in_reads1 * 100, 2)
    except ZeroDivisionError:
        reads_perc = 0
    print(f'{datetime.now().strftime("%H:%M:%S")} - {name}: {in_reads1:,} -> {out_reads1 + out_reads2:,} reads ({reads_perc}% passed).')

def python_quality_filtering(project_folder, file, settings_df):

    # Preprare output files
    # file = '/Users/tillmacher/Desktop/APSCALE_projects/test_dataset_apscale_nanopore/3_primer_trimming/data/Sample_3_trimmed.fastq.gz'
    input_file = Path(file)
    name = input_file.name.replace('.fastq.gz', '')
    output_folder_data = project_folder.joinpath('5_quality_filtering', 'data')
    # output files
    filtered_fasta = output_folder_data.joinpath(f'{name}_filtered.fasta')
    dereplicated_fasta = output_folder_data.joinpath(f'{name}_filtered_derep.fasta')

    # Collect required settings
    min_len = settings_df[settings_df['Category'] == 'minimum length']['Variable'].values.tolist()[0]
    max_len = settings_df[settings_df['Category'] == 'maximum length']['Variable'].values.tolist()[0]
    trunc_val = settings_df[settings_df['Category'] == 'minimum quality']['Variable'].values.tolist()[0]

    # Run python-based quality filtering
    reads_1 = 0
    total_reads = 0

    with gzip.open(input_file, "rt") as in_handle, open(filtered_fasta, "w") as out_handle:
        for record in SeqIO.parse(in_handle, "fastq"):
            total_reads += 1
            phred_scores = record.letter_annotations["phred_quality"]
            read_length = len(record.seq)
            if not phred_scores:
                continue
            if np.mean(phred_scores) >= trunc_val and min_len <= read_length <= max_len:
                fasta_record = SeqRecord(
                    Seq(str(record.seq)),
                    id=record.id,
                    description=""
                )
                SeqIO.write(fasta_record, out_handle, "fasta")
                reads_1 += 1

    # Run vsearch dereplication
    command = f"vsearch --threads 1 --derep_fulllength {filtered_fasta} --sizeout --relabel_sha1 --output {dereplicated_fasta}"
    process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    stdout, stderr = process.communicate()
    reads_2 = stderr.split()[-15]

    try:
        reads_perc = round(reads_1 / total_reads * 100, 2)
    except ZeroDivisionError:
        reads_perc = 0
    print(f'{datetime.now().strftime("%H:%M:%S")} - {name.replace("_trimmed", "")}: {total_reads:,} -> {reads_1:,} ({reads_perc}% passed) -> {int(reads_2):,} (dereplication)')

    if filtered_fasta.exists():
        os.remove(filtered_fasta)

def swarm_denoising(project_folder, file, settings_df):

    # Preprate output files
    input_file = Path(file)
    name = input_file.name.replace('_trimmed_filtered_derep.fasta', '')
    output_folder_data = project_folder.joinpath('6_denoising', 'data')
    cluster_file = output_folder_data.joinpath(f'{name}_clusters.fasta')

    # Collect required settings
    d_value = settings_df[settings_df['Category'] == 'd']['Variable'].values.tolist()[0]

    # Run swarm denoising
    command = f"swarm -d {d_value} --threads 1 -z --seeds {cluster_file} {input_file}"
    process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    stdout, stderr = process.communicate()

    # You can now use `stdout` and `stderr` as variables
    swarms = int(stderr.split()[-7])
    print(f'{datetime.now().strftime("%H:%M:%S")} - {name}: {swarms:,} swarms.')

def create_read_table(project_folder, settings_df):
    # Prepare output files
    swarm_files_path = project_folder.joinpath('6_denoising', 'data', '*_clusters.fasta')
    swarm_files = [Path(i) for i in glob.glob(str(swarm_files_path))]
    data = defaultdict(lambda: defaultdict(int))  # nested dict: hash -> sample -> size
    seq_dict = {}  # hash -> sequence

    # Collect required settings
    min_reads = settings_df[settings_df['Category'] == 'minimum reads']['Variable'].values.tolist()[0]

    # Parse files
    for file in sorted(swarm_files):
        sample = file.name.replace('_clusters.fasta', '')
        with open(file) as handle:
            for record in SeqIO.parse(handle, "fasta"):
                size = int(record.id.split(';')[1].replace('size=', ''))
                seq = str(record.seq)
                hash = hashlib.sha3_256(seq.encode("ascii")).hexdigest()
                data[hash][sample] += size
                seq_dict[hash] = seq

    # Create and prepare DataFrame
    df = pd.DataFrame.from_dict(data, orient='index').fillna(0).astype(int)
    df.insert(0, "Seq", df.index.map(seq_dict))

    # Filter low-abundance counts
    sample_cols = df.columns.difference(['Seq'])
    total_reads_before = df[sample_cols].sum().sum()
    n_swarms_before = len(df)
    df[sample_cols] = df[sample_cols].where(df[sample_cols] >= int(min_reads), 0)

    # Filter low-abundance ESVs
    df.index.name = "ID"
    df.reset_index(inplace=True)
    df['sum'] = df[sample_cols].sum(axis=1)
    df = df[df['sum'] != 0]
    df['sum'] = df[sample_cols].sum(axis=1)
    df = df.sort_values('sum', ascending=False)

    # Calculate stats
    total_reads_after = df['sum'].sum()
    removed = total_reads_before - total_reads_after
    n_swarms = len(df)
    swarms_discarded = n_swarms_before - n_swarms
    print(f'{datetime.now():%H:%M:%S} - Final read table contains {n_swarms:,} swarms accounting for {total_reads_after:,} reads.')
    print(f'{datetime.now():%H:%M:%S} - Discarded {swarms_discarded:,} swarms accounting for {removed:,} reads (<= {min_reads} reads).')

    # Final clean-up
    df.drop(columns='sum', inplace=True)

    # insert empty files
    for file in sorted(swarm_files):
        sample = file.name.replace('_clusters.fasta', '')
        if sample not in list(df.columns):
            df[sample] = [0] * len(df)

    # Collect name of the project
    project_name = project_folder.name.replace('_apscale_nanopore', '')

    # Write to files
    if df.shape[0] < 65000:
        excel_file = project_folder.joinpath('7_ESV_table', f'{project_name}_swarms.xlsx')
        df.to_excel(excel_file, index=False)
    parquet_file = project_folder.joinpath('7_ESV_table', f'{project_name}_swarms.parquet.snappy')
    df.to_parquet(parquet_file, compression='snappy')

    # Write sequences to fasta
    fasta_file = project_folder.joinpath('7_ESV_table', 'data', f'{project_name}_swarms.fasta')
    with open(fasta_file, "w") as output_handle:
        for hash, seq in df[['ID', 'Seq']].values.tolist():
            record = SeqRecord(Seq(seq), id=hash, description='')
            SeqIO.write(record, output_handle, "fasta")

    # Perform chimera detection
    nochimera_fasta = project_folder.joinpath('7_ESV_table', 'data', f'{project_name}_swarms_nochimera.fasta')
    command = f'vsearch --uchime3_denovo {fasta_file} --nonchimeras {nochimera_fasta}'
    process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    stdout, stderr = process.communicate()
    nochimera_total = stderr.split()[-12]
    nochimera_perc = stderr.split()[-11]
    print(f'{datetime.now().strftime("%H:%M:%S")} - Wrote {int(nochimera_total):,} {nochimera_perc} non-chimeras.')
    time.sleep(1)

def apscale_taxonomic_assignment(project_folder, settings_df):
    # Define files
    project_name = project_folder.name.replace('_apscale_nanopore', '')
    fasta_file = project_folder.joinpath('7_ESV_table', 'data', f'{project_name}_swarms_nochimera.fasta')
    results_folder = project_folder.joinpath('8_taxonomic_assignment', project_name)

    # Collect variables
    run_blastn = settings_df[settings_df['Category'] == 'apscale blast']['Variable'].values.tolist()[0]
    blastn_db = settings_df[settings_df['Category'] == 'apscale db']['Variable'].values.tolist()[0]

    # Run apscale blast
    if run_blastn == 'yes':
        try:
            shutil.rmtree(results_folder)
            os.makedirs(results_folder, exist_ok=True)
        except FileNotFoundError:
            pass
        command = f"apscale_blast -db {blastn_db} -q {fasta_file} -o {results_folder} -task megablast"
        process = subprocess.Popen(command, shell=True, text=True)
        process.wait()

def create_report(project_folder):
    # collect information
    project_name = project_folder.name.replace('_apscale_nanopore', '')
    read_table_file = project_folder.joinpath('7_ESV_table', f'{project_name}_swarms.parquet.snappy')
    read_table_df = pd.read_parquet(read_table_file).fillna('')
    taxonomy_table_file = project_folder.joinpath('8_taxonomic_assignment', project_name, f'{project_name}_taxonomy.xlsx')
    taxonomy_table_df = pd.read_excel(taxonomy_table_file).fillna('')

    # create folder to store batch results
    report_folder = project_folder.joinpath('9_nanopore_report')
    os.makedirs(report_folder, exist_ok=True)

    # Count previous runs
    previous_runs = glob.glob(str(report_folder.joinpath('batch_*.xlsx')))
    batch = len(previous_runs) + 1

    # Create output files
    output_xlsx = report_folder.joinpath(f'batch_{batch}.xlsx')
    output_pdf = report_folder.joinpath(f'batch_{batch}.pdf')
    output_html = report_folder.joinpath(f'batch_{batch}.html')

    # Collect stats
    n_species = len(set([i for i in taxonomy_table_df['Species'].values.tolist() if i != 'NoMatch' and i != '']))
    n_genera = len(set([i for i in taxonomy_table_df['Genus'].values.tolist() if i != 'NoMatch' and i != '']))
    n_families = len(set([i for i in taxonomy_table_df['Family'].values.tolist() if i != 'NoMatch' and i != '']))
    sample_df = read_table_df[read_table_df.columns.tolist()[2:]]
    n_reads_total = sample_df.sum().sum()
    mean_reads_per_sample = round(sample_df.sum().mean(), 2)
    sheet1 = pd.DataFrame([n_species, n_genera, n_families, n_reads_total, mean_reads_per_sample]).transpose()
    sheet1.columns = ['Mean species', 'Mean genera', 'Mean families', 'Reads (total)', 'Reads per sample (mean)']
    sheet2 = sample_df.describe().transpose()

    # Write to multiple sheets
    with pd.ExcelWriter(output_xlsx, engine="openpyxl") as writer:
        sheet1.to_excel(writer, sheet_name="Sheet1", index=False)
        sheet2.to_excel(writer, sheet_name="Sheet2", index=False)

def main():
    """
    APSCALE nanopore suite
    Command-line tool to process nanopore sequence data.
    """

    # Introductory message with usage examples
    message = """
    APSCALE nanopore command line tool - v0.0.1
    Example commands:
    $ apscale_nanopore create my_new_project
    $ apscale_nanopore run my_new_project

    """
    print(message)

    # Initialize main parser
    parser = argparse.ArgumentParser(description='APSCALE nanopore v0.0.1')
    subparsers = parser.add_subparsers(dest='command', required=True)

    # === Subparser: create ===
    create_parser = subparsers.add_parser('create', help='Create a new APSCALE nanopore project.')
    create_parser.add_argument('-p', '--project', type=str, required=True, help='Path to project.')

    # === Subparser: run ===
    run_parser = subparsers.add_parser('run', help='Run the APSCALE nanopore pipeline.')
    run_parser.add_argument('-p', '--project', type=str, required=True, help='Path to project.')
    run_parser.add_argument('-live', '--live_calling', action='store_true', help='Scan 1_raw_data for new batches.')
    run_parser.add_argument('-e1', type=str, help='Overwrite: allowed index demultiplexing errors.')
    run_parser.add_argument('-e2', type=str, help='Overwrite: allowed primer trimming errors.')
    run_parser.add_argument('-e3', type=str, help='Overwrite: allowed tag demultiplexing errors.')
    run_parser.add_argument('-t', type=str, help='Overwrite: quality truncation cutoff.')
    run_parser.add_argument('-minlen', type=str, help='Overwrite: minimum length.')
    run_parser.add_argument('-maxlen', type=str, help='Overwrite: maximum length.')
    run_parser.add_argument('-minq', type=str, help='Overwrite: minimum quality value.')
    run_parser.add_argument('-d', type=str, help="Overwrite: swarm's d value.")
    run_parser.add_argument('-minreads', type=str, help="Overwrite: Read filter threshold.")
    run_parser.add_argument('-step', type=str, help="Select step to re-run individually. "
                                                    "1:Index demultiplexing, 2:Primer trimming, "
                                                    "3:Tag demultiplexing, 4:Quality filtering, "
                                                    "5:Denoising, 6:ESV table, 7:Tax. assignment ")
    run_parser.add_argument('-steps', type=str, help="Select step from which to re-run all subsequent steps. "
                                                    "1:Index demultiplexing, 2:Primer trimming, "
                                                    "3:Tag demultiplexing, 4:Quality filtering, "
                                                    "5:Denoising, 6:ESV table, 7:Tax. assignment ")

    # Parse arguments
    args = parser.parse_args()

    # Create project
    if args.command == 'create':
        project_folder = Path(str(Path(args.project)) + '_apscale_nanopore')
        project_name = project_folder.name.replace('_apscale_nanopore', '')
        create_project(project_folder, project_name)

    # Run apscale
    elif args.command == 'run':

        # Collect step information
        all_steps = {"1": "Index demultiplexing", "2": "Tag demultiplexing", "3": "Primer trimming",
                     "4": "Quality filtering", "5": "Denoising", "6": "ESV table", "7": "Tax. assignment"}
        if args.step:
            steps = [all_steps[args.step]]
        elif args.steps:
            steps = [all_steps[str(i)] for i in range(int(args.steps), 8)]
        else:
            steps = list(all_steps.values())
        if steps == []:
            print('Error: Please choose a suitable step index!')
            return

        # Define folders and files
        project_folder = Path(args.project)
        project_name = project_folder.name.replace('_apscale_nanopore', '')
        settings_file = project_folder.joinpath(project_name + '_settings.xlsx')

        # Load settings dataframe
        if settings_file.exists():
            settings_df = pd.read_excel(settings_file, sheet_name='Settings').fillna('')
            demultiplexing_df = pd.read_excel(settings_file, sheet_name='Demultiplexing').fillna('')

            # Check if argument require to be adjusted
            if args.e1:
                index = settings_df[settings_df['Category'] == 'allowed errors index'].index[0]
                settings_df.loc[index, 'Variable'] = args.e1
                print(f'Adjusted value: Number of allowed index errors: {args.e1}')
            if args.e2:
                index = settings_df[settings_df['Category'] == 'allowed errors primer'].index[0]
                settings_df.loc[index, 'Variable'] = args.e2
                print(f'Adjusted value: Number of allowed primer errors: {args.e2}')
            if args.e3:
                index = settings_df[settings_df['Category'] == 'allowed errors tag'].index[0]
                settings_df.loc[index, 'Variable'] = args.e3
                print(f'Adjusted value: Number of allowed tag errors: {args.e3}')
            if args.t:
                index = settings_df[settings_df['Category'] == 'q_min'].index[0]
                settings_df.loc[index, 'Variable'] = args.t
                print(f'Adjusted value: Truncation cutoff: {args.t}')
            if args.minlen:
                index = settings_df[settings_df['Category'] == 'minimum length'].index[0]
                settings_df.loc[index, 'Variable'] = args.minlen
                print(f'Adjusted value: Minimum length: {args.minlen}')
            if args.maxlen:
                index = settings_df[settings_df['Category'] == 'maximum length'].index[0]
                settings_df.loc[index, 'Variable'] = args.maxlen
                print(f'Adjusted value: Maximum length: {args.maxlen}')
            if args.minq:
                index = settings_df[settings_df['Category'] == 'minimum quality'].index[0]
                settings_df.loc[index, 'Variable'] = args.minq
                print(f'Adjusted value: Maximum expected error: {args.minq}')
            if args.d:
                index = settings_df[settings_df['Category'] == 'd'].index[0]
                settings_df.loc[index, 'Variable'] = args.d
                print(f"Adjusted value: Swarm's d value: {args.d}")
            if args.minreads:
                index = settings_df[settings_df['Category'] == 'minimum reads'].index[0]
                settings_df.loc[index, 'Variable'] = args.minreads
                print(f"Adjusted value: Read filter threshold: {args.minreads}")

            # =======# Live processing #=======#
            print('')
            watch_folder(project_folder, settings_df, demultiplexing_df, args.live_calling, steps)

        else:
            print(settings_file)
            print('Error: Cannot find settings file!')

if __name__ == "__main__":
    main()