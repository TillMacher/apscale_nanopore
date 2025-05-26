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

# project_folder = Path('/Users/tillmacher/Desktop/APSCALE_projects/test_dataset_apscale_nanopore')
# settings_df = pd.read_excel('/Users/tillmacher/Desktop/APSCALE_projects/test_dataset_apscale_nanopore/test_dataset_settings.xlsx', sheet_name='Settings')
# demultiplexing_df = pd.read_excel('/Users/tillmacher/Desktop/APSCALE_projects/test_dataset_apscale_nanopore/test_dataset_settings.xlsx', sheet_name='Demultiplexing')

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
    sub_folders = ['1_raw_data', '2_demultiplexing', '3_quality_filtering', '4_denoising', '5_ESV_table', '6_taxonomic_assignment', '7_nanopore_report']
    for folder in sub_folders:
        folder_path = project_folder.joinpath(folder)
        os.makedirs(folder_path, exist_ok=True)
        data_folder_path = folder_path.joinpath('data')
        os.makedirs(data_folder_path, exist_ok=True)
        print(f'{datetime.now().strftime("%H:%M:%S")} - Created "{folder}" folder.')

    # Define path to settings file
    settings_file = project_folder.joinpath(project_name + '_settings.xlsx')

    # Create demultiplexing sheet
    cols = ['Forward index', 'Forward primer', 'Forward barcode', 'Reverse index', 'Reverse primer', 'Reverse barcode', 'ID']
    rows = [['', 'AAACTCGTGCCAGCCACC', 'CTGT', '', 'GGGTATCTAATCCCAGTTTG', 'GTCCTA', 'example_1_only_barcodes']]
    demultipexing_df_empty = pd.DataFrame(rows, columns=cols)

    # Create settings sheet
    cols = ['Step', 'Category', 'Variable', 'Comment']
    rows = [['General', 'cpu count', multiprocessing.cpu_count()-1, 'Number of cores to use'],
             ['demultiplexing',
              'allowed errors',
              2,
              'Allowed errors during in adapter sequence'],
            ['quality filtering',
             'truncation value',
             '20',
             'Reads below this length will be discarded'],
             ['quality filtering',
              'minimum length',
              '',
              'Reads below this length will be discarded'],
             ['quality filtering',
              'maximum length',
              '',
              'Reads above this length will be discarded'],
             ['quality filtering',
              'maxee',
              2,
              'Reads above this maxee value will be discarded'],
             ['swarm denoising', 'd', 1, 'Stringency of denoising'],
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

def watch_folder(project_folder, settings_df, demultiplexing_df, live_calling):

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

                    #=======# Demultiplexing #=======#
                    print(f'{datetime.now().strftime("%H:%M:%S")} - Starting cutadapt demultiplexing and primer trimming...')
                    cutadapt_demultiplexing(project_folder, main_file, settings_df, demultiplexing_df)
                    print(f'{datetime.now().strftime("%H:%M:%S")} - Finished cutadapt demultiplexing and primer trimming!')
                    print('')

                    #=======# Quality filtering #=======#
                    print(f'{datetime.now().strftime("%H:%M:%S")} - Starting vsearch quality filtering...')
                    fastq_files = glob.glob(str(project_folder.joinpath('2_demultiplexing', 'data', '*.fastq.gz')))
                    Parallel(n_jobs=cpu_count, backend='loky')(delayed(vsearch_quality_filtering)(project_folder, fastq_file, settings_df) for fastq_file in fastq_files)
                    print(f'{datetime.now().strftime("%H:%M:%S")} - Finished vsearch quality filtering!')
                    print('')

                    #=======# Denoising #=======#
                    print(f'{datetime.now().strftime("%H:%M:%S")} - Starting swarm denoising...')
                    fasta_files = glob.glob(str(project_folder.joinpath('3_quality_filtering', 'data', '*.fasta')))
                    Parallel(n_jobs=cpu_count, backend='loky')(delayed(swarm_denoising)(project_folder, fasta_file, settings_df) for fasta_file in fasta_files)
                    print(f'{datetime.now().strftime("%H:%M:%S")} - Finished swarm denoising...')
                    print('')

                    #=======# Read table #=======#
                    print(f'{datetime.now().strftime("%H:%M:%S")} - Starting to build read table...')
                    create_read_table(project_folder)
                    print(f'{datetime.now().strftime("%H:%M:%S")} - Finished building read table!')
                    print('')

                    #=======# Taxonomic assignment #=======#
                    print(f'{datetime.now().strftime("%H:%M:%S")} - Starting taxonomic assignment...')
                    apscale_taxonomic_assignment(project_folder, settings_df)
                    print(f'{datetime.now().strftime("%H:%M:%S")} - Finished building read table!')
                    print('')

                    # Move file to finish analysis for the file
                    new_file = Path(str(raw_tmp_folder.joinpath(name)) + '.fastq.gz')
                    shutil.move(main_file, new_file)
                    print(f'{datetime.now().strftime("%H:%M:%S")} - Moved {name}...')

                    # Finish file
                    print(f'{datetime.now().strftime("%H:%M:%S")} - Finished analysis for: {name}\n')
                    time.sleep(1)
                    i += 1

                # =======# Create report #=======#
                create_report(project_folder)

                print(f'{datetime.now().strftime("%H:%M:%S")} - Finished analysis for: {name} ({i}/{len(main_files)})')
                print('')

            if live_calling == False:
                break

    except KeyboardInterrupt:
        print('Stopping apscale nanopore live processing.')

def cutadapt_demultiplexing(project_folder, main_file, settings_df, demultiplexing_df):

    # Preprate output files
    # main_file = Path("/Users/tillmacher/Desktop/APSCALE_projects/test_dataset_apscale_nanopore/1_raw_data/backup/1_subset.fastq.gz")
    input_file = main_file
    name = input_file.name.replace('.fastq.gz', '')
    output_folder_tmp = project_folder.joinpath('2_demultiplexing', 'tmp')
    output_folder_data = project_folder.joinpath('2_demultiplexing', 'data')
    output_file = output_folder_tmp.joinpath("{name}.fastq")

    # Create tmp folder
    tmp_folder = project_folder.joinpath('2_demultiplexing', 'tmp')
    os.makedirs(tmp_folder, exist_ok=True)

    # Collect required settings
    number_of_errors = settings_df[settings_df['Category'] == 'allowed errors']['Variable'].values.tolist()[0]
    cpu_count = settings_df[settings_df['Category'] == 'cpu count']['Variable'].values.tolist()[0]

    # Run cutadapt demultiplexing
    g_args = []
    for _, row in demultiplexing_df.iterrows():
        # Create forward sequence
        fwd_seq = row['Forward index']
        # Create reverse sequence
        rvs_seq = row['Reverse index']
        # Combine to search sequence
        search_seq = f'{fwd_seq}...{rvs_seq}'
        g_args.extend(['-g', search_seq])

    # Run cutadapt demultiplexing and primer trimming
    command = f"cutadapt -e {number_of_errors} {' '.join(g_args)} --cores {cpu_count} -o {output_file} --discard-untrimmed --report=minimal {input_file}"
    process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    stdout, stderr = process.communicate()

    # You can now use `stdout` and `stderr` as variables
    in_reads = int(stdout.split()[11])
    out_reads = int(stdout.split()[-3])
    out_reads_perc = round(out_reads / in_reads * 100, 2)
    print(f'{datetime.now().strftime("%H:%M:%S")} - Finished demultiplexing of {name}: {out_reads_perc}% of reads kept.')

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

        # Remove tmp file
        os.remove(tmp_file)

def vsearch_quality_filtering(project_folder, file, settings_df):

    # Preprate output files
    # file = '/Users/tillmacher/Desktop/APSCALE_projects/test_dataset_apscale_nanopore/2_demultiplexing/data/Sample_1.fastq.gz'
    input_file = Path(file)
    name = input_file.name.replace('.fastq.gz', '')
    output_folder_data = project_folder.joinpath('3_quality_filtering', 'data')
    # output files
    filtered_fastq_len = output_folder_data.joinpath(f'{name}_filtered_len.fastq')
    filtered_fastq_trunc = output_folder_data.joinpath(f'{name}_filtered_trunc.fastq')
    filtered_fasta_maxee = output_folder_data.joinpath(f'{name}_filtered_maxee.fasta')
    dereplicated_fasta = output_folder_data.joinpath(f'{name}_filtered_derep.fasta')

    # Collect required settings
    min_len = settings_df[settings_df['Category'] == 'minimum length']['Variable'].values.tolist()[0]
    max_len = settings_df[settings_df['Category'] == 'maximum length']['Variable'].values.tolist()[0]
    trunc_val = settings_df[settings_df['Category'] == 'truncation value']['Variable'].values.tolist()[0]
    maxee = settings_df[settings_df['Category'] == 'maxee']['Variable'].values.tolist()[0]

    # Run vsearch quality filtering
    # Truncate
    command = f"vsearch --fastq_qmax 90 --fastq_filter {input_file} --fastq_truncqual {trunc_val}  --fastqout {filtered_fastq_trunc}"
    process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    stdout, stderr = process.communicate()
    reads_1 = stderr.split()[11]
    reads_0 = int(stderr.split()[-3]) + int(reads_1)

    # LENGTH
    command = f"vsearch --fastq_qmax 90 --fastq_filter {filtered_fastq_trunc} --fastq_minlen {min_len} --fastq_maxlen {max_len} --fastqout {filtered_fastq_len}"
    process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    stdout, stderr = process.communicate()
    reads_2 = stderr.split()[11]

    # Maxee_rate
    command = f"vsearch --fastq_qmax 90 --fastq_filter {filtered_fastq_len} --fastq_maxns 0 --fastq_maxee {maxee}  --fastaout {filtered_fasta_maxee}"
    process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    stdout, stderr = process.communicate()
    reads_3 = stderr.split()[11]

    # Run vsearch dereplication
    command = f"vsearch --derep_fulllength {filtered_fasta_maxee} --sizeout --relabel_sha1 --output {dereplicated_fasta}"
    process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    stdout, stderr = process.communicate()
    reads_4 = stderr.split()[-15]

    print(f'{datetime.now().strftime("%H:%M:%S")} - {name}: {reads_0} (in) -> {reads_1} (trunc) -> {reads_2} (len) -> {reads_3} (maxee) -> {reads_4} (derep)')

    if filtered_fastq_trunc.exists():
        os.remove(filtered_fastq_trunc)
    if filtered_fastq_len.exists():
        os.remove(filtered_fastq_len)
    if filtered_fasta_maxee.exists():
        os.remove(filtered_fasta_maxee)

def swarm_denoising(project_folder, file, settings_df):

    # Preprate output files
    input_file = Path(file)
    name = input_file.name.replace('_filtered_derep.fasta', '')
    output_folder_data = project_folder.joinpath('4_denoising', 'data')
    cluster_file = output_folder_data.joinpath(f'{name}_clusters.fasta')

    # Collect required settings
    d_value = settings_df[settings_df['Category'] == 'd']['Variable'].values.tolist()[0]

    # Run swarm denoising
    command = f"swarm -d {d_value} --threads 1 -z --seeds {cluster_file} {input_file}"
    process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    stdout, stderr = process.communicate()

    # You can now use `stdout` and `stderr` as variables
    swarms = int(stderr.split()[-7])
    print(f'{datetime.now().strftime("%H:%M:%S")} - {name}: {swarms} swarms.')

def create_read_table(project_folder):
    # Prepare output files
    swarm_files_path = project_folder.joinpath('4_denoising', 'data', '*_clusters.fasta')
    swarm_files = [Path(i) for i in glob.glob(str(swarm_files_path))]
    data = defaultdict(lambda: defaultdict(int))  # nested dict: hash -> sample -> size
    seq_dict = {}  # hash -> sequence

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

    # Create DataFrame
    df = pd.DataFrame.from_dict(data, orient='index')  # rows = hashes, columns = samples
    df.fillna(0, inplace=True)  # fill missing sample entries with 0
    df = df.astype(int)  # convert all counts to int
    df.insert(0, "Seq", df.index.map(seq_dict))  # add sequence column
    df.index.name = "ID"  # set index name
    df.reset_index(inplace=True)  # move hash to a column
    df['sum'] = df[df.columns.tolist()[2:]].sum(axis=1)
    df = df.sort_values('sum', ascending=False)
    # remove singletons
    total_reads_before = df['sum'].sum()
    df = df[df['sum'] > 1]
    total_reads_after = df['sum'].sum()
    singletons = total_reads_before - total_reads_after
    print(f'{datetime.now().strftime("%H:%M:%S")} - {total_reads_before} ESVs of which {singletons} were singletons.')
    df = df.drop(columns=['sum'])

    # insert empty files
    for file in sorted(swarm_files):
        sample = file.name.replace('_clusters.fasta', '')
        if sample not in list(df.columns):
            df[sample] = [0] * len(df)



    # Write to files
    if df.shape[0] < 65000:
        excel_file = project_folder.joinpath('5_esv_table', 'Swarms.xlsx')
        df.to_excel(excel_file, index=False)
    parquet_file = project_folder.joinpath('5_esv_table', 'Swarms.parquet.snappy')
    df.to_parquet(parquet_file, compression='snappy')

    # Write sequences to fasta
    fasta_file = project_folder.joinpath('5_esv_table', 'data', 'swarms.fasta')
    with open(fasta_file, "w") as output_handle:
        for hash, seq in seq_dict.items():
            record = SeqRecord(Seq(seq), id=hash, description='')
            SeqIO.write(record, output_handle, "fasta")

def apscale_taxonomic_assignment(project_folder, settings_df):
    # Define files
    fasta_file = project_folder.joinpath('5_esv_table', 'data', 'swarms.fasta')
    project_name = project_folder.name.replace('_apscale_nanopore', '')
    results_folder = project_folder.joinpath('6_taxonomic_assignment', project_name)

    # Collect variables
    run_blastn = settings_df[settings_df['Category'] == 'apscale blast']['Variable'].values.tolist()[0]
    blastn_db = settings_df[settings_df['Category'] == 'apscale db']['Variable'].values.tolist()[0]

    # Run apscale blast
    if run_blastn == 'yes':
        if results_folder.exists():
            shutil.rmtree(results_folder)
            os.makedirs(results_folder, exist_ok=True)
        command = f"apscale_blast -db {blastn_db} -q {fasta_file} -o {results_folder} -task megablast"
        process = subprocess.Popen(command, shell=True, text=True)
        process.wait()

def create_report(project_folder):
    print(f"Creating report")
    time.sleep(0.25)
    # Do your processing here...

def main():
    """
    APSCALE nanopore suite
    Command-line tool to process nanopore sequence data.
    """

    # Introductory message with usage examples
    message = """
    APSCALE nanopore command line tool - v0.0.1
    Example commands:
    $ apscale_nanopore --create_project
    $ apscale_blast --run_apscale

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
    run_parser.add_argument('-e', type=str, help='Overwrite: allowed errors.')
    run_parser.add_argument('-t', type=str, help='Overwrite: quality truncation cutoff.')
    run_parser.add_argument('-minlen', type=str, help='Overwrite: minimum length.')
    run_parser.add_argument('-maxlen', type=str, help='Overwrite: maximum length.')
    run_parser.add_argument('-maxee', type=str, help='Overwrite: maxee value.')
    run_parser.add_argument('-d', type=str, help="Overwrite: swarm's d value.")

    # Parse arguments
    args = parser.parse_args()

    # Create project
    if args.command == 'create':
        project_folder = Path(str(Path(args.project)) + '_apscale_nanopore')
        project_name = project_folder.name.replace('_apscale_nanopore', '')
        create_project(project_folder, project_name)

    # Run apscale
    elif args.command == 'run':
        # Define folders and files
        project_folder = Path(args.project)
        project_name = project_folder.name.replace('_apscale_nanopore', '')
        settings_file = project_folder.joinpath(project_name + '_settings.xlsx')

        # Load settings dataframe
        if settings_file.exists():
            settings_df = pd.read_excel(settings_file, sheet_name='Settings').fillna('')
            demultiplexing_df = pd.read_excel(settings_file, sheet_name='Demultiplexing').fillna('')

            # Check if argument require to be adjusted
            if args.e:
                index = settings_df[settings_df['Category'] == 'allowed errors'].index[0]
                settings_df.loc[index, 'Variable'] = args.e
                print(f'Adjusted value: Number of allowed errors: {args.e}')
            if args.t:
                index = settings_df[settings_df['Category'] == 'truncation value'].index[0]
                settings_df.loc[index, 'Variable'] = args.e
                print(f'Adjusted value: Truncation cutoff: {args.t}')
            if args.minlen:
                index = settings_df[settings_df['Category'] == 'minimum length'].index[0]
                settings_df.loc[index, 'Variable'] = args.minlen
                print(f'Adjusted value: Minimum length: {args.minlen}')
            if args.maxlen:
                index = settings_df[settings_df['Category'] == 'maximum length'].index[0]
                settings_df.loc[index, 'Variable'] = args.maxlen
                print(f'Adjusted value: Maximum length: {args.maxlen}')
            if args.maxee:
                index = settings_df[settings_df['Category'] == 'maxee'].index[0]
                settings_df.loc[index, 'Variable'] = args.maxee
                print(f'Adjusted value: Maximum expected error: {args.maxee}')
            if args.d:
                index = settings_df[settings_df['Category'] == 'd'].index[0]
                settings_df.loc[index, 'Variable'] = args.d
                print(f"Adjusted value: Swarm's d value: {args.d}")

            # =======# Live processing #=======#
            print('')
            watch_folder(project_folder, settings_df, demultiplexing_df, args.live_calling)

        else:
            print(settings_file)
            print('Error: Cannot find settings file!')

if __name__ == "__main__":
    main()