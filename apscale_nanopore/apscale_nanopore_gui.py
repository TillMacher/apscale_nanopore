import os
import shutil
import time
import streamlit as st
from streamlit_file_browser import st_file_browser
import importlib
from update_checker import update_check
from pathlib import Path
import glob
import pandas as pd
import subprocess
import multiprocessing
import platform
from ete3 import NCBITaxa
from distutils.util import strtobool
import webbrowser
import importlib.util
import glob, datetime, os
import pandas as pd
import numpy as np
from playwright.sync_api import sync_playwright
import zipfile
import importlib.metadata
import subprocess
import asyncio
import sys

# Apscale nanopore
from apscale_nanopore.__main__ import apscale_nanopore
from apscale_nanopore.__main__ import import_raw_data
from apscale_nanopore.__main__ import quality_report

# Apscale blast
from apscale_blast.a_blastn import main as a_blastn
from apscale_blast.b_filter import main as b_filter
from apscale_blast.__main__ import organism_filter as organism_filter


# help texts
n_cores_help = """
    All settings required by Apscale are configured within the project’s settings file. Each processing step has its
    own sheet in the document. The first tab you’ll see when opening the settings file is "0_general_settings". In this
    tab, you can define how many CPU cores Apscale should use for processing. By default, it uses the total number of
    available cores minus two.
  """

compression_level_help = """
    You can also set the gzip compression level here. Apscale compresses all output files to conserve disk space. The
    default compression level is 6, which is suitable for most use cases. If you need to save more space, you can
    increase it to 9—this will produce smaller files but may slow down processing.
  """

b_pe_merging_help = """ 
The first step performed by Apscale is merging paired-end reads using vsearch. The default settings are fairly relaxed
to merge the largest possible portion of reads, as quality filtering is handled in later steps.
                    """

c_primer_trimming_help = """ 
The next step performed by Apscale is primer trimming, which removes the primers used for target amplification since
they do not contain biologically relevant information.
                    """

d_quality_filtering_help = """ 
After primer trimming, Apscale performs quality filtering. This step filters out reads with an expected error higher
than the threshold defined in the settings, as well as sequences whose lengths fall outside the specified target range.
Typically, we use a tolerance of ±10 bases around the target length to allow for some biological variation while
removing artifacts such as primer dimers.
                    """

e_dereplication_help = """ 
Before running the modular workflow, the reads from all samples must be dereplicated. This step does not alter the data
itself but optimizes how it is stored. The output is a FASTA file containing all unique sequences with size annotations (e.g., >seq_1;size=100).
                    """

f_denoising_help = """ 
The denoising module performs sequence denoising using vsearch, processing each sample file individually. Pooling is
intentionally avoided to ensure that the resulting sequences remain independent of the overall dataset size. This
design choice guarantees that previously processed data remains unaffected when new samples are added and the project
is reanalyzed. During denoising, Apscale automatically assigns unique identifiers to each sequence using the SHA3-256
hashing algorithm.
                    """

f_denoising_threshold_help = """ 
Several threshold types are available to control which reads are considered for denoising.  By default, an absolute
threshold is applied (minsize = 4), meaning that only sequences with an abundance of four or more are retained —
effectively removing a substantial amount of low-abundance noise.  Alternatively, a relative threshold can be used to
retain only those sequences that represent a defined percentage of the sample’s total read count (e.g., 0.01%). Since
both absolute and relative thresholds are inherently arbitrary, we introduced a third option in version 4.0.0: power
law–based filtering. Read abundance distributions typically follow a power law, where a few sequences are highly
abundant (true biological signals) and many are rare (a mixture of real low-abundance taxa, sequencing noise, and
PCR artifacts). This filtering method fits a power law model to each sample’s read distribution and sets the threshold
at the point where the observed distribution deviates from the expected power law curve. The underlying assumption is
that this inflection marks a shift in the signal-to-noise ratio, with noise becoming dominant. This approach results
in a data-driven, rather than arbitrary, threshold for denoising.
                    """

g_swarm_clustering_help = """ 
Alternatively, files can be clustered individually using the Swarm algorithm with d=1 and the fastidious option enabled
by default (see https://github.com/torognes/swarm for details). We consider Swarm clustering as an alternative to denoising; when retaining
only the chosen centroid sequences, the results are generally quite similar. Additionally, the output from the
denoising module can be further clustered with Swarm if desired.
                    """

j_generate_read_table_help = """ 
This module generates the read table and performs threshold-based sequence grouping, similar to classical OTU
clustering. Apscale always outputs both sequences (ESVs) and sequence groups (OTUs). The read table is saved in
Parquet format and, if the dataset contains fewer than 1,000,000 distinct sequences, also in Excel format. Additionally,
this module creates a “read data store,” a DuckDB (https://duckdb.org/) database that contains comprehensive
information about sequences, groups, samples, and read counts. The read data store efficiently handles even very
large datasets—potentially billions of sequences—at high speed, without requiring the entire dataset to be loaded into
memory. This makes it especially useful for scaling up analyses.
                    """

def check_dependencies(tools=["cutadapt", "vsearch", "swarm", "blastn"]):
    missing = []
    for tool in tools:
        if shutil.which(tool) is None:
            missing.append(tool)
    if missing:
        missing_tools = ', '.join(missing)
        if len(missing) == 1:
            st.error(f'WARNING: The following tool is missing: **{missing_tools}**')
        else:
            st.error(f'WARNING: The following tools are missing: **{missing_tools}**')
        st.warning('⚠️ Please install all required tools, either manually or using the "apscale_installer".')

def check_package_update(package_name):
    # Get the currently installed version
    try:
        installed_version = importlib.metadata.version(package_name)
    except importlib.metadata.PackageNotFoundError:
        print(f"{package_name} is not installed.")
        return

    # Check for updates
    res = update_check(package_name, installed_version)
    return res

def get_package_versions():
    packages = ["apscale-nanopore","apscale_blast", "cutadapt"]

    for pkg in packages:
        try:
            version = importlib.metadata.version(pkg)
            st.write(f"{pkg}: {version}")
        except importlib.metadata.PackageNotFoundError:
            st.write(f"{pkg}: not installed")

    for tool in ['vsearch', 'swarm']:
        try:
            result = subprocess.run([tool, "--version"], capture_output=True, text=True, check=True)
            version = result.stderr.strip().split('\n')[0]
            st.write(f"{tool}: {version}")
        except:
            st.write(f"{tool}: not installed")

    try:
        result = subprocess.run(["blastn", "-version"], capture_output=True, text=True, check=True)
        version = result.stdout.strip().split('\n')[0]
        st.write(version)
    except:
        st.write(f"blastn: not installed")

def open_folder(folder_path):
    # Get the current operating system
    current_os = platform.system()

    # Open the folder based on the OS
    try:
        if current_os == "Windows":
            subprocess.Popen(f'explorer "{folder_path}"')
        elif current_os == "Darwin":  # macOS
            subprocess.Popen(['open', folder_path])
        else:  # Linux
            subprocess.Popen(['xdg-open', folder_path])
    except Exception as e:
        print(f"Failed to open folder: {e}")

def open_file(path: Path):
    path = str(path)
    if platform.system() == "Darwin":
        subprocess.run(["open", path])
    elif platform.system() == "Windows":
        os.startfile(path)
    else:
        subprocess.run(["xdg-open", path])

def create_project(project_folder, project_name, is_demultiplexed):
    """Create project folder structure and copy default settings file."""

    # Define subfolders
    sub_folders = [
        '1_raw_data', '2_index_demultiplexing', '3_primer_trimming',
        '4_quality_filtering', '5_clustering_denoising', '6_read_table',
        '7_taxonomic_assignment', '8_nanopore_report'
    ]

    # Create all subfolders with a nested 'data' folder
    for folder in sub_folders:
        (project_folder / folder / 'data').mkdir(parents=True, exist_ok=True)

    # Decide which default settings to use
    if is_demultiplexed == True:
        settings_file = 'default_settings_nd.xlsx'
    else:
        settings_file = 'default_settings.xlsx'
    src_settings = Path(__file__).resolve().parent / settings_file
    dst_settings = project_folder / f"{project_name}_settings.xlsx"

    shutil.copyfile(src_settings, dst_settings)

    # Log messages
    st.info(f'Created settings file: {dst_settings.name}')
    st.info(f'Copy your data into the "1_raw_data/data" folder.')

    # Open the copied settings file
    open_file(dst_settings)

def read_settings_file(settings_xlsx, settings_df, demultiplexing_df):
    for category, variable in settings_df[['Category', 'Variable']].values.tolist():
        # Normalize variable
        if isinstance(variable, str):
            var_lower = variable.strip().lower()
            if var_lower in ["true", "false"]:  # handle booleans
                variable = var_lower == "true"
        elif pd.isna(variable):  # handle NaN
            variable = ""
        else:
            variable = str(variable)
        st.session_state[category] = variable

def update_settings_file(settings_xlsx, settings_df, demultiplexing_df):
    print('')

    # Update Settings sheet with current session state values
    for category in settings_df['Category'].values.tolist():
        index = settings_df[settings_df['Category'] == category].index[0]
        value = st.session_state[category]

        # Try to convert to int, then float, otherwise keep as is
        if isinstance(value, str):
            try:
                if '.' in value:
                    num = float(value)
                    if num.is_integer():
                        value = int(num)
                    else:
                        value = num
                else:
                    value = int(value)
            except ValueError:
                pass  # keep original string if not numeric

        settings_df.loc[index, 'Variable'] = value
        print(f'Adjusted value: {category} -> {value} (type: {type(value).__name__})')

    # Save both sheets back into the Excel file
    with pd.ExcelWriter(settings_xlsx, engine='openpyxl') as writer:
        # Ensure Demultiplexing sheet is written first
        demultiplexing_df.to_excel(writer, sheet_name='Demultiplexing', index=False)
        settings_df.to_excel(writer, sheet_name='Settings', index=False)

    st.info(f"Settings were updated and saved to {settings_xlsx.name}")

def main():

    st.set_page_config(layout="wide")

    check_dependencies()

    # Sidebar inputs & outputs
    with st.sidebar:
        st.subheader("APSCALE projects")

        # read user_data.txt
        script_path = Path(__file__).resolve()
        user_data_txt = script_path.parent / '_user_data' / 'user_data.txt'

        default_value = ""
        if user_data_txt.exists():
            with open(user_data_txt, 'r', encoding='utf-8') as f:
                default_value = f.read().strip()  # strip removes newlines/spaces

        # Get user input
        path_to_projects = Path(st.text_input('Enter Path to APSCALE Projects', value=default_value))
        database_folder = path_to_projects / 'APSCALE_databases'

        if path_to_projects == Path('.'):
            st.write('Please select your APSCALE projects folder.')
        else:
            try:
                # Collect only folders that contain "_apscale" in their name
                if path_to_projects.exists() and path_to_projects.is_dir():
                    apscale_folders = [
                        f for f in path_to_projects.iterdir()
                        if f.is_dir() and "_apscale_nanopore" in f.name
                    ]

                    if st.button(label='Remember project folder', key='remember_project_folder',use_container_width=True):
                        script_path = Path(__file__).resolve()
                        user_data_txt = script_path.parent / '_user_data' / 'user_data.txt'
                        user_data_txt.parent.mkdir(parents=True, exist_ok=True)
                        with open(user_data_txt, 'w', encoding='utf-8') as f:
                            f.write(str(path_to_projects))
                        st.success(f"Saved project folder!")

                    if apscale_folders:
                        project_folder = st.selectbox(
                            "Select an APSCALE project folder:",
                            sorted(apscale_folders, key=lambda p: p.name.lower()),  # sort by name only
                            format_func=lambda p: p.name
                        )
                        project_name = project_folder.name.replace('_apscale_nanopore', '')
                        if st.button('Open project folder', use_container_width=True):
                            open_folder(project_folder)

                    else:
                        st.warning("No APSCALE project folders found here.")

                else:
                    st.error("The given path does not exist or is not a directory.")

                if not database_folder.exists():
                    st.error("Please first initialise the database folder.")
                    if st.button('Create database folder', use_container_width=True):
                        os.makedirs(database_folder, exist_ok=True)
                        st.rerun()

                st.markdown("---")
                st.subheader('Create new project')
                st.text_input(label='Enter name of new project', key='new_project_name')
                new_project_path = Path(str(path_to_projects / Path(st.session_state['new_project_name'])) + '_apscale_nanopore')
                new_project_name = new_project_path.name.replace('_apscale_nanopore', '')
                if new_project_path != path_to_projects:
                    st.selectbox("Is your data demultiplexed?", key='is_demultiplexed', options=[True, False])
                    if st.button('Create new project', use_container_width=True):
                        create_project(new_project_path, new_project_name, st.session_state['is_demultiplexed'])
                        st.success('Create new project folder!')

                st.markdown("---")
                st.subheader('APSCALE databases')
                n_databases = len(glob.glob(str(database_folder / '*')))
                st.info(f'{n_databases} databases available')

                if st.button('Open Database Hub', use_container_width=True):
                    webbrowser.open('https://seafile.rlp.net/d/474b9682a5cb4193a6ad/')
                if st.button('Download All Latest Databases', use_container_width=True):
                    public_url = 'https://seafile.rlp.net/d/474b9682a5cb4193a6ad/?p=%2FLatest&mode=list'
                    download_seafile_zip(public_url, database_folder)
                if st.button('Open APSCALE database folder', use_container_width=True):
                    open_folder(database_folder)

            except Exception as e:
                st.error(f"An error occurred: {e}")

        st.markdown("---")
        st.subheader("Refresh")

        if st.button("🔄 Refresh files and folders", use_container_width=True):
            st.session_state.pop("_last_settings_xlsx", None)
            st.success("Files and folders refreshed.")


    ####################################################################################################################
    if not path_to_projects.exists() and not path_to_projects.is_dir():
        st.info('Please select a project to continue!')
    elif 'project_folder' not in locals():
        st.info('Project folder is not selected yet!')
    else:
        ################################################################################################################
        # Raw data check
        st.header("APSCALE-Nanopore")

        raw_data = glob.glob(str(project_folder / '1_raw_data' / 'data' / '*.fastq*'))
        n_raw_data = len(raw_data)

        # Read Settings file
        settings_xlsx = project_folder / Path(str(project_name) + '_settings.xlsx')
        settings_df = pd.read_excel(settings_xlsx, sheet_name='Settings').fillna('')
        demultiplexing_df = pd.read_excel(settings_xlsx, sheet_name='Demultiplexing').fillna('')

        if st.session_state.get("_last_settings_xlsx") != settings_xlsx:
            st.session_state["_last_settings_xlsx"] = settings_xlsx
            read_settings_file(settings_xlsx, settings_df, demultiplexing_df)

        if n_raw_data == 0:
            st.info('No raw data files to demultiplex!')
        else:
            st.subheader('Nanopore Raw Data Processing')

            with st.expander("⚙️ General settings", expanded=True):
                col1, col2 = st.columns(2)
                with col1:
                    cpu_values = ['auto'] + [i for i in range(1, multiprocessing.cpu_count()+1)][::-1]
                    st.selectbox(label='CPU count', key='cpu count', options=cpu_values)
                with col2:
                    st.selectbox(label='APSCALE-Naonopore mode', key='nanopore mode', options=['linear', 'live base-calling'])

            with st.expander("🧬 Index Demultiplexing", expanded=True):
                st.selectbox(label='Skip demultiplexing', key='sd', options=[False, True])
                col1, col2 = st.columns(2)
                with col1:
                    st.text_input(label='Number of allowed errors', key='allowed errors index')
                with col2:
                    st.info('Please adjust the Demultiplexing sheet according to your dataset! (Only the first 10 rows are shown)')
                if st.session_state['sd'] == False:
                    st.dataframe(demultiplexing_df[['ID', 'Forward index 5-3', 'Reverse index 5-3']].head(10))
                    if st.button('Open demultiplexing sheet'):
                        open_file(settings_xlsx)

            with st.expander("✂️ Primer Trimming", expanded=True):
                col1, col2 = st.columns(2)
                with col1:
                    st.text_input(label='Number of allowed errors', key='allowed errors primer')
                with col2:
                    st.info('Please adjust the Demultiplexing sheet according to your dataset! (Only the first 10 rows are shown)')
                if st.session_state['sd'] == True:
                    st.dataframe(demultiplexing_df[['ID', 'Forward primer 5-3', 'Reverse primer 5-3', 'Minimum length','Maximum length']].head(10))
                    n_raw_files = len(glob.glob(str(project_folder / '1_raw_data' / 'data' / '*.fastq*')))
                    col3, col4 = st.columns(2)
                    with col3:
                        st.info(f'Found {n_raw_files} raw files.')
                    with col4:
                        if st.button('Import raw data and update "Demultiplexing" sheet'):
                            import_raw_data(project_folder, settings_xlsx, demultiplexing_df)
                else:
                    st.dataframe(demultiplexing_df[['ID', 'Forward primer 5-3', 'Reverse primer 5-3', 'Minimum length','Maximum length']].head(10))
                if st.button('Open demultiplexing sheet '):
                    open_file(settings_xlsx)

            with st.expander("📉 Quality Filtering", expanded=True):
                st.text_input(label='Minimum average PHRED quality', key='minimum quality')

            with st.expander("🌀 Clustering and Denoising", expanded=True):
                col1, col2 = st.columns(2)
                with col1:
                    st.selectbox(label='Clustering/denoising mode', key='mode', options=['ESVs', 'Swarms', 'denoised OTUs', 'Swarm OTUs'])
                with col2:
                    st.selectbox(label='Representative sequence', key='representative', options=['centroid', 'consensus'])
                col1, col2, col3 = st.columns(3)
                with col1:
                    st.text_input(label='Alpha value', key='alpha')
                with col2:
                    st.text_input(label='Percentage identity', key='percid')
                with col3:
                    st.text_input(label='SWARM d value', key='d')
                if st.session_state['mode'] == 'ESVs':
                    st.info(f'Using alpha value: {st.session_state["alpha"]}')
                if st.session_state['mode'] == 'Swarms':
                    st.info(f'Using SWARM d value: {st.session_state["d"]}')
                if st.session_state['mode'] == 'denoised OTUs':
                    st.info(f'Using alpha value: {st.session_state["alpha"]}')
                    st.info(f'Using percid value: {st.session_state["percid"]}')
                if st.session_state['mode'] == 'Swarm OTUs':
                    st.info(f'Using SWARM d value: {st.session_state["d"]}')
                    st.info(f'Using percid value: {st.session_state["percid"]}')

            with st.expander("📊 Read Table", expanded=True):
                col1, col2 = st.columns(2)
                with col1:
                    st.text_input(label='Read cutoff', key='minimum reads')
                with col2:
                    if float(st.session_state['minimum reads']) < 1:
                        st.info('Relative threshold detected.')
                        cutoff = float(st.session_state['minimum reads'])
                        st.info(f'Example: For 60,000 reads → cutoff = {60000 * cutoff:.0f} reads')
                    else:
                        st.info('Absolute threshold detected.')

            with st.expander("🔬 Taxonomic Assignment", expanded=True):
                st.selectbox(label='Perform taxonomic assignment', key='apscale blast', options=['Yes', 'No'])
                st.selectbox(label='Task', key='task', options=['megablast', 'blastn', 'dc-megablast'])
                available_databases = {Path(i).name:Path(i) for i in glob.glob(str(database_folder / '*'))}
                selected_db = st.selectbox(label='Database', key='selected_db', options=list(available_databases.keys()))

            if st.session_state['nanopore mode'] == 'live base-calling':
                st.warning('Live base-calling is currently only possible through the command line.')
            elif selected_db == None and st.session_state['apscale blast'] == 'Yes':
                st.warning('Please select a database!')
            else:
                st.markdown("---")
                st.subheader('Raw data analysis')

                col1, col2 = st.columns(2)

                with col1:
                    st.text('Select tasks to perform:')
                    all_steps = {"1": "Index demultiplexing", "2": "Primer trimming",
                                 "3": "Quality filtering", "4": "Clustering/denoising", "5": "Read table",
                                 "6": "Tax. assignment"}
                    selected_steps = {}
                    for key, label in all_steps.items():
                        selected_steps[key] = st.toggle(label, key=f"step_{key}", value=True)

                    # Example: get list of selected steps
                    active_steps = [all_steps[k] for k, v in selected_steps.items() if v]

                with col2:
                    if st.session_state['sd'] == True:
                        if st.session_state['sd'] is True:
                            st.warning(
                                '⚠️ Demultiplexing is disabled. Raw FASTQ files will be copied directly '
                                'to the folder **2_index_demultiplexing** without processing.'
                            )

                    st.info(
                        "✅ All required information is available. "
                        "Before starting raw data processing, please review the **Demultiplexing** sheet."
                    )
                    if st.button('Run Apscale-Nanopore'):
                        if st.session_state['apscale blast'] == 'Yes':
                            st.session_state['apscale db'] = available_databases[selected_db]
                        update_settings_file(settings_xlsx, settings_df, demultiplexing_df)
                        settings_df = pd.read_excel(settings_xlsx, sheet_name='Settings')
                        demultiplexing_df = pd.read_excel(settings_xlsx, sheet_name='Demultiplexing')
                        apscale_nanopore(project_folder, settings_df, demultiplexing_df, active_steps, st.session_state['sd'])

                st.markdown("---")
                st.subheader('Quality report')
                col1, col2 = st.columns(2)
                with col1:
                    st.selectbox(label='Select folder', key='qc_folder', options=["1_raw_data", "2_index_demultiplexing", "3_primer_trimming"])
                    if st.button(f'Generate report for "{st.session_state["qc_folder"]}"'):
                        quality_report(project_folder, st.session_state['qc_folder'])
                        st.success('Finished creating report!')
                with col2:
                    output_folder = project_folder / '8_nanopore_report' / 'data'
                    available_reports = {Path(i).name: Path(i) for i in glob.glob(str(output_folder / '*')) if Path(i).is_dir()}
                    files = None
                    if available_reports:
                        st.selectbox(label='Available reports', key='report_to_display', options=available_reports.keys())
                        files = glob.glob(str(available_reports[st.session_state['report_to_display']] / '*.html'))
                        plots = {
                            'Mean PHRED distribution': 'mean_phred_distribution',
                            'Quality by length-bin': 'quality_by_lengthbin',
                            'Quality vs. length': 'quality_vs_length',
                            'Read length distribution': 'read_length_distribution'
                        }
                        st.selectbox(label='Select plot to display', key='report_display', options=[None] + list(plots.keys()))
                if files and available_reports and st.session_state['report_display'] != None:
                    selected_files = [i for i in files if plots[st.session_state['report_display']] in i]
                    with open(selected_files[0], "r", encoding="utf-8") as f:
                        html_content = f.read()
                    st.components.v1.html(html_content, height=550, scrolling=True)

                st.markdown("---")
                st.subheader('Summary report')

                output_folder = project_folder / '8_nanopore_report' / 'main'
                available_reports = {Path(i).name: Path(i) for i in glob.glob(str(output_folder / '*.html'))}
                col1, col2 = st.columns(2)
                with col1:
                    st.toggle(label='Show summary report', key='display_summary', value=False)
                with col2:
                    st.info('The summary report is automatically created and updated when running APSCALE-Nanopore.')
                if available_reports and st.session_state['display_summary'] == True:
                    for file in available_reports.values():
                        if file.exists():
                            with open(file, "r", encoding="utf-8") as f:
                                html_content = f.read()
                            st.components.v1.html(html_content, height=550, scrolling=True)

        ############################################################################################################
        st.markdown("---")
        st.header("📚 Links and Tutorials")
        # GitHub Projects
        with st.expander("🔗 GitHub Repositories", expanded=False):
            st.markdown("""
            - [APSCALE-Nanopore](https://github.com/TillMacher/apscale_nanopore)  
            - [APSCALE-blast](https://github.com/TillMacher/apscale_blast)
            - [VSEARCH](https://github.com/torognes/vsearch)
            - [CUTADAPT](https://github.com/marcelm/cutadapt)
            - [SWARM](https://github.com/torognes/swarm)
            """)
        # Tutorials
        with st.expander("🎥 Video Tutorials", expanded=False):
            st.markdown("""
            - coming soon
            """)
        # Manual
        with st.expander("📖 Documentation", expanded=False):
            st.markdown("""
            - coming soon
            """)
        # Citations
        with st.expander("📑 Citations", expanded=False):
            st.markdown("""
            
            APSCALE-Nanopore         
            - coming soon
            
            CUTADAPT         
            - Martin, M. (2011). Cutadapt removes adapter sequences from high-throughput sequencing reads. EMBnet. Journal, 17(1), Article 1.
            
            VSEARCH
            - Rognes, T., Flouri, T., Nichols, B., Quince, C., & Mahé, F. (2016). VSEARCH: a versatile open source tool for metagenomics. PeerJ, 4, e2584.

            SWARM
            - Mahé, F., Czech, L., Stamatakis, A., Quince, C., de Vargas, C., Dunthorn, M., & Rognes, T. (2021). Swarm v3: Towards tera-scale amplicon clustering. Bioinformatics, 38(1), 267–269. https://doi.org/10.1093/bioinformatics/btab493

            """)
        with st.expander("📦 Package versions", expanded=False):
            get_package_versions()

if __name__ == "__main__":
    main()


