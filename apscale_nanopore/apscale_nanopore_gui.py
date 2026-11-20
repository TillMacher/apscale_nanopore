import os
import shutil
import time
import streamlit as st
import importlib
from update_checker import update_check
from update_checker import UpdateChecker
from pathlib import Path
import glob
import pandas as pd
import subprocess
import multiprocessing
import platform
from ete3 import NCBITaxa
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

def get_package_versions_disp(pkg):
    try:
        version = importlib.metadata.version(pkg)
        return version
    except importlib.metadata.PackageNotFoundError:
        return

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
    packages = ["apscale_nanopore", "apscale_blast", "cutadapt"]

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

def unzip_and_remove(zip_path, extract_to=None):
    zip_path = Path(zip_path)
    if extract_to is None:
        extract_to = zip_path.with_suffix("")  # folder named after zip
    extract_to = Path(extract_to)
    extract_to.mkdir(parents=True, exist_ok=True)

    with zipfile.ZipFile(zip_path, 'r') as zip_ref:
        zip_ref.extractall(extract_to)

    # remove the zip after extraction
    zip_path.unlink()

    print(f"Extracted to {extract_to} and removed {zip_path.name}")
    return extract_to

def flatten_zip_files(output_dir):
    output_dir = Path(output_dir)
    # Find all .zip files inside subfolders
    zip_files = list(output_dir.rglob("*.zip"))

    for zip_file in zip_files:
        target = output_dir / zip_file.name
        # Move to output_dir, overwrite if exists
        shutil.move(str(zip_file), str(target))
        print(f"Moved {zip_file} -> {target}")

    # Optional: clean up empty subfolders
    for folder in sorted(output_dir.rglob("*"), reverse=True):
        if folder.is_dir() and not any(folder.iterdir()):
            folder.rmdir()
            print(f"Removed empty folder {folder}")

if sys.platform.startswith("win"):
    asyncio.set_event_loop_policy(asyncio.WindowsProactorEventLoopPolicy())

def download_seafile_zip(public_url, output_dir="downloads"):
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    st.info("Downloading all latest database ZIP files. This may take a few minutes…")

    with sync_playwright() as p:
        browser = p.chromium.launch(headless=True)
        page = browser.new_page()
        page.goto(public_url)

        # Trigger the download inside the expect_download block
        with page.expect_download() as download_info:
            page.click("button:has-text('ZIP')")
        download = download_info.value

        # Save the downloaded file
        download_zip = output_dir / download.suggested_filename
        download.save_as(str(download_zip))
        st.success(f"Saved {download.suggested_filename}")
        browser.close()

    st.info('Extracting databases.')
    if download_zip.exists():
        unzip_and_remove(download_zip)
        flatten_zip_files(output_dir)

        all_zip_files = list(output_dir.glob("*.zip"))
        for zip_file in all_zip_files:
            unzip_and_remove(zip_file)

    st.success(f"Finished the extraction of {len(list(output_dir.glob('*')))} databases!")

########################################################################################################################
# Page Setup
st.set_page_config(
    page_title="APSCALE-Nanopore",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded"
)

@st.cache_data(ttl=86400, show_spinner=False)
def get_update_messages(packages):
    # The update check needs internet: run it only once per day (cached) instead of on every click
    messages = []
    for pkg in packages:
        try:
            installed_version = importlib.metadata.version(pkg)
            checker = UpdateChecker()
            result = checker.check(package_name=pkg, package_version=installed_version)
            if result:
                print(f'\n{result}\n$ pip install --upgrade {pkg}\n')
                messages.append(f'{result}\n\n$ pip install --upgrade {pkg}')
        except Exception:
            pass # offline or package not installed via pip
    return messages

def project_selected():
    # Pages that need a project show this hint
    if project_folder is not None:
        return True
    st.warning('Please select an APSCALE-Nanopore project in the sidebar first.')
    return False

########################################################################################################################
# Sidebar: Working directory & Project selection

# These variables are used by all pages
path_to_projects = Path('.')
database_folder = Path('.')
project_folder = None
project_name = ''

with st.sidebar:
    apscale_nanopore_version = get_package_versions_disp("apscale_nanopore")
    st.markdown(f"### APSCALE-Nanopore v{apscale_nanopore_version}" if apscale_nanopore_version else "### APSCALE-Nanopore")
    for message in get_update_messages(["apscale_nanopore", "apscale_blast"]):
        st.info(message)

    #####################################
    # Working directory (read user_data.txt)
    script_path = Path(__file__).resolve()
    user_data_txt = script_path.parent / '_user_data' / 'user_data.txt'

    default_value = ""
    if user_data_txt.exists():
        with open(user_data_txt, 'r', encoding='utf-8') as f:
            default_value = f.read().strip()  # strip removes newlines/spaces

    with st.expander('📂 Working directory', expanded=default_value == ''):
        path_to_projects = Path(st.text_input('Path to your APSCALE projects folder:', value=default_value, key='path_to_projects'))
        if st.button(label='💾 Remember project folder', key='remember_project_folder', width='stretch'):
            user_data_txt.parent.mkdir(parents=True, exist_ok=True)
            with open(user_data_txt, 'w', encoding='utf-8') as f:
                f.write(str(path_to_projects))
            st.success(f"Saved project folder!")

    database_folder = path_to_projects / 'APSCALE_databases'

    if path_to_projects == Path('.'):
        st.info('Please enter the path to your APSCALE projects folder.')

    elif not path_to_projects.is_dir():
        st.error("The given path does not exist or is not a directory.")

    else:
        if not database_folder.exists():
            st.warning("Please first initialise the database folder.")
            if st.button('Create database folder', width='stretch'):
                os.makedirs(database_folder, exist_ok=True)
                st.rerun()

        #####################################
        st.write(""" #### 🗂️ Project """)

        # Collect only folders that contain "_apscale_nanopore" in their name
        apscale_folders = [f for f in path_to_projects.iterdir() if f.is_dir() and "_apscale_nanopore" in f.name]

        if apscale_folders:
            project_folder = st.selectbox(
                "Active project:",
                sorted(apscale_folders, key=lambda p: p.name.lower()),  # sort by name only
                format_func=lambda p: p.name,
                key='active_project'
            )
            project_name = project_folder.name.replace('_apscale_nanopore', '')
        else:
            st.warning("No APSCALE-Nanopore project folders found here.")

        # Create a new project
        with st.expander('➕ Create new project'):
            st.text_input(label='Name of the new project:', key='new_project_name')
            st.selectbox("Is your data demultiplexed?", key='is_demultiplexed', options=[True, False])
            if st.button('⚙️ Create project', width='stretch'):
                if st.session_state['new_project_name'] != '':
                    new_project_name = st.session_state['new_project_name']
                    new_project_path = path_to_projects / f'{new_project_name}_apscale_nanopore'
                    if new_project_path.exists():
                        st.error('This project already exists!')
                    else:
                        create_project(new_project_path, new_project_name, st.session_state['is_demultiplexed'])
                        st.success('Created new project folder!')
                        st.rerun()
                else:
                    st.error('Please enter a valid project name!')

        if project_folder is not None:
            if st.button('📂 Open project folder', width='stretch'):
                open_folder(project_folder)

    #####################################
    st.write(""" #### 🔄 Refresh """)
    if st.button("🔄 Refresh files and folders", width='stretch'):
        st.session_state.pop("_last_settings_xlsx", None)
        st.success("Files and folders refreshed.")

########################################################################################################################
# PAGES
########################################################################################################################

########################################################################################################################
# Home
def page_home():
    st.title('🧬 APSCALE-Nanopore')
    st.markdown('Process your Nanopore DNA metabarcoding raw data from demultiplexing to taxonomic assignment. Use the **pages in the sidebar** to navigate between the steps.')

    check_dependencies()

    if project_folder is not None:
        n_raw = len(glob.glob(str(project_folder / '1_raw_data' / 'data' / '*.fastq*')))
        n_demultiplexed = len(glob.glob(str(project_folder / '2_index_demultiplexing' / 'data' / '*.fastq*')))
        n_fasta = len(glob.glob(str(project_folder / '6_read_table' / 'data' / '*.fasta')))
        n_databases = len(glob.glob(str(database_folder / '*')))

        st.write(f"### Active project: {project_folder.name}")
        col1, col2, col3, col4 = st.columns(4)
        col1.metric('Raw files', f'{n_raw:,}')
        col2.metric('Demultiplexed files', f'{n_demultiplexed:,}')
        col3.metric('Read table FASTA files', f'{n_fasta:,}')
        col4.metric('Databases', f'{n_databases:,}')
    else:
        st.info('👈 Start by entering your APSCALE projects folder in the sidebar and select a project.')

    st.markdown("""
    #### Workflow
    1. **🧬 APSCALE-Nanopore** – adjust the settings and run the raw data processing (incl. taxonomic assignment).
    2. **📈 Quality Report** – check the read quality of your raw, demultiplexed or primer-trimmed data.
    3. **📋 Summary Report** – view the summary report of your project.
    """)

    ############################################################################################################
    st.divider()
    st.markdown("## 📚 Links and Tutorials")
    col1, col2 = st.columns(2)
    with col1:
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
    with col2:
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

########################################################################################################################
# APSCALE-Nanopore (raw data processing)
def page_apscale_nanopore():
    st.title('🧬 APSCALE-Nanopore')
    st.markdown(
        "Adjust the settings of your project and run the raw data processing. Raw files are expected in `1_raw_data/data`. "
        "The settings are stored in the project's settings file and are **saved when you start the analysis**.")
    if not project_selected():
        return

    check_dependencies()

    # Databases
    available_databases = {Path(i).name: Path(i) for i in glob.glob(str(database_folder / '*')) if Path(i).is_dir()}

    with st.expander(f'🗄️ APSCALE databases ({len(available_databases)} available)', expanded=available_databases == {}):
        col1, col2, col3 = st.columns(3)
        with col1:
            if st.button('🌐 Open Database Hub', width='stretch'):
                webbrowser.open('https://seafile.rlp.net/d/c172d076de1e4c45b594/')
        with col2:
            if st.button('📥 Download All Latest Databases', width='stretch'):
                public_url = 'https://seafile.rlp.net/d/c172d076de1e4c45b594/?p=%2FLatest&mode=list'
                download_seafile_zip(public_url, database_folder)
        with col3:
            if st.button('📂 Open APSCALE database folder', width='stretch'):
                open_folder(database_folder)

    raw_data = glob.glob(str(project_folder / '1_raw_data' / 'data' / '*.fastq*'))
    n_raw_data = len(raw_data)
    if n_raw_data == 0:
        st.error(f'No raw data files were found! Please copy your .fastq files to "{project_folder}/1_raw_data/data"')
        return

    st.success(f'Found {n_raw_data} raw data files!')

    settings_xlsx = project_folder / Path(str(project_name) + '_settings.xlsx')
    if not settings_xlsx.exists():
        st.warning(f'⚠️ Could not find the settings file: "{settings_xlsx}"')
        return

    # Read Settings file
    settings_df = pd.read_excel(settings_xlsx, sheet_name='Settings').fillna('')
    demultiplexing_df = pd.read_excel(settings_xlsx, sheet_name='Demultiplexing').fillna('')

    # Streamlit forgets the widget values when you switch pages: then simply read the settings file again
    if st.session_state.get("_last_settings_xlsx") != settings_xlsx or 'minimum quality' not in st.session_state:
        st.session_state["_last_settings_xlsx"] = settings_xlsx
        read_settings_file(settings_xlsx, settings_df, demultiplexing_df)

    ############################################################################################################
    st.write('#### 1) Settings')
    col_left, col_right = st.columns(2)
    with col_left:
        with st.expander("⚙️ General settings", expanded=True):
            col1, col2 = st.columns(2)
            with col1:
                # Values are text, because the settings file is read as text (otherwise e.g. "8" is not found)
                cpu_values = ['auto'] + [str(i) for i in range(1, multiprocessing.cpu_count()+1)][::-1]
                if st.session_state.get('cpu count') not in cpu_values:
                    st.session_state['cpu count'] = 'auto'
                st.selectbox(label='CPU count', key='cpu count', options=cpu_values)
            with col2:
                st.selectbox(label='APSCALE-Nanopore mode', key='nanopore mode', options=['linear', 'live base-calling'])

        with st.expander("🧬 Index Demultiplexing", expanded=True):
            st.selectbox(label='Skip demultiplexing', key='sd', options=[False, True])
            col1, col2 = st.columns(2)
            with col1:
                st.text_input(label='Number of allowed errors', key='allowed errors index')
            with col2:
                st.info('Please adjust the Demultiplexing sheet according to your dataset! (Only the first 10 rows are shown)')
            if st.session_state['sd'] == False:
                st.dataframe(demultiplexing_df[['ID', 'Forward index 5-3', 'Reverse index 5-3']].head(10))
                if st.button('📄 Open demultiplexing sheet', width='stretch'):
                    open_file(settings_xlsx)

        with st.expander("✂️ Primer Trimming", expanded=True):
            col1, col2 = st.columns(2)
            with col1:
                st.text_input(label='Number of allowed errors', key='allowed errors primer')
            with col2:
                st.info('Please adjust the Demultiplexing sheet according to your dataset! (Only the first 10 rows are shown)')
            st.dataframe(demultiplexing_df[['ID', 'Forward primer 5-3', 'Reverse primer 5-3', 'Minimum length','Maximum length']].head(10))
            if st.session_state['sd'] == True:
                col3, col4 = st.columns(2)
                with col3:
                    st.info(f'Found {n_raw_data} raw files.')
                with col4:
                    if st.button('Import raw data and update "Demultiplexing" sheet', width='stretch'):
                        import_raw_data(project_folder, settings_xlsx, demultiplexing_df)
                        st.rerun() # show the updated sheet
            if st.button('📄 Open demultiplexing sheet ', width='stretch'):
                open_file(settings_xlsx)

    with col_right:
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
                st.info('Swarm always uses the seed (= centroid) as representative sequence.')
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
                try:
                    cutoff = float(st.session_state['minimum reads'])
                except ValueError:
                    cutoff = None
                    st.error('Please enter a number!')
                if cutoff is not None and cutoff < 1:
                    st.info('Relative threshold detected.')
                    st.info(f'Example: For 60,000 reads → cutoff = {60000 * cutoff:.0f} reads')
                elif cutoff is not None:
                    st.info('Absolute threshold detected.')

        with st.expander("🎯 Taxonomic Assignment", expanded=True):
            st.selectbox(label='Perform taxonomic assignment', key='apscale blast', options=['Yes', 'No'])
            st.selectbox(label='Task', key='task', options=['megablast', 'blastn', 'dc-megablast'])
            # Pre-select the database of the settings file
            database_names = sorted(available_databases.keys())
            last_db = Path(str(st.session_state.get('apscale db', ''))).name
            if 'selected_db' not in st.session_state and last_db in database_names:
                st.session_state['selected_db'] = last_db
            selected_db = st.selectbox(label='Database', key='selected_db', options=database_names)

    if st.button('📄 Open settings file', width='stretch'):
        open_file(settings_xlsx)

    ############################################################################################################
    st.write('#### 2) Run APSCALE-Nanopore')
    with st.container(border=True):
        if st.session_state['nanopore mode'] == 'live base-calling':
            st.warning('Live base-calling is currently only possible through the command line.')
            return
        if selected_db == None and st.session_state['apscale blast'] == 'Yes':
            st.warning('Please select a database! (see "🗄️ APSCALE databases" above)')
            return
        if selected_db != None and ' ' in str(available_databases[selected_db]) and st.session_state['apscale blast'] == 'Yes':
            st.error(
                "Invalid database path: it contains spaces.\n"
                "blastn cannot handle paths with whitespace.\n"
                "Please choose a path without spaces and try again.")
            return

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
                st.warning(
                    '⚠️ Demultiplexing is disabled. Raw FASTQ files will be copied directly '
                    'to the folder **2_index_demultiplexing** without processing.'
                )
            st.info(
                "✅ All required information is available. "
                "Before starting raw data processing, please review the **Demultiplexing** sheet."
            )

        if st.button('▶️ Run APSCALE-Nanopore', type='primary', width='stretch'):
            st.info('Starting APSCALE-Nanopore! Please refer to the terminal for live outputs!')
            if st.session_state['apscale blast'] == 'Yes':
                st.session_state['apscale db'] = available_databases[selected_db]
            update_settings_file(settings_xlsx, settings_df, demultiplexing_df)
            settings_df = pd.read_excel(settings_xlsx, sheet_name='Settings')
            demultiplexing_df = pd.read_excel(settings_xlsx, sheet_name='Demultiplexing')
            apscale_nanopore(project_folder, settings_df, demultiplexing_df, active_steps, st.session_state['sd'])
            st.success('Finished APSCALE-Nanopore!')

########################################################################################################################
# Quality Report
def page_quality_report():
    st.title('📈 Quality Report')
    st.markdown(
        "Create quality reports (PHRED scores and read lengths) for your raw, demultiplexed or primer-trimmed data. "
        "The reports are written to `8_nanopore_report/data`.")
    if not project_selected():
        return

    st.write('#### 1) Generate report')
    with st.container(border=True):
        st.selectbox(label='Select folder', key='qc_folder', options=["1_raw_data", "2_index_demultiplexing", "3_primer_trimming"])
        if st.button(f'📈 Generate report for "{st.session_state["qc_folder"]}"', type='primary', width='stretch'):
            quality_report(project_folder, st.session_state['qc_folder'])
            st.success('Finished creating report!')

    st.write('#### 2) Display report')
    output_folder = project_folder / '8_nanopore_report' / 'data'
    available_reports = {Path(i).name: Path(i) for i in glob.glob(str(output_folder / '*')) if Path(i).is_dir()}
    if not available_reports:
        st.info('No quality reports found. Please generate a report first.')
        return

    plots = {
        'Mean PHRED distribution': 'mean_phred_distribution',
        'Quality by length-bin': 'quality_by_lengthbin',
        'Quality vs. length': 'quality_vs_length',
        'Read length distribution': 'read_length_distribution'
    }
    col1, col2 = st.columns(2)
    with col1:
        st.selectbox(label='Available reports', key='report_to_display', options=available_reports.keys())
    with col2:
        st.selectbox(label='Select plot to display', key='report_display', options=[None] + list(plots.keys()))

    files = glob.glob(str(available_reports[st.session_state['report_to_display']] / '*.html'))
    if files and st.session_state['report_display'] != None:
        selected_files = [i for i in files if plots[st.session_state['report_display']] in i]
        if selected_files:
            with open(selected_files[0], "r", encoding="utf-8") as f:
                html_content = f.read()
            st.components.v1.html(html_content, height=550, scrolling=True)
        else:
            st.info('This plot is not available for the selected report.')

########################################################################################################################
# Summary Report
def page_summary_report():
    st.title('📋 Summary Report')
    st.markdown('The summary report is automatically created and updated when running APSCALE-Nanopore.')
    if not project_selected():
        return

    output_folder = project_folder / '8_nanopore_report' / 'main'
    available_reports = {Path(i).name: Path(i) for i in glob.glob(str(output_folder / '*.html'))}
    if not available_reports:
        st.info('No summary report found. Please run APSCALE-Nanopore first ("🧬 APSCALE-Nanopore").')
        return

    for file in available_reports.values():
        if file.exists():
            with open(file, "r", encoding="utf-8") as f:
                html_content = f.read()
            st.components.v1.html(html_content, height=800, scrolling=True)

########################################################################################################################
# Navigation

pages = {
    'Start': [
        st.Page(page_home, title='Home', icon='🏠', url_path='home', default=True),
    ],
    'Raw Data Processing': [
        st.Page(page_apscale_nanopore, title='APSCALE-Nanopore', icon='🧬', url_path='apscale_nanopore'),
    ],
    'Reports': [
        st.Page(page_quality_report, title='Quality Report', icon='📈', url_path='quality_report'),
        st.Page(page_summary_report, title='Summary Report', icon='📋', url_path='summary_report'),
    ],
}

pg = st.navigation(pages, position='sidebar', expanded=True)
pg.run()
