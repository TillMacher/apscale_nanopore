import pandas as pd
import plotly.graph_objects as go
import numpy as np
from pathlib import Path
import glob
from plotly.subplots import make_subplots

# === CONFIGURATION ===
BASE_PROJECT_PATH = Path('/Users/tillmacher/Desktop/APSCALE_projects/test_dataset_apscale_nanopore')
BENCHMARK_FOLDER = BASE_PROJECT_PATH / '9_benchmark'
REFERENCE_FILE = BASE_PROJECT_PATH / '9_benchmark_old/reference_taxon_table.xlsx'
existing_ids = [p.name for p in BENCHMARK_FOLDER.iterdir() if p.is_dir()]

# === INITIATE FINAL RESULTS CONTAINER ===
heatmap_df = pd.DataFrame()

# === MAIN LOOP FOR ALL IDS ===
for name in existing_ids:
    res_folder = BENCHMARK_FOLDER / name
    taxonomy_table = res_folder / f'{name}_taxonomy.xlsx'
    read_table = res_folder / f'{name}_read_table.xlsx'
    runtime_file = res_folder / f'{name}_time.xlsx'

    taxonomy_table_df = pd.read_excel(taxonomy_table).fillna('')
    read_table_df = pd.read_excel(read_table).fillna('')
    runtime = pd.read_excel(runtime_file).values.tolist()[0][0]

    samples = read_table_df.columns.tolist()[2:]
    nanopore_ids_dict, nanopore_species_dict = {}, {}

    for sample in samples:
        present_ids = read_table_df[read_table_df[sample] != 0]['ID'].tolist()
        present_species = sorted(set(
            taxonomy_table_df[taxonomy_table_df['unique ID'].isin(present_ids)]['Species']
        ))
        present_species = [i for i in present_species if i not in ['', 'NoMatch']]
        nanopore_ids_dict[sample] = present_ids
        nanopore_species_dict[sample] = present_species

    # Load reference table
    taxon_table_df = pd.read_excel(REFERENCE_FILE).fillna('')
    ref_ids_dict, ref_species_dict = {}, {}

    for sample in samples:
        present_ids = taxon_table_df[taxon_table_df[sample] != 0]['unique_ID'].tolist()
        present_species = sorted(set(
            taxon_table_df[taxon_table_df['unique_ID'].isin(present_ids)]['Species']
        ))
        present_species = [i for i in present_species if i not in ['', 'NoMatch']]
        ref_ids_dict[sample] = present_ids
        ref_species_dict[sample] = present_species

    # === METRICS COMPUTATION ===
    tests = ['ESVs', 'Species']
    stats_list = []

    for test in tests:
        for sample in samples:
            if test == 'ESVs':
                set_nano = set(nanopore_ids_dict[sample])
                set_ref = set(ref_ids_dict[sample])
            else:
                set_nano = set(nanopore_species_dict[sample])
                set_ref = set(ref_species_dict[sample])

            n_nanopore_only = len(set_nano - set_ref)
            n_shared = len(set_nano & set_ref)
            n_ref_only = len(set_ref - set_nano)

            consistency_score = (n_shared - n_nanopore_only) * n_shared
            consistency_efficiency_score = consistency_score / runtime if runtime > 0 else 0
            log_penalty_score = n_shared * np.log1p(n_shared / (n_nanopore_only + 1))
            purity_score = n_shared**2 / (n_shared + n_nanopore_only) if (n_shared + n_nanopore_only) > 0 else 0
            jaccard_index = n_shared / (n_shared + n_nanopore_only + n_ref_only) if (n_shared + n_nanopore_only + n_ref_only) > 0 else 0
            weighted_jaccard_score = jaccard_index * n_shared
            composite_score = (consistency_score * consistency_efficiency_score / runtime) if runtime > 0 else 0
            recall_score = n_shared / (n_shared + n_ref_only) if (n_shared + n_ref_only) > 0 else 0
            precision_score = n_shared / (n_shared + n_nanopore_only) if (n_shared + n_nanopore_only) > 0 else 0
            f1_score = (2 * recall_score * precision_score / (recall_score + precision_score)
                        if (recall_score + precision_score) > 0 else 0)

            stats_list.append([
                sample, test,
                n_nanopore_only, n_shared, n_ref_only, runtime,
                consistency_score, consistency_efficiency_score,
                log_penalty_score, purity_score,
                jaccard_index, weighted_jaccard_score, composite_score,
                recall_score, precision_score, f1_score
            ])

    # Save to dataframe and file
    cols = [
        "sample", "test",
        "n_nanopore_only", "n_shared", "n_ref_only", "runtime",
        "consistency_score", "consistency_efficiency_score",
        "log_penalty_score", "purity_score",
        "jaccard_index", "weighted_jaccard_score", "composite_score",
        "recall_score", "precision_score", "f1_score"
    ]
    res_df = pd.DataFrame(stats_list, columns=cols)
    res_file = res_folder / f'{name}_stats.xlsx'
    res_df.to_excel(res_file, index=False)

    score_names = cols[5:]
    for test in tests:
        sub_df = res_df[res_df['test'] == test]
        sub_df = sub_df[score_names].mean().reset_index()
        sub_df = pd.DataFrame([name, test] + sub_df[0].values.tolist()).T
        sub_df.columns = ['ID', 'test'] + score_names
        heatmap_df = pd.concat([heatmap_df, sub_df], ignore_index=True)

# === HEATMAP PLOTTING ===
for test in ['ESVs', 'Species']:
    sub_df = heatmap_df[heatmap_df['test'] == test].copy()
    sub_df = sub_df.sort_values(by=['f1_score', 'runtime', 'jaccard_index'])

    y_values = sub_df['ID'].tolist()
    scores = sub_df.columns[3:]

    fig = make_subplots(
        rows=1,
        cols=len(scores),
        shared_yaxes=True,
        horizontal_spacing=0.01
    )

    for i, score in enumerate(scores):
        z_values = sub_df[score].tolist()

        # Validate z-values
        if not all(np.isfinite(z_values)):
            print(f"Warning: Non-finite values found in score '{score}'. Replacing with 0.")
            z_values = [0 if not np.isfinite(v) else v for v in z_values]

        z = [[val] for val in z_values]

        cscale = 'Blues'
        zmin = zmax = None

        if score in ['jaccard_index', 'recall_score', 'precision_score', 'f1_score']:
            cscale = 'Greens'
            zmin, zmax = 0, 1
        elif score == 'runtime':
            cscale = 'Blues_r'

        fig.add_trace(
            go.Heatmap(
                z=z,
                x=[score],
                y=y_values,
                colorscale=cscale,
                zmin=zmin,
                zmax=zmax,
                showscale=(i == len(scores) - 1),
                colorbar=dict(title=score) if i == len(scores) - 1 else None
            ),
            row=1,
            col=i + 1
        )

    fig.update_layout(
        template='plotly_white',
        height=15 * len(y_values) if len(y_values) >= 5 else 500,
        width=100 * len(scores),
        showlegend=False
    )
    fig.update_yaxes(dtick='linear')
    fig.update_xaxes(tickangle=-45)

    fig.write_image(str(BENCHMARK_FOLDER / f'Scores_{test}.pdf'))
    fig.write_html(str(BENCHMARK_FOLDER / f'Scores_{test}.html'))