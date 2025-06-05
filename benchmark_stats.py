import pandas as pd
import plotly.graph_objects as go
import numpy as np

df = pd.read_excel('/Users/tillmacher/Desktop/APSCALE_projects/test_dataset_apscale_nanopore/9_benchmark/benchmark_results.xlsx')

tests = ['ESVs', 'Species']

for test in tests:
    indices = []
    for _, row in df.iterrows():
        id = row['id']
        n_new = row[f'{test}_only']           # Detected only by Nanopore
        n_shared = row[f'{test}_shared']      # Detected in both reference and Nanopore
        n_ref_only = row[f'{test}_ref']       # Detected only in reference
        time = row['elapsed_time']            # Time in seconds

        # 1. Consistency Score — favors shared values, penalizes new
        consistency_score = (n_shared - n_new) * n_shared

        # 2. Efficiency Score — same as above, normalized by time
        consistency_efficiency_score = consistency_score / time if time > 0 else 0

        # 3. Log-Penalty Score — penalizes new detections logarithmically
        log_penalty_score = n_shared * np.log1p(n_shared / (n_new + 1))

        # 4. Purity Score — how "clean" the detection is
        purity_score = n_shared**2 / (n_shared + n_new) if (n_shared + n_new) > 0 else 0

        # 5. Jaccard Index — proportion of overlap
        jaccard_index = n_shared / (n_shared + n_new + n_ref_only) if (n_shared + n_new + n_ref_only) > 0 else 0

        # 6. Weighted Jaccard Score — accounts for volume
        weighted_jaccard_score = jaccard_index * n_shared

        # 7. Composite Score — combining consistency and efficiency
        composite_score = consistency_score * consistency_efficiency_score / time if time > 0 else 0

        # 8. Normalized Recall Score — shared over all reference detections
        recall_score = n_shared / (n_shared + n_ref_only) if (n_shared + n_ref_only) > 0 else 0

        # 9. Detection Ratio Score — shared over all nanopore detections
        precision_score = n_shared / (n_shared + n_new) if (n_shared + n_new) > 0 else 0

        # 10. F1-style index — harmonic mean of recall and precision
        f1_score = (2 * recall_score * precision_score / (recall_score + precision_score)
                    if (recall_score + precision_score) > 0 else 0)

        indices.append([
            id,
            n_new,
            n_shared,
            n_ref_only,
            time,
            consistency_score,
            consistency_efficiency_score,
            log_penalty_score,
            purity_score,
            jaccard_index,
            weighted_jaccard_score,
            composite_score,
            recall_score,
            precision_score,
            f1_score
        ])


    cols = [
            "id",
            "n_new",
            "n_shared",
            "n_ref_only",
            "time",
            "consistency_score",
            "consistency_efficiency_score",
            "log_penalty_score",
            "purity_score",
            "jaccard_index",
            "weighted_jaccard_score",
            "composite_score",
            "recall_score",
            "precision_score",
            "f1_score"
        ]

    res_df = pd.DataFrame(indices, columns=cols)

    scores = cols[5:]
    ranks = {i:[] for i in res_df['id']}
    for score in scores:
        # Get ranking: {id: rank}
        ranked_ids = res_df.sort_values(score, ascending=False)['id'].values.tolist()
        sub_df = {j: i + 1 for i, j in enumerate(ranked_ids)}

        # Append rank for this score to each id
        for i in ranks:
            ranks[i].append(sub_df.get(i, None))

    ranks_df = pd.DataFrame(ranks).transpose()
    ranks_df_mean = ranks_df.median(axis=1).reset_index().sort_values(by=0)

    res_df.to_excel(f'/Users/tillmacher/Desktop/APSCALE_projects/test_dataset_apscale_nanopore/9_benchmark/benchmark_stats_{test}.xlsx', index=False)
    ranks_df_mean.to_excel(f'/Users/tillmacher/Desktop/APSCALE_projects/test_dataset_apscale_nanopore/9_benchmark/benchmark_ranks_{test}.xlsx', index=False)


    import plotly.graph_objects as go
    from plotly.subplots import make_subplots

    # Sort by F1 score and prepare data
    heatmap_df = res_df.copy().sort_values('f1_score', ascending=True)
    y_values = heatmap_df['id'].tolist()
    scores = heatmap_df.columns[5:]  # or use your defined `scores` list

    # Create subplots — one column per score
    fig = make_subplots(
        rows=1,
        cols=len(scores),
        shared_yaxes=True,
        horizontal_spacing=0.01,
    )

    for i, score in enumerate(scores):
        z = [[val] for val in heatmap_df[score]]  # shape (n, 1)

        # Default: no fixed range
        zmin = None
        zmax = None
        cscale = 'Blues'

        # Set fixed color scale range for selected metrics
        if score == 'jaccard_index' or score == 'recall_score' or score == 'precision_score' or score == 'f1_score':
            zmin = 0
            zmax = 1
            cscale = 'Greens'

        fig.add_trace(
            go.Heatmap(
                z=z,
                x=[score],
                y=y_values,
                colorscale=cscale,
                zmin=zmin,
                zmax=zmax,
                colorbar=dict(title=score) if i == len(scores) - 1 else None,
                showscale=(i == len(scores) - 1)
            ),
            row=1,
            col=i + 1
        )

    # Update layout
    fig.update_layout(
        template='plotly_white',
        height=15 * len(y_values),
        width=200 * len(scores),  # dynamically scale width
        showlegend=False
    )

    fig.update_yaxes(dtick='linear')

    # Save to file
    fig.write_image(f'/Users/tillmacher/Desktop/APSCALE_projects/test_dataset_apscale_nanopore/9_benchmark/Scores_{test}.pdf')
    fig.write_html(f'/Users/tillmacher/Desktop/APSCALE_projects/test_dataset_apscale_nanopore/9_benchmark/Scores_{test}.html')

