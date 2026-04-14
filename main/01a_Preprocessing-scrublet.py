import os
import scipy.io
import pandas as pd
import scrublet as scr

scrublet_dir = "scrublet"

counts_path = os.path.join(scrublet_dir, "counts.mtx")
barcodes_path = os.path.join(scrublet_dir, "barcodes.csv")
params_path = os.path.join(scrublet_dir, "scrublet_params.csv")
out_path = os.path.join(scrublet_dir, "scrublet_results.csv")

counts_matrix = scipy.io.mmread(counts_path).tocsr()
barcodes = pd.read_csv(barcodes_path)
params = pd.read_csv(params_path)

expected_doublet_rate = float(params["expected_doublet_rate"].iloc[0])

scrub = scr.Scrublet(counts_matrix, expected_doublet_rate=expected_doublet_rate)

doublet_scores, predicted_doublets = scrub.scrub_doublets(
    min_counts=2,
    min_cells=3,
    min_gene_variability_pctl=85,
    n_prin_comps=30
)

barcodes["scrublet_score"] = doublet_scores
barcodes["scrublet_call"] = predicted_doublets
barcodes.to_csv(out_path, index=False)

print(f"Wrote results to {out_path}")