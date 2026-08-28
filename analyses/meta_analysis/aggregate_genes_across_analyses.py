import pandas as pd
from pathlib import Path
from collections import defaultdict

dir_path = Path(__file__).parent

# Map each file to its gene column
sources = [
    ("IDA_condition_sensitive_tolerant.deseq2.results_filtered.tsv", "gene_id",   "\t"),
    ("LitReview_cross_study_DEG_comparison.csv",                      "LOC_ID",    ","),
    ("PADI_ds1-5_allDEGs_merged.csv",                                 "gene_id",   ","),
    ("6-gene-LOSO-panal_final_consensus_genes.csv",                   "gene",      ","),
    ("MutualInfo_top10genes.csv",                                      "gene_id",   ","),
    ("RandomForest_top_50_predictive_genes_for_classification.csv",    "gene_id",   ","),
]

# plain-text file, one gene per line
panel_genes = (dir_path / "6-gene-panal_final_panel_gene_list.txt").read_text().splitlines()
panel_genes = [g.strip() for g in panel_genes if g.strip()]

gene_analyses = defaultdict(list)

# Plain-text file
for gene in panel_genes:
    gene_analyses[gene].append("6-gene-panal_final_panel_gene_list")

# CSV/TSV files
for filename, gene_col, sep in sources:
    analysis_name = Path(filename).stem
    df = pd.read_csv(dir_path / filename, sep=sep, usecols=[gene_col], dtype=str)
    for gene in df[gene_col].dropna().unique():
        gene_analyses[gene.strip()].append(analysis_name)

# Build output table
rows = []
for gene, analyses in gene_analyses.items():
    rows.append({
        "gene_id": gene,
        "n_analyses": len(analyses),
        "analyses": "; ".join(analyses),
    })

result = (
    pd.DataFrame(rows)
    .sort_values(["n_analyses", "gene_id"], ascending=[False, True])
    .reset_index(drop=True)
)

out_path = dir_path / "all_genes_across_analyses.csv"
result.to_csv(out_path, index=False)
print(f"Wrote {len(result)} genes to {out_path}")
print(result["n_analyses"].value_counts().sort_index(ascending=False).to_string())
