"""
Create a combined biomarker table from all analyses.

Output: all_analyses_biomarkers.csv
Columns:
  gene_id     - gene identifier
  n_analyses  - number of analyses the gene was found in
  analyses    - semicolon-separated list of analysis names (from file names)
"""

import os
import urllib.request
import pandas as pd
from collections import defaultdict

BASE = os.path.dirname(os.path.abspath(__file__))
ANALYSES_DIR = os.path.join(BASE, "..")

# ---------------------------------------------------------------------------
# Define each gene-list source as (label, path_or_url, gene_col)
# ---------------------------------------------------------------------------
SOURCES = [
    (
        "top10genes",
        os.path.join(ANALYSES_DIR, "2025-01-03_RNAseq_all_AI_diffexp", "top10genes.csv"),
        "gene_id",
    ),
    (
        "final_panel_gene_list",
        os.path.join(
            ANALYSES_DIR,
            "Study1and5ThreeWay",
            "two_step_gene_expression_classifier",
            "final_panel_gene_list.txt",
        ),
        None,  # plain list, no header
    ),
    (
        "final_consensus_genes",
        os.path.join(
            ANALYSES_DIR,
            "Study1and5ThreeWay",
            "two_step_gene_expression_classifier",
            "results_08_Dec_2025",
            "final_consensus_genes.csv",
        ),
        "gene",
    ),
    (
        "top_50_predictive_genes_for_classification",
        "https://raw.githubusercontent.com/Resilience-Biomarkers-for-Aquaculture/"
        "project-predicting-phenotype/refs/heads/main/data/processed/"
        "top_50_predictive_genes_for_classification.csv",
        "gene_id",
    ),
    (
        "condition_sensitive_tolerant.deseq2.results_filtered",
        "https://gannet.fish.washington.edu/metacarcinus/USDA_MetaOmics/"
        "Cvirg_RNAseq/20260418_diffabund/allds/tables/differential/"
        "condition_sensitive_tolerant.deseq2.results_filtered.tsv",
        "gene_id",
    ),
    (
        "ds1-5_allDEGs_merged",
        os.path.join(ANALYSES_DIR, "20260418_diffAbund_allData", "ds1-5_allDEGs_merged.csv"),
        "gene_id",
    ),
    (
        "cross_study_DEG_comparison",
        os.path.join(ANALYSES_DIR, "LitReview", "cross_study_DEG_comparison.csv"),
        "LOC_ID",
    ),
]


def load_genes(label, path_or_url, gene_col):
    """Return a set of gene IDs for one source, or None if unavailable."""
    is_url = path_or_url.startswith("http")
    try:
        if is_url:
            req = urllib.request.urlopen(path_or_url, timeout=30)
            import io
            content = req.read().decode("utf-8")
            sep = "\t" if path_or_url.endswith(".tsv") else ","
            df = pd.read_csv(io.StringIO(content), sep=sep)
        else:
            sep = "\t" if path_or_url.endswith(".tsv") else ","
            if gene_col is None:
                # Plain list with no header
                df = pd.read_csv(path_or_url, header=None, names=["gene_id"])
                gene_col = "gene_id"
            else:
                df = pd.read_csv(path_or_url, sep=sep)
        genes = set(df[gene_col].dropna().astype(str).str.strip())
        print(f"  {label}: {len(genes)} genes loaded")
        return genes
    except Exception as exc:
        print(f"  WARNING – could not load {label}: {exc}")
        return None


def main():
    gene_to_analyses = defaultdict(list)

    for label, path_or_url, gene_col in SOURCES:
        genes = load_genes(label, path_or_url, gene_col)
        if genes is None:
            continue
        for gene in genes:
            gene_to_analyses[gene].append(label)

    rows = []
    for gene, analyses_list in gene_to_analyses.items():
        rows.append(
            {
                "gene_id": gene,
                "n_analyses": len(analyses_list),
                "analyses": "; ".join(sorted(analyses_list)),
            }
        )

    df_out = pd.DataFrame(rows).sort_values(
        ["n_analyses", "gene_id"], ascending=[False, True]
    ).reset_index(drop=True)

    out_path = os.path.join(BASE, "all_analyses_biomarkers.csv")
    df_out.to_csv(out_path, index=False)
    print(f"\nSaved {len(df_out)} genes to {out_path}")
    print(df_out["n_analyses"].value_counts().sort_index(ascending=False).to_string())


if __name__ == "__main__":
    main()
