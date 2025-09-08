import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os

# -----------------------------
# Inputs
# -----------------------------
VST_PATH = "DESEQ2_NORM_all.vst.tsv"
META_PATH = "../../data/differential_abundance_sheets/rnaseq_diffabundance_study1and5_samplesheet_filled.csv"
PANEL_PATH = "final_panel_gene_list.txt"

# -----------------------------
# Load data
# -----------------------------
vst = pd.read_csv(VST_PATH, sep="\t")
vst = vst.rename(columns={vst.columns[0]: "gene_id"}).set_index("gene_id")

meta = pd.read_csv(META_PATH)
# Clean header names
meta.columns = [c.replace("\ufeff", "").strip() for c in meta.columns]
# Normalize names
ren = {}
for c in meta.columns:
    cl = c.lower()
    if cl == "sample": ren[c] = "sample"
    elif cl == "condition": ren[c] = "condition"
    elif cl == "batch": ren[c] = "batch"
meta = meta.rename(columns=ren)

# Ensure required columns
assert {"sample","condition","batch"}.issubset(meta.columns)

meta["sample"] = meta["sample"].astype(str).str.strip()
meta["condition"] = meta["condition"].astype(str).str.strip().str.lower()
meta["batch"] = meta["batch"].astype(str).str.strip()

# Panel genes
genes = pd.read_csv(PANEL_PATH, header=None)[0].tolist()

# -----------------------------
# Align and subset
# -----------------------------
genes_present = [g for g in genes if g in vst.index]
expr_sub = vst.loc[genes_present].T
expr_sub.index.name = "sample"

df = expr_sub.join(meta.set_index("sample"), how="inner")

print(f"[INFO] Plotting for {len(genes_present)} genes across {df.shape[0]} samples.")
print(f"[INFO] Conditions: {sorted(df['condition'].unique())}")
print(f"[INFO] Batches: {sorted(df['batch'].unique())}")

# -----------------------------
# Output directory
# -----------------------------
OUTDIR = "./plots"
os.makedirs(OUTDIR, exist_ok=True)

# -----------------------------
# Helper: violin per gene
# -----------------------------
def save_violin_per_gene(df, gene, outpath):
    cond_order = sorted(df["condition"].unique())
    data = [df.loc[df["condition"] == c, gene].dropna().values for c in cond_order]

    plt.figure(figsize=(6, 4))
    # Violin plots
    plt.violinplot(data, showmeans=True, showextrema=True, showmedians=False)

    # Scatter overlay: split by batch
    batches = sorted(df["batch"].unique())
    x_positions = {c: i+1 for i, c in enumerate(cond_order)}
    rng = np.random.RandomState(42)
    for b in batches:
        sub = df[df["batch"] == b]
        for c in cond_order:
            yvals = sub.loc[sub["condition"] == c, gene].dropna().values
            if len(yvals) == 0:
                continue
            x = np.full_like(yvals, x_positions[c], dtype=float)
            x = x + (rng.rand(len(yvals)) - 0.5) * 0.15
            plt.scatter(x, yvals, s=20, alpha=0.7, linewidth=0.3)

    plt.xticks(range(1, len(cond_order)+1), cond_order)
    plt.xlabel("condition")
    plt.ylabel("VST expression")
    plt.title(f"{gene} by condition (points grouped by batch)")
    plt.tight_layout()
    plt.savefig(outpath, dpi=300)
    plt.close()

# -----------------------------
# Create violin plots
# -----------------------------
for g in genes_present:
    out = os.path.join(OUTDIR, f"violin_{g}.png")
    save_violin_per_gene(df, g, out)

# -----------------------------
# Heatmap + annotation bars
# -----------------------------
df_sorted = df.sort_values(["condition", "batch"])
H = df_sorted[genes_present].T.values
gene_labels = genes_present

# Heatmap
plt.figure(figsize=(max(8, len(df_sorted)*0.08), 4))
plt.imshow(H, aspect='auto')
plt.colorbar(label="VST expression")
plt.yticks(range(len(gene_labels)), gene_labels)
plt.xticks([])
plt.title("Final panel genes across samples (sorted by condition, batch)")
plt.tight_layout()
plt.savefig(os.path.join(OUTDIR, "heatmap_panel_genes.png"), dpi=300)
plt.close()

# Condition annotation bar
cond_codes, cond_uniques = pd.factorize(df_sorted["condition"], sort=True)
cond_bar = cond_codes[np.newaxis, :]
plt.figure(figsize=(max(8, len(df_sorted)*0.08), 1))
plt.imshow(cond_bar, aspect='auto')
plt.yticks([0], ["condition"])
plt.xticks([])
plt.title("Annotation: condition")
plt.tight_layout()
plt.savefig(os.path.join(OUTDIR, "annotation_condition.png"), dpi=300)
plt.close()

# Batch annotation bar
batch_codes, batch_uniques = pd.factorize(df_sorted["batch"], sort=True)
batch_bar = batch_codes[np.newaxis, :]
plt.figure(figsize=(max(8, len(df_sorted)*0.08), 1))
plt.imshow(batch_bar, aspect='auto')
plt.yticks([0], ["batch"])
plt.xticks([])
plt.title("Annotation: batch")
plt.tight_layout()
plt.savefig(os.path.join(OUTDIR, "annotation_batch.png"), dpi=300)
plt.close()

# Legend mapping
with open(os.path.join(OUTDIR, "annotation_legends.txt"), "w") as f:
    f.write("Condition legend:\n")
    for i, c in enumerate(cond_uniques):
        f.write(f"  code {i} -> {c}\n")
    f.write("Batch legend:\n")
    for i, b in enumerate(batch_uniques):
        f.write(f"  code {i} -> {b}\n")

print(f"Plots written to: {OUTDIR}")
