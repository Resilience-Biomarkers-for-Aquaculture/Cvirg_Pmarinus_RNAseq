#!/usr/bin/env python3
"""
Plot stability-selection probabilities for genes common to 3 CSVs
above a threshold in ALL three, with files living in separate subdirs.

Defaults are set for your current directory layout.
"""

from __future__ import annotations

import argparse
from pathlib import Path
import sys
import pandas as pd
import matplotlib.pyplot as plt


def _infer_cols(df: pd.DataFrame) -> tuple[str, str]:
    cols = [c.lower() for c in df.columns]
    lower_to_orig = {c.lower(): c for c in df.columns}

    gene_candidates = [c for c in cols if "gene" in c]
    prob_candidates = [c for c in cols if "prob" in c]

    if len(gene_candidates) != 1 or len(prob_candidates) != 1:
        raise ValueError(
            f"Could not uniquely infer columns. "
            f"gene_candidates={gene_candidates}, "
            f"prob_candidates={prob_candidates}, "
            f"all_cols={list(df.columns)}"
        )

    return (
        lower_to_orig[gene_candidates[0]],
        lower_to_orig[prob_candidates[0]],
    )


def load_fold_csv(csv_path: Path, fold_name: str) -> pd.DataFrame:
    df = pd.read_csv(csv_path)
    gene_col, prob_col = _infer_cols(df)

    out = df[[gene_col, prob_col]].copy()
    out = out.rename(columns={gene_col: "gene", prob_col: fold_name})

    # sanity / shape checks
    assert out.shape[1] == 2
    assert out["gene"].notna().all()
    assert out[fold_name].notna().all()

    out[fold_name] = pd.to_numeric(out[fold_name], errors="raise")

    # defensive deduplication
    if out["gene"].duplicated().any():
        out = out.groupby("gene", as_index=False)[fold_name].max()

    return out


def build_intersection(folds: dict[str, pd.DataFrame], cutoff: float) -> pd.DataFrame:
    fold_names = list(folds.keys())
    merged = folds[fold_names[0]]

    for fn in fold_names[1:]:
        merged = merged.merge(folds[fn], on="gene", how="inner")

    # all folds must be >= cutoff
    for fn in fold_names:
        merged = merged[merged[fn] >= cutoff]

    # at least one fold must be >= 0.2
    merged = merged[
        merged[fold_names].max(axis=1) >= 0.2
    ]

    return merged.sort_values("gene").reset_index(drop=True)


def plot_lines(df: pd.DataFrame, fold_order: list[str], cutoff: float, out_png: Path | None):
    plt.figure()

    for _, row in df.iterrows():
        plt.plot(
            fold_order,
            [row[fn] for fn in fold_order],
            marker="o",
            label=row["gene"],
        )

    plt.ylabel("Selection probability")
    plt.title(f"Genes common to all three folds (p ≥ {cutoff})")
    plt.legend(bbox_to_anchor=(1.05, 1), loc="upper left")
    plt.tight_layout()

    if out_png:
        out_png.parent.mkdir(parents=True, exist_ok=True)
        plt.savefig(out_png, dpi=200, bbox_inches="tight")
    else:
        plt.show()


def main(argv: list[str]) -> int:
    p = argparse.ArgumentParser()

    # ---- defaults set here ----
    p.add_argument("--base", type=Path, default=Path("results_08_Dec_2025"))
    p.add_argument("--subdir1", default="loso_P&S_2023")
    p.add_argument("--subdir2", default="loso_P&S_2020_2017")
    p.add_argument("--subdir3", default="loso_P&S_2020_2015")

    p.add_argument("--file1", default="stability_selection_P&S 2023.csv")
    p.add_argument("--file2", default="stability_selection_P&S 2020 2017.csv")
    p.add_argument("--file3", default="stability_selection_P&S 2020 2015.csv")

    p.add_argument("--name1", default="P&S_2023")
    p.add_argument("--name2", default="P&S_2020_2017")
    p.add_argument("--name3", default="P&S_2020_2015")

    p.add_argument("--cutoff", type=float, default=0.001)
    p.add_argument("--out", type=Path, default=Path("plots/common_genes_p001.png"))
    # ---------------------------

    args = p.parse_args(argv)

    csv1 = args.base / args.subdir1 / args.file1
    csv2 = args.base / args.subdir2 / args.file2
    csv3 = args.base / args.subdir3 / args.file3

    for path in (csv1, csv2, csv3):
        if not path.exists():
            raise FileNotFoundError(f"Missing file: {path}")

    folds = {
        args.name1: load_fold_csv(csv1, args.name1),
        args.name2: load_fold_csv(csv2, args.name2),
        args.name3: load_fold_csv(csv3, args.name3),
    }

    inter = build_intersection(folds, cutoff=args.cutoff)

    print(f"Common genes with p ≥ {args.cutoff}: {len(inter)}")
    if len(inter):
        print(inter[["gene", args.name1, args.name2, args.name3]].to_string(index=False))

    plot_lines(
        inter,
        fold_order=[args.name1, args.name2, args.name3],
        cutoff=args.cutoff,
        out_png=args.out,
    )

    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
