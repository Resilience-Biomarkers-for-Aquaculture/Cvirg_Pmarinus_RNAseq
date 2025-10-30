# This is from https://chatgpt.com/g/g-p-68dfc7876dc48191a04cbbe0be5c4097-minimal-biomarker/shared/c/68ffdd7a-b62c-8323-8807-a67dc0de3411?owner_user_id=user-AozXqg98DrkNBXzFTBoCy1ep
# It's for step 1.A.
import json
import numpy as np
import pandas as pd
from pathlib import Path
from sklearn.metrics import roc_auc_score, average_precision_score

# --- SANITY CHECKS ---
# We'll use scipy.stats in BCa; import once (minor speedup vs inside function)
import scipy.stats as st

rng = np.random.default_rng(42)


def brier_score(y, p):
    y = np.asarray(y, float); p = np.asarray(p, float)
    return np.mean((p - y) ** 2)


def _safe_metric(metric_name, y, p):
    # Compute metric with guards for degenerate cases
    y = np.asarray(y); p = np.asarray(p)
    if metric_name == "auroc":
        # Need both classes present
        if len(np.unique(y)) < 2:
            return np.nan
        return roc_auc_score(y, p)
    elif metric_name == "aupr":
        if (y == 1).sum() == 0 or (y == 0).sum() == 0:
            return np.nan
        return average_precision_score(y, p)
    elif metric_name == "brier":
        return brier_score(y, p)
    else:
        raise ValueError(metric_name)


def bca_ci_per_fold(y, p, metric_name, B=10000, alpha=0.05, rng=rng):
    """
    BCa CI for one metric in one fold. Returns (theta_hat, lo, hi, degenerate_skipped, total_attempted).
    y: 1D array of 0/1, p: 1D array of predicted probabilities
    """
    y = np.asarray(y); p = np.asarray(p); n = len(y)

    # 1) Observed
    theta_hat = _safe_metric(metric_name, y, p)

    # 2) Bootstrap distribution
    thetas = []
    b_count = 0
    degenerate = 0  # --- SANITY CHECKS ---
    attempts = 0    # --- SANITY CHECKS ---
    while b_count < B:
        idx = rng.integers(0, n, size=n)
        tb = _safe_metric(metric_name, y[idx], p[idx])
        attempts += 1
        if np.isnan(tb):
            degenerate += 1          # --- SANITY CHECKS ---
            continue  # resample again if degenerate
        thetas.append(tb)
        b_count += 1
    thetas = np.sort(np.array(thetas))

    # 3) Jackknife for acceleration 'a'
    # leave-one-out influence values
    jack = []
    for i in range(n):
        mask = np.ones(n, dtype=bool); mask[i] = False
        tj = _safe_metric(metric_name, y[mask], p[mask])
        if np.isnan(tj):  # if degenerate, approximate with observed
            tj = theta_hat
        jack.append(tj)
    jack = np.array(jack)
    jack_mean = jack.mean()
    num = np.sum((jack_mean - jack) ** 3)
    den = 6.0 * (np.sum((jack_mean - jack) ** 2) ** 1.5 + 1e-12)
    a = num / den if den != 0 else 0.0

    # 4) Bias-correction z0
    # proportion of bootstrap thetas less than observed
    prop = (thetas < theta_hat).mean()
    prop = np.clip(prop, 1e-6, 1 - 1e-6)

    # z = Phi^{-1}(prop)
    z0 = st.norm.ppf(prop)

    # 5) Adjusted quantiles
    z_alpha1 = st.norm.ppf(alpha / 2)
    z_alpha2 = st.norm.ppf(1 - alpha / 2)

    def adj_quantile(z):
        num = z0 + z
        den = 1 - a * (z0 + z)
        return st.norm.cdf(z0 + (z0 + z) / den)

    ql = adj_quantile(z_alpha1)
    qh = adj_quantile(z_alpha2)
    lo = np.quantile(thetas, ql)
    hi = np.quantile(thetas, qh)
    return float(theta_hat), float(lo), float(hi), int(degenerate), int(attempts)  # --- SANITY CHECKS ---


# --- SANITY CHECKS ---
def fold_class_balance(preds_csv):
    """Return (n_pos, n_neg, aupr_baseline)"""
    df = pd.read_csv(preds_csv)
    y = df["true_label"].values
    n_pos = int((y == 1).sum())
    n_neg = int((y == 0).sum())
    base_aupr = n_pos / (n_pos + n_neg) if (n_pos + n_neg) > 0 else np.nan
    return n_pos, n_neg, float(base_aupr)


def per_fold_bca_from_json(json_path, preds_key="preds_file", B=2000):
    with open(json_path) as f:
        info = json.load(f)
    preds_csv = info[preds_key]
    # prefix fold path dir of json_path
    preds_csv = str(Path(json_path).parent / preds_csv)
    df = pd.read_csv(preds_csv)
    y = df["true_label"].values
    p = df["pred_prob_tolerant"].values

    # --- SANITY CHECKS ---
    n_pos, n_neg, base_aupr = fold_class_balance(preds_csv)

    out = {}
    degenerate_info = {}  # metric -> (degenerate, attempts)
    for metr in ["auroc", "aupr", "brier"]:
        th, lo, hi, deg, att = bca_ci_per_fold(y, p, metr, B=B)
        out[metr] = {"point": th, "ci95": (lo, hi)}
        degenerate_info[metr] = (deg, att)

    out["n_test"] = int(len(y))
    out["fold"] = info.get("heldout", Path(json_path).stem)

    # --- SANITY CHECKS: JSON metric consistency printouts ---
    for metr in ["auroc", "aupr", "brier"]:
        if metr in info:
            orig = float(info[metr])
            now = float(out[metr]["point"])
            diff = abs(now - orig)
            print(f"[check:{out['fold']}] {metr} JSON={orig:.6f}  recomputed={now:.6f}  Δ={diff:.4g}")

    # --- SANITY CHECKS: print degenerate resamples per metric ---
    for metr in ["auroc", "aupr"]:
        deg, att = degenerate_info[metr]
        pct = 100.0 * deg / max(att, 1)
        print(f"[bootstrap:{out['fold']}] {metr} degenerate {deg}/{att} ({pct:.2f}%) skipped")

    # attach class balance
    out["n_pos"] = n_pos
    out["n_neg"] = n_neg
    out["aupr_baseline"] = base_aupr

    return out


root = Path("results")  # your LOSO root
rows = []
for jf in root.glob("loso_*/best_panel_summary_*.json"):
    res = per_fold_bca_from_json(jf, B=2000)  # bump to 10000 for final
    rows.append({
        "fold": res["fold"],
        "n_test": res["n_test"],
        "n_pos":  res["n_pos"],          # --- SANITY CHECKS ---
        "n_neg":  res["n_neg"],          # --- SANITY CHECKS ---
        "auroc": res["auroc"]["point"],
        "auroc_lo": res["auroc"]["ci95"][0],
        "auroc_hi": res["auroc"]["ci95"][1],
        "aupr":  res["aupr"]["point"],
        "aupr_lo": res["aupr"]["ci95"][0],
        "aupr_hi": res["aupr"]["ci95"][1],
        "aupr_baseline": res["aupr_baseline"],  # --- SANITY CHECKS ---
        "brier": res["brier"]["point"],
        "brier_lo": res["brier"]["ci95"][0],
        "brier_hi": res["brier"]["ci95"][1],
    })

df = pd.DataFrame(rows).sort_values("fold")

# size-weighted overall POINT estimates (not CIs)
w = df["n_test"] / df["n_test"].sum()
overall = {
    "auroc_weighted_mean": float((w * df["auroc"]).sum()),
    "aupr_weighted_mean":  float((w * df["aupr"]).sum()),
    "brier_weighted_mean": float((w * df["brier"]).sum()),
}

# Pretty print
cols = [
    "fold","n_test","n_pos","n_neg",
    "auroc","auroc_lo","auroc_hi",
    "aupr","aupr_lo","aupr_hi","aupr_baseline",
    "brier","brier_lo","brier_hi"
]
print("\n[Per-fold metrics with BCa 95% CIs and sanity checks]")
print(df[cols])

print("\n[Overall size-weighted point estimates]")
print(overall)
