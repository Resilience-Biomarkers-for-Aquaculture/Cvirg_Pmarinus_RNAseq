# Compute per-fold calibration slope/intercept from preds files
from pathlib import Path
import pandas as pd
import numpy as np
from sklearn.linear_model import LogisticRegression

rows=[]
for p in Path("results").glob("loso_*/preds_*.csv"):
    d = p.parent.name.replace("loso_","")
    df = pd.read_csv(p)
    y = df["true_label"].values.astype(int)
    p1 = df["pred_prob_tolerant"].clip(1e-6, 1-1e-6).values
    # Fit calibration model: logit(p) ~ a + b * score
    X = np.log(p1/(1-p1)).reshape(-1,1)
    # Use simple linear regression via closed form
    X1 = np.c_[np.ones_like(X), X]
    beta = np.linalg.lstsq(X1, np.log((y+0.5)/(1.5-y)), rcond=None)[0]  # Firth-ish smoothing
    rows.append({"batch": d, "cal_intercept": beta[0], "cal_slope": beta[1]})
print(pd.DataFrame(rows))
