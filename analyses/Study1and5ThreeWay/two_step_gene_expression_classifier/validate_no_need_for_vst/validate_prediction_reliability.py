import pandas as pd
from sklearn.metrics import roc_auc_score, average_precision_score
dfp = pd.read_csv("results/loso_all_predictions.csv")
print("Pooled AUROC:", roc_auc_score(dfp.true_label, dfp.pred_prob_tolerant))
print("Pooled AUPR :", average_precision_score(dfp.true_label, dfp.pred_prob_tolerant))
