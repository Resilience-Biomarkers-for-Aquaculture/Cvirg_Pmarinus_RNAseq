import pandas as pd, numpy as np
ints = pd.read_csv("results/loso_intercepts_by_fold.csv")
ints["p_tolerant_at_zero"] = 1/(1+np.exp(-ints["intercept"]))
print(ints)
