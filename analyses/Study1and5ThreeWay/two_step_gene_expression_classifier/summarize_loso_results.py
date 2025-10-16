import json, glob, pandas as pd, pathlib

records = []
for path in glob.glob("results/loso_*/best_panel_summary*.json"):
    batch = pathlib.Path(path).parent.name.replace("loso_", "")
    with open(path) as f:
        rec = json.load(f)
    rec["batch"] = batch
    records.append(rec)
df = pd.DataFrame(records)
print(df[["batch","n_test","auroc","aupr","brier"]])

print("\nMean metrics across folds:")
print(df[["auroc","aupr","brier"]].mean())

# Optional bar plot
try:
    import matplotlib.pyplot as plt
    df.plot(x="batch", y=["auroc","aupr"], kind="bar", ylim=(0,1))
    plt.title("LOSO performance per held-out batch")
    plt.ylabel("Score")
    plt.tight_layout()
    plt.show()
except Exception as e:
    print("Plot skipped:", e)
