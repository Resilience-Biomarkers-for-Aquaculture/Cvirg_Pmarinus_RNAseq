import pandas as pd
df = pd.read_csv('meta_ranked.csv')
col = df['model_used']
n = len(df)
for top in [10,25,50,100,250,500,1000,2500,5000,10000,n]:
    sub = col.iloc[:top]
    fe = (sub=='FE').sum()
    re = (sub=='RE').sum()
    print(f"Top {top:6d}: FE={fe:6d} ({fe/top*100:5.1f}%)  RE={re:6d} ({re/top*100:5.1f}%)")
