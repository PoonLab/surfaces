import pandas as pd

rcor_long = pd.read_csv ("/Users/sareh/Desktop/pearson_correlation_clean.csv")
rcor_wide = pd.pivot_table(rcor_long,index="gene1",columns='gene2', values='pearson_correlation',aggfunc=max)

print(rcor_long)
print(rcor_wide)

rcor_wide.to_csv("/Users/sareh/Desktop/rcor_wide.csv")
