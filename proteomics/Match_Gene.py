import pandas as pd

centroid = pd.read_csv('./ICA/Initial/Clean_proteomics_IC_centroid.txt', sep='\t', header=0, index_col=0)

disc = pd.read_excel('./Segundo TMT data to Krzysztof.xlsx', header=0).iloc[9:, -3:]

disc = disc.set_index('Accession')

J_centroid = centroid.join(disc, how='inner')

cols = list(J_centroid.columns)
cols = [cols[-1]] + cols[:-1]
cols = [cols[-1]] + cols[:-1]
J_centroid = J_centroid[cols]

J_centroid.to_csv('./ICA/Initial/J_Clean_proteomics_IC_centroid.csv', index=True, header=True)

