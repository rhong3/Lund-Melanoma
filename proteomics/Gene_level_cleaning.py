import pandas as pd

protein = pd.read_csv('./Clean_proteomics.tsv', sep='\t', header=0, index_col=0)

disc = pd.read_excel('./Segundo TMT data to Krzysztof.xlsx', header=0).iloc[9:, -3:]

disc = disc.set_index('Accession')

J_protein = protein.join(disc, how='inner')

cols = list(J_protein.columns)
cols = [cols[-1]] + cols[:-1]
cols = [cols[-1]] + cols[:-1]
J_centroid = J_protein[cols]

J_centroid = J_centroid.set_index('Gene Symbol')

J_centroid = J_centroid.drop('Description', axis=1)

J_centroid = J_centroid.groupby(J_centroid.index).mean()

J_centroid.to_csv('./J_Clean_proteomics.tsv', sep='\t', index=True, header=True)

