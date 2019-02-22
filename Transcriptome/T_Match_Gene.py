import pandas as pd

centroid = pd.read_csv('./ICA/Initial/Clean_transcriptome_IC_centroid.txt', sep='\t', header=0, index_col=0)

disc = pd.read_excel('./Transcriptomics_final.214.allMm.not.meancentered_rounded.xlsx', header=0, index_col=0).iloc[:,:7]

J_centroid = centroid.join(disc, how='inner')

cols = list(J_centroid.columns)
for i in range(7):
    cols = [cols[-1]] + cols[:-1]
J_centroid = J_centroid[cols]

J_centroid.to_csv('./ICA/Initial/J_Clean_transcriptome_IC_centroid.csv', index=True, header=True)

