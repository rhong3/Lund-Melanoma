import pandas as pd

RNA = pd.read_csv('Transcriptome/Clean_transcriptome.tsv', sep='\t', header=0, index_col=0)

disc = pd.read_excel('Transcriptome/Transcriptomics_final.214.allMm.not.meancentered_rounded.xlsx', header=0, index_col=0).iloc[:,:2]

J_centroid = RNA.join(disc, how='inner')

cols = list(J_centroid.columns)
for i in range(2):
    cols = [cols[-1]] + cols[:-1]
J_centroid = J_centroid[cols]

J_centroid = J_centroid.set_index('SYMBOL')

J_centroid = J_centroid.drop('PROBE_ID', axis=1)

J_centroid = J_centroid.groupby(J_centroid.index).mean()

J_centroid.to_csv('Transcriptome/J_Clean_transcriptome.tsv', sep='\t', index=True, header=True)

