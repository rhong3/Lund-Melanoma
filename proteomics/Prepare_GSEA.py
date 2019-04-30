import pandas as pd

ica = pd.read_csv("./ICA/Gene_level/J_Clean_proteomics_IC_centroid.txt", sep='\t', header=0, index_col=0)

wna_clinical = pd.read_csv('./wna_clinical.tsv', sep='\t', header=0, index_col=0)

proteomics = pd.read_csv('./J_Clean_proteomics.tsv', sep='\t', header=0, index_col=0)

ica = ica.sort_values(by=['77'])

wna_clinical = wna_clinical.sort_values(by=['local_Cutaneous', 'local_Lymph node', 'local_Subcutaneous', 'local_Visceral'])

ind = list(wna_clinical.index)

col = list(ica.index)

proteomics = proteomics[ind]

proteomics = proteomics.reindex(col)

print(proteomics)

proteomics.to_csv('./GSEA_proteomics.csv', index=True, header=True)

