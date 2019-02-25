import pandas as pd

ica = pd.read_csv("./ICA/Initial/Clean_proteomics_IC_centroid.txt", sep='\t', header=0, index_col=0)

wna_clinical = pd.read_csv('./wna_clinical.tsv', sep='\t', header=0, index_col=0)

proteomics = pd.read_csv('./Clean_proteomics.tsv', sep='\t', header=0, index_col=0)

ica = ica.sort_values(by=['62'], ascending=False)

wna_clinical = wna_clinical.sort_values(by=['5 year survival from primary diagnosis ((Date death-date prim diagn)-1825 days). Days differing from 5 years',
                                            '5 year survival from first metastasis ((Date death-date first met)-1825 days). Days differing from 5 years'])

ind = list(wna_clinical.index)

col = list(ica.index)

proteomics = proteomics[ind]

proteomics = proteomics.reindex(col)


print(proteomics)

proteomics.to_csv('./GSEA_proteomics.csv', index=True, header=True)

