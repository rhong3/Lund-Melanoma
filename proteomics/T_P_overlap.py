import pandas as pd

RNA = pd.read_csv('./Transcriptome/Clean_transcriptome.tsv', sep='\t', header=0, index_col=0)

protein = pd.read_csv('Clean_proteomics.tsv', sep='\t', header=0, index_col=0)

Rls = list(RNA.columns.values)

Pls = list(protein.columns.values)

PRls = list(set(Rls) & set(Pls))

R_clinical = pd.read_csv('./Transcriptome/T_wna_clinical.tsv' ,sep='\t', header=0, index_col=0)

P_clinical = pd.read_csv('wna_clinical.tsv' ,sep='\t', header=0, index_col=0)

RCls = list(R_clinical.index.values)

PCls = list(P_clinical.index.values)

PRCls = list(set(RCls) & set(PCls))

noin = protein.columns[~protein.columns.isin(Rls)]

print(noin)

