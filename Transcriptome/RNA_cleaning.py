import pandas as pd

RNA_raw = pd.read_excel('./Transcriptomics_final.214.allMm.not.meancentered_rounded.xlsx',
                        header=0, index_col=0)

RNA_raw = RNA_raw.drop(['PROBE_ID', 'SYMBOL', 'SEARCH_KEY', 'ILMN_GENE', 'CHROMOSOME',
                        'DEFINITION', 'SYNONYMS', 'undetected'], axis=1)

RNA_raw.columns = map(str.upper, RNA_raw.columns)

lst = list(RNA_raw.columns.values)

wna_clinical = pd.read_csv('../proteomics/wna_clinical.tsv', sep='\t', header=0, index_col=0)
wna_clinical = wna_clinical[wna_clinical.index.isin(lst)]
wona_clinical = pd.read_csv('../proteomics/wona_clinical.tsv', sep='\t', header=0, index_col=0)
wona_clinical = wona_clinical[wona_clinical.index.isin(lst)]

cls = list(wna_clinical.index.values)

RNA_raw = RNA_raw[RNA_raw.columns[RNA_raw.columns.isin(cls)]]

RNA_raw.to_csv('./Clean_transcriptome.tsv', sep='\t', index=True, header=True)
wna_clinical.to_csv('./T_wna_clinical.tsv', sep='\t', index=True, header=True)
wona_clinical.to_csv('./T_wona_clinical.tsv', sep='\t', index=True, header=True)

