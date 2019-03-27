import pandas as pd

# Proteomics data cleaning

proteomics = pd.read_csv('rm_Clean3.csv', header=0, index_col=0)

proteomics = proteomics.dropna(axis=0)

proteomics['MM692'] = (proteomics['MM692']+proteomics['MM692.1'])/2
proteomics['MM778'] = (proteomics['MM778']+proteomics['MM778.1'])/2
proteomics['MM807'] = (proteomics['MM807']+proteomics['MM807.1'])/2
proteomics['MM790'] = (proteomics['MM790']*77.76+proteomics['MM790.1']*78.1923076923077)/(77.76+78.1923076923077)*2
proteomics['MM814'] = proteomics['LG7438']

proteomics = proteomics.drop(['accession', 'MM807.1', 'MM790.1', 'LG7438', 'MM692.1', 'MM778.1'], axis=1)

lst = list(proteomics.columns.values)

# Clinical data cleaning
raw_clinical = pd.read_excel("Copy of ClinicalData_144samplesLundMM_Query2018-10-01 (2).xlsx", header=None)

raw_clinical = raw_clinical.drop(list(range(13, 24)), axis=1)
raw_clinical = raw_clinical.drop(31, axis=1)
raw_clinical = raw_clinical.rename(columns=raw_clinical.iloc[0])

clinical = raw_clinical.drop(0, axis=0)

clinical['sample'] = clinical['sample'].str.upper()

clinical.index = clinical['sample']

cls = list(clinical.index.values)

clinical = clinical.drop('sample', axis=1)
clinical = clinical.drop(index='MM1069', axis=0)

wona_clinical = clinical.dropna(0)
wona_clinical = pd.get_dummies(wona_clinical, columns=['stage', 'Alive 2016-12-05', 'local', 'clin.class.det',
                                                       'BRAF status', 'gender',
                                                       'Alive after 5 years from prim diagn (pos sign on days in Col T)',
                                                       'Alive after 5 years from first met (pos sign on days in Col V)'])
wona_clinical = wona_clinical.drop(['Alive after 5 years from prim diagn (pos sign on days in Col T)_False',
                                    'Alive after 5 years from first met (pos sign on days in Col V)_False'], axis=1)

wona_clinical = wona_clinical[wona_clinical.index.isin(lst)]
wona_clinical.to_csv('wona_clinical.tsv', sep='\t', index=True, header=True)


wna_clinical = pd.get_dummies(clinical, columns=['stage', 'Alive 2016-12-05', 'local', 'clin.class.det',
                                                       'BRAF status', 'gender',
                                                       'Alive after 5 years from prim diagn (pos sign on days in Col T)',
                                                       'Alive after 5 years from first met (pos sign on days in Col V)'],
                              dummy_na=True)
wna_clinical = wna_clinical.drop(['Alive after 5 years from prim diagn (pos sign on days in Col T)_False',
                                    'Alive after 5 years from first met (pos sign on days in Col V)_False'], axis=1)

wna_clinical = wna_clinical[wna_clinical.index.isin(lst)]

wna_clinical.to_csv('wna_clinical.tsv', sep='\t', index=True, header=True)


wona_clinical = wona_clinical[wona_clinical.index.isin(lst)]

proteomics = proteomics[proteomics.columns[proteomics.columns.isin(cls)]]

proteomics.to_csv('Clean_proteomics.tsv', sep='\t', index=True, header=True)
