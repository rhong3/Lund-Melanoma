import pandas as pd
import numpy as np

wna_clinical = pd.read_csv('./wna_clinical.tsv', sep='\t', header=0, index_col=0)

bins = [-1000, 0, 1000, 2000, 3000]

wna_clinical['Bin 5 year survival from primary diagnosis'] = \
    np.searchsorted(bins,
                    wna_clinical['5 year survival from primary diagnosis ((Date death-date prim diagn)-1825 days). Days differing from 5 years'].values)

wna_clinical['Bin 5 year survival from first metastasis'] = \
    np.searchsorted(bins,
                    wna_clinical['5 year survival from first metastasis ((Date death-date first met)-1825 days). Days differing from 5 years'].values)

wna_clinical = wna_clinical.drop('5 year survival from primary diagnosis ((Date death-date prim diagn)-1825 days). Days differing from 5 years', axis=1)
wna_clinical = wna_clinical.drop('5 year survival from first metastasis ((Date death-date first met)-1825 days). Days differing from 5 years', axis=1)

binsa = [40, 50, 60, 70, 80]

wna_clinical['bin age.diag'] = \
    np.searchsorted(binsa,
                    wna_clinical['age.diag'].values)

wna_clinical['bin prim.age.diag'] = \
    np.searchsorted(binsa,
                    wna_clinical['prim.age.diag'].values)

wna_clinical = wna_clinical.drop('age.diag', axis=1)
wna_clinical = wna_clinical.drop('prim.age.diag', axis=1)

binsb = [1, 2, 3, 4, 5]

wna_clinical['bin prim.breslow'] = \
    np.searchsorted(binsb,
                    wna_clinical['prim.breslow'].values)

wna_clinical = wna_clinical.drop('prim.breslow', axis=1)

binsc = [1000, 2000, 3000, 4000, 5000]

wna_clinical['bin days from primary diagnosis to death or censoring'] = \
    np.searchsorted(binsc,
                    wna_clinical['OS.DAYS     days from primary diagnosis to death or censoring'].values)

wna_clinical = wna_clinical.drop('OS.DAYS     days from primary diagnosis to death or censoring', axis=1)

binsd = [50, 200, 500, 1000, 2000]

wna_clinical['bin days from primary diagnosis to first metastasis'] = \
    np.searchsorted(binsd,
                    wna_clinical['DMFS.DAYS   days from primary diagnosis to first metastasis'].values)

wna_clinical = wna_clinical.drop('DMFS.DAYS   days from primary diagnosis to first metastasis', axis=1)

binse = [250, 500, 1000, 2000, 3500]

wna_clinical['bin days from first metastasis to death or censoring'] = \
    np.searchsorted(binse,
                    wna_clinical['DSS.DAYS.v1 days from first metastasis to death or censoring'].values)

wna_clinical['bin days from sample collection (surgery) to death or censoring'] = \
    np.searchsorted(binse,
                    wna_clinical['DSS.DAYS.v2   days from sample collection (surgery) to death or censoring'].values)

wna_clinical = wna_clinical.drop('DSS.DAYS.v1 days from first metastasis to death or censoring', axis=1)
wna_clinical = wna_clinical.drop('DSS.DAYS.v2   days from sample collection (surgery) to death or censoring', axis=1)

wna_clinical.to_csv('bin_wna_clinical.tsv', sep='\t', index=True, header=True)


