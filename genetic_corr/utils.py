import os
import pandas as pd
import numpy as np
# fdr correction
from statsmodels.stats.multitest import multipletests


def read_ldsc_res(subtypes, phenos, res_path):
    results = []
    header = None
    for subtype in subtypes:
        for i, p in enumerate(phenos):
            print('reading genetic correlation results', subtype, p)
            # file path
            file_path = os.path.join(res_path, f'{subtype}_{p}_h2.log')
            with open(file_path, 'r') as f:
                lines = f.readlines()
                if i == 0:
                    # line 61 as header
                    headline = lines[-5]
                    # remove '\n'
                    headline = headline.strip()
                # read line 62
                resline = lines[-4]
                # remove '\n'
                resline = resline.strip()

            # split
            header = headline.split()
            res = resline.split()
            # repalce 'NA' with np.nan
            res = [np.nan if x == 'NA' else x for x in res]
            results.append(res)

    # convert to dataframe, set column type: first two columns as string, others as float
    df = pd.DataFrame(results, columns=header)
    # set first two columns as string, others as float
    df.iloc[:, 0:2] = df.iloc[:, 0:2].astype(str)
    df.iloc[:, 2:] = df.iloc[:, 2:].astype(float)
    return df
