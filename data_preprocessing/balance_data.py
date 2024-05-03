import pandas as pd
from PSM import calc_smd, psm_norep_2
import numpy as np
from scipy import stats
import argparse


def main(args):
    dm_type = args.dm_type

    df_dm_biomarker = pd.read_csv('data/imputed/data_t2dm_biomarker_imputed.csv')
    # keep subject with genetic data
    # fam file for genetic data
    df_genetic = pd.read_csv('data/allchr.fam', sep=' ', header=None)
    df_genetic.columns = ['FID', 'IID', 'PID', 'MID', 'SEX', 'PHENOTYPE']

    df_dm_biomarker = df_dm_biomarker[df_dm_biomarker['Eid'].isin(df_genetic['IID'])]
    print('number of subjects: ', len(df_dm_biomarker))
    print(df_dm_biomarker['t2dm'].value_counts())

    # remove outliers for each biomarker, values > 5 standard deviation
    biomarkers = df_dm_biomarker.columns[11:].tolist()

    # remove subjects with insulin medication
    df_ins = pd.read_csv('data/medication/insulin.csv')
    eid_used_ins = df_ins[df_ins['insulin'] == 1]['Eid'].tolist()
    df_dm_biomarker = df_dm_biomarker[~df_dm_biomarker['Eid'].isin(eid_used_ins)]
    print('number of subjects after removing insulin medication: ', len(df_dm_biomarker))

    # remove outliers with values > 5 std
    for biomarker in biomarkers:
        mean = df_dm_biomarker[biomarker].mean()
        std = df_dm_biomarker[biomarker].std()
        df_dm_biomarker = df_dm_biomarker[abs(df_dm_biomarker[biomarker] - mean) < 5 * std]
        # print(f'{biomarker}: mean {mean}, std {std}, min {df_dm_biomarker[biomarker].min()}, '
             # f'max {df_dm_biomarker[biomarker].max()}')

    # keep dm_type
    df_time_diff = pd.read_csv('data/time/data_t2dm_time_diff.csv')
    # keep subjects with dm_type and controls
    Eids_dm_type = df_time_diff[df_time_diff['dm_type'].isin(['control', dm_type])]['Eid'].tolist()
    df_dm_biomarker = df_dm_biomarker[df_dm_biomarker['Eid'].isin(Eids_dm_type)]
    print('number of subjects after keeping dm_type: ', len(df_dm_biomarker))

    # remove controls with HbA1c >= 6.5%, glucose >= 7.0 mmol/L, as they are likely to be misclassified
    # HbA1c mmol/mol to %: HbA1c (%) = (HbA1c (mmol/mol) × 0.09148) + 2.152
    df_dm_biomarker['HbA1c_prec'] = (df_dm_biomarker['HbA1c'] * 0.09148) + 2.152
    df_dm_biomarker = df_dm_biomarker[~((df_dm_biomarker['t2dm'] == 0) & (df_dm_biomarker['HbA1c_prec'] >= 6.5))]
    df_dm_biomarker = df_dm_biomarker[~((df_dm_biomarker['t2dm'] == 0) & (df_dm_biomarker['Glucose'] >= 7.0))]

    print('number of subjects after removing outliers and likely diabetes subjects: ', len(df_dm_biomarker))
    print(df_dm_biomarker['t2dm'].value_counts())

    # df_new, df_balanced = psm_norep_2(df_dm_biomarker, 't2dm')
    #
    # df_new.to_csv('data/data_balanced/data_t2dm_all.csv', index=False)
    # df_balanced.to_csv(f'data/data_balanced/data_t2dm_{dm_type}_balanced_no_outliers.csv', index=False)
    #
    # for var in ['Sex', 'Age', 'Smoking_status', 'Drinking_status', 'Education', 'Income']:
    #     smd = calc_smd(df_new, var, 't2dm')
    #     print(f'Standardized Mean Difference for {var}: {smd}')
    #
    # print()
    #
    # for var in ['Sex', 'Age', 'Smoking_status', 'Drinking_status', 'Education', 'Income']:
    #     smd = calc_smd(df_balanced, var, 't2dm')
    #     print(f'Standardized Mean Difference for {var}: {smd}')


if __name__ == '__main__':
    # def args
    argparse = argparse.ArgumentParser()
    argparse.add_argument('--dm_type', type=str, default='pre',
                          help='type of diabetes, pre or post t2dm diagnosis')

    args = argparse.parse_args()
    main(args)