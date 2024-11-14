import os
from copy import copy
import pandas as pd
import numpy as np
from statsmodels.formula.api import logit
from scipy.optimize import linear_sum_assignment
import matplotlib.pyplot as plt
import seaborn as sns

import warnings
warnings.filterwarnings("ignore")


def calc_smd(df, variable, type='t1dm'):
    treated = df[df[type] == 1]
    control = df[df[type] == 0]
    mean_treated = treated[variable].mean()
    mean_control = control[variable].mean()
    std = df[variable].std()
    smd = np.abs(mean_treated - mean_control) / std
    return smd


def psm_norep_2(df, type='t1dm', covars=None):
    """
    :param df: the dataframe to be matched (including treatment variable, covariates)
    :param type: the treatment variable to be matched (or disease status)
    :param covars: the covariates to be used for matching
    :return: the original dataframe and the matched dataframe with control and treatment groups balanced
    """
    if covars is None:
        covars = ['Sex', 'Age', 'Smoking_status', 'Drinking_status', 'Education', 'Income']
    df = df.copy()
    # calculate propensity score by logistic regression for binary treatment
    if covars is None or len(covars) == 0:
        print('No covariates specified, cannot calculate propensity score')
        return None
    else:
        covars_formula = ' + '.join(covars)
    model = logit(f"{type} ~ {covars_formula}", df).fit()
    df['propensity_score'] = model.predict(df)
    # match samples
    df_1 = df[df[type] == 1].copy()
    df_0 = df[df[type] == 0].copy()
    score_matrix = abs(df_1['propensity_score'].values[:,None]-df_0['propensity_score'].values)
    row_ind,col_ind = linear_sum_assignment(score_matrix)
    return df, pd.concat([df_1.iloc[row_ind],df_0.iloc[col_ind]],axis=0)


#
# sns.histplot(data=df, x='propensity_score', hue='t1dm', element='step', stat='density', common_norm=False)
# plt.title('Propensity Score Distribution by Treatment Group')
# plt.show()
#
# sns.histplot(data=matched_df, x='propensity_score', hue='t1dm', element='step', stat='density', common_norm=False)
# plt.title('Propensity Score Distribution by Treatment Group')
# plt.show()
#
#
# for var in ['sex','alcohol','smoking','age','BMI','edu_score','income_score']:
#     smd = calc_smd(df, var, 't1dm')
#     print(f'Standardized Mean Difference for {var}: {smd}')
# print()
# for var in ['sex','alcohol','smoking','age','BMI','edu_score','income_score']:
#     smd = calc_smd(matched_df, var, 't1dm')
#     print(f'Standardized Mean Difference for {var}: {smd}')