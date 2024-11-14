import pandas as pd
import statsmodels.formula.api as smf


def zscore_regressed_out_covariates(df_data: pd.DataFrame, group, IDP_id: list, covariates: list):
    # make a copy of our dataframe (we don't want to overwrite our original data)
    df_data_copy = df_data.copy()

    # for each IDP, regress out the effect of covariates
    for idp in IDP_id:
        # create a formula for the regression
        formula = idp + ' ~ ' + ' + '.join(covariates)

        # fit the regression model with controls
        model = smf.ols(formula=formula,
                        data=df_data[df_data[group] == 0]).fit()

        # get the "predicted" values for all subjects based on the control model parameters
        predicted = model.predict(df_data[covariates + [idp]])

        # calculate our zscore: observed - predicted / SD of the control group residuals
        w_score = (df_data.loc[:, idp] - predicted) / model.resid.std()

        # save zscore back into our new (copied) dataframe
        df_data_copy.loc[:, idp] = w_score

    return df_data_copy


def normalize_by_controls(df_data: pd.DataFrame, inflammatory_disease: str, biomarkers: list):

    # control data
    df_control = df_data[df_data[inflammatory_disease] == 0]

    # copy data
    df_data_copy = df_data.copy()

    # for each IDP, calculate the mean and standard deviation of the control group
    for idp in biomarkers:
        df_control_mean = df_control[idp].mean()
        df_control_std = df_control[idp].std()
        # normalize data by controls
        df_data_copy[idp] = (df_data[idp] - df_control_mean) / df_control_std

    return df_data_copy


def biomarker_selection(threshold=0.1, num_biomarkers=None):
    """
    Select top biomarkers for modeling
    :param threshold: threshold of log2 fold change for biomarker selection, default is 0.1
    :param num_biomarkers: number of biomarkers to select, must be specified if threshold is not specified
    """
    # df_bio_tops = pd.read_excel('data/data.xlsx')
    df_bio_tops = pd.read_csv('data/biomarkers_log2fc_no_outlier.csv')

    if threshold is not None:
        biomarkers = df_bio_tops[df_bio_tops['abs_log2fc'] > threshold]['biomarker'].tolist()
        direction = df_bio_tops[df_bio_tops['abs_log2fc'] > threshold]['direction'].tolist()
        num_biomarkers = len(biomarkers)
    else:
        assert num_biomarkers is not None, "Please specify the number of biomarkers to select"
        assert 0 < num_biomarkers < len(df_bio_tops), ("Number of biomarkers to select "
                                                       "should be within the range of biomarkers")

        biomarkers = df_bio_tops.sort_values('abs_log2fc', ascending=False).head(num_biomarkers)['biomarker'].tolist()
        direction = df_bio_tops.sort_values('abs_log2fc', ascending=False).head(num_biomarkers)['direction'].tolist()

    print(f"Selected {len(biomarkers)} biomarkers for modeling:")
    # print with no quotes
    print(*biomarkers, sep=', ')

    return biomarkers, num_biomarkers, direction
