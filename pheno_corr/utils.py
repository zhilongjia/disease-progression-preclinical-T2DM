import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
from lifelines import CoxPHFitter
from lifelines.utils import concordance_index
from sklearn.model_selection import StratifiedKFold


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
    df_bio_tops = pd.read_csv('data/biomarkers_logreg_no_outlier.csv')

    if threshold is not None:
        biomarkers = df_bio_tops[df_bio_tops['abs_Beta'] > threshold]['Biomarker'].tolist()
        direction = df_bio_tops[df_bio_tops['abs_Beta'] > threshold]['direction'].tolist()
        num_biomarkers = len(biomarkers)
    else:
        assert num_biomarkers is not None, "Please specify the number of biomarkers to select"
        assert 0 < num_biomarkers < len(df_bio_tops), ("Number of biomarkers to select "
                                                       "should be within the range of biomarkers")

        biomarkers = df_bio_tops.sort_values('abs_Beta', ascending=False).head(num_biomarkers)['Biomarker'].tolist()
        direction = df_bio_tops.sort_values('abs_Beta', ascending=False).head(num_biomarkers)['direction'].tolist()

    print(f"Selected {len(biomarkers)} biomarkers for modeling:")
    # print with no quotes
    print(*biomarkers, sep=', ')

    return biomarkers, num_biomarkers, direction


# define a function that cut an array into x bins, and make sure the sum of the cut bins are the most equally distributed
# example: [1, 6, 5, 1, 1, 3, 3, 1], 3 bins, should be cut into [1, 6], [5, 1, 1], [3, 3, 1],
# so that the sum of each bin is the most equally distributed
# return the cut position of the array, in this case, [1, 4, 5]
def find_min_max_partition(arr, x):
    n = len(arr)
    # 前缀和，用于快速计算区间和
    prefix_sums = [0] * (n + 1)
    for i in range(n):
        prefix_sums[i + 1] = prefix_sums[i] + arr[i]

    # dp[i][j] 表示前 i 个元素分成 j 部分的最小的最大分组和
    dp = [[float('inf')] * (x + 1) for _ in range(n + 1)]
    # base case: 前 i 个元素分成 1 组就是它们的和
    for i in range(1, n + 1):
        dp[i][1] = prefix_sums[i]

    # 记录分割点
    partition = [[0] * (x + 1) for _ in range(n + 1)]

    for i in range(1, n + 1):
        for j in range(2, min(i, x) + 1):  # 分成至少两组
            for k in range(j - 1, i):  # 枚举最后一组的起点
                # 最大分组和的最小值
                cost = max(dp[k][j - 1], prefix_sums[i] - prefix_sums[k])
                if cost < dp[i][j]:
                    dp[i][j] = cost
                    partition[i][j] = k  # 记录分割点

    # 从 partition 回溯找到所有分割点
    splits = []
    current = n
    for i in range(x, 1, -1):
        splits.append(partition[current][i] - 1)
        current = partition[current][i]
    splits.reverse()

    # 重建解决方案
    parts = []
    last_idx = 0
    for split in splits + [n - 1]:
        parts.append(arr[last_idx:split + 1])
        last_idx = split + 1

    # sum of each part
    sum_parts = [sum(part) for part in parts]

    return parts, sum_parts, splits


def SAkfold(diseases, feats, df_bio, kf=5):
    # 5-fold cross-validation
    skf = StratifiedKFold(n_splits=kf, random_state=42, shuffle=True)
    res = np.zeros((len(diseases), kf))
    # use the first N components as features for prediction
    # feats = ['PC' + str(i) for i in range(N)]
    for i, dis in enumerate(diseases):
        df_sa_s1 = pd.read_csv(f'data/subtype1/survival_data_{dis}.csv')
        df_sa_s1['Subtype'] = np.where(df_sa_s1['Subtype'] == 1, 1, 0)
        df_sa_s2 = pd.read_csv(f'data/subtype2/survival_data_{dis}.csv')
        df_sa_s2['Subtype'] = np.where(df_sa_s2['Subtype'] == 1, 2, 0)
        # merge the two subtype data and remove duplicates
        df_sa = pd.concat([df_sa_s1, df_sa_s2], axis=0)
        df_sa = df_sa.drop_duplicates(subset=['Eid'])
        df_bio_x = df_bio[['Eid'] + feats]
        df = pd.merge(df_sa, df_bio_x, on='Eid', how='inner')
        fold = 0
        print('SA for', dis)
        for tr_idx, te_idx in skf.split(df[feats], df[dis]):
            print('fold: ', fold)
            df_tr, df_te = df.iloc[tr_idx], df.iloc[te_idx]
            cph = CoxPHFitter()
            df_tr = df_tr[['time', dis] + feats]
            cph.fit(df_tr, duration_col='time', event_col=dis)
            # predict the risk score
            df_te = df_te[['time', dis] + feats]
            y_pred = cph.predict_expectation(df_te)
            # calculate the concordance index
            c_index = concordance_index(df_te['time'], y_pred, df_te[dis])
            res[i, fold] = c_index
            print(dis, ': c-index on fold ', fold, ': ', c_index)
            fold += 1

    return res


def pred_risk(diseases, feats, df_bio, kf=5, subtype=None):
    skf = StratifiedKFold(n_splits=kf, random_state=42, shuffle=True)
    for i, dis in enumerate(diseases):
        df_sa_s1 = pd.read_csv(f'data/subtype1/survival_data_{dis}.csv')
        df_sa_s1['Subtype'] = np.where(df_sa_s1['Subtype'] == 1, 1, 0)
        df_sa_s2 = pd.read_csv(f'data/subtype2/survival_data_{dis}.csv')
        df_sa_s2['Subtype'] = np.where(df_sa_s2['Subtype'] == 1, 2, 0)
        # merge the two subtype data and remove duplicates
        df_sa = pd.concat([df_sa_s1, df_sa_s2], axis=0)
        df_sa = df_sa.drop_duplicates(subset=['Eid'])
        df_bio_x = df_bio[['Eid'] + feats]
        df = pd.merge(df_sa, df_bio_x, on='Eid', how='inner')
        fold = 0
        print('SA for', dis)
        pred_risks = []
        for tr_idx, te_idx in skf.split(df[feats], df[dis]):
            print('fold: ', fold)
            df_tr, df_te = df.iloc[tr_idx], df.iloc[te_idx]
            cph = CoxPHFitter()
            df_tr = df_tr[['time', dis] + feats]
            cph.fit(df_tr, duration_col='time', event_col=dis)
            # predict the risk score
            # df_te = df_te[['time', dis] + feats]
            # pred survival function at time 1-15 years
            y_pred = cph.predict_survival_function(df_te, times=np.arange(1, 16))
            y_pred = y_pred.T
            # to dataframe
            df_pred = pd.DataFrame(y_pred)
            df_pred['Eid'] = df_te['Eid'].values
            df_pred.columns = ['pred_risk_' + str(i) for i in range(1, 16)] + ['Eid']
            pred_risks.append(df_pred)
            fold += 1

        # save the predicted risk scores
        df_pred = pd.concat(pred_risks, axis=0)
        if subtype is not None:
            df_pred.to_csv(f'results/prediction/pred_risk/pred_risk_{dis}_{subtype}.csv', index=False)
        else:
            df_pred.to_csv(f'results/prediction/pred_risk/pred_risk_{dis}_baseline.csv', index=False)


def get_tf_positive(time, label, threshold, pred_prob, eval_time):
    tp = fp = 0
    num_eval = 0
    N = len(pred_prob)
    for i in range(N):
        # positive(died) at eval_time
        if time[i] <= eval_time and label[i] == 1:
            num_eval += 1
            if pred_prob[i] >= threshold:
                tp += 1
        # non-positive(survival) at eval_time
        elif time[i] > eval_time:
            num_eval += 1
            if pred_prob[i] >= threshold:
                fp += 1
    return tp, fp, num_eval


def cal_net_benefit(time, label, threshold, pred_prob, eval_time):
    if threshold == 0:
        threshold = 1e-8
    elif threshold == 1:
        threshold = 1. - 1e-8
    time = np.reshape(time, [-1])
    label = np.reshape(label, [-1])
    tp, fp, num_eval = get_tf_positive(time, label, threshold, pred_prob, eval_time=eval_time)
    if num_eval == 0:
        raise ValueError('Can not calculate net benefit')
    theta = threshold / (1 - threshold)
    res = tp * 1. / num_eval - fp * 1. / num_eval * theta
    return res


def get_km_scores(times, labels, fail_code=1, sort=False):
    """
    # estimate KM survival rate
    :param times: ndarray, shape(num_subject, ), event times or censoring times, shape,
    :param labels: ndarray, shape(num_subject, ), event labels
    :param fail_code: event_id, default 1, 0 for censored, 1 for event
    :param sort: whether sort by times, default False (we assume that the time is sorted in ascending order)
    :return:
    """
    N = len(times)
    times = np.reshape(times, [-1])
    labels = np.reshape(labels, [-1])
    # Sorting T and E in ascending order by T
    if sort:
        order = np.argsort(times)
        T = times[order]
        E = labels[order]
    else:
        T = times
        E = labels
    max_T = int(np.max(T)) + 1

    # calculate KM survival rate at time 0-T_max
    km_scores = np.ones(max_T)
    n_fail = 0
    n_rep = 0

    for i in range(N):

        if E[i] == fail_code:
            n_fail += 1

        if i < N - 1 and T[i] == T[i + 1]:
            n_rep += 1
            continue

        km_scores[int(T[i])] = 1. - n_fail / (N - i + n_rep)
        n_fail = 0
        n_rep = 0

    for i in range(1, max_T):
        km_scores[i] = km_scores[i - 1] * km_scores[i]

    return km_scores


if __name__ == '__main__':
    # test cut_array_to_bins
    arr = [1447, 754, 647, 668, 614, 529, 576, 483, 532, 576, 697, 613, 145, 219, 72, 317, 73, 101, 35, 92, 67, 91, 2, 1, 3, 8, 28, 3, 2, 5, 2]
    num_bins = 5
    index_bins = find_min_max_partition(arr, num_bins)

    print(index_bins)