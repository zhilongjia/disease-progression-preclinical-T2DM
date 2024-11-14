import os
import pandas as pd
import numpy as np
import pySuStaIn
from utils import zscore_regressed_out_covariates, biomarker_selection, normalize_by_controls
import argparse


def main(args):
    # top biomarkers
    # num_biomarkers = args.num_biomarkers

    thres = args.thres
    # dm type, pre or post t2dm
    dm_type = args.dm_type

    df_data = pd.read_csv(f'data/data_t2dm_{dm_type}_balanced_no_outliers.csv')

    biomarkers, num_biomarkers, direction = biomarker_selection(threshold=thres, dm_type=dm_type)

    if len(biomarkers) == 0:
        print(f"No biomarkers selected for threshold {thres}")
        return

    # Parameters
    N_startpoints = args.N_startpoints
    N_S_max = args.N_S_max
    N_iterations_MCMC = args.N_iterations_MCMC
    num_zvals = args.num_zvals

    print(f"Sustain modeling for type 2 diabetes")
    print(f"Number of t2dm: {len(df_data[df_data['t2dm'] == 1])}, "
          f"Number of controls: {len(df_data[df_data['t2dm'] == 0])}")

    # create an output path
    output_path = os.path.join(f'results', f't2dm_{dm_type}_thres{thres}')
    os.makedirs(output_path, exist_ok=True)

    SuStaInLabels = biomarkers

    print(f'{len(SuStaInLabels)} biomarkers for modeling sustain')

    # regress out covariates
    covariates = ['Sex', 'Age', 'Drinking_status', 'Smoking_status', 'Education', 'Income']

    # zscore for IDPs and biomarkers
    df_data_zscored = zscore_regressed_out_covariates(df_data, 't2dm', biomarkers, covariates)
    # df_data_zscored = normalize_by_controls(df_data, 't2dm', biomarkers)

    # min, max, mean, std for each biomarker
    for biomarker in SuStaInLabels:
        print(f"{biomarker} Min: {df_data_zscored[biomarker].min()}, Max: {df_data_zscored[biomarker].max()}, "
              f"Mean: {df_data_zscored[biomarker].mean()}, Std: {df_data_zscored[biomarker].std()}")

    # select data
    data = df_data_zscored[SuStaInLabels].values

    # apply direction
    data = data * np.array(direction)
    # data_controls = df_data_zscored[df_data_zscored['t2dm'] == 0][SuStaInLabels].values

    # setup dataset name
    dataset_name = f't2dm_{dm_type}'

    # setup z_vals and z_max
    z_max = np.quantile(data, 0.95, axis=0)

    # z_vals for each SuStaInLabels from 0 to z_max, num_zvals breakpoints
    z_vals = []
    for i in range(len(SuStaInLabels)):
        z_vals.append(np.linspace(0, z_max[i], num_zvals + 2)[1:num_zvals + 1])
    z_vals = np.array(z_vals)

    print('use_parallel_startpoints:', args.parallel_startpoints)

    print(z_vals.shape)
    # Initiate the SuStaIn object
    sustain_input = pySuStaIn.ZscoreSustain(
        data,
        z_vals,
        z_max,
        SuStaInLabels,
        N_startpoints,
        N_S_max,
        N_iterations_MCMC,
        output_path,
        dataset_name,
        use_parallel_startpoints=args.parallel_startpoints
    )

    # run pySuStaIn
    sustain_input.run_sustain_algorithm()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Sustain modeling for type 2 diabetes')
    parser.add_argument('--dm_type', type=str, default='pre', help='Type of diabetes, pre or post, default pre')
    parser.add_argument('--thres', type=float, default=0.3, help='Threshold for biomarker selection, default 0.3')
    parser.add_argument('--num_biomarkers', type=int, default=10, help='Number of biomarkers to select')
    parser.add_argument('--N_startpoints', type=int, default=25, help='Number of startpoints')
    parser.add_argument('--N_S_max', type=int, default=5, help='Maximum number of subtypes')
    parser.add_argument('--N_iterations_MCMC', type=int, default=int(1e5),
                        help='Number of iterations for MCMC')
    # num_zvals for each biomarker
    parser.add_argument('--num_zvals', type=int, default=2, help='Number of z values for each biomarker')
    parser.add_argument('--parallel_startpoints', action="store_false",
                        help='Whether or not to parallelize the maximum likelihood loop')
    args = parser.parse_args()
    main(args)
