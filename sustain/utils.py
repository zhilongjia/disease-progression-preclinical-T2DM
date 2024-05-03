import warnings

import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
from matplotlib import pyplot as plt
from pySuStaIn import AbstractSustain


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


def biomarker_selection(threshold=0.1, num_biomarkers=None, dm_type='pre'):
    """
    Select top biomarkers for modeling
    :param threshold: threshold of log2 fold change for biomarker selection, default is 0.1
    :param num_biomarkers: number of biomarkers to select, must be specified if threshold is not specified
    """
    # df_bio_tops = pd.read_excel('data/data.xlsx')
    df_bio_tops = pd.read_csv(f'data/logistic/biomarkers_{dm_type}_logreg_no_outlier.csv')

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


def plot_positional_var(samples_sequence, samples_f, n_samples, Z_vals, biomarker_labels=None,
                        ml_f_EM=None, cval=False, subtype_order=None, biomarker_order=None,
                        title_font_size=12, stage_font_size=10, stage_label='SuStaIn Stage',
                        stage_rot=0, stage_interval=1, label_font_size=10, label_rot=0, cmap="original",
                        biomarker_colours=None, figsize=None, subtype_titles=None, separate_subtypes=False,
                        save_path=None, save_kwargs={}):
    # Get the number of subtypes
    N_S = samples_sequence.shape[0]
    # Get the number of features/biomarkers
    N_bio = Z_vals.shape[0]
    # Check that the number of labels given match
    if biomarker_labels is not None:
        assert len(biomarker_labels) == N_bio
    # Set subtype order if not given
    if subtype_order is None:
        # Determine order if info given
        if ml_f_EM is not None:
            subtype_order = np.argsort(ml_f_EM)[::-1]
        # Otherwise determine order from samples_f
        else:
            subtype_order = np.argsort(np.mean(samples_f, 1))[::-1]
    elif isinstance(subtype_order, tuple):
        subtype_order = list(subtype_order)
    # Unravel the stage zscores from Z_vals
    stage_zscore = Z_vals.T.flatten()
    IX_select = np.nonzero(stage_zscore)[0]
    stage_zscore = stage_zscore[IX_select][None, :]
    # Get the z-scores and their number
    zvalues = np.unique(stage_zscore)
    N_z = len(zvalues)
    # Extract which biomarkers have which zscores/stages
    stage_biomarker_index = np.tile(np.arange(N_bio), (N_z,))
    stage_biomarker_index = stage_biomarker_index[IX_select]
    # Warn user of reordering if labels and order given
    if biomarker_labels is not None and biomarker_order is not None:
        warnings.warn(
            "Both labels and an order have been given. The labels will be reordered according to the given order!"
        )
    if biomarker_order is not None:
        # self._plot_biomarker_order is not suited to zscore version
        # Ignore for compatability, for now
        # One option is to reshape, sum position, and lowest->highest determines order
        if len(biomarker_order) > N_bio:
            biomarker_order = np.arange(N_bio)
    # Otherwise use default order
    else:
        biomarker_order = np.arange(N_bio)
    # If no labels given, set dummy defaults
    if biomarker_labels is None:
        biomarker_labels = [f"Biomarker {i}" for i in range(N_bio)]
    # Otherwise reorder according to given order (or not if not given)
    else:
        biomarker_labels = [biomarker_labels[i] for i in biomarker_order]
    # Check number of subtype titles is correct if given
    if subtype_titles is not None:
        assert len(subtype_titles) == N_S
    # Z-score colour definition
    if cmap == "original":
        # Hard-coded colours: hooray!
        colour_mat = np.array([[1, 0, 0], [0, 0, 1], [1, 0, 1], [0.5, 0, 1], [0, 1, 1], [0, 1, 0.5]])[:N_z]
        # We only have up to 5 default colours, so double-check
        if colour_mat.shape[0] > N_z:
            raise ValueError(f"Colours are only defined for {len(colour_mat)} z-scores!")
    else:
        raise NotImplementedError
    '''
    Note for future self/others: The use of any arbitrary colourmap is problematic, as when the same stage can have the same biomarker with different z-scores of different certainties, the colours need to mix in a visually informative way and there can be issues with RGB mixing/interpolation, particulary if there are >2 z-scores for the same biomarker at the same stage. It may be possible, but the end result may no longer be useful to look at.
    '''

    # Check biomarker label colours
    # If custom biomarker text colours are given
    if biomarker_colours is not None:
        biomarker_colours = AbstractSustain.check_biomarker_colours(
        biomarker_colours, biomarker_labels
    )
    # Default case of all-black colours
    # Unnecessary, but skips a check later
    else:
        biomarker_colours = {i:"black" for i in biomarker_labels}

    # Flag to plot subtypes separately
    if separate_subtypes:
        nrows, ncols = 1, 1
    else:
        # Determine number of rows and columns (rounded up)
        if N_S == 1:
            nrows, ncols = 1, 1
        elif N_S < 3:
            nrows, ncols = 1, N_S
        elif N_S < 7:
            nrows, ncols = 2, int(np.ceil(N_S / 2))
        else:
            nrows, ncols = 3, int(np.ceil(N_S / 3))
    # Total axes used to loop over
    total_axes = nrows * ncols
    # Create list of single figure object if not separated
    if separate_subtypes:
        subtype_loops = N_S
    else:
        subtype_loops = 1
    # Container for all figure objects
    figs = []
    # Loop over figures (only makes a diff if separate_subtypes=True)
    for i in range(subtype_loops):
        # Create the figure and axis for this subtype loop
        fig, axs = plt.subplots(nrows, ncols, figsize=figsize)
        figs.append(fig)
        # Loop over each axis
        for j in range(total_axes):
            # Normal functionality (all subtypes on one plot)
            if not separate_subtypes:
                i = j
            # Handle case of a single array
            if isinstance(axs, np.ndarray):
                ax = axs.flat[i]
            else:
                ax = axs
            # Check if i is superfluous
            if i not in range(N_S):
                ax.set_axis_off()
                continue

            this_samples_sequence = samples_sequence[subtype_order[i],:,:].T
            N = this_samples_sequence.shape[1]

            # Construct confusion matrix (vectorized)
            # We compare `this_samples_sequence` against each position
            # Sum each time it was observed at that point in the sequence
            # And normalize for number of samples/sequences
            confus_matrix = (this_samples_sequence==np.arange(N)[:, None, None]).sum(1) / this_samples_sequence.shape[0]

            # Define the confusion matrix to insert the colours
            # Use 1s to start with all white
            confus_matrix_c = np.ones((N_bio, N, 3))

            # Loop over each z-score event
            for j, z in enumerate(zvalues):
                # Determine which colours to alter
                # I.e. red (1,0,0) means removing green & blue channels
                # according to the certainty of red (representing z-score 1)
                alter_level = colour_mat[j] == 0
                # Extract the uncertainties for this z-score
                confus_matrix_zscore = confus_matrix[(stage_zscore==z)[0]]
                # Subtract the certainty for this colour
                confus_matrix_c[
                    np.ix_(
                        stage_biomarker_index[(stage_zscore==z)[0]], range(N),
                        alter_level
                    )
                ] -= np.tile(
                    confus_matrix_zscore.reshape((stage_zscore==z).sum(), N, 1),
                    (1, 1, alter_level.sum())
                )
            if subtype_titles is not None:
                title_i = subtype_titles[i]
            else:
                # Add axis title
                if cval == False:
                    temp_mean_f = np.mean(samples_f, 1)
                    # Shuffle vals according to subtype_order
                    # This defaults to previous method if custom order not given
                    vals = temp_mean_f[subtype_order]

                    if n_samples != np.inf:
                        title_i = f"Subtype {i+1} (f={vals[i]:.2f}, n={np.round(vals[i] * n_samples):n})"
                    else:
                        title_i = f"Subtype {i+1} (f={vals[i]:.2f})"
                else:
                    title_i = f"Subtype {i+1} cross-validated"

            # Plot the colourized matrix
            ax.imshow(
                confus_matrix_c[biomarker_order, :, :],
                interpolation='nearest'
            )
            # Add the xticks and labels
            stage_ticks = np.arange(0, N, stage_interval)
            ax.set_xticks(stage_ticks)
            ax.set_xticklabels(stage_ticks+1, fontsize=stage_font_size, rotation=stage_rot)
            # Add the yticks and labels
            ax.set_yticks(np.arange(N_bio))
            # Add biomarker labels to LHS of every row only
            if (i % ncols) == 0:
                ax.set_yticklabels(biomarker_labels, ha='right', fontsize=label_font_size, rotation=label_rot)
                # Set biomarker label colours
                for tick_label in ax.get_yticklabels():
                    tick_label.set_color(biomarker_colours[tick_label.get_text()])
            else:
                ax.set_yticklabels([])
            # Make the event label slightly bigger than the ticks
            ax.set_xlabel(stage_label, fontsize=stage_font_size+2)
            ax.set_title(title_i, fontsize=title_font_size)
        # Tighten up the figure
        fig.tight_layout()
        # Save if a path is given
        if save_path is not None:
            # Modify path for specific subtype if specified
            # Don't modify save_path!
            if separate_subtypes:
                save_name = f"{save_path}_subtype{i}"
            else:
                save_name = f"{save_path}_all-subtypes"
            # Handle file format, avoids issue with . in filenames
            if "format" in save_kwargs:
                file_format = save_kwargs.pop("format")
            # Default to png
            else:
                file_format = "png"
            # Save the figure, with additional kwargs
            fig.savefig(
                f"{save_name}.{file_format}",
                **save_kwargs
            )
    return figs, axs


def calculate_eFGR(Scr, age, sex, race):
    """
    Calculate eFGR using the CKD-EPI formula
    :param Scr: serum creatinine level (mg/dL)
    :param age: age
    :param sex: 0 for female, 1 for male
    :param race: 1 for black, 0 for non-black
    eFGR = 141 * min(Scr / κ, 1)^α * max(Scr / κ, 1)^-1.209 * 0.993^age * 1.018 (if female) * 1.159 (if black)
    where κ = 0.7 for female, 0.9 for male
    and α = -0.329 for female, -0.411 for male
    :return: eFGR
    """

    # define parameters
    kappa = np.where(sex == 0, 0.7, 0.9)
    alpha = np.where(sex == 0, -0.329, -0.411)
    sex_beta = np.where(sex == 0, 1.018, 1)
    beta = np.where(race == 1, 1.159, 1)

    eFGR = 141 * (np.minimum(Scr / kappa, 1) ** alpha) * (np.maximum(Scr / kappa, 1) ** -1.209) * 0.993 ** age
    eFGR *= sex_beta * beta

    return eFGR


if __name__ == '__main__':
    # test calculate_eFGR
    Scr = np.array([0.8, 0.9, 1.0, 0.6, 0.7, 0.8])
    np.random.seed(0)
    age = np.random.randint(50, 100, 6)
    sex = np.random.randint(0, 2, 6)
    race = np.random.randint(0, 2, 6)
    df = pd.DataFrame({'Scr': Scr, 'age': age, 'sex': sex, 'race': race})

    eFGR = calculate_eFGR(df['Scr'], df['age'], df['sex'], df['race'])
    print(eFGR)

