import os
from argparse import ArgumentParser

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import yaml
from dotenv import load_dotenv
from scipy import stats
from scipy.special import comb

# whether to make ribbon or spaghetti plot
ribbon = True

# set mode
allowed_options = ["production", "development"]
parser = ArgumentParser()
parser.add_argument(
    "-om", "--output-mode", dest="output_mode", default="development"
)
args = parser.parse_args()
output_mode = args.output_mode
if output_mode not in allowed_options:
    print(
        "%s not in allowed options, defaulting to development" % (output_mode)
    )
    output_mode = "development"

font = {"family": "Nimbus Roman", "weight": "normal", "size": 22}

mpl.rc("font", **font)

load_dotenv()
output_dir = os.path.expanduser(os.getenv("MAIN_OUTPUT_DIR"))

with open("products_config.yml", "r") as file:
    products_config = yaml.safe_load(file)

with open(os.path.join("analysis", "parameters.yml"), "r") as file:
    parameters_config = yaml.safe_load(file)

output_dir = os.path.join(output_dir, output_mode)

figure_path_out = os.path.join(
    output_dir,
    products_config["products"]["directory"],
)

# Ensure the figure directory exists
os.makedirs(
    figure_path_out,
    exist_ok=True,
)


def triangle_vl_fx(t, dp, tp, wp, wr):
    return (t <= tp) * (dp / wp) * (t - (tp - wp)) + (t > tp) * (
        dp - (dp / wr) * (t - tp)
    )


# taken straight from statsmodels.graphs source code
def banddepth(data, method="MBD"):
    """
    Calculate the band depth for a set of functional curves.

    Band depth is an order statistic for functional data (see `fboxplot`), with
    a higher band depth indicating larger "centrality".  In analog to scalar
    data, the functional curve with highest band depth is called the median
    curve, and the band made up from the first N/2 of N curves is the 50%
    central region.

    Parameters
    ----------
    data : ndarray
        The vectors of functions to create a functional boxplot from.
        The first axis is the function index, the second axis the one along
        which the function is defined.  So ``data[0, :]`` is the first
        functional curve.
    method : {'MBD', 'BD2'}, optional
        Whether to use the original band depth (with J=2) of [1]_ or the
        modified band depth.  See Notes for details.

    Returns
    -------
    ndarray
        Depth values for functional curves.

    Notes
    -----
    Functional band depth as an order statistic for functional data was
    proposed in [1]_ and applied to functional boxplots and bagplots in [2]_.

    The method 'BD2' checks for each curve whether it lies completely inside
    bands constructed from two curves.  All permutations of two curves in the
    set of curves are used, and the band depth is normalized to one.  Due to
    the complete curve having to fall within the band, this method yields a lot
    of ties.

    The method 'MBD' is similar to 'BD2', but checks the fraction of the curve
    falling within the bands.  It therefore generates very few ties.

    The algorithm uses the efficient implementation proposed in [3]_.

    References
    ----------
    .. [1] S. Lopez-Pintado and J. Romo, "On the Concept of Depth for
           Functional Data", Journal of the American Statistical Association,
           vol.  104, pp. 718-734, 2009.
    .. [2] Y. Sun and M.G. Genton, "Functional Boxplots", Journal of
           Computational and Graphical Statistics, vol. 20, pp. 1-19, 2011.
    .. [3] Y. Sun, M. G. Gentonb and D. W. Nychkac, "Exact fast computation
           of band depth for large functional datasets: How quickly can one
           million curves be ranked?", Journal for the Rapid Dissemination
           of Statistics Research, vol. 1, pp. 68-74, 2012.
    """
    n, p = data.shape
    rv = np.argsort(data, axis=0)
    rmat = np.argsort(rv, axis=0) + 1

    # band depth
    def _fbd2():
        down = np.min(rmat, axis=1) - 1
        up = n - np.max(rmat, axis=1)
        return (up * down + n - 1) / comb(n, 2)

    # modified band depth
    def _fmbd():
        down = rmat - 1
        up = n - rmat
        return ((np.sum(up * down, axis=1) / p) + n - 1) / comb(n, 2)

    if method == "BD2":
        depth = _fbd2()
    elif method == "MBD":
        depth = _fmbd()
    else:
        raise ValueError("Unknown input value for parameter `method`.")

    return depth


def fboxplot(data, xdata, ax, color, method="MBD", wfactor=1.5):
    data = np.asarray(data)
    if method not in ["MBD", "BD2"]:
        raise ValueError("Unknown value for parameter `method`.")
    else:
        depth = banddepth(data, method=method)

    # Inner area is 25%-75% region of band-depth ordered curves.
    ix_depth = np.argsort(depth)[::-1]
    median_curve = data[ix_depth[0], :]
    ix_IQR = data.shape[0] // 2
    lower = data[ix_depth[0:ix_IQR], :].min(axis=0)
    upper = data[ix_depth[0:ix_IQR], :].max(axis=0)
    # just thinking of if we instead wanted to get the
    # other quantiles, we want to write this in a generic
    # enough way

    # we want to allow for outlier detection to not be plotted
    # Determine region for outlier detection
    if wfactor is not None:
        inner_median = np.median(data[ix_depth[0:ix_IQR], :], axis=0)
        lower_fence = inner_median - (inner_median - lower) * wfactor
        upper_fence = inner_median + (upper - inner_median) * wfactor

        # find nonoutlier bound
        ix_nonout = []
        for ii in range(data.shape[0]):
            if not (
                np.any(data[ii, :] > upper_fence)
                or np.any(data[ii, :] < lower_fence)
            ):
                ix_nonout.append(ii)

        # Plot envelope of all non-outlying data
        lower_nonout = data[ix_nonout, :].min(axis=0)
        upper_nonout = data[ix_nonout, :].max(axis=0)
        ax.fill_between(
            xdata, lower_nonout, upper_nonout, color=color, alpha=0.1
        )

    # Plot central 50% region
    ax.fill_between(xdata, lower, upper, color=color, alpha=0.3)

    # Plot median curve
    ax.plot(xdata, median_curve, color="k", linewidth=8)


n_extracted_trajectories = 200

symp_categories = parameters_config["global_data_markers"][
    "symptom_categories"
]

max_si = parameters_config["scenario_parameters"]["max_si"]

minimum_day_considered = parameters_config["scenario_parameters"][
    "minimum_days_before_symptom_onset_considered"
]

empirical_values = pd.read_csv(
    os.path.join(
        output_dir,
        products_config["preprocessed_data"]["directory"],
        products_config["vl_model"]["stan_ready_preprocessed_data_output"],
    )
)

empirical_values_viral_load = empirical_values[
    [
        "symp_type_cat",
        "contiguous_id",
        "logVL",
        "days_since_symp_onset",
        "antigen",
        "culture",
    ]
]

empirical_values_symp_improvement = empirical_values[
    [
        "symp_duration",
        "symp_duration_censored",
        "symp_type_cat",
        "contiguous_id",
        "logVL",
        "days_since_symp_onset",
    ]
]

# remove individuals who have censored symp improv times.
empirical_values_symp_improvement = empirical_values_symp_improvement[
    empirical_values_symp_improvement["symp_duration_censored"] == 0
]

empirical_values_symp_improvement = empirical_values_symp_improvement[
    ~empirical_values_symp_improvement.contiguous_id.duplicated(keep="first")
]

empirical_values_symp_improvement_list = [
    empirical_values_symp_improvement.symp_duration[
        empirical_values_symp_improvement.symp_type_cat == i
    ].values
    for i in symp_categories
]

samples = pd.read_csv(
    os.path.join(
        output_dir,
        products_config["vl_model"]["stan_posterior_directory"],
        products_config["vl_model"]["posterior_individual"],
    )
)

tp_samples = [
    samples.tp[samples.symp_type_cat == i].values for i in symp_categories
]
wr_samples = [
    samples.wr[samples.symp_type_cat == i].values for i in symp_categories
]
dp_samples = [
    samples.dp[samples.symp_type_cat == i].values for i in symp_categories
]
wp_samples = [
    samples.wp[samples.symp_type_cat == i].values for i in symp_categories
]

weibull_samples = [
    samples.si_time[samples.symp_type_cat == i].values for i in symp_categories
]

dt = 0.1
ts = np.arange(minimum_day_considered, max_si, dt)

file_out = os.path.join(
    figure_path_out,
    "figure_manuscript_VL_symp_improvement_posteriors_gof",
)

colors = ["#1a1a1a", "#1a1a1a", "#1a1a1a"]

titles = [
    "Shortness of breath",
    "Fever or body aches",
    "Mild respiratory symptoms",
    "Non-specific symptoms",
]

vln_points = 8

fig = plt.figure(figsize=(30, 30 / 3), layout="constrained")
spec = fig.add_gridspec(
    ncols=len(symp_categories) - 1, nrows=2, height_ratios=[3, 1]
)

# don't include symp category 4 in main text!
for i in symp_categories[:-1]:
    # VL trajectory posteriors
    ax1 = fig.add_subplot(spec[0, i - 1])
    ax1.spines["top"].set_visible(False)
    ax1.spines["right"].set_visible(False)
    ax1.get_xaxis().tick_bottom()
    ax1.get_yaxis().tick_left()
    ax1.tick_params(axis="x", direction="out")
    ax1.tick_params(axis="y", direction="out")
    # offset the spines
    for spine in ax1.spines.values():
        spine.set_position(("outward", 5))
    # put the grid behind
    ax1.set_axisbelow(True)

    tp_samples_cat = tp_samples[i - 1]
    wr_samples_cat = wr_samples[i - 1]
    dp_samples_cat = dp_samples[i - 1]
    wp_samples_cat = wp_samples[i - 1]
    empirical_values_viral_load_cat = empirical_values_viral_load[
        empirical_values_viral_load.symp_type_cat == i
    ]

    if not ribbon:  # spaghetti plot
        random_samples = np.random.choice(
            len(dp_samples_cat),
            n_extracted_trajectories,
        )
        for sample in random_samples:
            ax1.plot(
                ts,
                triangle_vl_fx(
                    ts,
                    dp_samples_cat[sample],
                    tp_samples_cat[sample],
                    wp_samples_cat[sample],
                    wr_samples_cat[sample],
                ),
                color=colors[i - 1],
                linestyle="solid",
                alpha=0.2,
            )

        ax1.plot(
            ts * 0,
            0
            * triangle_vl_fx(
                ts,
                dp_samples_cat[sample],
                tp_samples_cat[sample],
                wp_samples_cat[sample],
                wr_samples_cat[sample],
            )
            - 2,  # so it gets cut off
            color="black",
            linestyle="solid",
            alpha=1,
            label="Modeled median logVL\n(posterior samples)",
        )

    else:  # ribbon
        samples_statistics_individual = samples[samples.symp_type_cat == i]
        logVL_central_all_curves = np.empty(
            (
                len(
                    samples_statistics_individual[
                        "study_participant_number_in_category"
                    ].unique()
                ),
                len(ts),
            )
        )
        logVL_median_all_curves = np.array(logVL_central_all_curves)
        for person in samples_statistics_individual[
            "study_participant_number_in_category"
        ].unique():
            temp_samples = samples_statistics_individual[
                samples_statistics_individual[
                    "study_participant_number_in_category"
                ]
                == person
            ]
            logVL_samples_per_person = np.empty(
                (
                    len(
                        samples_statistics_individual[
                            "posterior_draw_id"
                        ].unique()
                    ),
                    len(ts),
                )
            )
            for idx, t in enumerate(ts):
                temp_triangle_logVL = triangle_vl_fx(
                    t,
                    temp_samples.dp,
                    temp_samples.tp,
                    temp_samples.wp,
                    temp_samples.wr,
                )
                logVL_samples_per_person[:, idx] = temp_triangle_logVL
            depth = banddepth(logVL_samples_per_person)
            ix_depth = np.argsort(depth)[::-1]
            median_curve = logVL_samples_per_person[ix_depth[0], :]
            logVL_median_all_curves[person - 1, :] = median_curve

        fboxplot(logVL_median_all_curves, ts, ax1, colors[i - 1], wfactor=None)

        ax1.plot(
            ts * 0,
            0
            * triangle_vl_fx(
                ts,
                dp_samples_cat[0],
                tp_samples_cat[0],
                wp_samples_cat[0],
                wr_samples_cat[0],
            )
            - 2,  # so it gets cut off
            color="black",
            linestyle="solid",
            linewidth=4,
            alpha=1,
            label="Modeled median logVL\n(posterior samples)",
        )

    mask = (
        empirical_values_viral_load_cat.logVL
        != parameters_config["global_data_markers"]["skipped_test"]
    )
    logVL = empirical_values_viral_load_cat.logVL[mask]
    logVL[logVL == parameters_config["global_data_markers"]["lod"]] = 0.75
    logVL[
        logVL == parameters_config["global_data_markers"]["loq_lower"]
    ] = 1.75  # so they are all below 2
    ax1.scatter(
        empirical_values_viral_load_cat.days_since_symp_onset[mask]
        + np.random.normal(loc=0, scale=0.13, size=sum(mask)),
        logVL,
        color="grey",
        edgecolors="none",
        alpha=0.7,
    )

    ax1.scatter(
        empirical_values_viral_load_cat.days_since_symp_onset[mask],
        logVL * 0 - 2,
        color="grey",
        alpha=0.7,
        label="Measured logVL",
    )

    ax1.vlines(
        x=0,
        ymin=-2,
        ymax=14,
        linestyle="solid",
        color="black",
    )

    if i == 3:
        ax1.legend(frameon=True)

    ax1.text(
        -0.08,
        1.08,
        # A, B, C
        chr(ord("A") - 1 + i),
        transform=ax1.transAxes,
        fontsize=24,
        va="top",
        ha="right",
    )

    ax1.set_yticks(ticks=[3, 5, 7, 9])
    ax1.set_ylim([2.5, 9.5])
    ax1.set_xlim([minimum_day_considered, max_si])
    if i == 1:  # first symp category, leftmost plot
        ax1.set_ylabel(
            "RT-PCR measured viral load\n(log10 genome copies / mL)"
        )
    ax1.set_xlabel("Time since symptom onset (days)")

    ax1.set_title(titles[i - 1], pad=20, fontdict={"size": 32})

    # symptom improvement times subplot
    ax1 = fig.add_subplot(spec[1, i - 1])
    ax1.spines["top"].set_visible(False)
    ax1.spines["right"].set_visible(False)
    ax1.spines["left"].set_visible(False)
    ax1.get_xaxis().tick_bottom()
    ax1.get_yaxis().tick_left()
    ax1.tick_params(axis="x", direction="out")
    ax1.tick_params(axis="y", direction="out", length=0)
    # offset the spines
    for spine in ax1.spines.values():
        spine.set_position(("outward", 5))
    # put the grid behind
    ax1.set_axisbelow(True)

    violin_parts = ax1.violinplot(
        empirical_values_symp_improvement_list[i - 1],
        positions=[2],
        widths=0.7,
        showextrema=False,
        showmeans=False,
        showmedians=False,
        points=vln_points,
        vert=False,
    )

    for pc in violin_parts["bodies"]:
        pc.set_facecolor("grey")
        pc.set_alpha(0.7)

    bplot = ax1.boxplot(
        empirical_values_symp_improvement_list[i - 1],
        showmeans=True,
        positions=[2],
        whis=(25, 75),
        widths=0.5,
        patch_artist=True,
        showfliers=False,
        capprops=dict(linewidth=2),
        whiskerprops=dict(linewidth=2),
        meanprops=dict(color="k"),
        vert=False,
    )

    for patch in bplot["boxes"]:
        patch.set_facecolor("grey")
        patch.set_alpha(0.4)

    for item in ["means", "medians"]:
        tt1 = bplot[item]
        for tt2 in tt1:
            tt2.set_color("k")
            tt2.set_linewidth(2)
            if item == "means":
                tt2.set_linewidth(8)
                tt2.set_markerfacecolor("k")
                tt2.set_markersize(12)
                tt2.set_markeredgecolor("k")

    violin_parts = ax1.violinplot(
        weibull_samples[i - 1],
        positions=[1],
        showextrema=False,
        widths=0.7,
        showmeans=False,
        showmedians=False,
        points=vln_points,
        vert=False,
    )

    for pc in violin_parts["bodies"]:
        pc.set_facecolor("grey")
        pc.set_alpha(0.7)

    bplot = ax1.boxplot(
        weibull_samples[i - 1],
        showmeans=True,
        positions=[1],
        whis=(25, 75),
        widths=0.5,
        patch_artist=True,
        showfliers=False,
        capprops=dict(linewidth=2),
        whiskerprops=dict(linewidth=2),
        meanprops=dict(color="k"),
        vert=False,
    )

    for patch in bplot["boxes"]:
        patch.set_facecolor("grey")
        patch.set_alpha(0.4)

    for item in ["means", "medians"]:
        tt1 = bplot[item]
        for tt2 in tt1:
            tt2.set_color("k")
            tt2.set_linewidth(2)
            if item == "means":
                tt2.set_linewidth(8)
                tt2.set_markerfacecolor("k")
                tt2.set_markersize(12)
                tt2.set_markeredgecolor("k")

    ax1.vlines(
        x=0,
        ymin=-2,
        ymax=14,
        linestyle="solid",
        color="black",
    )

    ax1.set_ylim([0.5, 2.5])
    ax1.set_xticks(
        np.arange(0, parameters_config["scenario_parameters"]["max_si"], 5)
    )
    ax1.set_xlim(
        [
            minimum_day_considered,
            parameters_config["scenario_parameters"]["max_si"],
        ]
    )

    ax1.text(
        -0.08,
        1.08,
        # D, E, F
        chr(ord("A") + len(symp_categories) - 1 - 1 + i),
        transform=ax1.transAxes,
        fontsize=24,
        va="top",
        ha="right",
    )

    if i == 1:
        ax1.set_yticks([1, 2], ["Modeled", "Empirical"])
    else:
        ax1.set_yticks([])

    ax1.set_xlabel("Time to symptom improvement (days)")

    print(
        f"""symptom category {i} one-sided KS test (null hypothesis is
          that the empirical > modeled, which must be wrong because modeled has no censoring)
          is {stats.ks_2samp(empirical_values_symp_improvement_list[i - 1],
          weibull_samples[i - 1],
          alternative = "greater")}.
          Recall that we want the modeled > empirical, so p --> 0."""
    )

# Save the plot
fig.savefig(file_out + ".png", format="png", dpi=300, bbox_inches="tight")
fig.savefig(file_out + ".pdf", format="pdf", dpi=300, bbox_inches="tight")
