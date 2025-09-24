import os
from argparse import ArgumentParser

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import yaml
from dotenv import load_dotenv

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

plt.rcParams["figure.figsize"] = [12, 12]
font = {"family": "Nimbus Roman", "weight": "normal", "size": 32}

mpl.rc("font", **font)

load_dotenv()
output_dir = output_dir = os.path.expanduser(os.getenv("MAIN_OUTPUT_DIR"))

with open("products_config.yml", "r") as file:
    products_config = yaml.safe_load(file)

with open(os.path.join("analysis", "parameters.yml"), "r") as file:
    parameters_config = yaml.safe_load(file)

symp_categories = parameters_config["global_data_markers"][
    "symptom_categories"
]
symp_categories_reversed = symp_categories[::-1]

posterior_dir = os.path.join(
    output_dir,
    output_mode,
    products_config["vl_model"]["stan_posterior_directory"],
)

figure_path_out = os.path.join(
    output_dir,
    output_mode,
    products_config["products"]["directory"],
)

# Ensure the figure directory exists
os.makedirs(
    figure_path_out,
    exist_ok=True,
)

### plot of VL trajectory parameters

samples = pd.read_csv(
    os.path.join(
        posterior_dir,
        products_config["vl_model"]["posterior_pooled"],
    )
)

tp_samples = [
    samples.tp_mean[samples.symp_type_cat == i] for i in symp_categories
]
wr_samples = [
    samples.wr_mean[samples.symp_type_cat == i] for i in symp_categories
]
dp_samples = [
    samples.dp_mean[samples.symp_type_cat == i] for i in symp_categories
]
wp_samples = [
    samples.wp_mean[samples.symp_type_cat == i] for i in symp_categories
]

## box plots of clearance time, peak time by symptom category

file_out = os.path.join(
    figure_path_out,
    "figure_supplemental_box_plot_VL_pooled_params.png",
)

fig = plt.figure(figsize=(25, 24))
ax1 = fig.add_subplot(221)
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

bplot = ax1.boxplot(
    tp_samples,
    positions=symp_categories_reversed,
    showmeans=True,
    whis=(5, 95),
    patch_artist=True,
    showfliers=False,
    capprops=dict(linewidth=2),
    whiskerprops=dict(linewidth=2),
    meanprops=dict(color="k"),
    vert=False,
)
# quantiles = [[0.05, 0.95], [0.05, 0.95], [0.05, 0.95], [0.05, 0.95]])

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

ax1.set_yticks(
    symp_categories_reversed,
    [
        "1: Shortness\nof Breath",
        "2: Fever or body\naches",
        "3: Mild respiratory\nsymptoms",
        "4: Non-specific\nsymptoms",
    ],
)

ax1.set_xticks([-3, -2, -1, 0, 1, 2])
ax1.set_xlabel("Peak viral load timing\n(days since symptom onset)")
ax1.set_ylabel("Symptom Category", labelpad=20)

ax1.text(
    -0.50,
    1.08,
    "A",
    transform=ax1.transAxes,
    fontsize=24,
    va="top",
    ha="right",
)

ax1 = fig.add_subplot(222)
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

bplot = ax1.boxplot(
    wr_samples,
    positions=symp_categories_reversed,
    showmeans=True,
    whis=(5, 95),
    patch_artist=True,
    showfliers=False,
    capprops=dict(linewidth=2),
    whiskerprops=dict(linewidth=2),
    meanprops=dict(color="k"),
    vert=False,
)
# quantiles = [[0.05, 0.95], [0.05, 0.95], [0.05, 0.95], [0.05, 0.95]])

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

ax1.set_yticks(
    symp_categories_reversed,
    # [
    #     "1: Shortness\nof Breath",
    #     "2: Fever/\nBody aches",
    #     "3: Mild respiratory\nsymptoms",
    #     "4: Non-specific\nsymptoms",
    # ],
    ["", "", "", ""],
)

ax1.set_xticks([10, 15, 20, 25, 30])
ax1.set_xlim([7, 27])
ax1.set_xlabel("Clearance time (days)")
# ax1.set_ylabel("Symptom Category", labelpad = 20)

ax1.text(
    -0.08,
    1.08,
    "B",
    transform=ax1.transAxes,
    fontsize=24,
    va="top",
    ha="right",
)

ax1 = fig.add_subplot(223)
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

bplot = ax1.boxplot(
    wp_samples,
    positions=symp_categories_reversed,
    showmeans=True,
    whis=(5, 95),
    patch_artist=True,
    showfliers=False,
    capprops=dict(linewidth=2),
    whiskerprops=dict(linewidth=2),
    meanprops=dict(color="k"),
    vert=False,
)
# quantiles = [[0.05, 0.95], [0.05, 0.95], [0.05, 0.95], [0.05, 0.95]])

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

ax1.set_yticks(
    symp_categories_reversed,
    [
        "1: Shortness\nof Breath",
        "2: Fever or body\naches",
        "3: Mild respiratory\nsymptoms",
        "4: Non-specific\nsymptoms",
    ],
)

ax1.set_xticks([3, 6, 9, 12, 15])
ax1.set_xlim([2, 15])
ax1.set_xlabel("Proliferation time (days)")
ax1.set_ylabel("Symptom Category", labelpad=20)

ax1.text(
    -0.50,
    1.08,
    "C",
    transform=ax1.transAxes,
    fontsize=24,
    va="top",
    ha="right",
)

ax1 = fig.add_subplot(224)
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

bplot = ax1.boxplot(
    dp_samples,
    positions=symp_categories_reversed,
    showmeans=True,
    whis=(5, 95),
    patch_artist=True,
    showfliers=False,
    capprops=dict(linewidth=2),
    whiskerprops=dict(linewidth=2),
    meanprops=dict(color="k"),
    vert=False,
)
# quantiles = [[0.05, 0.95], [0.05, 0.95], [0.05, 0.95], [0.05, 0.95]])

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

ax1.set_yticks(
    symp_categories_reversed,
    # [
    #     "1: Shortness\nof Breath",
    #     "2: Fever/\nBody aches",
    #     "3: Mild respiratory\nsymptoms",
    #     "4: Non-specific\nsymptoms",
    # ],
    ["", "", "", ""],
)

ax1.set_xticks([2, 3, 4, 5, 6, 7, 8])
ax1.set_xlim([2, 8])
ax1.set_xlabel("Peak viral load (log IU / mL)")
# ax1.set_ylabel("Symptom Category", labelpad = 20)

ax1.text(
    -0.08,
    1.08,
    "D",
    transform=ax1.transAxes,
    fontsize=24,
    va="top",
    ha="right",
)

# Adjust layout
fig.tight_layout()

# Save the plot
fig.savefig(file_out, format="png", dpi=300, bbox_inches="tight")


def summary_statistics_on_posterior_params(
    variable_name, mathematical_name, symp_cat, running_dataframe_of_stats=None
):
    if running_dataframe_of_stats is None:
        running_dataframe_of_stats = pd.DataFrame(
            columns=[
                "pooled parameter name in stan code (following Kissler et al.)",
                "mathematical name of variable",
                "symptom category",
                "mean",
                "middle",
                "q2.5",
                "q97.5",
                "q5",
                "q95",
            ]
        )

    running_dataframe_of_stats = pd.concat(
        [
            running_dataframe_of_stats,
            pd.DataFrame(
                [
                    [
                        variable_name,
                        mathematical_name,
                        symp_cat,
                        np.mean(
                            samples[variable_name][
                                samples.symp_type_cat == symp_cat
                            ]
                        ),
                        np.median(
                            samples[variable_name][
                                samples.symp_type_cat == symp_cat
                            ]
                        ),
                        np.percentile(
                            samples[variable_name][
                                samples.symp_type_cat == symp_cat
                            ],
                            q=2.5,
                        ),
                        np.percentile(
                            samples[variable_name][
                                samples.symp_type_cat == symp_cat
                            ],
                            q=97.5,
                        ),
                        np.percentile(
                            samples[variable_name][
                                samples.symp_type_cat == symp_cat
                            ],
                            q=5,
                        ),
                        np.percentile(
                            samples[variable_name][
                                samples.symp_type_cat == symp_cat
                            ],
                            q=95,
                        ),
                    ]
                ],
                columns=running_dataframe_of_stats.columns,
            ),
        ],
        axis=0,
    )

    return running_dataframe_of_stats


vars = {
    "si_beta_0_exponentiated": "exp(lambda_0); baseline scale parameter for the Weibull symptom improvement distribution",
    "si_beta_wr": "lambda_c; Weibull scale association between the individual deviation in the clearance time and symptom improvement time",
    "si_shape": "alpha; shape parameter for the Weibull symptom improvement distribution",
    "wp_mean": "r; average proliferation period (days)",
    "dp_mean": "p; average peak viral load (log10 IU / mL)",
    "tp_mean": "t^p; average time of peak viral load since symptom onset (days)",
    "wr_mean": "c; average clearance time (days)",
    "wp_std": "sigma_r; standard deviation in proliferation period (days)",
    "dp_std": "sigma_p; standard deviation in peak viral load (log10 IU / mL)",
    "tp_std": "sigma_t^p; standard deviation in time of peak viral load since symptom onset (days)",
    "wr_std": "sigma_c; standard deviation in clearance time (days)",
    "antigen_50": "k_ant; logVL at which antigen positivity is 50% (log10 IU / mL)",
    "sigma_antigen": "sigma_ant; variability around k_ant in antigen positivity",
    "culture_50": "k_cul; logVL at which culture positivity is 50% at time of peak VL (log10 IU / mL)",
    "culture_beta": "b_cul; daily change in logVL needed to be 50% culture positive (log10 IU / mL / day)",
    "sigma_culture": "sigma_cul; variability around k_cul in culture positivity (log10 IU / mL)",
}

running_dataframe_of_stats = None
for symp_cat in symp_categories:
    for var in vars.keys():
        running_dataframe_of_stats = summary_statistics_on_posterior_params(
            var, vars[var], symp_cat, running_dataframe_of_stats
        )

table_path_out = os.path.join(
    output_dir,
    output_mode,
    products_config["products"]["directory"],
)

# Ensure the figure directory exists
os.makedirs(
    table_path_out,
    exist_ok=True,
)

running_dataframe_of_stats.to_csv(
    os.path.join(
        table_path_out,
        "table_supplemental_summary_statistics_VL_pooled_params.csv",
    ),
    index=False,
)
