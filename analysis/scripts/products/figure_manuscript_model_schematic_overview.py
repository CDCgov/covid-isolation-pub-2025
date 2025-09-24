import os
from argparse import ArgumentParser

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import yaml
from dotenv import load_dotenv
from scipy import special

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

with open("products_config.yml", "r") as file:
    products_config = yaml.safe_load(file)

with open(os.path.join("analysis", "parameters.yml"), "r") as file:
    parameters_config = yaml.safe_load(file)

font = {"family": "Nimbus Roman", "weight": "normal", "size": 22}

mpl.rc("font", **font)

load_dotenv()
output_dir = os.path.expanduser(os.getenv("MAIN_OUTPUT_DIR"))

stan_posterior_output_dir = os.path.join(
    output_dir,
    output_mode,
    products_config["vl_model"]["stan_posterior_directory"],
)

figure_path_out = os.path.join(
    output_dir,
    output_mode,
    products_config["products"]["directory"],
)

# Ensure the figures directory exists
os.makedirs(figure_path_out, exist_ok=True)


def triangle_vl_fx(t, dp, tp, wp, wr):
    return (t <= tp) * (dp / wp) * (t - (tp - wp)) + (t > tp) * (
        dp - (dp / wr) * (t - tp)
    )


max_si = parameters_config["scenario_parameters"]["max_si"]
minimum_day_considered = parameters_config["scenario_parameters"][
    "minimum_days_before_symptom_onset_considered"
]

samples = pd.read_csv(
    os.path.join(
        stan_posterior_output_dir,
        products_config["vl_model"]["posterior_pooled"],
    )
)

samples = samples[samples.symp_type_cat.isin(range(1, 4))]

symp_cat_weights = pd.read_csv(
    os.path.join(
        output_dir,
        output_mode,
        products_config["preprocessed_data"]["directory"],
        products_config["scenarios"]["symp_cat_weights"],
    )
)

samples = samples.merge(symp_cat_weights)

samples_by_cat = samples.groupby("symp_type_cat").mean().reset_index()

samples_by_cat["tp_mean"] *= samples_by_cat["wt"]
samples_by_cat["dp_mean"] *= samples_by_cat["wt"]
samples_by_cat["wp_mean"] *= samples_by_cat["wt"]
samples_by_cat["wr_mean"] *= samples_by_cat["wt"]
samples_by_cat["sigma"] *= samples_by_cat["wt"]
samples_by_cat["sigma_antigen"] *= samples_by_cat["wt"]
samples_by_cat["antigen_50"] *= samples_by_cat["wt"]
samples_by_cat["sigma_culture"] *= samples_by_cat["wt"]
samples_by_cat["culture_50"] *= samples_by_cat["wt"]
samples_by_cat["culture_beta"] *= samples_by_cat["wt"]

print(samples_by_cat)

samples_overall = samples_by_cat.sum(axis=0)

print(samples_overall)

dt = parameters_config["scenario_parameters"]["infectiousness_integration_dt"]
ts = np.arange(minimum_day_considered, max_si, dt)
tri_vl = triangle_vl_fx(
    ts,
    samples_overall["dp_mean"],
    samples_overall["tp_mean"],
    samples_overall["wp_mean"],
    samples_overall["wr_mean"],
)
tri_vl_samples = triangle_vl_fx(
    np.arange(minimum_day_considered, max_si, 1),
    samples_overall["dp_mean"],
    samples_overall["tp_mean"],
    samples_overall["wp_mean"],
    samples_overall["wr_mean"],
)

## schematic diagram of area under curve as transmission averted under proposed
## and current guidances for the two VL transformations to show comparisons

isolation_starting_time = 0
previous_isolation_end = 5
previous_postisolation_end = 10
updated_isolation_end = 3
updated_postisolation_end = updated_isolation_end + 5
lod = 1.5

tri = np.maximum(tri_vl, 1.5) - lod
tri = tri / (sum(tri) * dt)
sig = special.expit(
    (1 / samples_overall["sigma_antigen"])
    * (tri_vl - samples_overall["antigen_50"])
)
sig = sig / (sum(sig) * dt)
cul = special.expit(
    (1 / samples_overall["sigma_culture"])
    * (
        tri_vl
        - (
            samples_overall["culture_50"]
            + samples_overall["culture_beta"]
            * (ts - samples_overall["tp_mean"])
        )
    )
)
cul = cul / (sum(cul) * dt)

## example VL trajectory

file_out = os.path.join(
    figure_path_out, "figure_manuscript_model_schematic_overview"
)

# color scheme with inspiration from colorbrewer2 (qualitative paired palette)
vl_transformation_colors = {
    "logVL": "#a6bddb",
    "antigen": "#fb9a99",
    "culture": "#7fc97f",
}

# these are just colors from dark2 color palette
updated_previous_colors = {"previous": "#1b9e77", "updated": "#d95f02"}

# want two columns in the first row
# the second "row" is actually two rows, each half the height of the first row
# the second "rows" get divided into three subplots each
# since we use numpy slicing to get at the grid arrays, need
# to divide the space into 6 columns (3 * 2)
height = 20
fig = plt.figure(figsize=(20, height), layout="constrained")
spec = fig.add_gridspec(ncols=6, nrows=height)
ax1 = fig.add_subplot(spec[0:9, :3])
ax1.text(
    -0.1,
    1.08,
    "A",
    transform=ax1.transAxes,
    fontsize=24,
    va="top",
    ha="right",
)
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

np.random.seed(108)

ax1.scatter(
    np.arange(minimum_day_considered, max_si, 1),
    np.maximum(
        tri_vl_samples
        + np.random.normal(
            loc=0, scale=samples_overall["sigma"], size=tri_vl_samples.shape
        ),
        lod,
    ),
    color="black",
    label="RT-PCR measured",
)

ax1.plot(
    ts,
    tri_vl,
    linewidth=4,
    color=vl_transformation_colors["logVL"],
    label="Modeled underlying VL",
)

# add a line annotating the clearance time
ax1.errorbar(
    x=10,
    y=3.3,
    xerr=np.array(
        [[10 - samples_overall["tp_mean"]], [samples_overall["wr_mean"] - 10]]
    ),
    ecolor="grey",
    linewidth=2,
)
ax1.errorbar(
    x=samples_overall["tp_mean"], y=3.3, yerr=0.3, ecolor="grey", linewidth=2
)
ax1.errorbar(
    x=samples_overall["wr_mean"], y=3.3, yerr=0.3, ecolor="grey", linewidth=2
)
ax1.annotate(
    "Clearance time",
    fontstyle="italic",
    xy=((samples_overall["tp_mean"] + samples_overall["wr_mean"]) / 2, 3.4),
    horizontalalignment="center",
    c="grey",
)

arrow = mpl.patches.Arrow(
    samples_overall["wr_mean"] - 1, 3.65, 4, 0.75, color="grey", width=0.30
)
ax1.add_patch(arrow)
ax1.annotate(
    "Symptom improvement\ntime prediction",
    xy=(19, 4.5),
    horizontalalignment="center",
    c="grey",
)

ax1.hlines(
    y=lod,
    xmin=minimum_day_considered,
    xmax=max_si,
    linestyle="solid",
    linewidth=3,
    color="brown",
    label="RT-PCR limit of detection",
)

ax1.fill_between(
    x=np.arange(minimum_day_considered, max_si, 0.01),
    y1=np.arange(minimum_day_considered, max_si, 0.01) * 0 - 2,
    y2=np.arange(minimum_day_considered, max_si, 0.01) * 0 + lod,
    color="grey",
    alpha=0.7,
)

ax1.fill_between(
    x=np.arange(minimum_day_considered, max_si, 0.01),
    y1=np.arange(minimum_day_considered, max_si, 0.01) * 0 - 2,
    y2=np.arange(minimum_day_considered, max_si, 0.01) * 0 + lod,
    color="grey",
    alpha=0.4,
)

ax1.legend(frameon=True)

ax1.set_xlim([minimum_day_considered, max_si])
ax1.set_xticks(np.arange(minimum_day_considered, max_si, 5))
ax1.set_yticks(ticks=[3, 5, 7, 9])
ax1.set_ylim([1, 9.5])
ax1.set_xlim([minimum_day_considered, max_si])
ax1.set_ylabel("Viral load (log10 genome copies / mL)")
ax1.set_xlabel("Time since symptom onset (days)")

## schematic of transformations

ax1 = fig.add_subplot(spec[0:9, 3:])
ax1.text(
    -0.1,
    1.08,
    "B",
    transform=ax1.transAxes,
    fontsize=24,
    va="top",
    ha="right",
)
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

ax1.plot(
    ts,
    tri,
    linewidth=4,
    c=vl_transformation_colors["logVL"],
    label="Log VL",
)
ax1.plot(
    ts,
    sig,
    linewidth=4,
    c=vl_transformation_colors["antigen"],
    label="Antigen positivity",
)
ax1.plot(
    ts,
    cul,
    linewidth=4,
    c=vl_transformation_colors["culture"],
    label="Culture positivity",
)

ax1.set_xlim([minimum_day_considered, max_si])
ax1.set_xticks(np.arange(minimum_day_considered, max_si, 5))
ax1.legend(frameon=True)  # , loc=(0.45, 0.65)
ax1.set_yticks([])
ax1.set_xlim([minimum_day_considered, max_si])
ax1.set_ylabel("Estimated infectiousness")
ax1.set_xlabel("Time since symptom onset (days)")
ax1.set_ylim(bottom=0.0)

## schematic auc


def schematic_auc_infectiousness_averted(
    ax1,
    label1,
    isolation_start,
    isolation_end,
    postisolation_end,
    color,
    legend=False,
    updated=False,
    print_descriptor=False,
    describe_x=False,
):
    ax1.text(
        -0.1,
        1.08,
        label1,
        transform=ax1.transAxes,
        fontsize=24,
        va="top",
        ha="right",
    )
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

    ax1.fill_between(
        np.arange(isolation_start, isolation_end, dt),
        y1=tri[
            np.where(
                abs(ts - isolation_start) == min(abs(ts - isolation_start))
            )[0][0] : np.where(
                abs(ts - isolation_end) == min(abs(ts - isolation_end))
            )[
                0
            ][
                0
            ]
        ],
        alpha=0.8,
        color=color,
        label="Isolation",
    )

    ax1.plot(ts, tri, linewidth=4, c=color)

    ax1.fill_between(
        np.arange(isolation_end, postisolation_end, dt),
        y1=tri[
            np.where(abs(ts - isolation_end) == min(abs(ts - isolation_end)))[
                0
            ][0] : np.where(
                abs(ts - postisolation_end) == min(abs(ts - postisolation_end))
            )[
                0
            ][
                0
            ]
        ],
        alpha=0.4,
        color=color,
        label="Post-isolation",
    )

    if legend:
        ax1.legend(frameon=True)

    ax1.set_xlim([minimum_day_considered, max_si])
    ax1.set_xticks(np.arange(minimum_day_considered, max_si, 5))
    ax1.set_yticks([])
    ax1.set_ylim(bottom=0.0)
    if print_descriptor:
        if updated:
            ylabel = "Estimated infectiousness"
            ax1.annotate(
                "Updated guidance",
                xy=(-0.2, 0.5),
                xycoords="axes fraction",
                rotation=90,
                weight="bold",
                verticalalignment="center",
            )
        else:
            ylabel = "Estimated infectiousness"
            ax1.annotate(
                "Previous guidance",
                xy=(-0.2, 0.5),
                xycoords="axes fraction",
                rotation=90,
                weight="bold",
                verticalalignment="center",
            )
    else:
        ylabel = ""  # nothing, because shared y axis idea
    ax1.set_ylabel(ylabel)
    if describe_x:
        ax1.set_xlabel("Time since symptom onset (days)")


isolation_start = 0
pad = 20
# the upper row (spec[1, -]) is previous guidance
# we show three conditions -- (1) symp improv time of less than 5 days
# and not shortness of breath, (2) symp improv time of greater than 5 days
# and not shortness of breath,
# and (3) symp immprov time of greater than 5 days and shortness of breath

# previous -- shortness of breath (moderate illness)
ax1 = fig.add_subplot(spec[10:15, :2])
ax1.set_title("Shortness of breath/\nModerate illness", pad=pad)
schematic_auc_infectiousness_averted(
    ax1,
    label1="C",
    isolation_start=isolation_start,
    isolation_end=10,
    postisolation_end=10,
    color=updated_previous_colors["previous"],
    print_descriptor=True,
)

# previous -- symp improv time of greater than 5 days
ax1 = fig.add_subplot(spec[10:15, 2:4])
ax1.set_title("Symptom improvement time \n> 5 days", pad=pad)
schematic_auc_infectiousness_averted(
    ax1,
    label1="D",
    isolation_start=isolation_start,
    isolation_end=7,
    postisolation_end=10,
    color=updated_previous_colors["previous"],
)

# previous -- symp improv of less than five days
ax1 = fig.add_subplot(spec[10:15, 4:])
ax1.set_title("Symptom improvement time \n≤ 5 days", pad=pad)
schematic_auc_infectiousness_averted(
    ax1,
    label1="E",
    isolation_start=isolation_start,
    isolation_end=5,
    postisolation_end=10,
    color=updated_previous_colors["previous"],
    legend=True,
)

# updated -- shortness of breath and symp improv time of greater than 5 days
ax1 = fig.add_subplot(spec[15:, :2])
schematic_auc_infectiousness_averted(
    ax1,
    label1="F",
    isolation_start=isolation_start,
    isolation_end=8,
    postisolation_end=13,
    color=updated_previous_colors["updated"],
    print_descriptor=True,
    updated=True,
    describe_x=True,
)

# updated -- symp improv time of greater than 5 days
ax1 = fig.add_subplot(spec[15:, 2:4])
schematic_auc_infectiousness_averted(
    ax1,
    label1="G",
    isolation_start=isolation_start,
    isolation_end=7,
    postisolation_end=12,
    color=updated_previous_colors["updated"],
    updated=True,
    describe_x=True,
)

# updated -- symp improv time of less than 5 days
ax1 = fig.add_subplot(spec[15:, 4:])
schematic_auc_infectiousness_averted(
    ax1,
    label1="H",
    isolation_start=isolation_start,
    isolation_end=3,
    postisolation_end=8,
    color=updated_previous_colors["updated"],
    updated=True,
    legend=True,
    describe_x=True,
)

# Save the plot
fig.savefig(file_out + ".png", format="png", dpi=300, bbox_inches="tight")
fig.savefig(file_out + ".pdf", format="pdf", dpi=300, bbox_inches="tight")
