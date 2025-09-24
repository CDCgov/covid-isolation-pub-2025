import os
from argparse import ArgumentParser

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import yaml
from dotenv import load_dotenv

plt.rcParams["figure.figsize"] = [12, 12]
font = {"family": "Nimbus Roman", "weight": "normal", "size": 22}

mpl.rc("font", **font)

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

load_dotenv()
output_dir = output_dir = os.path.expanduser(os.getenv("MAIN_OUTPUT_DIR"))

with open("products_config.yml", "r") as file:
    products_config = yaml.safe_load(file)

with open(os.path.join("analysis", "parameters.yml"), "r") as file:
    parameters_config = yaml.safe_load(file)

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

file_out = os.path.join(
    figure_path_out,
    "figure_supplemental_example_VL_fits.png",
)

# Ensure the figure directory exists
os.makedirs(
    figure_path_out,
    exist_ok=True,
)


def plot_samples_data(
    ax1,
    data,
    samples,
    lod,
    loq,
    min_day,
    max_day,
    legend=False,
    loq_upper=False,
):
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
    for i in range(max(samples.sample_num)):
        temp = samples[samples.sample_num == i]
        ax1.plot(
            temp.time,
            temp.logVL,
            color=model_draws_color,
            linestyle="solid",
            alpha=0.2,
        )

    ax1.hlines(
        y=lod,
        xmin=min_day,
        xmax=max_day,
        linestyle="solid",
        linewidth=3,
        color=negative_color,
    )

    ax1.hlines(
        y=loq,
        xmin=min_day,
        xmax=max_day,
        linestyle="solid",
        linewidth=3,
        color=loq_color,
    )

    ax1.scatter(
        data.time[np.logical_and(data.antigen == 1, data.culture == 1)],
        data.logVL[np.logical_and(data.antigen == 1, data.culture == 1)]
        + (
            data.logVL[np.logical_and(data.antigen == 1, data.culture == 1)]
            == 0
        )
        * lod
        + (
            data.logVL[np.logical_and(data.antigen == 1, data.culture == 1)]
            == 1
        )
        * (loq - 1),  # minus 1 because the logVL value for these guys is 1
        marker="+",
        color="magenta",
        label="Antigen +, Culture +",
        s=size,
        zorder=10,
        clip_on=False,
    )

    ax1.scatter(
        data.time[np.logical_and(data.antigen == 1, data.culture == 0)],
        data.logVL[np.logical_and(data.antigen == 1, data.culture == 0)]
        + (
            data.logVL[np.logical_and(data.antigen == 1, data.culture == 0)]
            == 0
        )
        * lod
        + (
            data.logVL[np.logical_and(data.antigen == 1, data.culture == 0)]
            == 1
        )
        * (loq - 1),
        marker="+",
        color="green",
        label="Antigen +, Culture -",
        s=size,
        zorder=10,
        clip_on=False,
    )

    ax1.scatter(
        data.time[np.logical_and(data.antigen == 0, data.culture == 1)],
        data.logVL[np.logical_and(data.antigen == 0, data.culture == 1)]
        + (
            data.logVL[np.logical_and(data.antigen == 0, data.culture == 1)]
            == 0
        )
        * lod
        + (
            data.logVL[np.logical_and(data.antigen == 0, data.culture == 1)]
            == 1
        )
        * (loq - 1),
        marker="x",
        color="magenta",
        label="Antigen -, Culture +",
        s=size,
        zorder=10,
        clip_on=False,
    )

    ax1.scatter(
        data.time[np.logical_and(data.antigen == 0, data.culture == 0)],
        data.logVL[np.logical_and(data.antigen == 0, data.culture == 0)]
        + (
            data.logVL[np.logical_and(data.antigen == 0, data.culture == 0)]
            == 0
        )
        * lod
        + (
            data.logVL[np.logical_and(data.antigen == 0, data.culture == 0)]
            == 1
        )
        * (loq - 1),
        marker="x",
        color="green",
        label="Antigen -, Culture -",
        s=size,
        zorder=10,
        clip_on=False,
    )

    ax1.scatter(
        data.time[
            np.logical_and(
                data.culture != 0,
                np.logical_and(data.culture != 1, data.antigen == 0),
            )
        ],
        data.logVL[
            np.logical_and(
                data.culture != 0,
                np.logical_and(data.culture != 1, data.antigen == 0),
            )
        ]
        + (
            data.logVL[
                np.logical_and(
                    data.culture != 0,
                    np.logical_and(data.culture != 1, data.antigen == 0),
                )
            ]
            == 0
        )
        * lod
        + (
            data.logVL[
                np.logical_and(
                    data.culture != 0,
                    np.logical_and(data.culture != 1, data.antigen == 0),
                )
            ]
            == 1
        )
        * (loq - 1),
        marker="x",
        color="black",
        label="Antigen -, Culture NA",
        s=size,
        zorder=10,
        clip_on=False,
    )

    ax1.scatter(
        data.time[
            np.logical_and(
                data.culture != 0,
                np.logical_and(data.culture != 1, data.antigen == 1),
            )
        ],
        data.logVL[
            np.logical_and(
                data.culture != 0,
                np.logical_and(data.culture != 1, data.antigen == 1),
            )
        ]
        + (
            data.logVL[
                np.logical_and(
                    data.culture != 0,
                    np.logical_and(data.culture != 1, data.antigen == 1),
                )
            ]
            == 0
        )
        * lod
        + (
            data.logVL[
                np.logical_and(
                    data.culture != 0,
                    np.logical_and(data.culture != 1, data.antigen == 1),
                )
            ]
            == 1
        )
        * (loq - 1),
        marker="+",
        color="black",
        label="Antigen +, Culture NA",
        s=size,
        zorder=10,
        clip_on=False,
    )

    ax1.scatter(
        data.time[
            np.logical_and(
                data.antigen != 0,
                np.logical_and(data.antigen != 1, data.culture == 0),
            )
        ],
        data.logVL[
            np.logical_and(
                data.antigen != 0,
                np.logical_and(data.antigen != 1, data.culture == 0),
            )
        ]
        + (
            data.logVL[
                np.logical_and(
                    data.antigen != 0,
                    np.logical_and(data.antigen != 1, data.culture == 0),
                )
            ]
            == 0
        )
        * lod
        + (
            data.logVL[
                np.logical_and(
                    data.antigen != 0,
                    np.logical_and(data.antigen != 1, data.culture == 0),
                )
            ]
            == 1
        )
        * (loq - 1),
        marker="o",
        color="green",
        label="Antigen NA, Culture -",
        s=size,
        zorder=10,
        clip_on=False,
    )

    ax1.scatter(
        data.time[
            np.logical_and(
                data.antigen != 0,
                np.logical_and(data.antigen != 1, data.culture == 1),
            )
        ],
        data.logVL[
            np.logical_and(
                data.antigen != 0,
                np.logical_and(data.antigen != 1, data.culture == 1),
            )
        ]
        + (
            data.logVL[
                np.logical_and(
                    data.antigen != 0,
                    np.logical_and(data.antigen != 1, data.culture == 1),
                )
            ]
            == 0
        )
        * lod
        + (
            data.logVL[
                np.logical_and(
                    data.antigen != 0,
                    np.logical_and(data.antigen != 1, data.culture == 1),
                )
            ]
            == 1
        )
        * (loq - 1),
        marker="o",
        color="magenta",
        label="Antigen NA, Culture +",
        s=size,
        zorder=10,
        clip_on=False,
    )

    ax1.scatter(
        data.time[
            np.logical_and(
                np.logical_and(data.antigen != 0, data.antigen != 1),
                np.logical_and(data.culture != 0, data.culture != 1),
            )
        ],
        data.logVL[
            np.logical_and(
                np.logical_and(data.antigen != 0, data.antigen != 1),
                np.logical_and(data.culture != 0, data.culture != 1),
            )
        ]
        + (
            data.logVL[
                np.logical_and(
                    np.logical_and(data.antigen != 0, data.antigen != 1),
                    np.logical_and(data.culture != 0, data.culture != 1),
                )
            ]
            == 0
        )
        * lod
        + (
            data.logVL[
                np.logical_and(
                    np.logical_and(data.antigen != 0, data.antigen != 1),
                    np.logical_and(data.culture != 0, data.culture != 1),
                )
            ]
            == 1
        )
        * (loq - 1),
        marker="o",
        color="black",
        label="Antigen NA, Culture NA",
        s=size,
        zorder=10,
        clip_on=False,
    )

    ax1.fill_between(
        x=np.arange(min_day, max_day, 0.01),
        y1=np.arange(min_day, max_day, 0.01) * 0 - 1,
        y2=np.arange(min_day, max_day, 0.01) * 0 + lod,
        color="grey",
        alpha=0.4,
    )

    ax1.fill_between(
        x=np.arange(min_day, max_day, 0.01),
        y1=np.arange(min_day, max_day, 0.01) * 0 - lod,
        y2=np.arange(min_day, max_day, 0.01) * 0 + loq,
        color="grey",
        alpha=0.2,
    )

    if loq_upper:
        ax1.fill_between(
            x=np.arange(min_day, max_day, 0.01),
            y1=np.arange(min_day, max_day, 0.01) * 0 + 9,
            y2=np.arange(min_day, max_day, 0.01) * 0 + 13,
            color="grey",
            alpha=0.2,
        )
        ax1.text(
            -0.5, 9.6, "Above Upper RT-PCR LOQ", fontdict={"weight": "normal"}
        )

    ax1.set_yticks(ticks=[2, 4, 6, 8, 10, 12])
    ax1.set_yticks(ticks=[2, 4, 6, 8, 10, 12])

    ax1.text(
        -0.5, lod - 0.6, "Below RT-PCR LOD", fontdict={"weight": "normal"}
    )
    ax1.text(
        -0.5, loq - 0.6, "Below RT-PCR LOQ", fontdict={"weight": "normal"}
    )

    ax1.set_ylim([-1, 10.5])
    ax1.set_xlim([min_day, max_day])
    ax1.set_ylabel("Viral load (log10 IU / mL)")
    ax1.set_xlabel("Time since symptom onset (days)")

    if legend:
        ax1.legend()


# show four emblematic individuals
# picked from cat 2 because of good culture data there -- so hence the "2.csv"
vl_data = pd.read_csv(
    os.path.join(
        posterior_dir,
        "data_" + products_config["vl_model"]["run_name"] + "2.csv",
    )
)

vl_samples = pd.read_csv(
    os.path.join(
        posterior_dir,
        "samples_" + products_config["vl_model"]["run_name"] + "2.csv",
    )
)

size = 150

lod_inherent = parameters_config["INHERENT"]["lod"]
loq_inherent = parameters_config["INHERENT"]["loq_lower"]
lod_rvtn = parameters_config["RVTN"]["lod"]
loq_rvtn = parameters_config["RVTN"]["loq_lower"]

model_draws_color = "grey"  # "#3266ab"  # "#678dbf"
negative_color = "#cb8fe3"
loq_color = "orange"

fig = plt.figure(figsize=(20, 20))

if output_mode == "production":
    ax1 = fig.add_subplot(221)
    data = vl_data[vl_data.id == 3]
    samples = vl_samples[vl_samples.id == 3]
    plot_samples_data(
        ax1,
        data,
        samples,
        lod_inherent,
        loq_inherent,
        min_day=-10,
        max_day=50,
        loq_upper=True,
    )

    ax1.text(
        -0.08,
        1.08,
        "A",
        transform=ax1.transAxes,
        fontsize=24,
        va="top",
        ha="right",
    )

    ax1 = fig.add_subplot(222)
    data = vl_data[vl_data.id == 13]
    samples = vl_samples[vl_samples.id == 13]
    plot_samples_data(
        ax1, data, samples, lod_rvtn, loq_rvtn, min_day=-5, max_day=20
    )

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
    data = vl_data[vl_data.id == 5]
    samples = vl_samples[vl_samples.id == 5]
    plot_samples_data(
        ax1,
        data,
        samples,
        lod_inherent,
        loq_inherent,
        min_day=-5,
        max_day=90,
        legend=True,
        loq_upper=True,
    )

    ax1.text(
        -0.08,
        1.08,
        "C",
        transform=ax1.transAxes,
        fontsize=24,
        va="top",
        ha="right",
    )

    ax1 = fig.add_subplot(224)
    data = vl_data[vl_data.id == 81]
    samples = vl_samples[vl_samples.id == 81]
    plot_samples_data(
        ax1, data, samples, lod_rvtn, loq_rvtn, min_day=-5, max_day=15
    )

    ax1.text(
        -0.08,
        1.08,
        "D",
        transform=ax1.transAxes,
        fontsize=24,
        va="top",
        ha="right",
    )

else:
    print(
        "NOT making example VL plot because script is not running in production mode."
    )

# Adjust layout
fig.tight_layout()
# Save the plot
fig.savefig(file_out, format="png", dpi=300, bbox_inches="tight")
