import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def get_args(allowed_options, default_option="development"):
    """
    This function takes in a list of allowed options and a default option. It checks if the "--args" flag is present in the command-line arguments and retrieves the selected option. If no argument is provided or the argument value is invalid, it defaults to the default option. The function returns the selected option and a message that can be printed

    Parameters:
    - allowed_options (list)
    - default_option (str)

    Returns
    - selected_option (str)
    - printme (str)
    """

    selected_option = default_option
    printme = ""

    # ensure default option is allowed
    if default_option not in allowed_options:
        raise ValueError("Selected default option is not an allowed options")

    if "--args" in sys.argv:
        index = sys.argv.index("--args")
        if len(sys.argv) > index + 1:
            arg_value = sys.argv[index + 1]
            if arg_value in allowed_options:
                printme = "Selected option: " + arg_value
                selected_option = arg_value
            else:
                printme = (
                    "Invalid argument. Allowed options are: "
                    + ", ".join(allowed_options)
                )
                printme += "\nDefaulting to '" + default_option + "'"
        else:
            printme = (
                "No argument provided. Defaulting to '" + default_option + "'"
            )
    else:
        printme = (
            "No --args flag found. Defaulting to '" + default_option + "'"
        )

    return selected_option, printme


def add_symp_duration_columns(input_df):
    """
    This function adds 3 static (not timeseries) variables to the dataframe:
    symp_type_cat: categories 1-4 based on symptoms
    symp_duration: time to symptom improvement
    symp_duration_censored: whether or not time to symptom improvement is censored

    """
    df = input_df.copy()

    # Initialize columns
    df["symp_type_cat"] = np.nan
    df["symp_duration"] = np.nan
    df["symp_duration_censored"] = np.nan

    # Identify patients for each category
    # Category 1: Shortness of breath
    cat1_patients = df[df["symp_short_breath"] == 1]["id"].unique()

    # Category 2: Fever and/or body aches, excluding Category 1
    cat2_patients = df[
        ~(df["id"].isin(cat1_patients))
        & ((df["symp_fever"] == 1) | (df["symp_body_ache"] == 1))
    ]["id"].unique()

    # Category 3: Other specific symptoms, excluding Categories 1 and 2
    other_symps = [
        "symp_smell_loss",
        "symp_cough",
        "symp_sore_throat",
        "symp_nose",
        "symp_wheezing",
        "symp_chest_pain",
    ]
    cat3_patients = df[
        ~(df["id"].isin(cat1_patients))
        & ~(df["id"].isin(cat2_patients))
        & (df[other_symps].sum(axis=1) > 0)
    ]["id"].unique()

    # Category 4: Symptomatic but not in Categories 1, 2, or 3
    cat4_patients = df[
        ~(df["id"].isin(cat1_patients))
        & ~(df["id"].isin(cat2_patients))
        & ~(df["id"].isin(cat3_patients))
        & (df["symptomatic_ever"] == 1)
    ]["id"].unique()

    # Assign categories
    df.loc[df["id"].isin(cat1_patients), "symp_type_cat"] = 1
    df.loc[df["id"].isin(cat2_patients), "symp_type_cat"] = 2
    df.loc[df["id"].isin(cat3_patients), "symp_type_cat"] = 3
    df.loc[df["id"].isin(cat4_patients), "symp_type_cat"] = 4

    # Loop through all symptomatic people and calculate 'symp_duration' and 'symp_duration_censored' for each category
    for id_ in df[df["symptomatic_ever"] == 1]["id"].unique():
        df_id = df[df["id"] == id_].copy(deep=True)
        category = df_id["symp_type_cat"].iloc[0]

        ## note: 'final_day' refers to symptom improvement day
        ## note 2: 'final_day' was previously compared with the last survey day to determine
        ##         censorship, but due to fever edge cases, 'censor_final_day' was created.
        ##          'censor_final_day' equals 'final_day' except for those fever edge cases.

        # Remove days/rows with no symptom data (particularly relevant for INHERENT)
        df_id = df_id[pd.notnull(df_id["symp_count"])]
        df_id = df_id.reset_index()

        # Find the last day of data for this ID to handle censoring
        last_day_of_data = df_id[
            "days_since_symp_onset"
        ].max()  ## Does this produce any NA's?
        final_day_fever = 0
        final_day_ache = 0  # don't need to define final_day_short_breath ahead of time since always occurs for cat 1

        if category == 1:
            fever_ever = False

            # Check if the patient has a fever, and if so, update their final day of fever (+1)
            if not df_id[df_id["symp_fever"] == 1].empty:
                final_day_fever = (
                    df_id[df_id["symp_fever"] == 1][
                        "days_since_symp_onset"
                    ].max()
                    + 1  # Add one day because that would be the first day of symptoms improving
                )
                fever_ever = True

            final_day_short_breath = (
                df_id[df_id["symp_short_breath"] == 1][
                    "days_since_symp_onset"
                ].max()
                + 1
            )

            final_day = max(final_day_fever, final_day_short_breath)
            censor_final_day = final_day
            if (fever_ever) & (
                df_id[df_id["symp_fever"] == 1]["days_since_symp_onset"].max()
                == (final_day - 1)
            ):
                final_day += 1  # Additional day if fever is the last symptom

        elif category == 2:
            fever_ever = False
            if not df_id[df_id["symp_fever"] == 1].empty:
                final_day_fever = (
                    df_id[df_id["symp_fever"] == 1][
                        "days_since_symp_onset"
                    ].max()
                    + 1
                )
                fever_ever = True

            if not df_id[df_id["symp_body_ache"] == 1].empty:
                final_day_ache = (
                    df_id[df_id["symp_body_ache"] == 1][
                        "days_since_symp_onset"
                    ].max()
                    + 1
                )

            final_day = max(final_day_fever, final_day_ache)
            censor_final_day = final_day
            if (fever_ever) & (
                df_id[df_id["symp_fever"] == 1]["days_since_symp_onset"].max()
                == (final_day - 1)
            ):
                final_day += 1  # Additional day if fever is the last symptom

        elif category in [3, 4]:
            # Overall approach: We are going to trace through the symp count trajectory. First, we're going to find
            # the index of the first instance of the maximum symptom count (max_idx). We are going to then find
            # the corresponding days_since_symp_onset for that max_idx, and assign it to peak_symptoms_day. We are then
            # going to try stepping forward (max_idx +=1). If there is a row there, and if the value of symptom count
            # is the same (meaning that we are in a plateau), then we update peak_symptoms_day to reflect max_idx,
            # and then move forward in the 'while' loop until the plateau drops off

            # Calculate the peak symptoms day
            if not df_id["symp_count"].empty:
                # Find the index of the first occurrance of the maximum symp_count value
                max_idx = df_id["symp_count"].idxmax()
                max_symp_count = df_id["symp_count"].max()

                # step forward until the first symp_count decrease
                peak_symptoms_day = df_id.loc[max_idx, "days_since_symp_onset"]
                if df_id.shape[0] > 1:
                    df_id_afterpeak = df_id.loc[max_idx:, :]

                    df_postplateau = df_id_afterpeak.loc[
                        df_id_afterpeak["symp_count"] != max_symp_count, :
                    ]
                    if not df_postplateau.empty:
                        decrease_symp_indx = df_postplateau.index[0] - 1
                        peak_symptoms_day = df_id.loc[
                            decrease_symp_indx, "days_since_symp_onset"
                        ]

            else:
                print("ERROR in cat 3&4 calculation")
                break

            final_day = peak_symptoms_day + 1  # One day after peak symptoms
            censor_final_day = final_day

        # Assign final day to symp_duration, respecting censoring
        try:
            df.loc[df["id"] == id_, "symp_duration"] = min(
                final_day, last_day_of_data
            )

        except Exception:
            print(category)
            print(id_)
            break

        if (
            (category in [1, 2])
            & (censor_final_day < final_day)
            & (censor_final_day == last_day_of_data)
        ):  # special circumstance
            df.loc[df["id"] == id_, "symp_duration"] = final_day

        # Determine if the duration is censored
        df.loc[df["id"] == id_, "symp_duration_censored"] = int(
            censor_final_day > last_day_of_data
        )

    return df


def create_individual_patient_plot(
    df,
    yaxis1_var,
    time_var,
    title,
    additional_info,
    yaxis2_var=None,
    truefalse_vars=None,
):
    """
    This function takes in a timeseries symptom/testing dataframe for an individual patient and info on
    variables/descriptions and returns a plot

    df: dataframe for individual patient
    yaxis1_var: dictionary {'varname':'label'}
    time_var: string (variable name of time variable)
    title: string (plot title)
    additional_info: string (any other info to display in a box)
    yaxis2_var: (otional): dictionary {'varname':'label'}
    truefalse_vars (optional): dictionary {'varname1':'label1','varname2':'label2', etc.}

    returns: plot
    """
    fig, ax1 = plt.subplots()

    # Plotting yaxis1_var
    y1 = yaxis1_var[list(yaxis1_var.keys())[0]]
    ax1.plot(
        df[time_var],
        df[list(yaxis1_var.keys())[0]],
        "b-",
        label=y1,
        marker="o",
        linestyle="-",
    )
    ax1.set_xlabel(time_var)
    ax1.set_ylabel(y1, color="b")
    ax1.tick_params(axis="y", labelcolor="b")

    # If yaxis2_var is provided, plot on a secondary y-axis
    if yaxis2_var:
        ax2 = ax1.twinx()
        y2 = yaxis2_var[list(yaxis2_var.keys())[0]]
        ax2.plot(
            df[time_var],
            df[list(yaxis2_var.keys())[0]],
            "r-",
            label=y2,
            marker="o",
            linestyle="-",
        )
        ax2.set_ylabel(y2, color="r")
        ax2.tick_params(axis="y", labelcolor="r")

    # Plotting true/false variables if provided
    if truefalse_vars:
        for var, label in truefalse_vars.items():
            true_points = df[df[var] == 1]
            false_points = df[df[var] == 0]
            ax1.plot(
                true_points[time_var],
                true_points[list(yaxis1_var.keys())[0]],
                "g^",
                label=label + " True",
            )
            ax1.plot(
                false_points[time_var],
                false_points[list(yaxis1_var.keys())[0]],
                "rs",
                label=label + " False",
            )

    # Title and legend
    fig.suptitle(title)
    fig.legend(loc="upper left", bbox_to_anchor=(1.05, 1), borderaxespad=0.0)

    # Additional information box
    textstr = additional_info
    props = dict(boxstyle="round", facecolor="wheat", alpha=0.5)
    ax1.text(
        1.05,
        0.5,
        textstr,
        transform=ax1.transAxes,
        fontsize=10,
        verticalalignment="top",
        bbox=props,
    )

    # Layout adjustment
    fig.tight_layout(
        rect=[0, 0, 0.75, 1]
    )  # Adjust depending on the length of additional_info
    plt.show()
