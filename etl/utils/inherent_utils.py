import numpy as np
import pandas as pd
from utils import general_utils


def create_inherent_dataset(file_str: str):
    """
    Description: this function ingests the INHERENT (nursing home) response dataset (.xlsx),
    formats it, creates basic variables, and outputs a dataframe
    Input: file location (string)
    Output: pandas dataframe
    """

    # -------------------------------------------------------------------------
    # Load Data
    # -------------------------------------------------------------------------

    # Fill with your data location
    data_str = file_str

    # Testing Data
    df = pd.read_excel(data_str, sheet_name="Case list test results")

    # Staff Participant Info
    df_demo_staff = pd.read_excel(data_str, sheet_name="Consented staff list")
    df_demo_staff["type"] = "Staff"

    # Resident Participant Info
    df_demo_residents = pd.read_excel(
        data_str, sheet_name="NH resident enrollees"
    )
    df_demo_residents["type"] = "Resident"

    # Residents had their vaccination date stored as the actual date, not how long ago
    # But we still want vacc_ago for staff and not to drop it
    df_demo_residents["vacc_ago"] = np.nan

    df_demo = pd.concat([df_demo_staff, df_demo_residents], axis=0)

    # -------------------------------------------------------------------------
    # Get Symptomatic / Asymptomatic / Unknown Flag and Tier 2 Symptom Counts over Time
    # -------------------------------------------------------------------------
    # Symptomatic - 1
    # Asymptomatic - 0
    # Blank or Unknown - nan

    symp_indicators = np.array(
        [
            "Symptomatic, worsening",
            "Symptomatic, no change",
            "Symptomatic, improving",
        ]
    )
    symp_unknown = np.array(
        [
            "Symptom status unknown",
            "Unknown",
            "Valid skip",
            "NaN",
            "N/A, not Tier 1",
        ]
    )
    symp_cols = [
        "symp_t2_fever",
        "symp_t2_cough",
        "symp_t2_breath",
        "symp_t2_fatigue",
        "symp_t2_ache",
        "symp_t2_head",
        "symp_t2_taste",
        "symp_t2_throat",
        "symp_t2_nose",
        "symp_t2_nausea",
        "symp_t2_diarrhea",
        "symp_t2_appetite",
    ]

    sympT1_present = []
    sympT2_present = []
    sympT2_count = []
    for pat in df.study_id.unique():
        pat_df = df[df.study_id == pat]
        # get symptom data for Tier 1 and Tier 2
        tier1_traj = pat_df.symptom_tier1.fillna("NaN").values
        tier2_traj = pat_df[symp_cols].fillna("NaN").values

        for j in range(len(tier1_traj)):
            if str(tier1_traj[j]) in symp_indicators:
                sympT1_present.append(1)
            elif str(tier1_traj[j]) in symp_unknown:
                sympT1_present.append(np.nan)
            else:
                sympT1_present.append(0)

            if "Yes" in tier2_traj[j]:  # if any tier 2 symptom, symptomatic
                sympT2_present.append(1)
                sympT2_count.append(np.sum(tier2_traj[j] == "Yes"))
            elif (
                "No" in tier2_traj[j]
            ):  # if there are no symptoms and at least one no symptom survey was completed
                sympT2_present.append(0)
                sympT2_count.append(0)
            else:  # No symptom data present
                sympT2_present.append(np.nan)
                sympT2_count.append(np.nan)

    # Get flags to identify which symptom data to look at
    # We do this because if a patient is Tier 1 their Tier 2 symptoms will be NaN and vice versa
    tier2_flag = 1 * np.array([df.tier == 2, df.tier == 3]).any(axis=0)
    tier1_flag = 1 * (df.tier == 1).values
    index_flag = 1 * (df.tier == 0).values

    sympT1_present = np.array(sympT1_present)
    sympT2_present = np.array(sympT2_present)

    symp_T1 = np.zeros(sympT1_present.shape)
    symp_T1[tier1_flag == 1] = sympT1_present[tier1_flag == 1]

    symp_T2 = np.zeros(sympT2_present.shape)
    symp_T2[tier2_flag == 1] = sympT2_present[tier2_flag == 1]

    symp_T2[index_flag == 1] = sympT2_present[index_flag == 1]

    symp_present = symp_T1 + symp_T2

    # Adds a symptom flag for any symptom and specifically T1 and T2 symptoms
    # Also gives a count for # of tier 2 symtpoms
    df["symp_flag"] = symp_present
    df["sympT1_flag"] = sympT1_present
    df["sympT2_flag"] = sympT2_present
    df["sympT2_count"] = sympT2_count

    # -------------------------------------------------------------------------
    # Add key data to Case List Data Frame
    # -------------------------------------------------------------------------
    for col in [
        "age",
        "vacc_ago",
        "enr_date",
    ]:
        hold = []
        for pat in df.study_id.unique():
            pat_df = df[df.study_id == pat]
            pat_demo = df_demo[df_demo.study_id == pat]

            hold.append([pat_demo[col].values[0]] * len(pat_df))

        df[col] = np.concatenate(hold)

    # Modify age to be numeric
    df.age = pd.to_numeric(df.age.replace("90+", 90), errors="coerce")

    # -------------------------------------------------------------------------
    # Add time variables wrt certain events; these make plotting easier
    # -------------------------------------------------------------------------
    # delta      - time since first data collection
    # delta_peak - time since peak
    # delta_pos  - time since first positive
    # delta_symp - time since symptom onset
    # delta_ant  - time since antigen onset
    # delta_cult - time since culture onset
    # delta_pcr  - time since pcr onset

    # Add column that is only numeric Ct values
    df["Ct_numeric"] = pd.to_numeric(df.Ct, errors="coerce")
    df.colldate = pd.to_datetime(df.colldate)

    (
        delta,
        delta_peak,
        delta_pos,
        delta_symp,
        delta_ant,
        delta_cult,
        delta_pcr,
    ) = (
        [],
        [],
        [],
        [],
        [],
        [],
        [],
    )
    for pat in df.study_id.unique():
        pat_df = df[df.study_id == pat]

        # Gets time in days from first interaction
        pat_days = (
            pat_df.colldate.values - pat_df.colldate.values.min()
        ) / np.timedelta64(1, "D")
        # delta
        delta.append(pat_days)
        # delta_peak
        if pat_df.Ct_numeric.notna().any():
            delta_peak.append(
                pat_days - pat_days[np.nanargmin(pat_df.Ct_numeric)]
            )
        else:
            delta_peak.append(np.nan * np.zeros(len(pat_df)))
        # delta_pos
        first_pos_idx = np.concatenate(
            [
                np.where(pat_df.pcr == "Positive")[0],
                np.where(pat_df.antigen == "Positive")[0],
                np.where(pat_df.vculture == "Positive")[0],
            ]
        ).min()
        delta_pos.append(pat_days - pat_days[first_pos_idx])
        # delta_symp
        symp_idx = np.where(pat_df.symp_flag == 1)[0]
        if len(symp_idx) > 0:
            delta_symp.append(pat_days - pat_days[symp_idx[0]])
        else:
            delta_symp.append(np.nan * np.zeros(len(pat_df)))
        # delta_ant
        ant_idx = np.where(pat_df.antigen == "Positive")[0]
        if len(ant_idx) > 0:
            delta_ant.append(pat_days - pat_days[ant_idx[0]])
        else:
            delta_ant.append(np.nan * np.zeros(len(pat_df)))
        # delta_cult
        cult_idx = np.where(pat_df.vculture == "Positive")[0]
        if len(cult_idx) > 0:
            delta_cult.append(pat_days - pat_days[cult_idx[0]])
        else:
            delta_cult.append(np.nan * np.zeros(len(pat_df)))
        # pcr_ant
        pcr_idx = np.where(pat_df.pcr == "Positive")[0]
        if len(pcr_idx) > 0:
            delta_pcr.append(pat_days - pat_days[pcr_idx[0]])
        else:
            delta_pcr.append(np.nan * np.zeros(len(pat_df)))

    df["delta"] = np.concatenate(delta)
    df["delta_peak"] = np.concatenate(delta_peak)
    df["delta_pos"] = np.concatenate(delta_pos)
    df["delta_symp"] = np.concatenate(delta_symp)
    df["delta_ant"] = np.concatenate(delta_ant)
    df["delta_cult"] = np.concatenate(delta_cult)
    df["delta_pcr"] = np.concatenate(delta_pcr)

    df = df.sort_values(by=["study_id", "colldate"]).reset_index(drop=True)

    ## DROP DUPLICATES
    df.drop_duplicates(
        subset=["study_id", "colldate"], keep="last", inplace=True
    )

    return df


def process_inherent_dataset_for_STAN(df, parameters: dict):
    # -------------------------------------------------------------------------
    # Add/refactor variables
    # -------------------------------------------------------------------------

    # Create lists of symptomatic/asymptomatic/no symptom data individuals
    grouped = df.groupby("study_id", as_index=False).agg(
        {"symp_flag": "max", "sympT2_count": "max"}
    )
    symptomatic_list = grouped[grouped["symp_flag"] == 1]["study_id"].to_list()
    asymptomatic_list = grouped[grouped["symp_flag"] == 0][
        "study_id"
    ].to_list()
    # no_symptom_data_list = grouped[pd.isnull(grouped['symp_flag'])]['study_id'].to_list()

    # Add id variable
    id_map = {
        study_id: i + 1 for i, study_id in enumerate(df["study_id"].unique())
    }
    df["id"] = df["study_id"].map(id_map)

    # Add days since symptom onset
    df["days_since_symp_onset"] = df["delta_symp"]

    # Add age_cat
    age_under_18_list = list(df[df["age"] < 18].study_id.unique())
    age_over_18_under_50_list = list(
        df[np.logical_and(df["age"] >= 18, df["age"] < 50)].study_id.unique()
    )
    age_over_50_under_65_list = list(
        df[np.logical_and(df["age"] >= 50, df["age"] < 65)].study_id.unique()
    )
    age_65_plus_list = list(df[df["age"] >= 65].study_id.unique())
    conditions = [
        (df["study_id"].isin(age_under_18_list)),
        (df["study_id"].isin(age_over_18_under_50_list)),
        (df["study_id"].isin(age_over_50_under_65_list)),
        (df["study_id"].isin(age_65_plus_list)),
    ]
    values = [0, 1, 2, 3]
    df["age_cat"] = np.select(conditions, values, default=np.nan)

    # all INHERENT are after May 1
    # the latest INHERENT person is October 26th
    # XBB was still the dominant variant at that time
    # so just call this XBB era the year of 2023 basically
    df["cdc_period"] = "XBB etc: Jan 15-Oct 31, '23"  # match RVTN formatting

    # Add symptomatic_ever
    conditions = [
        (df["study_id"].isin(symptomatic_list)),
        (df["study_id"].isin(asymptomatic_list)),
    ]
    values = [1, 0]
    df["symptomatic_ever"] = np.select(conditions, values, default=np.nan)

    # Change antigen variable
    df["antigen"] = df["antigen"].replace(
        {"Positive": 1, "Not done": 100, "Inconclusive": 100, "Negative": 0}
    )

    # Change culture variable
    df["vculture"] = df["vculture"].replace(
        {
            "Positive": 1,
            "Will not be done": 100,
            "Not yet ordered": 100,
            "Inconclusive": 100,
            "Negative": 0,
        }
    )

    # Change PCR variable
    df["pcr"] = df["pcr"].replace(
        {
            "Positive": 1,
            "Not collected": 100,
            "Not tested": 100,
            "Not yet reported": 100,
            "Negative": 0,
        }
    )

    # change logVL variable
    conditions = [
        (df["pcr"] == 0),
        np.logical_and(
            (df["pcr"] == 1),
            (
                np.logical_or(
                    df["VL_status"] == "Less than 10,000 copies/mL",
                    df["VL_status"] == "VL run, not  detected",
                )
            ),
        ),
        np.logical_and(
            (df["pcr"] == 1),
            (df["VL_status"] == "More than 1 billion copies/mL"),
        ),
        np.logical_and((df["pcr"] == 1), (pd.notnull(df["logVL"]))),
        np.logical_or(
            np.logical_and((df["pcr"] != 1), (df["pcr"] != 0)),
            pd.isnull(df["logVL"]),
        ),
    ]
    # when using np.select, if multiple conditions are met, assigns value based on the first one
    # so just sticking with logVL available must be last condition
    values = [
        -10,
        0,
        parameters["global_data_markers"]["loq_upper"],
        df["logVL"],
        100,
    ]
    df["logVL"] = np.select(conditions, values, default=None)

    # check that the assignment worked as expected
    # all available logVL values should be between loq_lower and loq_upper
    pcr_positive_logVLs = df["logVL"][df["pcr"] == 1]
    # can enforce equality because we *know* that there are always LOQ
    # (both lower bound LOQ and upper bound LOQ) values in dataset
    assert (
        pcr_positive_logVLs.min()
        == parameters["global_data_markers"]["loq_lower"]
    )
    assert (
        pcr_positive_logVLs.max()
        == parameters["global_data_markers"]["loq_upper"]
        or
        # to account for null VL values
        pcr_positive_logVLs.max()
        == parameters["global_data_markers"]["skipped_test"]
    )
    available_logVLs = pcr_positive_logVLs[
        np.logical_and(
            np.logical_and(
                df["logVL"] != parameters["global_data_markers"]["loq_lower"],
                df["logVL"] != parameters["global_data_markers"]["loq_upper"],
            ),
            # to account for pcr positive but logVL null
            df["logVL"] != parameters["global_data_markers"]["skipped_test"],
        )
    ]
    # now we check that our bound for available logVLs and LOQ is correct
    assert available_logVLs.min() > parameters["INHERENT"]["loq_lower"]
    assert available_logVLs.max() < parameters["INHERENT"]["loq_upper"]

    # there is another way to go about this test -- that all VL_status "Within detectable range"
    # are within what we believe is the detectable range
    logVLs_within_detectable_range = df["logVL"][
        df["VL_status"] == "Within detectable range"
    ]
    assert (
        logVLs_within_detectable_range.min()
        > parameters["INHERENT"]["loq_lower"]
    )
    assert (
        logVLs_within_detectable_range.max()
        < parameters["INHERENT"]["loq_upper"]
    )

    # Refactor symptom variables
    symp_new_vars_dict = {
        "symp_t2_fever": "symp_fever",
        "symp_t2_cough": "symp_cough",
        "symp_t2_breath": "symp_short_breath",
        "symp_t2_fatigue": "symp_fatigue",
        "symp_t2_ache": "symp_body_ache",
        "symp_t2_head": "symp_headache",
        "symp_t2_taste": "symp_taste_loss",
        "symp_t2_throat": "symp_sore_throat",
        "symp_t2_nose": "symp_nose",
        "symp_t2_nausea": "symp_nausea",
        "symp_t2_diarrhea": "symp_diarrhea",
        "symp_t2_appetite": "symp_appetite",
        "vacc_ago": "months_since_last_covid_dose",
    }
    df.rename(columns=symp_new_vars_dict, inplace=True)

    # Refactor culture variable
    cult_new_vars_dict = {
        "vculture": "culture",
    }
    df.rename(columns=cult_new_vars_dict, inplace=True)

    rvtn_symp_vars = [
        "symp_smell_loss",
        "symp_congestion",
        "symp_wheezing",
        "symp_chest_pain",
    ]
    for var in rvtn_symp_vars:
        df[var] = np.nan

    df["symp_count"] = df["sympT2_count"]
    for symptom in symp_new_vars_dict.values():
        conditions = [df[symptom] == "Yes", df[symptom] == "No"]
        values = [1, 0]
        df[symptom] = np.select(conditions, values, default=np.nan)

    # add symptom duration columns
    df = general_utils.add_symp_duration_columns(df)

    # add logVL_quantifiable_ever
    cond_log_VL_quantified_set = set(
        df[(df["logVL"] > 1.5) & (df["logVL"] != 100)]["id"].unique()
    )
    df["logVL_quantifiable_ever"] = np.where(
        df["id"].isin(cond_log_VL_quantified_set), 1, 0
    )

    conditions = [df["type"] == "Staff", df["type"] == "Resident"]
    # The downstream pipeline expects `participant_type` to be a 3 or 4
    # to mark staff or residents -- this is our shorthand to make filtering
    # pre Stan easier
    values = [3, 4]
    df["participant_type"] = np.select(conditions, values, default=np.nan)

    # keep only desired variables
    keep = [
        "id",
        "study_id",
        "participant_type",
        "days_since_symp_onset",
        "logVL",
        "pcr",
        "antigen",
        "culture",
        "Ct_numeric",
        "symp_flag",
        "age",
        "age_cat",
        "cdc_period",
        "symptomatic_ever",
        "symp_count",
        "symp_duration",
        "symp_duration_censored",
        "logVL_quantifiable_ever",
        "symp_type_cat",
    ] + list(symp_new_vars_dict.values())
    df = df[keep]

    # fill in missing with -1
    df = df.fillna(-1)

    return df
