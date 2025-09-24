import numpy as np
import pandas as pd
from utils import general_utils


def create_rvtn_dataset(file_str_data: str, file_str_patients: str):
    """
    Description: this function ingests the RVTN (household) dataset (.xlsx),
    formats it, creates basic variables, and outputs a dataframe
    Input: file location (string) for main dataset and patient attributes
    Output: pandas dataframe
    """

    # -------------------------------------------------------------------------
    # Load Data
    # -------------------------------------------------------------------------

    data_str = file_str_data
    patient_str = file_str_patients

    df = pd.read_csv(data_str)
    df_patients = pd.read_csv(patient_str)

    # -------------------------------------------------------------------------
    # Add variables
    # -------------------------------------------------------------------------

    # symptom_count
    symptoms = [
        "fever",
        "cough",
        "sore_throat",
        "runny_nose",
        "congestion",
        "fatigue",
        "wheezing",
        "sob",
        "chest_pain",
        "smell",
        "headache",
        "abdominal_pain",
        "diarrhea",
        "vomiting",
        "aches",
    ]

    # inherent asks about congestion or runny nose as a single variable
    # rvtn asks about congestion and runny nose as separate variables
    # so we need to create a new variable for congestion/runny nose
    # call this variable "nose" because that is what inherent calls it
    df["nose"] = np.where(
        np.logical_or((df["congestion"] == 1), (df["runny_nose"] == 1)), 1, 0
    )

    # checks
    assert all(
        np.logical_or(
            df[df["nose"] == 1]["runny_nose"] == 1,
            df[df["nose"] == 1]["congestion"] == 1,
        )
    )
    assert all(
        np.logical_and(
            df[df["nose"] == 0]["runny_nose"] != 1,
            df[df["nose"] == 0]["congestion"] != 1,
        )
    )

    # remove congestion and runny nose from symptoms
    symptoms.remove("congestion")
    symptoms.remove("runny_nose")
    # instead add this new new variable
    symptoms.append("nose")

    df["symptom_count"] = df[symptoms].apply(lambda x: (x == 1).sum(), axis=1)

    # symp flag
    df["symp_flag"] = np.where(df["symptom_count"] > 0, 1, df["symptom_count"])

    # add asymptomatic variable to main dataset
    asymptomatic_ids = df_patients[df_patients["cdc_ever_symptoms"] == 0][
        "cdc_studyid"
    ].tolist()
    df["cdc_ever_symptoms"] = (
        df["cdc_studyid"].isin(asymptomatic_ids).apply(lambda x: 0 if x else 1)
    )

    # add time period of infection
    df["cdc_period"] = df["cdc_studyid"].map(
        df_patients.set_index("cdc_studyid")["cdc_period"].to_dict()
    )

    # add cdc_agegrp
    df["cdc_agegrp"] = df["cdc_studyid"].map(
        df_patients.set_index("cdc_studyid")["cdc_agegrp"].to_dict()
    )

    # add covid vaccine doses in the last twelve months
    df["days_lastcovdose2enroll"] = df["cdc_studyid"].map(
        df_patients.set_index("cdc_studyid")[
            "days_lastcovdose2enroll"
        ].to_dict()
    )

    # add cdc_incident
    df["cdc_incident"] = df["cdc_studyid"].map(
        df_patients.set_index("cdc_studyid")["cdc_incident"].to_dict()
    )

    # add participant type variable to main dataset
    index_ids = df_patients[df_patients["participant_type"] == 1][
        "cdc_studyid"
    ].tolist()
    df["participant_type"] = (
        df["cdc_studyid"].isin(index_ids).apply(lambda x: 1 if x else 2)
    )

    ## DROP DUPLICATES
    df.drop_duplicates(
        subset=["cdc_studyid", "days_since_onset"], keep="last", inplace=True
    )

    # Sort dataframe
    df = df.sort_values(
        by=["cdc_studyid", "days_since_symptom_onset"]
    ).reset_index(drop=True)

    return df


def process_rvtn_dataset_for_STAN(df):
    # -------------------------------------------------------------------------
    # Add/refactor variables
    # -------------------------------------------------------------------------

    # rename
    rename_dict = {
        "cdc_studyid": "id",
        "covqpcr_load": "logVL",
        "cdc_covpos_sample": "pcr",
        "cdc_antigen_pos": "antigen",
        "cdc_culture_pos": "culture",
        "days_since_symptom_onset": "days_since_symp_onset",
        "cdc_ever_symptoms": "symptomatic_ever",
        "cdc_incident": "incident",
        "symptom_count": "symp_count",
    }
    df.rename(columns=rename_dict, inplace=True)

    # fix symptomatic_ever variable
    max_symp_counts = df.groupby("id")["symp_count"].max().to_dict()
    df["symptomatic_ever"] = df["id"].map(
        lambda x: 1 if max_symp_counts.get(x, 0) > 0 else 0
    )

    # add in onset for asymptomatic patients
    df["days_since_symp_onset"] = np.where(
        df["symptomatic_ever"] == 0,
        df["days_since_onset"],
        df["days_since_symp_onset"],
    )

    # antigen refactor
    df["antigen"] = np.where(df["antigen"] == 99, 100, df["antigen"])

    # culture refactor (99s not actually observed in data)
    df["culture"] = np.where(df["culture"] == 99, 100, df["culture"])

    # pcr refactor
    df["pcr"] = np.where(df["pcr"] == 99, 100, df["pcr"])

    # logVL refactor
    conditions = [
        (df["pcr"] == 0),
        (df["pcr"] == 1) & (pd.isnull(df["logVL"])),
        (df["pcr"] == 1) & (pd.notnull(df["logVL"])),
        (df["pcr"] != 1) & (df["pcr"] != 0),
    ]
    values = [-10, 0, df["logVL"], 100]
    df["logVL"] = np.select(conditions, values, default=None)

    # age
    age_cat_dict = {
        "18-49": 1,
        "50-64": 2,
        "5-11": 0,
        "0-4": 0,
        "65+": 3,
        "12-17": 0,
    }
    df["age_cat"] = df["cdc_agegrp"].map(age_cat_dict)

    # adjust cdc period (time frame for when infections happened)
    # for xbb to include later times to be consistent with INHERENT
    df["cdc_period"] = df["cdc_period"].replace(
        "XBB etc: Jan 15-May 1, '23", "XBB etc: Jan 15-Oct 31, '23"
    )

    # convert days to month for time since last covid dose
    # on average, 30.437 days to a month
    # round rather than floor because that's what INHERENT did too
    df["months_since_last_covid_dose"] = -np.round(
        df["days_lastcovdose2enroll"] / 30.437
    )

    # study_id
    df["study_id"] = df["id"]

    # Refactor symptom variables
    symp_new_vars_dict = {
        "fever": "symp_fever",
        "cough": "symp_cough",
        "sob": "symp_short_breath",
        "fatigue": "symp_fatigue",
        "aches": "symp_body_ache",
        "headache": "symp_headache",
        "sore_throat": "symp_sore_throat",
        "nose": "symp_nose",
        "vomiting": "symp_vomiting",
        "chest_pain": "symp_chest_pain",
        "smell": "symp_smell_loss",
        "wheezing": "symp_wheezing",
    }
    df.rename(columns=symp_new_vars_dict, inplace=True)
    for symptom in symp_new_vars_dict.values():
        conditions = [df[symptom] == 1, df[symptom] == 0]
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

    # keep only desired variables
    keep = [
        "id",
        "study_id",
        "days_since_symp_onset",
        "logVL",
        "pcr",
        "antigen",
        "culture",
        "symp_flag",
        "symptomatic_ever",
        "participant_type",
        "symp_duration",
        "symp_duration_censored",
        "age_cat",
        "months_since_last_covid_dose",
        "cdc_period",
        "symp_count",
        "incident",
        "logVL_quantifiable_ever",
        "symp_type_cat",
    ] + list(symp_new_vars_dict.values())
    df = df[keep]

    # fill in missing with -1
    df = df.fillna(-1)

    return df
