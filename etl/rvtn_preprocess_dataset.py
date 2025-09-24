import os
import warnings
from argparse import ArgumentParser

import yaml
from dotenv import load_dotenv
from utils import rvtn_utils

# Set production mode from --args
allowed_options = ["production", "development"]
##output_mode, printme = general_utils.get_args(allowed_options)
##print(printme)

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

# Set up input and output paths
load_dotenv()
input_dir = os.path.expanduser(os.getenv("DATA_INPUT_DIR"))
output_dir = os.path.expanduser(os.getenv("MAIN_OUTPUT_DIR"))

with open("products_config.yml", "r") as file:
    products_config = yaml.safe_load(file)

output_folder = os.path.join(
    output_dir, output_mode, products_config["preprocessed_data"]["directory"]
)

lab_data_file = os.path.join(
    input_dir, products_config["input_data"]["RVTN_lab"]
)  # lab dataset (raw)
main_data_file = os.path.join(
    input_dir, products_config["input_data"]["RVTN"]
)  # main dataset (raw)
file_out_main = os.path.join(
    output_folder, products_config["preprocessed_data"]["RVTN"]
)  # main dataset (preprocessed)
file_out_patient_level = os.path.join(
    output_folder, products_config["preprocessed_data"]["RVTN_patient"]
)  # patient-level dataset (preprocessed)

# Ensure the output directory exists
os.makedirs(os.path.expanduser(output_folder), exist_ok=True)

# Suppress warnings
warnings.filterwarnings("ignore")

# Load data and process
df = rvtn_utils.create_rvtn_dataset(
    lab_data_file,
    main_data_file,
)

df = rvtn_utils.process_rvtn_dataset_for_STAN(df)

df_patient = df_patient = df.groupby("id", as_index=False).agg(
    {
        "study_id": "first",
        "symptomatic_ever": "first",
        "participant_type": "first",
        "symp_duration": "first",
        "symp_duration_censored": "first",
        "age_cat": "first",
        "incident": "first",
        "logVL_quantifiable_ever": "first",
        "symp_fever": "max",
        "symp_cough": "max",
        "symp_short_breath": "max",
        "symp_fatigue": "max",
        "symp_body_ache": "max",
        "symp_headache": "max",
        "symp_sore_throat": "max",
        "symp_nose": "max",
        "symp_vomiting": "max",
        "symp_chest_pain": "max",
        "symp_smell_loss": "max",
        "symp_wheezing": "max",
        "symp_type_cat": "max",
    }
)

# print(df.head(2))
print("######## RVTN ETL OUTPUT #########")
print(f"Output to {file_out_main}")
print(f"Output to {file_out_patient_level}")
print("")


# Output data
df.to_csv(
    file_out_main,
    index=False,
)

df_patient.to_csv(
    file_out_patient_level,
    index=False,
)
