import os
import warnings
from argparse import ArgumentParser

import yaml
from dotenv import load_dotenv
from utils import inherent_utils

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

with open(os.path.join("analysis", "parameters.yml"), "r") as file:
    parameters = yaml.safe_load(file)

output_folder = os.path.join(
    output_dir, output_mode, products_config["preprocessed_data"]["directory"]
)

file_in = os.path.join(input_dir, products_config["input_data"]["INHERENT"])
file_out = os.path.join(
    output_folder, products_config["preprocessed_data"]["INHERENT"]
)

# Ensure the output folder exists
os.makedirs(os.path.expanduser(output_folder), exist_ok=True)

# Suppress warnings
warnings.filterwarnings("ignore")

# Load data and process
df = inherent_utils.create_inherent_dataset(
    file_in,
)
df = inherent_utils.process_inherent_dataset_for_STAN(df, parameters)

# Output data
df.to_csv(
    file_out,
    index=False,
)

# print(df.head(2))
print("######## INHERENT ETL OUTPUT #########")
print(f"Output to {file_out}")
print("")
