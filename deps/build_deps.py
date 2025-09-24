import re
from argparse import ArgumentParser

import yaml


def parse_yaml(file_in):
    yaml_fh = open(file_in, "r")
    yaml_values = yaml.safe_load(yaml_fh)

    return yaml_values


def process_makefile(file_in, config):
    fh = open(file_in, "r")
    makefile_str = "".join(fh.readlines())
    fh.close()
    stan_processed_file = "stan_ready_preprocessed_data_output"

    replace_dict = {
        "PREPROCESS_DIR": config["preprocessed_data"]["directory"],
        "RVTN_DEPS": "%s %s"
        % (config["input_data"]["RVTN"], config["input_data"]["RVTN_lab"]),
        "INHERENT_DEPS": config["input_data"]["INHERENT"].replace(" ", "\\ "),
        "RVTN_FILE": config["preprocessed_data"]["RVTN"],
        "RVTN_PATIENT_FILE": config["preprocessed_data"]["RVTN_patient"],
        "INHERENT_FILE": config["preprocessed_data"]["INHERENT"],
        "STAN_PROCESSED_FILE": config["vl_model"][stan_processed_file],
        "STAN_POST_DIR": config["vl_model"]["stan_posterior_directory"],
        "STAN_FIT_DIR": config["vl_model"]["stan_fit_directory"],
        "STAN_POST_FILE": config["vl_model"]["posterior_individual"],
        "STAN_SAMPLE_FILE": config["vl_model"]["posterior_sampled"],
        "STAN_FIT_FILE": config["vl_model"]["run_name"],
        "SCENARIOS_DIR": config["scenarios"]["directory"],
        "GUIDANCE_FILE": config["scenarios"]["compare_guidance_df"],
        "OVERALL_GUIDANCE_FILE": config["scenarios"][
            "overall_compare_guidance_df"
        ],
        "PRODUCTS_PATH": config["products"]["directory"],
    }
    for k, v in replace_dict.items():
        tmp_make = re.sub(
            re.compile(r"^(%s :=).*$" % k, re.MULTILINE),
            "\\1 %s" % v,
            makefile_str,
        )
        makefile_str = tmp_make
    return makefile_str


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument(
        "-m", "--makefile", dest="make_file", default="deps/Makefile.in"
    )
    parser.add_argument(
        "-o", "--output", dest="outfile", default="Dependencies.mk"
    )
    parser.add_argument(
        "-y", "--yml", dest="yaml_file", default="products_config.yml"
    )

    args = parser.parse_args()
    products_dict = parse_yaml(file_in=args.yaml_file)

    makefile_out = process_makefile(
        file_in=args.make_file, config=products_dict
    )
    fh = open(args.outfile, "w")
    fh.write(makefile_out)
    fh.close()
