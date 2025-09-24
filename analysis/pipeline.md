# Pipeline

## Overview

This document details the steps to run all the files in this project (referred to as the "model pipeline").

Before proceeding, ensure that the set-up steps in `README.md` have been completed.

There are two potential workflows:
1. Run the Makefile, which will automatically run the scripts described below.
2. Run the model pipeline script by script in the order described below.

Regardless of the workflow, all scripts are intended to be run from project root. This means that your current working directory must be the file in which the project is stored (e.g., `covid-isolation-pub-2025`).

By default, the code will be run with the parameters in `analysis/scripts/parameters.yml` and in development mode.
The fraction of the dataset used to fit the model in development mode can be set with the `development_dataset_subset`
parameter, the MCMC settings/number of samples from the posterior can be set from the `stan_fit_parameters`,
and the number of synthetic individuals to simulate can be set with the `samples_per_symp_cat` parameter.

Most parameters except for those labeled as "prior", "scenario", or "sensitivity" parameters reflect
annotations of the data/predefined quantities about the data. To vary the parameters that impact model
results, the "prior" parameters define the priors used in the Stan model of disease natural history,
the "scenario" parameters define the parameters used in the main scenario analysis (i.e., simulated
interventions), and the "sensitivity" parameters specify the values of parameters that are
sensitivity tested.

## Run the pipeline via the Makefile

To run all steps in the Makefile, run the following line of code:

1. `make`

If running `make` is successful, the output files will be found in the path specified in `.env`. (For information on making a `.env` file, see `README.md`)

### Alternate runs of make

The above code runs the pipeline using development mode and will place outputs in a `development` folder. However, there are alternatives that a user can run.

#### Production

To run a production run, use `make PROD=1`. This will place outputs in a `production` folder.

#### Running pipeline substeps with make

**Note:** Each of these substeps will run both the specified files and any dependency files that have been updated more recently than the specified files.

1. To run the extract-transform-load (ETL) data preprocessing steps used in the pipeline (INHERENT + RVTN), run `make etl`.
    - `etl/inherent_preprocess_dataset.py`
    - `etl/rvtn_preprocess_dataset.py`
2. To run the model that fits viral load curves (in the file `triangle_vl_symp_hazard.stan`), run `make stan`.
    - `analysis/scripts/01_prepare_data_for_stan.R`
    - `analysis/scripts/02_fit_viral_kinetics_by_symp_cat.R`
    - `analysis/scripts/03_extract_parameters_from_stan_fit.R`
    - `analysis/scripts/04_subsample_extracted_parameters.R`
    - Note that this will also run ETL preprocessing steps if those files have changed.
3. Finally, to run different scenarios based on parameters listed in `analysis/parameters.yml`, run `make scenarios`.
    - `analysis/scripts/05_run_ind_scenarios.R`
    - `analysis/scripts/06_create_overall_category.R`
    - Note that this will also run the model that fits viral load curves and ETL preprocessing steps if those files have changed or have been affected by changes to files in [`isolation/R`](/isolation/R/).

#### Products

To specifically create the figures and tables used in the manuscript, run `make products`. By default this will run in development mode. To run in production mode run `make products PROD=1`.
- This will run all scripts in `analysis/scripts/products` and create a list of all output products in `deps/products_list`.
- Note that some figures/tables scripts produce `.Rout` and `.pyout` files that track whether the product is up-to-date. This means that if those figures/tables are deleted, the Makefile will consider the product up to date so long as the `.Rout` or `.pyout` file still exists and is up-to-date. These `. out` files will be cleared by `make clean_products`.

#### Clean up
Data output by the pipeline can be quickly and easily deleted by running `make clean_all`, which will delete everything output by the pipeline, except for the archives. In order to delete production run outputs run `make clean_all PROD=1`. Add `PROD=1` to any of the following commands to clean in production mode.

**Beware** that these commands will **delete everything** within the output folders. Do not use these commands if you store data unrelated to the pipeline within these folders, as this could lead to undesired data loss.

- In order to clean up files in the preprocessing folder run `make clean_preprocessed_data`.
- In order to clean up files in the stan fits and stan posterior folder run `make clean_stan_fits` and `make clean_stan_posteriors` respectively.
- In order to clean up all files in the scenarios folder, run `make clean_scenarios`.
- In order to clean up all files related to figures and tables, run `make clean_products`.

## Run the pipeline script-by-script (without using Makefile)

At a high level, the process consists of the following:

1. Run `etl/inherent_preprocess_dataset.py` (prepare raw data from INHERNT study)
2. Run `etl/rvtn_preprocess_dataset.py` (prepare raw data from RVTN study)
3. Run `analysis/scripts/01_prepare_data_for_stan.R` (prepare & combine raw data from both studies for input to Stan)
4. Run `analysis/scripts/02_fit_viral_kinetics_by_symp_cat.R` (Stan model fitting)
5. Run `analysis/scripts/03_extract_parameters_from_stan_fit.R` (parameter extraction from fitted Stan model)
6. Run `analysis/scripts/04_subsample_extracted_parameters.R` (takes a sample of parameters)
7. Run `analysis/scripts/05_run_ind_scenarios.R` (simulates the intervention)
8. Run `analysis/scripts/06_create_overall_category.R` (samples output to represent all symptom categories)
9. Run each script in `analysis/scripts/products` (generates figures and tables)

Individual Python scripts should be run by entering `poetry run python file/path/script_name.py` into the command line. In order to run in production mode, run `poetry run python file/path/script_name.py --output-mode production`.

Individual R scripts should be run by entering `Rscript file/path/script_name.R` into the command line. To run in production mode, run `Rscript file/path/script_name.R production`.

>[!IMPORTANT]
>For users interested in using the publicly available files to replicate portions of the analysis, the script-by-script approach is most straightforward. Specifically, users can run the intervention model (model of transmission potential reduced by isolation and post-isolation precautions) and sample from the output using the following commands:
>```
>Rscript analysis/scripts/05_create_overall_category.R production
>Rscript analysis/scripts/06_create_overall_category.R production
>```
>See below for a description of what these scripts do. Note, however, that the output from these scripts is already included in the publicly available files; indeed, these scripts will recreate and overwrite the downloaded files `intervention_params.csv`, `compare_guidance_df.csv`, and `overall_compare_guidance_df.csv`. Users interested solely in replicating figures and tables do not need to run this part of the pipeline.
>The earlier parts of the pipeline require non-public files and therefore cannot be directly replicated. However, we have created an artificial simulated dataset that can be used to demonstrate fitting the Stan model and the process of extracting parameters from the fitted model. (The code to create this dataset is included in `analysis/tests/test-stan_model.R` and the dataset is generated as a side effect of running `make isolationtest`.) A copy of the artificial simulated dataset is located in the `development` subdirectory of `covid-isolation-outputs` and therefore can be accessed when running in development mode. Specifically, users can run the Stan model and extract the parameters with the following commands:
>```
>Rscript analysis/scripts/02_fit_viral_kinetics_by_symp_cat.R
>Rscript analysis/scripts/03_extract_parameters_from_stan_fit.R
>Rscript analysis/scripts/04_subsample_extracted_parameters.R
>```

### Fitting the viral load and symptom improvement model

1. Run `analysis/scripts/01_prepare_data_for_stan.R`. This takes the data from the ETL phase and preprocesses it specifically for the Stan model.
     - This script uses `isolation/R/stan_preprocess.R` to process the joined INHERENT + RVTN data for Stan.
     - This use the function defined in `isolation/R/calculate_symp_cat_weights.R`, which calculates the proportion of study participants in each symptom category. These proportions serve as weights for later use in the analysis.

2. Run `analysis/scripts/02_fit_viral_kinetics_by_symp_cat.R`. This runs the viral kinetics model written in Stan for every symptom categories defined previously in the ETL phase.
    - This will run the _Stan model_ `isolation/inst/stan/triangle_vl_symp_hazard.stan`. It uses the function defined in `isolation/R/generate_stan_input_list.R` to turn the dataframe into a Stan ready data list. Using Stan's MCMC algorithm, the viral load model is fit and fitted parameters — both pooled and individual-level — are saved. Note that as part of the Stan model run, not only are parameters fit but simulated symptom improvement times are generated (based on the fitted parameters). This is useful because it means that all quantities necessary for assessing infectiousness and modeling isolation guidance are generated at once in the model.
    - This will use the function defined in `isolation/R/vl_fits_to_pdf.R` to export PDFs of the individual-level VL fits for visual inspection.
    - This script outputs the data as used in the model for later analysis and the VL fit samples.

3. Run `analysis/scripts/03_extract_parameters_from_stan_fit.R`.
    - This extracts the posterior samples — both for global parameters and individual-specific parameters — from the Stan fit for each of the four symptom categories using the function defined in `isolation/R/extract_posteriors_to_df.R`. This yields a large output file, where the number of rows of the global individual-level parameter samples is the number of study participants times the number of independent MCMC chains times the number of post-burn-in iterations in each chain. (Recall that the hierarchical model fits a set of parameters for the viral load "triangle" for each study participant, so each row is a particular draw of the individual-level parameters for a given individual.) The global parameter file contains as many rows as there are chains times post-burn-in MCMC samples.

4. Run `analysis/scripts/04_subsample_extracted_parameters.R`.
    - This script takes the output from `analysis/scripts/03_extract_parameters_from_stan_fit.R` and randomly draws a defined number of individual-level posterior parameters sets (i.e., rows from the output dataset) for each symptom category. At this point in the pipeline, we have exactly modeled the individuals present in the real data. Our parameters are specific to the particular people modeled, so for example each parameter sample is tied to a real person in the data. We seek to break this dependency, so we take a random sample of individual-level posterior parameter sets. Individual-level parameter sets are associated both with a person and a chain number (i.e., when in the MCMC chain that sample was drawn). By randomly sampling all parameter sets, agnostic to these assignments, we break the dependency of the parameter sets to real, specific individuals in the observed data. After this script, we only have parameter sets that can be thought of as representing simulated individuals who, collectively, are statistically representative of the study participants.

### Intervention simulation model: Isolation and post-isolation scenarios

1. Run `analysis/scripts/05_run_ind_scenarios.R`. This is the key script that takes as input samples from the posterior distribution (as selected by the previous script) and processes them to simulate interventions (i.e., the previous and updated guidance) and calculate for each simulated person with COVID-19 the proportion of transmission potential averted by the updated versus previous guidance. This script is entirely deterministic. (The algorithm coded in `isolation/R/simulate_isolation_current.R` that determines when a person isolates and takes post-isolation precautions under the previous guidance depends on probabilities of testing positive on an antigen test, but it returns expected values for the probability that a person will be in isolation or post-isolation precautions at a given time, so no Monte Carlo simulation is involved.)
    - This script uses parameters as defined in `analysis/parameters.yml` to create a data frame with one row of intervention scenario parameter levels for use in each analysis (the main analysis plus the sensitivity analyses). This data frame is saved as output.
    - This script loops through each simulated person with COVID-19 (i.e., each set of parameters subsampled from the hierarchical model posterior output) and uses the parameters to create vectors for the viral load, antigen positivity, and culture positivity over time for that person. These trajectories are re-used across analyses for each individual.
    - For each individual, the script also loops through all the analyses (main and sensitivity). This is where the interventions are simulated and the previous and updated guidance are compared. For each analysis, the output values for each simulated individual include the proportion of transmission potential averted under the previous and updated guidance and expected time spent in isolation and in post-isolation precautions under the previous and updated guidance, as well as the differences between the guidance (updated minus previous) for each of these outputs. This output is saved as `compare_guidance_df.csv`; this is the primary data set that is summarized, tabulated, and plotted to present the results of the analysis. The number of rows in this data frame equals the number simulated persons subsampled from each symptom category times the number of symptom categories times the number of analyses.
2. Run `analysis/scripts/06_create_overall_category.R`. This script takes as input the `compare_guidance_df.csv` and randomly samples from this data frame in proportion to the symptom category weights. The purpose is to create a synthetic "overall" category that averages across symptom categories in proportion to the prevalence of each symptom category among the study participants.

### Generate products (figures and tables)
Sixteen scripts that create the figures and tables associated with this analysis are located in `analysis/scripts/products`. Each one would need to be run individually to replicate all products.

>[!IMPORTANT]
>External users interested in replicating figures and tables using the publicly released files should consult the "Figures and tables" list in [`output_docs.md`](output_docs.md) to determine whether a given figure or table can be generated from the publicly available files and, if so, which script to run. Follow the instructions above for running individual scripts. Remember to run scripts in production mode.
