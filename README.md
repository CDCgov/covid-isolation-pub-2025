# covid-isolation-pub-2025

## Overview

This repository contains code necessary to reproduce a modeling analysis on the duration of isolation and post-isolation precautions for COVID-19. This analysis compares CDC's current guidance for [Preventing Spread of Respiratory Viruses When You're Sick](https://www.cdc.gov/respiratory-viruses/prevention/precautions-when-sick.html) ([updated](https://www.cdc.gov/ncird/whats-new/updated-respiratory-virus-guidance.html) on March 1, 2024) to the guidance for COVID-19 that existed previously.

An earlier version of this analysis is summarized at [Behind the Model: Estimating the Impact of Updated Isolation Guidance on COVID-19 Transmission](https://www.cdc.gov/cfa-behind-the-model/php/data-research/covid-isolation/index.html).

This analysis is being submitted for publication. It is subject to revision as part of the peer-review process.

## Getting Started

>[!IMPORTANT]
>This repository includes the code for the entire analytic pipeline, including the extract, transform, and load (ETL) process to prepare the data for analysis. Because the underlying data are not publicly available, external users will not be able to replicate the full pipeline. But as part of this release, we sharing a zip file that includes model output that can be used to replicate part of the pipeline and most figures and tables. We use notes like this one to alert external users to what portions of the analysis they can replicate using the model output files.

### Set up computing environment

This analysis was performed on Unbuntu 22.04.3 LTS running under Windows Subsystem for Linux version 2 (WSL2). Users may wish to consider using MacOS, Linux, or Windows Subsystem for Linux (WSL). Portions of the code presume that the analysis is being performed on a computer with at least 4 cores for parallel processing.

1. Create a `.env` file in the repository to specify the user's preferred input and output directories. The input directory is for the "raw" data files that are fed into the ETL part of the pipeline. The output directory is for all other files, including intermediary files that are output from a step in the pipeline and then input into one or more downstream steps. To create this file, create a copy of the [`.env.example`](.env.example) file that is named simply `.env` and then edit it as necessary. Use only directory names; do not leave comments on the same line as the directory path.

>[!IMPORTANT]
>Users interested in replicating portions of the analysis that do not require non-public participant-level data should download from the [release](https://github.com/cdcgov/covid-isolation-pub-2025/releases) the file `covid-isolation-outputs.zip` and unzip it to their preferred location and set the value for `MAIN_OUTPUT_DIR` in their `.env` file accordingly. It is important to maintain the file structure from the zip file. For example, if a user places the directory `covid-isolation-outputs` in their `D` drive, then the path of the file `VL_symp_cat_sampled.csv` would be `~/D/covid-isolation-outputs/production/viral_kinetics_posteriors/VL_symp_cat_sampled.csv` and the `.env` file would need to include the line `export MAIN_OUTPUT_DIR=~/D/covid-isolation-outputs`. The path specified for `DATA_INPUT_DIR` is only relevant when running the ETL portion of the pipeline and is therefore not applicable for most users. Documentation on the contents of `covid-isolation-outputs.zip` is available at [`analysis/output_docs.md`](analysis/output_docs.md).

2. Install specified versions of [Python](https://www.python.org/downloads/) and [R](https://cran.r-project.org/). This analysis was performed in Python 3.10 and R 4.4.1. These versions are specified in [`pyproject.toml`](pyproject.toml) and in [`renv.lock`](renv.lock).
3. Install and run [poetry](https://python-poetry.org/).
    - See poetry documentation for [installation instructions](https://python-poetry.org/docs/#installation).
    - Run `poetry install` with the working directory set to the directory for this project (e.g., `covid-isolation-pub-2025`). This will install the necessary dependencies and create a virtual environment. (This analysis uses the following Python libraries: `numpy`, `pandas`, `python-dotenv`, `matplotlib`, `pyyaml`, `openpyxl`, `scipy`, and dependencies thereof. Users should not need to perform manual installations as long as they use poetry.)
    - Sometimes users will need to activate the poetry environment. If modules are still not found after poetry install, running `$(poetry env activate)` or `source $(poetry env info --path)/bin/activate` will activate the environment.
4. Install and run [renv](https://rstudio.github.io/renv/).
    - Start an interactive `R` session and use the command `install.packages("renv")`.
    - In an interactive `R` session with the working directory set to the directory for this project (e.g., `covid-isolation-pub-2025`), run `renv::restore()`. This will install the necessary dependencies. This may take some time (e.g., 20 or more minutes). (This analysis uses the following R packages: `colorspace`, `devtools`, `DiagrammeR`, `DiscreteWeibull`, `dotenv`, `dplyr`, `ggplot2`, `htmltools`, `openxlsx`, `purrr`, `readr`, `renv`, `Rstan`, `testthat`, `tidyr`, `tidyselect`, `tibble`, `webshot`, `yaml`, and dependencies thereof. Users should not need to perform manual installations as long as they use renv.)
    - To verify successful installation of dependencies in accordance with the lockfile, the command `renv::status` in an interactive `R` session should return the statement `No issues found -- the project is in a consistent state.`.

### Specify parameters and intermediary output file names

Model parameters (including priors and parameter values for sensitivity analyses) are specified in the file [`analysis/parameters.yml`](analysis/parameters.yml). This file is used throughout the analytic pipeline. The names for each intermediary output file are stored in [`products_config.yml`](products_config.yml).

### Production vs Development

This repository can be run in two different modes: production mode and development mode. Development mode will run files on a subset of the full data used in production mode in order to allow the user to run analyses faster and test changes. The proportion of the full production data used during a development run can be set by changing the `development_dataset_subset` parameter in the [`analysis/parameters.yml`](analysis/pipeline.md) file.

Details on how to change between production and development are included in [`analysis/pipeline.md`](analysis/pipeline.md).

### Next Steps

See [`analysis/pipeline.md`](analysis/pipeline.md) for further instructions on running the analysis.

## Project Admin

Eric Q. Mooring (pgv5@cdc.gov) and Guido España (ukd0@cdc.gov)

## General Disclaimer

This repository was created for use by CDC programs to collaborate on public health related projects in support of the [CDC mission](https://www.cdc.gov/about/organization/mission.htm).  GitHub is not hosted by the CDC, but is a third party website used by CDC and its partners to share information and collaborate on software. CDC use of GitHub does not imply an endorsement of any one particular service, product, or enterprise.

## Public Domain Standard Notice

This repository constitutes a work of the United States Government and is not
subject to domestic copyright protection under 17 USC § 105. This repository is in
the public domain within the United States, and copyright and related rights in
the work worldwide are waived through the [CC0 1.0 Universal public domain dedication](https://creativecommons.org/publicdomain/zero/1.0/).
All contributions to this repository will be released under the CC0 dedication. By
submitting a pull request you are agreeing to comply with this waiver of
copyright interest.

## License Standard Notice

This repository is licensed under ASL v2 or later.

This source code in this repository is free: you can redistribute it and/or modify it under
the terms of the Apache Software License version 2, or (at your option) any
later version.

This source code in this repository is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE. See the Apache Software License for more details.

You should have received a copy of the Apache Software License along with this
program. If not, see <http://www.apache.org/licenses/LICENSE-2.0.html>

The source code forked from other open source projects will inherit its license.

## Privacy Standard Notice

This repository contains only non-sensitive, publicly available data and
information. All material and community participation is covered by the
[Disclaimer](https://github.com/CDCgov/template/blob/master/DISCLAIMER.md)
and [Code of Conduct](https://github.com/CDCgov/template/blob/master/code-of-conduct.md).
For more information about CDC's privacy policy, please visit [http://www.cdc.gov/other/privacy.html](https://www.cdc.gov/other/privacy.html).

## Contributing Standard Notice

Anyone is encouraged to contribute to the repository by [forking](https://help.github.com/articles/fork-a-repo)
and submitting a pull request. (If you are new to GitHub, you might start with a
[basic tutorial](https://help.github.com/articles/set-up-git).) By contributing
to this project, you grant a world-wide, royalty-free, perpetual, irrevocable,
non-exclusive, transferable license to all users under the terms of the
[Apache Software License v2](http://www.apache.org/licenses/LICENSE-2.0.html) or
later.

All comments, messages, pull requests, and other submissions received through
CDC including this GitHub page may be subject to applicable federal law, including but not limited to the Federal Records Act, and may be archived. Learn more at [http://www.cdc.gov/other/privacy.html](http://www.cdc.gov/other/privacy.html).

## Records Management Standard Notice

This repository is not a source of government records but is a copy to increase
collaboration and collaborative potential. All government records will be
published through the [CDC web site](http://www.cdc.gov).

## Additional Standard Notices

Please refer to [CDC's Template Repository](https://github.com/CDCgov/template)
for more information about [contributing to this repository](https://github.com/CDCgov/template/blob/master/CONTRIBUTING.md),
[public domain notices and disclaimers](https://github.com/CDCgov/template/blob/master/DISCLAIMER.md),
and [code of conduct](https://github.com/CDCgov/template/blob/master/code-of-conduct.md).
