##=========================================================#
## Authors: CMEI Team
## Isolation Guidance Model
##=========================================================#
include .env
include analysis/tests/test_parameters

PRODS_CONFIG := products_config.yml
DEPENDENCIES_MAKE := deps/Dependencies.mk
MAKEFILE_INPUT := deps/Makefile.in

all: etl stan scenarios products

python := poetry run python

$(shell $(python) deps/build_deps.py -y $(PRODS_CONFIG) -o $(DEPENDENCIES_MAKE) -m $(MAKEFILE_INPUT))

.PHONY: clean_products clean_preprocessed_data clean_stan_fits clean_stan_posteriors clean_scenarios clean_all

include $(DEPENDENCIES_MAKE)
##===========================================================#
## Specify output mode as make DEBUG=N from commandline
## PROD = 1 -> production
## PROD != 1 -> development
PROD = 0
OUT_MODE = development

ifeq ($(PROD), 1)
	OUT_MODE := production
else
	OUT_MODE := development
endif

FLAGS := --output-mode $(OUT_MODE)
RFLAGS := $(OUT_MODE)


OUTPUT_DIR := $(MAIN_OUTPUT_DIR)/$(OUT_MODE)
DATA_OUT_DIR := $(OUTPUT_DIR)/$(PREPROCESS_DIR)
ISOLATION_DIR:= isolation
SCRIPTS_DIR := analysis/scripts
FIGS_TABLES_DIR := analysis/products
##=========================================================#
## ETL ---------------------
##=========================================================#
ETL_DIR = etl
ETL_UTILS := $(wildcard $(ETL_DIR)/utils/*.py)


## output files
ETL_FILES := $(DATA_OUT_DIR)/$(RVTN_FILE)\
	$(DATA_OUT_DIR)/$(RVTN_PATIENT_FILE)\
	$(DATA_OUT_DIR)/$(INHERENT_FILE)


etl: $(DATA_OUT_DIR)/$(INHERENT_FILE) $(DATA_OUT_DIR)/$(RVTN_FILE)

$(DATA_OUT_DIR)/$(INHERENT_FILE): $(ETL_DIR)/inherent_preprocess_dataset.py $(DATA_INPUT_DIR)/$(INHERENT_DEPS) $(ETL_UTILS) #nolint: line_length_linter.
	@echo "Running Inherent Preprocess Dataset..." $(<)
	($(python) $(<) $(FLAGS))

$(DATA_OUT_DIR)/$(RVTN_FILE): $(ETL_DIR)/rvtn_preprocess_dataset.py $(addprefix $(DATA_INPUT_DIR)/, $(RVTN_DEPS)) $(ETL_UTILS) #nolint: line_length_linter.
	@echo "Running RVTN Preprocess Dataset..." $(<)
	($(python) $(<) $(FLAGS))

##===========================================================##
## Isolation package---------------
##===========================================================##
## Dependencies
RLIB_FILES := $(wildcard $(ISOLATION_DIR)/R/*.R) \
	$(wildcard $(ISOLATION_DIR)/inst/stan/*.stan)

## Target
isolationtest: analysis/tests/test-stan_model.R
	@echo $(RLIB_FILES)
	(Rscript $(<) --samples=$(samples) --pairs=$(pairs) --start-day=$(start_day))
##===========================================================##
## STAN---------------
##===========================================================##
## Dependencies
STAN_DEPS := analysis/parameters.yml $(PRODS_CONFIG) $(RLIB_FILES)
STAN_PROCESSED_DIR := $(OUTPUT_DIR)/$(PREPROCESS_DIR)

STAN_INPUT_DATA := $(STAN_PROCESSED_DIR)/$(STAN_PROCESSED_FILE)
STAN_FITS := $(OUTPUT_DIR)/$(STAN_FIT_DIR)/$(STAN_FIT_FILE)4.rds
STAN_POSTERIORS := $(OUTPUT_DIR)/$(STAN_POST_DIR)/$(STAN_POST_FILE)
STAN_SAMPLES := $(OUTPUT_DIR)/$(STAN_POST_DIR)/$(STAN_SAMPLE_FILE)

stan: $(STAN_INPUT_DATA) $(STAN_FITS) $(STAN_POSTERIORS) $(STAN_SAMPLES)

$(STAN_INPUT_DATA): $(SCRIPTS_DIR)/01_prepare_data_for_stan.R $(STAN_DEPS) $(ETL_FILES) #nolint: line_length_linter.
	@echo "Preparing data for stan..." $(<)
	(Rscript $(<) $(RFLAGS))

$(STAN_FITS): $(SCRIPTS_DIR)/02_fit_viral_kinetics_by_symp_cat.R $(STAN_DEPS) $(STAN_INPUT_DATA) #nolint: line_length_linter.
	@echo "Fitting stan model..." $(<)
	(Rscript $(<) $(RFLAGS))

$(STAN_POSTERIORS): $(SCRIPTS_DIR)/03_extract_parameters_from_stan_fit.R $(STAN_DEPS) $(STAN_FITS) #nolint: line_length_linter.
	@echo "Extracting stan params..." $(<)
	(Rscript $(<) $(RFLAGS))

$(STAN_SAMPLES): $(SCRIPTS_DIR)/04_subsample_extracted_parameters.R $(STAN_DEPS) $(STAN_POSTERIORS) #nolint: line_length_linter.
	@echo "Subsampling extracted parameters..." $(<)
	(Rscript $(<) $(RFLAGS))

##===========================================================##
## Scenarios---------------
##===========================================================##
## Dependencies
SCENARIO_DEPS := analysis/parameters.yml $(PRODS_CONFIG) $(RLIB_FILES)\
	$(STAN_POSTERIORS) $(STAN_SAMPLES)

SCENARIOS_OUT_DIR := $(OUTPUT_DIR)/$(SCENARIOS_DIR)
SCENARIOS_OUT_FILE := $(SCENARIOS_OUT_DIR)/$(GUIDANCE_FILE)
OVERALL_SCENARIOS_OUT_FILE := $(SCENARIOS_OUT_DIR)/$(OVERALL_GUIDANCE_FILE)

scenarios: $(SCENARIOS_OUT_FILE) $(OVERALL_SCENARIOS_OUT_FILE)

$(SCENARIOS_OUT_FILE): $(SCRIPTS_DIR)/05_run_ind_scenarios.R $(SCENARIO_DEPS)
	@echo "Running individual scenarios: " $(<)
	(Rscript $(<) $(RFLAGS))

$(OVERALL_SCENARIOS_OUT_FILE): $(SCRIPTS_DIR)/06_create_overall_category.R $(SCENARIO_DEPS) #nolint: line_length_linter.
	@echo "Aggregating results across all symptom categories: " $(<)
	(Rscript $(<) $(RFLAGS))

##===========================================================##
## Figures and tables---------------
##===========================================================##
#include deps/makefigs_regex.mk
include deps/dep_products.mk


##===========================================================##
## Clean up---------------
##===========================================================##
clean_all: clean_preprocessed_data clean_stan_fits clean_stan_posteriors clean_scenarios clean_products

clean_preprocessed_data:
	@echo "Clearing files in:" $(OUTPUT_DIR)/$(PREPROCESS_DIR)
	(rm -f --verbose $(OUTPUT_DIR)/$(PREPROCESS_DIR)/*.*)

clean_stan_fits:
	@echo "Clearing files in:" $(OUTPUT_DIR)/$(STAN_FIT_DIR)
	(rm -f --verbose $(OUTPUT_DIR)/$(STAN_FIT_DIR)/*.*)

clean_stan_posteriors:
	@echo "Clearing files in:" $(OUTPUT_DIR)/$(STAN_POST_DIR)
	(rm -f --verbose $(OUTPUT_DIR)/$(STAN_POST_DIR)/*.*)

clean_scenarios:
	@echo "Clearing files in:" $(SCENARIOS_OUT_DIR)
	(rm -f --verbose $(SCENARIOS_OUT_DIR)/*.*)
##===========================================================##
