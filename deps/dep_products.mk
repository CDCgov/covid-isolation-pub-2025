PRODUCT_DEPS := analysis/parameters.yml $(prods_config)\
	$(STAN_POSTERIORS) $(STAN_SAMPLES)\
	$(STAN_INPUT_DATA) $(SCENARIOS_OUT_FILE) $(OVERALL_SCENARIOS_OUT_FILE)

PRODUCTS_PATH := $(OUTPUT_DIR)/products
SCRIPT_PATH := analysis/scripts/products

deps_list := deps/products_list

$(shell find $(PRODUCTS_PATH) -maxdepth 1 -not -type d >> $(deps_list))
$(shell	sort $(deps_list) | uniq > $(deps_list)2)
$(shell mv $(deps_list)2 $(deps_list))

CURR_PRODS := $(basename $(notdir $(shell find $(PRODUCTS_PATH) -maxdepth 1 -not -type d)))
ALL_PRODS := $(notdir $(shell cat $(deps_list)))

## FIND DEPENDENCIES =========================#
PROD_SRC := $(notdir $(wildcard $(SCRIPT_PATH)/*.R) $(wildcard $(SCRIPT_PATH)/*.py))
NOTPRODS := $(filter-out $(basename $(ALL_PRODS)), $(basename $(PROD_SRC)))
INPRODS := $(filter $(basename $(ALL_PRODS)), $(basename $(PROD_SRC)))
DIFFPRODS := $(filter-out $(CURR_PRODS), $(filter-out $(INPRODS), $(basename $(ALL_PRODS))))

TARGET_FIGS := $(filter $(addsuffix .png, $(INPRODS)) $(addsuffix .pdf, $(INPRODS)),  $(ALL_PRODS))
TARGET_TABLES := $(filter $(addsuffix .csv, $(INPRODS)) $(addsuffix .xlsx, $(INPRODS)), $(ALL_PRODS))

TARGET_R :=  $(patsubst %.R, %.Rout, $(filter $(addsuffix .R, $(NOTPRODS)), $(PROD_SRC)))
TARGET_PYTHON := $(patsubst %.py, %.pyout, $(filter $(addsuffix .py, $(NOTPRODS)), $(PROD_SRC)))

ifneq ($(strip $(DIFFPRODS)), )
$(shell rm $(addprefix $(SCRIPT_PATH)/, $(TARGET_R) $(TARGET_PYTHON)))
endif

## BUILD FIGURES AND TABLES ==================#
products: $(addprefix $(PRODUCTS_PATH)/, $(TARGET_FIGS)) \
	$(addprefix $(PRODUCTS_PATH)/, $(TARGET_TABLES)) \
	$(addprefix $(SCRIPT_PATH)/, $(TARGET_R)) \
	$(addprefix $(SCRIPT_PATH)/, $(TARGET_PYTHON))

$(addprefix $(PRODUCTS_PATH)/, %.pdf %.png %.csv %.xlsx): $(SCRIPT_PATH)/%.R $(PRODUCT_DEPS)
	(Rscript $< $(RFLAGS))

$(addprefix $(PRODUCTS_PATH)/, %.pdf %.png %.csv %.xlsx): $(SCRIPT_PATH)/%.py $(PRODUCT_DEPS)
	($(python) $< $(FLAGS))

$(addprefix $(SCRIPT_PATH)/, %.Rout): $(SCRIPT_PATH)/%.R
	(Rscript $< $(RFLAGS) > $@)

$(addprefix $(SCRIPT_PATH)/, %.pyout): $(SCRIPT_PATH)/%.py
	($(python) $< $(FLAGS) > $@)

productecho:
	@echo "Products not generated yet: " $(DIFFPRODS)


##===========================================================##
## Clean---------------
##===========================================================##
clean_products:
	@echo "Clearing products in:" $(PRODUCTS_PATH)
	(rm -f --verbose $(PRODUCTS_PATH)/*.pdf)
	(rm -f --verbose $(PRODUCTS_PATH)/*.png)
	(rm -f --verbose $(PRODUCTS_PATH)/*.csv)
	(rm -f --verbose $(PRODUCTS_PATH)/*.xlsx)
	@echo "Removing .Routs and .pyouts in:" $(SCRIPT_PATH)
	(rm -f --verbose $(SCRIPT_PATH)/*.Rout)
	(rm -f --verbose $(SCRIPT_PATH)/*.pyout)
	@echo "Removing products_list"
	(rm -f --verbose deps/products_list)
