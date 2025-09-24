#' Preprocess data to input into Stan model
#'
#' Data processing and associated checks to prepare it for the Stan model
#'
#' @param data_to_prep INHERENT/RVTN data to prep for stan
#' @param parameters_config a config specifying special values of this dataset
#' @param dataset_name data set being analyzed here's name to match config
#' @param at_least_one_test only include days for which there is >= 1 test?
#'
#' @return a list that contains 1) processed data and 2) a mermaid diagram
#' @export
stan_preprocess <- function(data_to_prep,
                            parameters_config,
                            dataset_name,
                            at_least_one_test = TRUE) {
  # the parameters_config must have a section that describes
  # config parameters specific to this data set
  stopifnot(dataset_name %in% names(parameters_config))
  dataset_config <- parameters_config[[dataset_name]]
  global_data_markers <- parameters_config$global_data_markers

  mermaid_diagram <- "graph TB\n"

  # filter ETL data for Stan
  mermaid_diagram <- paste0(
    mermaid_diagram,
    "A[",
    length(unique(data_to_prep$id)),
    " symptomatic individuals in ",
    dataset_name,
    "]-->"
  )
  # filter to make sure everyone who is symptomatic also has symp counts
  max_symp_count_expd_zero <- data_to_prep |>
    dplyr::group_by(id) |>
    dplyr::summarise(max_symp_count = max(symp_count)) |>
    dplyr::filter(max_symp_count == 0)

  data_to_prep <- data_to_prep |>
    dplyr::filter(!(id %in% max_symp_count_expd_zero$id))
  mermaid_diagram <- paste0(
    mermaid_diagram,
    "B[",
    length(unique(data_to_prep$id)),
    " symptomatic individuals with recorded symptoms in ",
    dataset_name,
    "]\n"
  )

  data_to_prep <- data_to_prep |>
    # filter to only include relevant participant type
    dplyr::filter(participant_type %in% dataset_config$participant_type_keep)
  if (dataset_name == "INHERENT") {
    participant_type_str <- " -- nursing home staff -- "
  } else if (dataset_name == "RVTN") {
    participant_type_str <- " -- household contacts -- "
  }
  mermaid_diagram <- paste0(
    mermaid_diagram,
    "B-->C[",
    length(unique(data_to_prep$id)),
    " relevant individuals",
    participant_type_str,
    "in ",
    dataset_name,
    "]\n"
  )
  data_to_prep <- data_to_prep |>
    # filter to only include people whose logVL is quantified
    # at least once during infection
    dplyr::filter(logVL_quantifiable_ever == 1)
  mermaid_diagram <- paste0(
    mermaid_diagram,
    "C-->D[",
    length(unique(data_to_prep$id)),
    " individuals with >= 1 quantifiable logVL in ",
    dataset_name,
    "]\n"
  )

  # now add how the people shake out by category
  symp_cat_table <- as.data.frame(table(factor(
    unique(data_to_prep[
      ,
      c(
        "id",
        "symp_type_cat"
      )
    ])$symp_type_cat,
    levels = parameters_config$global_data_markers$symptom_categories
  )))
  for (symp_type_cat in symp_cat_table$Var1) {
    mermaid_diagram <- paste0(
      mermaid_diagram,
      "D-->",
      LETTERS[4 + as.integer(symp_type_cat)], # 4 bc we've used A, B, C, D
      "[",
      symp_cat_table[symp_cat_table$Var1 == symp_type_cat, "Freq"],
      " individuals in Symptom Category ",
      symp_type_cat,
      " in ", dataset_name,
      "]\n"
    )
  }

  # there should be no individuals who are not symptomatic
  # meaning that everybody has a symptom duration
  stopifnot(data_to_prep$symp_duration > 0)
  stopifnot(data_to_prep$logVL_quantifiable_ever != 0)
  stopifnot(
    data_to_prep$participant_type %in% dataset_config$participant_type_keep
  )

  # truncate si times greater than max_si
  # consists of two parts -- making sure no si times are greater than max
  # but also that no censored times are above that threshold because then
  # how could they be observed? recall no density above max_si
  # first print out so the user can see when/how often this happens
  cases_where_symp_duration_greater_than_max_si <- data_to_prep |> # nolint: object_length_linter, line_length_linter.
    dplyr::filter(symp_duration >= parameters_config$scenario_parameters$max_si)
  print(paste(
    "There are",
    length(unique(cases_where_symp_duration_greater_than_max_si$id)),
    "individuals who have symptom improvement times greater
  than the max time set in config"
  ))
  if (length(unique(cases_where_symp_duration_greater_than_max_si$id)) > 0) {
    print("The IDs and symptom improvement times are:")
    cases_where_symp_duration_greater_than_max_si <- cases_where_symp_duration_greater_than_max_si |> # nolint: object_length_linter, line_length_linter.
      dplyr::select(id, symp_duration, symp_duration_censored) |>
      dplyr::group_by(id) |>
      dplyr::summarise(
        symp_duration = max(symp_duration),
        symp_duration_censored = max(symp_duration_censored)
      )
    print(cases_where_symp_duration_greater_than_max_si)
    print(paste(
      "Recall the max symptom improvement time set in config is",
      parameters_config$scenario_parameters$max_si
    ))
    print("For people with a symptom improvement time of greater than the max,
  their symptom improvement times will be truncated
  to the max for modeling purposes")
  }
  data_to_prep <- data_to_prep |>
    dplyr::mutate(
      symp_duration = dplyr::case_when(
        symp_duration > parameters_config$scenario_parameters$max_si ~
          parameters_config$scenario_parameters$max_si,
        TRUE ~ symp_duration
      ),
      symp_duration_censored = dplyr::case_when(
        # need a >= here because symp_duration == max_si also can't be censored
        symp_duration >= parameters_config$scenario_parameters$max_si ~ 0,
        TRUE ~ symp_duration_censored
      )
    )

  # all symp improv times are less than or equal to max_si
  stopifnot(
    data_to_prep$symp_duration <=
      parameters_config$scenario_parameters$max_si
  )

  # all censored symp improv times are strictly less than max_si
  stopifnot(dplyr::filter(
    data_to_prep,
    symp_duration_censored == 1
  )$symp_duration <
    parameters_config$scenario_parameters$max_si)

  # change variable values to be consistent across data sets bc we concatenate
  # the datasets together eventually
  data_to_prep <- data_to_prep |>
    # parameters_config$global_data_markers$skipped_test is the
    # value for a "skipped test" that we adopt (meaning no test that day)
    # we differentiate between skipped test and "inconclusive" or the like
    # skipped/inconclusive tests are set to `antigen_missing_flat`
    dplyr::mutate(antigen = dplyr::case_when(
      antigen == global_data_markers$antigen_missing ~
        global_data_markers$skipped_test,
      TRUE ~ antigen
    )) |>
    dplyr::mutate(culture = dplyr::case_when(
      culture == global_data_markers$culture_missing ~
        global_data_markers$skipped_test,
      TRUE ~ culture
    )) |>
    # keep the LOD for this dataset in a column
    dplyr::mutate(lod = dataset_config$lod) |>
    # keep the LOQ for this dataset in a column
    dplyr::mutate(loq_lower = dataset_config$loq_lower)
  # keep the loq_upper for this dataset in a column
  # note that if there is no loq_upper in that datasets config,
  # R will assign this a null
  # but, Stan doesn't know how to read in a null value
  # so we have to reassign the null as something -- very large
  loq_upper <- dataset_config$loq_upper
  if (is.null(loq_upper)) {
    loq_upper <- global_data_markers$skipped_test
  }
  data_to_prep <- data_to_prep |>
    dplyr::mutate(loq_upper = loq_upper) |>
    # change any of the values of logVL that are the LOQ
    # to a global LOQ marker
    # this is just changing the marker to be consistent across all datasets
    # so it is easier to check for one marker in Stan
    # the actual LOQ is kept in the line above, and that is used in the interval
    # censoring function in stan
    dplyr::mutate(logVL = dplyr::case_when(
      logVL == dataset_config$loq_lower_marker ~ global_data_markers$loq_lower,
      TRUE ~ logVL
    ))

  # when antigen_missing != skipped_test (which it should never),
  # we never expect there to be any cases of `antigen_missing` in `antigen`
  if (global_data_markers$skipped_test != global_data_markers$antigen_missing) {
    stopifnot(data_to_prep$antigen != global_data_markers$antigen_missing)
  }
  # when culture_missing != skipped_test (which it should never),
  # we never expect there to be any cases of `culture_missing` in `culture`
  if (global_data_markers$skipped_test != global_data_markers$culture_missing) {
    stopifnot(data_to_prep$culture != global_data_markers$culture_missing)
  }
  # loq column should be all the loq value we set
  stopifnot(data_to_prep$loq_lower == dataset_config$loq_lower)
  # same idea as above with the antigen test
  if (dataset_config$loq_lower_marker != global_data_markers$loq_lower) {
    stopifnot(data_to_prep$logVL != dataset_config$loq_lower_marker)
  }

  # grab those days for which there is at least one test available
  # this basically reduces the amount of unnecessary data given to stan model
  # to improve its iteration speed/making sure it only gets the info it needs
  if (at_least_one_test) {
    data_to_prep <- data_to_prep |>
      dplyr::mutate(pcr_test = logVL == global_data_markers$skipped_test) |>
      dplyr::mutate(
        antigen_test = antigen ==
          global_data_markers$skipped_test
      ) |>
      dplyr::mutate(
        culture_assay = culture ==
          global_data_markers$skipped_test
      ) |>
      dplyr::mutate(no_test = pcr_test & antigen_test & culture_assay) |>
      dplyr::mutate(at_least_one_test = !no_test) |>
      dplyr::filter(at_least_one_test)
  }
  # there should be no entries that have test results as entirely missing
  # so, if we groupby id and look at each individual's logVL and antigen,
  # they should not be all equal to `missing_data_flag`
  check_avail_tests <- data_to_prep |>
    dplyr::group_by(id) |>
    dplyr::summarise(
      count_avail = sum(
        logVL != global_data_markers$skipped_test
      ) + sum(
        antigen != global_data_markers$skipped_test
      )
    ) |>
    dplyr::ungroup()
  stopifnot(check_avail_tests$count_avail > 0)

  # select only the variables we need.
  data_to_prep <- data_to_prep |>
    dplyr::select(
      id,
      days_since_symp_onset,
      logVL,
      antigen,
      culture,
      symp_duration,
      symp_duration_censored,
      lod,
      loq_lower,
      loq_upper,
      symp_type_cat,
      age_cat,
      cdc_period,
      months_since_last_covid_dose,
      symp_count
    )

  data_to_prep$contiguous_id <- dplyr::consecutive_id(data_to_prep$id)

  return(list(
    data_to_prep = data_to_prep,
    mermaid_diagram = mermaid_diagram
  ))
}
