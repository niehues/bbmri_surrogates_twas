################################################################################
## Title:         plan.R
## Description:   Definition of analysis workflow using drake
## Author:        Anna Niehues
## Date created:  2020-11-20
## Email:         anna.niehues@radboudumc.nl 
################################################################################
## Notes:
## 
################################################################################

plan <- drake::drake_plan(

  # check if output directory exists
  outdir = check_dir(outdir_path),

  # model parameters (comparisons, models and covariates, cohorts)
  model_parameters = target(
    as.list(models_list[[model_id]]),
    transform = map(
      model_id = !!seq_len(num_models),
      comparison = !!lapply(models_list, function(x) {x$comparison}),
      model = !!lapply(models_list, function(x) {x$model}),
      cohort = !!lapply(models_list, function(x) {x$cohort}),
      .id = c(comparison, model, cohort))),
  
  # groups of models
  model_parameters_by_comparison_by_cohort = target(
    model_parameters,
    transform = combine(model_parameters, 
                        .by = c(comparison, cohort), 
                        .id = c(comparison, cohort))),
  
  # read predicted values - EpiSmokEr
  epismoker_values = target(
    read_predicted_values(epismoker_results_paths[[model_parameters[["cohort"]]]]),
    transform = map(model_parameters, 
                    .id = c(comparison, model, cohort))),

  epismoker_variables = target(
    unique(unlist(lapply(list(epismoker_values), colnames))),
    transform = combine(epismoker_values, .id = FALSE)),
  
  # read predicted values - metabolic scores
  surrogates_values = target(
    read_predicted_values(
      metabolomics_surrogates_paths[[model_parameters[["cohort"]]]]),
    transform = map(model_parameters, 
                    .id = c(comparison, model, cohort))),

  surrogates_variables = target(
    unique(unlist(lapply(list(surrogates_values), colnames))),
    transform = combine(surrogates_values, .id = FALSE)),
  
  # combine predicted values (EpiSmokEr and metabolic scores)
  predicted_values = target(
    combine_metadata(list(epismoker_values, surrogates_values)),
    transform = map(epismoker_values, surrogates_values, 
                    .id = c(comparison, model, cohort))),

  # load and preprocess RNA seq data
  # either from BBMRIomics 
  rnaseq_counts = target(
    load_and_prepare(model_parameters[["cohort"]], predicted_values),
    transform = map(model_parameters, predicted_values, 
                    .id = c(comparison, model, cohort))),

  # filter RNA seq data (per comparison based on metadata availability)
  counts_filtered = target(
    filter_rna(rnaseq_counts, model_parameters[["comparison_variables"]]),
    transform = map(rnaseq_counts, model_parameters,
                    .id = c(comparison, model, cohort))),

  # # define design
  design = target(
    create_design(counts_filtered, model_parameters[["covariates"]]),
    transform = map(counts_filtered, model_parameters,
                    .id = c(comparison, model, cohort))),

  # normalize counts
  counts_normalized = target(
    normalize_rnaseq_data(counts_filtered, design,
                          check_dir(file.path(outdir, model_parameters[["comparison"]])), 
                          paste0(model_parameters[["cohort"]], 
                                         model_parameters[["model"]], 
                                         sep = "_")),
    transform = map(counts_filtered, design, model_parameters,
                    .id = c(comparison, model, cohort))),

  # run TWAS for all alternative models - entire data set
  TWAS = target(
    run_twas(counts_normalized, design),
    transform = map(counts_normalized, design,
                    .id = c(comparison, model, cohort))),
  
  TWAS_file = target(
    save(TWAS,
         file = file.path(outdir, comparison, 
                          paste0("TWAS_", model, "_", cohort, ".RData", sep = ""))),
    transform = map(TWAS, .id = c(comparison, model, cohort))),
  
  # volcano plot
  volcano_plot = target(
    plot_volcano(TWAS,
                 paste0(model_parameters[["cohort"]], ", ", 
                        model_parameters[["model"]]),
                 check_dir(file.path(outdir, model_parameters[["comparison"]]))),
    transform = map(TWAS, 
                    model_parameters, 
                    .id = c(comparison, model, cohort))),

  # TWAS group by cohort
  betas_by_comparison_by_cohort = target(
    lapply(list(TWAS), function(x) {x[["efit"]][["coefficients"]]}),
    transform = combine(TWAS,
                        .by = c(comparison, cohort),
                        .id = c(comparison, cohort))),
  betas_bacon_by_comparison_by_cohort = target(
    lapply(list(TWAS), function(x) {x[["results"]][["eff.size.bacon"]]}),
    transform = combine(TWAS,
                        .by = c(comparison, cohort),
                        .id = c(comparison, cohort))),

  # correlation analysis of betas
  beta_cor_plots = target(
    correlation_analysis(
      betas_by_comparison_by_cohort, 
      model_parameters_by_comparison_by_cohort,
      cohort,
      check_dir(file.path(
        outdir, model_parameters_by_comparison_by_cohort[[1]][["comparison"]]))),
    transform = map(betas_by_comparison_by_cohort, 
                    model_parameters_by_comparison_by_cohort, 
                    .id = c(comparison, cohort))),
  
  # bias and inflation
  bias_inflation = target(
    calc_bias_inflation(TWAS,
                        model_parameters[["comparison"]],
                        model_parameters[["cohort"]],
                        model_parameters[["model"]],
                        alpha),
    transform = map(TWAS, model_parameters,
                    .id = c(comparison, model, cohort))),
 
  # TWAS group by model 
  TWAS_by_comparison_by_model = target(
    setNames(lapply(list(TWAS), function(x) {x[["results"]]}),
             lapply(list(model_parameters), function(x) {x[["cohort"]]})),
    transform = combine(TWAS, model_parameters,
                        .by = c(comparison, model), 
                        .id = c(comparison, model))),
  
  # meta-analysis of all cohorts
  all_meta_analysis = target(
    meta_analysis(
      lapply(
        TWAS_by_comparison_by_model, 
        function(twas) {twas[Reduce(
          intersect, 
          lapply(TWAS_by_comparison_by_model, 
                 function(twas2) {rownames(twas2)})
          ),]})),
    transform = map(TWAS_by_comparison_by_model,
                    .id = c(comparison, model))),
  
  # LOO TWAS meta-analysis -
  loo_meta_analysis = target(
    loo_analysis(TWAS_by_comparison_by_model, comparison, model, alpha,
                 check_dir(file.path(outdir, comparison))),
    transform = map(TWAS_by_comparison_by_model,
                      # rep_cohort = !!models_json[["comparison"]][["cohorts"]],
                      .id = c(comparison, model))),

  # combine all LOO meta-analyses by model
  loo_meta_analysis_all = target(
   do.call(rbind, list(loo_meta_analysis)),
   transform = combine(loo_meta_analysis, .id = FALSE)),
  
  meta_GSEA = target(
    perform_gsea(all_meta_analysis, metaanalysis = TRUE),
    transform = map(all_meta_analysis,
                    .id = c(comparison, model))),
  
  meta_GSEA_by_comparison = target(
      setNames(list(meta_GSEA),
               unique(lapply(list(model_parameters), 
                      function(x) {x[["model"]]}))),
      transform = combine(meta_GSEA, 
                          model_parameters,
                          .by = comparison,
                          .id = comparison)),
  
  GSEA = target(
    perform_gsea(TWAS$results, metaanalysis = FALSE),
    transform = map(TWAS,
                    .id = c(comparison, model, cohort))),

  # combine TWAS by c(comparison, cohort) to compare results
  GSEA_by_comparison_by_cohort = target(
    setNames(list(GSEA),
             lapply(list(model_parameters), function(x) {x[["model"]]})),
    transform = combine(GSEA, model_parameters,
                        .by = c(comparison, cohort),
                        .id = c(comparison, cohort))),

  # Render markdown report
  # TWAS
  report = rmarkdown::render(
    knitr_in("report.Rmd"),
    output_file = file.path(
      outdir, paste0("report.html")),
    quiet = TRUE),
  # # 
  # # GSEA - detailed results
  report_gsea = rmarkdown::render(
    knitr_in("report_gsea.Rmd"),
    output_file = file.path(
      outdir, paste0("report_gsea.html")),
    quiet = TRUE),
  
  report_gsea_supplementary = rmarkdown::render(
    knitr_in("report_gsea_supplementary.Rmd"),
    output_file = file.path(
      outdir, paste0("report_gsea_supplementary.html")),
    quiet = TRUE),
  
  report_model_parameters = rmarkdown::render(
    knitr_in("report_parameters.Rmd"),
    output_file = file.path(
      outdir, paste0("report_parameters.html")),
    quiet = TRUE),
  
  trace = TRUE
)
