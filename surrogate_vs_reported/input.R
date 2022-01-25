################################################################################
## Title:         input.R
## Description:   Definition of file paths and workflow parameters
## Author:        Anna Niehues
## Date created:  2020-11-20
## Email:         anna.niehues@radboudumc.nl
################################################################################
## Notes:
## 
################################################################################

outdir_path <- here::here("output")

### definition of covariates and models
models_json <- rjson::fromJSON(file = "parameters.json")
# names(models_json) <- lapply(models_json, function(x) {x[["comparison"]]})
models_list <- list()
num_models <- 0
for (comparison in models_json) {
  comparison_variables <- unique(unlist(lapply(
    comparison[["models"]], function(x) {x})))
  for (model_name in names(comparison$models)) {
    for (cohort in comparison$cohorts) {
      parameter_set <- list(
        comparison = comparison$comparison,
        model = model_name,
        covariates = comparison$models[[model_name]],
        comparison_variables = comparison_variables,
        cohort = cohort
      )
      num_models <- num_models + 1
      models_list[[num_models]] <- parameter_set
    }
  }
}

# alpha level for p-value
alpha = 0.05

### input files - predicted values
# EpiSmokEr results
#"","SampleName","smokingScore","methylationScore","logOdds_CS","logOdds_FS",
#"logOdds_NS","probs_CS","probs_FS","probs_NS","PredictedSmokingStatus","reported"
epismoker_results_paths <- lapply(
  list("NTR" = "smoking_status_EpiSmokEr_NTR.csv",
       "LL" = "smoking_status_EpiSmokEr_LL.csv",
       "LLS" = "smoking_status_EpiSmokEr_LLS.csv",
       "RS" = "smoking_status_EpiSmokEr_RS.csv"),
  function(x) {file.path(
    here::here("~/researchdrive/aniehues/epismoker_210320"), x)})
# Metabolic scores from Daniele 
#"visit_id","bios_id","uuid",s_sex","s_diabetes","s_lipidmed",
#"s_blood_pressure_lowering_med","s_current_smoking","s_metabolic_syndrome",
#"s_alcohol_consumption","s_high_age","s_middle_age","s_low_age","s_obesity",
#"s_high_hscrp","s_high_triglycerides","s_high_ldl_chol","s_low_hdlchol",
#"s_high_totchol","s_low_eGFR","s_low_wbc","s_low_hgb"
metabolomics_surrogates_paths <- lapply(
  list("NTR" = "surrogates_VUNTR.csv",
       "LL" = "surrogates_LL.csv",
       "LLS" = "surrogates_LLS_PAROFF.csv",
       "RS" = "surrogates_RS.csv"),
  function(x) {file.path(
    here::here("~/researchdrive/RSC_BIOS/Users/dbizzarri/Surrogates_BIOS"), x)})

