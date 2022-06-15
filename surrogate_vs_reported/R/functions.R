################################################################################
## Title:         functions.R
## Description:   Definition of functions called in analysis plan
## Author:        Anna Niehues
## Date created:  2020-11-20
## Email:         anna.niehues@radboudumc.nl
################################################################################
## Notes:
##
################################################################################

#' Check if directory exists and create it if necessary.
#'
#' @param directory Path to directory.
#' @return Path to directory.
check_dir <- function(directory) {
  if (!dir.exists(directory)) {
    dir.create(directory, recursive = TRUE)
  }
  return(directory)
}


#' Load data from BBMRIomics
#'
#' @param data_set_name Name of BBMRIomics data set
load_data_set <- function(cohort) {
  require(BBMRIomics)
  message("[load data] ", cohort)
  tmpenv <- new.env()
  if (cohort == "RS") {
    BBMRIomics::bbmri.data(rnaSeqData_ReadCounts_RS_Freeze2_unrelated,
                           envir = tmpenv)
  } else if (cohort == "LL") {
    BBMRIomics::bbmri.data(rnaSeqData_ReadCounts_LL_Freeze2_unrelated,
                           envir = tmpenv)
  } else if (cohort == "LLS") {
    BBMRIomics::bbmri.data(rnaSeqData_ReadCounts_LLS_Freeze2_unrelated,
                           envir = tmpenv)
  } else if (cohort == "NTR") {
    BBMRIomics::bbmri.data(rnaSeqData_ReadCounts_NTR_Freeze2_unrelated,
                           envir = tmpenv)
  }
  counts <- tmpenv[["counts"]]
  return(counts)
}

# TODO create test data for running workflow when BBMRIomics is not available
create_test_data <- function(){
  require(SummarizedExperiment)
  nrows <- 100
  ncols <- 6
  counts <- matrix(runif(nrows * ncols, 1, 1e4), nrows)
  rowRanges <- GRanges(rep(c("chr1", "chr2"), c(50, 150)),
                       IRanges(floor(runif(nrows, 1e5, 1e6)), width = 100),
                       strand = sample(c("+", "-"), nrows, TRUE),
                       feature_id = sprintf("ID%03d", 1:nrows))
  colData <- DataFrame(Treatment = rep(c("ChIP", "Input"), 3),
                       row.names = LETTERS[1:6])

  SummarizedExperiment(assays = list(counts = counts),
                       rowRanges = rowRanges, colData = colData)
}


#'Read predicted values from csv files
#'"","SampleName","smokingScore","methylationScore","logOdds_CS","logOdds_FS","logOdds_NS",
#'"probs_CS","probs_FS","probs_NS","PredictedSmokingStatus","reported"
#'
read_predicted_values <- function(path_to_csv) {
  df <- read.csv(path_to_csv)
  if ("uuid" %in% colnames(df)) {
    df$ID <- as.character(df$uuid)
  } else if ("SampleName" %in% colnames(df)) {
    df$ID <- as.character(df$SampleName)
  }
  # remove sample that don't have a BIOS ID (uuid)
  df <- df[!is.na(df$ID), ]
  # make sure that IDs are unique
  stopifnot(length(df$ID) == length(unique(df$ID)))
  # set IDs as rownames
  rownames(df) <- df$ID
  # rename factor levels for PredictedSmokingStatus
  if ("PredictedSmokingStatus" %in% colnames(df)) {
    df$PredictedSmokingStatus <- plyr::revalue(
      df$PredictedSmokingStatus,
      c("Current Smoker" = "current smoker",
        "Former Smoker" = "former-smoker",
        "Never Smoker" = "non-smoker"))
  }
  return(df)
}


#' Combine additional metadata into single data frame
#'
combine_metadata <- function(values_df_list) {
  # samples present in all additional data sets
  samples_intersect <- Reduce(
    intersect, lapply(values_df_list, function(x) {rownames(x)}))
  # subset data frames
  values_df_list <- lapply(values_df_list, function(x) {x[samples_intersect,]})
  # check that rownames are the same
  stopifnot(rownames(values_df_list[[1]]) == rownames(values_df_list[[2]]))
  # merge
  combined_df <- do.call(cbind, values_df_list)
  return(combined_df)
}


#' Ordered quantile normalization of features
#'
#' @param col_data Data frame with features
normalize_covariates <- function(col_data) {
  require(bestNormalize)
  vars <- c(
    "s_low_eGFR",
    "triglycerides",
    "s_high_triglycerides",
    "totchol",
    "s_high_totchol",
    "ldlchol",
    "s_high_ldl_chol",
    "hdlchol",
    "s_low_hdlchol",
    "sex_male",
    "s_sex",
    "lipidmed_statins",
    "s_lipidmed",
    "hscrp",
    "s_high_hscrp",
    "s_diabetes",
    "s_metabolic_syndrome",
    "s_blood_pressure_lowering_med",
    "sampling_age",
    "s_high_age",
    "bmi",
    "s_obesity",
    "hgb",
    "s_low_hgb",
    "wbc",
    "s_low_wbc",
    "smoking_current",
    "probs_CS",
    "s_current_smoking",
    "s_alcohol_consumption")
  # ordered quantile normaliztion
  # col_data.orq <- apply(col_data[, vars], 2, function(x) {
  #   orderNorm(as.numeric(x))$x.t})
  # colnames(col_data.orq) <- sapply(vars, function(x) {paste0(x, "_orq")})
  # z-score transformation
  col_data.z <- apply(col_data[, vars], 2, function(x) {
    scale(as.numeric(x))})
  colnames(col_data.z) <- sapply(vars, function(x) {paste0(x, "_z")})
  # ln and z-score transformation
  # col_data.zln <- apply(col_data[, vars], 2, function(x) {
  #   scale(log(as.numeric(x)+1))}) # plus 1 to deal with zero values
  # colnames(col_data.zln) <- sapply(vars, function(x) {paste0(x, "_zln")})
  # return(do.call(cbind, list(col_data.orq, col_data.z, col_data.zln)))
  return(col_data.z)
}


#' Filter and normalize RNA-Seq data
#'
#' @param counts A summarizedExperiment object of RNA-Seq read counts
#' @return Filtered and normalized RNA-Seq read counts (summarizedExperiment)
#' @examples
#' filter_rna(counts)
filter_rna <- function(counts, variables) {
  require(SummarizedExperiment)
  require(edgeR)
  # ------------------------------------
  # Filtering: drop samples (RNA-Seq) with missing covariates
  # ------------------------------------
  ## determine samples with missing covariates
  message("[filter] missing covariates")
  samples2drop <- apply(as.data.frame(colData(counts))[, variables], 1, anyNA)
  ## drop samples
  counts <- counts[, !samples2drop]
  message("Dropped ", table(samples2drop)["TRUE"],
          " samples with missing covariates, ", table(samples2drop)["FALSE"],
          " samples remaining.")
  # ------------------------------------
  # Filtering: drop samples that have many missing features (>10%)
  # ------------------------------------
  message("[filter] samples with missing features")
  min_features <- 0.1 * dim(counts)[1]
  samples2keep <- colSums(is.na(assays(counts)$data)) <= min_features
  counts <- counts[, samples2keep]
  message("Dropped ", table(samples2keep)["FALSE"],
          " samples with >10% missing features, ", table(samples2keep)["TRUE"],
          " samples remaining.")
  # ------------------------------------
  # Filtering: drop lowly expressed genes
  # ------------------------------------
  min_count <- 1
  min_samples <- 0.1 * sum(samples2keep)
  features2keep <- rowSums(assays(counts)$data > min_count) > min_samples
  counts <- counts[features2keep,]
  message("Dropped ", table(features2keep)["FALSE"],
          " features with in count<", min_count, " in >10% of samples, ",
          table(features2keep)["TRUE"], " features remaining.")
  return(counts)
}


load_and_prepare <- function(cohort,
                             predicted_values,
                             nonBIOS_test = FALSE) {
  require(SummarizedExperiment)
  if (!nonBIOS_test) { # load BBMRI BIOS data
    counts <- load_data_set(cohort)
  } else { # dummy data to test workflow when not working in the BIOS VM
    counts <- create_test_data()
  }
  # make subset - remove samples with predicted values which not present in the omics data
  samples_not_predicted <- setdiff(rownames(colData(counts)),
                                   rownames(predicted_values))
  tmp <- predicted_values[samples_not_predicted,]
  rownames(tmp) <- samples_not_predicted
  predicted_values <- rbind(predicted_values, tmp)

  # add predicted feature values to counts data
  colData(counts) <- cbind(colData(counts),
                           predicted_values[rownames(colData(counts)),])
  # numeric smoking variable
  colData(counts)$smoking_current <- NA
  colData(counts)$smoking_current[colData(counts)$smoking == "current smoker"] <- 1
  colData(counts)$smoking_current[colData(counts)$smoking %in% c("former-smoker", "non-smoker")] <- 0
  colData(counts)$smoking_current <- as.numeric(colData(counts)$smoking_current)
  # numeric sex
  colData(counts)$sex_male <- NA
  colData(counts)$sex_male[colData(counts)$sex == "male"] <- 1
  colData(counts)$sex_male[colData(counts)$sex == "female"] <- 0
  colData(counts)$sex_male <- as.numeric(colData(counts)$sex_male)
  # numeric lipidmed
  colData(counts)$lipidmed_statins <- NA
  colData(counts)$lipidmed_statins[colData(counts)$lipidmed == "statins"] <- 1
  colData(counts)$lipidmed_statins[colData(counts)$lipidmed %in% c("no", "yes, but no statins")] <- 0
  # calculate BMI
  colData(counts)$bmi <- apply(colData(counts), 1, function(x) {
    as.numeric(x[["weight"]]) / (as.numeric(x[["height"]])/100)^2 })
  # normalize values
  colData(counts) <- cbind(colData(counts),
                           normalize_covariates(colData(counts)))
  return(counts)
}


#' Select subset of data where predicted and reported covariate values are similar
get_training_data <- function(metadata, covars2compare, cohorts) {
  require(bestNormalize)
  lapply(cohorts, function(cohort) {
    m <- metadata[metadata$biobank_id == cohort, unlist(covars2compare)]
    covars2compare.bin <- setNames(lapply(names(covars2compare), function(x) {
      binarize(m[, covars2compare[[x]]], location_measure = "median")$x.t}), names(covars2compare))
    covars2compare.bin$match <- covars2compare.bin$reported == covars2compare.bin$predicted
    covars2compare.bin <- as.data.frame(covars2compare.bin, row.names = rownames(m))
  })
}


### TWAS

#' Make formula based on response variable and covariates
#'
#' @param covariates A character vector defining all covariates.
#' @return Formula object with symbolic model formula.
#' @example
#' my_formula <- make_formula('sampling_age', c('sex', 'smoking'))
make_formula <- function(covariates) {
  model_formula <- as.formula(paste("~",
                                    paste(covariates,
                                          sep = "",
                                          collapse = "+"),
                                    sep = ""))
  return(model_formula)
}


#' Create design for summarizedExperiments object for given response variable
#' and covariates
#'
#' @param se A summarized experiments object.
#' @param covariates A character vector with covariates.
#' @return The model design matrix.
#' @example
#' my_design <- create_design(counts, 'sampling_age', c('sex', 'smoking'))
create_design <- function(se, covariates) {
  require(SummarizedExperiment)
  message("[TWAS] create design")
  design <- as.data.frame(colData(se))[, c(covariates), drop = FALSE]
  # design matrix for regression-like model with formula and data
  model_formula <- make_formula(covariates)
  design <- model.matrix(model_formula, data = design)
  return(design)
}


#' Normalize RNA-seq data (CPM)
#'
#' @param counts_se A summarizedExperiment object with read counts as log counts-per-million
#' @param design Model design matrix
normalize_rnaseq_data <- function(counts_se, design, plot_dir, tag) {
  require(edgeR)
  require(limma)
  require(SummarizedExperiment)
  require(ggplot2)
  ## convert counts to Digital Gene Expression data - class
  counts_dge <- DGEList(counts = assays(counts_se)$data)
  # ------------------------------------
  # Normalization of RNA-Seq read counts to account for different library sizes
  # ------------------------------------
  message("[normalize] trimmed mean of M-values (TMM)")
  counts_dge <- calcNormFactors(counts_dge) # trimmed mean of M-values (TMM)
  # stored in DGEList$samples$norm.factors
  # ------------------------------------
  # normalize using voom
  # ------------------------------------
  message("[TWAS] normalize")
  pdf(file.path(plot_dir, paste0(tag, "_voom_mean-var.pdf")))
  data_voom <- voom(counts_dge, design, plot = TRUE)
  dev.off()
  return(data_voom)
  # save mean-variance trand - voom plot
  # gg_v <- ggplot(data.frame(x = data_voom$voom.xy$x,
  #                           y = data_voom$voom.xy$y,
  #                           linex = data_voom$voom.line$x,
  #                           liney = data_voom$voom.line$y)) +
  #   geom_point(aes(x = x, y = y), size = 0.5) +
  #   geom_path(aes(x = linex, y = liney), color = "red")

  # ggsave(file.path(plot_dir, paste0(tag, "_voom_mean-var.png")), gg_v)
  # assumes normal distribution of logCPM values (may not be the case, see above)
}



#' TWAS - fit linear model for each gene (least squares) for a quantitative
#' trait (e.g. age)
#'
#' @param data_voom Voom normalized RNA-seq data (output from `normalize_rnaseq_data()`)
#' @param design Model design matrix
#' @return A list with results
#' @example
#' run_twas(counts_se, design)
run_twas <- function(data_voom, design) {
  require(limma)
  # ------------------------------------
  # Perform linear modeling
  # ------------------------------------
  message("[TWAS] fit model")
  fitLogCPM <- lmFit(data_voom, design)
  # ------------------------------------
  # Statistics
  # ------------------------------------
  message("[TWAS] empirical Bayes")
  efit <- eBayes(fitLogCPM) # Empirical Bayes Statistics for Differential Expression
  # NOTE: response variable (e.g. age) at position 2 in coefficients
  message("[TWAS] stats")
  results <- data.frame(eff.size = efit$coefficients[, 2],
                        std.err = sqrt(efit$s2.post) * efit$stdev.unscaled[, 2])
  results$tstat <- results$eff.size / results$std.err # t-statistics
  results$pval <- 2 * pnorm(-abs(results$tstat)) # ~efit$p.value
  results$z <- qnorm(1 - results$pval)
  message("[TWAS] bacon")
  results <- bacon_adjusted_TWAS(results)
  # adjust p-values for multiple testing
  results$padj <- p.adjust(results$pval, method = "fdr")
  results$padj.bacon <- p.adjust(results$pval.bacon, method = "fdr")
  # rank genes
  message("[rank genes]")
  results$ranking <- -log10(results$padj) * sign(results$eff.size)
  results$ranking.bacon <- -log10(results$padj.bacon) * sign(results$eff.size.bacon)
  # add gene symbols
  results <- add_gene_symbols(results)
  # clear memory
  gc()
  return(list(results = results, efit = efit))
}


#' Select significant genes based on given alpha, corrected for multiple testing
#' return results sorted by p-value (lowest first)
get_signif_genes <- function(se, alpha){
  signif <- se$results[se$results$padj < alpha,]
  signif_genes <- as.character(unlist(rownames(signif[order(signif$padj),])))
  return(signif_genes)
}


#' Perform correlation analysis of regression coefficients
#' #standardized coefficients (beta coeeficients)
#' Two or more models
#'
correlation_analysis <- function(betas, model_parameters, cohort, outdir) {
  require(reshape2)
  require(ggplot2)
  require(gtools)
  plot_cor <- function(model1_idx, model2_idx, betas,
                       model1, model2, cohort, outdir) {
    cor_m <- cor(betas[[model1_idx]], betas[[model2_idx]], method = "pearson")
    cor_m <- melt(cor_m)
    gg_p <- ggplot(data = cor_m) +
      geom_tile(aes(x = Var1, y = Var2, fill = value), color = "white") +
      scale_fill_gradient2(
        low = "red", high = "blue", mid = "white",
        midpoint = 0, limit = c(-1,1), space = "Lab",
        name="Pearson\nCorrelation") +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
      ggtitle(paste0("Correlation of regression coefficients, ", cohort)) +
      xlab(model1) +
      ylab(model2)
    ggsave(file.path(outdir,
                     paste0("beta_cor_", cohort, "_",
                            model1, "_", model2, ".png")), gg_p)
    return(gg_p)
  }
  plot_coef <- function(model1_idx, model2_idx, betas, covar_num,
                        model1, model2, cohort, outdir) {
    df <- as.data.frame(cbind(betas[[model1_idx]][, covar_num],
                              betas[[model2_idx]][, covar_num]))
    labels <- c(names(as.data.frame(betas[[model1_idx]]))[[covar_num]],
                names(as.data.frame(betas[[model2_idx]]))[[covar_num]])
    names(df) <- c("model1", "model2")
    gg_p <- ggplot(data = df, aes(x = model1, y = model2)) +
      geom_point() +
      geom_smooth(method = lm) +
      ggtitle(paste0("Regression coefficients, ", cohort)) +
      xlab(paste0(c(model1, labels[[1]]), collapse = ", ")) +
      ylab(paste0(c(model2, labels[[2]]), collapse = ", "))
    ggsave(file.path(outdir,
                     paste0("beta_coef_", cohort, "_",
                            model1, "_", model2, ".png")), gg_p)
    return(gg_p)
  }
  model_names <- unlist(lapply(model_parameters, function(x) x[["model"]]))
  print(model_names)
  if (length(model_names) > 1) {
    combos <- gtools::combinations(
      length(model_names), 2, seq_len(length(model_names)))
    my_plots <- list()
    for (i in seq(1, dim(combos)[1])) {
      model1_idx <- combos[i,1]
      model2_idx <- combos[i,2]
      model1 <- model_names[[model1_idx]]
      model2 <- model_names[[model2_idx]]
      my_plots <- c(
        my_plots,
        setNames(list(
          plot_cor(model1_idx, model2_idx, betas,
                   model1, model2, cohort, outdir)),
          c(paste0(c("cor", model1, model2), collapse = "_"))))
      covar_num <- 2
      my_plots <- c(
        my_plots,
        setNames(list(
          plot_coef(model1_idx, model2_idx, betas, covar_num,
                    model1, model2, cohort, outdir)),
          c(paste0(c("coef", covar_num, model1, model2), collapse = "_"))))
    }
  } else {my_plots = NULL}
  return(my_plots)
}


#' Bias and inflation of test statistics
#'
#' van Iterson, M., van Zwet, E.W., the BIOS Consortium. et al.
#' Controlling bias and inflation in epigenome- and transcriptome-wide
#' association studies using the empirical null distribution.
#' Genome Biol 18, 19 (2017). https://doi.org/10.1186/s13059-016-1131-9
#'
#'
bacon_adjusted_TWAS <- function(TWASresults) {
  require(bacon)
  bc <- bacon(teststatistics = NULL,
              effectsizes = TWASresults$eff.size,
              standarderrors = TWASresults$std.err)
  TWASresults$pval.bacon <- bacon::pval(bc, corrected = TRUE)
  TWASresults$eff.size.bacon = bacon::es(bc, corrected = TRUE)
  TWASresults$std.err.bacon = bacon::se(bc, corrected = TRUE)
  return(TWASresults)
}

calc_bias_inflation <- function(TWAS, comparison, cohort, model, alpha) {
  require(bacon)
  bc <- bacon(teststatistics = NULL,
              effectsizes = TWAS$results$eff.size,
              standarderrors = TWAS$results$std.err)

  results <- list(
    bias = bacon::bias(bc),
    inflation = bacon::inflation(bc),
    num_significant = sum(TWAS$results$padj < alpha),
    num_significant.bacon = sum(TWAS$results$padj.bacon < alpha),
    mean_abseffsize = mean(abs(TWAS$results$eff.size)),
    mean_abseffsize.bacon = mean(abs(TWAS$results$eff.size.bacon)),
    comparison = comparison,
    cohort = cohort,
    model = model)
  return(results)
}

# TODO currently not used for further analysis
select_model_aic <- function(counts_normalized){
  require(limma)
  #require(SummarizedExperiment)
  # samples <- Reduce(intersect, lapply(counts_normalized, function(x) {
    # rownames(x$design) }))
  v <- counts_normalized[[1]]$E
  designs <- lapply(counts_normalized, function(x) {x$design})
  aic_select <- selectModel(v, #assays(counts)$data,
                            designlist = designs,
                            criterion = "aic")
  return(aic_select)
}

aic_table <- function(aic_selection) {
  df_all <- do.call(rbind, lapply(names(aic_selection), function(x) {
    df <- as.data.frame(table(aic_selection[[x]]$pref))
    df$cohort <- c(x)
    print(df)
    }))
  colnames(df_all) <- c("Model", "Number.of.genes", "Cohort")
  df_all
}

#'Meta analysis
#' As performed in the X_omics summer school
#' 2019-07 Integrative X-omics Analyses Empowering Personalized Healthcare/BBMRIomics_Intro.html
#' and implemented in
#' https://www.bioconductor.org/packages/release/bioc/vignettes/bacon/inst/doc/bacon.html
#' see https://rdrr.io/bioc/bacon/src/R/BaconMethods.R
#'
meta_analysis <- function(results){

  message("[meta] reducing")
  results <- do.call("cbind", results)

  message("[meta] perform fixed-effect meta-analysis using bacon-corrected effect sizes")
  ES <- results[, grepl("eff.size.bacon", colnames(results))]
  SE <- results[, grepl("std.err.bacon", colnames(results))]
  W <- 1/SE^2
  V <- 1/rowSums(W)
  TS <- rowSums(ES*W)*V

  results$eff.size.meta <- TS
  results$std.err.meta <- sqrt(V)
  results$tstat.meta <- TS/sqrt(V)
  results$pval.meta <- 2*pnorm(-abs(results$tstat.meta))
  results$padj.meta <- p.adjust(results$pval.meta, method = "bonf")

  invisible(results[order(results$padj.meta),])
}


#' Perform meta analysis
#' Leave-one-cohort-out method as described in Rooij et al. 2019, genome biology, 20: 235
loo_analysis <- function(twas_results,
                         comparison_name,
                         model_name,
                         alpha = 0.05,
                         outdir) {
  if (length(twas_results) < 3) {
    return(warning("Require >= 3 cohorts for metaanalysis - skipping metaanalysis"))
  } else {
    stopifnot(length(twas_results) > 2)
    results <- list()
    idx <- 0
    # intersection of transcripts found in all cohorts
    transcript_intersect <- Reduce(
      intersect,
      lapply(twas_results, function(twas) {rownames(twas)}))
    twas_results <- lapply(twas_results, function(twas) {twas[transcript_intersect,]})
    pdf(file = file.path(outdir, paste0("bacon_metaanalysis_", model_name, ".pdf")))
    for (replication_cohort in names(twas_results)) {
      meta_cohorts <- names(twas_results)[
        names(twas_results) != replication_cohort]
      # metaanalysis of studies
      meta_results <- meta_analysis(twas_results[meta_cohorts])
        # lapply(twas_results[meta_cohorts],
        #        function(twas) {twas[transcript_intersect,]}))
      meta_results_signif <- rownames(
        meta_results[meta_results$padj.meta < alpha,])
      # replication
      replicated <- intersect(
        meta_results_signif,
        rownames(twas_results[[replication_cohort]][
          twas_results[[replication_cohort]]$padj.bacon < alpha,]))
      metaanalyzed <- length(meta_results_signif)
      replicated <- length(replicated)
      rankscore.cor <- cor(
        -log10(meta_results$padj.meta)*sign(meta_results$eff.size.meta),
        -log10(twas_results[[replication_cohort]]$padj.bacon
               )*sign(twas_results[[replication_cohort]]$eff.size.bacon),
                      method = "pearson")
      eff.size.cor <-cor(meta_results$eff.size.meta,
                     twas_results[[replication_cohort]]$eff.size.bacon,
                     method = "pearson")

      idx <- idx + 1
      results[[idx]] <- as.data.frame(list(
        meta.cohorts = paste(unlist(meta_cohorts), collapse = "_"),
        meta.analyzed = metaanalyzed,
        replicated = replicated,
        eff.size.cor = eff.size.cor,
        rankscore.cor = rankscore.cor,

        comparison = comparison_name,
        model = as.character(model_name),
        replication.cohort = replication_cohort))
      results[[idx]]$percentage_replicated <- replicated / metaanalyzed
    }
    dev.off()
    return(do.call(rbind, results))
  }
}


#' Volcano plot
#'
#'
plot_volcano <- function(TWAS, plot_title, plot_dir) {
  require(ggplot2)
  require(ggrepel)
  # add labels for top 25 genes
  TWAS$results$label <- NA
  # top_genes <- rownames(limma::topTable(TWAS$efit, number = 25, coef = 2))
  top_genes <- rownames(head(TWAS$results[order(TWAS$results$padj.bacon),], n = 25))
  # TWAS$results$tmporder <- abs(-log10(TWAS$results$padj)*TWAS$results$eff.size)
  TWAS$results$label[rownames(TWAS$results) %in% top_genes
                     ] <- TWAS$results$SYMBOL[rownames(TWAS$results) %in% top_genes]
  # volcano plot
  thresholds = c(0.05, 0.01)
  ggp <- ggplot() +
    geom_point(
      aes(x = eff.size.bacon, y = -log10(padj.bacon),
          col = abs(-log10(padj.bacon)*eff.size.bacon)),
      data = TWAS$results) +
    scale_colour_viridis_c(direction = -1) +
    # horizontal lines at different adjusted p-value thersholds
    geom_hline(aes(
      lty = sapply(thresholds, function(x) paste0("p = ", x)),
      yintercept = sapply(thresholds, function(x) -log10(x)))) +
    labs(lty = "threshold") +
    geom_text_repel(
      aes(x = eff.size.bacon, y = -log10(padj.bacon),
          col = abs(-log10(padj.bacon)*eff.size.bacon), label = label),
      data = TWAS$results) +
    ggtitle(plot_title)
  # save plot
  ggsave(file.path(
    plot_dir,
    paste0("TWAS_volcanoplot_", gsub("[,. ]+", "_", plot_title), ".png")), ggp)
  return(ggp)
}


add_gene_symbols <- function(TWASresults) {
  # add gene IDs
  message("[adding symbols]")
  library(org.Hs.eg.db)

  symbol <- select(org.Hs.eg.db, keys = rownames(TWASresults),
                   keytype = "ENSEMBL", # The ensembl ID as indicated by ensembl
                   columns = c("SYMBOL", # The official gene symbol
                               "ENTREZID" # Entrez gene Identifiers
                   ))
  TWASresults$ENTREZID <- symbol$ENTREZID[!duplicated(symbol$ENSEMBL)]
  TWASresults$SYMBOL <- symbol$SYMBOL[!duplicated(symbol$ENSEMBL)]
  return(TWASresults)
}


perform_gsea <- function(TWASresults, metaanalysis = FALSE){
  # compute serially; otherwise error with BiocParallel in fgsea
  library(BiocParallel)
  register(MulticoreParam(1, log = TRUE))
  if("dplyr" %in% (.packages())){
    detach("package:dplyr", unload = TRUE) # causes problems with fgsea
  }
  require(fgsea) # https://bioconductor.org/packages/release/bioc/html/fgsea.html
  require(reactome.db)
  # for testing
  # TWAS <- readd(TWAS_89L_hdlchol_z_s_low_hdlchol_z_NTR)
  # ranks <- TWAS$results[genes, ]$ranking.bacon # from rank_genes()
  # note that using adjusted p-values for ranking yields zero values and causes fgsea to hang
  if (metaanalysis == FALSE) { # data from TWAS
    TWASresults$rankstats <- -log10(TWASresults$pval.bacon) * sign(TWASresults$eff.size.bacon)
  } else { # data from meta analysis of TWAS
    TWASresults$rankstats <- -log10(TWASresults$pval.meta) * sign(TWASresults$eff.size.meta)
    TWASresults$ENTREZID <- TWASresults[, grepl("ENTREZID", colnames(TWASresults))][,1]
  }
  # some -log10 pvalues are Inf which causes error in fgsea; use max. rank value +1 instead
  TWASresults$rankstats[is.infinite(TWASresults$rankstats)] <- max(
    TWASresults$rankstats[is.finite(TWASresults$rankstats)]) + 1
  # filter genes by presence of ENTREZID
  genes <- !is.na(TWASresults$ENTREZID) & !is.na(TWASresults$rankstats)
  # rank stats
  ranks <- TWASresults[genes, ]$rankstats
  names(ranks) <- TWASresults[genes, ]$ENTREZID
  # get pathways from reactome for ranked genes
  pathways <- reactomePathways(names(ranks))
  message("[GSEA] pathway enrichment")
  set.seed(123)
  fgsea_result <- fgsea::fgseaMultilevel(
    pathways = pathways,
    stats = ranks,
    minSize = 15,
    maxSize = Inf,
    eps = 1e-50, # boundary for p-value calculation
    nproc = 1, # number of workers; note that setting this to 1 might make register(MulticoreParam(1, log = TRUE)) unnecessary
    nPermSimple = 10000)
   # in some cases (GSEA_89L_hdlchol_z_s_low_hdlchol_z_NTR) fgsea hangs; see also https://github.com/ctlab/fgsea/issues/57
  # z <- replicate(1e5, calcGseaStat(ranks, sample.int(length(ranks), 500)))
  # hist(z) # check for skewing
  # summary(ranks)
  result <- list(fgsea = fgsea_result,
                 pathways = pathways,
                 ranks = ranks)
  gc()
  return(result)
}


write_parameters_table <- function(models_list,
                                   epismoker_variables,
                                   surrogates_variables,
                                   outfile) {
  get_covariatetype <- function(outcome_variable) {
    if (outcome_variable %in% epismoker_variables) {
      covariate_type <- "Epigenetic surrogate"
    } else if (outcome_variable %in% surrogates_variables) {
      covariate_type <- "Metabolic surrogate"
    } else {
      covariate_type <- "Reported variable"
    }
    return(covariate_type)
  }
  rows <- lapply(models_list, function(x) {
    row <- list(
      `Comparison name` = x$comparison,
      `Model name` = x$model,
      `Outcome variable type` = get_covariatetype(x$model),
      `Model covariates` = paste0(x$covariates, collapse = ", "),
      Cohort = x$cohort)})
  df <- data.table::rbindlist(rows)
  write.csv(df, file = outfile)
}
