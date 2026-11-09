# Repeated mixed-effects models for microbiome/clinical subgroup analyses.
# Keep numeric extraction positions consistent with the original manuscript code.

run_subgroup_associations <- function(response, predictor, output_file,
                                      change_data, baseline_data) {
  groups <- list(
    list(column = "sex", label = "sex", n = 2L, test = "coefficient"),
    list(column = "age_gp", label = "age group", n = 2L, test = "coefficient"),
    list(column = "diabmeds_gp", label = "# glucose-control medications used", n = 3L, test = "likelihood_ratio"),
    list(column = "cholmeds_gp", label = "# cholestrol medications used", n = 3L, test = "likelihood_ratio"),
    list(column = "htnmeds_gp", label = "# anti-hypertension medications used", n = 3L, test = "likelihood_ratio"),
    list(column = "metformin", label = "Metformin use", n = 2L, test = "coefficient"),
    list(column = "glp1", label = "GLP1 receptor agonist use", n = 2L, test = "coefficient"),
    list(column = "sglt2", label = "SGLT-2 inhibitors use", n = 2L, test = "coefficient")
  )

  fit_model <- function(group = NULL, include_interaction = TRUE) {
    term <- predictor
    if (!is.null(group)) {
      operator <- if (include_interaction) "*" else "+"
      term <- paste(predictor, operator, group)
    }
    formula <- stats::as.formula(paste(response, "~", term, "+ visit + (1|pptid)"))
    lmerTest::lmer(formula, data = change_data)
  }

  subgroup_results <- lapply(groups, function(spec) {
    model <- fit_model(spec$column)
    if (spec$test == "likelihood_ratio") {
      reduced_model <- fit_model(spec$column, include_interaction = FALSE)
      test_table <- as.data.frame(stats::anova(reduced_model, model))
      p_value <- test_table[2, 8]
    } else {
      p_value <- coef(summary(model))[7, 5]
    }

    trends <- as.data.frame(emmeans::emtrends(
      model,
      specs = stats::as.formula(paste0("~", spec$column)),
      var = predictor
    ))
    rows <- seq_len(spec$n)
    data.frame(
      Variable = rep(spec$label, spec$n),
      Group = trends[rows, 1],
      N = as.vector(table(baseline_data[[spec$column]])[rows]),
      Est = trends[rows, 2],
      lb = trends[rows, 5],
      ub = trends[rows, 6],
      p = rep(p_value, spec$n)
    )
  })

  marginal_model <- fit_model()
  marginal_coef <- coef(summary(marginal_model))
  marginal_ci <- stats::confint(marginal_model)
  subgroup_results[[length(subgroup_results) + 1L]] <- data.frame(
    Variable = "marginal", Group = NA, N = NA,
    Est = marginal_coef[2, 1],
    lb = marginal_ci[4, 1],
    ub = marginal_ci[4, 2],
    p = NA
  )

  result <- do.call(rbind, subgroup_results)
  utils::write.table(result, output_file, sep = "\t", quote = FALSE,
                     row.names = FALSE, na = "")
  result
}

format_subgroup_summary <- function(df) {
  dplyr::mutate(df, print = sprintf("%.2f(%.2f,%.2f)", Est, lb, ub))
}

# Baseline and change-from-baseline associations with the eight prespecified features.
summarize_feature_associations <- function(response, analysis, data, output_file,
                                           compact_prausnitzii_label = FALSE) {
  stopifnot(analysis %in% c("baseline", "change"))
  feature_names <- c("dash_index", "shannon", "butyrate_producer",
                     "Faecalibacterium.prausnitzii", "But_log", "EtfA_log",
                     "EtfB_log", "Acetyl_CoA_log")
  if (analysis == "baseline") {
    feature_labels <- c(
      "DASH4D microbe index", "Shannon", "log2(Butyrate producer)",
      if (compact_prausnitzii_label) "log2(F.prausnitzii)" else "log2(F. prausnitzii)",
      "log2(But gene per 1 million reads)", "log2(EtfA gene per 1 million reads)",
      "log2(EtfB gene per 1 million reads)", "log2(Acetyl-CoA pathway per 1 million reads)"
    )
  } else {
    feature_labels <- c("DASH4D microbe index", "Shannon", "log2(Butyrate producer)",
                        "log2(F. prausnitzii)", "But gene", "EtfA gene", "EtfB gene",
                        "Acetyl-CoA pathway")
  }

  results <- Map(function(predictor, feature_label) {
    covariates <- if (analysis == "baseline") {
      "sex + age_gp + diabmeds_gp + cholmeds_gp + htnmeds_gp + metformin + glp1 + sglt2"
    } else {
      "visit + (1|pptid)"
    }
    formula <- stats::as.formula(paste0(response, " ~ scale(", predictor, ") + ", covariates))
    model <- if (analysis == "baseline") stats::lm(formula, data = data) else lmerTest::lmer(formula, data = data)
    coefficients <- coef(summary(model))
    confidence_interval <- stats::confint(model)
    coefficient_row <- 2L
    ci_row <- if (analysis == "baseline") 2L else 4L
    p_column <- if (analysis == "baseline") 4L else 5L
    data.frame(Feature = feature_label,
               est = coefficients[coefficient_row, 1],
               lb = confidence_interval[ci_row, 1],
               ub = confidence_interval[ci_row, 2],
               pval = coefficients[coefficient_row, p_column])
  }, feature_names, feature_labels)

  result <- do.call(rbind, results)
  utils::write.table(result, output_file, sep = "\t", quote = FALSE,
                     row.names = FALSE, na = "")
  result
}
