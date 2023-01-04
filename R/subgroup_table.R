#' Generate a tibble meant for a forest plot describing subgroup analyses
#'
#' This function performs data preparation for subgroup analyses using a forest
#' plot. The goal is to describe a binary outcome in two treatment groups
#' stratified by pre-specified subgroups with three different data:
#'  * absolute and relative frequencies of favorable outcomes in each subgroup
#'   stratum grouped by treatment group,
#'  * odds ratios and 95% confidence intervals of the treatment group effect
#'  calculated on each subgroup stratum separately by binary logistic regression
#'  models (`glm(family=binomial)`) of the outcome, and
#'  * a p-value of the interaction between the subgroup and the treatment group
#'  as calculated by a type III ANOVA on the binary logistic regression model
#'  on the full study population.
#'
#' @param data a data frame containing individual data of the study population.
#' @param formula_overall a string specifying a formula for a binary logistic.
#' regression model of the overall treatment effect (usually without subgroup covariates).
#' This overall formula is only necessary if `OR_interact = FALSE`.
#' @param formulas a named list of strings specifying formulas for binary logistic
#' regression models of the overall treatment effect adjusted for subgroup levels
#' including an interaction term between treatment group and subgroup variable.
#' These models are used for calculating the p-value of the interaction. The list
#' entries must be named by the subgroup variable names.
#' @param formulas_subgroups a named list of strings specifying formulas for binary logistic
#' regression models of the overall treatment effect (usually without subgroup covariates).
#' These models are used for calculating the odds ratios and its 95% confidence
#' intervals within each subgroup stratum. The list entries must be named by the
#' subgroup variable names. This parameter is only necessary if
#' `OR_interact = false`.
#' @param group a string specifying the name of the treatment group variable in
#' `data`. The group variable must be a factor.
#' @param outcome a string specifying the name of the outcome variable in `data`.
#' The outcome variable must be a factor.
#' @param outcome_positive a string specifying the name favourable outcome level.
#' @param OR_interact a string specifying whether you want to calculate the odds
#' ratio of the treatment groups effects in the subgroup strata (`OR_interact = FALSE` with
#' models for odds ratio supplied by `formulas_subgroups`) or the odds ratios
#' ratio of the interaction terms (`OR_interact = TRUE` with
#' models for odds ratio supplied by `formulas`).
#' @return A tibble ready for input into `forestplot::forestplot()`.
#'
#' @export
#' @examples
#' library(dplyr)
#'
#' # Data preparation
#' data(knee, package="catdata")
#'
#' knee %>% mutate(Th = factor(Th),
#'                 Sex = factor(Sex, levels = c(0, 1), labels = c("male", "female")),
#'                 Age = cut(Age, breaks = c(min(Age), 30, 40, max(Age)), include.lowest = TRUE),
#'                 R4_improved = factor(if_else(R4 < R1, 1, 0), levels = c(0, 1), labels = c("no", "yes"))) ->
#'   knee
#' formula_overall = "R4_improved ~ Th + R1"
#'
#' # Generate table
#' knee %>%
#'   subgroup_table_binomial(data = .,
#'                           formula_overall,
#'                           formulas = list(
#'                             Sex = "R4_improved ~ Th*Sex + R1",
#'                             Age = "R4_improved ~ Th*Age + R1"
#'                           ),
#'                           formulas_subgroups = list(
#'                             Sex = formula_overall,
#'                             Age = formula_overall
#'                           ),
#'                           group = "Th",
#'                           outcome = "R4_improved",
#'                           outcome_positive = "yes") ->
#'   res
#'
#' # Add a header
#' header <- tibble(subgroup = c("", "", "Subgroup"),
#'                  group1 = c("Improvement", "after ten days", "of placebo"),
#'                  group2 = c("Improvement", "after ten days", "of treatment"),
#'                  OR = c("", "Odds ratio", "[95% CI]"),
#'                  summary = TRUE,
#'                  p = c("", "", "p"))
#' res <- bind_rows(header,
#'                  res)
#'
#' # Visualize using forestplot
#' library(forestplot)
#' res %>%
#'   forestplot(labeltext = c(subgroup, group1, group2, OR, p),
#'              align = c("l", "r", "r", "r", "c"),
#'              is.summary = summary,
#'              graph.pos = ncol(.) - 4,
#'              hrzl_lines = TRUE,
#'              clip = c(1/sqrt(2)^7, 64),
#'              xlab = "Odds ratio (95% CI)",
#'              graphwidth = unit(50, "mm"),
#'              colgap = unit(2, "mm"),
#'              lineheight = unit(4.5, "mm"),
#'              line.margin = unit(0.7, "mm"),
#'              xlog = TRUE,
#'              col = fpColors(box = "royalblue",
#'                             line = "darkblue",
#'                             summary = "royalblue"),
#'              ci.vertices = TRUE,
#'              txt_gp = fpTxtGp(summary = gpar(fontsize = 8.7, cex=1),
#'                               label = gpar(fontsize = 7.7, cex=1),
#'                               xlab = gpar(fontsize = 7.7, cex=0.8),
#'                               ticks =gpar(fontsize = 7.7, cex=0.8))) ->
#'   plot
#' plot
#' ## Not run:
#' library(Cairo)
#' Cairo::CairoPNG("knee.png", width = 500, height = 300, units="pt", dpi=310)
#' plot
#' dev.off()
#' ## End(**Not run**)
subgroup_table_binomial <- function(data,
                           formula_overall = NULL,
                           formulas,
                           formulas_subgroups,
                           group,
                           outcome,
                           outcome_positive,
                           OR_interact = FALSE){
  # which level is the positive_outcome?
  outcome_positive_pos = which(levels(data[[outcome]]) == outcome_positive)
  if(outcome_positive_pos != 2){
    warning("The positive outcome is the first level of the factor, not
             the second. This may lead to unintuitive odds ratios. We recommend
             changing the order of the outcome levels.")
  }
  # Variables names of subgroup variables
  subgroup_vars <- names(formulas)
  # Labels for subgroup variables
  subgroup_labs <- character(length = length(subgroup_vars))
  subgroup_labs = subgroup_vars
  names(subgroup_labs) <- subgroup_vars
  for(subgroup in subgroup_vars){
    lab = attr(data[[subgroup]], "label")
    if(!is.null(lab)){
      subgroup_labs[[subgroup]] = lab
    }
  }
  # Treatment group levels
  if(!is.factor(data[[group]])){
    stop(paste("The group variable", group, "must be a factor variable."))
  }
  group_levels = levels(data[[group]])
  # Name of the second group
  group2_name = paste0(group, group_levels[2])
  # Data frame with results
  res <- data.frame()
  # Adding results for overall treatment effect if the odds ratios of treatment
  # effects in subgroups are shown
  if(!OR_interact){
    # Logistic regression model for overall treatment effect
    mod_overall = glm(formula = formula_overall,
                      family = binomial,
                      data=data)
    # Type III ANOVA is not necessary as no interaction terms are given
    p_str = rd_p(summary(mod_overall)$coefficients[group2_name,"Pr(>|z|)"])
    # Frequency table of positive outcome
    group_t = table(data[[group]])
    outcome_per_group_t =
      table(data[[group]], data[[outcome]])[, outcome_positive_pos]
    outcome_per_group_freq_t =
      matrix(paste0(outcome_per_group_t, "/",
                    group_t, " (",
                    rd_pct(outcome_per_group_t/group_t), "%)"),
             dim(group_t))
    # Confidence interval
    ci = exp(confint.default(mod_overall))
    # Odds ratio of treatment group in the model calculated on the subgroup
    OR_num = unname(exp(mod_overall$coefficients[group2_name]))
    OR_string = paste0(rd(exp(mod_overall$coefficients[group2_name])),
                       " [",rd(ci[group2_name,1]), ",", rd(ci[group2_name,2]), "]")
    # Combine results
    res_overall <- tibble(mean  =  OR_num,
                          lower =  ci[group2_name,1],
                          upper =  ci[group2_name,2],
                          subgroup = "Overall",
                          group1 = outcome_per_group_freq_t[1,],
                          group2 = outcome_per_group_freq_t[2,],
                          OR = OR_string,
                          summary = FALSE,
                          p = p_str
    )
    res <- bind_rows(res,
                     res_overall)

  }
  # Adding results for single subgroups
  for(subgroup in subgroup_vars){
    if(!is.factor(data[[subgroup]])){
      stop(paste("The subgroup variable", subgroup, "must be a factor variable."))
    }
    # Contrasts for mod_interact
    contrasts = list(subgroup=contr.sum, group=contr.sum)
    names(contrasts) = c(subgroup, group)
    # Logistic regression model with interaction term
    mod_interact = glm(formula = formulas[subgroup][[1]],
                       family = binomial,
                       contrasts=contrasts,
                       data=data)
    # Subgroup levels
    subgroup_levels = levels(data[[subgroup[[1]]]])
    # p-value is calculated by type III ANOVA
    aov_mod = car::Anova(mod_interact,
                    type=3)
    interact_Anova = paste0(subgroup, ":", group)
    p = aov_mod[interact_Anova,"Pr(>Chisq)"]
    # The name of the interaction term may be interchanged.
    if(is.na(p)){
      interact_Anova = paste0(group, ":", subgroup)
      p = aov_mod[interact_Anova,"Pr(>Chisq)"]
    }
    p_str = rd_p(p)
    # Frequency table of positive outcome
    subgroup_t = table(data[[subgroup]], data[[group]])
    if(0 %in% rowSums(subgroup_t)){
      stop(paste("The factor variable", subgroup, "contains unused levels. Consider using droplevels()."))
    }
    outcome_per_subgroup_t =
      table(data[[subgroup]], data[[group]], data[[outcome]])[, , outcome_positive_pos]
    outcome_per_subgroup_freq_t =
      matrix(paste0(outcome_per_subgroup_t, "/",
                    subgroup_t, " (",
                    rd_pct(outcome_per_subgroup_t/subgroup_t), "%)"),
             dim(subgroup_t))
    # Confidence intervals and ORs are either calculated as interaction terms
    # of the ANOVA model or as the treatment effects of separate models calculated
    # on the specific subgroups
    OR_num = numeric(length(subgroup_levels))
    ci_lower = numeric(length(subgroup_levels))
    ci_upper = numeric(length(subgroup_levels))
    OR_string = character(length(subgroup_levels))

    if(OR_interact){
      # Name of interaction term
      interact_mod = paste0(subgroup, 1:(length(subgroup_levels)-1), ":", group, "1")
      # Confidence interval(s)
      ci = exp(confint.default(mod_interact))
      ci_lower = unname(c(NA, ci[interact_mod,1]))
      ci_upper = unname(c(NA, ci[interact_mod,2]))
      # Odds ratio(s) of interaction term(s)
      OR_num = unname(c(1, exp(mod_interact$coefficients[interact_mod])))
      OR_string = as.character(rd(c(1, exp(mod_interact$coefficients[interact_mod]))))
    } else{
      k = 0
      for(level in subgroup_levels){
        k = k + 1
        # Calculate model on subgroup
        mod_subgroup = glm(formula = formulas_subgroups[subgroup][[1]],
                           family = binomial,
                           data=data[data[[subgroup]] == level & !is.na(data[[subgroup]]), ])
        # Confidence interval(s)
        ci = exp(confint.default(mod_subgroup))
        ci_lower[k] = ci[group2_name,1]
        ci_upper[k] = ci[group2_name,2]
        # Odds ratio of treatment group in the model calculated on the subgroup
        OR_num[k] = unname(exp(mod_subgroup$coefficients[group2_name]))
        OR_string[k] = paste0(rd(exp(mod_subgroup$coefficients[group2_name])),
                              " [",rd(ci[group2_name,1]), ",", rd(ci[group2_name,2]), "]")
      }
    }

    # Combine results
    res_subgroup <- tibble(mean  =  c(NA, OR_num),
                           lower =  c(NA, ci_lower),
                           upper =  c(NA, ci_upper),
                           subgroup = c(subgroup_labs[[subgroup]], paste(" ", subgroup_levels)),
                           group1 = c("", outcome_per_subgroup_freq_t[,1]),
                           group2 = c("", outcome_per_subgroup_freq_t[,2]),
                           OR = c("", OR_string),
                           summary = FALSE,
                           p = c(p_str, rep("", length(subgroup_levels)))
    )
    res <- bind_rows(res,
                     res_subgroup)
  }
  return(res)
}
