
<!-- README.md is generated from README.Rmd. Please edit that file -->

# forestplotutil

<!-- badges: start -->
<!-- badges: end -->

forestplotutil helps generate a data frame with subgroup analyses using
binomial logistic regression models. This data frame is ready for input
into `forestplot::forestplot`.

## Installation

You can install the development version of forestplotutil from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("LukasDSauer/forestplotutil")
```

## Example

This is a basic example on how to generate a subgroup analysis plot
using forestplotutil and forestplot.

``` r
library(forestplotutil)
library(forestplot)
#> Lade nötiges Paket: grid
#> Lade nötiges Paket: checkmate
#> Lade nötiges Paket: abind
library(dplyr)
#> 
#> Attache Paket: 'dplyr'
#> Die folgenden Objekte sind maskiert von 'package:stats':
#> 
#>     filter, lag
#> Die folgenden Objekte sind maskiert von 'package:base':
#> 
#>     intersect, setdiff, setequal, union

# Data preparation
data(knee, package="catdata")

knee %>% mutate(Th = factor(Th),
                Sex = factor(Sex, levels = c(0, 1), labels = c("male", "female")),
                Age = cut(Age, breaks = c(min(Age), 30, 40, max(Age)), include.lowest = TRUE),
                R4_improved = factor(if_else(R4 < R1, 1, 0))) ->
  knee
formula_overall = "R4_improved ~ Th + R1"

# Generate table
knee %>%
  subgroup_table_binomial(data = .,
                          formula_overall,
                          formulas = list(
                            Sex = "R4_improved ~ Th*Sex + R1",
                            Age = "R4_improved ~ Th*Age + R1"
                          ),
                          formulas_subgroups = list(
                            Sex = formula_overall,
                            Age = formula_overall
                          ),
                          group = "Th",
                          outcome = "R4_improved",
                          outcome_positive = 1) ->
  res

# Add a header
header <- tibble(subgroup = c("", "", "Subgroup"),
                 group1 = c("Improvement", "after ten days", "of placebo"),
                 group2 = c("Improvement", "after ten days", "of treatment"),
                 OR = c("", "Odds ratio", "[95% CI]"),
                 summary = TRUE,
                 p = c("", "", "p"))
res <- bind_rows(header,
                 res)

# Visualize using forestplot
library(forestplot)
res %>%
  forestplot(labeltext = c(subgroup, group1, group2, OR, p),
             align = c("l", "r", "r", "r", "c"),
             is.summary = summary,
             graph.pos = ncol(.) - 4,
             hrzl_lines = TRUE,
             clip = c(1/sqrt(2)^7, 64),
             xlab = "Odds ratio (95% CI)",
             graphwidth = unit(50, "mm"),
             colgap = unit(2, "mm"),
             lineheight = unit(4.5, "mm"),
             line.margin = unit(0.7, "mm"),
             xlog = TRUE,
             col = fpColors(box = "royalblue",
                            line = "darkblue",
                            summary = "royalblue"),
             ci.vertices = TRUE,
             txt_gp = fpTxtGp(summary = gpar(fontsize = 8.7, cex=1),
                              label = gpar(fontsize = 7.7, cex=1),
                              xlab = gpar(fontsize = 7.7, cex=0.8),
                              ticks =gpar(fontsize = 7.7, cex=0.8)))
```

<img src="man/figures/README-example-1.png" width="100%" />
