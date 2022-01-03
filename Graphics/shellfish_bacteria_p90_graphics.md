Revised Graphics Summarizing Shellfish Bacteria Data
================
Curtis C. Bohlen, Casco Bay Estuary Partnership.
02/17/2021

-   [Introduction](#introduction)
-   [Relevant Standards](#relevant-standards)
-   [Load Libraries](#load-libraries)
-   [Load Data](#load-data)
    -   [Main Data](#main-data)
        -   [Address Censored Data](#address-censored-data)
        -   [Remove NAs](#remove-nas)
        -   [Remove Sites not in Region](#remove-sites-not-in-region)
    -   [Weather Data](#weather-data)
        -   [Incorporate Weather Data](#incorporate-weather-data)
    -   [Summary Statistics Dataframe](#summary-statistics-dataframe)
-   [Critical levels (Reminder)](#critical-levels-reminder)
-   [Graphics for 238 Stations](#graphics-for-238-stations)
    -   [Bootstrap Confidence Interval
        Function](#bootstrap-confidence-interval-function)
    -   [Confidence Intervals for the Geometric
        Mean](#confidence-intervals-for-the-geometric-mean)
    -   [Confidence Intervals for the 90th
        Percentile](#confidence-intervals-for-the-90th-percentile)
    -   [Combining Results](#combining-results)
    -   [Geometric Mean Plot](#geometric-mean-plot)
    -   [p90 Plot](#p90-plot)
-   [Combined Graphics; Horizontal
    Layout](#combined-graphics-horizontal-layout)
    -   [Assemble Long Data](#assemble-long-data)
    -   [Base Plot](#base-plot)
        -   [Mimicing Graphic As Modified by the Graphic
            Designer](#mimicing-graphic-as-modified-by-the-graphic-designer)
        -   [Alternate Annotations](#alternate-annotations)
    -   [Shapes with outlines Plot](#shapes-with-outlines-plot)
        -   [Alternate Annotations](#alternate-annotations-1)
        -   [Annotation Dataframe](#annotation-dataframe)
-   [Combined Graphics: Vertical
    Layout](#combined-graphics-vertical-layout)
    -   [One Panel](#one-panel)
    -   [Two Panels](#two-panels)
-   [Years Plot](#years-plot)
-   [Frequency of Exceedences](#frequency-of-exceedences)

<img
    src="https://www.cascobayestuary.org/wp-content/uploads/2014/04/logo_sm.jpg"
    style="position:absolute;top:10px;right:50px;" />

# Introduction

Exploratory analysis highlights the extreme skewness of the distribution
of bacteria data, both here with the shellfish-related data collected by
DMR, and with the data related to recreational beaches, collected by
towns, and managed by DEP.

Here we present graphical summaries of the data, with an emphasis on
reporting observed quantities that relate to regulatory thresholds,
especially geometric means and 90th percentiles related to whether sites
can be classified as “Approved” for harvest.

The motivation for this notebook came from comments from reviewers of
our draft State of Casco Bay chapter that our rosy analysis of shellfish
bacteria data held for geometric mean thresholds, but less so for the
“p90” or 90th percentile thresholds. That means we needed to find ways
to depict both geometric mean and p90 values in at least some graphics.

# Relevant Standards

| Growing Area Classification | Activity Allowed                                                          | Geometric mean FC/100ml | 90th Percentile (P90) FC/100ml |
|-----------------------------|---------------------------------------------------------------------------|-------------------------|--------------------------------|
| Approved                    | Harvesting allowed                                                        | ≤ 14                    | ≤ 31                           |
| Conditionally Approved      | Harvesting allowed except during specified conditions                     | ≤ 14 in open status     | ≤ 31 in open status            |
| Restricted                  | Depuration harvesting or relay only                                       | ≤ 88 and &gt;15         | ≤ 163 and &gt;31               |
| Conditionally Restricted    | Depuration harvesting or relay allowed except during specified conditions | ≤ 88 in open status     | ≤ 163 in open status           |
| Prohibited                  | Aquaculture seed production only                                          | &gt;88                  | &gt;163                        |

# Load Libraries

``` r
library(readr)
library(tidyverse)      # Loads another `select()`
#> -- Attaching packages --------------------------------------- tidyverse 1.3.1 --
#> v ggplot2 3.3.5     v dplyr   1.0.7
#> v tibble  3.1.6     v stringr 1.4.0
#> v tidyr   1.1.4     v forcats 0.5.1
#> v purrr   0.3.4
#> -- Conflicts ------------------------------------------ tidyverse_conflicts() --
#> x dplyr::filter() masks stats::filter()
#> x dplyr::lag()    masks stats::lag()

library(emmeans)        # For marginal means
library(mgcv)           # For GAMs, here used principally for hierarchical models
#> Loading required package: nlme
#> 
#> Attaching package: 'nlme'
#> The following object is masked from 'package:dplyr':
#> 
#>     collapse
#> This is mgcv 1.8-38. For overview type 'help("mgcv-package")'.

library(CBEPgraphics)
load_cbep_fonts()
theme_set(theme_cbep())

library(LCensMeans)
```

# Load Data

## Main Data

``` r
sibfldnm <- 'Data'
parent <- dirname(getwd())
sibling <- file.path(parent,sibfldnm)

dir.create(file.path(getwd(), 'figures'), showWarnings = FALSE)
```

``` r
fl1<- "Shellfish data 2015 2018.csv"
path <- file.path(sibling, fl1)

coli_data <- read_csv(path, 
    col_types = cols(SDate = col_date(format = "%Y-%m-%d"), 
        SDateTime = col_datetime(format = "%Y-%m-%dT%H:%M:%SZ"), # Note Format!
        STime = col_time(format = "%H:%M:%S"))) %>%
  mutate_at(c(6:7), factor) %>%
  mutate(Class = factor(Class, levels = c( 'A', 'CA', 'CR',
                                           'R', 'P', 'X' ))) %>%
  mutate(Tide = factor(Tide, levels = c("L", "LF", "F", "HF",
                                        "H", "HE", "E", "LE"))) %>%
  mutate(DOY = as.numeric(format(SDate, format = '%j')),
         Month = as.numeric(format(SDate, format = '%m'))) %>%
  mutate(Month = factor(Month, levels = 1:12, labels = month.abb))
```

### Address Censored Data

We calculate a estimated conditional mean to replace the (left) censored
values. The algorithm is not entirely appropriate, as it assumes
lognormal distribution, and our data are closer to Pareto-distributed.
Still, it handles the non-detects on a more rational basis than the
usual conventions.

Second, we calculate a version of the data where non-detects are
replaced by half the value of the detection limit. However, we plan to
use the LOG of coliform counts in Gamma GLM models, which require
response variables to be strictly positive. The most common Reporting
Limit in these data is `RL == 2`. Half of that is 1.0, and
`log(1.0) == 0`. Consequently, we replace all values 1.0 with 1.1, as
log(1.1) is positive, and thus can be modeled by a suitable `gam()`
model.

In this notebook, we do not use either of these altered versions of the
data, but instead rely only on results calculated assuming “non-detects”
equal the reporting limit. This (obviously) over estimates the true
(unobserved) values of the non-detects, but it is consistent,
convenient, and transparent.

``` r
coli_data <- coli_data %>%
  mutate(ColiVal_ml = sub_cmeans(ColiVal, LCFlag)) %>%
  mutate(ColiVal_hf = if_else(LCFlag, ColiVal/2, ColiVal),
         ColiVal_hf = if_else(ColiVal_hf == 1, 1.1, ColiVal_hf))
```

### Remove NAs

``` r
coli_data <- coli_data %>%
  filter (! is.na(ColiVal))
```

### Remove Sites not in Region

We have some data that was selected for stations outside of Casco Bay.
To be  
careful, we remove sampling data for any site in th two adjacent Growing
Areas, “WH” and “WM”.

``` r
coli_data <- coli_data %>%
  filter(GROW_AREA != 'WH' & GROW_AREA != "WM") %>%
  mutate(GROW_AREA = fct_drop(GROW_AREA),
         Station = factor(Station))
```

## Weather Data

``` r
sibfldnm    <- 'Data'
parent      <- dirname(getwd())
sibling     <- file.path(parent,sibfldnm)

fn <- "Portland_Jetport_2015-2019.csv"
fpath <- file.path(sibling, fn)

weather_data <- read_csv(fpath, 
 col_types = cols(station = col_skip())) %>%
  select( ! starts_with('W')) %>%
  rename(sdate = date,
         Precip=PRCP,
         MaxT = TMAX,
         MinT= TMIN,
         AvgT = TAVG,
         Snow = SNOW,
         SnowD = SNWD) %>%
  mutate(sdate = as.Date(sdate, format = '%m/%d/%Y'))
```

``` r
weather_data <- weather_data %>%
  arrange(sdate) %>%
  
  select(sdate, Precip, AvgT, MaxT) %>%
  mutate(AvgT = AvgT / 10,
         MaxT = MaxT / 10,
         Precip = Precip / 10,
         Precip_d1 = dplyr::lag(Precip,1),
         Precip_d2 = dplyr::lag(Precip,2),
         log1precip    = log1p(Precip), 
         log1precip_d1 = log1p(Precip_d1),
         log1precip_d2 = log1p(Precip_d2),
         log1precip_2   = log1p(Precip_d1 + Precip_d2),
         log1precip_3   = log1p(Precip + Precip_d1 + Precip_d2)) %>%
  rename_with(tolower)
```

### Incorporate Weather Data

``` r
coli_data <- coli_data %>%
  left_join(weather_data, by = c('SDate' = 'sdate'))
```

## Summary Statistics Dataframe

We only work with the raw observed summary statistics here.

``` r
sum_data <- coli_data %>%
  mutate(logcoli = log(ColiVal),
         logcoli2 = log(ColiVal_ml)) %>%
  group_by(Station) %>%
  summarize(mean1 = mean(ColiVal),
            median1 = median(ColiVal),
            iqr1 = IQR(ColiVal),
            p901 = quantile(ColiVal, 0.9),
            meanlog1 = mean(logcoli, na.rm = TRUE),
            sdlog1 = sd(logcoli, na.rm = TRUE),
            nlog1 = sum(! is.na(logcoli)),
            selog1 = sdlog1/sqrt(nlog1),
            gmean1 = exp(meanlog1),
            U_CI1 = exp(meanlog1 + 1.96 * selog1),
            L_CI1 = exp(meanlog1 - 1.96 * selog1)) %>%
  mutate(Station = fct_reorder(Station, gmean1))
```

# Critical levels (Reminder)

Geometric Mean include:  
 &lt;  = 14 and  &lt;  = 88  
and for the p90  
 &lt; 31 and  &lt;  = 16

# Graphics for 238 Stations

## Bootstrap Confidence Interval Function

This is a general function, so needs to be passed log transformed data
to produce geometric means. note that we can substitute another function
for the mean to generate confidence intervals for other statistics, here
the 90th percentile.

``` r
boot_one <- function (dat, fun = "mean", sz = 1000, width = 0.95) {
  
  low <- (1 - width)/2
  hi <- 1 - low

  vals <- numeric(sz)
  for (i in 1:sz) {
    vals[i] <- eval(call(fun, sample(dat, length(dat), replace = TRUE)))
  }
  return (quantile(vals, probs = c(low, hi)))
}
```

## Confidence Intervals for the Geometric Mean

We need to first calculate confidence intervals on a log scale, then
build a tibble and back transform them.

``` r
gm_ci <- tapply(log(coli_data$ColiVal), coli_data$Station, boot_one)

# Convert to data frame (and then tibble...)
# This is convenient because as_tibble() drops the row names,
# which we want to keep.
gm_ci <- as.data.frame(do.call(rbind, gm_ci)) %>%
  rename(gm_lower1 = `2.5%`, gm_upper1 = `97.5%`) %>%
  rownames_to_column('Station')

# Back Transform
gm_ci <- gm_ci %>%
  mutate(gm_lower1 = exp(gm_lower1),
         gm_upper1 = exp(gm_upper1))
```

## Confidence Intervals for the 90th Percentile

Because we use `eval()` and
call()`inside the`boot\_one()`function, we need to pass the function we want to bootstrap as a string. We can't pass in an anonymous function.  The function`call()`assembles a call object (unevaluated).  It's first argument must be a character string.  Then`eval()\`
evaluates the call, seeking the named function among function
identifiers in the current environment.

All of that could be addressed with more advanced R programming, such as
quoting function parameters or passing additional parameters to the call
object using R’s ellipsis operator (`...`). But for our current purpose,
it is far simpler to write a named function rather than revise and
generalize the `boot_one()` function. Besides, there are good packages
to support bootstrapping available. If we needed a more capable
bootstrap function, we would have used the `boot` package.

``` r
p90 <- function(.x) quantile(.x, 0.9)
```

This takes a while to run because calculating percentiles is harder than
calculating the mean.

``` r
p90_ci <- tapply(coli_data$ColiVal, coli_data$Station,
               function(d) boot_one(d, 'p90'))

# Convert to data frame (and then tibble...) 
p90_ci <- as.data.frame(do.call(rbind, p90_ci)) %>%
  rename(p90_lower1 = `2.5%`, p90_upper1 = `97.5%`) %>%
  rownames_to_column('Station')
```

## Combining Results

We add results to the summary data. (Because this uses `left_join()`,
rerunning it without deleting the old versions will generate errors in
later steps.)

``` r
sum_data <-  sum_data %>% 
  left_join(gm_ci, by = 'Station') %>%
  left_join(p90_ci, by = 'Station') %>%
  mutate(Station = fct_reorder(Station, gmean1))
```

## Geometric Mean Plot

``` r
plt <- ggplot(sum_data, aes(gmean1, Station)) + 
  geom_pointrange(aes(xmin = gm_lower1, xmax = gm_upper1),
                  color = cbep_colors()[6],
                  size = .2) +
  scale_x_log10(breaks = c(1,3,10,30, 100)) +
  
  xlab('Geometric Mean Fecal Coliforms\n(CFU / 100ml)') +

  ylab('Location') +
  
  geom_vline(xintercept = 14, lty = 2) +
  annotate('text', y = 30, x = 16, label = "14 CFU",
           size = 3, hjust = 0, angle = 270) +
  
  theme_cbep(base_size = 12) + 
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line  = element_line(size = 0.5, 
                                  color = 'gray85'),
        panel.grid.major.x = element_line(size = 0.5, 
                                          color = 'gray85', 
                                          linetype = 2)) 
plt
```

<img src="shellfish_bacteria_p90_graphics_files/figure-gfm/station_bootstrap_gm_graphics-1.png" style="display: block; margin: auto;" />

``` r
ggsave('figures/stations_gm_bootstrap.pdf', device = cairo_pdf, 
       width = 3.25, height = 5)
```

## p90 Plot

``` r
plt <- ggplot(sum_data, aes(p901, Station)) + 
  geom_pointrange(aes(xmin = p90_lower1, xmax = p90_upper1),
                  color = cbep_colors()[4],
                  size = .2) +
  scale_x_log10(breaks = c(1,3,10,30, 100)) +
  
  xlab('90th Percentile Fecal Coliforms\n(CFU / 100ml)') +

  ylab('Location') +
  
  geom_vline(xintercept = 31, lty = 2) +
  annotate('text', y = 30, x = 40, label = "31 CFU",
           size = 3, hjust = 0, angle = 270) +
  
  theme_cbep(base_size = 12) + 
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line  = element_line(size = 0.5, 
                                  color = 'gray85'),
        panel.grid.major.x = element_line(size = 0.5, 
                                          color = 'gray85', 
                                          linetype = 2)) 
plt
```

<img src="shellfish_bacteria_p90_graphics_files/figure-gfm/station_bootstrap_p90_graphics-1.png" style="display: block; margin: auto;" />

``` r
ggsave('figures/stations_p90_bootstrap.pdf', device = cairo_pdf, 
       width = 3.25, height = 5)
```

# Combined Graphics; Horizontal Layout

## Assemble Long Data

``` r
long_dat <- sum_data %>%
  select(Station, gmean1, gm_lower1, gm_upper1, p901, p90_lower1, p90_upper1) %>%
  rename(gm_value = gmean1, p90_value = p901) %>%
  rename_with( ~sub('1', '', .x)) %>%
  pivot_longer(gm_value:p90_upper,
               names_to = c('parameter', 'type'), 
               names_sep = '_') %>%
  pivot_wider(c(Station, parameter), names_from = type, values_from = value) %>%
  mutate(parameter = factor(parameter, 
                            levels = c('p90', 'gm'), 
                            labels = c('90th Percentile', 'Geometric Mean'))) %>%
  arrange(parameter, Station)
```

## Base Plot

``` r
plt <- ggplot(long_dat, aes(Station, value)) + 

  geom_linerange(aes(ymin = lower, ymax = upper, color = parameter),
                 alpha = 0.25,
                 size = .1) +
  geom_point(aes(color = parameter), size = 0.75) + 
               
  scale_y_log10(labels = scales::comma) +
  
  scale_color_manual(name = '', values = cbep_colors()[c(6,4)]) +

  ylab('Fecal Coliforms (CFU / 100ml)')+

  xlab('Sampling Sites Around Casco Bay') +
  expand_limits(x = 240) +  # this ensures the top dot is not cut off
  
  theme_cbep(base_size = 7) + 
  theme(legend.position = c(.2, 0.8),
        legend.text = element_text(size = 7),
        legend.key.height = unit(1, 'lines'),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 7),
        axis.line  = element_line(size = 0.5, 
                                  color = 'gray85'),
        panel.grid.major.y = element_line(size = 0.5, 
                                          color = 'gray85', 
                                          linetype = 2)) +
  guides(color = guide_legend(override.aes = list(lty = 0, size = 2)))
```

### Mimicing Graphic As Modified by the Graphic Designer

``` r
plt +
   geom_hline(yintercept = 32, 
             lty = 2, color = 'gray25') +
  geom_hline(yintercept = 14, 
             lty = 2, color = 'gray25') +
  
  annotate('text', x = 3, y = 40, label = '31 CFU',
           hjust = 0,
           size = 1.75, 
           color = cbep_colors()[6]) +
  annotate('text', x = 3, y = 18, label = '14 CFU',
           hjust = 0,
           size = 1.75, 
           color = cbep_colors()[4])
```

<img src="shellfish_bacteria_p90_graphics_files/figure-gfm/bootstrap__graphic_1-1.png" style="display: block; margin: auto;" />

``` r
ggsave('figures/stations_both_one_revised.pdf', device = cairo_pdf, 
       width = 4, height = 3)
```

### Alternate Annotations

``` r
plt +
   geom_hline(yintercept = 32, 
             lty = 2, color = 'gray25') +
  geom_hline(yintercept = 14, 
             lty = 2, color = 'gray25') +
  
  annotate('text', x = 3, y = 40, label = 'DMR P90 threshold, 31 CFU',
           hjust = 0,
           size = 1.75, 
           color = cbep_colors()[6]) +
  annotate('text', x = 3, y = 17.5, label = 'DMR Geometric Mean threshold, 14 CFU',
           hjust = 0,
           size = 1.75, 
           color = cbep_colors()[4])
```

<img src="shellfish_bacteria_p90_graphics_files/figure-gfm/bootstrap_graphic_2-1.png" style="display: block; margin: auto;" />

``` r
ggsave('figures/stations_both_one_revised_alt.pdf', device = cairo_pdf, 
       width = 4, height = 3)
```

## Shapes with outlines Plot

These plots look fairly poor here in the Markdown document, but the PDFs
work O.K. The goal was to make the markers standout a bit more

``` r
plt <- ggplot(long_dat, aes(Station, value)) + 

  geom_linerange(aes(ymin = lower, ymax = upper, color = parameter),
                 alpha = 0.25,
                 size = .1) +
  geom_point(aes(fill = parameter), color = 'gray60', 
             stroke = 0.25,  size = 1, shape = 21) + 
               
  scale_y_log10(labels = scales::comma) +
  
  scale_color_manual(name = '', values = cbep_colors()[c(6,4)]) +
  scale_fill_manual(name = '', values = cbep_colors()[c(6,4)]) +

  ylab('Fecal Coliforms (CFU / 100ml)')+

  xlab('Sampling Sites Around Casco Bay') +
  expand_limits(x = 240) +  # this ensures the top dot is not cut off
  
  theme_cbep(base_size = 7) + 
  theme(legend.position = c(.2, 0.8),
        legend.text = element_text(size = 7),
        legend.key.height = unit(1, 'lines'),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 7),
        axis.line  = element_line(size = 0.5, 
                                  color = 'gray85'),
        panel.grid.major.y = element_line(size = 0.5, 
                                          color = 'gray85', 
                                          linetype = 2)) +
  guides(color = guide_legend(override.aes = list(lty = 0, size = 2)))
```

``` r
plt +
  geom_hline(yintercept = 32, 
             lty = 2, color = 'gray25') +
  geom_hline(yintercept = 14, 
             lty = 2, color = 'gray25') +
  
  annotate('text', x = 3, y = 40, label = '31 CFU',
           hjust = 0,
           size = 1.75, 
           color = cbep_colors()[6]) +
  annotate('text', x = 3, y = 18, label = '14 CFU',
           hjust = 0,
           size = 1.75, 
           color = cbep_colors()[4])
```

<img src="shellfish_bacteria_p90_graphics_files/figure-gfm/bootstrap__graphic_3-1.png" style="display: block; margin: auto;" />

``` r
ggsave('figures/stations_both_one_revised_outlines.pdf', device = cairo_pdf, 
       width = 4, height = 3)
```

### Alternate Annotations

``` r
plt +
   geom_hline(yintercept = 32, 
             lty = 2, color = 'gray25') +
  geom_hline(yintercept = 14, 
             lty = 2, color = 'gray25') +
  
  annotate('text', x = 3, y = 40, label = 'DMR P90 threshold, 31 CFU',
           hjust = 0,
           size = 1.75, 
           color = cbep_colors()[6]) +
  annotate('text', x = 3, y = 17.5, label = 'DMR Geometric Mean threshold, 14 CFU',
           hjust = 0,
           size = 1.75, 
           color = cbep_colors()[4])
```

<img src="shellfish_bacteria_p90_graphics_files/figure-gfm/bootstrap_graphic_4-1.png" style="display: block; margin: auto;" />

``` r
ggsave('figures/stations_both_one_revised_alt_outlines.pdf', device = cairo_pdf, 
       width = 4, height = 3)
```

### Annotation Dataframe

We create a data frame that uses the same factor names, so we can direct
different annotations to each panel

``` r
annot_df <- tibble(parameter = factor(c(1,2), 
                                      labels = c('Geometric Mean',
                                                 '90th Percentile')),
                   threshold = c(14, 31),
                   txt = c('14 CFU', '31 CFU'),
                   adjust = c(3, 7.5),
                   adjust_no_bars = c(2, 4.5),
                   adjust_two = c(5, 12),
                   adjust_two_no_bars = c(3.5, 8))
```

# Combined Graphics: Vertical Layout

## One Panel

``` r
ggplot(long_dat, aes(value, Station)) + 
  geom_linerange(aes(xmin = lower, xmax = upper, color = parameter),
                 alpha = 0.5,
                 size = .1) +
  geom_point(aes(color = parameter), size = 0.75) + 
  scale_x_log10() +
  scale_color_manual(name = '', values = cbep_colors()[c(4,6)]) +
  scale_fill_manual(name = '', values = cbep_colors()[c(4,6)]) +
  xlab('Fecal Coliforms\n(CFU / 100ml)')+
  ylab('Location') +
  expand_limits(y = 240) +  # this ensures the top dot is not cut off
  geom_vline(data = annot_df,
             mapping = aes(xintercept = threshold), 
             lty = 2, color = 'gray25') +
  geom_text(data = annot_df, y = 30, 
             mapping = aes(x = threshold + adjust, 
                           label = txt),
             size = 3, hjust = 0, angle = 270) +
  theme_cbep(base_size = 12) + 
  theme(legend.position = c(.7, .25),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line  = element_line(size = 0.5, 
                                  color = 'gray85'),
        panel.grid.major.x = element_line(size = 0.5, 
                                          color = 'gray85', 
                                          linetype = 2)) +
  guides(color = guide_legend(override.aes = list(lty = 0, size = 2)))
```

<img src="shellfish_bacteria_p90_graphics_files/figure-gfm/station_bootstrap_one_graphics-1.png" style="display: block; margin: auto;" />

``` r
ggsave('figures/stations_both_one.pdf', device = cairo_pdf, 
       width = 4, height = 5)
```

## Two Panels

``` r
ggplot(long_dat, aes(value, Station)) + 
  geom_linerange(aes(xmin = lower, xmax = upper, color = parameter),
                 alpha = 0.5,
                 size = .1) +
  geom_point(aes(color = parameter), size = 0.75) + 
  scale_x_log10() +
  facet_wrap('parameter', nrow = 1) +
  scale_color_manual(values = cbep_colors()[c(4,6)]) +
  scale_fill_manual(values = cbep_colors()[c(4,6)]) +
  xlab('Fecal Coliforms\n(CFU / 100ml)') +
  ylab('Location') +
  geom_vline(data = annot_df,
             mapping = aes(xintercept = threshold), 
             lty = 2, color = 'gray25') +
  geom_text(data = annot_df, y = 30, 
             mapping = aes(x = threshold + adjust_two, label = txt),
             size = 3, hjust = 0, angle = 270) +
  theme_cbep(base_size = 12) + 
  theme(legend.position = 'None',
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line  = element_line(size = 0.5, 
                                  color = 'gray85'),
        panel.grid.major.x = element_line(size = 0.5, 
                                          color = 'gray85', 
                                          linetype = 2)) 
```

<img src="shellfish_bacteria_p90_graphics_files/figure-gfm/station_bootstrap_combined_graphics-1.png" style="display: block; margin: auto;" />

``` r
ggsave('figures/stations_both_two.pdf', device = cairo_pdf, 
       width = 5, height = 5)
```

# Years Plot

Technically, this produces a plot with summary statistics calculated on
the log-transformed data. That means `stat_sumamry()` correctly displays
the geometric mean. However the `quantile()` function finds (and plots)
the 90th percentile of the LOG of the data, rather than the log
transform of the 90th percentile of the untransformed data. Given the
abundance of data here, that makes no practical difference in a graphic.
Nevertheless, we commented out the lines in this code that would add the
90th percentiles to the plot.

``` r
plt <- ggplot(coli_data, aes(YEAR, ColiVal)) +
  geom_jitter(aes(color = LCFlag), alpha = 0.25, height = 0.05, width = 0.4) +
  ## We use the MEAN here because `stat_summary()` works on data after
  ## applying the transformation to the y axis, thus implicitly calculating the
  ## geometric mean.
  stat_summary(fun = mean, 
               fill = 'red', shape = 22) +
 # stat_summary(fun = ~ quantile(.x, 0.9),
 #              fill = 'orange', shape = 23) +
  
  scale_color_manual(values = cbep_colors(), name = '',
                     labels = c('Observed', 'Below Detection')) +
  
  
  xlab('') +
  ylab('Fecal Coliforms\n(CFU / 100ml)') +

  scale_y_log10() +

  theme_cbep(base_size = 12) +
  theme(legend.position = "bottom") +
  
  guides(color = guide_legend(override.aes = list(alpha = c(0.5,0.75),
                                                  size = 3) ))
```

``` r
 plt +  
  geom_hline(yintercept = 14, lty = 2) +
  annotate('text', x = 2020, y = 17, 
           size = 3, hjust = .75, label = "14 CFU") +
  
  geom_hline(yintercept = 31, lty = 2) +
  annotate('text', x = 2020, y = 40, 
           size = 3, hjust = .75, label = "31 CFU") +
  scale_x_continuous(breaks = c(2015, 2017, 2019))
#> Warning: Removed 5 rows containing missing values (geom_segment).
```

<img src="shellfish_bacteria_p90_graphics_files/figure-gfm/years_add_ref_lines-1.png" style="display: block; margin: auto;" />

``` r
ggsave('figures/years.pdf', device = cairo_pdf, 
       width = 5, height = 4)
#> Warning: Removed 5 rows containing missing values (geom_segment).
```

# Frequency of Exceedences

``` r
sum_data %>%
  summarize(exceeds_p90 = sum(p901 >= 31),
            exceeds_gm = sum(gmean1 > 14))
#> # A tibble: 1 x 2
#>   exceeds_p90 exceeds_gm
#>         <int>      <int>
#> 1          47          1
```

So just about 20% of sites fail the p90 “Accepted” threshold.
