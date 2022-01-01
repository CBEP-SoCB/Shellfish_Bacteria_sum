Frequency of Exceedences, Shellfish Bacteria Data
================
Curtis C. Bohlen, Casco Bay Estuary Partnership.
02/17/2021

-   [Introduction](#introduction)
    -   [Growing Area Classification
        Standards](#growing-area-classification-standards)
        -   [Maine State Class SB Waters
            Standards](#maine-state-class-sb-waters-standards)
-   [Load Libraries](#load-libraries)
-   [Load Data](#load-data)
    -   [Main Data](#main-data)
        -   [Remove NAs](#remove-nas)
    -   [Weather Data](#weather-data)
    -   [Incorporate Weather Data](#incorporate-weather-data)
    -   [Remove Sites not in Region](#remove-sites-not-in-region)
    -   [Calculate Indicator Variables](#calculate-indicator-variables)
    -   [Raw Station Probabilitiies](#raw-station-probabilitiies)
    -   [Add Imperviousness Data](#add-imperviousness-data)
-   [Export Data for GIS](#export-data-for-gis)
-   [Exploratory Graphics](#exploratory-graphics)
-   [Utility Functions](#utility-functions)
-   [Cleanup](#cleanup)
-   [Modeling](#modeling)
    -   [Binomial Models](#binomial-models)
        -   [Station Only](#station-only)
        -   [Region Model](#region-model)
        -   [Season Model](#season-model)
        -   [Rainfall Model](#rainfall-model)
        -   [Imperviousness Models](#imperviousness-models)
    -   [Proportional Odds Models](#proportional-odds-models)
        -   [Summary of the Data](#summary-of-the-data)
        -   [Region Models](#region-models)
        -   [Using MASS](#using-mass)
        -   [Rainfall Models](#rainfall-models)
-   [Further Questions](#further-questions)

<img
    src="https://www.cascobayestuary.org/wp-content/uploads/2014/04/logo_sm.jpg"
    style="position:absolute;top:10px;right:50px;" />

# Introduction

The bacteria data is highly skewed. We have found no ready way to model
a geometric mean for these data in a robust way. The data appears
distributed close to a Pareto distribution, which is highly skewed, with
a heavy right tail. This distribution is likely to be difficult to
model, so we can not readily estimate parameters, moments and other
summary statistics.

One possibility would be to run some sort of (generalized) linear model
on the log of data, calculate means and medians, and back transform. We
explore linear model approaches in another R Notebook.

Here we take another approach, which is to focus on binomial,
quasi-binomial, and multinomial proportional odds models to estimate
probabilities of exceeding different regulatory thresholds.

## Growing Area Classification Standards

| Growing Area Classification | Activity Allowed                                                          | Geometric mean FC/100ml | 90th Percentile (P90) FC/100ml |
|-----------------------------|---------------------------------------------------------------------------|-------------------------|--------------------------------|
| Approved                    | Harvesting allowed                                                        | ≤ 14                    | ≤ 31                           |
| Conditionally Approved      | Harvesting allowed except during specified conditions                     | ≤ 14 in open status     | ≤ 31 in open status            |
| Restricted                  | Depuration harvesting or relay only                                       | ≤ 88 and &gt;15         | ≤ 163 and &gt; 31              |
| Conditionally Restricted    | Depuration harvesting or relay allowed except during specified conditions | ≤ 88 in open status     | ≤ 163 in open status           |
| Prohibited                  | Aquaculture seed production only                                          | &gt;88                  | &gt;163                        |

So, critical levels for Geometric Mean include:

-   *G**M* ≤ 14: Approved, or Open status at Conditionally Approved
    sites

-   *G**M* ≤ 88: Depuration harvesting or Relay Only

-   *G**M* &gt; 88 : Prohibited

And for the p90:

-   *P*90 &lt; 31 Approved or Open status at Conditionally Approved

-   *P*90 ≤ 163 Depuration harvesting or Relay Only

-   *P*90 &gt; 163 Prohibited

### Maine State Class SB Waters Standards

Maine’s water quality criteria includes an additional standard, which
applies only indirectly to these data:  
&gt; the number of enterococcus bacteria in these waters may not exceed
a geometric mean of 8 CFU per 100 milliliters in any 90-day interval or
54 CFU per 100 milliliters in more than 10% of the samples in any 90-day
interval.

38 M.R.S. §465-B(2)(B)

A “90 day interval” might apply to a summer’s worth of data, but in most
years that will only represent a handful of observations at each site,
so it is of limited value. More seriously, the standard is written in
terms of “enterococcus” bacteria, not the fecal coliform data used by
DMR.

# Load Libraries

``` r
library(MASS)   # Load before tidyverse because it has a select() function
library(mgcv)   # For GAMs and GAMMs; used her for seasonal smoothers
#> Warning: package 'mgcv' was built under R version 4.0.5
#> Loading required package: nlme
#> This is mgcv 1.8-38. For overview type 'help("mgcv-package")'.
library(tidyverse)
#> Warning: package 'tidyverse' was built under R version 4.0.5
#> -- Attaching packages --------------------------------------- tidyverse 1.3.1 --
#> v ggplot2 3.3.5     v purrr   0.3.4
#> v tibble  3.1.6     v dplyr   1.0.7
#> v tidyr   1.1.4     v stringr 1.4.0
#> v readr   2.1.0     v forcats 0.5.1
#> Warning: package 'ggplot2' was built under R version 4.0.5
#> Warning: package 'tidyr' was built under R version 4.0.5
#> Warning: package 'dplyr' was built under R version 4.0.5
#> Warning: package 'forcats' was built under R version 4.0.5
#> -- Conflicts ------------------------------------------ tidyverse_conflicts() --
#> x dplyr::collapse() masks nlme::collapse()
#> x dplyr::filter()   masks stats::filter()
#> x dplyr::lag()      masks stats::lag()
#> x dplyr::select()   masks MASS::select()

library(readr)
library(GGally)
#> Warning: package 'GGally' was built under R version 4.0.5
#> Registered S3 method overwritten by 'GGally':
#>   method from   
#>   +.gg   ggplot2

library(emmeans)   # For marginal means
#> Warning: package 'emmeans' was built under R version 4.0.5
#> 
#> Attaching package: 'emmeans'
#> The following object is masked from 'package:GGally':
#> 
#>     pigs
library(mblm)      # for the Thiel-Sen estimators -- not really successful here
library(VGAM)      # For Pareto GLMs and estimation.
#> Loading required package: stats4
#> Loading required package: splines
#> 
#> Attaching package: 'VGAM'
#> The following object is masked from 'package:tidyr':
#> 
#>     fill
#> The following object is masked from 'package:mgcv':
#> 
#>     s
library(mgcv)      # For GAMs, here used principally for hierarchical models

library(CBEPgraphics)
load_cbep_fonts()
theme_set(theme_cbep())
```

# Load Data

## Main Data

``` r
sibfldnm <- 'Derived_Data'
parent <- dirname(getwd())
sibling <- file.path(parent,sibfldnm)
fl1<- "Shellfish data 2015 2018.csv"
path <- file.path(sibling, fl1)

coli_data <- read_csv(path, 
  col_types = cols(SDate = col_date(format = "%Y-%m-%d"), 
        SDateTime = col_datetime(format = "%Y-%m-%dT%H:%M:%SZ"), # Note Format!
        STime = col_time(format = "%H:%M:%S"))) %>%
  mutate(across(Station:WDIR, factor)) %>%
  mutate(Class = factor(Class, levels = c( 'A', 'CA', 'CR',
                                           'R', 'P', 'X' ))) %>%
  mutate(Tide = factor(Tide, levels = c("L", "LF", "F", "HF",
                                        "H", "HE", "E", "LE"))) %>%
  mutate(DOY = as.numeric(format(SDate, format = '%j')),
         Month = as.numeric(format(SDate, format = '%m'))) %>%
  mutate(Month = factor(Month, levels = 1:12, labels = month.abb)) %>%
  rename_with(tolower)
```

### Remove NAs

``` r
coli_data <- coli_data %>%
  filter (! is.na(colival))
```

## Weather Data

We simplify the weather data somewhat.

``` r
sibfldnm    <- 'Original_Data'
parent      <- dirname(getwd())
sibling     <- file.path(parent,sibfldnm)

fn <- "Portland_Jetport_2015-2019.csv"
fpath <- file.path(sibling, fn)

weather_data <- read_csv(fpath, 
 col_types = cols(AWNDattr = col_skip(), 
        FMTM = col_skip(), FMTMattr = col_skip(), 
        PGTM = col_skip(), PGTMattr = col_skip(),
        PRCPattr = col_character(), SNOWattr = col_character(), 
        SNWD = col_skip(), SNWDattr = col_skip(),
        TAVG = col_number(), TAVGattr = col_character(), 
        TMIN = col_number(), TMINattr = col_character(), 
        TMAX = col_number(), TMAXattr = col_character(), 
        station = col_skip())) %>%
  select( ! starts_with('W')) %>%
  select(! ends_with('attr')) %>%
  rename(sdate = date,
         Precip=PRCP,
         MaxT = TMAX,
         MinT= TMIN,
         AvgT = TAVG,
         Snow = SNOW) %>%
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
         Log1Precip    = log1p(Precip), 
         Log1Precip_d1 = log1p(Precip_d1),
         Log1Precip_d2 = log1p(Precip_d2),
         Log1Precip_2   = log1p(Precip_d1 + Precip_d2),
         Log1Precip_3   = log1p(Precip + Precip_d1 + Precip_d2)) %>%
  rename_with(tolower)
```

## Incorporate Weather Data

``` r
coli_data <- coli_data %>%
  left_join(weather_data, by = 'sdate')
```

## Remove Sites not in Region

We have some data that was selected for stations outside of Casco Bay.
To be  
careful, we remove sampling data for any site in the two adjacent
Growing Areas, “WH” and “WM”.

``` r
coli_data <- coli_data %>%
  filter(grow_area != 'WH' & grow_area != "WM") %>%
  mutate(grow_area = fct_drop(grow_area))
```

## Calculate Indicator Variables

Here we calculate variables that indicate whether each sample exceeds
our four thresholds. Then they are combined to produce three multinomial
ordered factors.

The key thresholds are these:

``` r
coli_limits <- list(open = 0,   gmrelay=14,    p90relay=31, 
                    gmclosed=88, p90closed=163, high= 50000)
```

We create a data frame containing a number of binomial and multinomial
responses.

-   The first four are binary responses, showing whether each individual
    sample meets or fails to meet a single standard. Since the standards
    are written to apply to long-term records, that is not 100%
    legitimate, especially for the geometric mean standards, but it is
    what we have to work with. We can apply the p90 standard by ensuring
    that the probability of exceeding that threshold is less than 0.90.

-   The next three variables are multinomial responses that classify
    variables into a sequence of of ordered categories

``` r
freq_data <- coli_data %>%
  mutate(gm_open   = colival <= coli_limits$gmrelay,
         gm_relay  = colival <= coli_limits$gmclosed,
         p90_open  = colival <= coli_limits$p90relay,
         p90_relay = colival <= coli_limits$p90closed) %>%
  mutate(all_lvls = cut(colival, coli_limits,
                    labels = c('open', 'gm_relay', 'p90_relay', 
                               'gm_closed', 'p90_closed'),
                    ordered_result = TRUE)) %>%
  mutate(p90_lvls = cut(colival, coli_limits[c(1,3,5,6)], 
                    labels = c('p90open', 'p90relay', 'p90closed'),
                    ordered_result = TRUE)) %>%
  mutate(gm_lvls = cut(colival, coli_limits[c(1,2,4,6)], 
                      labels = c('gmopen', 'gmrelay', 'gmclosed'),
                      ordered_result = TRUE))
```

``` r
freq_data <- freq_data %>% 
  select(-coliscore, -rawcoli, -lcflag, -rcflag, -colival)
```

## Raw Station Probabilitiies

We want to be able to compare the results of modeling to observed
relative frequencies of meeting or exceeding standards, so we want a
simple data frame containing those probabilities. We structured our
binomial observations with TRUE = open, so these are probabilities of
meeting standards.

``` r
rawprobs <- freq_data %>%
  group_by(station) %>%
  summarize(grow_area = first(grow_area),
            p_gm_open  = sum(gm_open)/sum(! is.na(gm_open)),
            p_p90_open = sum(p90_open)/sum(! is.na(p90_open)),
            p_gm_relay  = sum(gm_relay)/sum(! is.na(gm_relay)),
            p_p90_relay = sum(p90_relay)/sum(! is.na(p90_relay)))

lowerfun <- function(data, mapping){
  ggplot(data = data, mapping = mapping)+
    geom_point() +
    geom_abline(slope = 1, intercept = 0) +
    scale_x_continuous(limits = c(0,1)) +
    scale_y_continuous(limits = c(0,1))
  }  

rawprobs %>%
  select(-station, - grow_area) %>%
  ggpairs(lower = list(continuous = wrap(lowerfun)),
          progress = FALSE) +
  theme_cbep(base_size = 10)
```

<img src="frequency-analysis_files/figure-gfm/calculate_observed_frequencies-1.png" style="display: block; margin: auto;" />

Which looks pretty good. Recall that the reason the points are all above
the 1:1 line is that it is impossible to exceed the higher standard
without also exceeding the lower threshold. Note the relatively high
correlations.

That one extreme value probably corresponds to the one very high
geometric mean site we already observed.

## Add Imperviousness Data

``` r
sibfldnm    <- 'Derived_Data'
parent      <- dirname(getwd())
sibling     <- file.path(parent,sibfldnm)

fn <- "station_imperviousness.csv"
fpath <- file.path(sibling, fn)

imperv_data <- read_csv(fpath) %>%
  select(Station, pct_1000, pct_500, pct_100, pct_l_1000, pct_l_500, pct_l_100)
#> Rows: 259 Columns: 17
#> -- Column specification --------------------------------------------------------
#> Delimiter: ","
#> chr  (2): Station, Grow_Area
#> dbl (15): OBJECTID, Lat, Long, land_1000, land_500, land_100, imperv_1000, i...
#> 
#> i Use `spec()` to retrieve the full column specification for this data.
#> i Specify the column types or set `show_col_types = FALSE` to quiet this message.
```

``` r
freq_data <- freq_data %>%
  mutate(s = as.character(station)) %>%
  left_join(imperv_data, by = c('s' = 'Station')) %>%
  select (-s)

rawprobs<- rawprobs %>%
  mutate(s = as.character(station)) %>%
  left_join(imperv_data, by = c('s' = 'Station')) %>%
  select (-s)

rm(imperv_data)
```

# Export Data for GIS

``` r
freq_data %>%
  write_csv('shellfish_exceeds_data.csv')
```

# Exploratory Graphics

``` r
freq_data %>% select(month, all_lvls) %>%
  ggpairs(aes(month, all_lvls)) 
```

<img src="frequency-analysis_files/figure-gfm/simple_pairsplot-1.png" style="display: block; margin: auto;" />

1.  Most observations are from summer months  
2.  Most observations are lower than all our cut points  
3.  Probability of exceedences appears slightly higher in summer months.

``` r
freq_data %>%
  count(month, all_lvls) %>%
  ggplot(aes(month, n, fill=all_lvls)) +
   geom_bar(stat = "identity", position='fill') +
  theme_minimal()
```

<img src="frequency-analysis_files/figure-gfm/freq_bar_all-1.png" style="display: block; margin: auto;" />

The following is the equivalent of 1- p(anything better than P90\_relay)
in the prior graphic.

``` r
freq_data %>%
  count(month, p90_open) %>%
  ggplot(aes(month, n, fill= p90_open )) +
   geom_bar(stat = "identity", position='fill') +
  theme_minimal()
```

<img src="frequency-analysis_files/figure-gfm/freq_bar_simple-1.png" style="display: block; margin: auto;" />

# Utility Functions

The coefficients of a binomial or multinomial GLM are actually the logit
(log odds) of the probabilities we are interested in, so a couple of
utility functions may come in handy. These are probably not the most
numerically stable versions of these functions for extreme values, but
they work.

``` r
logit <- function(p) {log(p)-log(1-p)}
inv_logit <- function(x) {1/(1+exp(-x))}
```

# Cleanup

We remove unnecessary data frames to free memory for calculations.

``` r
rm(coli_data, coli_limits, weather_data)
```

# Modeling

We can explore three different modeling strategies: 1. A binomial GLM
for exceedences of the gm\_mean and p90 thresholds. 2. A multinomial
proportional odds model using the `polr()` function from `MASS` or
`vglm()` from `VGAM`.

Ideally, we should account for different sampling histories to provide
robust comparable estimates of site probabilities.

## Binomial Models

We focus on the probability of meeting the lower p90 threshold
(`p90_open`). That threshold has potential consequences. If the
probability of violating that standard is high enough (&gt; 90%) the
site would have to be approved only for relay. Stations failed the
higher `p90_relay` standard so rarely that we run into problems with
modeling.

It is less clear how we would relate the probability of violating the
Geometric Mean standards to our models.

### Station Only

We start with a simple model looking only at stations. The probability
of failing a standard may be an appropriate way to symbolize stations in
GIS.

``` r
system.time(p90_open_glm_1 <- glm(p90_open  ~ station,
             data = freq_data,
             family=binomial(link=logit)))
#>    user  system elapsed 
#>    8.70    0.09    8.80
```

``` r
anova(p90_open_glm_1, test='LRT')
#> Analysis of Deviance Table
#> 
#> Model: binomial, link: logit
#> 
#> Response: p90_open
#> 
#> Terms added sequentially (first to last)
#> 
#> 
#>          Df Deviance Resid. Df Resid. Dev  Pr(>Chi)    
#> NULL                      9398     4632.6              
#> station 237   661.34      9161     3971.3 < 2.2e-16 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

By the LRT, differences among station are highly significant.

``` r
p_res <- summary(emmeans(p90_open_glm_1, "station", 
                        type = 'response')) %>%
  mutate(station = fct_reorder(station, prob)) %>%
  arrange(station)


plot(p_res) +
  
  theme_cbep(base_size = 12) + 
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line  = element_line(size = 0.5, 
                                  color = 'gray85'),
        panel.grid.major.x = element_line(size = 0.5, 
                                          color = 'gray85', 
                                          linetype = 2)) +
  
  ylab('Station') +
  xlab('Probability of Meeting\nP90 Standard' )
```

<img src="frequency-analysis_files/figure-gfm/station_glm_graphic-1.png" style="display: block; margin: auto;" />

Many sites show unstable standard errors, as they were never observed
with a sample that failed the standard, triggering a Hauke-Donner
effect. That makes estimation of parameters and standard errors
effectively impossible for those sites. We can treat them as P = 1.0 for
practical purposes in a model like this.

``` r
tmp <- rawprobs %>%
  left_join(p_res, by = 'station')

ggplot(tmp, aes(p_p90_open, prob)) +
    geom_pointrange(mapping = aes(ymin = asymp.LCL, ymax =asymp.UCL)) +
    geom_abline(intercept = 0, slope = 1) +
    scale_x_continuous(limits = c(.5,1)) +
    scale_y_continuous(limits = c(.5,1))
#> Warning: Removed 66 rows containing missing values (geom_segment).
```

<img src="frequency-analysis_files/figure-gfm/observed_vs_predicted-1.png" style="display: block; margin: auto;" />

In a simple model, the predictions match the observed relative
frequencies and the standard errors behave as expected (wider near the
middle of a binomial distribution, narrower and asymmetrical near the
limits).

### Region Model

For more complex models, we treat Stations as random factors, using the
`gam()` function, which can fit random factors using a smoothing term
with basis designated as `bs = 're'`. This helps protect against making
claims on the basis of the specific stations from which data are
available.

``` r
system.time(
  p90_open_gam_regions <- gam(p90_open  ~ grow_area + s(station, bs = 're'),
             data = freq_data,
             family=binomial(link=logit))
  )
#>    user  system elapsed 
#>   51.22    0.62   51.85
```

``` r
anova(p90_open_gam_regions, test='LRT')
#> Warning in anova.gam(p90_open_gam_regions, test = "LRT"): test argument ignored
#> 
#> Family: binomial 
#> Link function: logit 
#> 
#> Formula:
#> p90_open ~ grow_area + s(station, bs = "re")
#> 
#> Parametric Terms:
#>           df Chi.sq  p-value
#> grow_area  3  28.51 2.84e-06
#> 
#> Approximate significance of smooth terms:
#>            edf Ref.df Chi.sq p-value
#> s(station) 127    234  304.5  <2e-16
```

By the LRT, differences among regions are highly significant, as
expected.

``` r
p_res <- summary(emmeans(p90_open_gam_regions, "grow_area", 
                         nesting = " station %in% grow_area",
                        type = 'response'))

plot(p_res) +
  xlab('Region') +
  ylab('Probability of Meeting\nP90 Standard' ) +
  coord_flip()
```

<img src="frequency-analysis_files/figure-gfm/region_draft_graphic-1.png" style="display: block; margin: auto;" />

#### Draft Graphic

``` r
plt <- ggplot(p_res, aes(grow_area, prob)) + 
  geom_jitter(data = rawprobs,
              mapping = aes(y = p_p90_open),
              height = 0, width = 0.4, color = cbep_colors()[5] ) +
  geom_pointrange(aes(ymin = lower.CL, ymax = upper.CL),
                  fill = cbep_colors()[2], size = .75, shape = 23) +
  geom_hline(yintercept = 0.9, color = 'gray85', lty =2, size = 1) +
  
  ylab(expression(atop('Frequency of Meeting', 
                           'P90 ' ~ 'Fecal Coliforms' ~ ' Standard'))) +
  xlab('Maine DMR Growing Area') +

  theme_cbep(base_size = 12)
plt
```

<img src="frequency-analysis_files/figure-gfm/plot_regions_emm-1.png" style="display: block; margin: auto;" />

### Season Model

``` r
system.time(p90_open_gam_months <- gam(p90_open ~ month + s(station, bs = 're'),
             data = freq_data,
             family=binomial(link=logit)))
#>    user  system elapsed 
#>   57.74    0.51   58.25
```

``` r
anova(p90_open_gam_months)
#> 
#> Family: binomial 
#> Link function: logit 
#> 
#> Formula:
#> p90_open ~ month + s(station, bs = "re")
#> 
#> Parametric Terms:
#>       df Chi.sq p-value
#> month 11  201.4  <2e-16
#> 
#> Approximate significance of smooth terms:
#>              edf Ref.df Chi.sq p-value
#> s(station) 158.9  237.0  448.1  <2e-16
```

We see strong evidence that time of year matters, as expected from our
analysis of geometric means.

``` r
mm <- emmeans(p90_open_gam_months, "month", type = 'response')
mms <- summary(mm)  # Summary has class dataframe, making access and display easier.
```

#### Draft Graphic

We have a problem here, as Stations are not contained “within” months,
so sites are not uniquely attributable to them. Accordingly, we can not
use the strategy of the prior graphic to show both Stations and Months.

``` r
plt <- ggplot(mms, aes(month, prob)) + 
  geom_pointrange(aes(ymin = lower.CL, ymax = upper.CL),
                  fill = cbep_colors()[2], size = .75, shape = 23) +
  geom_line(aes(x = as.numeric(month)), color =  cbep_colors()[3]) +
  geom_hline(yintercept = 0.9, color = 'gray85', lty =2, size = 1) +
  
  xlab('') + 
  ylab(expression(atop('Probability of Meeting', 
                           'P90 ' ~ 'Fecal Coliforms' ~ ' Standard'))) +
  ylim(0.75, 1.0) +
  theme_cbep(base_size = 12)
plt
```

<img src="frequency-analysis_files/figure-gfm/plot_months_emm-1.png" style="display: block; margin: auto;" />

### Rainfall Model

Again, we artificially reduce the number of degrees of freedom
associated with each smoother to make them less wiggly.

``` r
system.time(
  p90_open_gam_rain <- gam(p90_open ~ s(log1precip, k = 5) + 
                             s(log1precip_d1, k = 5) + 
                             s(station, bs = 're'),
                           data = freq_data,
                           family=binomial(link=logit)))
#>    user  system elapsed 
#>   90.27    1.01   91.32
```

``` r
gam.check(p90_open_gam_rain)
```

<img src="frequency-analysis_files/figure-gfm/check_rainfall_gam-1.png" style="display: block; margin: auto;" />

    #> 
    #> Method: UBRE   Optimizer: outer newton
    #> full convergence after 6 iterations.
    #> Gradient range [1.169244e-08,9.880935e-07]
    #> (score -0.5832307 & scale 1).
    #> Hessian positive definite, eigenvalue range [1.813532e-06,0.002113119].
    #> Model rank =  247 / 247 
    #> 
    #> Basis dimension (k) checking results. Low p-value (k-index<1) may
    #> indicate that k is too low, especially if edf is close to k'.
    #> 
    #>                      k'    edf k-index p-value  
    #> s(log1precip)      4.00   3.97    0.92    0.12  
    #> s(log1precip_d1)   4.00   3.98    0.91    0.04 *
    #> s(station)       238.00 153.80      NA      NA  
    #> ---
    #> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
anova(p90_open_gam_rain)
#> 
#> Family: binomial 
#> Link function: logit 
#> 
#> Formula:
#> p90_open ~ s(log1precip, k = 5) + s(log1precip_d1, k = 5) + s(station, 
#>     bs = "re")
#> 
#> Approximate significance of smooth terms:
#>                      edf  Ref.df Chi.sq p-value
#> s(log1precip)      3.966   3.999  203.4  <2e-16
#> s(log1precip_d1)   3.977   3.999  299.0  <2e-16
#> s(station)       153.796 237.000  436.7  <2e-16
```

We see strong evidence that rainfall matters, as expected from our
analysis of geometric means. Both current day rainfall and prior day
rainfall matter. Relationship of the linear predictors with rainfall are
close to linear, except at the extremes.

``` r
summary(p90_open_gam_rain)
#> 
#> Family: binomial 
#> Link function: logit 
#> 
#> Formula:
#> p90_open ~ s(log1precip, k = 5) + s(log1precip_d1, k = 5) + s(station, 
#>     bs = "re")
#> 
#> Parametric coefficients:
#>             Estimate Std. Error z value Pr(>|z|)    
#> (Intercept)  3.23710    0.09335   34.68   <2e-16 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> Approximate significance of smooth terms:
#>                      edf  Ref.df Chi.sq p-value    
#> s(log1precip)      3.966   3.999  203.4  <2e-16 ***
#> s(log1precip_d1)   3.977   3.999  299.0  <2e-16 ***
#> s(station)       153.796 237.000  436.7  <2e-16 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> R-sq.(adj) =  0.161   Deviance explained = 22.5%
#> UBRE = -0.58323  Scale est. = 1         n = 9397
```

``` r
plot(p90_open_gam_rain)
```

<img src="frequency-analysis_files/figure-gfm/plot_rainfall_gam-1.png" style="display: block; margin: auto;" /><img src="frequency-analysis_files/figure-gfm/plot_rainfall_gam-2.png" style="display: block; margin: auto;" /><img src="frequency-analysis_files/figure-gfm/plot_rainfall_gam-3.png" style="display: block; margin: auto;" />

### Imperviousness Models

Is risk of failing standards related to local land use, as indexed by
imperviousness? We examine value of percent imperviousness as a
predictor of risk.

#### Exploratory Graphic

We display the horizontal axis using a square root transform, which is
may stabilize the variation of a percentage.

``` r
rawprobs %>%
  select(p_p90_open, starts_with('pct_')) %>%
  pivot_longer(starts_with('pct_'), 
               names_to = 'category', 
               values_to = 'value') %>%
  ggplot(aes(value, p_p90_open)) +
  geom_point(alpha = 0.25) +
  geom_smooth(method = 'lm') +
  #scale_x_sqrt() +
  facet_wrap(~category, nrow = 2, scales = 'free_x')
#> `geom_smooth()` using formula 'y ~ x'
#> Warning: Removed 12 rows containing non-finite values (stat_smooth).
#> Warning: Removed 12 rows containing missing values (geom_point).
```

<img src="frequency-analysis_files/figure-gfm/imperv_expl_grap-1.png" style="display: block; margin: auto;" />

While those almost all suggest a decline in probability of meeting
standards with imperviousness, the relationships rest on the large
number of sites that never fail the standard at low impervious cover.

Intellectually, I prefer the land based IC values. There is no great way
to decide which range makes the most sense. A full km (1000 m) feels
like a long way on the narrow peninsulas of the Eastern Bay.

#### Models

Here, we need to model each station, since risk is a station by station
property. We could treat Stations as either random or fixed factors, but
it is more intellectually honest to acknowledge that we want to draw
lessons about a larger population of possible sampling locations. We
treat Stations as random factors.

We run all three distances, using land-based percent cover estimates.
First, with untransformed percent cover values.

``` r
system.time(imp_100 <- gam(p90_open ~ pct_l_100 + s(station, bs = 're'),
             data = freq_data,
             family=binomial(link=logit)))
#>    user  system elapsed 
#>   41.92    0.50   42.42
system.time(imp_500 <- gam(p90_open ~ pct_l_500 + s(station, bs = 're'),
             data = freq_data,
             family=binomial(link=logit)))
#>    user  system elapsed 
#>   44.48    0.66   45.17
system.time(imp_1000 <- gam(p90_open ~ pct_l_1000 + s(station, bs = 're'),
             data = freq_data,
             family=binomial(link=logit)))
#>    user  system elapsed 
#>   48.17    0.63   48.81
```

``` r
anova(imp_100)
#> 
#> Family: binomial 
#> Link function: logit 
#> 
#> Formula:
#> p90_open ~ pct_l_100 + s(station, bs = "re")
#> 
#> Parametric Terms:
#>           df Chi.sq p-value
#> pct_l_100  1  0.634   0.426
#> 
#> Approximate significance of smooth terms:
#>              edf Ref.df Chi.sq p-value
#> s(station) 135.5  224.0  372.6  <2e-16
cat('\n\n')
anova(imp_500)
#> 
#> Family: binomial 
#> Link function: logit 
#> 
#> Formula:
#> p90_open ~ pct_l_500 + s(station, bs = "re")
#> 
#> Parametric Terms:
#>           df Chi.sq p-value
#> pct_l_500  1  0.605   0.437
#> 
#> Approximate significance of smooth terms:
#>              edf Ref.df Chi.sq p-value
#> s(station) 146.1  236.0    409  <2e-16
cat('\n\n')
anova(imp_1000)
#> 
#> Family: binomial 
#> Link function: logit 
#> 
#> Formula:
#> p90_open ~ pct_l_1000 + s(station, bs = "re")
#> 
#> Parametric Terms:
#>            df Chi.sq p-value
#> pct_l_1000  1   6.01  0.0142
#> 
#> Approximate significance of smooth terms:
#>              edf Ref.df Chi.sq p-value
#> s(station) 141.9  236.0  390.7  <2e-16
```

The only version that shows a statistically robust response is the 1000
meter model. The others are not significant.

Now we run similar models on the square root of percent cover.

``` r
system.time(imp_100_t <- gam(p90_open ~ sqrt(pct_l_100) + s(station, bs = 're'),
             data = freq_data,
             family=binomial(link=logit)))
#>    user  system elapsed 
#>   35.81    0.53   36.34
system.time(imp_500_t <- gam(p90_open ~ sqrt(pct_l_500) + s(station, bs = 're'),
             data = freq_data,
             family=binomial(link=logit)))
#>    user  system elapsed 
#>   43.91    0.53   44.45
system.time(imp_1000_t <- gam(p90_open ~ sqrt(pct_l_1000) + s(station, bs = 're'),
             data = freq_data,
             family=binomial(link=logit)))
#>    user  system elapsed 
#>   53.55    0.57   54.13
```

``` r
anova(imp_100_t)
#> 
#> Family: binomial 
#> Link function: logit 
#> 
#> Formula:
#> p90_open ~ sqrt(pct_l_100) + s(station, bs = "re")
#> 
#> Parametric Terms:
#>                 df Chi.sq p-value
#> sqrt(pct_l_100)  1  0.131   0.718
#> 
#> Approximate significance of smooth terms:
#>              edf Ref.df Chi.sq p-value
#> s(station) 135.7  224.0  373.5  <2e-16
cat('\n\n')
anova(imp_500_t)
#> 
#> Family: binomial 
#> Link function: logit 
#> 
#> Formula:
#> p90_open ~ sqrt(pct_l_500) + s(station, bs = "re")
#> 
#> Parametric Terms:
#>                 df Chi.sq p-value
#> sqrt(pct_l_500)  1  1.345   0.246
#> 
#> Approximate significance of smooth terms:
#>              edf Ref.df Chi.sq p-value
#> s(station) 145.7  236.0  408.3  <2e-16
cat('\n\n')
anova(imp_1000_t)
#> 
#> Family: binomial 
#> Link function: logit 
#> 
#> Formula:
#> p90_open ~ sqrt(pct_l_1000) + s(station, bs = "re")
#> 
#> Parametric Terms:
#>                  df Chi.sq p-value
#> sqrt(pct_l_1000)  1  4.729  0.0297
#> 
#> Approximate significance of smooth terms:
#>              edf Ref.df Chi.sq p-value
#> s(station) 142.9  236.0  394.3  <2e-16
```

Again, only the 1000 meter version shows significant patterns with land
use. The square root transform makes almost no difference regarding the
adequacy of the fit.

``` r
anova(imp_1000, imp_1000_t)
#> Analysis of Deviance Table
#> 
#> Model 1: p90_open ~ pct_l_1000 + s(station, bs = "re")
#> Model 2: p90_open ~ sqrt(pct_l_1000) + s(station, bs = "re")
#>   Resid. Df Resid. Dev      Df Deviance
#> 1    9203.8     4104.9                 
#> 2    9202.9     4103.6 0.84634   1.2469
```

We continue with the untransformed model, as easier to explain to our
audience.

``` r
summary(imp_1000)
#> 
#> Family: binomial 
#> Link function: logit 
#> 
#> Formula:
#> p90_open ~ pct_l_1000 + s(station, bs = "re")
#> 
#> Parametric coefficients:
#>             Estimate Std. Error z value Pr(>|z|)    
#> (Intercept)  3.27367    0.18323  17.867   <2e-16 ***
#> pct_l_1000  -0.05958    0.02430  -2.451   0.0142 *  
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> Approximate significance of smooth terms:
#>              edf Ref.df Chi.sq p-value    
#> s(station) 141.9    236  390.7  <2e-16 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> R-sq.(adj) =  0.0501   Deviance explained = 11.4%
#> UBRE = -0.53266  Scale est. = 1         n = 9399
```

We extract regularly spaced predictions with `emmeans()`. This
simplifies averaging across stations, and correctly incorporates the
station by station errors.

``` r
mm <- emmeans(imp_1000, "pct_l_1000",
              at = list(pct_l_1000 = seq(0, 20, 0.5)), type = 'response')
mms <- summary(mm)  # Summary has class dataframe, making access and display easier.
```

#### Graphic

``` r
rawprobs %>%
  select(pct_l_1000, p_p90_open) %>%
  ggplot(aes(pct_l_1000, p_p90_open)) +
  geom_point(alpha = 0.25, color = cbep_colors()[4]) +

  geom_line(data = mms, mapping = aes(x = pct_l_1000, y = prob)) +
  xlab('Percent Impervious, 1000 m') +
  ylab('Prob. Meeting P90 Standard') +
  theme_cbep(base_size= 12)
```

<img src="frequency-analysis_files/figure-gfm/imperv_graphic-1.png" style="display: block; margin: auto;" />

So, the relationship is statistically significant, but not especially
striking or robust. There is a lot of scatter, and the trend is driven
mostly by (a) the frequency of sites that never see bad water quality,
and (b) a handful of sites with high imperviousness.

#### Cleanup

``` r
rm(imperv_data,
   tmp, plt, mm, mms, p_res, 
   p90_open_glm_1, p90_open_gam_months, p90_open_gam_regions,
   p90_open_gam_rain, imp_100, imp_100_t, imp_500, imp_500_t, imp_1000, imp_1000_t,
   lowerfun)
#> Warning in rm(imperv_data, tmp, plt, mm, mms, p_res, p90_open_glm_1,
#> p90_open_gam_months, : object 'imperv_data' not found
```

## Proportional Odds Models

We explore proportional odds models to see if they provide any further
insight.

### Summary of the Data

``` r
summ <- as_tibble(ftable(xtabs(~ grow_area + all_lvls, data = freq_data))) %>%
  group_by(grow_area) %>%
  mutate(Prop = Freq / sum(Freq)) %>%
  ungroup()
summ
#> # A tibble: 20 x 4
#>    grow_area all_lvls    Freq    Prop
#>    <fct>     <fct>      <int>   <dbl>
#>  1 WI        open        1505 0.812  
#>  2 WJ        open        2694 0.908  
#>  3 WK        open        1814 0.877  
#>  4 WL        open        2287 0.911  
#>  5 WI        gm_relay     131 0.0707 
#>  6 WJ        gm_relay     130 0.0438 
#>  7 WK        gm_relay     106 0.0512 
#>  8 WL        gm_relay     100 0.0398 
#>  9 WI        p90_relay    122 0.0658 
#> 10 WJ        p90_relay     84 0.0283 
#> 11 WK        p90_relay     80 0.0387 
#> 12 WL        p90_relay     71 0.0283 
#> 13 WI        gm_closed     50 0.0270 
#> 14 WJ        gm_closed     24 0.00809
#> 15 WK        gm_closed     35 0.0169 
#> 16 WL        gm_closed     28 0.0112 
#> 17 WI        p90_closed    45 0.0243 
#> 18 WJ        p90_closed    35 0.0118 
#> 19 WK        p90_closed    34 0.0164 
#> 20 WL        p90_closed    24 0.00956
```

``` r
ggplot(summ, aes(x = grow_area, fill = all_lvls, y = Prop)) +
  geom_col(size = .75, position = 'dodge')  +
  theme_cbep(base_size = 12)
```

<img src="frequency-analysis_files/figure-gfm/hist_levels_by_region-1.png" style="display: block; margin: auto;" />

### Region Models

#### Using `VGAM`

VGAM does not allow for random effects. A model fitting station %in%
grow\_area bogged down with high memory usage.

##### Basic Model

The parameter `parallel = TRUE` in the family function sets this up as a
proportional odds model. The cutpoints for all regions are proportional.

``` r
system.time(
  pom_region <- vglm(all_lvls ~ grow_area,
                     data = freq_data,
                     family = cumulative(link = "logitlink",
                                         parallel = TRUE,
                                         reverse = FALSE))
)
#>    user  system elapsed 
#>    0.20    0.01    0.22
```

``` r
summary(pom_region)
#> 
#> Call:
#> vglm(formula = all_lvls ~ grow_area, family = cumulative(link = "logitlink", 
#>     parallel = TRUE, reverse = FALSE), data = freq_data)
#> 
#> Coefficients: 
#>               Estimate Std. Error z value Pr(>|z|)    
#> (Intercept):1  1.45814    0.05910  24.671  < 2e-16 ***
#> (Intercept):2  2.07309    0.06396  32.413  < 2e-16 ***
#> (Intercept):3  2.95141    0.07802  37.830  < 2e-16 ***
#> (Intercept):4  3.65817    0.09836  37.191  < 2e-16 ***
#> grow_areaWJ    0.83526    0.08672   9.632  < 2e-16 ***
#> grow_areaWK    0.50100    0.08900   5.629 1.81e-08 ***
#> grow_areaWL    0.87123    0.09165   9.506  < 2e-16 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> Names of linear predictors: logitlink(P[Y<=1]), logitlink(P[Y<=2]), 
#> logitlink(P[Y<=3]), logitlink(P[Y<=4])
#> 
#> Residual deviance: 9405.086 on 37589 degrees of freedom
#> 
#> Log-likelihood: -4702.543 on 37589 degrees of freedom
#> 
#> Number of Fisher scoring iterations: 4 
#> 
#> Warning: Hauck-Donner effect detected in the following estimate(s):
#> '(Intercept):4'
#> 
#> 
#> Exponentiated coefficients:
#> grow_areaWJ grow_areaWK grow_areaWL 
#>    2.305407    1.650378    2.389837
```

##### Model Framework

The estimated model can be written as:

$$
\\begin{aligned}
logit( \\text{open \| gm\_relay}) = logit( \\hat{P}(Y &lt;= 1)) =  1.46 + \\beta\_i \\\\
logit( \\text{gm\_relay \| p90\_relay}) = logit( \\hat{P}(Y &lt;= 2)) =  2.07 + \\beta\_i \\\\
logit( \\text{p90\_relay \| gm\_closed}) = logit( \\hat{P}(Y &lt;= 3)) =  2.95 + \\beta\_i \\\\
logit( \\text{gm\_closed \| p90\_closed}) = logit( \\hat{P}(Y &lt;= 4)) =  3.66 + \\beta\_i \\\\
\\end{aligned}
$$

Where the *β*<sub>*i*</sub> are log odds between the base case (here
“WI”) and the other cases. It is hard to interpret those coefficients in
the absence of the other coefficients.

Each of the “cutpoints” (shown numerically) is the log odds of being
below the threshold between levels. There for the probability of having
better water quality.

##### Alternate (Non-proportional) Model

``` r
system.time(
  om_region<- vglm(all_lvls ~ grow_area,
                     data = freq_data,
                     family = cumulative(link = "logitlink",
                                         parallel = FALSE,
                                         reverse = FALSE))
)
#>    user  system elapsed 
#>    0.23    0.00    0.23
```

##### Compare Models

``` r
AIC(pom_region)
#> [1] 9419.086
AIC(om_region)
#> [1] 9429.989
```

The proportional odds model is better by AIC, as it uses fewer
parameters, and does nearly as good a job fitting the data.

##### Interpretation

We can look at predicted (conditional) probabilities. We work through
the logic of the models, to clarify the structure of the model
predictions.

First, we calculate the log odds ratios. This is the linear predictor
generated by the model. These are log odds of being below each
“cutpoint.”

``` r
pp <- predict(pom_region, 
              newdata = data.frame(grow_area = c('WI', 'WJ', 'WK', 'WL')))
pp
#>   logitlink(P[Y<=1]) logitlink(P[Y<=2]) logitlink(P[Y<=3]) logitlink(P[Y<=4])
#> 1           1.458137           2.073090           2.951410           3.658174
#> 2           2.293394           2.908347           3.786668           4.493431
#> 3           1.959141           2.574094           3.452415           4.159178
#> 4           2.329362           2.944315           3.822636           4.529399
```

Then we convert from log odds to probabilities of being below each
threshold.

``` r
pp <- inv_logit(pp)
colnames(pp) <- c('p[<=1]', 'p[<=2]', 'p[<=3]', 'p[<=5]')
rownames(pp) <- c('WI', 'WJ', 'WK', 'WL')
pp
#>       p[<=1]    p[<=2]    p[<=3]    p[<=5]
#> WI 0.8112475 0.8882600 0.9503301 0.9748683
#> WJ 0.9083284 0.9482575 0.9778316 0.9889415
#> WK 0.8764400 0.9291756 0.9693031 0.9846199
#> WL 0.9112797 0.9499941 0.9785980 0.9893280
```

And finally, we find the differences between successive threshold
probabilities to generate probabilities of observing each outcome.

``` r
pp <- cbind(pp, 1)
ppp <- pp
colnames(ppp) <- c('p = 1', 'p = 2', 'p = 3', 'p = 4', 'p = 5')
for (i in 2:5)
  ppp[,i] <- pp[,i] - pp[,i-1]

ppp
#>        p = 1      p = 2      p = 3      p = 4      p = 5
#> WI 0.8112475 0.07701248 0.06207013 0.02453823 0.02513166
#> WJ 0.9083284 0.03992907 0.02957405 0.01110989 0.01105855
#> WK 0.8764400 0.05273563 0.04012749 0.01531678 0.01538014
#> WL 0.9112797 0.03871435 0.02860388 0.01072998 0.01067203
```

You can spit out the same type of prediction of probabilities from
`VGAM` using parameter `type = 'response'`.

We use that to look at results from the complete odds model, which are
nearly identical, as suggested by the very similar AIC values we
calculated before.

``` r
pp <- predict(om_region, 
              newdata = data.frame(grow_area = c('WI', 'WJ', 'WK', 'WL')),
              type = 'response')
pp
#>        open   gm_relay  p90_relay   gm_closed  p90_closed
#> 1 0.8121964 0.07069617 0.06583918 0.026983271 0.024284943
#> 2 0.9079879 0.04381530 0.02831143 0.008088979 0.011796427
#> 3 0.8767521 0.05123248 0.03866602 0.016916385 0.016433059
#> 4 0.9111554 0.03984064 0.02828685 0.011155378 0.009561753
```

### Using MASS

We can fit the same model using the `polr()` function from `MASS`. The
primary difference are:  
1. `polr()` can not fit the non-proportional odds model we fit as an
alternate in `vglm()`.  
2. `polr()` parameterizes the model in a different way, so the model
parameters are the negative of what was generated in `vglm()`.  
3. `MASS` provides the option of producing predicted probabilities with
the `type = 'p'` parameter to the
‘predict.polr()`function, rather than the`type =
’response’`parameter used by`VGAM`.   4.`vglm()\` has many more models
and alternate forms for presenting the data, making it both more
flexible, and easier to misapply.

``` r
system.time(
  polr_region <- polr(all_lvls ~ grow_area,
                    data = freq_data, Hess = TRUE,
                    method = "logistic")
)
#>    user  system elapsed 
#>    0.22    0.00    0.22
```

``` r
summary(polr_region)
#> Call:
#> polr(formula = all_lvls ~ grow_area, data = freq_data, Hess = TRUE, 
#>     method = "logistic")
#> 
#> Coefficients:
#>               Value Std. Error t value
#> grow_areaWJ -0.8352    0.08673  -9.630
#> grow_areaWK -0.5010    0.08916  -5.619
#> grow_areaWL -0.8712    0.09170  -9.501
#> 
#> Intercepts:
#>                      Value   Std. Error t value
#> open|gm_relay         1.4581  0.0593    24.6091
#> gm_relay|p90_relay    2.0730  0.0641    32.3481
#> p90_relay|gm_closed   2.9514  0.0781    37.8073
#> gm_closed|p90_closed  3.6581  0.0984    37.1814
#> 
#> Residual Deviance: 9405.086 
#> AIC: 9419.086
```

MASS provides an option of predicting the probabilities with the
`type = 'p'` parameter to the ’predict.polr()\` function.

``` r
pp <- predict(polr_region, 
              newdata = data.frame(grow_area = c('WI', 'WJ', 'WK', 'WL')),
              type = 'p')
pp
#>        open   gm_relay  p90_relay  gm_closed p90_closed
#> 1 0.8112418 0.07701365 0.06207237 0.02453906 0.02513311
#> 2 0.9083222 0.03993122 0.02957626 0.01111070 0.01105961
#> 3 0.8764401 0.05273503 0.04012766 0.01531676 0.01538046
#> 4 0.9112786 0.03871444 0.02860444 0.01073014 0.01067242
```

### Rainfall Models

We continue by examining models of the impact of precipitation \#\#\#\#
Full model

``` r
system.time(
  pom_rain <- vglm(all_lvls ~ log1precip + log1precip_d1,
                     data = freq_data,
                     family = cumulative(link = "logitlink",
                                         parallel = TRUE,
                                         reverse = FALSE))
)
#>    user  system elapsed 
#>    0.22    0.00    0.22
```

``` r
summary(pom_rain)
#> 
#> Call:
#> vglm(formula = all_lvls ~ log1precip + log1precip_d1, family = cumulative(link = "logitlink", 
#>     parallel = TRUE, reverse = FALSE), data = freq_data)
#> 
#> Coefficients: 
#>               Estimate Std. Error z value Pr(>|z|)    
#> (Intercept):1  2.66924    0.04900   54.48   <2e-16 ***
#> (Intercept):2  3.31354    0.05682   58.31   <2e-16 ***
#> (Intercept):3  4.22463    0.07393   57.14   <2e-16 ***
#> (Intercept):4  4.94574    0.09581   51.62   <2e-16 ***
#> log1precip    -0.39271    0.02628  -14.94   <2e-16 ***
#> log1precip_d1 -0.43803    0.02631  -16.65   <2e-16 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> Names of linear predictors: logitlink(P[Y<=1]), logitlink(P[Y<=2]), 
#> logitlink(P[Y<=3]), logitlink(P[Y<=4])
#> 
#> Residual deviance: 9045.134 on 37582 degrees of freedom
#> 
#> Log-likelihood: -4522.567 on 37582 degrees of freedom
#> 
#> Number of Fisher scoring iterations: 5 
#> 
#> Warning: Hauck-Donner effect detected in the following estimate(s):
#> '(Intercept):1', '(Intercept):3', '(Intercept):4'
#> 
#> 
#> Exponentiated coefficients:
#>    log1precip log1precip_d1 
#>     0.6752261     0.6453073
```

#### Simplified Model

``` r
system.time(
  pom_rain_2 <- vglm(all_lvls ~ log1p(precip + precip_d1),
                     data = freq_data,
                     family = cumulative(link = "logitlink",
                                         parallel = TRUE,
                                         reverse = FALSE))
)
#>    user  system elapsed 
#>     0.2     0.0     0.2
```

``` r
AIC(pom_rain)
#> [1] 9057.134
AIC(pom_rain_2)
#> [1] 9155.063
```

So the full model is preferable.

``` r
summary(pom_rain)
#> 
#> Call:
#> vglm(formula = all_lvls ~ log1precip + log1precip_d1, family = cumulative(link = "logitlink", 
#>     parallel = TRUE, reverse = FALSE), data = freq_data)
#> 
#> Coefficients: 
#>               Estimate Std. Error z value Pr(>|z|)    
#> (Intercept):1  2.66924    0.04900   54.48   <2e-16 ***
#> (Intercept):2  3.31354    0.05682   58.31   <2e-16 ***
#> (Intercept):3  4.22463    0.07393   57.14   <2e-16 ***
#> (Intercept):4  4.94574    0.09581   51.62   <2e-16 ***
#> log1precip    -0.39271    0.02628  -14.94   <2e-16 ***
#> log1precip_d1 -0.43803    0.02631  -16.65   <2e-16 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> Names of linear predictors: logitlink(P[Y<=1]), logitlink(P[Y<=2]), 
#> logitlink(P[Y<=3]), logitlink(P[Y<=4])
#> 
#> Residual deviance: 9045.134 on 37582 degrees of freedom
#> 
#> Log-likelihood: -4522.567 on 37582 degrees of freedom
#> 
#> Number of Fisher scoring iterations: 5 
#> 
#> Warning: Hauck-Donner effect detected in the following estimate(s):
#> '(Intercept):1', '(Intercept):3', '(Intercept):4'
#> 
#> 
#> Exponentiated coefficients:
#>    log1precip log1precip_d1 
#>     0.6752261     0.6453073
```

The easiest way to interpret these results is in terms of how the
probability that a sample meets all standards changes with rainfall. We
calculate probabilities on a grid of rainfall amounts.

``` r
rain_df <- tibble(precip = rep(seq(0,2, 0.25),5),
                      log1precip = log1p(precip),
                      precip_d1 = rep(seq(0,2,0.5),each = 9),
                      log1precip_d1 = log1p(precip_d1))

pp <- predict(pom_rain, 
              newdata = rain_df,
              type = 'response')
rain_df <- rain_df %>% cbind(pp)
```

``` r
rain_df %>%
  
   ggplot(aes(x = precip, y = open, color = factor(precip_d1))) +  
   geom_line() +
   
   scale_y_log10(limits = c(.75,1)) +
  
   
   theme_cbep(base_size = 12) +
   theme(axis.text.x = element_text(angle = 90)) +
  
   labs(color = "Yesterday's\nPrecipitation\n(in)",
       y = 'Probability a Sample Meets All Criteria',
       x = "Today's Precipitation")
```

<img src="frequency-analysis_files/figure-gfm/unnamed-chunk-18-1.png" style="display: block; margin: auto;" />

Or, to put it another way, (and not graphically), rainfall drops the
probability of meeting all criteria from 0.935187 to 0.8527812, by the
same token, it increases the risk of any single observation exceeding
all the thresholds from 0.0070634 to 0.0174111.

``` r
max(rain_df$open)
#> [1] 0.935187
min(rain_df$open)
#> [1] 0.8527812
min(rain_df$p90_closed)
#> [1] 0.0070634
max(rain_df$p90_closed)
#> [1] 0.0174111
```

Expressed in odds, about one on fourteen samples fails the lowest
standard when there has been no recent rain, which after a couple of
days of heavy rain, that climbs to about one in six. Conversely, only
about one in 141 samples fails the highest (P90) standard without recent
rain, which more than doubles to about 1 in 56 after heavy rain.

``` r
max(rain_df$open)/ (1- max(rain_df$open))
#> [1] 14.42901
min(rain_df$open)/ (1- min(rain_df$open))
#> [1] 5.792613

# Since our highest standard is an "exceeds" probability, we calculate the
# inverse odds for interpretation purposes.
(1- max(rain_df$p90_closed))/ max(rain_df$p90_closed)
#> [1] 56.43464
(1- min(rain_df$p90_closed))/min(rain_df$p90_closed)
#> [1] 140.5749
```

Of course, those odds are across all samples, and conditions at some
individual station will be worse. Conversely, many station never fail
standards.

# Further Questions

More complex models (region by rainfall, site by rainfall, weighted
rainfall models, etc.) are possible. Such models could be either
binomial or proportional odds models. However, further analyses are
unlikely to qualitatively alter our understanding of patterns.

The one additional question that could be interesting to examine would
be to look at whether the impact of rainfall differs from Station to
Station or Growing Area to Growing Area.
