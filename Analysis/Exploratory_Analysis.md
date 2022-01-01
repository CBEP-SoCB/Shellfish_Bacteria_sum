Exploratory Data Analysis of Shellfish Sanitation Program Data
================
Curtis C. Bohlen, Casco Bay Estuary Partnership.
11/14/2020

<img
    src="https://www.cascobayestuary.org/wp-content/uploads/2014/04/logo_sm.jpg"
    style="position:absolute;top:10px;right:50px;" />

# Load Libraries

``` r
library(readr)
library(tidyverse)
#> Warning: package 'tidyverse' was built under R version 4.0.5
#> -- Attaching packages --------------------------------------- tidyverse 1.3.1 --
#> v ggplot2 3.3.5     v dplyr   1.0.7
#> v tibble  3.1.6     v stringr 1.4.0
#> v tidyr   1.1.4     v forcats 0.5.1
#> v purrr   0.3.4
#> Warning: package 'ggplot2' was built under R version 4.0.5
#> Warning: package 'tidyr' was built under R version 4.0.5
#> Warning: package 'dplyr' was built under R version 4.0.5
#> Warning: package 'forcats' was built under R version 4.0.5
#> -- Conflicts ------------------------------------------ tidyverse_conflicts() --
#> x dplyr::filter() masks stats::filter()
#> x dplyr::lag()    masks stats::lag()
library(corrr)  # Used for correlate(), which produces a data frame
library(GGally)
#> Warning: package 'GGally' was built under R version 4.0.5
#> Registered S3 method overwritten by 'GGally':
#>   method from   
#>   +.gg   ggplot2

library(CBEPgraphics)
load_cbep_fonts()
theme_set(theme_cbep())
```

# Load Data

## Folder References

``` r
sibfldnm <- 'Derived_Data'
parent <- dirname(getwd())
sibling <- file.path(parent,sibfldnm)
fn = 'Shellfish data 2015 2018.csv'
path <- file.path(sibling, fn)

#dir.create(file.path(getwd(), 'figures'), showWarnings = FALSE)
```

## Read Data

``` r
coli_data <- read_csv(path)
#> Rows: 10130 Columns: 19
#> -- Column specification --------------------------------------------------------
#> Delimiter: ","
#> chr  (8): Station, GROW_AREA, OpenClosed, WDIR, Tide, Class, CATEGORY, RawColi
#> dbl  (6): ...1, Temp, Sal, ColiScore, YEAR, ColiVal
#> lgl  (2): LCFlag, RCFlag
#> dttm (1): SDateTime
#> date (1): SDate
#> time (1): STime
#> 
#> i Use `spec()` to retrieve the full column specification for this data.
#> i Specify the column types or set `show_col_types = FALSE` to quiet this message.
```

## Data Preparation

### Convert to Factors

``` r
coli_data <- coli_data %>%
  mutate_at(4:8, factor) %>%
  mutate(Class = factor(Class,levels = c( 'A', 'CA', 'CR',
                                         'R', 'P', 'X' ))) %>%
  mutate(Tide = factor(Tide, levels = c("L", "LF", "F", "HF",
                                        "H", "HE", "E", "LE"))) %>%
  mutate(DOY = as.numeric(format(SDate, format = '%j')),
         Month = as.numeric(format(SDate, format = '%m')))
```

# Exploratory Data Analysis

## Missing Data

``` r
with(coli_data, length(ColiVal))
#> [1] 10130
with(coli_data, sum(is.na(ColiVal)))
#> [1] 684
```

What does it mean to have missing values in this context? That is not
entirely clear, but these appear to be samples that are recorded with
all site information, but no sample-related information, so my guess is,
these represent samples that were scheduled but not actually collected.

That suggests we should simply drop these rows as uninformative. We hold
off on doing that for now, as we want to explore the structure of the
data.

## What values are represented in our data?

(Note that `ColiVal` has addressed inconsistent handling of non-detects,
by rebuilding numeric data from the source. `ColiVal` includes censored
values. Censoring is shown with `LCFlag` and `RCFlag`).

``` r
test <- coli_data %>%
  select(ColiVal, RCFlag, LCFlag ) %>%
  group_by(factor(ColiVal)) %>%
  select(-ColiVal) %>%
  rename(ColiVal =`factor(ColiVal)`) %>%
  summarize(some_censored = as.logical(sum(any(RCFlag | LCFlag))),
            .groups = 'drop') %>%
  mutate(row = row_number(),
         ColiVal = as.numeric(as.character(ColiVal)))
  
ggplot(test, aes(row, ColiVal)) +
  geom_point(aes(color = some_censored))
#> Warning: Removed 1 rows containing missing values (geom_point).
```

<img src="Exploratory_Analysis_files/figure-gfm/unnamed-chunk-6-1.png" style="display: block; margin: auto;" />

### Low Values

``` r
test$ColiVal[which(test$ColiVal <= 10)]
#>  [1]  2.0  2.8  3.6  4.0  5.5  6.0  6.9  7.0  7.3  7.7  8.0  9.0  9.1 10.0
```

### Non-integer Values

Non-integer values include the following:

``` r
test$ColiVal[which(test$ColiVal != as.integer(test$ColiVal))]
#>  [1]  2.8  3.6  5.5  6.9  7.3  7.7  9.1 10.1 10.9 12.7 14.5 16.4 23.6 27.3 30.9
#> [16] 32.7 34.5 36.4 38.2 41.8
```

These are discrete, low, non-integer values. Presumably these reflect
the possible numerical values derived from the methods used to estimate
number of colony forming units in each water sample.

## Exploratory Graphics

``` r
ggplot(coli_data, aes(SDateTime, ColiVal)) + 
  geom_point(aes(color=LCFlag | RCFlag)) +
  scale_y_log10()
#> Warning: Removed 684 rows containing missing values (geom_point).
```

<img src="Exploratory_Analysis_files/figure-gfm/unnamed-chunk-9-1.png" style="display: block; margin: auto;" />

This shows us: 1. The discrete of values observed 2. The importance of
censoring 3. Possible coding errors where censored values were perhaps
not consistently coded in the original data. 4. Data is highly skewed,
so that even the log of the data remains skewed. A lognormal density is
not appropriate for these data.

We need to return to how we interpreted the raw data, in particular
understanding how we ended up with nominal uncensored values at the
censored value of 1.9 and 2.0. (This is consistent with a pattern
observed in the related Beaches data, where we had an excess of
observations at the lower reporting limit).

That may reflect the way we determined whether a sample was censored –
by looking for a less than sign (“&lt;”) in the original data. We could
have alternatively looked at whether the raw col and coliscore values
were identical.

So let’s go to the the trouble of reshaping our data the way we really
want.

``` r
coli_data %>%
  filter(! is.na(ColiVal)) %>%
  group_by(GROW_AREA, YEAR) %>%
  summarize(gmColi = exp(mean(log(ColiVal))),
            .groups = 'drop') %>%
  ggplot(aes(YEAR, gmColi, color=GROW_AREA)) + geom_line(lwd=2)
```

<img src="Exploratory_Analysis_files/figure-gfm/unnamed-chunk-10-1.png" style="display: block; margin: auto;" />

So, geometric means vary year to year and Growing Area to Growing Area.
We obviously need indications of variability to raw any conclusions, but
this points in useful directions.

## Histograms

First, a histogram of the log of the fecal coliform numbers.

``` r
ggplot(coli_data, aes(log10(ColiVal))) +
  geom_histogram(aes(color=RCFlag | LCFlag), bins = 100)
#> Warning: Removed 684 rows containing non-finite values (stat_bin).
```

<img src="Exploratory_Analysis_files/figure-gfm/unnamed-chunk-11-1.png" style="display: block; margin: auto;" />

And then a log-log histogram. A Pareto random variable produces a linear
relation on in a log-log histogram.

``` r
ggplot(coli_data, aes(ColiVal)) +
  geom_histogram(aes(color=RCFlag | LCFlag), bins = 100) +
  scale_x_log10() + scale_y_log10()
#> Warning: Removed 684 rows containing non-finite values (stat_bin).
#> Warning: Transformation introduced infinite values in continuous y-axis
#> Warning: Removed 115 rows containing missing values (geom_bar).
```

<img src="Exploratory_Analysis_files/figure-gfm/unnamed-chunk-12-1.png" style="display: block; margin: auto;" />

As suggested before, this is a strongly skewed, heavy-tailed
distribution. So, Gamma, Exponential, and perhaps Pareto distributions
might work, but that nearly linear decline int eh histogram suggests a
Pareto-style distribution (with added complications due to censoring).

## Rank Correlations

``` r
coli_data %>% 
  select(Temp, Sal, ColiVal, YEAR) %>%
  correlate(method = 'spearman')
#> 
#> Correlation method: 'spearman'
#> Missing treated using: 'pairwise.complete.obs'
#> # A tibble: 4 x 5
#>   term       Temp     Sal ColiVal    YEAR
#>   <chr>     <dbl>   <dbl>   <dbl>   <dbl>
#> 1 Temp    NA       0.219   0.192   0.0294
#> 2 Sal      0.219  NA      -0.0801  0.168 
#> 3 ColiVal  0.192  -0.0801 NA       0.0117
#> 4 YEAR     0.0294  0.168   0.0117 NA
```

So, no strong linear rank correlations. Correlations are slightly lower
using Pearson correlations. Weak correlations with Temperature and
Salinity and salinity and fecal coliform are likely meaningful. Also,
there may be non-linear relationships at play. The salinity to year
correlation is unexpected….

``` r
coli_data %>% 
  select(Temp, Sal, ColiVal, YEAR) %>%
  mutate(ColiVal = log10(ColiVal)) %>%
  ggpairs(lower = list(continuous=wrap('smooth_loess', color='red')),
          progress=FALSE)
#> Warning: Removed 829 rows containing non-finite values (stat_density).
#> Warning in ggally_statistic(data = data, mapping = mapping, na.rm = na.rm, :
#> Removed 868 rows containing missing values
#> Warning in ggally_statistic(data = data, mapping = mapping, na.rm = na.rm, :
#> Removed 849 rows containing missing values
#> Warning in ggally_statistic(data = data, mapping = mapping, na.rm = na.rm, :
#> Removed 829 rows containing missing values
#> Warning: Removed 868 rows containing non-finite values (stat_smooth).
#> Warning: Removed 868 rows containing missing values (geom_point).
#> Warning: Removed 707 rows containing non-finite values (stat_density).
#> Warning in ggally_statistic(data = data, mapping = mapping, na.rm = na.rm, :
#> Removed 708 rows containing missing values
#> Warning in ggally_statistic(data = data, mapping = mapping, na.rm = na.rm, :
#> Removed 707 rows containing missing values
#> Warning: Removed 849 rows containing non-finite values (stat_smooth).
#> Warning: Removed 849 rows containing missing values (geom_point).
#> Warning: Removed 708 rows containing non-finite values (stat_smooth).
#> Warning: Removed 708 rows containing missing values (geom_point).
#> Warning: Removed 684 rows containing non-finite values (stat_density).
#> Warning in ggally_statistic(data = data, mapping = mapping, na.rm = na.rm, :
#> Removed 684 rows containing missing values
#> Warning: Removed 829 rows containing non-finite values (stat_smooth).
#> Warning: Removed 829 rows containing missing values (geom_point).
#> Warning: Removed 707 rows containing non-finite values (stat_smooth).
#> Warning: Removed 707 rows containing missing values (geom_point).
#> Warning: Removed 684 rows containing non-finite values (stat_smooth).
#> Warning: Removed 684 rows containing missing values (geom_point).
```

<img src="Exploratory_Analysis_files/figure-gfm/unnamed-chunk-14-1.png" style="display: block; margin: auto;" />

So, data are very highly scattered, obscuring patterns even when the
fecal coliform counts are log transformed. It appears fecal coliform
Levels are weakly related to salinity and temperature, which are not
especially correlated themselves. The highly skewed fecal coliform data
is problematic.

``` r
coli_data %>%
  select(DOY, Month, YEAR, ColiVal) %>%
  mutate(ColiVal = log10(ColiVal)) %>%
  ggpairs(lower = list(continuous=wrap('smooth_loess', color='red')), progress=FALSE) 
#> Warning in ggally_statistic(data = data, mapping = mapping, na.rm = na.rm, :
#> Removed 684 rows containing missing values

#> Warning in ggally_statistic(data = data, mapping = mapping, na.rm = na.rm, :
#> Removed 684 rows containing missing values

#> Warning in ggally_statistic(data = data, mapping = mapping, na.rm = na.rm, :
#> Removed 684 rows containing missing values
#> Warning: Removed 684 rows containing non-finite values (stat_smooth).
#> Warning: Removed 684 rows containing missing values (geom_point).
#> Warning: Removed 684 rows containing non-finite values (stat_smooth).
#> Warning: Removed 684 rows containing missing values (geom_point).
#> Warning: Removed 684 rows containing non-finite values (stat_smooth).
#> Warning: Removed 684 rows containing missing values (geom_point).
#> Warning: Removed 684 rows containing non-finite values (stat_density).
```

<img src="Exploratory_Analysis_files/figure-gfm/unnamed-chunk-15-1.png" style="display: block; margin: auto;" />

No strong patterns with time of year, month, or year, except weak
relationships to time of year. Obvious artifacts due to the discrete
nature of values of the fecal coliform data at low values. Fewer
late-year observations from 2019, so we probably need to either model
time of year or remove 2019 from the data entirely.

# Statistical Considerations

The highly skewed and heavy-tailed nature of the distribution suggests
that any form of linear model or generalized linear model to predict the
observed values will be difficult.

The data is distributed approximately Pareto, so perhaps it could be
modeled with a Pareto GLM (over a fixed support of perhaps 0:inf, or
1:inf). The package VGAM supports such Pareto GLMs.

It may be more appropriate to transform observations into exceedances of
one or more thresholds, and analyze the probability of exceedances with
a binomial or quasi-binomial GLM.

A final alternative may be to use a one-dimensional PERMANOVA to analyze
the observed data. Using euclidean distance metric, this is equivalent
to an ANOVA analysis, but with standard errors and significance levels
driven by permutation tests rather than assumptions of normality. I am
not certain how well that would work with these data’s heavy tails.

If we are willing to forego the modeling sophistication possible with
linear models, we might be able to test for differences in medians by
classes using a rank-based procedure, such as Wilcoxon / Kruskal-Wallis,
etc.

## Why are we doing this analysis anyway?

Selection of methods here hinges on our analytic goals. To some extent,
those goals are unclear, which complicates planning of the analysis.
Results of analysis will end up as one or two graphics in State of the
Bay. The most important graphic will be a map of some sort of levels of
risk or concern. Another graphic could be used to show off interesting
relationships, such has between levels of contamination and time of
year, temperature, or salinity. A third possibility is that results
could be used to create some sort of geospatial analysis looking for
correlations with nearby land use.

In the State of the Bay Report we want to convey certain messages: 1.
Some sites are more likely to have high bacteria counts than others –
and we’d like to show that as a map. 2. We’d like to be able to give
some explanation for some of those differences. 3. We may want to show
relationships to other predictors, including time of year, rainfall,
salinity, temperature, land use etc. 4. We’d like to see whether sites
falling under different classifications have expected differences in
probability of extreme values.

So, our core need is for some sort of summary statistic for individual
sites that we can use to show geographic patterns. The second need is
for models that allow us to explore the impact of both categorical and
quantitative predictors.

Here are the criteria used for establishing shellfish harvest area

# Relevant Standards

## Growing Area Classification Standards

| Growing Area Classification | Activity Allowed                                                          | Geometric mean FC/100ml | 90th Percentile (P90) FC/100ml |
|-----------------------------|---------------------------------------------------------------------------|-------------------------|--------------------------------|
| Approved                    | Harvesting allowed                                                        | ≤ 14                    | ≤ 31                           |
| Conditionally Approved      | Harvesting allowed except during specified conditions                     | ≤ 14 in open status     | ≤ 31 in open status            |
| Restricted                  | Depuration harvesting or relay only                                       | ≤ 88 and &gt;15         | ≤ 163 and &gt;31               |
| Conditionally Restricted    | Depuration harvesting or relay allowed except during specified conditions | ≤ 88 in open status     | ≤ 163 in open status           |
| Prohibited                  | Aquaculture seed production only                                          | &gt;88                  | &gt;163                        |

So, critical levels for Geometric Mean include:  &lt;  = 14 and
 &lt;  = 88 and for the p90  &lt; 31 and  &lt;  = 163

## Maine State Class SB Waters Standards

Maine’s water quality criteria includes an additional standard, which
applies only indirectly to these data:  
&gt; the number of enterococcus bacteria in these waters may not exceed
a geometric mean of 8 CFU per 100 milliliters in any 90-day interval or
54 CFU per 100 milliliters in more than 10% of the samples in any 90-day
interval.

38 M.R.S. §465-B(2)(B)

A “90 day interval” might apply to a summer’s worth of data, but in most
years that will only represent a handful of observations at each site.
Also note that this standard is written in terms of “enterococci”, not
“fecal coliform” or “coliformes”.

## Evaluation

We can readily calculate a geometric mean for each site – which is the
basis on which shellfish area closures are calculated – but we can not
readily model geometric means for this distribution. The most straight
forward way to do that would be to analyze the log of counts, but even
the log of the raw data is highly skewed, and heavy tailed.
