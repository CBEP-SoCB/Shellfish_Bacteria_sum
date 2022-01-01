Examining Distribution of Shellfish Pathogen Data
================
Curtis C. Bohlen, Casco Bay Estuary Partnership.
02/10/2021

-   [Introduction](#introduction)
-   [Load Libraries](#load-libraries)
-   [Load Data](#load-data)
    -   [Address Censored Data](#address-censored-data)
    -   [Remove NAs](#remove-nas)
-   [Dealing with Skewed Data](#dealing-with-skewed-data)
    -   [Analytic Graphics](#analytic-graphics)
    -   [Gamma Distribution](#gamma-distribution)
    -   [Pareto Distribution](#pareto-distribution)
        -   [What is the Pareto
            Distribution?](#what-is-the-pareto-distribution)
        -   [Means and Variances](#means-and-variances)
        -   [Pareto Distribution Functions in
            R](#pareto-distribution-functions-in-r)
        -   [Effect of Parameters on Shape of the
            Distribution](#effect-of-parameters-on-shape-of-the-distribution)
        -   [Fitting the Pareto
            Distribution](#fitting-the-pareto-distribution)
        -   [Compare Fits](#compare-fits)
        -   [Fitting a Pareto Distribution to Censored
            Data](#fitting-a-pareto-distribution-to-censored-data)
        -   [Means and Variances](#means-and-variances-1)
-   [Quantile matching](#quantile-matching-1)
-   [Maximum Likelihood with Censored
    Data](#maximum-likelihood-with-censored-data)
-   [Conclusions](#conclusions)

<img
    src="https://www.cascobayestuary.org/wp-content/uploads/2014/04/logo_sm.jpg"
    style="position:absolute;top:10px;right:50px;" />

# Introduction

As part of exploratory data analysis, and specifically to provide
insight into how best to handle non-detects and censored values, we want
to understand the distribution of observations in the fecal coliforms
shellfish bacteria data. Here we take a look at overall data
distributions, to consider possible modeling strategies. Since the data
distribution may also be shaped by the distribution of predictors, this
analysis is exploratory only, and any distributional assumptions need to
be confirmed in the model context.

# Load Libraries

``` r
library(readr)
library(fitdistrplus)  # For cullen-fray graph etc.
#> Warning: package 'fitdistrplus' was built under R version 4.0.5
#> Loading required package: MASS
#> Loading required package: survival
library(actuar)        # For a particular version of Pareto Distribution fxns
#> Warning: package 'actuar' was built under R version 4.0.5
#> 
#> Attaching package: 'actuar'
#> The following object is masked from 'package:grDevices':
#> 
#>     cm
library(tidyverse)  # Loads another `select()`
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
#> x dplyr::select() masks MASS::select()

library(emmeans)   # For marginal means
#> Warning: package 'emmeans' was built under R version 4.0.5
#> 
#> Attaching package: 'emmeans'
#> The following object is masked from 'package:actuar':
#> 
#>     emm

library(mblm)      # for the Thiel-Sen estimators -- not really successful here

library(VGAM)      # For Pareto GLMs and estimation.
#> Loading required package: stats4
#> Loading required package: splines
#> 
#> Attaching package: 'VGAM'
#> The following object is masked from 'package:tidyr':
#> 
#>     fill
#> The following objects are masked from 'package:actuar':
#> 
#>     dgumbel, dlgamma, dpareto, pgumbel, plgamma, ppareto, qgumbel,
#>     qlgamma, qpareto, rgumbel, rlgamma, rpareto

library(GGally)
#> Warning: package 'GGally' was built under R version 4.0.5
#> Registered S3 method overwritten by 'GGally':
#>   method from   
#>   +.gg   ggplot2
#> 
#> Attaching package: 'GGally'
#> The following object is masked from 'package:emmeans':
#> 
#>     pigs

library(CBEPgraphics)
load_cbep_fonts()
theme_set(theme_cbep())

library(LCensMeans)
```

# Load Data

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
  mutate_at(c(4:8), factor) %>%
  mutate(Class = factor(Class, levels = c( 'A', 'CA', 'CR',
                                           'R', 'P', 'X' ))) %>%
  mutate(Tide = factor(Tide, levels = c("L", "LF", "F", "HF",
                                        "H", "HE", "E", "LE"))) %>%
  mutate(DOY = as.numeric(format(SDate, format = '%j')),
         Month = as.numeric(format(SDate, format = '%m')))
```

## Address Censored Data

Right censored values are sufficiently rare as to be relatively
unimportant we leave the undressed, but address left censored values.
Interval Censored values add complexity, but are unlikely in our setting
to sharply alter qualitative conclusions, so we chose not to address
them.

``` r
coli_data <- coli_data %>%
  mutate(ColiVal_2 = sub_cmeans(ColiVal, LCFlag))
```

``` r
ggplot(coli_data, aes(x = ColiVal, y = ColiVal_2, color = LCFlag)) +
  geom_point() +
  xlim(0,25) +
  ylim(0,25)
#> Warning: Removed 1405 rows containing missing values (geom_point).
```

<img src="dist_bacteria_files/figure-gfm/unnamed-chunk-4-1.png" style="display: block; margin: auto;" />

Almost all censored values were at 2.

``` r
coli_data %>%
  filter(ColiVal == 2, LCFlag) %>%
  pull(ColiVal_2) %>%
  summary
#>    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#>  0.5491  0.6000  0.6116  0.6115  0.6229  0.6736
```

So, our (lognormal) based estimator for censored values estimates a
non-detect at somewhere around 0.61 CF per 100 ml.

## Remove NAs

``` r
coli_data <- coli_data %>%
  filter (! is.na(ColiVal))
```

# Dealing with Skewed Data

The data is so highly skewed that even the logs are skewed, which
suggests we may have a hard time finding a parametric way of analyzing
these data.

Given the non-detects and censored upper values, we can consider 1.
Fitting data to a parametric family of distributions, accounting for
censoring  
2. Robust or resistant statistical methods  
3. Non-parametric methods

## Analytic Graphics

The package ‘fitdistrplus’ provides functions for fitting standard as
well as non-standard distributions both with and without censored data.
We test it first without accounting for censoring.

``` r
test_dat <- coli_data %>%
  filter(! is.na(ColiVal)) %>%
  select(ColiVal, LCFlag, RCFlag)
plotdist(test_dat$ColiVal, histo = TRUE, demp = TRUE)
```

<img src="dist_bacteria_files/figure-gfm/unnamed-chunk-7-1.png" style="display: block; margin: auto;" />

The package contains and interesting function that provides a graphical
representation of data in terms of kurtosis and skewness. The output is
helpful for evaluating potential alternative data distributions.

``` r
descdist(test_dat$ColiVal, boot = 1000)
```

<img src="dist_bacteria_files/figure-gfm/unnamed-chunk-8-1.png" style="display: block; margin: auto;" />

    #> summary statistics
    #> ------
    #> min:  2   max:  1600 
    #> median:  2 
    #> mean:  17.01091 
    #> estimated sd:  99.07828 
    #> estimated skewness:  12.25917 
    #> estimated kurtosis:  172.719

That plot confirms that the data has heavier tails than either lognormal
or gamma distributions A Beta distribution is not appropriate, since
Beta variates are bounded between zero and one. We need a highly skewed,
heavy-tailed distribution with a support on the positive real line.

``` r
descdist(log(test_dat$ColiVal), boot = 1000)
```

<img src="dist_bacteria_files/figure-gfm/unnamed-chunk-9-1.png" style="display: block; margin: auto;" />

    #> summary statistics
    #> ------
    #> min:  0.6931472   max:  7.377759 
    #> median:  0.6931472 
    #> mean:  1.281953 
    #> estimated sd:  1.134261 
    #> estimated skewness:  2.346337 
    #> estimated kurtosis:  8.856334

Even log-transformed data is problematic to model. It is significantly
more skew than an exponential or gamma distribution.

## Gamma Distribution

Gamma can be fit by the method of moments, with parameters Shape = *α*
and Scale = *θ*. The mean and variance are as follows:
*μ* = *α**θ*
*σ*<sup>2</sup> = *α**θ*<sup>2</sup>

The estimates of the parameters therefore are:
*θ* = *α*<sup>2</sup>/*μ*

*α* = *μ*/*θ* = *μ*<sup>2</sup>/*σ*<sup>2</sup>

The empirical moments are:

``` r
(a <- with(coli_data, list(mean=mean(ColiVal_2, na.rm=TRUE),
                           var=var(ColiVal_2, na.rm=TRUE))))
#> $mean
#> [1] 16.24422
#> 
#> $var
#> [1] 9839.914
```

So, by the method of moments, we can examine a gamma distribution with
the following parameters:

``` r
(theta = a$var/a$mean)
#> [1] 605.7487
(alpha = a$mean/theta)
#> [1] 0.02681676
```

## Pareto Distribution

### What is the Pareto Distribution?

The Pareto distribution is an extreme value distribution. The Pareto
distribution comes in several flavors. The one we are focusing on is a
three parameter family, with Location, Scale and Shape parameters. The
first is often treated as known in an estimation context. The latter two
must be positive. This version of the Pareto distribution is often
described as the “Pareto(II)” distribution, or when Location == 0, as
the “Lomax” distribution. Wikipedia includes a discussion of this
distribution under “Generalized Pareto Distribution”, but with different
parameterization.

The “classic” Pareto Distribution (a two parameter version of
Pareto(II), with `Location == Scale`) ends up looking linear on a
log-log plot. Other versions are close to linear. This property (shared
with the Gamma distribution) suggests the Pareto may be a good bet for
the bacteria data, which looks similar (if you ignore censoring).

The Pareto distribution can be fit several ways. First and second
moments have relatively simple closed form values, so parameters could
be estimated using either the method of moments, or maximum likelihood
methods.

### Means and Variances

#### Mean

$$E(x) = \\frac{\\sigma}{\\alpha - 1} + \\mu $$

Where parameters are: \* Location = *μ* (The minimum of the
distribution; 0 for a “Lomax” distribution)  
\* Shape = *α* and  
\* Scale = /*s**i**g**m**a*

Mean goes UP with the scale parameter, and DOWN with the shape
parameter.

#### Variance

The variance is a little trickier. Different parameterizations from
different packages and references are confusing. It appears that (using
our notation so far) variance is:

$$ Var(x) = E((x -E(x))^2) = \\frac{\\sigma^2 \\alpha} {(\\alpha - 1)^2(\\alpha - 2))} \\quad (\\alpha &gt; 1)  $$
Variance is undefined for *α* &lt; 1 and infinite for *α* &lt; 2.

Generally, variance goes UP with the scale parameter, and DOWN with the
shape parameter.

### Pareto Distribution Functions in R

The Pareto distribution is not available in vanilla R, although it is
available in many packages. See
(<https://cran.r-project.org/web/views/Distributions.html>). We examined
versions from `actuar`, `EnvStats` and `VGAM`. The implementations are
slightly different. Confusingly, the names of the parameters of the
distribution are inconsistent from package to package, and even within
the `VGAM` package. We use `VGAM` as it includes nice modeling functions
for what amount to Pareto GLM models (although the mathematical details
differ).

We mostly work with the distribution functions from VGAM (note the
default values):

`dparetoI(x, scale = 1, shape = 1, log = FALSE)`

`dparetoII(x, location = 0, scale = 1, shape = 1, log = FALSE)`

Unfortunately, these two functions use the ‘scale’ parameter somewhat
differently. The ‘scale’ parameter from `dparetoI()` is related to the
‘location’ parameter in `dparetoII()`. In each case, it represents the
MINIMUM value of the distribution. In estimation, the minimum value is
often assumed known.

The reason the terminology changes is that the equivalent `dparetoII()`
function to a `dparetoI()` example has `location == scale`, so
`dparetoI()` needs one fewer parameters, and the choice of which name to
carry in the reduced function was essentially arbitrary.

Here’s an example showing the correspondence. Note the linear decline on
a log-log plot. That linear decline is only an identity for a pareto(I)
distribution.

``` r
ggplot() +
  scale_x_log10(limits = c(0.1,100)) +
  scale_y_log10() +
  geom_function(fun = dparetoI,
                args = list(scale = 2, shape = .5 )) +
geom_function(fun = dparetoII,
                args = list(location = 2, scale = 2, shape = .5 ), 
              color = 'red',
              lty = 2, size = 1)
```

<img src="dist_bacteria_files/figure-gfm/unnamed-chunk-12-1.png" style="display: block; margin: auto;" />

The important distinction is to realize that `paretoII()` has more
flexibility in shape, and in particular the relationship between mean
and weight of the extreme tails.

### Effect of Parameters on Shape of the Distribution

We can get a feel for the implications of shape and scale parameters by
plotting Pareto densities.

``` r
ggplot() +
  scale_x_log10(limits = c(1,1000)) +
  #scale_y_log10() +

  theme_cbep(base_size = 10) +
  
  geom_function(fun = dparetoII,
                 args = list(location = 0, scale = 1, shape = 1),
                 color = 'red') +
  geom_function(fun = dparetoII,
                 args = list(location = 0, scale = 2, shape = 1),
                 color = 'orange') +
  geom_function(fun = dparetoII,
                 args = list(location = 0, scale = 5, shape = 1),
                 color = 'yellow') +
  geom_function(fun = dparetoII,
                 args = list(location = 0, scale = 10, shape = 1),
                 color = 'black')
```

<img src="dist_bacteria_files/figure-gfm/pareto_scale_demo-1.png" style="display: block; margin: auto;" />

So, as SCALE goes up, low values drop, mid and higher values rise.

``` r
ggplot() +
  scale_x_log10(limits = c(1,1000)) +
  #scale_y_log10() +

  theme_cbep(base_size = 10) +
  
  geom_function(fun = dparetoII,
                 args = list(location = 0, scale = 5, shape = .1),
                 color = 'red') +
  geom_function(fun = dparetoII,
                 args = list(location = 0, scale = 5, shape = 1),
                 color = 'orange') +
  geom_function(fun = dparetoII,
                 args = list(location = 0, scale = 5, shape = 2),
                 color = 'yellow') +
  geom_function(fun = dparetoII,
                 args = list(location = 0,  scale = 5, shape = 5),
                 color = 'black')
```

<img src="dist_bacteria_files/figure-gfm/pareto_shape_demo-1.png" style="display: block; margin: auto;" />

The SHAPE parameter has to do with how closely the distribution hugs
towards the “origin” – more correctly the location parameter, which is
the minimum.

Low shape spreads things out a bit more. High shape pulls things in
closely.

A classic Pareto distribution often has `location = scale = 1` has
`location = 1`.

We will often use `location = 0`, and not constrain the scale parameter.
In principal, the value of concentration of bacteria has support on the
entire positive real number line. It is only because of the limitations
of our data that we have a positive minimum.

### Fitting the Pareto Distribution

#### Fitting Based on Moments

Some comments on closed-form maximum likelihood estimates based on
sample moments can be found here:
(<https://stats.stackexchange.com/questions/27426/how-do-i-fit-a-set-of-data-to-a-pareto-distribution-in-r>)

Here is the simple function to estimate parameters for a Pareto(I)
distribution. (for non-censored data). It is from the StackExchange
answer referenced above. The first parameter here is a minimum value,
the second, a shape parameter.

``` r
pareto.MLE <- function(X)
{
   X = X[! is.na(X)]
   n <- length(X)
   m <- min(X, na.rm=TRUE)
   a <- n/sum(log(X)-log(m), na.rm=TRUE)
   return( c(m,a)) 
}

(moment_parms <- pareto.MLE(coli_data$ColiVal_2))
#> [1] 0.5490743 0.8126995
```

#### Fitting the Pareto with `VGAM`

We can also fit a Pareto(II) distribution using `VGAM`. We are fitting a
Pareto distribution with minimum value `location = 0.5` to be as
consistent as possible with the prior. Note that the PARAMETERS are not
very similar at all, which point to a challenge to interpreting modeling
results later. We also note that selection of the minimum value has a
large effect on fitted parameters.

``` r
paretofit = vglm(ColiVal_2~ 1, paretoII(location = 0.5) , data = coli_data, maxit = 100)
vgam_parms <- exp(coef(paretofit))
names(vgam_parms) <- c('Scale', 'Shape')
vgam_parms
#>     Scale     Shape 
#> 0.1422401 0.4866275
```

Given scaling issues, it is difficult to compare these three model fits
to data. We (arbitrarily) divide the observed density by 100 to better
line up the models with the observed data. All are reasonable, none are
obviously correct.

#### Fitting with the `fitdistrplus` Package

The `fitdist()` function provides four methods for fitting a
distribution : “mle”, for ‘maximum likelihood estimation’, “mme” for
‘moment matching estimation’, “qme” for ‘quantile matching estimation’
and “mge” for ‘maximum goodness-of-fit estimation’. Not all appear to
work with these data, and these versions of the Pareto distribution
functions, but it is not clear why.

##### MLE Estimation (equivalent to the VGAM Default)

``` r
f_mle <- fitdist(coli_data$ColiVal_2, "paretoII",
                 method='mle',
                 start = list(shape=1, scale= 1),
                 fix.arg = list(location = 0.5))
#> Warning in fitdist(coli_data$ColiVal_2, "paretoII", method = "mle", start =
#> list(shape = 1, : The dparetoII function should return a zero-length vector when
#> input has length zero
#> Warning in fitdist(coli_data$ColiVal_2, "paretoII", method = "mle", start =
#> list(shape = 1, : The pparetoII function should return a zero-length vector when
#> input has length zero
cat('\n')
(fd_parms <- f_mle$estimate[c(2,1)])
#>     scale     shape 
#> 0.1422302 0.4865680
```

##### Quantile Matching

Quantile matching also works, with a list of quantiles to be matched.
Since we have two parameters to fit, we can ask for fits at two
quantiles. It’s not obvious statistically what quantiles make the most
sense, but the shellfish sanitation program is interested in the
geometric mean and the p90, Here we could go with the median and the
90th percentile.

But note : (1) This does not work with `ColiVal_2`, only with `ColiVal`,
and (2) resulting parameter estimates look quite different.

``` r
f_qme <- fitdist(test_dat$ColiVal, "paretoII",
                 method='qme', 
                 probs=c(.50, .90),
              start = list(shape=1, scale=1),
              fix.arg = list(location = 0.5))
#> Warning in fitdist(test_dat$ColiVal, "paretoII", method = "qme", probs =
#> c(0.5, : The dparetoII function should return a zero-length vector when input
#> has length zero
#> Warning in fitdist(test_dat$ColiVal, "paretoII", method = "qme", probs =
#> c(0.5, : The pparetoII function should return a zero-length vector when input
#> has length zero
#> Warning in fitdist(test_dat$ColiVal, "paretoII", method = "qme", probs =
#> c(0.5, : The qparetoII function should return a zero-length vector when input
#> has length zero
cat('\n')
(fd_parms_alt <- f_qme$estimate[c(2,1)])
#>     scale     shape 
#> 1.1340503 0.8225913
```

The parameters LOOK quite different from the prior estimate, which is
troubling.

### Compare Fits

``` r
ggplot(coli_data, aes(ColiVal_2)) +

  # geom_histogram(data = coli_data, 
  #               mapping=aes(x=ColiVal_2, y=..density../50), 
  #               fill = NA, color = 'gray20') +
  scale_x_log10(limits = c(.1,100) ) +
 # scale_y_log10(limits = c(10^-4, 2)) +

  theme_cbep(base_size = 10) +
  
  geom_function(fun = dgamma,
                args = list(shape=alpha, scale=theta),
                color = 'orange') + 
  geom_function(fun = dparetoII,
                args = list(location = moment_parms[[1]], scale = moment_parms[[1]],
                             shape = moment_parms[[2]]),
                color = 'red') +
  geom_function(fun = dparetoII,
                args = list(location = .5, scale = vgam_parms[[1]],
                             shape = vgam_parms[[2]]),
                color = 'purple', lwd = 2) +
  geom_function(fun = dparetoII,
                args = list(location = 0.5, scale = fd_parms[[1]],
                             shape = fd_parms[[2]]),
                color = 'yellow', lwd = 1) +
  geom_function(fun = dparetoII,
                args = list(location = 0.5, scale = fd_parms_alt[[1]],
                             shape = fd_parms_alt[[2]]),
                color = 'green')
```

<img src="dist_bacteria_files/figure-gfm/plot_pareto_fit-1.png" style="display: block; margin: auto;" />

All the different Pareto distribution fits are similar. Two are
effectively identical. One lesson is that apparent differences in
parameters may have little effect on the real distribution of
probability mass density.

The gamma distribution (orange) is lower, and has significant mass at
low (below 1) values.

### Fitting a Pareto Distribution to Censored Data

The FAQ file for the `fitdistplus` files says the following about
fitting censored data with the `fitdistcens` function):

> Censored data must be represented in the package by a dataframe of two
> columns respectively named left and right, describing each observed
> value as an interval. The left column contains either NA for left
> censored observations, the left bound of the interval for interval
> censored observations, or the observed value for non-censored
> observations. The right column contains either NA for right censored
> observations, the right bound of the interval for interval censored
> observations, or the observed value for non-censored observations.

The function uses maximum likelihood to estimate parameters.

So, we need to address several cases here for censoring: 1. Values below
2.0 2. Discrete nature of values reported by the sampling method,
especially at lower values 3. Right censored values that exceed the
method upper limit

``` r
cens_data <- coli_data %>%
  select(ColiVal, LCFlag, RCFlag) %>%
  filter(! is.na(ColiVal)) %>%
  mutate(left  = ifelse(LCFlag, NA, ColiVal)) %>%
  mutate(right = ifelse(RCFlag, NA, ColiVal)) %>%
  select(left,right)
```

We use the `pareto2` family of functions from `actuar`, rather than the
`paretoII()` functions from `VGAM`, because `actuar` includes a defined
moment function, `mpareto2()`, and quantile function, `qpareto2()`.
Those functions are required for several of the alternative fitting
algorithms in `fitdistplus`.

``` r
ff <- fitdistcens(data.frame(cens_data), "pareto2",
                 start = list(shape=1, scale= 1),
                 fix.arg = list(min = 0.5))
cat('\n')
(cens_parms <- ff$estimate[c(2,1)])
#>     scale     shape 
#> 0.7063829 0.7401517
```

The parameters are quite different from what we observed fitting to our
data without addressing non-detects directly. (We DID fit to data with
non-detects replaced by estimated conditional means, but those estimated
conditional means are based on assuming a lognormal distribution, and
thus probably overestimate the conditional means.)

How different are our estimates with and without explicit modeling of
non-detects?

``` r
ggplot(coli_data, aes(ColiVal_2)) +

  # geom_histogram(data = coli_data, 
  #               mapping=aes(x=ColiVal_2, y=..density../50), 
  #               fill = NA, color = 'gray20') +
  scale_x_log10(limits = c(.1,100) ) +
  # scale_y_log10(limits = c(10^-4, 3)) +

  theme_cbep(base_size = 10) +

  geom_function(fun = dparetoII,
                 args = list(location = 0.5, scale = fd_parms[[1]],
                             shape = fd_parms[[2]]),
                color = 'blue', lwd = 1) +
  geom_function(fun = dparetoII,
                args = list(location = 0.5, scale = fd_parms_alt[[1]],
                             shape = fd_parms_alt[[2]]),
                color = 'purple') +
  geom_function(fun = dparetoII,
                 args = list(location = 0.5, scale = cens_parms[[1]],
                             shape =cens_parms[[2]]),
                color = 'red')
```

<img src="dist_bacteria_files/figure-gfm/plot_pareto_fit_2-1.png" style="display: block; margin: auto;" />

So, explicit modeling of non-detects flattens out the curve a fair
amount, with the largest effects reducing probability at low values and
increasing them at slightly higher values. The revised curve looks a LOT
like the curve fit to the median and 90th percentile.

### Means and Variances

These don’t work very well for these estimated sets of parameters
because estimated parameters are low enough so both the mean and the
variance of the fitted Pareto distribution are unbounded.

``` r
est_pareto_mean <- function(scale, shape, location){
  if_else(shape > 1,
          location + (scale / (shape - 1)),
          NaN)
}
est_pareto_variance <- function(scale, shape, location){
  if_else(shape > 1, 
          if_else(shape>2,
                  scale^2 * shape / ((shape - 1)^2 * (shape - 2)),
                  Inf),
          NaN)
}
```

#### Maximum Likelihood

``` r
fd_parms
#>     scale     shape 
#> 0.1422302 0.4865680
est_pareto_mean(fd_parms[[1]], fd_parms[[2]], 0.5)
#> [1] NaN
est_pareto_variance(fd_parms[[1]], fd_parms[[2]], 0.5)
#> [1] NaN
```

# Quantile matching

``` r
fd_parms_alt
#>     scale     shape 
#> 1.1340503 0.8225913
est_pareto_mean(fd_parms_alt[[1]], fd_parms_alt[[2]], 0.5)
#> [1] NaN
est_pareto_variance(fd_parms_alt[[1]], fd_parms_alt[[2]], 0.5)
#> [1] NaN
```

# Maximum Likelihood with Censored Data

``` r
cens_parms
#>     scale     shape 
#> 0.7063829 0.7401517
est_pareto_mean(cens_parms[[1]], cens_parms[[2]], 0.5)
#> [1] NaN
est_pareto_variance(cens_parms[[1]], cens_parms[[2]], 0.5)
#> [1] NaN
```

# Conclusions

While it is clear that the bacteria data has much in common with a
Pareto distribution, procedures for fitting that distribution to data
are convoluted, and results difficult to interpret, as there is no
simple relationship between parameters and simple summary statistics.
Further, fitted parameters are strongly influenced by analytic choices,
such as selecting which version of the Pareto distribution to fit, and
setting the Location parameter.
