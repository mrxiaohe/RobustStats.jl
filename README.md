RobustStats
==============


[![Build Status](https://travis-ci.org/maximsch2/RobustStats.jl.svg?branch=master)](https://travis-ci.org/maximsch2/RobustStats.jl)

This package contains a variety of functions from the field robust statistical methods. Many are estimators of location or dispersion; others estimate the standard error or the confidence intervals for the location or dispresion estimators, generally computed by the bootstrap method.

Many functions in this package are based on the R package WRS ([an R-Forge repository](https://r-forge.r-project.org/projects/wrs/)) by [Rand Wilcox](http://dornsife.usc.edu/cf/labs/wilcox/wilcox-faculty-display.cfm). Others were contributed by users as needed. [References](#References) to the statistics literature can be found below.

This package requires `Compat`, `Rmath`, `Dataframes`, and `Distributions`. They can be installed automatically, or by invoking `Pkg.add("packagename")`.

## Estimators

### Location estimators:
* `tmean(x, tr=0.2)`  - Trimmed mean: mean of data with the lowest and highest fraction `tr` of values omitted.
* `winmean(x, tr=0.2)`- Winsorized mean: mean of data with the lowest and highest fraction `tr` of values squashed to the 20%ile or 80%ile value, respectively.
* `tauloc(x)`         - Tau measure of location by Yohai and Zamar.
* `onestep(x)`        - One-step M-estimator of location using Huber's ψ
* `mom(x)`            - Modified one-step M-estimator of location (MOM)
* `bisquareWM(x)`     - Mean with weights given by the bisquare rho function.
* `huberWM(x)`        - Mean with weights given by Huber's rho function.
* `trimean(x)`        - Tukey's trimean, the average of the median and the midhinge.

### Dispersion estimators:
* `winvar(x, tr=0.2)` - Winsorized variance.
* `wincov(x, y, tr=0.2)` - Winsorized covariance.
* `pbvar(x)`          - Percentage bend midvariance.
* `bivar(x)`          - Biweight midvariance.
* `tauvar(x)`         - Tau measure of scale by Yohai and Zamar.
* `iqrn(x)`           - Normalized inter-quartile range (normalized to equal σ for Gaussians).
* `shorthrange(x)`    - Length of the shortest closed interval containing at least half the data.
* `scaleQ(x)`         - Normalized Rousseeuw & Croux Q statistic, from the 25%ile of all 2-point distances.
* `scaleS(x)`         - Normalized Rousseeuw & Croux S statistic, from the median of the median of all 2-point distances.
* `shorthrange!(x)`, `scaleQ!(x)`, and `scaleS!(x)` are non-copying (that is, `x`-modifying) forms of the above.

### Confidence interval or standard error estimates:
* `trimse(x)` - Standard error of the trimmed mean.
* `trimci(x)` - Confidence interval for the trimmed mean.
* `msmedse(x)` - Standard error of the median.
* `binomci(s,n)` - Binomial confidence interval (Pratt's method).
* `acbinomci(s,n)` - Binomial confidence interval (Agresti-Coull method).
* `sint(x)`  - Confidence interval for the median (with optional p-value).
* `momci(x)` - Confidence interval of the modified one-step M-estimator of location (MOM).
* `trimpb(x)` - Confidence interval for trimmed mean.
* `pcorb(x)` - Confidence intervale for Pearson's correlation coefficient.
* `yuend` - Compare the trimmed means of two dependent random variables.
* `bootstrapci(x, est=f)` - Compute a confidence interval for estimator `f(x)` by bootstrap methods.
* `bootstrapse(x, est=f)` - Compute a standard error of estimator `f(x)` by bootstrap methods.

### Utility functions:
* `winval(x, tr=0.2)`         - Return a Winsorized copy of the data.
* `idealf(x)`                 - Ideal fourths, interpolated 1st and 3rd quartiles.
* `outbox(x)`                 - Outlier detection.
* `hpsi(x)`                   - Huber's ψ function.
* `contam_randn`              - Contaminated normal distribution (generates random deviates).
* `_weightedhighmedian(x)`    - Weighted median (breaks ties by rounding up). Used in scaleQ.

### Recommendations:
For location, consider the `bisquareWM` with k=3.9σ, if you can make any reasonable guess as to the "Gaussian-like width" σ (see dispersion estimators for this).  If not, `trimean` is a good second choice, though less efficient. Also, though the author personally has no experience with them, `tauloc`, `onestep`, and `mom` might be useful.

For dispersion, the `scaleS` is a good general choice, though `scaleQ` is very efficient for nearly Gaussian data.  The MAD is the most robust though less efficient.  If scaleS doesn't work, then shorthrange is a good second choice.

The first reference on scaleQ and scaleS (below) is a lengthy discussion of the tradeoffs among scaleQ, scaleS, shortest half, and median absolute deviation (MAD, see BaseStats.mad for Julia implementation). All four have the virtue of having the maximum possible breakdown point, 50%. This means that replacing up to 50% of the data with unbounded bad values leaves the statistic still bounded. The efficiency of Q is better than S and S is better than MAD (for Gaussian distributions), and the influence of a single bad point and the bias due to a fraction of bad points is only slightly larger on Q or S than on MAD. Unlike MAD, the other three do not implicitly assume a symmetric distribution.

To choose between Q and S, the authors note that Q has higher statistical efficiency, but S is typically twice as fast to compute and has lower gross-error sensitivity. An interesting advantage of Q over the others is that its influence function is continuous. For a rough idea about the efficiency, the large-N limit of the standardized variance of each quantity is 2.722 for MAD, 1.714 for S, and 1.216 for Q, relative to 1.000 for the standard deviation (given Gaussian data). The paper gives the ratios for Cauchy and exponential distributions, too; the efficiency advantages of Q are less for Cauchy than for the other distributions.


## Examples

    #Set up a sample dataset:
    x=[1.672064, 0.7876588, 0.317322, 0.9721646, 0.4004206, 1.665123, 3.059971, 0.09459603, 1.27424, 3.522148,
       0.8211308, 1.328767, 2.825956, 0.1102891, 0.06314285, 2.59152, 8.624108, 0.6516885, 5.770285, 0.5154299]

    julia> mean(x)     #the mean of this dataset
    1.853401259

#### `tmean`: trimmed mean

    julia> tmean(x)            #20% trimming by default
    1.2921802666666669

    julia> tmean(x, tr=0)      #no trimming; the same as the output of mean()
    1.853401259

    julia> tmean(x, tr=0.3)    #30% trimming
    1.1466045875000002

    julia> tmean(x, tr=0.5)    #50% trimming, which gives you the median of the dataset.
    1.1232023


#### `winval`: winsorize data

That is, return a copy of the input array, with the extreme low or high values
replaced by the lowest or highest non-extreme value, repectively. The fraction
considered extreme can be between 0 and 0.5, with 0.2 as the default.

    julia> winval(x)           #20% winsorization; can be changed via the named argument `tr`.
    20-element Any Array:
     1.67206
     0.787659
     0.400421
     0.972165
     ...
     0.651689
     2.82596
     0.51543


#### `winmean`, `winvar`, `wincov`: winsorized mean, variance, and covariance

    julia> winmean(x)          #20% winsorization; can be changed via the named argument `tr`.
    1.4205834800000001
    julia> winvar(x)
    0.998659015947531
    julia> wincov(x, x)
    0.998659015947531
    julia> wincov(x, x.^2)
    3.2819238397424004


#### `trimse`: estimated standard error of the trimmed mean

    julia> trimse(x)           #20% winsorization; can be changed via the named argument `tr`.
    0.3724280347984342


#### `trimci`: (1-α) confidence interval for the trimmed mean
Can be used for paired groups if `x` consists of the difference scores of two paired groups.

    julia> trimci(x)                 #20% winsorization; can be changed via the named argument `tr`.
    (1-α) confidence interval for the trimmed mean

    Degrees of freedom:   11
    Estimate:             1.292180
    Statistic:            3.469611
    Confidence interval:  0.472472       2.111889
    p value:              0.005244



#### `idealf`: the ideal fourths:

Returns `(q1,q3)`, the 1st and 3rd quartiles. These will be a weighted sum of
the values that bracket the exact quartiles, analogous to how we handle the
median of an even-length array.

    julia> idealf(x)
    (0.4483411416666667,2.7282743333333332)


#### `pbvar`: percentage bend midvariance

A robust estimator of scale (dispersion).
See [NIST ITL webpage](http://www.itl.nist.gov/div898/software/dataplot/refman2/auxillar/pbendmv.htm) for more.

    julia> pbvar(x)
    2.0009575278957623


#### `bivar`: biweight midvariance

A robust estimator of scale (dispersion).
See [NIST ITL webpage](http://www.itl.nist.gov/div898/software/dataplot/refman2/auxillar/biwmidv.htm) for more.

    julia> bivar(x)
    1.5885279811329132


#### `tauloc`, `tauvar`: tau measure of location and scale

Robust estimators of location and scale, with breakdown points of 50%.

See Yohai and Zamar _JASA_, vol 83 (1988), pp 406-413 and  Maronna and Zamar _Technometrics_, vol 44 (2002), pp. 307-317.

    julia> tauloc(x)       #the named argument `cval` is 4.5 by default.
    1.2696652567510853
    julia> tauvar(x)
    1.53008203090696


#### `outbox`: outlier detection

Use a modified boxplot rule based on the ideal fourths; when the named argument `mbox` is set to `true`, a modification of the boxplot rule suggested by Carling (2000) is used.

    julia> outbox(x)
    Outlier detection method using
    the ideal-fourths based boxplot rule

    Outlier ID:         17
    Outlier value:      8.62411
    Number of outliers: 1
    Non-outlier ID:     1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 18, 19, 20


#### `msmedse`: Standard error of the median
Return the standard error of the median, computed through the method recommended
by McKean and Schrader (1984).

    julia> msmedse(x)
    0.4708261134886094


#### `binomci()`, `acbinomci()`: Binomial confidence interval

Compute the (1-α) confidence interval for p, the binomial probability of success, given
`s` successes in `n` trials. Instead of `s` and `n`, can use a vector `x` whose values are all
0 and 1, recording failure/success one trial at a time. Returns an object.

`binomci` uses Pratt's method;
`acbinomci` uses a generalization of the Agresti-Coull method that was studied by Brown, Cai, & DasGupta.

    julia> binomci(2, 10)           # # of success and # of total trials are provided. By default alpha=.05
    p_hat:               0.2000
    confidence interval: 0.0274   0.5562
    Sample size          10


    julia> trials=[1, 0, 1, 1, 0, 0, 1, 0, 0, 1, 1, 0]
    julia> binomci(trials, alpha=0.01)    #trial results are provided in array form consisting of 1's and 0's.
     p_hat:               0.5000
     confidence interval: 0.1768   0.8495
     Sample size          12

    julia> acbinomci(2, 10)           # # of success and # of total trials are provided. By default alpha=.05
    p_hat:               0.2000
    confidence interval: 0.0459   0.5206
    Sample size          10



#### `sint()`
Compute the confidence interval for the median. Optionally, uses the
Hettmansperger-Sheather interpolation method to also estimate a p-value.

    julia> sint(x)
    Confidence interval for the median

     Confidence interval:  0.547483       2.375232

    julia> sint(x, 0.6)
    Confidence interval for the median with p-val

     Confidence interval:  0.547483       2.375232
     p value:              0.071861



#### `hpsi`
Compute Huber's ψ. The default bending constant is 1.28.

    julia> hpsi(x)
    20-element Array{Float64,1}:
    1.28
    0.787659
    0.317322
    0.972165
    0.400421
    ...


#### `onestep`
Compute one-step M-estimator of location using Huber's ψ. The default bending constant is 1.28.

    julia> onestep(x)
    1.3058109021286803



#### `bootstrapci`, `bootstrapse`
Compute a bootstrap, (1-α) confidence interval (`bootstrapci`) or a standard error (`bootstrapse`) for the measure of location corresponding to the argument `est`. By default, the median is used. Default α=0.05.

    julia> ci = bootstrapci(x, est=onestep, nullvalue=0.6)
     Estimate:             1.305811
     Confidence interval:  0.687723       2.259071
     p value:              0.026000


    julia> se = bootstrapse(x, est=onestep)
    0.41956761772722817




#### `mom` and `mom!`
Returns a modified one-step M-estimator of location (MOM), which is the unweighted
mean of all values not more than (bend times the `mad(x)`) away from the data
median.

    julia> mom(x)
    1.2596462322222222



#### `momci`
Compute the bootstrap (1-α) confidence interval for the MOM-estimator of location
based on Huber's ψ.  Default α=0.05.

    julia> momci(x, seed=2, nboot=2000, nullvalue=0.6)
    Estimate:             1.259646
    Confidence interval:  0.504223       2.120979
    p value:              0.131000




#### `contam_randn`
Create contaminated normal distributions. Most values will by from a N(0,1) zero-mean
unit-variance normal distribution. A fraction `epsilon` of all values will have `k`
times the standard devation of the others. Default: `epsilon=0.1` and `k=10`.


    julia> srand(1);
    julia> std(contam_randn(2000))
    3.516722458797104


#### `trimpb`
Compute a (1-α) confidence interval for a trimmed mean by bootstrap methods.

    julia> trimpb(x, nullvalue=0.75)
     Estimate:             1.292180
     Confidence interval:  0.690539       2.196381
     p value:              0.086000


#### `pcorb`
Compute a .95 confidence interval for Pearson's correlation coefficient. This function uses an adjusted percentile bootstrap method that gives good results when the error term is heteroscedastic.

    julia> pcorb(x, x.^5)
     Estimate:             0.802639
     Confidence interval:  0.683700       0.963478


#### `yuend`
Compare the trimmed means of two dependent random variables using the data in x and y. The default amount of trimming is 20%.

    julia> srand(3)
    julia> y2 = randn(20)+3;
    julia> yuend(x, y2)

    Comparing the trimmed means of two dependent variables.

    Sample size:          20
    Degrees of freedom:   11
    Estimate:            -1.547776
    Standard error:       0.460304
    Statistic:           -3.362507
    Confidence interval: -2.560898      -0.534653
    p value:              0.006336


### Unmaintained functions
See `UNMAINTAINED.md` for information about functions that the maintainers have not yet
understood but also not yet deleted entirely.


## References
* Percentage bend and related estimators come from L.H. Shoemaker and T.P. Hettmansperger ["Robust estimates and tests for the one- and two-sample scale models"](https://doi.org/10.1093/biomet/69.1.47) in _Biometrika_ Vol 69 (1982) pp. 47-53.

* Tau measures of location and scale are from V.J. Yohai and R.H. Zamar ["High Breakdown-Point Estimates of Regression by Means of the Minimization of an Efficient Scale"](http://doi/10.1080/01621459.1988.10478611) in _J. American Statistical Assoc._ vol 83 (1988) pp. 406-413.

* The `outbox(..., mbox=true)` modification was suggested in K. Carling, ["Resistant outlier rules and the non-Gaussian case"](http://dx.doi.org/10.1016/S0167-9473(99)00057-2) in _Computational Statistics and Data Analysis_ vol 33 (2000), pp. 249-258. doi:10.1016/S0167-9473(99)00057-2

* The estimate of the standard error of the median, `msmedse(x)`, is computed by the method of J.W. McKean and
R.M. Schrader, ["A comparison of methods for studentizing the sample median"](http://dx.doi.org/10.1080/03610918408812413) in _Communications in Statistics: Simulation and Computation_ vol 13 (1984) pp. 751-773. doi:10.1080/03610918408812413

* For Pratt's method of computing binomial confidence intervals, see J.W. Pratt (1968) ["A normal approximation for binomial, F, Beta, and other common, related tail probabilities, II"](http://dx.doi.org/10.1080/01621459.1968.10480939) _J. American Statistical Assoc._, vol 63, pp. 1457- 1483, doi:10.1080/01621459.1968.10480939.  Also R.G. Newcombe ["Confidence Intervals for a binomial proportion"](http://dx.doi.org/10.1002/sim.4780131209) _Stat. in Medicine_ vol 13 (1994) pp 1283-1285, doi:10.1002/sim.4780131209.

* For the Agresti-Coull method of computing binomial confidence intervals, see L.D. Brown, T.T. Cai, & A. DasGupta ["Confidence Intervals for a Binomial Proportion and Asymptotic Expansions"](http://www.jstor.org/stable/2700007) in _Annals of Statistics_, vol 30 (2002), pp. 160-201.

* Shortest Half-range comes from P.J. Rousseeuw and A.M. Leroy, ["A Robust Scale Estimator Based on the Shortest Half"](http://onlinelibrary.wiley.com/doi/10.1111/j.1467-9574.1988.tb01224.x/abstract) in _Statistica Neerlandica_ Vol 42 (1988), pp. 103-116. doi:10.1111/j.1467-9574.1988.tb01224.x . See also R.D. Martin and R. H. Zamar, ["Bias-Robust Estimation of Scale"](http://projecteuclid.org/euclid.aos/1176349161)  in _Annals of Statistics_ Vol 21 (1993) pp. 991-1017.  doi:10.1214/aoe/1176349161

* Scale-Q and Scale-S statistics are described in P.J. Rousseeuw and C. Croux ["Alternatives to the Median Absolute Deviation"](http://www.jstor.org/stable/2291267) in _J. American Statistical Assoc._ Vo 88 (1993) pp 1273-1283. The time-efficient algorithms for computing them appear in C. Croux and P.J. Rousseeuw, ["Time-Efficient Algorithms for Two Highly Robust Estimators of Scale"](http://link.springer.com/chapter/10.1007/978-3-662-26811-7_58) in _Computational Statistics, Vol I_ (1992), Y. Dodge and J. Whittaker editors, Heidelberg, Physica-Verlag, pp 411-428. If link fails, see ftp://ftp.win.ua.ac.be/pub/preprints/92/Timeff92.pdf
