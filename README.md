RobustStats
==============


[![Build Status](https://travis-ci.org/maximsch2/RobustStats.jl.svg?branch=master)](https://travis-ci.org/maximsch2/RobustStats.jl)

Most functions in this file are robust statistical methods based on the R package WRS ([an R-Forge repository](https://r-forge.r-project.org/projects/wrs/)) by [Rand Wilcox](http://dornsife.usc.edu/cf/labs/wilcox/wilcox-faculty-display.cfm). Only a handful of functions are included at this point. Others were contributed by users as needed. [References](#References) can be found below.

This package requires `Compat`, `Rmath`, `Dataframes`, and `Distributions`. They can be installed automatically, or by invoking `Pkg.add("packagename")`.


##Examples

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


#### `trimci`: (1-alpha) confidence interval for the trimmed mean
Can be used for paired groups if `x` consists of the difference scores of two paired groups.

    julia> trimci(x)                 #20% winsorization; can be changed via the named argument `tr`.
    1-alpha confidence interval for the trimmed mean

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



####17. `akerd`:
Compute adaptive kernel density estimate for univariate data (See Silverman, 1986)

    julia> akerd(x, title="Lognormal Distribution", xlab="x", ylab="Density")

![plot](http://img824.imageshack.us/img824/4894/sme4.png)

    julia> srand(10)
    julia> x3=rnorm(100, 1, 2);
    julia> akerd(x3, title="Normal Distribution; mu=1, sd=2", xlab="x", ylab="Density", color="red", plottype="dash")

![plot](http://img547.imageshack.us/img547/4920/zb8j.png)


####18. `indirectTest`:
This function is adapted from Andrew Hayes' SPSS macro, which evaluates indirect/mediation effects via percentile bootstrap and the Sobel Test.

    julia> srand(1)
    julia> m = randn(20);        #Mediator
    julia> srand(2)
    julia> y = randn(20) + 2.0;  #Outcome variable

    julia> indirectTest(y, x, m)   #5000 bootstrap samples by default
    TESTS OF INDIRECT EFFECT

    Sample size: 20
    Number of bootstrap samples: 5000

    DIRECT AND TOTAL EFFECTS
              Estimate  Std.Error     t value   Pr(>|t|)
    b(YX):    0.125248   0.087729    1.427664   0.170508
    b(MX):    0.082156   0.106198    0.773611   0.449202
    b(YM.X):  0.140675   0.089359    1.574264   0.133852
    b(YX.M): -0.187775   0.195111   -0.962402   0.349338

    INDIRECT EFFECT AND SIGNIFICANCE USING NORMAL DISTRIBUTION
              Estimate  Std.Error     z score      CI lo      CI hi  Pr(>|z|)
    Sobel:   -0.015427   0.019161   -0.805126  -0.052981   0.022128  0.420747

    BOOTSTRAP RESULTS OF INDIRECT EFFECT
              Estimate  Std.Error       CI lo      CI hi    P value
    Effect:  -0.048598   0.084280   -0.298647   0.019394   0.420800

if we add `plotit=true`, we get a kernel density plot of the effects derived from bootstrap samples:

![plot](http://img405.imageshack.us/img405/7603/wdiu.png)




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



#### `bootstrapci()`
Compute a bootstrap, (1-α) confidence interval for the measure of location corresponding to the argument `est`. By default, a one-step M-estimator of location based on Huber's ψ is used. The default number of bootstrap samples is `nboot=2000`. `nv` is the target value used when computing a p-value. Default α=0.05.

    julia> bootstrapci(x)
    Estimate:             1.305811
    Confidence interval:  0.687723       2.259071
    p value:             < 10e-16




#### `mom` and `mom!`
Returns a modified one-step M-estimator of location (MOM), which is the unweighted
mean of all values not more than (bend times the `mad(x)`) away from the data
median.

    julia> mom(x)
    1.2596462322222222



#### `momci`
Compute the bootstrap (1-α) confidence interval for the MOM-estimator of location
based on Huber's ψ.  Default α=0.05.

    julia> momci(x, seed=2, nboot=2000, nv=0.6)
    Estimate:             1.259646
    Confidence interval:  0.504223       2.120979
    p value:              0.131000




#### `contam_randn()`
Create contaminated normal distributions. Most values will by from a N(0,1) zero-mean
unit-variance normal distribution. A fraction `epsilon` of all values will have `k`
times the standard devation of the others. Default: `epsilon=0.1` and `k=10`.


    julia> srand(1);
    julia> std(contam_randn(2000))
    3.516722458797104


####29. `trimpb()`
Compute a 1-alpha confidence interval for a trimmed mean. The default number of bootstrap samples is nboot=2000. win is the amount of Winsorizing before bootstrapping when WIN=true.

plotit=true gives a plot of the bootstrap values
* pop=1 results in the expected frequency curve.
* pop=2 kernel density estimate    NOT IMPLEMENTED
* pop=3 boxplot                    NOT IMPLEMENTED
* pop=4 stem-and-leaf              NOT IMPLEMENTED
* pop=5 histogram
* pop=6 adaptive kernel density estimate.

fr controls the amount of smoothing when plotting the bootstrap values via the function rdplot. fr=NaN means the function will use fr=.8 (When plotting bivariate data, rdplot uses fr=.6 by default.)


    julia> trimpb(x, plotit=true, pop=6)
    Compute a 1-alpha confidence interval for a trimmed mean using the bootstrap percentile method.

     Confidence interval:  0.708236       2.253489
     p value:             < 10e-16

![plot](http://img585.imageshack.us/img585/9997/iua.png)




####30. `trimcibt()`
Compute a 1-alpha confidence interval for the trimmed mean using a bootstrap percentile t method. The default amount of trimming is tr=.2. side=true,  indicates the symmetric two-sided method. side=false yields an equal-tailed confidence interval. NOTE: p.value is reported when side=true only.

    julia> trimcibt(x, nboot=5000, plotit=true)
    Bootstrap .95 confidence interval for the trimmed mean
    using a bootstrap percentile t method

     Estimate:             1.292180
     Statistic:            3.469611
     Confidence interval:  0.292162       2.292199
     p value:              0.022600

![plot](http://imageshack.us/a/img844/1970/34o.png)




####31. `pcorb()`
Compute a .95 confidence interval for Pearson's correlation coefficient. This function uses an adjusted percentile bootstrap method that gives good results when the error term is heteroscedastic.

    julia> pcorb(x, y)
    Estimate:             0.318931
    Confidence interval: -0.106467       0.663678


####32. `yuend()`
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


####33. `t1way()`
A heteroscedastic one-way ANOVA for trimmed means using a generalization of Welch's method. When `tr=0`, the function conducts a heteroscedastic 1-way ANOVA without trimming.

There are two ways to specify data for the function:
* A two dimensional array where each column represents a group. Hence, a 10 X 3 array means that a one way ANOVA will be performed on 3 groups. This also means equal sample sizes amongst the groups.
* A vector containing the outcome variable and another vector containing group information.

Two examples are shown below to demonstrate the two ways of specifying data:

    #Data in two dimensional array form
    #Prepare data
    julia> srand(12)
    julia> m2 = reshape(sort(randn(30)), 10, 3);

    julia> t1way(m2)
    Heteroscedastic one-way ANOVA for trimmed means
    using a generalization of Welch's method.

    Sample size:          10   10   10
    Degrees of freedom:   2.00   8.25
    Statistic:            20.955146
    p value:              0.000583


    #Data in vector form
    julia> srand(12)
    julia> m3 = sort(randn(30));
    julia> group = rep(1:3, [8,12,10]);   #Unequal sample sizes: n1 = 8, n2 = 12, n3 = 10

    julia> t1way(m3, group)
    Heteroscedastic one-way ANOVA for trimmed means
    using a generalization of Welch's method.

    Sample size:          8   12   10
    Degrees of freedom:   2.00   8.11
    Statistic:            28.990510
    p value:              0.000202

### Unmaintained functions
The functions in this section are undocumented. If you know what any of them mean and can help us to document, fix, and/or improve the code, please contact the maintainers.


#### `stein1`, `stein2`: Two steps of Stein's method

    julia> stein1(x, 1)
    41
    julia> srand(1)         #set seed
    julia> x2=rnorm(21);    #suppose additional 21 data points were collected.
    julia> stein2(x, x2)
    df:        19
    estimate:  0.8441066102423266
    confint:   [-3.65043, 5.33865]
    statistic: 2.516970580783766
    crit:      2.093024054408309
    pval:      0.020975586092444765

#### `stein1_tr()`, `stein2_tr()`
Extension of Stein's method based on the trimmed mean.

    julia> stein1_tr(x, 0.2)
    89

    julia> stein2_tr(x, x2)
    Extension of the 2nd stage of Stein's method based on the trimmed mean

    Statistic:            3.154993
    Critical value:       2.200985


This function is able to handle multiple dependent groups. Suppose that the original dataset contained 4 dependent groups and the sample size is 5 (`xold`), and we collected more data (`xnew`).

    julia> srand(2)
    julia> xnew = rand(6, 4)
    julia> xold = reshape(x, 5, 4)
    julia> stein2_tr(xold, xnew)

    Extension of the 2nd stage of Stein's method based on the trimmed mean

     Statistic             Group 1  Group 2  Statistic
                                 1        2  -0.933245
                                 1        3   0.291014
                                 1        4  -0.949618
                                 2        3   1.212577
                                 2        4  -0.608510
                                 3        4  -0.649426
     Critical value:       10.885867




## Estimators from other sources
###Location estimators:
* `bisquareWM(x)`     - Mean with weights given by the bisquare rho function.
* `huberWM(x)`        - Mean with weights given by Huber's rho function.
* `trimean(x)`        - Tukey's trimean, the average of the median and the midhinge.

### Dispersion estimators:
* `shorthrange(x)`    - Length of the shortest closed interval containing at least half the data.
* `scaleQ(x)`         - Normalized Rousseeuw & Croux Q statistic, from the 25%ile of all 2-point distances.
* `scaleS(x)`         - Normalized Rousseeuw & Croux S statistic, from the median of the median of all 2-point distances.
* `shorthrange!(x)`, `scaleQ!(x)`, and `scaleS!(x)` are non-copying (that is, `x`-modifying) forms of the above.

### Utility functions:
* `_weightedhighmedian(x)`    - Weighted median (breaks ties by rounding up). Used in scaleQ.

### Recommendations:
For location, consider the `bisquareWM` with k=3.9*sigma, if you can make any reasonable guess as to the "Gaussian-like width" sigma (see dispersion estimators for this).  If not, trimean is a good second choice, though less efficient.

For dispersion, the `scaleS` is a good general choice, though `scaleQ` is very efficient for nearly Gaussian data.  The MAD is the most robust though less efficient.  If scaleS doesn't work, then shorthrange is a good second choice.

The first reference on scaleQ and scaleS (below) is a lengthy discussion of the tradeoffs among scaleQ, scaleS, shortest half, and median absolute deviation (MAD, see BaseStats.mad for Julia implementation). All four have the virtue of having the maximum possible breakdown point, 50%. This means that replacing up to 50% of the data with unbounded bad values leaves the statistic still bounded. The efficiency of Q is better than S and S is better than MAD (for Gaussian distributions), and the influence of a single bad point and the bias due to a fraction of bad points is only slightly larger on Q or S than on MAD. Unlike MAD, the other three do not implicitly assume a symmetric distribution.

To choose between Q and S, the authors note that Q has higher statistical efficiency, but S is typically twice as fast to compute and has lower gross-error sensitivity. An interesting advantage of Q over the others is that its influence function is continuous. For a rough idea about the efficiency, the large-N limit of the standardized variance of each quantity is 2.722 for MAD, 1.714 for S, and 1.216 for Q, relative to 1.000 for the standard deviation (given Gaussian data). The paper gives the ratios for Cauchy and exponential distributions, too; the efficiency advantages of Q are less for Cauchy than for the other distributions.


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
