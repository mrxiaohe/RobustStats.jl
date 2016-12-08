RobustStats
==============


Most functions in this file are robust statistical methods based on the R package WRS ([R-Forge repository](https://r-forge.r-project.org/projects/wrs/)) by [Rand Wilcox](http://dornsife.usc.edu/cf/labs/wilcox/wilcox-faculty-display.cfm). Only a handful of functions are included at this point. More will be added soon.


In order to use the functions, `Rmath`, `Stats`, `Dataframes`, `Distributions`, and `Winston` are required. They can be installed by invoking `Pkg.add("packagename")`.


##Examples

    #Set up a sample dataset:
    x=[1.672064, 0.7876588, 0.317322, 0.9721646, 0.4004206, 1.665123, 3.059971, 0.09459603, 1.27424, 3.522148, 
       0.8211308, 1.328767, 2.825956, 0.1102891, 0.06314285, 2.59152, 8.624108, 0.6516885, 5.770285, 0.5154299]

    julia> mean(x)     #the mean of this dataset
    1.853401259
    
####1. `tmean`: trimmed mean
    
    julia> tmean(x)            #20% trimming by default
    1.2921802666666669
    
    julia> tmean(x, tr=0)      #no trimming; the same as the output of mean()
    1.853401259
    
    julia> tmean(x, tr=0.3)    #30% trimming
    1.1466045875000002
    
    julia> tmean(x, tr=0.5)    #50% trimming, which gives you the median of the dataset.
    1.1232023
    
    
####2. `winval`: winsorize data
    
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


####3. `winmean`: winsorized mean

    julia> winmean(x)          #20% winsorization; can be changed via the named argument `tr`.
    1.4205834800000001


####4. `winvar`: winsorized variance

    julia> winvar(x)           #20% winsorization; can be changed via the named argument `tr`.
    0.998659015947531


####5. `trimse`: estimated standard error of the gamma trimmed mean

    julia> trimse(x)           #20% winsorization; can be changed via the named argument `tr`.       
    0.3724280347984342


####6. `trimci`: 1-alpha confidence interval for the trimmed mean
Can be used for paired groups if `x` consists of the difference scores of two paired groups.
  
    julia> trimci(x)                 #20% winsorization; can be changed via the named argument `tr`.       
    1-alpha confidence interval for the trimmed mean

    Degrees of freedom:   11
    Estimate:             1.292180
    Statistic:            3.469611
    Confidence interval:  0.472472       2.111889
    p value:              0.005244



####7. `stein1`: Stein's method

    julia> stein1(x, 1)     
    41


####8. `stein2`: the second stage of Stein's method.
 
    julia> srand(1)         #set seed
    julia> x2=rnorm(21);    #suppose additional 21 data points were collected.
    julia> stein2(x, x2)
    df:        19
    estimate:  0.8441066102423266
    confint:   [-3.65043, 5.33865]
    statistic: 2.516970580783766
    crit:      2.093024054408309
    pval:      0.020975586092444765


####9. `idealf`: the ideal fourths:

    julia> idealf(x)
    The ideal fourths
    Lower quartile:  0.448341
    Upper quartile:  2.728274


####10. `pbvar`: percentage bend midvariance

    julia> pbvar(x)
    2.0009575278957623
    
    
####11. `bivar`: biweight midvariance

    julia> bivar(x)
    1.588527039652727
    
    
####12. `tauloc`: tau measure of location 

(see Yohai and Zamar (JASA, 83, 406-413))

    julia> tauloc(x)       #the named argument `cval` is 4.5 by default.
    1.2696674487629664


####13. `tauvar`: the tau measure of scale 

(see Yohai and Zamar (JASA, 1988, 83, 406-413) and  Maronna and Zamar (Technometrics, 2002, 44, 307-317))

    julia> tauvar(x)
    1.5300804969271078


####14. `outbox`: outlier detection 

using a modified boxplot rule based on the ideal fourths; when the named argument `mbox` is set to `true`, a modification of the boxplot rule suggested by Carling (2000) is used

    julia> outbox(x)     
    Outlier detection method using 
    the ideal-fourths based boxplot rule

    Outlier ID:         17                                                        
    Outlier value:      8.62411
    Number of outliers: 1
    Non-outlier ID:     1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 18, 19, 20
    

####15. `msmedse`:
Computing standard error of the median using a method recommended by McKean and Shrader (1984)

    julia> msmedse(x)
    0.4708261134886094


####16. `binomci`: 
Compute a 1-alpha confidence interval for p, the probability of success for a binomial dist. using Pratt's method y is a vector of 1's and 0's, x is the number of successes observed among n trials.

    julia> binomci(2, 10)           #when number of success and number of total trials are provided. By default alpha=.05
    phat:                0.2000
    confidence interval: 0.0274   0.5562
    Sample size          10

    
    julia> trials=[1, 0, 1, 1, 0, 0, 1, 0, 0, 1, 1, 0]
    julia> binomci(trials, 0.01)    #trial results are provided in array form consisting of 1's and 0's. Alpha=.01
    phat:    0.5
    confint: [0.176811, 0.849496]
    n:       12
    

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




####18. `sint()`
Computing the confidence interval for the median.

    julia> sint(x)
    Confidence interval for the median

     Confidence interval:  0.547483       2.375232
     


####19. `acbinomci()`
Compute a 1-alpha confidence interval for p, the probability of success for a binomial dist using a generalization of the Agresti-Coull  method that was studied by Brown, Cai DasGupta.

    julia> acbinomci(2, 10)           #when number of success and number of total trials are provided. By default alpha=.05
    phat:                0.2000
    confidence interval: 0.0459   0.5206
    Sample size          10



####21. `stein1_tr()`
Extension of Stein's method based on the trimmed mean.

    julia> stein1_tr(x, 0.2)
    89



####22. `stein2_tr()`
Extension of the 2nd stage of Stein's method based on the trimmed mean.

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
     

####23. `sintv2()`
Confidence interval for the median using the Hettmansperger-Sheather interpolation method.

    julia> sintv2(x)
    Confidence interval for the median
    using the Hettmansperger-Sheather interpolation method

     Confidence interval:  0.547483       2.375232
     p value:              0.000100



####24. `onestep()`
Compute one-step M-estimator of location using Huber's Psi. The default bending constant is 1.28.

    julia> onestep(x)
    1.384070801414857
    



####25. `onesampb()`
Compute a bootstrap, .95 confidence interval for the measure of location corresponding to the argument est. By default, a one-step M-estimator of location based on Huber's Psi is used. The default number of bootstrap samples is nboot=2000. nv=0 when computing a p-value.

    julia> onesampb(x)
    Estimate:             1.384071
    Confidence interval:  0.760978       2.378872
    p value:             < 10e-16
    



####26. `mom()`
Compute MOM-estimator of location. The default bending constant is 2.24

    julia> mom(x)
    1.2596462322222222
    


####27. `momci()`
Compute the bootstrap .95 confidence interval for the MOM-estimator of location based on Huber's Psi.

    julia> momci(x, nboot=2)
    Bootstrap .95 confidence interval for the MOM-
    estimator of location based on Huber's Psi

     Confidence interval:  0.652212       1.684510
     
     

####28. `cnorm()`
Create contaminated normal distributions.

    julia> srand(1)
    julia> c = cnorm(2000)
    julia> akerd(c)

![plot](http://img818.imageshack.us/img818/7345/8uq8.png)




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


####32. `t1way()`
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
    



