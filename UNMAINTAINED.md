## Unmaintained Functions

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

     #### `akerd`:
     Compute adaptive kernel density estimate for univariate data (See Silverman, 1986)

         julia> akerd(x, title="Lognormal Distribution", xlab="x", ylab="Density")


         julia> srand(10)
         julia> x3=rnorm(100, 1, 2);
         julia> akerd(x3, title="Normal Distribution; mu=1, sd=2", xlab="x", ylab="Density", color="red", plottype="dash")



     #### `indirectTest`:
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
