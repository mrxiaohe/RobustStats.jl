"""`tmean(x; tr=0.2)`

Trimmed mean of real-valued array `x`.

Find the mean of `x`, omitting the lowest and highest `tr` fraction of the data.
This requires `0 <= tr <= 0.5`. The amount of trimming defaults to `tr=0.2`.
"""
function tmean{S <: Real}(x::AbstractArray{S}; tr::Real=0.2)
    tmean!(copy(x), tr=tr)
end


"""`tmean!(x; tr=0.2)`

Trimmed mean of real-valued array `x`, which sorts the vector `x` in place.

Find the mean of `x`, omitting the lowest and highest `tr` fraction of the data.
This requires `0 <= tr <= 0.5`. The trimming fraction defaults to `tr=0.2`.
"""
function tmean!{S <: Real}(x::AbstractArray{S}; tr::Real=0.2)
    if tr < 0 || tr > 0.5
        error("tr cannot be smaller than 0 or larger than 0.5")
    elseif tr == 0
        return mean(x)
    elseif tr == .5
        return median!(x)
    else
        n   = length(x)
        lo  = floor(Int64, n*tr)+1
        hi  = n+1-lo
        return mean(sort!(x)[lo:hi])
    end
end


"""`winval(x; tr=0.2)`

Winsorize real-valued array `x`.

Return a copy of `x` in which extreme values (that is, the lowest and highest
fraction `tr` of the data) are replaced by the lowest or highest non-extreme
value, as appropriate. The trimming fraction defaults to `tr=0.2`.
"""
function winval{S <: Real}(x::AbstractArray{S}; tr::Real=0.2)
    const n = length(x)
    xcopy   = sort(x)
    ibot    = floor(Int64, tr*n)+1
    itop    = n-ibot+1
    xbot, xtop = xcopy[ibot], xcopy[itop]
    return  [x[i]<=xbot ? xbot : (x[i]>=xtop ? xtop : x[i]) for i=1:n]
end

"""`winmean(x; tr=0.2)`

Winsorized mean of real-valued array `x`.

See `winval` for what Winsorizing (clipping) signifies.
"""
winmean{S <: Real}(x::AbstractArray{S}; tr=0.2) = mean(winval(x, tr=tr))

"""`winvar(x; tr=0.2)`

Winsorized variance of real-valued array `x`.

See `winval` for what Winsorizing (clipping) signifies.
"""
winvar{S <: Real}(x::AbstractArray{S}; tr=0.2) = var(winval(x, tr=tr))

"""`winstd(x; tr=0.2)`

Winsorized standard deviation of real-valued array `x`.

See `winval` for what Winsorizing (clipping) signifies.
"""
winstd{S <: Real}(x::AbstractArray{S}; tr=0.2) = std(winval(x, tr=tr))


"""`wincov(x, y; tr=0.2)`

Compute the Winsorized covariance between `x` and `y`.

See `winval` for what Winsorizing (clipping) signifies.
"""
function wincov{S <: Real, T <: Real}(x::AbstractArray{S}, y::AbstractArray{T}; tr::Real=0.2)
    xvec = winval(x, tr=tr)
    yvec = winval(y, tr=tr)
    wcov = cov(xvec, yvec)
end



"""`trimse(x; tr=0.2)`

Estimated standard error of the mean for Winsorized real-valued array `x`.

See `winval` for what Winsorizing (clipping) signifies.
"""
trimse{S <: Real}(x::AbstractArray{S}; tr::Real=0.2) =
    sqrt(winvar(x,tr=tr))/((1-2tr)*sqrt(length(x)))

"""`trimci(x; tr=0.2, alpha=0.05, ...)`

Compute a (1-α) confidence interval for the trimmed mean.

Returns a `RobustStats.testOutput` object.
"""
function trimci{S <: Real}(x::AbstractArray{S}; tr::Real=0.2, alpha::Real=0.05, nullvalue::Real=0, method=true)
    se  = trimse(x, tr=tr)
    n   = length(x)
    df::Int64   = n-2*floor(tr*n)-1
    estimate    = tmean(x, tr=tr)
    confint     = [estimate-Rmath.qt(1-alpha/2, df)*se,
                   estimate+Rmath.qt(1-alpha/2, df)*se]
    statistic   = (estimate-nullvalue)/se
    pval        = 2*(1-Rmath.pt(abs(statistic),df))
    METHOD      = method ? "(1-α) confidence interval for the trimmed mean\n": nothing
    output           = testOutput()
    output.method    = METHOD
    output.df        = df
    output.estimate  = estimate
    output.ci        = confint
    output.statistic = statistic
    output.p         = pval
    return output
end


"""`idealf(x)`

Compute the ideal fourths (interpolated quartiles) of real-valued array `x`.

Returns a tuple of (1st_quartile, 3rd_quartile)
"""
function idealf{S <: Real}(x::AbstractArray{S})
    y       = sort(x)
    n       = length(x)
    j       = floor(Int64, n/4+5/12) # 25%ile is in [y[j], y[j+1]]
    k       = n-j+1        # 75%ile is in [y[k],y[k-1]]
    g       = n/4+5/12 - j   # weighting for the two data surrounding quartiles.
    (1-g).*y[j]+g.*y[j+1], (1-g).*y[k]+g.*y[k-1]
end

"""`pbvar(x; beta=0.2)`

Return the percentage bend midvariance of real-valued array `x`, a robust, efficient
measure of scale (dispersion). Lower values of beta increase efficiency but reduce
robustness.
This requires `0 <= beta <= 0.5`. The trimming fraction defaults to `beta=0.2`.
"""
function pbvar{S <: Real}(x::AbstractArray{S}; beta::Real=0.2)
    const n = length(x)
    med = median(x)
    absdev = abs.(x-med)
    sort!(absdev)

    m = floor(Int64, (1-beta)*n+0.5)
    ω = absdev[m]
    if ω <= 0   # At least a fraction (1-beta) of all values are identical
        return 0.0
    end

    z = 0.0
    counter = 0
    for i = 1:n
        ψ = absdev[i]/ω
        if abs(ψ) >= 1.0
            z += 1.0
        else
            z += ψ^2
            counter += 1
        end
    end
    n*(ω^2)*z/(counter^2)
end


"""`bivar(x; beta=0.2)`

Return the biweight midvariance of real-valued array `x`, a robust, efficient
measure of scale (dispersion). Lower values of beta increase efficiency but reduce
robustness.
This requires `0 <= beta <= 0.5`. The trimming fraction defaults to `beta=0.2`.
"""
function bivar{S <: Real}(x::AbstractArray{S})
    const n = length(x)
    med = median(x)
    MAD = mad(x)
    q = Rmath.qnorm(0.75)
    top = bot = 0.0
    for i = 1:n
        u = abs(x[i]-med)./(9.*q.*MAD)
        if u<1.0
            top += n*(x[i]-med)*(x[i]-med)*(1-u*u).^4
            bot += (1-u*u)*(1-5*u*u)
        end
    end
    top/(bot^2)
end


"""`tauloc(x; cval=4.5)`

Return the tau measure of location of real-valued array `x`, a robust, efficient
estimator.
"""
function tauloc{S <: Real}(x::AbstractArray{S}; cval::Real=4.5)
    const n = length(x)
    med = median(x)
    s = Rmath.qnorm(0.75)*mad(x)
    Wnom = Wden = 0.0
    for i in 1:n
        y = (x[i]-med)/s
        temp = (1.0-(y/cval)^2)^2
        if abs(temp) <= cval
            Wnom += temp*x[i]
            Wden += temp
        end
    end
    Wnom/Wden
end


"""`tauvar(x; cval=3.0)`

Return the tau measure of dispersion of real-valued array `x`, a robust, efficient
estimator.
"""
function tauvar{S <: Real}(x::AbstractArray{S}; cval::Real=3.0)
    const n = length(x)
    s     = Rmath.qnorm(0.75)*mad(x)
    tloc  = tauloc(x)
    W     = 0.0
    cval2 = cval*cval
    [W    += min(((x[i]-tloc)/s)*((x[i]-tloc)/s), cval2) for i=1:n]
    s*s*W/n
end


"""`outbox(x; mbox::Bool=false, ...)`

Use a modified boxplot rule based on the ideal fourths (`idealf`). When the named argument
`mbox` is set to true, a modification of the boxplot rule suggested by Carling (2000) is used.

Returns an object with vectors `keepid` and `outid` giving the kept/rejected element numbers,
`nout` (the number of rejected elements), and `outval`, an array of the outlier values.
"""
function outbox{S <: Real}(x::AbstractArray{S}; mbox::Bool=false, gval::Real=NaN, method::Bool=true)
    const n = length(x)
    lower_quartile, upper_quartile = idealf(x)
    IQR = upper_quartile-lower_quartile
    cl = cu = 0.0
    if mbox
        if isnan(gval)
            gval=(17.63*n-23.64)/(7.74*n-3.71)
        end
        cl = median(x) - gval*IQR
        cu = median(x) + gval*IQR
    elseif !mbox
        if isnan(gval)
            gval=1.5
        end
        cl = lower_quartile - gval*IQR
        cu = upper_quartile + gval*IQR
    end
    flag = (x.<cl) .| (x.>cu)
    vec = 1:n
    outid  = vec[flag]
    keepid = vec[!flag]
    outval = x[flag]
    nout = length(outid)
    if method && !mbox
        METHOD = "Outlier detection method using \nthe ideal-fourths based boxplot rule\n"
    elseif method && mbox
        METHOD = "Outlier detection method using \nthe ideal-fourths based boxplot rule\n(using the modification suggested by Carling (2000))\n"
    else
        METHOD = nothing
    end
    outOutput(outid, keepid, outval, nout, METHOD)
end


"""`msmedse(x)`

Return the standard error of the median, computed through the method recommended
by McKean and Sshrader (1984)."""
function msmedse{S <: Real}(x::AbstractArray{S})
    const n = length(x)
    y = sort(x)
    if duplicated(y)
        warn("Tied values detected. Estimate of standard error might be highly inaccurate, even with n large")
    end
    q995 = Rmath.qnorm(.995)
    av::Int = round((n+1)/2 - q995*sqrt(n/4))
    if av == 0
        av = 1
    end
    top::Int = n-av+1
    abs((y[top]-y[av])/(2q995))
end


"""`binomci(s, n; alpha=0.05)`

Compute the (1-α) confidence interval for p, the binomial probability of success, given
`s` successes in `n` trials. Returns an object with components `p_hat` (the observed
fraction of successes) and `confint=[lo,hi]` (the confidence interval). The computation
uses Pratt's method.

Can also use `binomci(x; alpha=0.05)`, where x is an array consisting only of 0s
and 1s. It's equivalent to `binomci(sum(x), length(x), alpha=alpha)`."""
function binomci(s::Int, n::Int; alpha::Real=0.05)
    if s > n
        error("binomci requires s≤n (no more successes than trials)")
    elseif s < 0
        error("binomci requires s≥0")
    elseif n <= 1
        error("binomci requires n≥2 (at least 2 trials)")
    end
    p_hat=s/n
    if s == 0
        upper = 1.0-alpha.^(1/n)
        return binomciOutput(p_hat, [0,upper], n)
    elseif s == 1
        lower = 1-(1-alpha/2).^(1/n)
        upper = 1-(alpha/2).^(1/n)
        return binomciOutput(p_hat, [lower, upper], n)
    elseif s == (n-1)
        lower = (alpha/2).^(1/n)
        upper = (1-alpha/2).^(1/n)
        return binomciOutput(p_hat, [lower, upper], n)
    elseif s == n
        lower = alpha.^(1/n)
        upper = 1
        return binomciOutput(p_hat, [lower, upper], n)
    end

    z     = Rmath.qnorm(1-alpha/2)
    A     = ((s+1)/(n-s))*((s+1)/(n-s))
    B     = 81.*(s+1)*(n-s)-9.*n-8
    C     = (0-3)*z*sqrt(9.*(s+1)*(n-s)*(9*n+5-z^2)+n+1)
    D     = 81.*(s+1)^2-9.*(s+1)*(2+z^2)+1
    E     = 1+A*((B+C)/D)^3
    upper = 1/E

    A     = (s/(n-s-1))*(s/(n-s-1))
    B     = 81.*s*(n-s-1)-9.*n-8
    C     = 3.*z*sqrt(9.*s*(n-s-1)*(9.*n+5-z^2)+n+1)
    D     = 81.*s^2-9.*s*(2+z^2)+1
    E     = 1+A*((B+C)/D)^3
    lower = 1/E
    binomciOutput(p_hat, [lower, upper], n)
end


function binomci(x::Vector{Int}; alpha::Real=0.05)
    for i = 1:length(x)
        if x[i]<0 || x[i] > 1
            error("x vector must contain only values 0 or 1.")
        end
    end
    binomci(sum(x), length(x), alpha=alpha)
end



"""`acbinomci(s, n; alpha=0.05)`

Compute the (1-α) confidence interval for p, the binomial probability of success, given
`s` successes in `n` trials. Returns an object with components `p_hat` (the observed
fraction of successes) and `confint=[lo,hi]` (the confidence interval). The computation
uses a generalization of the Agresti-Coull  method that was studied by Brown, Cai, & DasGupta.

Can also use `acbinomci(x; alpha=0.05)`, where `x` is an array consisting only of 0s
and 1s. It's equivalent to `acbinomci(sum(x), length(x), alpha=alpha)`."""
function acbinomci(s::Int, n::Int; alpha::Real=0.05)
    if s > n
        error("acbinomci requires s≤n (no more successes than trials)")
    elseif s < 0
        error("acbinomci requires s≥0")
    elseif n <= 1
        error("acbinomci requires n≥2 (at least 2 trials)")
    end
    p_hat=s/n

    if s == 0
        upper = 1.0-alpha.^(1/n)
        return binomciOutput(p_hat, [0, upper], n)
    elseif s == 1
        lower = 1-(1-alpha/2)^(1/n)
        upper = 1-(alpha/2)^(1/n)
        return binomciOutput(p_hat, [lower, upper], n)
    elseif s == (n-1)
        lower = (alpha/2)^(1/n)
        upper = (1-alpha/2)^(1/n)
        return binomciOutput(p_hat, [lower, upper], n)
    elseif s == n
        lower = alpha^(1/n)
        upper = 1
        return binomciOutput(p_hat, [lower, upper], n)
    end

    cr    = Rmath.qnorm(1-alpha/2)
    ntil  = n+cr^2
    ptil  = (s+cr^2/2)/ntil
    lower = ptil-cr*sqrt(ptil*(1-ptil)/ntil)
    upper = ptil+cr*sqrt(ptil*(1-ptil)/ntil)
    binomciOutput(p_hat, [lower, upper], n)
end

function acbinomci(x::Vector{Int}; alpha::Real=0.05)
    for i = 1:length(x)
        if x[i]<0 || x[i] > 1
            error("x vector must contain only values 0 or 1.")
        end
    end
    acbinomci(sum(x), length(x), alpha=alpha)
end


"""`_estimate_dispersion(x)`

Estimate dispersion by the following methods. Return the first value that gives
a non-zero dispersion. Each are normalized to 1.0 for Gaussian distributions:

1. Normalized median absolute deviation `mad`,
1. Normalized inter-quartile range `iqrn`,
1. Normalized winsorized variance `winvar`."""
function _estimate_dispersion{S <: Real}(x::AbstractArray{S})
    m =  mad(x)
    m > 0 && return m

    m = iqrn(x)
    m > 0 && return m

    m =  sqrt(winvar(x)./0.4129)
    m > 0 && return m

    error("All measures of dispersion are equal to 0")
end



"""`sint(x; alpha=.05)`
`sint(x, testmedian; alpha=.05)`

Compute the (1-α) confidence interval for the median. In the second form,
use the Hettmansperger and Sheather interpolation method to estimate a p-value
for the `testmedian`."""
function sint{S <: Real}(x::AbstractArray{S}; alpha::Real=0.05, method::Bool=true)
    const n = length(x)
    k = Int(Rmath.qbinom(alpha/2.0, n, 0.5))
    gk = Rmath.pbinom(n-k, n, .5) - Rmath.pbinom(k-1, n, .5)
    if gk < (1 - alpha)
        k = k - 1
        gk = Rmath.pbinom(n-k, n, .5) - Rmath.pbinom(k-1, n, .5)
    end
    gkp1 = Rmath.pbinom(n-k-1, n, .5) - Rmath.pbinom(k, n, .5)
    kp = k + 1

    xsort=sort(x)
    nmk = n-k
    nmkp = nmk+1
    ival = (gk-1+alpha)/(gk-gkp1)
    lam = ((n-k)*ival)/(k+(n-2k)*ival)
    low = lam*xsort[kp]+(1-lam)*xsort[k]
    hi = lam*xsort[nmk]+(1-lam)*xsort[nmkp]
    if method
        METHOD="Confidence interval for the median\n"
        if duplicated(x)
            METHOD *= "Duplicate values detected; hdpb() might have more power\n"
        end
    else
        METHOD=nothing
    end
    output=testOutput()
    output.method=METHOD
    output.ci=[low, hi]
    output
end


function sint{S <: Real}(x::AbstractArray{S}, testmedian;
    alpha::Real=0.05, method::Bool=true)
    ci = sint(x, alpha=alpha, method=false).ci
    med = median(x)
    cichoice = testmedian<med ? 1 : 2

    # Find the pvalue that excludes testmedian by binary search.
    minloga = -8.0
    maxloga = -0.001
    ciA = sint(x, alpha=exp(minloga)).ci[cichoice]-testmedian
    ciB = sint(x, alpha=exp(maxloga)).ci[cichoice]-testmedian
    if ciA*ciB > 0
        if ciB*(med-testmedian)<0
            pval = 1.0
        else
            pval = 0.0
        end
    else
        while (maxloga-minloga > .0001)
            newloga = (maxloga+minloga)/2
            newci = sint(x, alpha=exp(newloga)).ci[cichoice]-testmedian
            if newci*ciB >= 0
                ciB = newci
                maxloga = newloga
            else
                ciA = newci
                minloga = newloga
            end
        end
        pval = exp((maxloga+minloga)/2.0)
    end
    if method
        METHOD="Confidence interval for the median with p-val.\n"
        if duplicated(x)
            METHOD *= "Duplicate values detected; hdpb() might have more power\n"
        end
    else
        METHOD=nothing
    end
    output = testOutput()
    output.method = METHOD
    output.ci     = ci
    output.p      = pval
    output
end


"""`hpsi(x, bend=1.28)`

Evaluate Huber's ψ function for each value in the vector `x`.
ψ(x) = max( min(x,bend), -bend)."""
function hpsi{S <: Real}(x::AbstractArray{S}, bend::Real=1.28)
    ψ = Array(x)
    ψ[x .> bend] = bend
    ψ[x .< -bend] = -bend
    ψ
end


"""`onestep(x, bend=1.28)`

Compute one-step M-estimator of location using Huber's ψ."""
function onestep{S <: Real}(x::AbstractArray{S}, bend::Real=1.28)
    MED = median(x)
    MAD = mad(x)
    y = (x-MED)/MAD
    A = sum(hpsi(y, bend))
    B = sum(abs.(y) .<= bend)
    return MED + MAD*A/B
end

"""`bootstrapci(x; est=onestep, alpha=0.05, nboot=2000, nullvalue=NaN)`

Compute a (1-α) confidence interval for the location-estimator function `est`
using a bootstrap calculation. The default estimator is `onestep`. If `nullvalue` is
given, it is the target value used when computing a p-value.
"""
function bootstrapci{S <: Real}(x::AbstractArray{S}; est::Function=onestep,
    alpha::Real=0.05, nboot::Integer=2000, seed=2, nullvalue::Real=NaN)
    if isa(seed, Int)
        srand(seed)
    elseif seed
        srand(2)
    end
    const n = length(x)
    bvec = zeros(nboot)
    for i = 1:nboot
        randid=rand(1:n, n)
        bvec[i]=est(x[randid])
    end
    low::Int = round((alpha/2)*nboot) + 1
    up = nboot-low + 1
    sort!(bvec)

    pv = NaN
    if nullvalue != NaN
        pv = mean(bvec.>nullvalue)+0.5*mean(bvec.==nullvalue)
        pv = 2min(pv, 1-pv)
    end
    estimate = est(x)
    output = testOutput()
    output.estimate = estimate
    output.ci = [bvec[low], bvec[up]]
    output.p = pv
    output
end

"""`bootstrapse(x; est=median, alpha=0.05, nboot=2000)`

Compute the standard error of the location-estimator function `est`
using a bootstrap calculation. The default estimator is `median`.
"""
function bootstrapse{S <: Real}(x::AbstractArray{S};
        nboot::Integer=1000, est::Function=median, seed=2)
    if isa(seed, Int)
        srand(seed)
    elseif seed
        srand(2)
    end
    const n = length(x)
    bvec = zeros(nboot)
    for i = 1:nboot
        randid=rand(1:n, n)
        bvec[i]=est(x[randid])
    end
    std(bvec)
end


"""`mom(x; bend=2.24)`

Returns a modified one-step M-estimator of location (MOM), which is the unweighted
mean of all values not more than (bend times the `mad(x)`) away from the data
median.
"""
function mom{S <: Real}(x::AbstractArray{S}; bend::Real=2.24)
    mom!(copy(x), bend=bend)
end

"""`mom!(x)`

Like `mom`, but will sort the input vector."""
function mom!{S <: Real}(x::AbstractArray{S}; bend::Real=2.24)
    const n = length(x)
    med = median!(x)
    MAD = mad(x)
    not_extreme = abs.(x-med) .<= bend*MAD
    mean(x[not_extreme])
end


"""`momci(x; bend=2.24, alpha=0.05, nboot=2000)`

Compute a bootstrap, (1-α) confidence interval for the MOM-estimator of location based on Huber's ψ.
The default number of bootstrap resamplings is nboot=2000."""
function momci{S <: Real}(x::AbstractArray{S}; bend::Real=2.24, alpha::Real=0.05,
    nboot::Integer=2000, seed=2, nullvalue::Real=NaN)
    estimator(z) = mom!(z, bend=bend)
    bootstrapci(copy(x), est=estimator, alpha=alpha, nboot=nboot, seed=seed, nullvalue=nullvalue)
end

"""`contam_randn([T=Float64], n; epsilon=0.1, k=10.0)`

Contaminated normal distribution N(0,1). (That is, with μ=0, σ=1.) A fraction `epsilon` of
values will be N(0,`k`)."""
function contam_randn(T::Type, n::Integer; epsilon::Real=0.1, k::Real=10)
    k <= 0 && error("k > 0 is required")
    epsilon > 1 || epsilon < 0 && error("0 ≤ epsilon ≤ 1 is required")
    output = randn(T, n)
    contaminated = rand(n) .< epsilon
    output[contaminated] *= k
    return output
end

contam_randn(n::Integer; epsilon::Real=0.1, k::Real=10) =
    contam_randn(Float64, n, epsilon=epsilon, k=k)



"""`trimpb(x; tr=0.2, alpha=0.05, nboot=2000, win=false, nullvalue=0.0)`

Compute a (1-α) confidence interval for a trimmed mean with a trimming fraction of `tr`.

Use `nboot` bootstrap samples and `alpha` for α.

If `win` is a real number, it is the amount of Winsorizing before bootstrapping.
If `win` is true, then use 10% Winsorizing. If `win` is false, no Winsorizing is done.

The p-value is for the hypothesis that trimmed mean equals `nullvalue`.
"""
function trimpb{S <: Real}(x::AbstractArray{S}; tr::Real=0.2, alpha::Real=0.05, nboot::Integer=2000,
                win=false, nullvalue::Real=0.0, seed=2)
    if isa(win, Bool) && win
        win = 0.1
    end
    if win > tr
        error("trimpb() requires that the amount of Winsorizing ≤ the amount of trimming.")
    end
    wx = winval(x, tr=win)

    estimator(x) = tmean!(x, tr=tr)
    bootstrapci(wx, est=estimator, alpha=alpha, nboot=nboot, seed=seed, nullvalue=nullvalue)
end



"""`procb(x, y; seed=2)`

Compute a (1-α) confidence interval for Pearson's correlation coefficient.

This function uses an adjusted percentile bootstrap method that
gives good results when the error term is heteroscedastic.
"""
function pcorb{S <: Real, T <: Real}(x::AbstractArray{S}, y::AbstractArray{T}; seed=2)
   if isa(seed, Bool)
        seed && srand(2)
    else
        srand(seed)
    end
    const n = length(x)
    # Wow. Every number in this function is totally magic.
    bvec=zeros(Float64, 599)
    for i=1:599
        randid=rand(1:n, n)
        tempx = x[randid]
        tempy = y[randid]
        bvec[i]=cor(tempx, tempy)
    end
    if n >= 250
        ilow, ihi = 15, 584
    elseif n >= 180
        ilow, ihi = 14, 585
    elseif n >= 80
        ilow, ihi = 11, 588
    elseif n >= 40
        ilow, ihi = 8, 592
    else
        ilow, ihi = 7, 593
    end
    sort!(bvec)
    output = testOutput()
    output.estimate = cor(x, y)
    output.ci = [bvec[ilow], bvec[ihi]]
    output
end



"""`yuend(x,y; tr=0.2, alpha=0.05)`

Compare the trimmed means of two dependent random variables `x` and `y`.
The default amount of trimming `tr` is 20%.

A (1-α) confidence interval for the difference of trimmed mean of `x` minus
the trimmed mean of `y` is computed and returned in `yuend.ci`.
The significance level is returned in `yuend.siglevel`.
"""
function yuend{S <: Real, T <: Real}(x::AbstractArray{S}, y::AbstractArray{T};
        tr::Real=0.2, alpha::Real=0.05, method::Bool=true)
    const n = length(x)
    if n != length(y)
        error("`x` and `y` must agree in length")
    end
    h1::Integer = n - 2*floor(tr*n)
    q1 = (n - 1)*winvar(x, tr=tr)
    q2 = (n - 1)*winvar(y, tr=tr)
    q3 = (n - 1)*wincov(x, y, tr=tr)
    df = h1 - 1
    se = sqrt((q1 + q2 - 2q3)/(h1*(h1-1)))
    crit = Rmath.qt(1 - alpha/2, df)
    meandif = tmean(x, tr=tr) - tmean(y, tr=tr)
    confint = [meandif - crit*se, meandif + crit*se]
    test = meandif/se
    p = 2*(1 - Rmath.pt(abs(test), df))
    if method
        METHOD="Comparing the trimmed means of two dependent variables.\n"
    else
        METHOD=nothing
    end
    output = testOutput()
    output.method = METHOD
    output.ci = confint
    output.p = p
    output.estimate = meandif
    output.se = se
    output.statistic = test
    output.n = n
    output.df = df
    output
end

function pbos{S <: Real}(x::AbstractArray{S}; beta::Real=0.2)
    temp    = sort( abs.( x - median(x) ))
    nval    = length( x )
    omhatid::Integer = floor( (1 - beta)*nval )
    omhatx  = temp[ omhatid ]
    psi     = ( x - median(x) )./ omhatx
    i1      = length(psi[ psi .< -1 ])
    i2      = length(psi[ psi .> 1 ])
    sx      = 0.0
    [ sx += psi[i] < -1 ? 0 : [ psi[i] > 1 ? 0 : x[i] ]  for i=1:nval ]
    return ( sx  + omhatx * (i2 - i1))/(nval - i1 - i2)
end


#Compute the percentage bend correlation between x and y
#beta is the bending constant for omega sub N.
function pbcor{S <: Real, T <: Real}(x::AbstractArray{S}, y::AbstractArray{T}; beta::Real=0.2)
    nval = length(x)
    if length(y) != nval
        error("x and y do not agree in length.")
    end
    temp    = sort( abs.( x - median(x) ))
    omhatid::Integer = floor( (1 - beta)*nval )
    omhatx  = temp[ omhatid ]
    temp    = sort( abs.( y - median(y) ))
    omhaty  = temp[ omhatid ]
    a       = (x .- pbos(x, beta=beta) )./omhatx
    b       = (y .- pbos(y, beta=beta) )./omhaty
    for i = 1:nval
        if a[i] < -1
            a[i] = -1
        elseif a[i] > 1
            a[i] = 1
        end
        if b[i] < -1
            b[i] = -1
        elseif b[i] > 1
            b[i] = 1
        end
    end
    Pbcor   = sum( a.*b )/sqrt(sum( a.*a ) * sum( b.*b ))
    test    = Pbcor*sqrt( ( nval - 2 )/( 1 - Pbcor*Pbcor ) )
    sig     = 2*( 1 - Rmath.pt(abs(test), nval-2))

    METHOD=nothing

    output = testOutput()
    output.method = METHOD
    output.p = sig
    output.estimate = Pbcor
    output.statistic = test
    output.n = nval
    output.df = nval - 2
    output
end
