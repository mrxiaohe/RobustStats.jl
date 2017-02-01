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
    METHOD      = method ? "1-alpha confidence interval for the trimmed mean\n": nothing
    output           = testOutput()
    output.method    = METHOD
    output.df        = df
    output.estimate  = estimate
    output.ci        = confint
    output.statistic = statistic
    output.p         = pval
    return output
end


# We have no idea what this does.
function stein1{S <: Real}(x::AbstractArray{S}, delta::Real; alpha::Real=0.05,
    pow::Real=0.8, oneside::Bool=false, n=nothing, variance=nothing)
    delta = abs(delta)
    if n==nothing; n=length(x); end
    if variance==nothing; variance = var(x); end
    if !oneside
        alpha = alpha/2
    end
    df       = n-1
    d        = (delta/(Rmath.qt(pow, df)-Rmath.qt(alpha, df)))*(delta/(Rmath.qt(pow, df)-Rmath.qt(alpha, df)))
    N        = max(n, floor(Int64, variance/d)+1)
    return N
end


# We have no idea what this does.
function stein2{S <: Real, T <: Real}(x1::AbstractArray{S}, x2::AbstractArray{T}; mu0::Real=0.0, alpha::Real=0.05, method::Bool=true)
    n       = length(x1)
    df      = n-1
    N       = n+length(x2)
    xbar    = mean([x1, x2])
    test    = sqrt(N)*(xbar-mu0)/std(x1)
    crit    = Rmath.qt(1-alpha/2, df)
    lo, hi  = xbar-crit*std(x1), xbar+crit*std(x1)
    sig     = 2*(1-Rmath.pt(test, df))
    METHOD  = method?"The 2nd stage of Stein's method":nothing
    output  = testOutput()
    output.method    = METHOD
    output.df        = df
    output.estimate  = xbar
    output.ci        = [lo, hi]
    output.statistic = test
    output.crit      = crit
    output.p         = sig
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
    absdev = abs(x-med)
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
        if (abs(ψ) >= 1.0)
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
        if (u<1.0)
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
    flag = (x.<cl) | (x.>cu)
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


#Compute a 1-alpha confidence interval for p, the probability of success for a binomial dist. using Pratt's method
#y is a vector of 1's and 0's, x is the number of successes observed among n trials.
"""`binomci(s, n; alpha=0.05)`

Compute the (1-α) confidence interval for p, the binomial probability of success, given
`s` successes in `n` trials. Returns an object with components `p_hat` (the observed
fraction of successes) and `confint=[lo,hi]` (the confidence interval). The computation
uses Pratt's method.

Can also use `binomci(x; alpha=0.05)`, where x is an array consisting only of 0s
and 1s. It's equivalent to `binomci(sum(x), length(x), alpha=alpha)`.
"""
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
    end
    if s == 1
        lower = 1-(1-alpha/2).^(1/n)
        upper = 1-(alpha/2).^(1/n)
        return binomciOutput(p_hat, [lower, upper], n)
    end
    if s == (n-1)
        lower = (alpha/2).^(1/n)
        upper = (1-alpha/2).^(1/n)
        return binomciOutput(p_hat, [lower, upper], n)
    end
    if s == n
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

#Compute adaptive kernel density estimate for univariate data
function akerd{S <: Real}(x::AbstractArray{S}; hval::Real=NaN, aval::Real=0.5, op::Integer=1,
               fr::Real=0.8, pts=NaN, plotit=true, xlab="",
               ylab="", title="", color="black")
    if isnan(pts)
       pts = x[:]
    end
    pts  = sort!(pts)

    if op == 1
        m = _estimate_dispersion(x)
        fhat  = rdplot(x, pts=pts, plotit=false, fr=fr)
    # elseif op == 2
    #     init = kde(x) ## NOT DEFINED!
    #     fhat = init.density
    #     x    = init.x
    end
    const n = length(x)
    if isnan(hval)
        A = min(std(x), iqrn(x))
        if A==0.0; A = winstd(x)/0.64; end
        hval = 1.06*A/n^0.2
    end
    gm = 0.0
    gm_int = 0
    const nfhat = length(fhat)
    for i = 1:nfhat
        if fhat[i] > 0.0
            gm += log(fhat[i])
            gm_int += 1
        end
    end
    gm = exp(gm/gm_int)
    alam = (fhat/gm).^(-aval)
    dhat = akerd_loop(x, pts, hval, alam)
    if plotit
        plot(pts, dhat, color=color)
        plt[:title](title)
        plt[:xlabel](xlab)
        plt[:ylabel](ylab)
    end
    dhat
end

#Expected frequency curve. fr controls amount of smoothing, theta is the azimuthal direction and
#phi the colatitude
function rdplot{S <: Real}(x::AbstractArray{S}; fr::Real=NaN, pts=NaN,
                           plotit=true, title="", xlab="", ylab="", color="black")
    if isnan(fr); fr = 0.8; end
    if isnan(pts); pts = x[:];end
    rmd = [sum(near(x, pts[i], fr))*1.0 for i=1:length(pts)]
    rmd /= length(x)
    MAD = mad(x)
    if MAD != 0.0
        rmd /= 2fr*MAD
    end
    if plotit
        index = sortperm(pts);
        clf()
        plot(pts[index], rmd[index], color=color)
        plt[:title](title)
        plt[:xlabel](xlab)
        plt[:ylabel](ylab)
    end
    rmd
end


"""`near(x, pt, fr=1.0)`

Determine which values in `x` are near `pt`. Return a BitArray giving whether
each value of `x` is within `fr*m` of `pt`, where `m` is the dispersion measure
returned by `_estimate_dispersion(x)`"""
function near{S <: Real}(x::AbstractArray{S}, pt::Real, fr::Real=1.0)
    m = _estimate_dispersion(x)
    return abs(x-pt) .<= fr*m
end


#Computing the confidence interval for the median using the
#Hettmansperger and Sheather method.
function sint{S <: Real}(x::AbstractArray{S}; alpha::Real=0.05, method::Bool=true)
    n  = length(x)
    k  = Rmath.qbinom(alpha/2.0, n, 0.5)
    gk = Rmath.pbinom(length(x)-k, n, .5) - Rmath.pbinom(k-1, n, .5)
    if gk >= (1 - alpha)
        gkp1 = Rmath.pbinom(n-k-1, n, .5) - Rmath.pbinom(k, n, .5)
        kp   = k + 1
    else
        k    = k - 1
        gk   = Rmath.pbinom(n-k, n, .5)   - Rmath.pbinom(k-1, n, .5)
        gkp1 = Rmath.pbinom(n-k-1, n, .5) - Rmath.pbinom(k, n, .5)
        kp   = k+1
    end
    xsort=sort(x)
    nmk  = n-k
    nmkp = nmk+1
    ival = (gk-1+alpha)/(gk-gkp1)
    lam  = ((n-k)*ival)/(k+(n-2*k)*ival)
    low  = lam*xsort[kp]+(1-lam)*xsort[k]
    hi   = lam*xsort[nmk]+(1-lam)*xsort[nmkp]
    if method
        if duplicated(x)
            METHOD="Confidence interval for the median\nDuplicate values detected; hdpb() might have more power\n"
        else
            METHOD="Confidence interval for the median\n"
        end
    else
        METHOD=nothing
    end
    output=testOutput()
    output.method=METHOD
    output.ci=[low, hi]
    output
end


#Computing a 1-alpha confidence interval for p, the probability of
#success for a binomial distribution, using a generalization of the
#Agresti-Coull  method that was studied by Brown, Cai DasGupta
#(Annals of Statistics, 2002, 30, 160-201.)
function acbinomci(x::Int, n::Int; alpha::Real=0.05)
    if x > n
        error("x must be smaller than or equal to n")
    elseif x < 0
        error("x cannot be negative")
    elseif n == 1
        error("Something is wrong: number of observations is only 1")
    end
    if x != n && x != 0
        cr    = Rmath.qnorm(1-alpha/2)
        ntil  = n+cr.*cr
        ptil  = (x+cr.*cr./2)./ntil
        lower = ptil-cr.*sqrt(ptil.*(1-ptil)./ntil)
        upper = ptil+cr.*sqrt(ptil.*(1-ptil)./ntil)
    end
    if x == 0
        lower = 0.0
        upper = 1.0-alpha.^(1/n)
    end
    if x == 1
        lower = 1-(1-alpha/2).^(1/n)
        upper = 1-(alpha/2).^(1/n)
    end
    if x == (n-1)
        lower = (alpha/2).^(1/n)
        upper = (1-alpha/2).^(1/n)
    end
    if x == n
        lower = alpha.^(1/n)
        upper = 1
    end
    phat = x/n
    binomciOutput(phat, [lower, upper], n)
end
acbinomci(x::Vector{Int}; alpha::Real=0.05)=
sum([x[i]!=0 && x[i]!=1 for i=1:length(x)])>0?error("Data must be 0's and 1's"):acbinomci(length(x.==1), length(x), alpha=alpha)


#Extension of Stein's method based on the trimmed mean.
function stein1_tr{S <: Real}(x::AbstractArray{S}, del::Real; alpha::Real=0.05, pow::Real=0.8, tr::Real=0.2)
    if tr <0 || tr >= .5
        error("Argument tr must be between 0 and .5")
    end
    n =length(x)
    g::Int    = floor(Int64, tr*n)
    df::Int   = n-2*g-1
    t1        = Rmath.qt(pow, df)
    t2        = Rmath.qt(alpha./2, df)
    dv        = (del/(t1-t2))*(del/(t1-t2))
    nvec::Int = floor(Int64, trimse(x, tr=tr)./dv)+1
    n > nvec ? n : nvec
end

function stein1_tr{S <: Real}(x::Array{S, 2}, del::Real; alpha::Real=0.05, pow::Real=0.8, tr::Real=0.2)
    if tr<0||tr>=.5
        error("Argument tr must be between 0 and .5")
    end
    if size(x, 2)==1 || size(x, 1)==1
        temp=zeros(length(x))
        temp[:]=x
        return stein1_tr(temp, del, alpha=alpha, pow=pow, tr=tr)
    else
        J     = size(x,2)
        ic    = 0
        ntest::Integer = (J*J-J)/2
        n     = size(x, 1)
        g::Integer = floor(Int64, tr*n)
        df    = n-2g-1
        tdiff = Rmath.qt(pow, df)-Rmath.qt(alpha/2.0/ntest, df)
        dv    = (del/tdiff)*(del/tdiff)
        N     = zeros(Int, ntest)
        if ntest > 1
            for i = 1:J
                for j = 1:J
                    if i < j
                        N[ic+=1] = floor(Int64, trimse(x[:,i]-x[:,j], tr=tr)/dv)+1
                    end
                end
            end
            return max(vcat(N, n))
        elseif ntest == 1
            return floor(Int64, trimse(x[:,1], tr=tr)/dv)+1
        end
    end
end



#Extension of the second stage of Stein's method when performing all pairwise comparisons
#among J dependent groups.
function stein2_tr{S <: Real, T <: Real}(x::AbstractArray{S}, y::AbstractArray{T}; alpha::Real=0.05, tr::Real=0.2, method::Bool=true)
    if tr < 0 || tr >= .5
        error("Argument tr must be between 0 and .5")
    end
    m       = vcat(x, y)
    METHOD  = method? "Extension of the 2nd stage of Stein's method based on the trimmed mean\n":nothing
    output  = testOutput()
    output.method    = METHOD
    output.statistic = sqrt(length(m))*(1-2tr)*tmean(m, tr=tr)/winstd(m)
    output.crit         = Rmath.qt(1-alpha/2, floor(length(x)-2floor(tr*length(x))-1))
    return output
end

function stein2_tr(x::Array, y::Array; alpha::Real=0.05, tr::Real=.2, method::Bool=true)
    if tr < 0 || tr >= .5
        error("Argument tr must be between 0 and .5")
    end
    if size(x, 2) != size(y, 2)
        error("Number of variables in x not equal to number of variables in y")
    end
    if size(x, 2) == 1 || size(x, 1) ==1
        return stein2_tr(x[:], del, alpha=alpha, pow=pow, tr=tr, method=method)
    else
        g::Int  = floor(tr*size(x,1))
        df::Int = size(x,1)-2.*g-1
        ic::Int = 0
        J::Int  = size(x,2)
        ntest::Int = (J.*J-J)./2
        m       = vcat(x, y)
        nm      = size(m, 1)
        test=DataFrame([Integer, Integer, Float64], ["grp1", "grp2", "test_stat"], ntest)
        if ntest > 1
            for i=1:J
                for j=1:J
                    if i<j
                        temp  = m[:,i]-m[:,j]
                        test[ic+=1, 1], test[ic, 2] = i, j
                        test[ic, 3] = sqrt(nm)*(1-2tr)*tmean(temp, tr=tr)/winstd(temp)
                    end
                end
            end
        else ntest == 1
            test[1,1], test[1, 2] = 1, 2
            test[1, 3] = test[ic, 3] = sqrt(nm)*(1-2tr)*tmean(m[:,1], tr=tr)/winstd(m[:,1])
        end
        if method
            METHOD="Extension of the 2nd stage of Stein's method based on the trimmed mean\n"
        else
            METHOD=nothing
        end
        output = testOutput()
        output.method    = METHOD
        output.statistic = test
        output.crit      = Rmath.qt(1-alpha/2./ntest, df)
        return output
    end
end


#Compute 1-alpha confidence interval for the median using the Hettmansperger-Sheather interpolation method
function sintv2{S <: Real}(x::AbstractArray{S}; alpha::Real=0.05, nullval::Real=0.0, method::Bool=true)
    ci      = sint(x, alpha=alpha, method=false).ci
    alph    = [1:99]./100
    irem    = zeros(Int, 1)
    for i = 1:99
        irem  = i
        chkit = sint(x, alpha=alph[i], method=false).ci
        if chkit[1].>nullval || chkit[2].<nullval
            break
        end
    end
    pval=irem/100
    if pval.<=0.1
        iup  = (irem+1)/100
        alph = seq(0.001, iup, 0.001)
        for i = 1:length(alph)
            pval  = alph[i]
            chkit = sint(x, alpha=alph[i], method=false).ci
            if chkit[1] > nullval || chkit[2] < nullval
                break
            end
        end
    end
    if pval <= 0.001
        alph = seq(.0001,.001,.0001)
        for i = 1:length(alph)
            pval  = alph[i]
            chkit = sint(x, alpha=alph[i], method=false).ci
            if chkit[1] > nullval || chkit[2] < nullval
                break
            end
        end
    end
    if method
        if duplicated(x)
            METHOD="Confidence interval for the median\nusing the Hettmansperger-Sheather interpolation method\nDuplicate values detected; hdpb() might have more power\n"
        else
            METHOD="Confidence interval for the median\nusing the Hettmansperger-Sheather interpolation method\n"
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
    B = sum(abs(y) .<= bend)
    return MED + MAD*A/B
end

#Compute a bootstrap, .95 confidence interval for the
#measure of location corresponding to the argument est.
#By default, a one-step
#M-estimator of location based on Huber's ψ is used.
#The default number of bootstrap samples is nboot=500
#
#nv=null value when  computing a p-value

function onesampb{S <: Real}(x::AbstractArray{S}; est::Function=onestep,
    alpha::Real=0.05, nboot::Integer=2000, seed=2, nv::Real=0)
    if isa(seed, Int)
        srand(seed)
    elseif seed
        srand(2)
    end
    n        = length(x)
    bvec     = zeros(nboot)
    temp     = zeros(n)
    randid   = rand(1:n, n*nboot)
    for i = 1:nboot
        [temp[j] = x[randid[(i-1)*n+j]] for j=1:n]
        bvec[i]  = est(temp)
    end
    low      = round((alpha./2).*nboot) + 1
    up       = nboot-low + 1
    pv       = mean(bvec.>nv)+0.5.*mean(bvec.==nv)
    pv       = 2.*min(pv, 1-pv)
    bvec     = sort!(bvec)
    estimate = est(x)
    output   = testOutput()
    output.estimate = estimate
    output.ci       = [bvec[low], bvec[up]]
    output.p        = pv
    output
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
    not_extreme = abs(x-med) .<= bend*MAD
    mean(x[not_extreme])
end


"""`momci(x; bend=2.24, alpha=0.05, nboot=2000)`

Compute a bootstrap, (1-α) confidence interval for the MOM-estimator of location based on Huber's ψ.
The default number of bootstrap resamplings is nboot=2000."""
function momci{S <: Real}(x::AbstractArray{S}; bend::Real=2.24, alpha::Real=0.05,
    nboot::Integer=2000, seed=2, method=true)
    if isa(seed, Int)
        srand(seed)
    elseif seed
        srand(2)
    end
    const n = length(x)
    bvec = zeros(nboot)
    temp = zeros(n)
    randid=rand(1:n, n, nboot)
    for i = 1:nboot
        for j = 1:n
            temp[j] = x[randid[j,i]]
        end
        bvec[i]=mom!(temp, bend=bend)
    end
    low::Int = round((alpha/2)*nboot) + 1
    up = nboot-low + 1
    sort!(bvec)
    output = testOutput()
    output.ci = [bvec[low], bvec[up]]
    if method
        x = 1-alpha
        output.method = "Bootstrap $(x) confidence interval for the MOM-estimator of location based on Huber's ψ.\n"
    end
    output
end

#Contaminated normal distribution
function cnorm(n::Integer; epsilon::Real=0.1, k::Real=10)
    if epsilon > 1
        error("epsilon must be less than or equal to 1")
    elseif epsilon < 0
        error("epsilon must be greater than or equal to 0")
    end
    if k <= 0
        error("k must be greater than 0")
    end
    output  = zeros(n)
    epsilondiff = 1.0 - epsilon
    [output[i] = rand() > epsilondiff ? k*randn():randn() for i=1:n]
    return output
end


#   Compute a 1-alpha confidence interval for
#   a trimmed mean.
#
#   The default number of bootstrap samples is nboot=2000
#
#   win is the amount of Winsorizing before bootstrapping
#   when WIN=T.
#
#   Missing values are automatically removed.
#
#  nv is null value. That test hypothesis trimmed mean equals nv
#
#  plotit=TRUE gives a plot of the bootstrap values
#  pop=1 results in the expected frequency curve.
#  pop=2 kernel density estimate    NOT IMPLEMENTED
#  pop=3 boxplot                    NOT IMPLEMENTED
#  pop=4 stem-and-leaf              NOT IMPLEMENTED
#  pop=5 histogram
#  pop=6 adaptive kernel density estimate.
#
#  fr controls the amount of smoothing when plotting the bootstrap values
#  via the function rdplot. fr=NA means the function will use fr=.8
#  (When plotting bivariate data, rdplot uses fr=.6 by default.)

function trimpb{S <: Real}(x::AbstractArray{S}; tr::Real=0.2, alpha::Real=0.05, nboot::Integer=2000,
                win=false, plotit::Bool=false, pop::Int=1, nullval::Real=0.0,
                xlab="X", ylab="Density", fr=nothing, seed=2,
                method::Bool=true)
    if isa(win, Bool)
        if win
            x=winval(x, tr=0.1)
        end
    elseif win>tr
        error("The amount of Winsorizing must be <= to the amount of trimming")
    else
        x=winval(x, tr=win)
    end
    crit=alpha/2.0
    icl=round(crit*nboot)+1
    icu=nboot-icl
    if isa(seed, Bool)
        if seed
            srand(2)
        end
    else
        srand(seed)
    end
    n=length(x)
    bvec=zeros(nboot)
    randid=rand(1:n, n*nboot)
    for i=1:nboot
        temp=zeros(n)
        for j=1:n
            temp[j]=x[randid[(i-1)*n+j]]
        end
        bvec[i]=tmean!(temp, tr=tr)
    end
    bvec=sort!(bvec)
    pval1=pval2=0.0
    for i=1:nboot
        if bvec[i]<nullval
            pval1+=1./nboot
        end
        if bvec[i]==nullval
            pval2+=0.5/nboot
        end
    end
    pval=2*min(pval1+pval2, 1-pval1-pval2)
    ci=[bvec[icl], bvec[icu]]
    if method
        METHOD::String="Compute a 1-alpha confidence interval for a trimmed mean using the bootstrap percentile method.\n"
    else
        METHOD=nothing
    end
    if plotit
        if pop==1
            if fr==nothing
                rdplot(bvec, fr=0.6, xlab=xlab, ylab=ylab)
            else
                rdplot(bvec, fr=fr, xlab=xlab, ylab=ylab)
            end
        elseif pop==5
            p=FramedPlot()
            add(p, Histogram(hist(bvec)[2], 2))
            if xlab!=nothing
                setattr(p, "xlabel", xlab)
            end
            if ylab!=nothing
                setattr(p, "ylabel", ylab)
            end
            Winston.tk(p)
        elseif pop==6
            akerd(bvec, xlab=xlab, ylab=ylab)
        end
    end
    output = testOutput()
    output.method = METHOD
    output.ci     = ci
    output.p      = pval
    output
end



#  Compute a 1-alpha confidence interval for the trimmed mean
#  using a bootstrap percentile t method.
#
#  The default amount of trimming is tr=.2
#  side=T, for true,  indicates the symmetric two-sided method
#
#
#  Side=F yields an equal-tailed confidence interval
#
#
#  NOTE: p.value is reported when side=T only.
#

function trimcibt{S <: Real}(x::AbstractArray{S}; tr::Real=0.2, alpha::Real=0.05, nboot::Integer=2000, side::Bool=true,
                             plotit::Bool=false, op::Integer=1, nullval::Real=0, seed=2, method::Bool=true)
    if isa(seed, Bool)
        if seed
            srand(2)
        end
    else
        srand(seed)
    end
    n=length(x)
    test=(tmean(x, tr=tr)-nullval)./trimse(x, tr=tr)
    randid=rand(1:n, n*nboot)
    tempout=trimcibt_loop(x, n, nboot, tr, side, randid, test)
    if side
        tval=tempout[1]
        pval=tempout[2]
    end
    icrit=floor((1-alpha)*nboot+0.5)
    ibot=round(alpha*nboot/2)+1
    itop=nboot-ibot-1
    ci=zeros(2)
    if method && side
        METHOD="Bootstrap .95 confidence interval for the trimmed mean\nusing a bootstrap percentile t method\n"
    elseif method && !side
        METHOD="Bootstrap .95 confidence interval for the trimmed mean\nusing a bootstrap percentile t method\n[NOTE: p value is computed only when side=true]\n"
    else
        METHOD=nothing
    end
    if !side
        if plotit
            if op==1
                akerd(tempout)
            elseif op==2
                rdplot(tempout)
            end
        end
        ci[1]=tmean(x, tr=tr)-tempout[itop]*trimse(x, tr=tr)
        ci[2]=tmean(x, tr=tr)-tempout[ibot]*trimse(x, tr=tr)
        output = testOutput()
        output.method = METHOD
        output.estimate = tmean(x, tr=tr)
        output.ci = ci
        output.statistic = test
        return output
    else
        if plotit
            if op==1
                akerd(tval)
            elseif op==2
                rdplot(tval)
            end
        end
        ci[1]=tmean(x, tr=tr)-tval[icrit]*trimse(x, tr=tr)
        ci[2]=tmean(x, tr=tr)+tval[icrit]*trimse(x, tr=tr)
        output = testOutput()
        output.method = METHOD
        output.estimate = tmean(x, tr=tr)
        output.ci = ci
        output.statistic = test
        output.p = pval
        return output
    end
end

#   Compute bootstrap estimate of the standard error of the
#   estimator est
#   The default number of bootstrap samples is nboot=1000

function bootse{S <: Real}(x::AbstractArray{S}; nboot::Integer=1000, est::Function=median, seed=2)
   if isa(seed, Bool)
        if seed
            srand(2)
        end
    else
        srand(seed)
    end
    n=length(x)
    temp=zeros(n)
    bvec=zeros(nboot)
    randid=rand(1:n, n*nboot)
    for i=1:(nboot*n)
        if (i%n)!=0
            temp[i%n]=x[randid[i]]
        else
            temp[n]=x[randid[i]]
            bvec[div(i, n)]=est(temp)

        end
    end
    return std(bvec)
end

#   Compute a .95 confidence interval for Pearson's correlation coefficient.
#
#   This function uses an adjusted percentile bootstrap method that
#   gives good results when the error term is heteroscedastic.
function pcorb{S <: Real, T <: Real}(x::AbstractArray{S}, y::AbstractArray{T}; seed=2, plotit::Bool=false)
   if isa(seed, Bool)
        if seed
            srand(2)
        end
    else
        srand(seed)
    end
    n=length(x)
    randid=rand(1:n, n*599)
    tempx=zeros(n)
    tempy=zeros(n)
    bvec=zeros(599)
    for i=1:(599*n)
        if (i%n)!=0
            tempx[i%n], tempy[i%n]=x[randid[i]], y[randid[i]]
        else
            tempx[n], tempy[n]=x[randid[i]], y[randid[i]]
            bvec[div(i, n)]=cor!(tempx, tempy)
        end
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
    bvec=sort!(bvec)
    if plotit
        akerd(bvec, title="Distribution of Bootstrap Pearson Correlation Coefficients")
    end
    r=cor!(x, y)
    output = testOutput()
    output.estimate = r
    output.ci = [bvec[ilow], bvec[ihi]]
    output
end



# Test the hypothesis of independence between x and y by
# testing the hypothesis that the regression surface is a horizontal plane.
# Stute et al. (1998, JASA, 93, 141-149).
#
#  flag=1 gives Kolmogorov-Smirnov test statistic
#  flag=2 gives the Cramer-von Mises test statistic
#  flag=3 causes both test statistics to be reported.
#
#  tr=0 results in the Cramer-von Mises test statistic when flag=2
#      With tr>0, a trimmed version of the test statistic is used.
#
function indt(x::AbstractArray, y::AbstractArray; flag::Int=1, nboot::Int=599,
              tr::Float64=0.0, seed=2, method=true)
    if ndims(x)==1
        x=reshape(x,length(x),1)
    end
    if length(findin(flag, 1:3))==0
        error("flag must be set to 1, 2, or 3")
    end
    n=size(x,1)
    np=size(x,2)
    y=reshape(y, length(y), 1)
    if length(y)!=n
        error("Incondistent dimensions of x and y: number of x must match number of y")
    end
    mflag=indt_mflag(x)
    yhat=mean(y)
    res=zeros(n)
    [res[i]=y[i]-yhat for i=1:n]
    if isa(seed, Bool)
        if seed
            srand(2)
        end
    else
        srand(seed)
    end
    data=(rand(nboot, n)-0.5).*sqrt(12)
    rvalb=zeros(n, nboot)
    const sqrtn=sqrt(n)
    [rvalb[:,i]=regts1(data[i,:], yhat, res, mflag, x, 0) for i=1:nboot]
    [[rvalb[i,j]=abs(rvalb[i,j])/sqrtn for i=1:n] for j=1:nboot]

    dstatb=zeros(nboot)
    [dstatb[i]=max(rvalb[:,i]) for i=1:nboot]

    [[rvalb[i,j]=rvalb[i,j].*rvalb[i,j] for i=1:n] for j=1:nboot]
    wstatb=mean(rvalb, 1)
    rval=regts1(fill(1.0, n), yhat, res, mflag, x, tr)./sqrtn

    dstat=pval_d=wstat=pval_w=nothing
    if flag==1 || flag==3
        [rval[i]=abs(rval[i]) for i=1:n]
        dstat=max(rval)
        pval_d=0.0
        pval_d=sum(dstat.>=dstatb)./nboot
        pval_d=1.0-pval_d
    end
    if flag==2 || flag==3
        [rval[i]=rval[i].*rval[i] for i=1:n]
        wstat=tmean(rval, tr=tr)
        pval_w=0.0
        pval_w=sum(wstat.>=wstatb)./nboot
        pval_w=1.0-pval_w
    end
    if method
        METHOD::String="Test whether x and y are independent by testing the hypothesis\nthat the regression surface is a horizontal plane.\n"
    else
        METHOD=nothing
    end
    indtOutput(
        METHOD,
        dstat,
        pval_d,
        wstat,
        pval_w,
        flag
        )
end


function indirectTest{S <: Real, T <: Real, W <: Real}(dv::AbstractArray{S}, iv::AbstractArray{T}, m::Vector{W};
            nboot::Integer=5000, alpha::Real=0.05, scale::Bool=false, seed=2, plotit::Bool=false)
    if isa(seed, Bool)
        if seed
            srand(2)
        end
    else
        srand(seed)
    end
    n = length(iv)
    randid  = rand(1:n, n*nboot)
    bvec    = sort!(bootindirect(iv, dv, m, nboot))
    bbar    = mean(bvec)
    bootci  = [bvec[round(alpha*nboot/2)], bvec[nboot - round(alpha*nboot/2) + 1]]
    bootest = mean(bvec)
    bootse  = std(bvec)
    p       = mean(bvec .<0) + 0.5*mean(bvec .==0)
    p       = 2*min(p, 1-p)
    data    = DataFrame(iv, m, dv)
    regfit1 = coeftable(lm( :(x3 ~ x1     ), data))
    regfit2 = coeftable(lm( :(x2 ~ x1     ), data))
    regfit3 = coeftable(lm( :(x3 ~ x1 + x2), data))
    regfit  = rbind(regfit1[2,:], regfit2[2,:], regfit3[2:3,:])

    estimate  = regfit2[2,1]*regfit3[3,1]
    sobel_se  = sqrt(regfit[4, 1]*regfit[4, 1]*regfit[2, 2]*regfit[2, 2]+
                     regfit[2, 1]*regfit[2, 1]*regfit[4, 2]*regfit[4, 2]+
                     regfit[2, 2]*regfit[2, 2]*regfit[4, 2]*regfit[4, 2])
    sobel     = DataFrame(estimate,
                          sobel_se,
                          estimate/sobel_se,
                          estimate - Rmath.qnorm(.975)*sobel_se,
                          estimate + Rmath.qnorm(.975)*sobel_se,
                          2*(1-Rmath.pnorm(abs(estimate/sobel_se))))
    #return nboot, n, regfit, sobel, bootest, bootci, p
    if plotit
        dens = kde(bvec)
        plot(dens.x, dens.density)
    end
    output = indirectTestOutput()
    output.nboot   = nboot
    output.n       = n
    output.sobel   = sobel
    output.regfit  = regfit
    output.bootest = bootest
    output.bootci  = bootci
    output.p       = p
    output.bootse  = bootse
    output
end


#   Compute the Winsorized correlation between x and y.
#
#   tr is the amount of Winsorization
#   This function also returns the Winsorized covariance
function wincor{S <: Real, T <: Real}(x::AbstractArray{S}, y::AbstractArray{T}; tr::Real=0.2)
    n = length(x)
    if n != length(y)
        error("`x` and `y` must agree in length")
    end
    g::Integer = floor(tr*n)
    xvec = winval(x, tr=tr)
    yvec = winval(y, tr=tr)
    wcor = cor(xvec, yvec)
    wcov = cov(xvec, yvec)

    if sum(x.==y) != n
        test = wcor*sqrt((n - 2)/(1 - wcor*wcor))
        sig  = 2*(1 - Rmath.pt(abs(test), n-2*g-2))
        return wcor, wcov, sig, n
    else
        return wcor, wcov, n
    end
end



#
#  Compare the trimmed means of two dependent random variables
#  using the data in x and y.
#  The default amount of trimming is 20%
#
#  Missing values (values stored as NA) are not allowed.
#
#  A confidence interval for the trimmed mean of x minus the
#  the trimmed mean of y is computed and returned in yuend$ci.
#  The significance level is returned in yuend$siglevel
function yuend{S <: Real, T <: Real}(x::AbstractArray{S}, y::AbstractArray{T}; tr::Real=0.2, alpha::Real=0.05, method::Bool=true)
    n = length(x)
    if n != length(y)
        error("`x` and `y` must agree in length")
    end
    h1::Integer = n - 2*floor(tr*n)
    q1 = (n - 1)*winvar(x, tr=tr)
    q2 = (n - 1)*winvar(y, tr=tr)
    q3 = (n - 1)*wincor(x, y, tr=tr)[2]
    df = h1 - 1
    se = sqrt((q1 + q2 - 2*q3)/(h1*(h1-1)))
    crit = Rmath.qt(1 - alpha/2, df)
    dif = tmean(x, tr=tr) - tmean(y, tr=tr)
    confint = [dif - crit*se, dif + crit*se]
    test = dif/se
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
    output.estimate = dif
    output.se = se
    output.statistic = test
    output.n = n
    output.df = df
    output
end

#  A heteroscedastic one-way ANOVA for trimmed means using a generalization of Welch's method.

function t1way{S <: Real}(x::Array{S, 2}; tr::Real=0.2, method::Bool=true)
    n = size(x, 1)
    g = [1:size(x, 2)]
    grp = rep(g, rep(n, size(x, 2)))
    x = x[:]
    t1waycore(x, grp, tr, method)
end

function t1way{S <: Real}(x::AbstractArray{S}, grp::Vector; tr::Real=0.2, method::Bool=true)
    g = unique(grp)
    grpcopy = [find(g.==grp[i])[1] for i=1:length(grp)]
    t1waycore(x, grpcopy, tr, method)
end

function pbos{S <: Real}(x::AbstractArray{S}; beta::Real=0.2)
    temp    = sort( abs( x - median(x) ))
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
    temp    = sort( abs( x - median(x) ))
    omhatid::Integer = floor( (1 - beta)*nval )
    omhatx  = temp[ omhatid ]
    temp    = sort( abs( y - median(y) ))
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

function outer{S <: Real, T <: Real}(x::AbstractArray{S}, y::AbstractArray{T}, f::Function)
    nx      = length(x)
    ny      = length(y)
    output  = zeros(nx, ny)
    for i = 1:ny
        for j = 1:nx
            output[j, i] = f(x[j], y[i])
        end
    end
    return output
end

function hd{S <: Real}(x::AbstractArray{S}; q::Real=0.5)
    #Compute the Theil-Sen regression estimator.
    # Only a single predictor is allowed in this version
    const n = length(x)
    m1   = ( n + 1 )*q
    m2   = ( n + 1 )*(1 - q)
    vec1 = [1:n]./n
    vec2 = ([1:n] - 1)./n
    w    = Rmath.pbeta( vec1, m1, m2 ) - Rmath.pbeta( vec2,  m1, m2 )
    return sum( w.*sort(x) )
end

function tsp1reg{S <: Real, T <: Real}(x::AbstractArray{S}, y::AbstractArray{T}, HD::Bool)
    order = sortperm( x )
    xsort = x[ order ]
    ysort = y[ order ]
    vec1  = outer( ysort, ysort, - )
    vec2  = outer( xsort, xsort, - )
    v1    = vec1[ vec2 .> 0 ]
    v2    = vec2[ vec2 .> 0 ]
    b1    = median!( v1./v2 )
    b0    = 0.0
    if !HD
        b0 = median( y ) - b1 * median( x )
    else
        b0 = hd( y ) - b1*hd( x )
    end
    return [b0, b1]
end

function tsreg_coef(mf::ModelFrame, HD::Bool, iter::Integer)
    y     = vector( model_response( mf ) )
    x     = ModelMatrix(mf).m[:,2:end]
    np, n = size(x, 2), size(x, 1)

    #ONE PREDICTOR
    if np == 1
        coef = tsp1reg( x[:], y, false )
    #    coef[1], coef[2] =  output.intercept, output.slope
    #    res = output.res
    else
    #MULTIPLE PREDICTORS
        coef_temp = zeros( np )
        for p = 1:np
            coef_temp[p] = tsp1reg( x[:,p], y, false )[2]
        end
        res = y - x*coef_temp
        if !HD
            b0 = median!( res )
        else
            b0 = hd( res )
        end
        #r = zeros( n )
        #coef_old = coef_temp[:]
        for i = 1:iter
            for p = 1:np
                r = y - x*coef_temp - b0 + coef_temp[p].*x[:,p]
                coef_temp[p] = tsp1reg(x[:,p], r, false)[2]
            end
            if !HD
                b0 = median!( y - x*coef_temp )
            else
                b0 = hd( y - x*coef_temp )
            end
            coef_old = coef_temp[:]
        end
        coef = [b0, coef_temp]
    end
    return coef
end



function tsreg(formula::Expr, dataframe::AbstractDataFrame; varfun::Function=pbvar, corfun::Function=pbcor, HD::Bool=false, iter::Integer=10)
#  Compute Theil-Sen regression estimator
#  Gauss-Seidel algorithm is used when there is more than one predictor
    mf = ModelFrame( formula, dataframe )

    #WILL ADD THE ABILITY TO DO MULTIPLE REGRESSION
    output = regOut()
    coef   = zeros(2)
    res    = zeros(n)
    if np == 1
        output = tsp1reg( vector( mf[2] ), mr )
        coef[1], coef[2] =  output.intercept, output.slope
    else
        stop("Only 1 predictor is allowed.")
    end

    res  = temp1.res
    yhat = y - res

    epow   = pbvar(yhat)/pbvar(y)
    if epow >= 1
        epow = sqrt( pbcor(yhat, y).estimate )
    end
    stre   = sqrt(epow)
    output = DataFrame()
    output["b0"] = temp1.intercept
    output["b1"] = temp1.slope
    output["strength of assoc."] = stre
    output["explanatory power"]  = epow

    if np == 1
        lm_coef = coef(lm( formula, dataframe ))
        p     = FramedPlot()
        pts   = Points( vector(mf[2]), mr , "type", "dot")
        s1    = Slope( lm_coef[2], (0, lm_coef[1]), "type", "solid", "color", "blue")
        s2    = Slope( temp1.slope, (0, temp1.intercept), "type", "solid", "color", "red")
        add(p, pts, s1, s2)
        Winston.display(p)
    end
    output
end
