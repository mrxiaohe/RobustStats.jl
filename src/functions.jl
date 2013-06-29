#Trimmed mean; amount of trimming defaults to 20% or 0.2
function tmean{S <: Real}(x::Vector{S}; tr::Real=0.2)
    if tr < 0 || tr > 0.5
        error("tr cannot be smaller than 0 or larger than 0.5")
    elseif tr == 0
        return mean(x)
    elseif tr == .5
        return median(x, checknan=false)
    else 
        n   = length(x)
        lo  = floor(n*tr)+1
        hi  = n+1-lo
        return mean(sort!(copy(x))[lo:hi])
    end
end
#tmean{S <: Real}(x::DataVector{S}; tr::Real=0.2)=tmean(removeNA(x), tr=tr)

function tmean!{S <: Real}(x::Vector{S}; tr::Real=0.2)
    if tr < 0 || tr > 0.5
        error("tr cannot be smaller than 0 or larger than 0.5")
    elseif tr == 0
        return mean(x)
    elseif tr == .5
        return median!(x, checknan=false)
    else 
        n   = length(x)
        lo  = floor(n*tr)+1
        hi  = n+1-lo
        return mean(sort!(x)[lo:hi])
    end
end


#Winsorize data
function winval{S <: Real}(x::Vector{S}; tr::Real=0.2)
    xcopy   = copy(x)
    n       = length(x)
    xcopy   = sort!(xcopy)
    ibot    = floor(tr*n)+1
    itop    = n-ibot+1
    xbot, xtop = xcopy[ibot], xcopy[itop]
    return  [x[i]<=xbot?xbot:(x[i]>=xtop?xtop:x[i]) for i=1:n]
end
winval{S <: Real}(x::DataVector{S}; tr::Real=0.2)=winval(removeNA(x), tr=tr)

#Winsorized mean
winmean{S <: Real}(x::Vector{S}; tr=0.2)=mean(winval(x, tr=tr))
#winmean{S <: Real}(x::DataVector{S}; tr::Real=0.2)=mean(winval(removeNA(x), tr=tr))

#Winsorized variance
winvar{S <: Real}(x::Vector{S}; tr=0.2)=var(winval(x, tr=tr))
#winvar{S <: Real}(x::DataVector{S}; tr::Real=0.2)=var(winval(remove(x), tr=tr))

#winsorized standard deviation
winstd{S <: Real}(x::Vector{S}; tr=0.2)=std(winval(x, tr=tr))
#winstd{S <: Real}(x::DataVector{S}; tr::Real=0.2)=std(winval(remove(x), tr=tr))


#Estimate the standard error of the gamma trimmed mean
trimse{S <: Real}(x::Vector{S}; tr::Real=0.2)=
    sqrt(winvar(x,tr=tr))/((1-2tr)*sqrt(length(x)))
#trimse{S <: Real}(x::DataVector{S}; tr::Real=0.2)=trimse(removeNA(x), tr=tr)

#Compute a 1-alpha confidence interval for the trimmed mean
function trimci{S <: Real}(x::Vector{S}; tr::Real=0.2, alpha::Real=0.05, nullvalue::Real=0, method=true)
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
#trimci{S <: Real}(x::DataVector{S}; tr=0.2, alpha=0.05, nullvalue=0, method=true)=
    trimci(removeNA(x), tr=tr, alpha=alpha, nullvalue=nullvalue, method=method)


#Stein's method
function stein1{S <: Real}(x::Vector{S}, delta::Real; alpha::Real=0.05, pow::Real=0.8, oneside::Bool=false, n=nothing, variance=nothing)
    delta    = abs(delta)
    n        = n==nothing?length(x):n 
    variance = variance==nothing?var(x):variance 
    df       = n-1
    alpha    = !oneside?alpha/2:alpha
    d        = (delta/(Rmath.qt(pow, df)-Rmath.qt(alpha, df)))*(delta/(Rmath.qt(pow, df)-Rmath.qt(alpha, df)))
    N        = max([n, floor(variance/d)+1])
    return integer(N)
end
#stein1{S <: Real}(x::DataVector{S}, delta::Real; alpha::Real=0.05, pow::Real=0.8, oneside::Bool=false, n=nothing, variance=nothing)=
 #   stein1(removeNA(x), delta, alpha=alpha, pow=pow, oneside=oneside, n=n, variance=variance)

    
function stein2{S <: Real, T <: Real}(x1::Vector{S}, x2::Vector{T}; mu0::Real=0.0, alpha::Real=0.05, method::Bool=true)
    n       = length(x1)
    df      = n-1
    N       = n+length(x2)
    xbar    = mean([x1, x2])
    test    = sqrt(N)*(xbar-mu0)/std(x1)
    crit    = Rmath.qt(1-alpha/2, df)
    lo, hi  = xbar-crit*std(x1), xbar+crit*std(x1)
    sig     = 2*(1-Rmath.pt(test, df))
    METHOD  = method?"The 2nd stage of Stein's method":nothing
    output  =testOutput()
    output.method    = METHOD
    output.df        = df
    output.estimate  = xbar
    output.ci        = [lo, hi]
    output.statistic = test
    output.crit      = crit
    output.p         = sig
    return output
end
#stein2{S <: Real, T <: Real}(x1::DataVector{S}, x2::DataVector{T}; mu0::Real=0.0, alpha::Real=0.05, method::Bool=true)=
 #   stein2(removeNA(x1), removeNA(x2), mu0=mu0, alpha=alpha, method=method)


#Ideal fourths
function idealf{S <: Real}(x::Vector{S}; method::Bool=true)
    n       = length(x)
    j       = integer(floor(n/4+5/12))
    y       = sort(x)
    g       = n/4-j+5/12
    k       = n-j+1
    METHOD  = method?"The ideal fourths":nothing
    idealfOutput((1-g).*y[j]+g.*y[j+1], (1-g).*y[k]+g.*y[k-1], METHOD)
end
#idealf{S <: Real}(x::DataVector{S}; method::Bool=true)=idealf(removeNA(x), method=method)

#Percentage bend midvariance
function pbvar{S <: Real}(x::Vector{S}; beta::Real=0.2)
    output  = 0.0
    const n = length(x)
    w       = zeros( n)
    med     = median(x, checknan=false)
    w       = [abs(x[i]-med) for i=1:n]
    w       = sort!(w)
    m       = floor((1-beta).*n+0.5)
    z       = 0.0
    omega   = w[m]
    counter = 0
    if omega > 0
        for i = 1:n
            y = (x[i]-med)/omega
            z += y>1.0?1.0:(y< -1.0?1.0:y*y)
            counter+=abs(y) < 1?1:0
        end
        output = n*(omega*omega)*z/(counter*counter)
    end
    return output
end
#pbvar{S <: Real}(x::DataVector{S}; beta::Real=0.2)=pbvar(removeNA(x), beta=beta)


#Biweight midvariance
function bivar{S <: Real}(x::Vector{S})
    m   = median(x, checknan=false)
    n   = length(x)
    u   = 0.0
    MAD = mad(x)
    av  = zeros( n)
    top = bot=0.0
    q   = Rmath.qnorm(0.75)
    for i = 1:n
        u     = abs(x[i]-m)./(9.*q.*MAD)
        av[i] = u<1.0?1:0
        top   += (n*av[i]*(x[i]-m)*(x[i]-m)*(1-u*u).^4)
        bot   += av[i]*(1-u*u)*(1-5*u*u)
    end
    top/(bot*bot)
end
#bivar{S <: Real}(x::DataVector{S})=bivar(removeNA(x))


#Compute the tau measure of location as described in Yohai and Zamar (JASA, 1988, 83, 406-413)
function tauloc{S <: Real}(x::Vector{S}; cval::Real=4.5)
    s     = Rmath.qnorm(0.75)*mad(x)
    med   = median(x, checknan=false)
    n     = length(x)
    Wnom  = Wden = 0.0
    for i in 1:n
        y = (x[i]-med)/s
        temp = (1.0-y*y/cval/cval)*(1.0-y*y/cval/cval)
        if abs(temp) <= cval
            Wnom += temp*x[i]
            Wden += temp
        end
    end
    return Wnom/Wden
end
#tauloc{S <: Real}(x::DataVector{S}; cval::Real=4.5)=tauloc(removeNA(x), cval=cval)


#Compute the tau measure of scale as described in Yohai and Zamar (JASA, 1988, 83, 406-413)
function tauvar{S <: Real}(x::Vector{S}; cval::Real=3.0)
    x     = x[:]
    s     = Rmath.qnorm(0.75)*mad(x)
    tloc  = tauloc(x)
    n     = length(x)
    W     = 0.0
    cval2 = cval*cval
    [W    += min(((x[i]-tloc)/s)*((x[i]-tloc)/s), cval2) for i=1:n]
    s*s*W/n
end
#tauvar{S <: Real}(x::DataVector{S}; cval::Real=3.0)=tauvar(removeNA(x), cval=cval)


#Outlier detection method using the boxplot rule. The ideal fourths are
#used to estimate the quartiles
function outbox{S <: Real}(x::Vector{S}; mbox::Bool=false, gval::Real=NaN, method::Bool=true)
    n    = length(x)
    temp = idealf(x, method=false)
    cl   = 0.0
    cu   = 0.0
    if mbox
        if isnan(gval)
            gval=(17.63*n-23.64)/(7.74*n-3.71)
        end
        cl = median(x, checknan=false) - gval*(temp.upper_quartile-temp.lower_quartile)
        cu = median(x, checknan=false) + gval*(temp.upper_quartile-temp.lower_quartile)
    elseif !mbox
        if isnan(gval)
            gval=1.5
        end
        cl = temp.lower_quartile - gval*(temp.upper_quartile-temp.lower_quartile)
        cu = temp.upper_quartile + gval*(temp.upper_quartile-temp.lower_quartile)
    end
    flag=falses(n)
    for i=1:n
        if (x[i].<cl) || (x[i].>cu)
            flag[i]=true
        end
    end
    vec=1:n
    outid  = sum(flag)>0?vec[flag]:[]
    keepid = vec[!flag]
    outval = x[flag]
    nout   = length(outid)
    if method && !mbox
        METHOD = "Outlier detection method using \nthe ideal-fourths based boxplot rule\n"
    elseif method && mbox
        METHOD = "Outlier detection method using \nthe ideal-fourths based boxplot rule\n(using the modification suggested by Carling (2000))\n"
    else 
        METHOD = nothing
    end
    outOutput(outid, keepid, outval, nout, METHOD)
end
#outbox{S <: Real}(x::DataVector{S}; mbox=false, gval::Real=NaN, method=true)=outbox(removeNA(x), mbox=mbox, gval=gval, method=method)

#Compute adaptive kernel density estimate for univariate data
function akerd{S <: Real}(x::Vector{S}; hval::Real=NaN, aval::Real=0.5, op::Integer=1, 
               fr::Real=0.8, pyhat::Bool=false, pts=NaN, plotit=true, xlab="", 
               ylab="", title="", color="black", plottype="solid")
    xsort =  sort(x)
    if op == 1
        m =  mad(x)
        if m == 0.0
            temp = idealf(x, method=false)
            m = (temp.upper_quartile-temp.lower_quartile)/(Rmath.qnorm(0.75)-Rmath.qnorm(0.25))
        end
        if m  == 0.0
            m =  sqrt(winvar(x)./0.4129)
        end
        if m  == 0.0
            error("All measures of dispersion are equal to 0")
        end
        fhat  = rdplot(x, pyhat=true, plotit=false, fr=fr, pyhat=true)
        if m  >  0.0
            [fhat[i] = fhat[i]/(2.0*fr*m) for i=1:length(fhat)]
        end
    end
    if op == 2
        init = kde(x)
        fhat = init.density
        x    = init.x
    end
    n = length(x)
    if isnan(hval)
        sig  = std(x)
        temp = idealf(x, method=false)
        iqr  = (temp.upper_quartile-temp.lower_quartile)/1.34
        A    = min(sig, iqr)
        A    = A==0.0? winstd(x)/0.64:A
        hval = 1.06*A/n^0.2
    end
    gm = 0.0
    gm_int = 0
    nfhat = length(fhat)
    for i = 1:nfhat
        if fhat[i] > 0.0
            gm += log(fhat[i])
            gm_int += 1
        end
    end
    gm = exp(gm/gm_int)
    alam = zeros(n)
    [alam[i] = (fhat[i]/gm)^(0-aval) for i=1:nfhat]
    if isnan(pts)
        pts = x[:]
    end
    pts  = sort!(pts)
    dhat = akerd_loop(x, pts, hval, alam)
    if plotit
        p = FramedPlot()
        add(p, Curve(pts, dhat, "color", color, "type", plottype))
        if title != nothing
            setattr(p, "title", title)
        end
        if xlab != nothing
            setattr(p, "xlabel", xlab)
        end
        if ylab != nothing
            setattr(p, "ylabel", ylab)
        end
        Winston.tk(p)
    end
    if pyhat
        return dhat
    end
end
#akerd{S <: Real}(x::DataVector{S}; hval::Real=NaN, aval::Real=0.5, 
 #                op::Integer=1, fr::Real=0.8, pyhat::Bool=false, 
  #               pts=NaN, plotit=true, xlab="", ylab="", title="", 
   #              color="black", plottype="solid")=akerd(removeNA(x), 
    #             hval=hval, aval=aval, op=op, fr=fr, pyhat=pyhat, pts=pts, 
     #            plotit=plotit, xlab=xlab, ylab=ylab, title=title, color=color, 
      #           plottype=plottype)

#Expected frequency curve. fr controls amount of smoothing, theta is the azimuthal direction and
#phi the colatitude
function rdplot{S <: Real}(x::Vector{S}; fr::Real=NaN, plotit=true, 
                           theta=50, phi=25, expand=0.5, pyhat=false, pts=NaN,
                           title="", xlab="", ylab="", color="black", plottype="solid")
    x = x[:]
    if isnan(fr)
        fr = 0.8
    end
    if isnan(pts)
        pts = x[:]
    end
    nx   = length(x)
    npts = length(pts)
    rmd  = [sum(near(x, pts[i], fr))*1.0 for i=1:npts]
    MAD  = mad!(x)
    if MAD != 0.0
        [rmd[i] = rmd[i]/(2*fr*MAD) for i=1:npts]
    end
    [rmd[i] = rmd[i]/nx i=1:npts]
    if plotit
        index = sortperm(pts);
        p     = FramedPlot()
        add(p, Curve(pts[index], rmd[index], "color", color, "type", plottype))
        if title != nothing
            setattr(p, "title", title)
        end
        if xlab != nothing
            setattr(p, "xlabel", xlab)
        end
        if ylab != nothing
            setattr(p, "ylabel", ylab)
        end
        Winston.tk(p)
    end
    if pyhat 
        return rmd
    end
end
#rdplot{S <: Real}(x::DataVector{S}; fr::Real=NaN, plotit=true, theta=50, phi=25, expand=0.5, 
 #                 pyhat=false, pts=NaN, title="", xlab="", ylab="", color="black", plottype="solid")=
#rdplot(removeNA(x), fr=fr, plotit=plotit, theta=theta, phi=phi, expand=expand, 
 #   pyhat=pyhat, pts=pts, title=title, xlab=xlab, ylab=ylab, color=color, plottype=plottype)


#Determine which values in x are near pt
function near{S <: Real}(x::Vector{S}, pt::Real, fr::Real=1.0)
    m = mad!(x)
    if m == 0.0
        temp = idealf(x, method=false)
        m = (temp.upper_quartile-temp.lower_quartile)/(Rmath.qnorm(0.75)-Rmath.qnorm(.25))
    end
    if m == 0.0 
        m = sqrt(winvar(x)/0.4129)
    end
    if m == 0.0
        error("All measures of dispersion are equal to 0")
    end
    n = length(x)
    dflag = falses(n)
    fr_m = fr*m
    #[dflag[i]=abs(x[i]-pt) <= fr_m ?true:false for i=1:n]
    for i = 1:n
        if abs(x[i]-pt) <= fr_m
            dflag[i] = true
        end
    end

    return dflag
end
#near{S <: Real}(x::DataVector{S}, pt::Real, fr::Real=1.0)=near(removeNA(x), pt, fr)


#Computing the confidence interval for the median using the 
#Hettmansperger and Sheather method.
function sint{S <: Real}(x::Vector{S}; alpha::Real=0.05, method::Bool=true)
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
#sint{S <: Real}(x::DataVector{S}; alpha::Float64=0.05, method=true)=sint(removeNA(x), alpha=alpha, method=method)


#Computing standard error of the median using a method recommended by McKean and Shrader (1984)
function msmedse{S <: Real}(x::Vector{S})
    if duplicated(x)
        warn("Tied values detected. Estimate of standard error might be highly inaccurate, even with n large")
    end
    y = sort(x)
    n = length(y)
    av::Int = round((n+1)./2-Rmath.qnorm(.995).*sqrt(n/4))
    if av == 0
        av = 1
    end
    top::Int = n-av+1
    abs((y[top]-y[av])/(2*Rmath.qnorm(.995)))
end
#msmedse{S <: Real}(x::DataVector{S})=msmedse(removeNA(x))


#Compute a 1-alpha confidence interval for p, the probability of success for a binomial dist. using Pratt's method
#y is a vector of 1's and 0's, x is the number of successes observed among n trials.

function binomci(x::Int, n::Int; alpha::Real=0.05)
    if x > n
        error("x must be smaller than or equal to n")
    elseif x < 0
        error("x cannot be negative")
    elseif n == 1
        error("Something is wrong: number of observations is only 1")
    end
    if x != n && x != 0
        z     = Rmath.qnorm(1-alpha/2)
        A     = ((x+1)/(n-x)).*((x+1)/(n-x))
        B     = 81.*(x+1).*(n-x)-9.*n-8
        C     = (0-3).*z.*sqrt(9.*(x+1).*(n-x).*(9*n+5-z.*z)+n+1)
        D     = 81.*(x+1).*(x+1)-9.*(x+1).*(2+z*z)+1
        E     = 1+A.*((B+C)./D).*((B+C)./D).*((B+C)./D)
        upper = 1/E
        A     = (x./(n-x-1)).*(x./(n-x-1))
        B     = 81.*x.*(n-x-1)-9.*n-8
        C     = 3.*z.*sqrt(9.*x.*(n-x-1).*(9.*n+5-z.*z)+n+1)
        D     = 81.*x.*x-9.*x.*(2+z.*z)+1
        E     = 1+A.*((B+C)./D).^3
        lower = 1/E
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
    phat=x/n
    binomciOutput(phat, [lower, upper], n)
end
binomci(x::Vector{Int}; alpha::Real=0.05)=
sum([x[i]!=0 && x[i]!=1 for i=1:length(x)])>0?error("Data must be 0's and 1's"):binomci(length(x.==1), length(x), alpha=alpha)
#binomci(x::DataVector{Int}; alpha::Real=0.05)=
#sum([x[i]!=0 && x[i]!=1 for i=1:length(x)])>0?error("Data must be 0's or 1's"):binomci(removeNA(x), alpha=alpha)


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
#acbinomci(x::DataVector{Int}; alpha::Real=0.05)=
#sum([x[i]!=0 && x[i]!=1 for i=1:length(x)])>0?error("Data must be 0's or 1's"):acbinomci(removeNA(x), alpha=alpha)


#Extension of Stein's method based on the trimmed mean.
function stein1_tr{S <: Real}(x::Vector{S}, del::Real; alpha::Real=0.05, pow::Real=0.8, tr::Real=0.2)
    if tr <0 || tr >= .5
        error("Argument tr must be between 0 and .5")
    end
    n =length(x)
    g::Int    = floor(tr*n)
    df::Int   = n-2*g-1
    t1        = Rmath.qt(pow, df)
    t2        = Rmath.qt(alpha./2, df)
    dv        = (del/(t1-t2))*(del/(t1-t2))
    nvec::Int = floor(trimse(x, tr=tr)./dv)+1
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
        g::Integer = floor(tr*n)
        df    = n-2g-1
        tdiff = Rmath.qt(pow, df)-Rmath.qt(alpha/2.0/ntest, df)
        dv    = (del/tdiff)*(del/tdiff)
        N     = zeros(Int, ntest)
        if ntest > 1
            for i = 1:J
                for j = 1:J
                    if i < j
                        N[ic+=1] = floor(trimse(x[:,i]-x[:,j], tr=tr)/dv)+1
                    end
                end
            end
            return max(vcat(N, n))
        elseif ntest == 1
            return int(floor(trimse(x[:,1], tr=tr)/dv)+1)
        end
    end
end
#stein1_tr{S <: Real}(x::DataArray{S, 2}, del::Real; alpha::Real=0.05, pow::Real=0.8, tr::Real=0.2)=
 #   stein1_tr(deleterows!(x), del, alpha=alpha, pow=pow, tr=tr)



#Extension of the second stage of Stein's method when performing all pairwise comparisons
#among J dependent groups.
function stein2_tr{S <: Real, T <: Real}(x::Vector{S}, y::Vector{T}; alpha::Real=0.05, tr::Real=0.2, method::Bool=true)
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

function stein2_tr{S <: Real, T <: Real}(x::Array{S}, y::Array{T}; alpha::Real=0.05, tr::Real=.2, method::Bool=true)
    if tr <0 || tr >= .5
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
        test=DataFrame({Integer, Integer, Float64}, ["grp1", "grp2", "test_stat"], ntest)
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

#stein2_tr{S <: Real, T <: Real}(x::DataVector{S}, y::DataVector{S}; alpha::Real=0.05, tr::Real=0.2, method::Bool=true)=
#    stein2_tr(removeNA(x), removeNA(y), alpha=alpha, tr=tr, method=method)
#stein2_tr{S <: Real, T <: Real}(x::Array{Int}, y::Array{Int}; alpha=0.05, tr=0.2, method=true)=
#    stein2_tr(convert(Array{Float64}, x), convert(Array{Float64}, y), alpha=alpha, tr=tr, method=method)


#Compute 1-alpha confidence interval for the median using the Hettmansperger-Sheather interpolation method
function sintv2{S <: Real}(x::Vector{S}; alpha::Real=0.05, nullval::Real=0.0, method::Bool=true)
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

#sintv2{S <: Real}(x::DataVector{S}; alpha::Real=0.05, nullval::Real=0.0, method::Bool=true)=
 #   sintv2(removeNA(x), alpha=alpha, nullval=nullval, method=method)

#Evaluate Huber's Psi function for each value in the vector x.
#The bending constant defaults to 1.28
hpsi{S <: Real}(x::Vector{S}, bend::Real=1.28)=[abs(x[i])<=bend?x[i]:bend*sign(x[i]) for i=1:length(x)]
#hpsi{S <: Real}(x::DataVector{S}, bend::Real=1.28)= hpsi(removeNA(x), bend=bend)


#Compute one-step M-estimator of location using Huber's Psi.
#The default bending constant is 1.28
function onestep{S <: Real}(x::Vector{S}, bend::Real=1.28)
    MED = median(x, checknan=false)
    MAD = mad(x)
    n   = length(x)
    y   = [(x[i]-MED)/MAD for i = 1:n]
    A   = sum(hpsi(y, bend))
    B   = 0
    [B  += abs(y[i]) <= bend ? 1:0 for i=1:n]
    return MED + MAD*A/B
end
#onestep{S <: Real}(x::DataVector{S}, bend::Real=1.28)=onestep(removeNA(x), bend=bend)

#Compute a bootstrap, .95 confidence interval for the
#measure of location corresponding to the argument est.
#By default, a one-step
#M-estimator of location based on Huber's Psi is used.
#The default number of bootstrap samples is nboot=500
#
#nv=null value when  computing a p-value

function onesampb{S <: Real}(x::Vector{S}; est::Function=onestep, alpha::Real=0.05, nboot::Integer=2000, seed::Union(Int, Bool)=2, nv::Real=0)
    if iseltype(seed, Int)
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
#onesampb{S <: Real}(x::DataVector{S}; est::Function=onestep, alpha::Real=0.05, nboot::Integer=2000, seed::Union(Int, Bool)=2, nv::Real=0)=
 #   onesampb(removeNA(x), est=est, alpha=alpha, nboot=nboot, seed=seed, nv=nv)



#MOM-estimator of location. The default bending constant is 2.24
function mom{S <: Real}(x::Vector{S}; bend::Real=2.24)
    n    = length(x)
    flag = trues(n)
    MED  = median(x, checknan=false)
    MAD  = mad(x)
    for i = 1:n
        if x[i] > MED+bend*MAD || x[i] < MED-bend*MAD
            flag[i]=false
        end
    end
    mean(x[flag])
end
#mom{S <: Real}(x::DataVector{S}; bend=2.24)=mom(removeNA(x), bend=bend)

function mom!{S <: Real}(x::Vector{S}; bend::Real=2.24)
    n    = length(x)
    flag = trues(n)
    MED  = median!(x, checknan=false)
    MAD  = mad!(x)
    for i = 1:n
        if x[i] > MED+bend*MAD || x[i] < MED-bend*MAD
            flag[i]=false
        end
    end
    mean(x[flag])
end



#Compute a bootstrap, .95 confidence interval for the
#MOM-estimator of location based on Huber's Psi.
#The default number of bootstrap samples is nboot=500
function momci{S <: Real}(x::Vector{S}; bend::Real=2.24, alpha::Real=0.05, nboot::Integer=2000, seed::Union(Int, Bool)=2, method=true)
    if iseltype(seed, Int)
        srand(seed)
    elseif seed
        srand(2)
    end
    n    = length(x)
    bvec = zeros(nboot)
    temp = zeros(n)
    randid=rand(1:n, n*nboot)
    for i = 1:nboot
        [temp[j] = x[randid[(i-1)*n+j]] for j=1:n]
        bvec[i]=mom!(temp, bend=bend)
    end
    low  = round((alpha./2).*nboot) + 1
    up   = nboot-low + 1
    bvec = sort!(bvec)
    METHOD = method? "Bootstrap .95 confidence interval for the MOM-\nestimator of location based on Huber's Psi\n":nothing
    output = testOutput()
    output.method = METHOD
    output.ci     = [bvec[low], bvec[up]]
    output
end
#momci{S <: Real}(x::DataVector{S}; bend::Real=2.24, alpha::Real=0.05, nboot::Integer=2000, seed::Union(Int, Bool)=2, method=true)=
 #   momci(removeNA(x), bend=bend, alpha=alpha, nboot=nboot, seed=seed, method=method)


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

function trimpb{S <: Real}(x::Vector{S}; tr::Real=0.2, alpha::Real=0.05, nboot::Integer=2000, 
                win::Union(Bool, Real)=false, plotit::Bool=false, pop::Int=1, nullval::Real=0.0, 
                xlab="X", ylab="Density", fr::Union(Nothing, Real)=nothing, seed::Union(Bool, Integer)=2, 
                method::Bool=true)
    if iseltype(win, Bool)
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
    if iseltype(seed, Bool)
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

#trimpb{S <: Real}(x::DataVector{S}; tr::Real=0.2, alpha::Real=0.05, nboot::Integer=2000, 
 #               win::Union(Bool, Real)=false, plotit::Bool=false, pop::Int=1, nullval::Real=0.0, 
  #              xlab="X", ylab="Density", fr::Union(Nothing, Real)=nothing, seed::Union(Bool, Integer)=2, 
   #             method::Bool=true)=trimpb(removeNA(x), tr=tr, alpha=alpha, nboot=nboot, win=win, plotit=plotit, 
    #            pop=pop, nullval=nullval, xlab=xlab, ylab=ylab, fr=fr, seed=seed, method=method)



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

function trimcibt{S <: Real}(x::Vector{S}; tr::Real=0.2, alpha::Real=0.05, nboot::Integer=2000, side::Bool=true, 
                             plotit::Bool=false, op::Integer=1, nullval::Real=0, seed::Union(Bool, Int)=2, method::Bool=true)
    if iseltype(seed, Bool)
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
#trimcibt{S <: Real}(x::DataVector{S}; tr::Real=0.2, alpha::Real=0.05, nboot::Integer=2000, side::Bool=true, 
 #                   plotit::Bool=false, op::Integer=1, nullval::Real=0, seed::Union(Bool, Int)=2, method::Bool=true) = 
  #                  trimcibt(removeNA(x), tr=tr, alpha=alpha, nboot=nboot, side=side, plotit=plotit, op=op, nullval=nullval, 
   #                 seed=seed, method=method)

#   Compute bootstrap estimate of the standard error of the
#   estimator est
#   The default number of bootstrap samples is nboot=1000

function bootse{S <: Real}(x::Vector{S}; nboot::Integer=1000, est::Function=median, seed::Union(Bool, Integer)=2)
   if iseltype(seed, Bool)
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

#bootse{S <: Real}(x::DataVector{S}; nboot::Integer=1000, 
 #                 est::Function=median, seed::Union(Bool, Integer)=2)=
  #                bootse(removeNA(x), nboot=nboot, est=est, seed=seed)


#   Compute a .95 confidence interval for Pearson's correlation coefficient.
#
#   This function uses an adjusted percentile bootstrap method that
#   gives good results when the error term is heteroscedastic.
function pcorb{S <: Real, T <: Real}(x::Vector{S}, y::Vector{T}; seed::Union(Bool, Integer)=2, plotit::Bool=false)
   if iseltype(seed, Bool)
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
#pcorb{S <: Real, T <: Real}(x::DataVector{S}, y::DataVector{T}; nboot::Integer=2000, seed::Union(Bool, Integer)=2)
#= pcorb(removeNAVector{Int}, y::Vector{Float64}; nboot=2000, seed::Union(Bool, Int)=2)=
    


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
              tr::Float64=0.0, seed::Union(Bool, Int)=2, method=true)
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
    if iseltype(seed, Bool)
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


function indirectTest{S <: Real, T <: Real, W <: Real}(dv::Vector{S}, iv::Vector{T}, m::Vector{W}; 
            nboot::Integer=5000, alpha::Real=0.05, scale::Bool=false, seed::Union(Bool, Int)=2, plotit::Bool=false)
    if iseltype(seed, Bool)
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
function wincor{S <: Real, T <: Real}(x::Vector{S}, y::Vector{T}; tr::Real=0.2)
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
function yuend{S <: Real, T <: Real}(x::Vector{S}, y::Vector{T}; tr::Real=0.2, alpha::Real=0.05, method::Bool=true)
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




# Data=[-1.0915999  0.7010442  0.2540099 -1.1499019  1.5375046  1.5362423 -1.2846310  1.1551006 -0.1107871  0.1731488  0.1824683 -0.9044951  0.9024446 -1.1056716 -0.3281121]











