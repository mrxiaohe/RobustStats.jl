using PyPlot

type indtOutput
    method
    dstat
    pval_d
    wstat
    pval_w
    flag
end

type indirectTestOutput
    nboot
    n
    regfit
    sobel
    bootest
    bootci
    bootse
    p
end

indirectTestOutput()=indirectTestOutput(nothing,
                                        nothing,
                                        nothing,
                                        nothing,
                                        nothing,
                                        nothing,
                                        nothing,
                                        nothing)


function show(io::IO, object::indtOutput)
    if object.method!=nothing
        println(io, "$(object.method)")
    end
    if object.flag!=2
            @printf("Kolmogorov-Smirnov test statistic:  % .6f\n", object.dstat)
        if object.pval_d >=1
            @printf("Kolmogorov-Smirnov p value:          1.0\n")
        elseif object.pval_d < 10e-16
            @printf("Kolmogorov-Smirnov p value:         < 10e-16\n")
        else
            @printf("Kolmogorov-Smirnov p value:         % .6f\n", object.pval_d)
        end
    end
    if object.flag!=1
            @printf("Cramer-von Mises test statistic:    % .6f\n", object.wstat)
        if object.pval_w >=1
            @printf("Cramer-von Mises p value:            1.0\n")
        elseif object.pval_w < 10e-16
            @printf("Cramer-von Mises p value:           < 10e-16\n")
        else
            @printf("Cramer-von Mises p value:           % .6f\n", object.pval_w)
        end
    end
end


function show(io::IO, object::indirectTestOutput)
    @printf("\nTESTS OF INDIRECT EFFECT\n\n")
    @printf("Sample size:  %d\n", object.n)
    @printf("Number of bootstrap samples: %d\n\n", object.nboot)
    @printf("DIRECT AND TOTAL EFFECTS\n")
    @printf("          % 8s  %9s    %8s   %8s\n", "Estimate", "Std.Error", "t value", "Pr(>|t|)")
    @printf("b(YX):   % 8.6f   %8.6f   % 8.6f   %8.6f\n",
            object.regfit[1,1], object.regfit[1,2], object.regfit[1,3], object.regfit[1,4])
    @printf("b(MX):   % 8.6f   %8.6f   % 8.6f   %8.6f\n",
            object.regfit[2,1], object.regfit[2,2], object.regfit[2,3], object.regfit[2,4])
    @printf("b(YM.X): % 8.6f   %8.6f   % 8.6f   %8.6f\n",
            object.regfit[3,1], object.regfit[3,2], object.regfit[3,3], object.regfit[3,4])
    @printf("b(YX.M): % 8.6f   %8.6f   % 8.6f   %8.6f\n\n",
            object.regfit[4,1], object.regfit[4,2], object.regfit[4,3], object.regfit[4,4])

    @printf("INDIRECT EFFECT AND SIGNIFICANCE USING NORMAL DISTRIBUTION\n")
    @printf("          % 8s  %9s    %8s   %8s   %8s  %8s \n",
            "Estimate", "Std.Error", "z score", "CI lo", "CI hi", "Pr(>|z|)")
    @printf("Sobel:   % 8.6f   %8.6f   % 8.6f  % 8.6f  % 8.6f  %8.6f\n\n",
            object.sobel[1,1], object.sobel[1,2], object.sobel[1,3],
            object.sobel[1,4], object.sobel[1,5], object.sobel[1,6])

    @printf("BOOTSTRAP RESULTS OF INDIRECT EFFECT\n")
    @printf("          % 8s  %9s    %8s   %8s   %8s \n",
            "Estimate", "Std.Error", "CI lo", "CI hi", "P value")
    @printf("Effect:  % 8.6f   %8.6f   % 8.6f  % 8.6f   %8.6f\n",
            object.bootest, object.bootse, object.bootci[1], object.bootci[2], object.p)
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
function stein2{S <: Real, T <: Real}(x1::AbstractArray{S}, x2::AbstractArray{T};
    mu0::Real=0.0, alpha::Real=0.05, method::Bool=true)
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



#Compute adaptive kernel density estimate for univariate data
function akerd{S <: Real}(x::AbstractArray{S}; hval::Real=NaN, aval::Real=0.5,
               fr::Real=0.8, pts=NaN,
               plotit=true, xlab="", ylab="", title="", color="black")
    if isnan(pts)
       pts = x[:]
    end
    pts  = sort!(pts)

    m = _estimate_dispersion(x)
    fhat  = rdplot(x, pts=pts, plotit=false, fr=fr)
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

function akerd_loop(x, pts, hval, alam)
    npts  = length(pts)
    dhat  = zeros(npts)
    sqrt5 = sqrt(5)
    n     = length(x)
    temp  = zeros(n)
    for i = 1:npts
        epan_alam_hval = 0.0
        for j = 1:n
            temp[j] = (pts[i]-x[j])/(hval*alam[j])
            if abs(temp[j]) < sqrt5
                epan_alam_hval += (0.75*(1-0.2*temp[j]*temp[j])/sqrt5)/alam[j]
            end
        end
        dhat[i]=epan_alam_hval/npts
    end
    return dhat./hval
end



#Expected frequency curve. fr controls amount of smoothing, theta is the azimuthal direction and
#phi the colatitude
function rdplot{S <: Real}(x::AbstractArray{S}; fr::Real=NaN, pts=NaN,
                           plotit=true, title="", xlab="", ylab="", color="black")
    if fr == NaN; fr = 0.8; end
    if pts == NaN; pts = x[:];end
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


function regts1(vstar, yhat, res, mflag, x, tr)
    const n=length(res)
    ystar=zeros(n)
    [ystar[i]=yhat+res[i].*vstar[i] for i=1:n]
    rval=zeros(n)
    ystar_bar=tmean(ystar, tr=tr)
    rval=zeros(n)
    [rval[i]= sum(ystar[mflag[:,i]]-ystar_bar) for i=1:n]
    return rval
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




function bootindirect{S <: Real, T <: Real, W <: Real}(x::Vector{S}, y::Vector{T}, m::Vector{W}, nboot::Integer)
    n = length(x)
    tempx  = hcat(rep(1.0, n), zeros(n))
    tempy  = zeros(n)
    tempm  = zeros(n)
    tempmx = hcat(rep(1.0, n), zeros(Real, n), zeros(Real, n))
    bvec   = zeros(nboot)
    randid = rand(1:n, n*nboot)
    for i = 1:nboot
        for j = 1:n
            tempx[j,2] = tempmx[j,2] = x[randid[(i-1)*n + j]]
            tempy[j]   = y[randid[(i-1)*n + j]]
            tempm[j]   = tempmx[j,3] = m[randid[(i-1)*n + j]]
        end
        bvec[i] = (inv(tempx'tempx)'*(tempx'tempm))[2]*(inv(tempmx'tempmx)'*(tempmx'tempy))[3]
    end
    bvec
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



function t1waycore(args...)
    x, grp, tr, method = args[1], args[2], args[3], args[4]
    g = unique(grp)
    J = length(g)
    h, w, xbar, n = zeros(J), zeros(J), zeros(J), zeros(J)

    for j = 1:J
        n[j] = sum(grp.==g[j])
        temp = x[grp.==g[j]]
        h[j] = n[j] - 2*floor(tr*n[j])
        w[j] = h[j]*(h[j] - 1)/((n[j] -1)*winvar(temp, tr=tr))
        xbar[j] = tmean(temp, tr=tr)
    end
    u = sum(w)
    xtil = sum(w .* xbar)/u
    A = sum(w .* (xbar - xtil).*(xbar - xtil))/(J - 1)
    B = 2*(J - 2)*sum((1 - w./u).*(1 - w./u)./(h - 1))/(J*J - 1)
    statistic = A/(B + 1)
    nu1 = J - 1
    nu2 = 1/(3 * B/2/(J-2))
    p   = 1 - Rmath.pf(statistic, nu1, nu2)
    if method
        METHOD="Heteroscedastic one-way ANOVA for trimmed means\nusing a generalization of Welch's method.\n"
    else
        METHOD=nothing
    end
    output    = testOutput()
    output.n  = n
    output.df = [nu1, nu2]
    output.statistic = statistic
    output.p  = p
    output.method = METHOD
    output
end


#  Compute a (1-α) confidence interval for the trimmed mean
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
                             nullvalue::Real=0, seed=2, method::Bool=true)
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
        METHOD="Bootstrap (1-α) confidence interval for the trimmed mean\nusing a bootstrap percentile t method\n"
    elseif method && !side
        METHOD="Bootstrap (1-α) confidence interval for the trimmed mean\nusing a bootstrap percentile t method\n[NOTE: p value is computed only when side=true]\n"
    else
        METHOD=nothing
    end
    if !side
        ci[1]=tmean(x, tr=tr)-tempout[itop]*trimse(x, tr=tr)
        ci[2]=tmean(x, tr=tr)-tempout[ibot]*trimse(x, tr=tr)
        output = testOutput()
        output.method = METHOD
        output.estimate = tmean(x, tr=tr)
        output.ci = ci
        output.statistic = test
        return output
    else
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


function trimcibt_loop(x, n, nboot, tr, side, randid, test)
    xbar=tmean(x, tr=tr)
    temp=zeros(n)
    tval=zeros(nboot)
    if side
        test=abs(test)
        pval=0.0
        for i=1:(nboot*n)
            if (i%n)!=0
                temp[i%n]=x[randid[i]]-xbar
            else
                temp[n]=x[randid[i]]-xbar
                tval[div(i, n)]=abs(tmean(temp, tr=tr)./trimse(temp, tr=tr))
                pval += tval[div(i, n)]>=test? 1/nboot :0.0
            end
        end
        return sort!(tval), pval
    else
        for i=1:(nboot*n)
            if (i%n)!=0
                temp[i%n]=x[randid[i]]-xbar
            else
                temp[n]=x[randid[i]]-xbar
                tval[div(i, n)]=tmean(temp, tr=tr)./trimse(temp, tr=tr)
            end
        end
        return sort!(tval)
    end
end
