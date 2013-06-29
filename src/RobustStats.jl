using DataFrames
using Rmath
using Distributions
using Stats
using GLM
using Winston

module robust
using DataFrames
using Rmath
using Winston
using Distributions
using GLM
using Stats
import Base.show
import DataFrames.complete_cases
import DataFrames.deleterows!


export 
    outOutput,
    idealfOutput,
    testOutput,
    tmean,
    winval,
    winmean,
    winvar,
    trimse,
    trimci,
    stein1,
    stein2,
    idealf,
    pbvar,
    bivar,
    tauloc,
    tauvar,
    outbox,
    akerd,
    akerd_C,
    rdplot,
    sint,
    binomci,
    acbinomci,
    near, 
    msmedse,
    stein1_tr,
    stein2_tr,
    sintv2,
    seq,
    cnorm,
    hpsi,
    onestep,
    onesampb,
    mom,
    momci,
    trimpb,
    trimcibt,
    indt,
    pcorb,
    bootse,
    indirectTest, 
    yuend
    

type outOutput
    outid
    keepid
    outval
    nout
    method
end

type indtOutput
    method
    dstat
    pval_d
    wstat
    pval_w
    flag
end

type idealfOutput
    lower_quartile
    upper_quartile
    method
end

type testOutput
    method
    h0
    n
    df
    estimate
    se
    statistic
    crit
    ci
    p
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

testOutput()=testOutput(
                        nothing,
                        nothing,
                        nothing,
                        nothing,
                        nothing,
                        nothing,
                        nothing,
                        nothing,
                        nothing,
                        nothing)

function show(io::IO, object::testOutput)
    if object.method != nothing
        @printf("%s\n", object.method)
    end
    if object.h0 != nothing
        @printf(" Alternative hypothesis:   % s\n", object.h0)
    end 
    if object.n != nothing
        @printf(" Sample size:         % d\n", object.n)
    end
    if object.df != nothing
        if iseltype(object.df, Integer)
        @printf(" Degrees of freedom:   %d\n", object.df)
        elseif iseltype(object.df, FloatingPoint)
        @printf(" Degrees of freedom:  % .2f\n", object.df)
        end
    end
    if object.estimate != nothing
        @printf(" Estimate:            % .6f\n", object.estimate)
    end
    if object.se != nothing
        @printf(" Standard error:      % .6f\n", object.se)
    end
    if object.statistic != nothing
        if ndims(object.statistic)==0
            @printf(" Statistic:           % .6f\n", object.statistic)
        else ndims(object.statistic)==2
            @printf(" Statistic:            % 7s  % 7s  % 9s\n", "Group 1", "Group 2", "Statistic")
            for i = 1:size(object.statistic, 1)
                @printf("                       % 7d  % 7d  %  7.6f\n", 
                        object.statistic[i,1], object.statistic[i,2], object.statistic[i,3])
            end
        end
    end
    if object.crit != nothing
        @printf(" Critical value:      % .6f\n", object.crit)
    end
    if object.ci != nothing
        if iseltype(object.ci, Integer)
            @printf(" Confidence interval: % d      % d\n", object.ci[1], object.ci[2])
        else 
            @printf(" Confidence interval: % .6f      % .6f\n", object.ci[1], object.ci[2])
        end
    end 
    if object.p != nothing
        if object.p < 10e-16
            @printf(" p value:             < 10e-16\n")
        elseif object.p < 1 
            @printf(" p value:             % .6f\n", object.p)
        else 
            @printf(" p value:              1.0")
        end      
    end
end    

type binomciOutput
    phat::FloatingPoint
    confint::Array
    n::Int
end

function show(io::IO, object::outOutput)
    if object.method!=nothing
        println(io, "$(object.method)")
    end
    @printf("Outlier ID:         ")
    if length(object.outid) == 0 
        @printf("\n")
    else 
        for i = 1:length(object.outid)
            if i < length(object.outid)
                @printf("%d, ", object.outid[i])
            else 
                @printf("%d\n", object.outid[i])
            end
        end 
    end 
    @printf("Outlier value:      ")
    if length(object.outid) == 0 
        @printf("\n")
    else 
        if iseltype(object.outval, Integer)
            for i = 1:length(object.outval)
                if i < length(object.outval)
                    @printf("%d, ", object.outval[i])
                else 
                    @printf("%d\n", object.outval[i])
                end
            end         
        else 
            for i=1:length(object.outval)
                if i < length(object.outval)
                    @printf("%.5f, ", object.outval[i])
                else 
                    @printf("%.5f\n", object.outval[i])
                end
            end
        end
    end
    @printf("Number of outliers: %d\n", object.nout)
    @printf("Non-outlier ID:     ")
    for i=1:length(object.keepid)
        if i < length(object.keepid)
            @printf("%d, ", object.keepid[i])
        else 
            @printf("%d\n", object.keepid[i])
        end
    end
end    

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


function show(io::IO, object::idealfOutput)
    if object.method!=nothing
        println(io, "$(object.method)")
    end
    @printf("Lower quartile:  %.6f\n", object.lower_quartile) 
    @printf("Upper quartile:  %.6f\n", object.upper_quartile)
end    

function show(io::IO, object::binomciOutput)
   @printf(" phat:                %.4f\n", object.phat)
   @printf(" confidence interval: %.4f   %.4f\n", object.confint[1], object.confint[2])
   @printf(" Sample size          %d\n", object.n)
end    


function complete_cases(x::DataArray)
    xna=sum(!isna(x), 2)
    ind=Int[]
    n=size(x,1)
    ncols=size(x,2)
    vec=1:n
    for i=1:n
        if xna[i]==ncols
            push!(ind, i)
        end
    end
    return vec[ind]
end

function deleterows!(x::DataArray)
    keep_inds::Array{Int}=complete_cases(x)
    output=DataArray(Any, length(keep_inds), size(x,2))
    for i=1:size(x,2)
        for j=1:length(keep_inds)
            output[j,i]=x[keep_inds[j],i]
        end
    end
    return output
end

#Finding duplicated values
function duplicated{S <: Real}(x::Vector{S})
    xsort=sort(x)
    dup=false
    for i=1:(length(xsort)-1)
        if xsort[i]==xsort[i+1]
            dup=true
            break
        end
    end
    return dup
end
#duplicated(x::Vector{Int})=duplicated(convert(Vector{Float64}, x))
#duplicated(x::DataVector)=duplicated(removeNA(x))

function seq(Start::Number, End::Number, increment::Number)
    delta=End-Start
    if delta.<=0
        error("`start` must be larger than `end`")
    else
        n::Int=floor(delta/increment)+1
        output=zeros(n)
        for i=1:n
            output[i]=(i-1).*increment+Start
        end
    end
    return output
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


#Faster Pearson's r for 
function cor!(x, y)
    n=length(x)
    xbar=mean(x); ybar=mean(y)
    top=bot1=bot2=0.0
    for i=1:n
        top+=(x[i]-xbar).*(y[i]-ybar)
        bot1+=(x[i]-xbar)*(x[i]-xbar)
        bot2+=(y[i]-ybar)*(y[i]-ybar)
    end
    top./sqrt(bot1.*bot2)
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

function indt_mflag(x)
    n=size(x, 1)
    np=size(x, 2)
    output=falses(n, n)
    for i=1:n
        for j=1:n
            total=0
            for k=1:np
                total+=x[i,k]<=x[j,k]?1:0
            end
            if total==np
                output[i,j]=true
            end
        end
    end
    return output
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

function mad!(v::AbstractArray)
    v=copy(v)
    center = median!(v, checknan=false)
    for i in 1:length(v)
        v[i] = abs(v[i]-center)
    end
    1.4826 * median!(v, checknan=false)
end

########################################################


include("functions.jl")
include("data.jl")
end


