


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

"""`duplicated(x)` returns whether `x` contains any duplicate values"""
function duplicated(x::AbstractArray)
    xsort = sort(x)
    for i = 1:(length(xsort)-1)
        if xsort[i] == xsort[i+1]
            return true
        end
    end
    false
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
