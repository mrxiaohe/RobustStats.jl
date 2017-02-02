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
