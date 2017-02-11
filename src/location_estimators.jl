using StatsBase

"""`bisquareWM(x, k, center, tol)`

Return the bisquare weighted mean, a location estimator.
The bisquare weighted mean of the data `x` with a k-value of `k`.
A sensible choice of `k` is 3 to 5 times the rms width or 1.3 to 2 times the
full width at half max of a peak.  For strictly Gaussian data, the choices of
k= 3.14, 3.88, and 4.68 times sigma will be 80%, 90%, and 95% efficient.
Data values a distance of more than `k` from the weighted mean are given no
weight.

`center` is used as an initial guess at the weighted mean; if omitted,
then the data median will be used.

The answer is found iteratively, revised until it changes by less than `tol`.  If
`tol` is omitted, then `tol` will use 1e-5 times the median absolute
deviation of `x` about its median."""
function bisquareWM{T <: Real}(x::AbstractArray{T}, k::Real, center::Real, tol::Real)
    for _iteration in 1:100
        weights = (1.0-((x-center)/k).^2).^2
        weights[abs.(x-center).>k] = 0.0
        newcenter = dot(weights, x)/sum(weights)
        if abs(newcenter - center)<tol
            return newcenter
        end
        center = newcenter
    end
    error("bisquareWM used too many iterations.\nConsider using higher `tol` or better `center`, or change to trimean(x).")
end


bisquareWM{T<:Real}(x::AbstractArray{T}, k::Real, center::Real) =
        bisquareWM(x, k, center, 1e-5*StatsBase.mad(float(x)))
bisquareWM{T<:Real}(x::AbstractArray{T}, k::Real) = bisquareWM(x, k, median(x))


"""huberWM(x, k, center, tol)`

Return Huber's weighted mean, a location estimator.
Huber's weighted mean of the data `x` with a k-value of `k`.
A sensible choice of `k` is 1 to 1.5 times the rms width or 0.4 to 0.6 times the
full width at half max of a peak.  For strictly Gaussian data, the choices of
k=1.0 and 1.4 sigma give ...
Data values a distance of more than `k` from the weighted mean are given no
weight.

`center` is used as an initial guess at the weighted mean.
If `center` is None, then the data median will be used.

The answer is found iteratively, revised until it changes by less than `tol`.  If
`tol` is None (the default), then `tol` will use 1e-5 times the median absolute
deviation of `x` about its median."""
function huberWM{T <: Real}(x::AbstractArray{T}, k::Real, center::Real, tol::Real)
    for _iteration = 1:100
        weights = float(k)./abs.(x-center)
        weights[weights.>1.0] = 1.0
        newcenter = dot(weights, x)/sum(weights)
        if abs(newcenter - center) < tol
            return newcenter
        end
        center = newcenter
    end
    error("huberWM used too many iterations.\nConsider using higher `tol` or better `center`, or change to trimean(x).")
end

huberWM{T<:Real}(x::AbstractArray{T}, k::Real, center::Real) =
        huberWM(x, k, center, 1e-5*StatsBase.mad(float(x)))
huberWM{T<:Real}(x::AbstractArray{T}, k::Real) = huberWM(x, k, median(x))



"""`trimean(x)`

Return Tukey's trimean, a location estimator, for a data set `x`.

If (q1,q2,q3) are the quartiles (i.e., the 25%ile, median, and 75 %ile),
the trimean is (q1+q3)/4 + q2/2. Uses `idealf(x)` for finding quartiles.
"""
function trimean{T <: Real}(x::AbstractArray{T})
    mean([mean(idealf(x)), median(x)])
end
