


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
