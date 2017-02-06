


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
