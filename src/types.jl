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

testOutput()=testOutput(nothing,
                        nothing,
                        nothing,
                        nothing,
                        nothing,
                        nothing,
                        nothing,
                        nothing,
                        nothing,
                        nothing)