type outOutput
    outid
    keepid
    outval
    nout
    method
end

type binomciOutput
    p_hat::Number
    confint::Array
    n::Int
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
