function show(io::IO, object::testOutput)
    if object.method != nothing
        @printf("%s\n", object.method)
    end
    if object.h0 != nothing
        @printf(" Alternative hypothesis:   % s\n", object.h0)
    end
    if object.n != nothing
        @printf(" Sample size:         ")
        for i = 1:length(object.n)
            @printf("% d  ", object.n[i])
        end
        @printf("\n")
    end
    if object.df != nothing
            @printf(" Degrees of freedom:  ")
        for i = 1:length(object.df)
            if isa(object.df[i], Integer)
                @printf(" %d  ", object.df[i])
            elseif isa(object.df, Number)
                @printf("% .2f  ", object.df[i])
            end
        end
       @printf("\n")
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
        if isa(object.ci, Integer)
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
        if isa(object.outval, Integer)
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

function show(io::IO, object::binomciOutput)
   @printf(" p_hat:               %.4f\n", object.p_hat)
   @printf(" confidence interval: %.4f   %.4f\n", object.confint[1], object.confint[2])
   @printf(" Sample size          %d\n", object.n)
end
