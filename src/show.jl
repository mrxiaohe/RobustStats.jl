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

type binomciOutput
    phat::Number
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

