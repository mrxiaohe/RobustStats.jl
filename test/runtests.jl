using Base.Test
using RobustStats

x=[1.672064, 0.7876588, 0.317322, 0.9721646, 0.4004206, 1.665123, 3.059971, 0.09459603, 1.27424, 3.522148,
   0.8211308, 1.328767, 2.825956, 0.1102891, 0.06314285, 2.59152, 8.624108, 0.6516885, 5.770285, 0.5154299]
y = sort(x)

# Tests are based on examples from README

@test mean(x) ≈ 1.853401259
@test tmean(x) ≈ 1.2921802666666669
@test tmean(x, tr=0) ≈ 1.853401259
@test tmean(x, tr=0.3) ≈ 1.1466045875000002
@test tmean(x, tr=0.5) ≈ 1.1232023

# Winsorized data
winval(x)      # just verify that it doesn't crash
@test winmean(x) ≈ 1.4205834800000001
@test winvar(x) ≈ 0.998659015947531
@test RobustStats.winstd(x) ≈ 0.9993292830431474
@test trimse(x) ≈ 0.3724280347984342
out = trimci(x)
@test out.estimate ≈ 1.2921802666666669
@test out.p ≈ 0.005243565819244678

@test stein1(x, 1) == 41

q1,q3 = idealf(x)
@assert length(x) == 20
@test q1 ≈ y[5]*7/12+y[6]*5/12
@test q3 ≈ y[15]*5/12+y[16]*7/12

@test pbvar(x) ≈ 2.0009575278957623

#srand(2)
#y = randn(20) + 2.0
y = [2.73962,1.25549,1.39149,0.276543,1.32438,2.55665,1.13842,2.51809,3.24823,3.16911,1.90698,0.707597,1.98813,2.57913,1.31969,1.73992,2.2905,1.63649,3.24779,2.27172]
# 31. pcorb()
corval = pcorb(x, y)
@test_approx_eq_eps corval.estimate 0.318931 1e-6
@test_approx_eq_eps corval.ci[1] -0.106467 1e-2
@test_approx_eq_eps corval.ci[2] 0.663678 1e-2

#srand(3)
#y2 = randn(20)+3;
y2 = [4.19156,0.480267,5.07481,2.02675,2.89839,1.45749,3.10079,2.99803,4.00879,3.84422,4.15807,2.52484,2.75539,3.07187,2.48719,1.41976,3.51166,2.29273,1.8339,2.56761]

# 32. yuend()
res = yuend(x, y2)
@test_approx_eq_eps res.estimate -1.547776 1e-6
@test_approx_eq_eps res.se 0.460304 1e-6
@test_approx_eq_eps res.p 0.006336 1e-4

# Basic tests of the weighted means
a = collect(-2:2)
@test_approx_eq_eps bisquareWM(a,3,.1,1e-5) 0.0 1e-4
@test_approx_eq_eps huberWM(a,3,.1,1e-5) 0.0 1e-4
push!(a, 97)
@test_approx_eq_eps bisquareWM(a,3,.1,1e-5) 0.0 1e-4
push!(a, 98)
@test_approx_eq_eps bisquareWM(a,3,.1,1e-5) 0.0 1e-4
push!(a, 99)
@test_approx_eq_eps bisquareWM(a,3,.1,1e-5) 0.0 1e-4
push!(a, 98)
@test_approx_eq_eps bisquareWM(a,3,.1,1e-5) 0.0 1e-4
append!(a, [98,98,98])
@test_approx_eq_eps bisquareWM(a,3,97,1e-4) 98.0 1e-3

@test_throws ErrorException bisquareWM(collect(0:9), 0.1)

# Test 10 sets of 100 N(0,1) random numbers
# Add up to 8 values of much larger tailiness
# Expected mean is within [-.2, +.2] (2-sigma), so test for
# being in the [-1, +1] range, even with added large values.
for i = 1:10
    r = randn(100)
    for j = 1:5
        b = bisquareWM(r, 4, 0, 0.01)
        h = huberWM(r, 1.3, 0, 0.01)
        @test_approx_eq_eps b 0 1
        @test_approx_eq_eps h 0 1
        @test_approx_eq_eps b h 1.5
        append!(r, randn(2)*100)
    end
end


# Trimean: exact quartiles
a = collect(0:12)
@test_approx_eq trimean(a) 6.0
@test_approx_eq trimean(a.^2) (3^2+9^2)/4.+6^2/2.

# Trimean: exact median, inexact 1st and 3rd quartiles
a = collect(0:10)
@test_approx_eq trimean(a) 5.0
@test_approx_eq trimean(a.^2) (6.5+56.5)/4. + 25.0/2

# Trimean: inexact quartiles
a = collect(0:9)
@test_approx_eq trimean(a) 4.5
@test_approx_eq trimean(a.^2) (4*.75+9*.25 + 36*.25+49*.75)/4. + (16*.5+25*.5)/2


# Shortest half-range
a = [10,3,5,6,6.5,7,8,0,13]  # Odd # of values. 9 values, with "half" containing 5
@test_approx_eq shorthrange(a) 3
@test_approx_eq a[1] 10 # Be sure it didn't rearrange a
@test_approx_eq shorthrange!(a) 3
@test_approx_eq a[1] 0  # Should have sorted a
a = [0,4.5,6,7,8,10,14,99]  # Even # of values. 8 values, with each "half" containing 5
@test_approx_eq shorthrange(a) 10-4.5

# Weighted high median
whm = RobustStats._weightedhighmedian
whm! = RobustStats._weightedhighmedian!

# Check that whm and whm! give correct and equal answers.
function _verify_whm{T<:Real,U<:Integer}(a::Vector{T}, wts::Vector{U})
	answer = whm(a,wts)
	wlow = whigh = wexact =  0
	for i=1:length(a)
		if a[i] < answer
			wlow += wts[i]
		elseif a[i] > answer
			whigh += wts[i]
		else
			wexact += wts[i]
		end
	end
	wtotal = wlow + wexact + whigh
	@test 2*wlow <= wtotal && 2*whigh < wtotal

	answermangled = whm!(a,wts)
	@test answer == answermangled
end


_verify_whm(collect(1:5), [1,1,1,1,1])
_verify_whm(collect(1:5), [1,2,3,4,5])
_verify_whm(collect(1:5), [1,2,3,4,9])
_verify_whm(collect(1:5), [1,2,3,4,10])
_verify_whm(collect(1:5), [2,1,1,1,1])
_verify_whm(collect(1:5), [1,1,1,2,1])

_verify_whm(collect(5:-1:1), [1,1,1,1,1])
_verify_whm(collect(5:-1:1), [1,2,3,4,5])
_verify_whm(collect(5:-1:1), [1,2,3,4,9])
_verify_whm(collect(5:-1:1), [1,2,3,4,10])
_verify_whm(collect(5:-1:1), [1,2,3,4,11])
_verify_whm(collect(5:-1:1), [2,1,1,1,1])
_verify_whm(collect(5:-1:1), [1,1,1,2,1])
_verify_whm([1,4,2,5,3,6], [1,4,2,5,3,6])
_verify_whm([1,4,2,5,3,6], [1,4,2,5,3,5])
_verify_whm([1,4,2,5,3,6], [1,4,2,5,3,4])

datalengths = [10,11,1000,1111]
wttypes = [Int8, Int16, Int32, Int64, UInt8, UInt16, UInt32, UInt64]
for N in datalengths
	for wttype in wttypes
		a = randn(N)
        w = rand(wttype, N)
		# Careful! Can't have sum of weights overflow, and can't have negative weights
		for k=1:N
			if w[k]<0
				w[k] = 1
			end
			w[k] %= 8192
		end
		result = whm(a,w)
		_verify_whm(a, w)
		# println("Success on size $(N) and type $(eltype(w))")
	end
end
@test_throws ArgumentError whm(collect(1:5),collect(1:4))


@test_approx_eq RobustStats._slow_scaleQ([1,2,3,4,5,10]) scaleQ!([1,2,3,4,5,10])
@test_approx_eq RobustStats._slow_scaleS([1,2,3,4,5,10]) scaleS!([1,2,3,4,5,10])
@test_approx_eq RobustStats._slow_scaleQ([1,2,3,4,5,10.5]) scaleQ!([1,2,3,4,5,10.5])
@test_approx_eq RobustStats._slow_scaleS([1,2,3,4,5,10.5]) scaleS!([1,2,3,4,5,10.5])

NTESTS = 10
for N in [2,3,4,5,6,7,8,9,10,11,12,15,20,25,50,100,151,200,225,250,299,350,399,500]
    for _=1:NTESTS
        a = randn(N)
        Q = RobustStats._slow_scaleQ(a)
        S = RobustStats._slow_scaleS(a)
        @test_approx_eq Q scaleQ(a)
        @test_approx_eq Q scaleQ!(copy(a))
        @test_approx_eq S scaleS(a)
        @test_approx_eq S scaleS!(a)
        @test_approx_eq S scaleS!(a) # scaleS! sorts the array only, so result is unchanged
        @test_approx_eq Q scaleQ!(a)
        @test_approx_eq Q scaleQ!(a) # scaleQ! sorts the array only, so result is unchanged
    end
    println("Success on $NTESTS scaleQ and scaleS tests with size $(N)")
end
