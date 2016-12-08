using Base.Test
using RobustStats

x=[1.672064, 0.7876588, 0.317322, 0.9721646, 0.4004206, 1.665123, 3.059971, 0.09459603, 1.27424, 3.522148, 
   0.8211308, 1.328767, 2.825956, 0.1102891, 0.06314285, 2.59152, 8.624108, 0.6516885, 5.770285, 0.5154299]

# Tests are based on examples from README
#1
@test_approx_eq mean(x) 1.853401259
@test_approx_eq tmean(x) 1.2921802666666669
@test_approx_eq tmean(x, tr=0) 1.853401259
@test_approx_eq tmean(x, tr=0.3) 1.1466045875000002
@test_approx_eq tmean(x, tr=0.5) 1.1232023

# 2
winval(x) # doesn't crash

# 3 
@test_approx_eq winmean(x) 1.4205834800000001

# 4
@test_approx_eq winvar(x) 0.998659015947531

# 5. trimse: estimated standard error of the gamma trimmed mean
@test_approx_eq trimse(x)   0.3724280347984342

# 7. stein1: Stein's method
@test stein1(x, 1) == 41



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






