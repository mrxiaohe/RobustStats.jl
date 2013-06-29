#Datasets used for testing the functions

x=[1.672064, 0.7876588, 0.317322, 0.9721646, 0.4004206, 1.665123, 3.059971, 0.09459603, 1.27424, 3.522148, 
   0.8211308, 1.328767, 2.825956, 0.1102891, 0.06314285, 2.59152, 8.624108, 0.6516885, 5.770285, 0.5154299]

srand(1)
x2=rnorm(21)
srand(10)
x3=rnorm(100, 1, 2);

srand(1)
m = randn(20);        
srand(2)             
y = randn(20) + 2.0;  

srand(2)
xnew = rand(6, 4)
xold = reshape(x, 5, 4)

srand(1)
c = cnorm(2000)

srand(3)
y2 = randn(20)+3; 