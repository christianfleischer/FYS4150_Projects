from pylab import *
import numpy as np

table10 = np.loadtxt('tridiag_solution_n10.dat', skiprows=0)
x10 = table10[:,0]
v10 = table10[:,1]

tableLU10 = np.loadtxt('LU_solution_n10.dat', skiprows=0)
vLU10 = tableLU10[:,1]

table100 = np.loadtxt('tridiag_solution_n100.dat', skiprows=0)
x100 = table100[:,0]
v100 = table100[:,1]

tableLU100 = np.loadtxt('LU_solution_n100.dat', skiprows=0)
vLU100 = tableLU100[:,1]

table1000 = np.loadtxt('tridiag_solution_n1000.dat', skiprows=0)
x1000 = table1000[:,0]
v1000 = table1000[:,1]

tableLU1000 = np.loadtxt('LU_solution_n1000.dat', skiprows=0)
vLU1000 = tableLU1000[:,1]

exact10 = 1 - (1-exp(-10))*x10 - exp(-10*x10)
exact100 = 1 - (1-exp(-10))*x100 - exp(-10*x100)
exact1000 = 1 - (1-exp(-10))*x1000 - exp(-10*x1000)

figure(1)

subplot(3,1,1)
plot(x10, v10,
    x10, exact10)
title('Convergence of tridiag solution to exact solution for increasing n')
ylabel('u(x)')
legend(['Numerical, n=10', 'Exact'] ,prop={'size':10})

subplot(3,1,2)
plot(x100, v100,
    x100, exact100)
ylabel('u(x)')
legend(['Numerical, n=100', 'Exact'] ,prop={'size':10})

subplot(3,1,3)
plot(x1000, v1000,
    x1000, exact1000)
xlabel('x')
ylabel('u(x)')
legend(['Numerical, n=1000', 'Exact'] ,prop={'size':10})

figure(2)

subplot(3,1,1)
plot(x10, v10,
    x10, vLU10)
title('Tridiagonal vs. LU Decomp.')
ylabel('u(x)')
legend(['Tri, n=10', 'LU, n=10'] ,prop={'size':10})

subplot(3,1,2)
plot(x100, v100,
    x100, vLU100)
ylabel('u(x)')
legend(['Tri, n=100', 'LU, n=100'] ,prop={'size':10})

subplot(3,1,3)
plot(x1000, v1000,
    x1000, vLU1000)
ylabel('u(x)')
xlabel('x')
legend(['Tri, n=1000', 'LU, n=1000'] ,prop={'size':10})

show()
