from pylab import *
import numpy as np

table0_01wr = np.loadtxt('omega0_01withrepulsion.dat', skiprows=0)
x = table0_01wr[:,0]
u0_01wr = table0_01wr[:,1]

table0_01nr = np.loadtxt('omega0_01norepulsion.dat', skiprows=0)
u0_01nr = table0_01nr[:,1]

table0_5wr = np.loadtxt('omega0_5withrepulsion.dat', skiprows=0)
u0_5wr = table0_5wr[:,1]

table0_5nr = np.loadtxt('omega0_5norepulsion.dat', skiprows=0)
u0_5nr = table0_5nr[:,1]

table1wr = np.loadtxt('omega1withrepulsion.dat', skiprows=0)
u1wr = table1wr[:,1]

table1nr = np.loadtxt('omega1norepulsion.dat', skiprows=0)
u1nr = table1nr[:,1]

table5wr = np.loadtxt('omega5withrepulsion.dat', skiprows=0)
u5wr = table5wr[:,1]

table5nr = np.loadtxt('omega5norepulsion.dat', skiprows=0)
u5nr = table5nr[:,1]

figure(1)

subplot(4,1,1)
plot(x, u0_01wr**2,
    x, u0_01nr**2)
title(r'Probability density $|u(\rho)|^2$ for various $\omega_r$ with N=200 and $\rho_{max}$=5')
ylabel(r'$|u(\rho)|^2$')
legend([r'$\omega_r$=0.01, With repulsion', '$\omega_r$=0.01, Without repulsion'] ,prop={'size':10})

subplot(4,1,2)
plot(x, u0_5wr**2,
    x, u0_5nr**2)
ylabel(r'$|u(\rho)|^2$')
legend([r'$\omega_r$=0.5, With repulsion', '$\omega_r$=0.5, Without repulsion'] ,prop={'size':10})

subplot(4,1,3)
plot(x, u1wr**2,
    x, u1nr**2)
ylabel(r'$|u(\rho)|^2$')
legend([r'$\omega_r$=1, With repulsion', '$\omega_r$=1, Without repulsion'] ,prop={'size':10})

subplot(4,1,4)
plot(x, u5wr**2,
    x, u5nr**2)
xlabel(r'$\rho$')
ylabel(r'$|u(\rho)|^2$')
legend([r'$\omega_r$=5, With repulsion', '$\omega_r$=5, Without repulsion'] ,prop={'size':10})

show()
