import sys
import numpy as np
import matplotlib.pyplot as plt

def ReadDataFromFile(filename):
    
    #Read data from file
    data = np.loadtxt(filename, unpack=True)
    
    v = data[0]     #Array of solution at time t.
    x = data[1]     #Array of x-values.
    t = data[2]	    #Point in time.
    
    return v, x, t[0]
    
def PlotData():
    times = ["0.100", "1.000"]
    dx = ["0.010", "0.100"]
    solvers = ["FE", "BE", "CN"]
    
    #Set up figures
    plt.figure(1)
    plt.title('Foward Euler\ndx = %s, dt = 2e-5' % dx[0])
    #plt.xlim(0,1)
    #plt.ylim(0,1)
    plt.xlabel('x', fontsize="large")
    plt.ylabel('u(x, t)', fontsize="large")
    
    plt.figure(2)
    plt.title('Foward Euler\ndx = %s, dt = 2e-5' % dx[1])
    #plt.xlim(0,1)
    #plt.ylim(0,1)
    plt.xlabel('x', fontsize="large")
    plt.ylabel('u(x, t)', fontsize="large")
    
    plt.figure(3)
    plt.title('Backward Euler\ndx = %s, dt = 2e-5' % dx[0])
    #plt.xlim(0,1)
    #plt.ylim(0,1)
    plt.xlabel('x', fontsize="large")
    plt.ylabel('u(x, t)', fontsize="large")
    
    plt.figure(4)
    plt.title('Backward Euler\ndx = %s, dt = 2e-5' % dx[1])
    #plt.xlim(0,1)
    #plt.ylim(0,1)
    plt.xlabel('x', fontsize="large")
    plt.ylabel('u(x, t)', fontsize="large")
    
    plt.figure(5)
    plt.title('Crank-Nicolson\ndx = %s, dt = 2e-5' % dx[0])
    #plt.xlim(0,1)
    #plt.ylim(0,1)
    plt.xlabel('x', fontsize="large")
    plt.ylabel('u(x, t)', fontsize="large")
    
    plt.figure(6)
    plt.title('Crank-Nicolson\ndx = %s, dt = 2e-5' % dx[1])
    #plt.xlim(0,1)
    #plt.ylim(0,1)
    plt.xlabel('x', fontsize="large")
    plt.ylabel('u(x, t)', fontsize="large")
    
    for time in times:
        for s in solvers:
        
            #Read data from file
            filename1 = 'Data/%s_dx%s_t%s.txt' % (s, dx[0], time)
            v1, x1, t1 = ReadDataFromFile(filename1)
            filename2 = 'Data/%s_dx%s_t%s.txt' % (s, dx[1], time)
            v2, x2, t2 = ReadDataFromFile(filename2)
            
            #Fill figures
            if s == "FE":
                plt.figure(1)
                plt.plot(x1, v1, label="t = %s" %t1)
                plt.legend()
                #plt.savefig("Plots/FE_dx001.png")
                
                plt.figure(2)
                plt.plot(x2, v2, label="t = %s" %t2)
                plt.legend()
                #plt.savefig("Plots/FE_dx01.png")
                
            if s == "BE":
                plt.figure(3)
                plt.plot(x1, v1, label="t = %s" %t1)
                plt.legend()
                #plt.savefig("Plots/BE_dx001.png")
        
                plt.figure(4)
                plt.plot(x2, v2, label="t = %s" %t2)
                plt.legend()
                #plt.savefig("Plots/BE_dx01.png")

            if s == "CN":
                plt.figure(5)
                plt.plot(x1, v1, label="t = %s" %t1)
                plt.legend()
                #plt.savefig("Plots/CN_dx001.png")
        
                plt.figure(6)
                plt.plot(x2, v2, label="t = %s" %t2)
                plt.legend()
                #plt.savefig("Plots/CN_dx01.png")

    plt.show()
    
def CompareWithAnalytical():
    
    def Analytical(x, t):
        sum_n = 0
        N = 200
        for n in range(1, N):
            sum_n += 1./n*np.sin(n*np.pi*x)*np.exp(-(n*np.pi)**2*t) 
        Analytical = 1 - x -2/np.pi*sum_n
        return Analytical
    
    times = ["0.100"]
    dx = ["0.010"]
    solvers = ["FE", "BE", "CN"]
    
    #Set up figures
    plt.figure(1)
    plt.title('Abs Error\ndx = %s, t = %s, dt = 2e-5' % (dx[0], times[0]))
    #plt.xlim(0,1)
    #plt.ylim(0,1)
    plt.xlabel('x', fontsize="large")
    plt.ylabel('Abs Error', fontsize="large")
    
    plt.figure(2)
    plt.title('Rel Error\ndx = %s, t = %s, dt = 2e-5' % (dx[0], times[0]))
    #plt.xlim(0,1)
    #plt.ylim(0,1)
    plt.xlabel('x', fontsize="large")
    plt.ylabel('Rel Error', fontsize="large")
    
    fig=1
    for time in times:
        for d in dx:
        
            #Read data from file
            filename1 = 'Data/%s_dx%s_t%s.txt' % (solvers[0], d, time)
            vFE, x, t = ReadDataFromFile(filename1)
            filename2 = 'Data/%s_dx%s_t%s.txt' % (solvers[1], d, time)
            vBE, x, t = ReadDataFromFile(filename2)
            filename3 = 'Data/%s_dx%s_t%s.txt' % (solvers[2], d, time)
            vCN, x, t = ReadDataFromFile(filename3)

            vA = Analytical(x, t)
            
            AbsErrorFE = abs(vA-vFE)
            AbsErrorBE = abs(vA-vBE)
            AbsErrorCN = abs(vA-vCN)
            
            RelErrorFE = AbsErrorFE/vA
            RelErrorBE = AbsErrorBE/vA
            RelErrorCN = AbsErrorCN/vA
            
            plt.figure(fig)
            plt.plot(x, AbsErrorFE, label="Forward Euler")
            plt.plot(x, AbsErrorBE, label="Backward Euler")
            plt.plot(x, AbsErrorCN, label="Crank-Nicolson")
            plt.legend(loc="best")
            #plt.savefig("Plots/AbsError.png")
            fig += 1
            plt.figure(fig)
            plt.plot(x[:-1], RelErrorFE[:-1], label="Forward Euler")
            plt.plot(x[:-1], RelErrorBE[:-1], label="Backward Euler")
            plt.plot(x[:-1], RelErrorCN[:-1], label="Crank-Nicolson")
            plt.legend(loc="best")
            #plt.savefig("Plots/RelError.png")
            fig +=1
    plt.show()

def PlotMonteCarlo():
    #Read data from file
    times = ["0.010", "0.100", "1.000", "10.000"]
    
    plt.figure(1)
    plt.title('Monte Carlo, Gaussian step (idum=-1)\nNumWalkers = 1e5, dt = 1e-5')
    #plt.xlim(0,1)
    #plt.ylim(0,1)
    plt.xlabel('x', fontsize="large")
    plt.ylabel('Normalized MC count(x, t)', fontsize="large")
    
    for t in times:
        filename = "Data/MC_Gauss_t%s.txt" %t
        data = np.loadtxt(filename, unpack=True)
    
        v = data[0]
        x = data[1]
        
        label = "t=%s" %t
        if t == "10.000":
            plt.plot(x, v, "--", label=label)
        else:
            plt.plot(x, v, label=label)
        plt.legend()
    #plt.savefig("Plots/MC_Gauss.png")
    plt.show()
    
    
def PlotUnstable():
    times = ["1.000"]
    dx = ["0.100"]
    solvers = ["FE", "BE", "CN"]
    
    #Set up figures
    plt.figure(1)
    plt.title('dx = %s, dt = 1e-4' % dx[0])
    #plt.xlim(0,1)
    #plt.ylim(0,1)
    plt.xlabel('x', fontsize="large")
    plt.ylabel('u(x, t)', fontsize="large")
    
    for s in solvers:
        
        #Read data from file
        filename1 = 'Data/%s_dx%s_t%s_unstable.txt' % (s, dx[0], times[0])
        v, x, t = ReadDataFromFile(filename1)
            
        #Fill figures
        plt.figure(1)
        if s=="CN":
            plt.plot(x, v, "--", label="Method: %s" %s)
        else:
            plt.plot(x, v, label="Method: %s" %s)
        plt.legend()
    #plt.savefig("Plots/Unstable.png")
    plt.show()

if __name__ == "__main__":
    #PlotData()
    #CompareWithAnalytical()
    PlotMonteCarlo()
    #PlotUnstable()
