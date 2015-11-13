import sys
import numpy as np
import matplotlib.pyplot as plt

def ReadDataFromFile(filename, function):
    
    #Read data from file
    data = np.loadtxt(filename, unpack=True)
    
    T = data[0]                 #Temperature [kT/J].
    CyclesMC = data[1]          #Total amount of Monte Carlo cycles.
    AcceptedConfigs = data[2]	#Toatal amount of accepted configurations.
    expectE = data[3]           #Mean energy.
    expectAbsM = data[4]        #Mean magnetization. (absolute value) 
    varianceE = data[5]         #Variance of energy.
    CV = data[6]                #Specific heat.
    chi = data[7] 		   	    #Susceptibility.
    
    #Return only what the calling function needs:
    if function == 'PlotMeanMC':
        return CyclesMC, expectE, expectAbsM
    elif function == 'PlotAcceptedConfigsMC':
        return CyclesMC, AcceptedConfigs
    elif function == 'CalculateProb':
        return expectE, varianceE
    elif function == 'PlotMeanTemp':
        return T, expectE, expectAbsM, CV, chi
    elif function == 'CalculateCriticalTemp':
        return T, chi


    
def PlotMeanMC():
    
    Temperatures = [1.0, 2.4]
    L = 20                 #Size of latice dimensions.
    InitStates = [0, 1]    #0 for all up, 1 for random.
    
    #Set up figures
    plt.figure(1)
    plt.title('Mean Energy, L = %d, T = %.1f' % (L,Temperatures[0]))
    plt.xlim(0,5000)
    plt.ylim(-2.1,-0.6)
    plt.xlabel('Monte Carlo cycles')
    plt.ylabel(r'$\langle E \rangle$')
    
    plt.figure(2)
    plt.title('Mean Magnetization, L = %d, T = %.1f' % (L,Temperatures[0]))
    plt.xlim(0,25000)
    plt.ylim(0.2,1.2)
    plt.xlabel('Monte Carlo cycles')
    plt.ylabel(r'$\langle |M| \rangle$')
    
    plt.figure(3)
    plt.title('Mean Energy, L = %d, T = %.1f' % (L,Temperatures[1]))
    plt.xlim(0,25000)
    #plt.ylim()
    plt.xlabel('Monte Carlo cycles')
    plt.ylabel(r'$\langle E \rangle$')
    
    plt.figure(4)
    plt.title('Mean Magnetization, L = %d, T = %.1f' % (L,Temperatures[1]))
    #plt.xlim()
    #plt.ylim()
    plt.xlabel('Monte Carlo cycles')
    plt.ylabel(r'$\langle |M| \rangle$')
    
    for Temp in Temperatures:
        for R in InitStates:
        
            #Read data from file
            filename = 'Data/ExpValsMC_T%.1f_L%d_R%d.txt' % (Temp, L, R)
            CyclesMC, expectE, expectAbsM = ReadDataFromFile(filename, function='PlotMeanMC')
        
            if R == 0:
                label = 'Initial Config: All spins up'
            else:
                label = 'Initial Config: Random'
        
            #Fill figures
            if Temp == 1.0:
                plt.figure(1)
                plt.plot(CyclesMC, expectE, label=label)
                plt.legend()
        
                plt.figure(2)
                plt.plot(CyclesMC, expectAbsM, label=label)
                plt.legend()
                
            elif Temp == 2.4:
                plt.figure(3)
                plt.plot(CyclesMC, expectE, label=label)
                plt.legend()
        
                plt.figure(4)
                plt.plot(CyclesMC, expectAbsM, label=label)
                plt.legend()

    plt.show()

def PlotAcceptedConfigsMC():

    Temperatures = [1.0, 2.0, 2.4]
    L = 20              #Size of latice dimensions.
    InitState = 0       #All spins up as initial state.
    
    #Set up figure
    plt.figure(5)
    plt.title('Accepted Configs, L = %d' %L)
    plt.xlim(0,2000)
    plt.ylim(-1,122)
    plt.xlabel('Monte Carlo cycles')
    plt.ylabel('Accepted Configs')
    
    for Temp in Temperatures:
        
        #Read data from file
        filename = 'Data/ExpValsMC_T%.1f_L%d_R%dA.txt' % (Temp, L, InitState)
        CyclesMC, AcceptedConfigs = ReadDataFromFile(filename, function='PlotAcceptedConfigsMC')
        
        #Fill figures
        label = 'T = %.1f' %Temp
        plt.figure(5)
        plt.plot(CyclesMC, AcceptedConfigs, label=label)
        plt.legend(loc='lower right')
        
    plt.show()
    
def PlotProb():
    
    Temperatures = [1.0, 2.4]
    L = 20                  #Size of latice dimensions.
    InitState = 0           #All spins up as initial state.
    SliceStart = 15000      #Use values after equilibrium 
                            #state has been reached.
    
    for Temp in Temperatures:
        
        #Read data from files
        filenameExpVals = 'Data/ExpValsMC_T%.1f_L%d_R%d.txt' % (Temp, L, InitState)
        filenameEnergies = 'Data/Energies_T%.1f_L%d_R%d.txt' % (Temp, L, InitState)
        
        expectE, varianceE = ReadDataFromFile(filenameExpVals, function='CalculateProb')
        Energies = np.loadtxt(filenameEnergies, unpack=True)
        
        #Use same amount of decimals for both 
        #datasets so we can count later:
        expectE = np.round(expectE, decimals=2)
        Energies = np.round(Energies, decimals=2)
        #Use values after equilibrium state has been reached
        Energies = Energies[SliceStart:]
        
        #Plot probability distribution in a histogram:
        hist, bin_edges = np.histogram(Energies, normed=1)
        width = bin_edges[1] - bin_edges[0]
        plt.bar(bin_edges[:-1], hist*width, width)
        plt.title(r'Probability distribution, L = %d, T = %.1f, $\langle E\rangle$ = %.2f, $\sigma^2_E$ = %.2f' % (L, Temp, expectE[-1], varianceE[-1]))
        plt.xlabel(r'$E$')
        plt.ylabel(r'$P(E)$')
        plt.show()
    
        print "T = %.2f" % Temp
        print "<E> = %.2f" % expectE[-1]
        print "Var(E) = %.2f" % varianceE[-1]
        #The sum of the probabilities for all E should be 1:
        assert np.round(width*hist.sum()) == 1.0
        print "Sum = %.1f\n" %(width*hist.sum())
        
def PlotMeanTemp():
    
    Ls = [20, 40, 60, 80]   #Sizes of latice dimensions.
    InitState = 0           #All spins up as initial state.
    
    #Set up figures
    plt.figure(6)
    plt.title('Mean Energy as Function of Temperature')
    plt.xlabel(r'$T \hspace{1} [kT/J]$')
    plt.ylabel(r'$\langle E \rangle$')
    
    plt.figure(7)
    plt.title('Mean Magnetization as Function of Temperature')
    plt.xlabel(r'$T \hspace{1} [kT/J]$')
    plt.ylabel(r'$\langle |M| \rangle$')
    
    plt.figure(8)
    plt.title('Specific Heat as Function of Temperature')
    plt.xlabel(r'$T \hspace{1} [kT/J]$')
    plt.ylabel(r'$C_V$')
    
    plt.figure(9)
    plt.title('Susceptibility as Function of Temperature')
    plt.xlabel(r'$T \hspace{1} [kT/J]$')
    plt.ylabel(r'$\chi$')
    
    for L in Ls:
    
        #Read data from file:
        filename = 'Data/ExpValsTemp_L%d_R%d.txt' %(L, InitState)
        T, expectE, expectAbsM, CV, chi = ReadDataFromFile(filename, function='PlotMeanTemp')
        
        #Fill figures
        label = 'L = %d' %L
        plt.figure(6)
        plt.plot(T, expectE, label=label)
        plt.legend(loc='best')
        
        plt.figure(7)
        plt.plot(T, expectAbsM, label=label)
        plt.legend(loc='best')
        
        plt.figure(8)
        plt.plot(T, CV, label=label)
        plt.legend(loc='best')
    
        plt.figure(9)
        plt.plot(T, chi, label=label)
        plt.legend(loc='best')
        
    plt.show()
    
def CalculateCriticalTemp():

    Ls = [20, 40, 60, 80]    #Sizes of latice dimensions.
    InitState = 0            #All spins up as initial state.
    
    exactTC = 2.269
    #Constants:
    nu = 1.
    a = 1.
    
    print "exact TC = %.3f" %exactTC
    
    for L in Ls:
    
        #Read data from file
        filename = 'Data/ExpValsTemp_L%d_R%d.txt' %(L, InitState)
        T, chi = ReadDataFromFile(filename, function='CalculateCriticalTemp')
        
        #Estimating the critical temperature 
        #based on the numerically calculated 
        #values for the susceptibility
        index = np.argmax(chi)
        Tmax = T[index]
        TC = Tmax - a*L**(-1./nu)
        
        print "L = %d, TC = %.3f" %(L, TC)
        
    
if __name__ == '__main__':
    PlotMeanMC()
    PlotAcceptedConfigsMC()
    PlotProb()
    PlotMeanTemp()
    CalculateCriticalTemp()
