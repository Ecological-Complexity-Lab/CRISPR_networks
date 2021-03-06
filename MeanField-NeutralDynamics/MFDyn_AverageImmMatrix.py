

# Created by by Sergio A. Alcala-Corona 
# Creative Commons License: CC BY-NC-SA 4.0 (2019)

import sys
import numpy as np
import matplotlib.pyplot as plt

from scipy.integrate import odeint


###########################################################################################


# This is the simplified mean-field model (from Childs et. al.) to be integrated
def ChildsModel(y, t, phi, K, p, q, M, b):
    N,V = y
    return [N*(1-(N/K)) - ((1-q)*(1-M) + p*M)*phi*N*V, b*phi*((1-q)*(1-M) + p*M)*N*V - (phi*N + 0.1)*V]



# The function PlotNullDynamics make the numerical integration and plot the resulting phase portrait for each M value
def PlotNullDynamics(M,pos,intialConditios):

    ax = plt.subplot(pos) 

    plt.ylabel('$Virus$',fontsize=BIGGER_SIZE)  #  Putting some labels for the corresponding plot
    plt.xlim([1000, 350000])                    #  setting the limits for the plot
    plt.ylim([0, 50000000])


# Numerical integration of the model
    for y0 in intialConditios:
        tspan = np.linspace(0, 1000, 10000) 
        ys = odeint(ChildsModel, y0, tspan, args=(phi,K, p, q, M, b))  
        plt.plot(ys[:,0], ys[:,1], 'b-') # Plotting the phase portrait trajectories
        
    plt.plot(label='Virus$<M_{ij}> = %s$'%M, fontsize=10)    # plotting axes labels
        
    N = 0.1/((((1-q)*(1-M) + p*M)*b - 1)*phi)   # computing the corresponding nullcline for bacterial population
    plt.axvline(x=N)                            # and plotting it
    #print N

    x = np.linspace(1000, 350000,100000)         # computing the corresponding isolcline for viral population
    y = (1-(x/K))/(((1-q)*(1-M) + p*M)*phi)
    plt.plot(x, y, '-r')                        # and plotting it



    # plotting the trivial equilibrium points
    a = 3
    m = 0.1
    plt.plot(0+a,0+a,'go',markersize=7)
    plt.plot(K,0,'go',markersize=7)


    # computting and plotting the C3 equilibrium point
    C3 = (1/(phi*((1-q)*(1-M) + p*M))*(1-(m/(phi*K*(b*((1-q)*(1-M) + p*M)-1)))))
    plt.plot(N, C3,'mo',markersize=5)
   
    plt.text(275000,45000000, '$<M_{ij}> = %s$'%M,fontsize=MEDIUM_SIZE)  # labeling each plot
    
    
    #plt.show()



###########################################################################################

# Setting the parameters for the simplified model (from Childs et. al.)
phi = 0.0000001
K = 316227.766016838
p = 0.00001
q = 0.00001
b = 50

# Getting four M values from the command line
M1 = float(sys.argv[1])     # M1 = 0.35
M2 = float(sys.argv[2])     # M2 = 0.7
M3 = float(sys.argv[3])     # M3 = 0.85
M4 = float(sys.argv[4])     # M4 = 0.9167


# Some variables to control the size of labels in the plots
SMALL_SIZE = 4
MEDIUM_SIZE = 11.5
BIGGER_SIZE = 18


# Making the frame of the general plot
plt.subplots_adjust(top=0.975, bottom=0.055, left=0.035, right=0.980, wspace=0.05, hspace=0.1)
plt.rc('axes', labelsize=BIGGER_SIZE)    # fontsize of the x and y labels


# Some set of initial conditions. An array containing several initial points (N, V) in the phase portrait, from which the dynamics will be integrated numerically
intialConditios = [[40000,1000000],[140000,10000000],[240000,20000000],[340000,1000000],[340000,10000000],[340000,20000000],[340000,30000000],[340000,40000000]]


# Making a set of plots for the different values of and positions in the lo panel
# The function PlotNullDynamics make the numerical integration and plot the resulting phase portrait for each M value
for M,pos in [(M1,221),(M2,222),(M3,223),(M4,224)]:
    PlotNullDynamics(M,pos,intialConditios)

plt.show()
#plt.tight_layout()  
#plt.savefig("Fig_PhaseSpace_MF.png",dpi=500)
#plt.close() 


