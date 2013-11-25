# Peterson noise model

import numpy as np
import matplotlib.pyplot as plt
from math import log10, log,  exp,  pow,  pi

def peterson(unit, plot):   
    """
    get and plot Peterson's 1993 new high/new low noise model.
    unit: string, 'ACC','VEL' or 'DIS
    plot: Boolean, plot or not
    
    output:
    period axis of NHNM, NHNM, period axis of NLNM, NLNM
    
    """
    
     
    #- NHNM ====================================================================================================
    A=[-108.73,-150.34,-122.31,-116.85,-108.48,-74.66,0.66,-93.37,73.54,-151.52,-206.66]
    B=[-17.23, -80.50, -23.87, 32.51, 18.08, -32.95, -127.18, -22.42, -162.98, 10.01, 31.63]
    P=[0.1, 0.22, 0.32, 0.80, 3.80, 4.60, 6.30, 7.90, 15.40, 20.0, 354.80, 100000]
    
    
    nhnm=np.zeros(len(P))
    nhnm_test=np.zeros(len(P))
    for i in range(len(P)-1):
        nhnm[i]=A[i]+B[i]*log10(P[i])
        nhnm_test[i]=10*log10(exp(0.1*A[i])*pow(P[i], (0.1*B[i])/log(10.0)))
    nhnm[-1]=A[-1]+B[-1]*log10(P[-1])
    nhnm_test[-1]=10*log10(exp(0.1*A[-1])*pow(P[-1], (0.1*B[-1])/log(10.0)))
    
    #- NLNM ====================================================================================================
    A1=[-162.36, -166.7, -170.00, -166.40, -168.60, -159.98, -141.10, -71.36, -97.26, -132.18, -205.27, -37.65, -114.37, -160.58, -187.50, -216.47, -185.00, -168.34, -217.43, -258.28, -346.88 ]
    B1=[5.64, 0.00, -8.30, 28.90, 52.48, 29.81, 0.00, -99.77, -66.49, -31.57, 36.16, -104.33, -47.10, -16.28, 0.00, 15.70, 0.00, -7.61, 11.90, 26.60, 48.75]
    P1=[0.1, 0.17, 0.40, 0.80, 1.24, 2.40, 4.30, 5.00, 6.00, 10.00, 12.00, 15.60, 21.90, 31.60, 45.00, 70.00, 101.00, 154.00, 328.00, 600.00, 10000.00, 100000.00]
    
    nlnm=np.zeros(len(P1))
    nlnm_test=np.zeros(len(P1))
    for i in range(len(P1)-1):
        nlnm[i]=A1[i]+B1[i]*log10(P1[i])
        nlnm_test[i]=10*log10(exp(0.1*A1[i])*pow(P1[i], (0.1*B1[i])/log(10.0)))
    nlnm[-1]=A1[-1]+B1[-1]*log10(P1[-1])
    nlnm_test[-1]=10*log10(exp(0.1*A1[-1])*pow(P1[-1], (0.1*B1[-1])/log(10.0)))
    
   
    
    if unit=='ACC':
        
        if plot:
            plt.semilogx(P, nhnm, linewidth=2.0)
            plt.semilogx(P1, nlnm)
            plt.grid()
            plt.show()
            #plt.semilogx(P, nhnm_test)
            #plt.semilogx(P1, nlnm_test)
            #plt.grid()
            #plt.show()
        
        return P, nhnm, P1, nlnm
        
        
    elif unit=='VEL':
        for i in range(0, len(P)-1):
            nhnm[i]=nhnm[i]+20*log10(P[i]/(2*pi))
            nlnm[i]=nlnm[i]+20*log10(P1[i]/(2*pi)) 
            
        if plot:
            plt.semilogx(P, nhnm, linewidth=2.0)
            plt.semilogx(P1, nlnm)
            plt.grid()
            plt.show()
           
            
        return P, nhnm, P1, nlnm
    
    elif unit=='DIS':
        for i in range(0, len(P)-1):
            nhnm[i]=nhnm[i]+20*log10(P[i]**2/(2*pi)**2)
            nlnm[i]=nlnm[i]+20*log10(P1[i]**2/(2*pi)**2)
              
        if plot:
            plt.semilogx(P, nhnm, linewidth=2.0)
            plt.semilogx(P1, nlnm)
            plt.grid()
            plt.show()
            
            
        return P, nhnm, P1, nlnm
     
