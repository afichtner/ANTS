from cmath import sqrt, exp, pi

import scipy as sp
import numpy as np
import matplotlib.pyplot as plt

import warnings

def getkernel():
    
    xmin = -20000.0
    xmax = 20000.
    dx = 40.
    
    ymin = -6000.
    ymax = 6000.
    dy = 40.
    
    #fmin = 0.005
    #fmax = 0.0051
    #df = 0.005
    f = 0.01
    
    x1 = -2000.
    y1 = 0.
    
    x2 = 2000.
    y2 = 0.
    
    c = 3.0
    rho = 3.0
    Q = 100.
    
    sigma = 100.0
    
    #if fmin < (3.*c)/min((xmax-xmin),(ymax-ymin)):
    #    msg = 'Longest wavelengths fit into model domain less than three times. Consider using higher fmin.'
    #    warnings.warn(msg)
    #    
    #if fmax > c/(2.*min(dx,dy)):
    #    msg = 'Spatial sampling is not sufficient for shortest wavelength. Consider using lower frequencies or decrease dx, dy.'
    #

    
    #K = kernel_g2d_example(x=[xmin,xmax,dx],y=[ymin,ymax,dy],f=[fmin,fmax,df], recx1=x1,recy1=y1,recx2=x2,recy2=y2,c=c,rho=rho,sigma=sigma)
    K = kernel_g2d_example(x=[xmin,xmax,dx],y=[ymin,ymax,dy],f=f, recx1=x1,recy1=y1,recx2=x2,recy2=y2,c=c,rho=rho,sigma=sigma)
    
    #- Write to file so you can compare with matlab
    ofid_re = open('examplekernel_re.txt','w')
    ofid_im = open('examplekernel_im.txt','w')
    
    for i in np.arange(0,len(K[:,1])):
        for j in np.arange(0,len(K[1,:])):
            
            ofid_re.write('%9.6f ' %np.real(K[i,j]))
            ofid_im.write('%9.6f ' %np.imag(K[i,j]))
            
        ofid_re.write('\n')
        ofid_im.write('\n')
    ofid_re.close()
    ofid_im.close()
    
#==============================================================================
#==============================================================================
#==============================================================================


def kernel_g2d_example(x=[-20000.0e3,20000.0e3,40.0e3],y=[-6000.0e3,6000.0e3,40.0e3],\
    f=0.01,recx1=-2000.0e3,recy1=0.0,recx2=2000.0e3,recy2=0.0,\
    c=3000.0,rho=3000.0,Q=100.0,sigma=100.0,saveeps=False,savepng=True,filter=True):
    
    """
    Compute 2-D analytic noise source kernels. 

    x:                          x coordinates in km
    y:                          y coordinates in km
    f:                          frequency in Hz
    recx1, recx2, recy1, recy2: receiver coordinates in km
    c:                          phase velocity
    Q:                          Q
    sigma:                      width of the Gaussian measurement window in s
    filter:                     if True, the kernels are spatially filtered to half the wavelength of the waves 0.5*c/f.
    
    """
    
    
    #========================================================================
    #= Make a grid :)
    #========================================================================
    
    x_range=np.arange(x[0],x[1]+0.5*x[2],x[2])
    y_range=np.arange(y[0],y[1]+0.5*y[2],y[2])
    
    omega=2.0*np.pi*f
      
    xv,yv=np.meshgrid(x_range,y_range)
    
    #========================================================================
    #= Compute kernel
    #========================================================================

    K=kern_g2d(xv,yv,recx1,recy1,recx2,recy2,omega,c,rho,Q,sigma)

    K=filter_kern2d(x_range,y_range,K,0.5*c/f)
    
    #========================================================================
    #= Plot kernel
    #========================================================================

    cmap = plt.get_cmap('RdBu')
    plt.pcolormesh(xv,yv,K,cmap=cmap,shading='interp')
    plt.axis('image')
    plt.clim(-np.max(np.max(K))*0.5,np.max(np.max(K))*0.5)
    plt.colorbar()
    plt.title('Noise source kernel $[m^2 N^{-2} s^{-1}]$')

    #========================================================================
    #= Save figure
    #========================================================================
    
    figname = ("TEST/%4.3ffHz_%gkms-1_%gkgm3-1_spacing%gkm" %(f,c,rho,recx2-recx1))
    if saveeps == True:
        plt.savefig(figname+'.eps', format='eps', dpi=200)
    if savepng == True:
        plt.savefig(figname+'.png', format='png', dpi=200)
    
    plt.show()

    #========================================================================
    #= Compute integral in y direction
    #========================================================================

    K_line=np.sum(K,axis=0)
            
    plt.plot(x_range,K_line,'k',linewidth=1.4)
    plt.plot(recx1,0.5*(max(K_line)-min(K_line))+min(K_line),'rv',markersize=10)
    plt.plot(recx2,0.5*(max(K_line)-min(K_line))+min(K_line),'rv',markersize=10)
    plt.xlabel('Distance (m)')
    plt.ylabel('Summed Kernel')
        
    figname = ("TEST/Line%4.3fHz_%gkms-1_%gkgm3-1_spacing%gkm" %(f,c,rho,recx2-recx1))
            
    if saveeps == True:
        plt.savefig(figname+'.eps', format='eps', dpi=200)
    if savepng == True: 
        plt.savefig(figname+'.png', format='png', dpi=200)
        
    plt.show()
    
    return K

#==============================================================================
#==============================================================================
#==============================================================================

def get_humkernel(stationdist,numseg):
    print 'Nuts and kernels.'
    
    
#==============================================================================
#==============================================================================
#==============================================================================

def kern_g2d(x,y,x1,y1,x2,y2,omega,c,rho,Q,sigma):
    """
    Compute noise correlation source kernel. Input: see above.
    """
         
    rec_dist=np.sqrt((x1-x2)**2+(y1-y2)**2)

    K=np.zeros(np.shape(x),dtype=np.complex)
    
    #- Green's functions
    G1=green2d(x,y,x1,y1,omega,c,rho,Q)
    G2=green2d(x,y,x2,y2,omega,c,rho,Q)
    #- Adjoint source
    f=adj_src(x,y,x1,y1,x2,y2,omega,c,rho,Q,rec_dist,sigma)
    K+=np.multiply(np.multiply(G1,np.conj(G2)),f)

    return 2.0*np.real(K)


#==============================================================================
#==============================================================================
#==============================================================================

def green2d(x,y,xr,yr,omega,c,rho,Q):
    
    """
    Analytic 2-D Green's function for a homogeneous plane.
    
    
    Input:
    
    :type x,y: numpy array
    :param x,y: Source coordinate grid
    
    :type xr,yr: float
    :param xr,yr: Receiver coordinates
    
    :type omega: float
    :param omega: Circular frequency
    
    :type c: float
    :param c: Propagation velocity
    
    :type rho: float
    :param rho: Density of medium
    
    :type Q: float
    :param Q: Quality factor of the medium
    
    
    Output:
    
    :type G: numpy array
    :param G: 2-D analytic Green's function received in point x,y at frequency omega
    
    """

    dx = x[1,2]-x[1,1]
    dy = y[2,1]-y[1,1]
    
    radius=np.sqrt((x-xr)**2+(y-yr)**2)
    radius+=(dx+dy)/4.0

    A=(1.0/(4.0*1.j*rho*c**2))*np.sqrt(2.0*c/(np.pi*radius*np.abs(omega)))
    G=A*np.exp(-1.j*(omega*radius/c-np.pi/4.0))
    G=G*np.exp(-np.abs(omega)*radius/(2.0*c*Q))

    return G
    
#==============================================================================
#==============================================================================
#==============================================================================

def adj_src(x,y,x1,y1,x2,y2,omega,c,rho,Q,rec_dist,sigma):
    """
    Compute adjoint source function. The measurement is the logarithmic energy ratio
    of the causal and anticausal parts of the correlation function.
    """

    # Frequency axis:
    omega_min=-5.0*omega
    omega_max=5.0*omega
    domega=omega/10.0

    ws = np.arange(-5.0*omega,5.0*omega+omega/10.0,omega/10.0)*np.pi*2.0
    
    # Adjoint source terms
    f1=0.
    f2=0.

    # Scaling
    E=0.

    # March through the frequencies to compute convolution and scaling terms
    for w in ws:
    
        if np.abs(w)>0.1*domega:

            #- correlation function for frequency w
            corr=corr_g2d(x,y,x1,y1,x2,y2,w,c,rho,Q)

            #- time-shifted Gaussian
            G=1.0/(2*sigma*np.sqrt(np.pi))
            G=G*np.exp(-0.25*(omega-w)**2*sigma**2)

            #- convolution product
            f1=f1+corr*np.exp(-1.j*(omega-w)*rec_dist/c)*G;
            f2=f2+corr*np.exp(1.j*(omega-w)*rec_dist/c)*G;

            #- scaling
            E=E+(corr*np.exp(-0.5*(w**2)*(sigma**2)))**2


    f=(f1-f2)/E
    f=np.conj(f)

    return f

    

def corr_g2d(x,y,x1,y1,x2,y2,w,c,rho,Q):
    
    """
    Obtain a cross-correlation from 2-D Green's functions.
    
    """
    
    dx = x[1,2]-x[1,1]
    dy = y[2,1]-y[1,1]
    
    
    G1 = green2d(x,y,x1,y1,w,c,rho,Q)
    G2 = green2d(x,y,x2,y2,w,c,rho,Q)
    
    correlation_re = np.zeros(np.shape(G1),dtype=np.float128)
    correlation_im = np.zeros(np.shape(G1),dtype=np.float128)
    correlation_re += np.multiply(np.real(G1),np.real(G2))+np.multiply(np.imag(G1),np.imag(G2))
    correlation_im += np.multiply(np.imag(G1),np.real(G2))-np.multiply(np.real(G1),np.imag(G2))
    
    correlation = np.sum(correlation_re+1.j*correlation_im)*dx*dy

    return correlation
    
    
def filter_kern2d(x,y,K,sigma):
    """
    Apply a Gaussian moving window filter
    
    x: 1-D x axis
    y: 1-D y axis
    K: 2-D kernel, dimensions len(y)xlen(x)
    sigma: standard deviation, half width of Gaussian window
    returns filtered version of kernels, dimensions len(y)xlen(x)
    """
    
    dx = x[2]-x[1]
    dy = y[2]-y[1]
    
    K_f = np.zeros((len(y),len(x)),dtype=np.complex)
    
    #- y direction
    
    for j in np.arange(0,len(y)):
        
        #- Gaussian window
        g = (1./(np.sqrt(2.*np.pi)*sigma)) * np.exp(-np.power((y-y[j]),2) /(2.*sigma**2))
        #- Window is normalized by area under the window
        N = np.sum(g) * dy
        
        #- Convolution of this Gauss window with Kernel
        for i in np.arange(0,len(x)):
            K_f[j,i] = np.sum(K[:,i]*g)*dy/N
    
    #- x direction
            
    for i in np.arange(0,len(x)):
        
        #- Gaussian window
        g = (1./(np.sqrt(2.*np.pi)*sigma)) * np.exp(-np.power((x-x[i]),2) /(2.*sigma**2))
        #- Window is normalized by area under the window
        N = np.sum(g) * dx
        
        #- Convolution of this Gauss window with Kernel
        for j in np.arange(0,len(y)):
            K_f[j,i] = np.sum(K_f[j,:]*g)*dx/N
            
    return K_f
        
    
    
    
    
    
    