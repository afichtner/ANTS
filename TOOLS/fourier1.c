#include <csubs.h>
#include <math.h>
#include <stdio.h>
#define SWAP(a,b) tempr=(a);(a)=(b);(b)=tempr
#define pi 3.14159265

void fourier1(float *data, int nn, int isign)
{
    int n, mmax, m, j, istep, i;
    double wtemp, wr, wpr, wpi, wi, theta;
    double tempr, tempi;
    
    /* Bit reversal of the input */
    n=nn<<1;
    j=1;
    for (i=1;i<n;i+=2) {
        if (j>i) {
            SWAP(data[j],data[i]);
            SWAP(data[j+1],data[i+1]);
        }
        
        m=n>>1;
        while (m>=2 && j>m) {
            j -=m;
            m>>=1;
        }
        
        j+=m;
        
    }
    
    /* Danielson-Lanczos loop */
    /*Don't ask me, it s complicated.*/
    mmax=2;
    while (n>mmax) {
        istep=mmax << 1;
        theta=isign*(6.28318530717959/mmax);
        wtemp=sin(0.5*theta);
        wpr=-2.0*wtemp*wtemp;
        wpi=sin(theta);
        wr=1.0;
        wi=0.0;
        
        for (m=1;m<mmax;m+=2) {
            for (i=m;i<=n;i+=istep){
                j=i+mmax;
                tempr=wr*data[j]-wi*data[j+1];
                tempi=wr*data[j+1]+wi*data[j];
                data[j]=data[i]-tempr;
                data[j+1]=data[i+1]-tempi;
                data[i]+=tempr;
                data[i+1] +=tempi;
            }
            
            wr=(wtemp=wr)*wpr-wi*wpi+wr;
            wi=wi*wpr+wtemp*wpi+wi;
            
        }
        mmax=istep;
    }
            
    
}



void tworealfft(float *data1, float *data2, float *fft1, float *fft2, int n)
{   printf("Entering 2 real fft routine.\n");
    /* Declaration */
    void fourier1(float *data, int nn, int isign);
    
    int nn3, nn2, jj, j;
    float rep, rem, aip, aim;
    nn2=2+n+n;
    nn3=1+nn2;
  
    
    for (j=1,jj=2; j<=n; j++,jj+=2) {
        fft1[jj-1]=data1[j];
        fft1[jj]=data2[j];
    }
   
    fourier1(fft1,n,1);
    printf("Successfully calculated 2 real fft.\n");
    fft2[1]=fft1[2];
    fft1[2]=fft2[2]=0.0;
    
    for (j=3;j<=n+1;j+=2) {
        
        rep=0.5*(fft1[j]+fft1[nn2-j]);
        rem=0.5*(fft1[j]-fft1[nn2-j]);
        aip=0.5*(fft1[j+1]+fft1[nn3-j]);
        aim=0.5*(fft1[j+1]-fft1[nn3-j]);
        
        fft1[j]=rep;
        fft1[j+1]=aim;
        fft1[nn2-j]=rep;
        fft1[nn3-j]=-aim;
       
        fft2[j]=aip;
        fft2[j+1]=-rem;
        fft2[nn2-j]=aip; 
        fft2[nn3-j]=rem;
        
    }
    printf("Done with Fourier transforms.\n");
}

