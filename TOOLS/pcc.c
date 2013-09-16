/*--------------------------------------------------------------------
# Filename: testpcc.c
#  Purpose: phase cross correlation in time domain
#   Author: Hansruedi Maurer
#  Changes: Laura Ermert
#  Copyright?
#---------------------------------------------------------------------*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <memory.h>
#include <csubs.h>
#define pi 3.14159265

/*--------------------------------------------------------------------
# Main 
# To be included: 
#--All in one array
#--nu
#-Testing: auto - phase - correlation
--------------------------------------------------------------------*/
void main(void){
    void asignal(float *,float *, float *, float *, int);
    
    FILE *fh_t, *fh_u, *fh_v, *fh_pcc;
    
    int a, b,i1,i2, ndat=1024*512,len,n2=2*ndat;
    float *tra1, *tra2;
    float *atra1, *atra2;
    float *re1, *im1;
    float *re2, *im2;
    float abs1, abs2;
    float am=1.0,f=0.001,cyc=10.0,t[ndat+1],u[ndat+1],dt=1.0/(20.0*f);
    printf("dt=%g\n",dt);
    int maxlag=1000/dt;
    double *pccf;
        
    printf("Starting...\n");
    printf("maxlag = %d\n", maxlag);
    /* Make the example signal*/
    for (a=1;a<=ndat;a++){
        t[a]=a*dt;
        u[a]=am*(1 -cos(2*pi*f*t[a]/cyc))*cos(2*pi*f*t[a]);}
    
    fh_t=fopen("test_pcc_t.txt","w");
    fh_u=fopen("test_pcc_u.txt","w");
   
     for (a=1;a<=ndat;a++){
        fprintf(fh_t,"%g\n",t[a]);
        fprintf(fh_u,"%g\n",u[a]);
        }
        
    fclose(fh_t);
    fclose(fh_u);
    
    printf("Example signal done.\n");
        
    /* Allocate space for the traces */
    tra1 = (float *)malloc((size_t) (ndat+1)* sizeof(float));
    if (tra1 == NULL) 
        {
            fprintf(stderr,"\nMemory allocation error!\n");
            exit(EXIT_FAILURE);
        }
    tra2 = (float *)malloc((size_t) (ndat+1)* sizeof(float));
    if (tra2 == NULL) 
        {
            fprintf(stderr,"\nMemory allocation error!\n");
            exit(EXIT_FAILURE);
        }
    atra1 = (float *)calloc((size_t) n2+1, sizeof(float));
    if (atra1 == NULL) 
        {
            fprintf(stderr,"\nMemory allocation error!\n");
            exit(EXIT_FAILURE);
        }
    atra2 = (float *)calloc((size_t) n2+1, sizeof(float));
    if (atra2 == NULL) 
        {
            fprintf(stderr,"\nMemory allocation error!\n");
            exit(EXIT_FAILURE);
        }
    re1 = (float *)malloc((size_t) (ndat+1)* sizeof(float));
    if (re1== NULL) 
        {
            fprintf(stderr,"\nMemory allocation error!\n");
            exit(EXIT_FAILURE);
        }
    re2 = (float *)malloc((size_t) (ndat+1)*sizeof(float));
    if (re2 == NULL) 
        {
            fprintf(stderr,"\nMemory allocation error!\n");
            exit(EXIT_FAILURE);
        }
    im1 = (float *)malloc((size_t) (ndat+1)*sizeof(float));
    if (im1 == NULL) 
        {
            fprintf(stderr,"\nMemory allocation error!\n");
            exit(EXIT_FAILURE);
        }
    im2 = (float *)malloc((size_t) (ndat+1)* sizeof(float));
    if (im2 == NULL) 
        {
            fprintf(stderr,"\nMemory allocation error!\n");
            exit(EXIT_FAILURE);
        }
        
    pccf = (double *)malloc((size_t) (2*maxlag+1)* sizeof(double));
    if (pccf == NULL) 
        {
            fprintf(stderr,"\nMemory allocation error!\n");
            exit(EXIT_FAILURE);
        }
        
    printf("Traces allocated.\n");

    /* Determing the maximum 'usable' window and if there are enough samples*/
    len=ndat-2*maxlag;
   
    /* If there aren't enough samples, go back.*/
    if (len <= 0)
        {
            printf("Cannot cross-correlate as the number of samples is not sufficient for the maximum lag. Choose a smaller maximum lag and try again.");
            exit(EXIT_FAILURE);
        }	
    
    
    /* If there are enough samples: Get going! */
    else
    {
        /* For now, no mean removal and no normalization. This should be taken care off already by preprocessing! */
        /* Same for the quality control of the trace (can reintroduce it one day)*/
        
        /*----------------------------------------------------------------------------------------------------------------------
        #-First, copy the traces to a new location 
        #-We want unit offset arrays because this is what the FFT code works with
        -----------------------------------------------------------------------------------------------------------------------*/
       for (a=1;a<=ndat;a++)
            {
                tra1[a]=u[a];
                tra2[a]=u[a];
                }
        
     printf("Example trace copied to new traces.\n");  
        
        /*-----------------------------------------------------------------------------------------------------------------------
        #-Second, need to obtain analytic signal of the traces
        #-This is handled by asignal subroutine
        #-The resulting analytic signals are written to atra1,atra2; these arrays contain 2*ndat numbers according
        #-to the scheme: Re1 Im1 Re2 Im1......Re ndat Im ndat
        -----------------------------------------------------------------------------------------------------------------------*/
        
        asignal(tra1,tra2, atra1, atra2, ndat);
                printf("Analytic signal determined.\n"); 
               
                /*for now, since this is less complicated: Split atra arrays up into their real and imaginary parts. 
                Later on, this can be replaced by one single array, which will be more convenient but less readable.*/
        for (a=1;a<=ndat;a++){
            
            re1[a-1]=atra1[2*a-1];
            im1[a-1]=atra1[2*a];
            re2[a-1]=atra2[2*a-1];
            im2[a-1]=atra2[2*a];
        }
                
        printf("Obtained real and imaginary parts of analytic signal.\n"); } 
        /*----------------------------------------------------------------------------------------------------------------------
        #-Third, obtain the phase cross correlation from the analytic signal real and im. parts
        -----------------------------------------------------------------------------------------------------------------------*/     
        for (a=-maxlag;a<=maxlag;a++){
            i1=a+maxlag;
            i2=a+ndat-maxlag;
            pccf[i1]=0.0;
            
            for (b=i1;b<=i2;b++)
            {
                abs1=sqrt(re1[b]*re1[b]+im1[b]*im1[b]);
                re1[b]/=abs1;
                im1[b]/=abs1;
                
                abs2=sqrt(re2[b-a]*re2[b-a]+im2[b-a]*im2[b-a]);
                re2[b-a]/=abs2;
                im2[b-a]/=abs2;
                
            pccf[i1]+=sqrt((re1[b]+re2[b-a])*(re1[b]+re2[b-a])+(im1[b]+im2[b-a])*(im1[b]+im2[b-a]))-sqrt((re1[b]-re2[b-a])*(re1[b]-re2[b-a])+(im1[b]-im2[b-a])*(im1[b]-im2[b-a]));
            }       
               
            pccf[i1]/=(2*len);  
        }
        
    printf("Phase cross correlation calculated.\n");
        

    fh_pcc=fopen("test_pcc.txt","w");
        for (a=0;a<=2*maxlag;a++){
            fprintf(fh_pcc,"%g\n",pccf[a]);
        }
    fclose(fh_pcc);
       
    free((char *)tra1);
    free((char *)tra2);
    free((char *)atra1);
    free((char *)atra2);
} 




/*--------------------------------------------------------------------
# Obtain analytic signal
--------------------------------------------------------------------*/
/* ----------------------------------------------------------------------------------------- */
/* ----------------------------------------------------------------------------------------- */
/* This is a short script to get the analytic signal in C. FFT is based on the code provided in
'Numerical Recipes in C'. Precision is down to float again (had it double, but overlying 
cross correlation routine uses float.) The FT routine for transforming two 
real arrays at once is used. For IFFT the single-signal routine has to be used as the 
resulting array is complex.*/
/* ----------------------------------------------------------------------------------------- */
/* ----------------------------------------------------------------------------------------- */

void asignal(float *tra1,float *tra2, float *atra1, float *atra2, int ndata){
    /*Declaration*/
    int i;
    /*
    /* The analytic signal can be obtained by calculating IFFT(FFT(signal)*sign(signal))
    I am not sure here: This should be 2*FFT for the first half of 
    the FFT and 0 for the second half? */
    printf("Entering asignal routine.\n");
    tworealfft(tra1,tra2,atra1,atra2,ndata);
    
   for (i=1;i<=ndata;i++){
        atra1[i]*=2.0;
        atra2[i]*=2.0;
        
    }
    
    for (i=ndata+1;i<=2*ndata;i++){
        atra1[i]=0.0;
        atra2[i]=0.0;
        
    }
    

    /*Inverse FT yields the analytic signal*/
   fourier1(atra1,ndata,-1);
    fourier1(atra2,ndata,-1);
}



