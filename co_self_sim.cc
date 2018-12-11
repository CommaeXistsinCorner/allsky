#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "ana_wid_sim.h"
#include <iostream>
#include <fstream>
using namespace std;

float     n0;
float     n;
double    eff[NZEAZ];

int main(int argc,char *argv[])
{
//   fstream fout("out.txt",ios::out|ios::binary);

  FILE                     *fp_in,*fp_out;
  FILE                     *fp_log;
  int                       i,j,isr;
  int                      k,k0; 
  double                    hour;
  double                    max=0.;
  double                    n_b=0.;
  int                      ismon=1;
  char                     filename[200];
   
  if(argc<2) { printf("%s  input(.dat)  [1 for monitor mode]\n",argv[0]);exit(0);}
  if(argc>=3)  ismon=atoi(argv[2]);
  strcpy(filename,argv[1]);  
  if((fp_in=fopen(strcat(argv[1],".dat"),"rb"))==NULL) {printf("cannot open input file\n"); exit(0);}
  strcpy(argv[1],filename);
  if((fp_out=fopen(strcat(argv[1],"c.dat"),"wb"))==NULL) {printf("cannot create output file\n"); exit(0);}
  strcpy(argv[1],filename);
  if((fp_log=fopen(strcat(argv[1],".log"),"a+"))==NULL)  {printf("cannot create log file\n");exit(0);}

  fprintf(fp_log,"\n\n//////////co_self.exe//////////\n");

  for(int tmp_i=0,tmp_k=0; tmp_i<NZE; tmp_i++)  {
    NAZ0[tmp_i] = ceil( 360.0 * sin((0.5+tmp_i)*hori_sys_width*PI/180.0) / hori_sys_width );
    NAZ[tmp_i] = tmp_k;
    tmp_k += NAZ0[tmp_i];
  }

  for(k=0;k<NZEAZ;k++) eff[k]=0.;

  for(isr=0;isr<NLST;isr++) {
    if(ismon==1) printf("reading data!   %3d%% \r",int(isr*100/NLST));
    for(k=0;k<NZEAZ;k++)      {
      fread(&n0,sizeof(n0),1,fp_in);
//if(!(n0<1e-5&&n0>-1e-5)) printf("%f\n",n0);
      eff[k]+=n0;
    }
  }
    
  k=0;
  for(i=0;i<NZE;i++)      { 
    if(ismon==1) printf("calibrating!   %3d%% \r",int(i*100/NZE));
    n_b=0.;
    k0=k;
    for(j=0;j<NAZ0[i];j++) {
      n_b+=eff[k];
      k++;
    }
    n_b/=NAZ0[i];
    k=k0;
    for(j=0;j<NAZ0[i];j++) {
    if((eff[k]<1e-5)||(n_b<1e-5)) {eff[k]=1;}
    else{
      eff[k]=eff[k]/n_b;
      }
      k++;
    }
    //  printf("\n");
  }
    // exit(0);
    ///////////////////////eff calculate////////////////////////
   fseek(fp_in,0L,0);

   for(isr=0;isr<NLST;isr++) {  
     if(ismon==1) printf("correcting data!   %3d%% \r",int(isr*100/NLST));
     for(k=0;k<NZEAZ;k++)      {
       fread(&n0,sizeof(n0),1,fp_in);
       n=n0/eff[k];
       if(n<0.0) n=0.0;
       if(isnan(n)!=0||isinf(n)!=0) n=0.0;
//       	  fprintf(stderr,"%d %f %f\n",n0,n,phi_cor(&iflag,&theta,&phi));
       fwrite(&n,sizeof(n),1,fp_out);
//       cout<<n<<endl;
       max+=n;
     }                   			          
   }

   printf("\n%f\tmax = %f\n",n,max);
   fprintf(fp_log,"%f\tmax = %f\n",n,max);  max = 0.;

   fclose(fp_log);
   fclose(fp_in);
   fclose(fp_out);
   if(ismon==0)  {strcpy(argv[1],filename);remove(strcat(argv[1],".dat"));}
   return 0;
}
