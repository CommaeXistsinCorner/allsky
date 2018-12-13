#include         <math.h>
#include        <stdlib.h>
#include        <stdio.h>
#include        "TROOT.h"
#include        "TFile.h"
#include        "TTree.h"
#include        "TH2.h"
#include        "TH1.h"
#include        "ana_wid_sim.h"
//#include        "point_source_wcda.cc"
#include        "spectra_sim.cc"
#include        "Astro.c"
double normlize(int norm_f,int clean_flag);
//void htoe2(double zen,double azim,double lst,double *ra,double *dec);

static float   n[NZEAZ];
static int     p_ra[NZEAZ],p_dec[NZEAZ];  //for point
//static float   n[NLST][NZEAZ];
//static double  n_a[NLST][NZE],n_b;
static double  beta[NRA][NDEC];
static double  a[NRA][NDEC],b[NRA][NDEC];
double  ki2,kia,kib,ki20=0.,ki21=0.;
double sig;
//for smooth

int main(int argc,char *argv[])
{   
  FILE          *fp_out,*fp_in;
  FILE          *fp_na,*fp_log;
  int           isr;
  int           i, j, k, k0, t;
  double        zen,azim,ra,dec,lst;
  int           i_ra,j_dec;  
  double        max=0.;
  int           iter=25;
  double        n_b=0.0;
  float         tmp_fp_na;
  int           ismon=1;
  char          filename[200];

  if(argc<2) { 
    printf("%s  input(c.dat)  [iter=25]  [1 for monitor mode]\n",argv[0]);
    exit(0);
  }  
  if(argc>=3) {iter=atoi(argv[2]);printf("iter=%d\n",iter);} 
  if(argc>=4)  ismon=atoi(argv[3]);

  strcpy(filename,argv[1]);
  if((fp_in=fopen(strcat(argv[1],"c.dat"),"rb"))==NULL) {printf("cannot open input file\n"); exit(0);}
  strcpy(argv[1],filename);
  if((fp_out=fopen(strcat(argv[1],"d.dat"),"wb+"))==NULL) {printf("cannot create output file\n"); exit(0);}
  strcpy(argv[1],filename);
  if((fp_log=fopen(strcat(argv[1],".log"),"a+"))==NULL)  {printf("cannot open log file\n");exit(0);}
  strcpy(argv[1],filename);
  if((fp_na=fopen(strcat(argv[1],"_tmp_na"),"wb+"))==NULL)  {printf("cannot create temp file\n");exit(0);}

  fprintf(fp_log,"\n\n//////////ana_wid_sim.exe//////////\n");

  for(int tmp_i=0,tmp_k=0; tmp_i<NZE; tmp_i++)  {
    NAZ0[tmp_i] = ceil( 360.0 * sin((0.5+tmp_i)*hori_sys_width*PI/180.0) / hori_sys_width );
    NAZ[tmp_i] = tmp_k;
    tmp_k += NAZ0[tmp_i];
  }

  ////////////intitial the data the point ///////////////////////////////////////
  for(i=0;i<NRA;i++) { 
    for(j=0;j<NDEC;j++) {
      beta[i][j]=1.;
      a[i][j]=0.;
      b[i][j]=0.;
    }
  }
  ////initial the point int isr=0 
  k=0;lst=0.5*equa_sys_width;
  for(i=0;i<NZE;i++)   {
    zen=(i+0.5)*hori_sys_width; 
    for(j=0;j<NAZ0[i];j++)   	{ 
      azim=(j+0.5)*360./NAZ0[i];
      htoe2(zen,azim,lst,&ra,&dec);
      p_ra[k]=int(ra/equa_sys_width);p_dec[k++]=int((dec-DEC_MIN)/equa_sys_width);
    }
  }
  if(ismon==1)  printf("initiation done\n");

  /////////////////////////iteraion work////////////////////////////////////////
  for(t=1;t<iter+1;t++)   { 
    ki2=0.;
    ////read from the data 
    fseek(fp_in,0,SEEK_SET);
    fseek(fp_na,0,SEEK_SET);
    for(tmp_fp_na=0.0,isr=0;isr<NLST;isr++) {
      for(i=0;i<NZE;i++) {
        fwrite(&tmp_fp_na,sizeof(tmp_fp_na),1,fp_na);
      }
    }

    //iteration
    for(isr=0;isr<NLST;isr++)     {

      if(ismon==1) printf("%d/%d    chi square calculating!   %3d%% \r",t,iter,int(isr*100/NLST));
      k=0;
      lst=(isr+0.5)*360.0/NLST;
      
//fprintf(stderr,"ok\n");
      for(i=0;i<NZE;i++) 	         {
        tmp_fp_na=0.0;
	k0=k; 
	fread(n,sizeof(float)*NAZ0[i],1,fp_in);

        fseek(fp_na,(isr*NZE+i)*sizeof(tmp_fp_na),SEEK_SET);
        fread(&tmp_fp_na,sizeof(tmp_fp_na),1,fp_na);
	for(j=0;j<NAZ0[i];j++)       { 
	  i_ra=(p_ra[k]+isr)<NLST?(p_ra[k]+isr):(p_ra[k]+isr-NLST);
	  j_dec=p_dec[k];
          tmp_fp_na += n[j]*beta[i_ra][j_dec];
	  k++;
//	  if(beta[i_ra][j_dec]==0) {printf("%d %d %d %f\n",NDEC,i_ra,j_dec,beta[i_ra][j_dec]);exit(0);}
	}	
        fseek(fp_na,-1*sizeof(tmp_fp_na),SEEK_CUR);
        fwrite(&tmp_fp_na,sizeof(tmp_fp_na),1,fp_na);	

	k=k0;
	for(j=0;j<NAZ0[i];j++)       {
	  i_ra=(p_ra[k]+isr)<NLST?(p_ra[k]+isr):(p_ra[k]+isr-NLST);
	  j_dec=p_dec[k];
	  n_b=(tmp_fp_na-n[j]*beta[i_ra][j_dec])/(NAZ0[i]-1); 
	  b[i_ra][j_dec]+=n_b;
	  if(t==1) a[i_ra][j_dec]+=n[j];
	  if(n[j]!=0)    {
	    kia=pow((n[j]*beta[i_ra][j_dec]-n_b),2);
	    kib=n[j]*(beta[i_ra][j_dec]*beta[i_ra][j_dec]);
	    ki2+=kia/kib;
//     printf("Now t=%d\t  ki2=%.20f\n",t,ki2);
//printf("%f %f %f %f %f %d %d\n",kia,kib,ki2,n[j],beta[i_ra][j_dec],i_ra,j_dec);
//printf("%d %d %d\n",isr,i,j);
            }
	  k++;
	}    
      }
//fprintf(stderr,"ok\n");
//fprintf(stderr,"%d %d %d %d\n",isr,i,i_ra,j_dec);
  }
    
//fprintf(stderr,"%d %d %d %d\n",isr,i,i_ra,j_dec);
    // printf("at the end t=%d\t %f\t%.20f\n",t,beta[358][40],ki2);
    //  if(t>2&&fabs(ki20-ki21)<fabs(ki21-ki2)) break;
    ki20=ki21;
    ki21=ki2;
    printf("Now t=%d\t  ki2=%.20f\n",t,ki2);
    fprintf(fp_log,"Now t=%d\t  ki2=%.20f\n",t,ki2);
//printf("%f %f %f\n",kia,kib,ki2);
    normlize(1,t<iter?0:1);
//    if(t==iter) break;
    if(t>9&&fabs(ki20-ki21)/ki20<fabs(ki21-ki2)/ki21) break;
  }
  /////////////////////////////////////////////////////////

  //normlize(1);
    if(ismon==1) printf("begin to fill the root \n");  
    strcpy(argv[1],filename);
    TFile hfile(strcat(argv[1],"c.root"), "RECREATE");
    TH2D *h8021 = new TH2D("h8021","Intens",NRA,0.,360.,NDEC,DEC_MIN,DEC_MAX);
    TH2D *h8020 = new TH2D("h8020","Signi",NRA,0.,360.,NDEC,DEC_MIN,DEC_MAX);
    TH2D *hnon = new TH2D("hnon","N_on",NRA,0.,360.,NDEC,DEC_MIN,DEC_MAX);
    TH2D *hnoff = new TH2D("hnoff","N_off",NRA,0.,360.,NDEC,DEC_MIN,DEC_MAX);
    TH1D *h10 = new TH1D("h10","Sigin. dis",140,-30.,30.);
    TH1D *hra = new TH1D("hra_in","RA dis.",NRA,0.,360.);
    TH1D *hdec = new TH1D("hdec_in","DEC dis.",NDEC,DEC_MIN,DEC_MAX);
    TH1D *hsig = new TH1D("hsig_in","Signi",3000,-500,500);
    normlize(3,1);
    double ps_sig,ps_bg,sig_max=0.;
    double ps_ra,ps_dec;
      for(i=0;i<NRA;i++)	       {
        if(ismon==1) printf("RA-DEC map!   %3d%% \r",int(i*100/NRA));
        for(j=0;j<NDEC;j++)             {
    //		   printf("%f\n",a[j][i]);
          fwrite(&a[i][j],8,1,fp_out);
          fwrite(&b[i][j],8,1,fp_out);
          max+=a[i][j];
          sig=(1-beta[i][j])*sqrt(pow(a[i][j],3)/pow(b[i][j],2));
          if(isnan(sig)!=0||isinf(sig)!=0) sig=0;
          if(sig_max < sig && fabs((i+0.5)*equa_sys_width-CRAB[0])<2 && fabs((j+0.5)*equa_sys_width+DEC_MIN-CRAB[1])<2)
          {
            ps_ra = (i+0.5)*equa_sys_width;
            ps_dec = (j+0.5)*equa_sys_width+DEC_MIN;
            sig_max = sig;
            ps_sig = a[i][j];
            ps_bg = b[i][j];
          }
//             sig_max=sig_max>sig?sig_max:sig;
//printf("%f %f\n",sig_max,sig);
          h8020->Fill((i+0.5)*equa_sys_width,(j+0.5)*equa_sys_width+DEC_MIN,sig); 
	  h8021->Fill((i+0.5)*equa_sys_width,(j+0.5)*equa_sys_width+DEC_MIN,b[i][j]==0?1:a[i][j]/b[i][j]);
          hnon->Fill((i+0.5)*equa_sys_width,(j+0.5)*equa_sys_width+DEC_MIN,a[i][j]-b[i][j]);
          hnoff->Fill((i+0.5)*equa_sys_width,(j+0.5)*equa_sys_width+DEC_MIN,b[i][j]);
	  h10->Fill(sig);
	  hsig->Fill(sig);
          hra->Fill((i+0.5)*equa_sys_width,a[i][j]);
          hdec->Fill((j+0.5)*equa_sys_width+DEC_MIN,a[i][j]);		        
        }
      }
        printf("\nsig_max=%f\n ps_sig=%f\n ps_bg=%f\n",sig_max,ps_sig,ps_bg);
        printf("sig_max(ra,dec) = ( %3f , %3f )\n",ps_ra,ps_dec);
printf("max=%f\n",max);  
        fprintf(fp_log,"sig_max=%f\n ps_sig=%f\n ps_bg=%f\n",sig_max,ps_sig,ps_bg);
        fprintf(fp_log,"sig_max(ra,dec) = ( %3f , %3f )\n",ps_ra,ps_dec);
	      fprintf(fp_log,"max=%f\n",max);  
        hfile.Write();
	      hfile.Close();

        fclose(fp_in);
        fclose(fp_log);
        fclose(fp_na);
        fclose(fp_out);

        if(ismon==0)
        {
          strcpy(argv[1],filename);
          remove(strcat(argv[1],"_tmp_na"));
          strcpy(argv[1],filename);
//           remove(strcat(argv[1],"c.dat"));
        }
 	 return 0;
}


double normlize(int norm_f,int clean_flag)
{
  int i,j;
  ///////////////////////global normal//////////////////
  if(norm_f==0) {
    double I_c0=0.,I_c0_err=0.;
    
    for(i=0;i<NDEC;i++)     {
      for(j=0;j<NRA;j++)     {
	if(a[j][i]==0) beta[j][i]=1.;
	else   beta[j][i]=b[j][i]/a[j][i];

	I_c0+=1./beta[j][i]*a[j][i]*beta[j][i]*beta[j][i];
	I_c0_err+=a[j][i]*beta[j][i]*beta[j][i];
      }
    }
    I_c0/=I_c0_err;
    for(i=0;i<NDEC;i++)     {
      for(j=0;j<NRA;j++)     {
	beta[j][i]*= I_c0;
      }
    }

    printf("normlize=%f \n",I_c0);
  }
  
  else if(norm_f==1) {
    double beta_c=0.;
    for(i=0;i<NDEC;i++)    {
      for(j=0;j<NRA;j++)       {
	if(a[j][i]==0) beta[j][i]=1.;
	else   beta[j][i]=b[j][i]/a[j][i];
	beta_c+=1./beta[j][i];
      }
    }
    beta_c/=(NRA*NDEC);
    for(i=0;i<NDEC;i++)    {
      for(j=0;j<NRA;j++)         {
	beta[j][i]=beta[j][i]*beta_c;
      }
    }
  }
  
///////////////////////Dec. normal//////////////////
  else if(norm_f==2) {
    double I_c[NDEC],I_c_err[NDEC];
    for(i=0;i<NDEC;i++)    {
      I_c[i]=0.;I_c_err[i]=0.;
      for(j=0;j<NRA;j++)       {
	if(a[j][i]==0) beta[j][i]=1.;
	else   beta[j][i]=b[j][i]/a[j][i];
	I_c[i]+=1./beta[j][i]*a[j][i]*beta[j][i]*beta[j][i];
	I_c_err[i]+=a[j][i]*beta[j][i]*beta[j][i];
	
      }
      I_c[i]/=I_c_err[i];
    }
    for(i=0;i<NDEC;i++)    {
      for(j=0;j<NRA;j++)         {
	beta[j][i]*=I_c[i];
      }
    }
  }

  else if(norm_f==3) {
    double beta_c[NDEC];
    for(i=0;i<NDEC;i++)    {
      beta_c[i]=0.;
      for(j=0;j<NRA;j++)       {
	if(a[j][i]==0) beta[j][i]=1.;
	else   beta[j][i]=b[j][i]/a[j][i];
	beta_c[i]+=1./beta[j][i];
      }
    }
    
    //beta_c=beta_c/(NRA*NDEC);
    for(i=0;i<NDEC;i++)    {
      beta_c[i]/=NRA;
      for(j=0;j<NRA;j++)         {
	beta[j][i]=beta[j][i]*beta_c[i];
      }
    }
  }
  

  //////////clear////////////////////
  if(clean_flag==0) {
    for(i=0;i<NDEC;i++)       {
      for(j=0;j<NRA;j++)        {
	b[j][i]=0.;
      }
    } 
  }
  else if (clean_flag==1) {
    for(i=0;i<NDEC;i++)       {
      for(j=0;j<NRA;j++)        {
	b[j][i]=beta[j][i]*a[j][i];
      }
    } 
  }

  return 0;
}

