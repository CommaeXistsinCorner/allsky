#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "spectra_sim.cc"
#include <TCanvas.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH2D.h>
#include <TMath.h>
#include <TF1.h>
#include <TH2.h>
#include <TString.h>
#include <iostream>
#include <fstream>

int main(int argc, char *argv[])
{
  char  filename[200];
  FILE  *fp_in,*fp_log;

  if(argc<2) { printf("%s  output(d.root)\n", argv[0]); exit(0); }

  strcpy(filename,argv[1]);
  if((fp_in=fopen(strcat(argv[1],"d.root"),"rb"))==NULL) { printf("cannot open input file\n"); exit(0);}
  strcpy(argv[1],filename);
  if((fp_log=fopen(strcat(argv[1],".log"),"a+"))==NULL) { printf("cannot open input file\n"); exit(0);}

  fprintf(fp_log,"\n\n////////////// identify_sources.exe ///////////////\n");

  strcpy(argv[1],filename);
  TFile file(strcat(argv[1],"d.root"));
  TH1D *h = (TH1D*)file.Get("h8020");

  int PS_num = sizeof(PS_TEVCAT)/(sizeof(double)*7);
  double ra,dec;
  for(int i_PS=0; i_PS<PS_num; i_PS++)
  {
    double ra = PS_TEVCAT[i_PS][0];
    double dec = PS_TEVCAT[i_PS][1];
    int n_global = h->FindBin(ra,dec);
    int imid,jmid,kmid;
    h->GetBinXYZ(n_global,imid,jmid,kmid);
    int imin = imid-4;
    int imax = imid+4;
    int jmin = jmid-4;
    int jmax = jmid+4;
    double sigma_max=0.;
    for(int i = imin; i <= imax; i++)
    {
      for(int j = jmin; j <= jmax; j++)
      {
        double sigma = h->GetBinContent(i,j);
        if(sigma_max < sigma) sigma_max = sigma;
      }
    }
    printf("( i_PS , sigma_max ) = ( %3d , %8.8f )\n", i_PS, sigma_max);
    fprintf(fp_log, "( i_PS , sigma_max ) = ( %3d , %8.8f )\n", i_PS, sigma_max);
  }

  file.Close();
  fclose(fp_in);
  fclose(fp_log);
  return 0;
}
