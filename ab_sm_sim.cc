#include "stdio.h"
#include "math.h"
#include <stdlib.h>
#include "TFile.h"
#include "TH2D.h"
#include "TH1D.h"
#include "ana_wid_sim.h"
#include "Astro.c"
#include "detector_sim.cc"
#include "spectra_sim.cc"
//double intens_b[ra_bin],intens_a[ra_bin];
void sig_correct(int i0, int j0, double r_smooth, double grid_width, double *sig, double *intens, double *Nsig, double *Nbg, int *Nsmthbin);
void wcda_ang_maxsigma(double *theta_maxsigma, double *maxsigma, double *keepratio_ga, double *keepratio_cr, double *num_ga, double *num_cr);
//void htoe2(double zen,double azim,double lst,double *ra,double *dec);
double a[NRA][NDEC], b[NRA][NDEC], beta[NRA][NDEC];
double intens_b[NRA], intens_a[NRA];
double signi[NRA][NDEC], intens[NRA][NDEC];
double a_tmp, b_tmp;
double num_sig, num_bg, num_sig_max, num_bg_max;
int num_bin, num_bin_max;
int smth[200000][2];
int smth_num;
int normalize(int norm_f);
double wcda_ang(double ep)
{
  int i = 0;
  double theta;

  ep = pow(10.0, ep - 12.0);

  if (ep <= exb[0])
    theta = angoptb[0];
  else if (ep >= exb[nvvb - 1])
    theta = angoptb[nvvb - 1];
  else
  {
    while (ep > exb[i])
      i++;
    theta = (angoptb[i] - angoptb[i - 1]) / (exb[i] - exb[i - 1]) * (ep - exb[i - 1]) + angoptb[i - 1];
  }
  if (theta < 0 || isnan(theta) != 0)
    theta = 0.;
  return theta;
}

int main(int argc, char *argv[])
{
  FILE *fp_in;
  FILE *fp_log;
  int isr;
  int i, j, k, k0, t;
  double zen, azim, ra, dec, lst;
  double sig, max0 = 0.;
  double sig_max;
  double ps_ra, ps_dec;
  double max = 0.;
  double reso1, reso2 = 5;
  double ediff;
  int t_reuse = 1;
  int ismon = 1;
  char filename[200];
  double maxsigma, keepratio_ga, keepratio_cr, num_ga, num_cr;

  if (argc < 3)
  {
    printf("%s  input(d.dat)  ediff  [reso1]  [reso2]  [t_reuse]  [1 for monitor mode]  (use log10(ediff))  set Reso1=0 for default\n", argv[0]);
    exit(0);
  }
  ediff = atof(argv[2]);
  //  reso1=1.5852*wcda_ang(ediff);
  wcda_ang_maxsigma(&reso1, &maxsigma, &keepratio_ga, &keepratio_cr, &num_ga, &num_cr);
  if (argc >= 4 && atof(argv[3]) > 1e-4)
    reso1 = atof(argv[3]);
  if (argc >= 5)
    reso2 = atof(argv[4]);
  if (argc >= 6)
    t_reuse = atoi(argv[5]);
  if (argc == 7)
    ismon = atoi(argv[6]);
  //  reso=atof(argv[2]);
  strcpy(filename, argv[1]);
  if ((fp_in = fopen(strcat(argv[1], "d.dat"), "rb")) == NULL)
  {
    printf("cannot open input file.\n");
    exit(0);
  }
  strcpy(argv[1], filename);
  if ((fp_log = fopen(strcat(argv[1], ".log"), "a+")) == NULL)
  {
    printf("cannot open log file.\n");
    exit(0);
  }

  fprintf(fp_log, "\n\n//////////ab_sm_sim.exe//////////\n");
  printf("reso %g %g     t_reuse %d\n", reso1, reso2, t_reuse);
  fprintf(fp_log, "reso %g %g     t_reuse %d\n", reso1, reso2, t_reuse);

  for (i = 0; i < NRA; i++)
  {
    for (j = 0; j < NDEC; j++)
    {
      a[i][j] = 0.;
      b[i][j] = 0.;
    }
  }

  for (t = 0; t < t_reuse; t++)
  {
    fseek(fp_in, 0, SEEK_SET);
    for (i = 0; i < NRA; i++)
    {
      if (ismon == 1)
        printf("%d/%d    initiating!   %3d%% \r", t, t_reuse, int(i * 100 / NRA));
      for (j = 0; j < NDEC; j++)
      {
        fread(&a_tmp, 8, 1, fp_in);
        fread(&b_tmp, 8, 1, fp_in);
        max = max + a_tmp;
        a[i][j] += a_tmp;
        b[i][j] += b_tmp;
      }
    }
  }
  fclose(fp_in);
  printf("max= %f\n", max);
  fprintf(fp_log, "max= %f\n", max);
  for (i = 0; i < NRA; i++)
  {
    for (j = 0; j < NDEC; j++)
    {
      beta[i][j] = b[i][j] / a[i][j];
      //printf("%f\n",beta[i][j]);
    }
  }

  normlize(1);
  strcpy(argv[1], filename);
  TFile hfile(strcat(argv[1], "d.root"), "RECREATE");
  if (ismon == 1)
    printf("ROOT file part 1 is preparing \n");
  TH2D *h8021 = new TH2D("h8021", "Intens", NRA, 0., 360., NDEC, DEC_MIN, DEC_MAX);
  TH2D *h8020 = new TH2D("h8020", "Signi", NRA, 0., 360., NDEC, DEC_MIN, DEC_MAX);
  //  TH2D *hnon = new TH2D("hnon","N_on",NRA,0.,360.,NDEC,DEC_MIN,DEC_MAX);
  //  TH2D *hnoff = new TH2D("hnon","N_off",NRA,0.,360.,NDEC,DEC_MIN,DEC_MAX);
  TH1D *h10 = new TH1D("h10", "Sigin. dis", 140, -30., 30.);
  TH1D *hra = new TH1D("hra", "RA. dis", NRA, 0., 360.);
  TH1D *hdec = new TH1D("hdec", "DEC. dis", NDEC, DEC_MIN, DEC_MAX);
  TH1D *hsig = new TH1D("hsig_in", "Signi", 3000, -500, 500);

  sig_max = 0.;
  for (i = 0; i < NRA; i++)
  {
    intens_b[i] = 0.;
    intens_a[i] = 0.;
  }
  for (i = 0; i < NRA; i++)
  {
    for (j = 0; j < NDEC; j++)
    {
      intens_b[i] += b[i][j];
      intens_a[i] += a[i][j];
    }
  }

  for (i = 0; i < NRA; i++)
  {
    hra->SetBinContent(i + 1, intens_a[i] / intens_b[i]);
    hra->SetBinError(i + 1, sqrt(intens_a[i]) / intens_b[i]);
  }

  for (j = 0; j < NDEC; j++)
  { // for reusing the selected grid when DEC is same.
    if (ismon == 1)
      printf("point source smooth set!   %3d%%\r", int(j * 100 / NDEC));
    for (i = 0; i < NRA; i++)
    {
      hdec->Fill((j + 0.5) * equa_sys_width + DEC_MIN, a[i][j]);
      sig_correct(i, j, reso1, equa_sys_width, &signi[i][j], &intens[i][j], &num_sig, &num_bg, &num_bin);
      if (isnan(signi[i][j]) != 0 || isinf(signi[i][j]) != 0)
        signi[i][j] = 0;
      if (sig_max < signi[i][j] && fabs((i + 0.5) * equa_sys_width - CRAB[0]) < 2 && fabs((j + 0.5) * equa_sys_width + DEC_MIN - CRAB[1]) < 2)
      {
        ps_ra = (i + 0.5) * equa_sys_width;
        ps_dec = (j + 0.5) * equa_sys_width + DEC_MIN;
        sig_max = signi[i][j];
        num_sig_max = num_sig;
        num_bg_max = num_bg;
        num_bin_max = num_bin;
      }
      h8020->Fill((i + 0.5) * equa_sys_width, (j + 0.5) * equa_sys_width + DEC_MIN, signi[i][j]);
      h10->Fill(signi[i][j]);
      hsig->Fill(signi[i][j]);
      h8021->SetBinContent(i + 1, j + 1, intens[i][j]);
      h8021->SetBinError(i + 1, j + 1, (intens[i][j] - 1) / signi[i][j]);
    }
    //printf("%d\n",smth_num);
  }
  printf("\nsig_max = %f    (ra,dec) = ( %3f , %3f )\n", sig_max, ps_ra, ps_dec);
  printf("(Nsig,Nbg) = ( %3f , %3f )\n", num_sig_max, num_bg_max);
  printf("Nsmthbin = %d\n", num_bin_max);

  fprintf(fp_log, "sig_max = %f    (ra,dec) = ( %3f , %3f )\n", sig_max, ps_ra, ps_dec);
  fprintf(fp_log, "(Nsig,Nbg) = ( %3f , %3f )\n", num_sig_max, num_bg_max);
  fprintf(fp_log, "Nsmthbin = %d\n", num_bin_max);
  ///ADD LARGE SMOOTH//////////////////////////////////////////////////////
  if (ismon == 1)
    printf("ROOT file part 2 is preparing \n");
  TH2D *h7021 = new TH2D("h7021", "Intens", NRA, 0., 360., NDEC, DEC_MIN, DEC_MAX);
  TH2D *h7020 = new TH2D("h7020", "Signi", NRA, 0., 360., NDEC, DEC_MIN, DEC_MAX);
  TH2D *h5s = new TH2D("h5s", "Signi", NRA, 0., 360., NDEC, DEC_MIN, DEC_MAX);
  for (j = 0; j < NDEC; j++)
  {
    if (ismon == 1)
      printf("large area smooth set!   %3d%%\r", int(j * 100 / NDEC));
    for (i = 0; i < NRA; i++)
    {
      sig_correct(i, j, reso2, equa_sys_width, &signi[i][j], &intens[i][j], &num_sig, &num_bg, &num_bin);
      h7020->Fill((i + 0.5) * equa_sys_width, (j + 0.5) * equa_sys_width + DEC_MIN, signi[i][j]);
      h7021->Fill((i + 0.5) * equa_sys_width, (j + 0.5) * equa_sys_width + DEC_MIN, intens[i][j]);
      if (fabs(intens[i][j] - 1) / signi[i][j] < 0.0002)
        h5s->Fill((i + 0.5) * equa_sys_width, (j + 0.5) * equa_sys_width + DEC_MIN, intens[i][j]);
    }
  }

  ///////////////mid structure////////////////////////////////////////////
  if (ismon == 1)
    printf("ROOT file part 3 is preparing \n");
  TH2D *h6021 = new TH2D("h6021", "Intens", NRA, 0., 360., NDEC, DEC_MIN, DEC_MAX);
  TH2D *h6020 = new TH2D("h6020", "Signi", NRA, 0., 360., NDEC, DEC_MIN, DEC_MAX);
  double intensi0, err0;
  for (i = 0; i < NRA; i++)
  {
    if (ismon == 1)
      printf("ps - area smooth  set!   %3d%%\r", int(i * 100 / NRA));
    for (j = 0; j < NDEC; j++)
    {
      intensi0 = h8021->GetBinContent(i + 1, j + 1) - h7021->GetBinContent(i + 1, j + 1);
      err0 = (h8021->GetBinContent(i + 1, j + 1) - 1) / h8020->GetBinContent(i + 1, j + 1);
      //  printf("intenis %g\n",err0);
      h6021->SetBinContent(i + 1, j + 1, intensi0 + 1);
      h6020->SetBinContent(i + 1, j + 1, intensi0 / err0);
      h6021->SetBinError(i + 1, j + 1, err0);
    }
  }

  hfile.Write();
  hfile.Close();
  fclose(fp_log);
  if (ismon == 0)
  {
    strcpy(argv[1], filename);
    remove(strcat(argv[1], "d.dat"));
  }
  return 0;
}

void sig_correct(int i0, int j0, double r_smooth, double grid_width, double *sig, double *intens, double *Nsig, double *Nbg, int *Nsmthbin)
{
  int i, j, k;
  double Nsig_tot = 0., Nbg_tot = 0.;
  if (i0 == 0)
  {
    double ra0, dec0, ra, dec;
    smth_num = 0;
    ra0 = (i0 + 0.5) * grid_width;
    dec0 = (j0 + 0.5) * grid_width + DEC_MIN;
    for (i = 0; i < NRA; i++)
    {
      ra = (i + 0.5) * grid_width;
      for (j = 0; j < NDEC; j++)
      {
        dec = (j + 0.5) * grid_width + DEC_MIN;
        if (r_smooth < distance_equatorial(ra, dec, ra0, dec0))
          continue;
        smth[smth_num][0] = i;
        smth[smth_num][1] = j;
        smth_num++;
      }
    }
  }
  for (k = 0; k < smth_num; k++)
  {
    i = smth[k][0] + i0 < NRA ? (smth[k][0] + i0) : (smth[k][0] + i0 - NRA);
    j = smth[k][1];
    Nsig_tot += a[i][j] - b[i][j];
    Nbg_tot += b[i][j];
  }
  *Nsig = Nsig_tot;
  *Nbg = Nbg_tot;
  *Nsmthbin = smth_num;
  *sig = Nsig_tot / sqrt(Nbg_tot);
  *intens = 1. + Nsig_tot / Nbg_tot;
}

int normalize(int norm_f)
{
  int i, j;
  ///////////////////////global normal//////////////////
  if (norm_f == 0)
  {
    double I_c0 = 0., I_c0_err = 0.;

    for (i = 0; i < NDEC; i++)
    {
      for (j = 0; j < NRA; j++)
      {
        if (a[j][i] == 0)
          beta[j][i] = 1.;
        else
          beta[j][i] = b[j][i] / a[j][i];

        I_c0 += 1. / beta[j][i] * a[j][i] * beta[j][i] * beta[j][i];
        I_c0_err += a[j][i] * beta[j][i] * beta[j][i];
      }
    }
    I_c0 /= I_c0_err;
    for (i = 0; i < NDEC; i++)
    {
      for (j = 0; j < NRA; j++)
      {
        beta[j][i] *= I_c0;
      }
    }

    //printf("normlize=%f \n",I_c0);
  }
  ///////////////////////Dec. normal//////////////////
  else if (norm_f == 2)
  {
    double I_c[NDEC], I_c_err[NDEC];
    for (i = 0; i < NDEC; i++)
    {
      I_c[i] = 0.;
      I_c_err[i] = 0.;
      for (j = 0; j < NRA; j++)
      {
        if (a[j][i] == 0)
          beta[j][i] = 1.;
        else
          beta[j][i] = b[j][i] / a[j][i];
        I_c[i] += 1. / beta[j][i] * a[j][i] * beta[j][i] * beta[j][i];
        I_c_err[i] += a[j][i] * beta[j][i] * beta[j][i];
      }
      I_c[i] /= I_c_err[i];
    }
    for (i = 0; i < NDEC; i++)
    {
      for (j = 0; j < NRA; j++)
      {
        beta[j][i] *= I_c[i];
        b[j][i] = beta[j][i] * a[j][i];
      }
    }
  }

  else if (norm_f == 1)
  {
    double beta_c[NDEC];
    for (i = 0; i < NDEC; i++)
    {
      beta_c[i] = 0.;
      for (j = 0; j < NRA; j++)
      {
        if (a[j][i] == 0)
          beta[j][i] = 1.;
        else
          beta[j][i] = b[j][i] / a[j][i];
        beta_c[i] += 1. / beta[j][i];
      }
    }

    //beta_c=beta_c/(NRA*NDEC);
    for (i = 0; i < NDEC; i++)
    {
      beta_c[i] /= NRA;
      for (j = 0; j < NRA; j++)
      {
        beta[j][i] = beta[j][i] * beta_c[i];
        b[j][i] = beta[j][i] * a[j][i];
      }
      //printf("%f\n",beta_c[i]);
    }
    printf("yes\n");
  }

  return 0;
}

TFile *file_ang = new TFile("angres.root");

void wcda_ang_maxsigma(double *theta_maxsigma, double *maxsigma, double *keepratio_ga, double *keepratio_cr, double *num_ga, double *num_cr)
{
  double sigma[700];
  double n_ga, n_cr, n_cr_tot;
  double n_tmp = 0.;
  int n_maxsigma;

  n_cr_tot = ((TH1F *)file_ang->Get("hpol_cr_smth"))->Integral(1, 700);
  for (int i = 1; i <= 700; i++)
  {
    n_ga = ((TH1F *)file_ang->Get("hpol_smth"))->Integral(1, i);
    n_cr = n_cr_tot * (1. - cos(i * 0.01 * D2R)) / (1. - cos(7. * D2R));
    ;
    sigma[i - 1] = n_ga / sqrt(n_cr);
    if (sigma[i - 1] < n_tmp)
      continue;
    n_tmp = sigma[i - 1];
    n_maxsigma = i;
  }

  *theta_maxsigma = n_maxsigma * 0.01;
  *maxsigma = n_tmp;
  *num_ga = ((TH1F *)file_ang->Get("hpol_smth"))->Integral(1, n_maxsigma);
  *num_cr = ((TH1F *)file_ang->Get("hpol_cr_smth"))->Integral(1, n_maxsigma);
  *keepratio_ga = (((TH1F *)file_ang->Get("hpol_smth"))->Integral(1, n_maxsigma)) / (((TH1F *)file_ang->Get("hpol_smth"))->Integral(1, 700));
  *keepratio_cr = (((TH1F *)file_ang->Get("hpol_cr_smth"))->Integral(1, n_maxsigma)) / (((TH1F *)file_ang->Get("hpol_cr_smth"))->Integral(1, 700));
}
