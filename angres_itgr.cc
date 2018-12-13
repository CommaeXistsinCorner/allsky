#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "Astro.c"
#include "detector_sim.cc"
#include <TCanvas.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TMath.h>
#include <TF1.h>
#include <TString.h>
#include <TF2.h>
#include <TH2F.h>

//double bigaussian(double *x, double *par);
//double gaussian(double *x, double *par);
double gaussian_polar_coordinates(double *x, double *par);

void error_proj(double theta, double phi, double thetaF, double phiF, double *kx, double *ky)
{
  double distance;
  double direction;

  distance = distance_horizontal(theta, phi, thetaF, phiF);
  direction = direction_horizontal(theta, phi, thetaF, phiF);

  *kx = sin(distance * deg_rad) * sin(direction * deg_rad) * rad_deg;
  *ky = sin(distance * deg_rad) * cos(direction * deg_rad) * rad_deg;
}

int main(int argc, char *argv[])
{
  float theta, phi, thetaF, phiF;
  double kx, ky;
  //  float e_le,e_ue,zen_le,zen_ue,dAngle_ue;
  float zenistar, cAngle, dAngle, vcxPE[9], weight;
  int nTop;
  double nTop_le, nTop_ue;
  double emin, copt;
  //  dAngle_ue = 5.;
  int smth_t = 1;
  int nopic = 0;
  double day_ga = 2.e6;
  double day_cr = 390.;

  if (argc < 4)
  {
    printf("%s  angres.root  emin  [smooth_times]  (use log10(emin))\n", argv[0]);
    exit(0);
  }

  emin = atof(argv[2]);
  if (argc == 4)
    smth_t = atoi(argv[3]);

  emin = pow(10.0, emin - 12.0);

  if (emin <= exb[0])
  {
    nTop_le = npmtb[0];
    copt = coptb[0];
  }
  else if (emin >= exb[nvvb - 1])
  {
    nTop_le = npmtb[0];
    copt = coptb[nvvb - 1];
  }
  else
  {
    int i = 0;
    while (emin + 1e-4 > exb[i])
      i++;
    nTop_le = npmtb[i - 1];
    copt = coptb[i - 1];
  }
  printf("(nTop,copt) = ( %f , %f )\n", nTop_le, copt);
  nTop_le = 128;
  nTop_ue = 10000;
  copt = 14.4;

  TString filename = Form("%s", argv[1]);
  TFile *f_ga = new TFile("Rn_gamma.root");
  TTree *t1ga = (TTree *)f_ga->Get("t1");
  t1ga->SetBranchAddress("nTop", &nTop);
  t1ga->SetBranchAddress("weight", &weight);
  t1ga->SetBranchAddress("thetaF", &thetaF);
  //  t1->SetBranchAddress("energy",&energy);
  t1ga->SetBranchAddress("theta", &theta);
  t1ga->SetBranchAddress("phi", &phi);
  t1ga->SetBranchAddress("phiF", &phiF);
  t1ga->SetBranchAddress("zenistar", &zenistar);
  t1ga->SetBranchAddress("dAngle", &dAngle);
  t1ga->SetBranchAddress("vcxPE", &vcxPE);

  TFile *f_cr = new TFile("Rn_nuclei.root");
  TTree *t1cr = (TTree *)f_cr->Get("t1");
  t1cr->SetBranchAddress("nTop", &nTop);
  t1cr->SetBranchAddress("weight", &weight);
  t1cr->SetBranchAddress("thetaF", &thetaF);
  t1cr->SetBranchAddress("zenistar", &zenistar);
  t1cr->SetBranchAddress("cAngle", &cAngle);
  t1cr->SetBranchAddress("vcxPE", &vcxPE);

  int n_ga = t1ga->GetEntries();
  int n_cr = t1cr->GetEntries();

  TFile *fa = new TFile(filename, "recreate");
  TH2F *h = new TH2F("h", "kx_ky", 1001, -5.005, 5.005, 1001, -5.005, 5.005);
  TH1F *hx = new TH1F("hx", "kx", 201, -5.025, 5.025);
  TH1F *hy = new TH1F("hy", "ky", 201, -5.025, 5.025);
  TH1F *hpol = new TH1F("hpol", "theta", 700, 0, 7);
  TH1F *hpol_smth = new TH1F("hpol_smth", "theta", 700, 0, 7);
  TH1F *hpol_rho = new TH1F("hpol_rho", "theta", 700, 0, 7);
  TH1F *hpol_cr_smth = new TH1F("hpol_cr_smth", "theta", 700, 0, 7);
  TH1F *hsigma = new TH1F("hsigma", "theta", 700, 0, 7);

  for (int i = 0; i < n_ga; i++)
  {
    t1ga->GetEntry(i);
    if (nTop < nTop_le || nTop >= nTop_ue || thetaF < 0 || zenistar > 60. || nTop / (vcxPE[6] + 0.01) < copt || dAngle > 7.557)
      continue;
    error_proj(theta, phi, thetaF, phiF, &kx, &ky);
    h->Fill(kx, ky, weight);
    hx->Fill(kx, weight);
    hy->Fill(ky, weight);
    hpol->Fill(dAngle, weight);
    hpol_smth->Fill(dAngle, weight / day_ga);
  }

  double cal_factor = hpol->Integral();
  hpol->Scale(1. / cal_factor);

  for (int i = 0; i < n_cr; i++)
  {
    t1cr->GetEntry(i);
    if (nTop < nTop_le || nTop >= nTop_ue || thetaF < 0 || zenistar > 60. || nTop / (vcxPE[6] + 0.01) < copt)
      continue;
    hpol_cr_smth->Fill(cAngle, weight / day_cr);
  }

  hpol->Smooth(smth_t);
  hpol_smth->Smooth(smth_t);
  hpol_cr_smth->Smooth(smth_t);
  h->Smooth(smth_t);
  hx->Smooth(smth_t);
  hy->Smooth(smth_t);

  for (int i = 1; i <= 700; i++)
  {
    double a = hpol->GetBinContent(i);
    a = a / (2. * PI * (cos((i - 1) * 0.01 * D2R) - cos(i * 0.01 * D2R)) * R2D * R2D);
    hpol_rho->SetBinContent(i, a);
  }
  hpol_rho->Smooth(smth_t);

  TCanvas *c1 = new TCanvas("c1");
  hpol->Draw();
  hpol_rho->Draw();

  //  hpol->Fit("pol9");

  //  TF1 *fgpol = new TF1("fgpol",bigaussian,-5,5,6);

  fa->Write();
  return 0;
  //  c1->SaveAs(filename);
}

double gaussian_polar_coordinates(double *x, double *par)
{
  return par[0] * exp(-x[0] * x[0] / (2. * par[1] * par[1])) * x[0];
}

//double bigaussian(double *x, double *par)
//{
//  return gaussian(x,par) + gaussian(x,&par[3]);
//}
//
//double gaussian(double *x, double *par)
//{
//  return (par[0]/par[2])*(exp(-((x[0]-par[1])*(x[0]-par[1]))/(2.*par[2]*par[2])));
//}
