#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "Astro.c"
#include "spectra_sim.cc"
#include "detector_sim.cc"
#include <TCanvas.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TMath.h>
#include <TF1.h>
#include <TH2.h>
#include <TString.h>

#include <time.h>
#include <iostream>
#include <fstream>

double coverage(double zenith);

int main(int argc,char *argv[])
{
  double emid,nTop_le,nTop_ue,copt;
  int smth_t = 1;

  if(argc<4) {
    printf("%s  effarea.root  emid  smooth_times  (use log10(emid))\n",argv[0]);
    exit(0);
  }
  emid = atof(argv[2]);
  if(argc==4) smth_t=atoi(argv[3]);

  emid = pow(10.0, emid - 12.0);

  if (emid<=exb[0]) {
    nTop_le = npmtb[0];
    copt = coptb[0];
  }
  else if (emid>=exb[nvvb-1]) {
    nTop_le = npmtb[0];
    copt =  coptb[nvvb-1];
  }
  else {
    int i=0;
    while (emid + 1e-4 > exb[i]) i++;
    nTop_le = npmtb[i-1];
    copt = coptb[i-1];
  }

  

printf("(nTop,copt) = ( %f , %f )\n",nTop_le,copt);
nTop_le = 128;
nTop_ue = 1000;
copt = 14.4;
double test1=0;
double test2=0;
double test3=0;
double test4=0;

  double dura_w[18000] = {0};
  int d_zen;
  double ZEN,AZI;
  double lst;
  double step_lst=360./(86400*1000);

  TString filename = Form("%s",argv[1]);

  for(lst=0.;lst<360.;lst+=step_lst)
  {
    equator_horizon_lst_ybj( lst, CRAB[0], CRAB[1], &ZEN, &AZI);
    d_zen = floor(100. * ZEN);
    dura_w[d_zen] += 1./(86400*1000);
  }

//  for(int i=0;i<18000;i++) Printf("%f\n",dura_w[i]);

  float energy,theta,thetaF,dAngle,zenistar,vcxPE[9],weight;
  int nTop;

  TFile *f_ga = new TFile("Rn_gamma.root");
  TTree *t1ga = (TTree*)f_ga->Get("t1");
  t1ga->SetBranchAddress("energy",&energy);
  t1ga->SetBranchAddress("theta",&theta);
  t1ga->SetBranchAddress("thetaF",&thetaF);
  t1ga->SetBranchAddress("dAngle",&dAngle);
  t1ga->SetBranchAddress("zenistar",&zenistar);
  t1ga->SetBranchAddress("vcxPE",&vcxPE);
  t1ga->SetBranchAddress("weight",&weight);
  t1ga->SetBranchAddress("nTop",&nTop);

  TFile *f_cr = new TFile("Rn_nuclei.root");
  TTree *t1cr = (TTree*)f_cr->Get("t1");
  t1cr->SetBranchAddress("energy",&energy);
  t1cr->SetBranchAddress("theta",&theta);
  t1cr->SetBranchAddress("thetaF",&thetaF);
  t1cr->SetBranchAddress("dAngle",&dAngle);
  t1cr->SetBranchAddress("zenistar",&zenistar);
  t1cr->SetBranchAddress("vcxPE",&vcxPE);
  t1cr->SetBranchAddress("weight",&weight);
  t1cr->SetBranchAddress("nTop",&nTop);

  double weight1, weight2, weight1_itgr, weight2_itgr;
  double zen_le=0.;
  double zen_ue=60.;
  double e_le=1.;
  double e_ue=5.;
  double zen_bin_width = 0.5;
  double e_bin_width = 0.04;
  int zen_bin_num = (zen_ue - zen_le)/zen_bin_width;
  int e_bin_num = (e_ue - e_le)/e_bin_width;
  double t_zen[zen_bin_num];
  double expo_zen[zen_bin_num];

  TFile *effarea = new TFile(filename,"recreate");
  TH2F *harea_p0 = new TH2F("harea_p0","effarea - energy x zenith",e_bin_num,e_le,e_ue,zen_bin_num,zen_le,zen_ue);
  TH2F *hratio_p0 = new TH2F("hratio_p0","ratio - energy x zenith",e_bin_num,e_le,e_ue,zen_bin_num,zen_le,zen_ue);
  TH2F *harea_p1 = new TH2F("harea_p1","effarea - energy x zenith",e_bin_num,e_le,e_ue,zen_bin_num,zen_le,zen_ue);
  TH2F *hratio_p1 = new TH2F("hratio_p1","ratio - energy x zenith",e_bin_num,e_le,e_ue,zen_bin_num,zen_le,zen_ue);
  TH1F *he_ga = new TH1F("he_ga","energy ga",e_bin_num,e_le,e_ue);
  TH1F *he_cr = new TH1F("he_cr","energy cr",e_bin_num,e_le,e_ue);
  TH1F *harea_p0_59 = new TH1F("harea_p0_59","effarea - zen",60,zen_le,zen_ue);
  TH1F *harea_p1_59 = new TH1F("harea_p1_59","effarea - zen",60,zen_le,zen_ue);

  double day_ga = 2.e6;
  double day_cr = 390.;
//  double copt=15.73;
  double d_ang = 15.;
  int ang_bin_num = int(100. * d_ang);

  double sumw;
  double A_eff;

  for(int i=0;i<zen_bin_num;i++)
  {
    t_zen[i] = 0.;
    for(d_zen = floor(i * zen_bin_width * 100.);d_zen < floor((i+1) * zen_bin_width * 100.);d_zen++)
    {
      t_zen[i] += dura_w[d_zen];
    }
    t_zen[i] *= day_ga * 86400.;
  }

  int n = t1ga->GetEntries();
  for(int i=0;i<n;i++)
  {
    t1ga->GetEntry(i);
test1=test1+weight;
    if(!(nTop>=nTop_le&&nTop<nTop_ue&&thetaF>0&&zenistar<=60)) continue;
test2=test2+weight;
    weight1 = (weight/t_zen[int(theta/zen_bin_width)])/CRAB_itgr_intens_HEGRA(int((9+log10(energy))/e_bin_width)*e_bin_width,int((9+log10(energy))/e_bin_width +1)*e_bin_width);
    if(isnan(weight1)!=0||isinf(weight1)!=0) weight1=0.;
    harea_p0->Fill(log10(energy),theta,weight1);
    if(log10(energy)>3.36&&log10(energy)<3.4)    harea_p0_59->Fill(theta,weight1);
    if(!(nTop/(vcxPE[6]+0.01)>copt&&dAngle<7.557)) continue;
test3=test3+weight;
    weight2 = (weight/t_zen[int(theta/zen_bin_width)])/CRAB_itgr_intens_HEGRA(int((9+log10(energy))/e_bin_width)*e_bin_width,int((9+log10(energy))/e_bin_width +1)*e_bin_width);
    if(isnan(weight2)!=0||isinf(weight2)!=0) weight2=0.;
    hratio_p0->Fill(log10(energy),theta,weight2);
    he_ga->Fill(log10(energy),weight);
  }
printf("test1,2,3 = %f,%f,%f\n",test1,test2,test3);
  for(int i = 0; i <= (e_bin_num+2) * (zen_bin_num+2); i++)  {
    weight1 = harea_p0->GetBinContent(i);
    weight2 = hratio_p0->GetBinContent(i);
    if(isnan(weight2/weight1)!=0||isinf(weight2/weight1)!=0)  {
      weight2 = 0.;
    }
    else {
    weight2 = weight2/weight1;
    }
    hratio_p0->SetBinContent(i,weight2);
  }

//  for(int i=0;i<zen_bin_num;i++)
//  {
//    t_zen[i] *= day_cr/day_ga;
//  }

  for(int i=0;i<zen_bin_num;i++)
//  for(int i=238;i<239;i++)
  {
    double d_zen_min;
    double d_zen_max;
    double min,max;
    double S_arc;
    double a,b,c,A;

    expo_zen[i] = 0.;

    min = (0>(i * zen_bin_width - 7.) * 100.) ? 0 : (i * zen_bin_width - 7.) * 100.;
    max = (60 * 100. < ((i+1) * zen_bin_width + 7.) * 100.) ? 60 * 100. : ((i+1) * zen_bin_width + 7.) * 100.;
//    d_zen_min = round(max(0 , (i * zen_bin_width - 7.) * 100.));
//    d_zen_max = round(min(60 * 100. , ((i+1) * zen_bin_width + 7.) * 100.));
    d_zen_min = min;
    d_zen_max = max;
    for(d_zen = round(d_zen_min);d_zen <= round(d_zen_max);d_zen++)
    {
      if(abs(d_zen/100. - (i+0.5) * zen_bin_width) > 7. + 0.5 * zen_bin_width - 1e-5)
      {
        S_arc = 0.;
      }
      else if(abs(d_zen/100. - (i+0.5) * zen_bin_width) < 7. - 0.5 * zen_bin_width + 1e-5)
      {
        if((d_zen/100.) + (i+0.5) * zen_bin_width < 7. + 1e-5)
        {
          A = PI;
        }
        else
        {
          a = 7. * D2R;
          b = (i+0.5) * zen_bin_width * D2R;
          c = (d_zen/100.) * D2R;
          A = acos((cos(a) - cos(b) * cos(c))/(sin(b) * sin(c)));
        }
        S_arc = 2. * A * sin((i+0.5) * zen_bin_width * D2R) * zen_bin_width * D2R;
      }
      else
      {
        if((d_zen/100.) < ((i+0.5) * zen_bin_width))
        {
          a = 7. * D2R;
          b = i * zen_bin_width * D2R;
          c = (d_zen/100.) * D2R;
          A = acos((cos(a) - cos(b) * cos(c))/(sin(b) * sin(c)));
          S_arc = A * sin(i * zen_bin_width * D2R) * zen_bin_width * D2R;
        }
        else
        {
          a = 7. * D2R;
          b = (i+1.) * zen_bin_width * D2R;
          c = (d_zen/100.) * D2R;
          A = acos((cos(a) - cos(b) * cos(c))/(sin(b) * sin(c)));
          S_arc = A * sin((i+1.) * zen_bin_width * D2R) * zen_bin_width * D2R;
        }
      }
      if(isnan(S_arc)!=0||isinf(S_arc)!=0) S_arc=0.;
      expo_zen[i] += dura_w[d_zen] * S_arc;
    }
    expo_zen[i] *= day_cr * 86400.;
//    printf("%d %f\n",i,expo_zen[i]);
  }

  n = t1cr->GetEntries();
  weight1_itgr = 0.;
  weight2_itgr = 0.;
  double tmp_w;
test1=0;
test2=0;
test3=0;
  for(int i=0;i<n;i++)
  {
    t1cr->GetEntry(i);
test1=test1+weight;
    if(!(nTop>=nTop_le&&nTop<nTop_ue&&thetaF>0&&zenistar<60.)) continue;
test2=test2+weight;
    weight1 = (weight/expo_zen[int(theta/zen_bin_width)])/(CR_itgr_intens_all(int((9+log10(energy))/e_bin_width)*e_bin_width,int((9+log10(energy))/e_bin_width +1)*e_bin_width));
//printf("111 %f %f %f\n",theta,expo_zen[int(theta/zen_bin_width)],weight1);
    if(isnan(weight1)!=0||isinf(weight1)!=0) weight1=0.;
    harea_p1->Fill(log10(energy),theta,weight1);
    if(log10(energy)>3.36&&log10(energy)<3.4)    harea_p1_59->Fill(theta,weight1);
    tmp_w = weight/expo_zen[int(theta/zen_bin_width)];
    if(isnan(tmp_w)!=0||isinf(tmp_w)!=0) tmp_w=0.;
    weight1_itgr += tmp_w;
    if(!(nTop/(vcxPE[6]+0.01)>copt)) continue;
test3=test3+weight;
    tmp_w = weight/expo_zen[int(theta/zen_bin_width)];
    if(isnan(tmp_w)!=0||isinf(tmp_w)!=0) tmp_w=0.;
    weight2_itgr += tmp_w;
    he_cr->Fill(log10(energy),weight);
//printf("222 %f %f\n",weight1_itgr,weight2_itgr);
  }
printf("test1,2,3 = %f,%f,%f\n",test1,test2,test3);
  if(isnan(weight2_itgr/weight1_itgr)!=0||isinf(weight2_itgr/weight1_itgr)!=0)  {
    weight2_itgr = 0.;
  }
  else {
  weight2_itgr = weight2_itgr/weight1_itgr;
  }

  for(int i = 0; i <= (e_bin_num+2) * (zen_bin_num+2); i++)  {
    hratio_p1->SetBinContent(i,weight2_itgr);
  }

//  TString cutstr = "";
//
//  cutstr = Form("weight*((nTop/vcxPE[6]>%f)&&dAngle<%f&&zenistar>=%f&&zenistar<%f&&thetaF>=0.&&energy>=%f&&energy<%f)",copt,d_ang,zen_le,zen_ue,e_le,e_ue);
//  TFile *f_ga = new TFile("/afs/ihep.ac.cn/users/c/changxc/eos/sci_code/simdata/Rn_gamma.root");
//  TTree *t1ga;
//  t1ga = (TTree*)f_ga->Get("t1");
//  TH1F *hdAngle_ga = new TH1F("hdAngle_ga","dAngle dist.",ang_bin_num,0.,d_ang);
//  t1ga->Draw("theta:energy>>harea_p0",cutstr);

//  sumw = harea_p0->GetSumOfWeights();

//  printf("%e\n",sumw);
  for(int tmp_t=0; tmp_t<smth_t; tmp_t++)
  {
    harea_p0->Smooth(1,"k5b");
    harea_p1->Smooth(1,"k5b");
    hratio_p0->Smooth(1,"k5b");
    hratio_p1->Smooth(1,"k5b");
  }
  effarea->Write();
  return 0;
}


double coverage(double zenith)
{
//8.08 is the minimum of zenith of CRAB track, 67.93 is the maximum; 7 is the radii of observing window.
  double round_over_arch;
//round_over_arch is a factor which should be multiplied to the effective area at a fixed theta.
  if(zenith>=(8.08+7.) && zenith<=(67.93-7.))
  {
    round_over_arch = 1.;
  }
  else if(zenith<=(8.08-7.) || zenith>=(67.93+7.))
  {
    round_over_arch = 0.;
  }
  else if(zenith>(8.08-7.) && zenith<(8.08+7.))
  {
    round_over_arch = PI*7.*7./(7.*7.*acos((8.08-zenith)/7.)-(8.08-zenith)*sqrt(7.*7.-(8.08-zenith)*(8.08-zenith)));
  }
  else if(zenith>(67.93-7.) && zenith<(67.93+7.))
  {
    round_over_arch = PI*7.*7./(7.*7.*acos((zenith-67.93)/7.)-(zenith-67.93)*sqrt(7.*7.-(zenith-67.93)*(zenith-67.93)));
  }
  else
  {
    round_over_arch = 0.;
  }
  return round_over_arch;
}

