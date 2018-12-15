#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TH2.h"
#include "TH1.h"
#include "ana_wid_sim.h"
#include "Astro.c"
#include "detector_sim.cc"
#include "spectra_sim.cc"
//#include "angres_func.cc"
//#include "point_source_wcda.cc"
//#include "detector_wcda_zenith_cal.cpp"
//#define Tibet_Lngi      90.522                   /* Longitude of Tibet (deg.) */
//#define Tibet_Lati      30.102                   /* Latitude of Tibet  (deg.) */
//#define mjd_st 56789.0    /* a arbitrary virtual start time to calculate the position of point sources   */
//static  float   n[NLST][NZEAZ];
//static  float   nb_tmp[NLST][NZE];
//memset(n,0.0,sizeof(n));
//memset(nb_tmp,0.0,sizeof(nb_tmp));

///////////// function definations!!!//////////
int rand_possion(double lambda);

double rand_gauss(double mu, double sigma);

double bigaussian_mean(double theta, double sigma);

double bigaussian(double theta, double sigma);

void e_median(double *e_ga, double *e_cr);

double wcda_eff_area_ga(double eg, double zen);

double wcda_eff_area_cr(double ep, double zen);

double wcda_filter_ga(double eg, double zen);

double wcda_filter_cr(double ep, double zen);

double wcda_ang(double ep);

double wcda_theta_dist(double theta);

double wcda_theta_dist1(double theta);

void wcda_ang_maxsigma(double day_num, double *theta_maxsigma, double *maxsigma, double *keepratio_ga, double *keepratio_cr, double *num_ga, double *num_cr);
///////////// function definations end!!!/////

int main(int argc, char *argv[])
{
  //printf("%d,%d",NLST,NZEAZ);
  int k0;
  double evnum_p[NZE];
  double evnum_ga[NZE];
  double area_nze[NZE];

  //const int NZEAZ = nzeaz();

  //static  float   n[NLST][NZEAZ];
  //static  float   nb_tmp[NLST][NZE];
  //static double grid_ra[NZEAZ],grid_dec[NZEAZ];
  //3 lines above are tradition path to create arrays.

  //memset(n,0.0,sizeof(n));
  //memset(nb_tmp,0.0,sizeof(nb_tmp));

  //  float  (*n)[NZEAZ] = new float[NLST][NZEAZ];
  //  float  (*nb_tmp)[NZEAZ] = new float[NLST][NZEAZ];
  //  memset(nb_tmp,0.0,sizeof(nb_tmp));

  FILE *fp_out, *fp_log;
  FILE *fp_nb, *fp_grid_ra_dec;
  int isr, k;
  double cr_max = 0.;
  double eTeV, nTop_le, copt;
  double emin_ga = 10., emin_cr = 10., emax_ga = 14., emax_cr = 14., emid_ga, emid_cr;
  double days;
  int i, j;
  double lst, RA, DEC, ZEN, AZI;
  float tmp_stream1, tmp_stream2;
  double zen, azim, ra, dec;
  int ismon = 1, isroot = 1;
  //  char         *tmp_n,*tmp_nb,*tmp_grid_ra_dec;
  int N_int, n_bin;
  double sum_int;
  double bin_h, bin_area;
  //  double       grid_sr;
  char filename[200];

  if (argc < 5)
  {
    printf("%s  output(.dat)  emid_ga  emax_ga  days  [1 for .root]  [1 for monitor mode]  (use log10(emid_ga/emax_ga))  set emax_ga=0 for default: 14.\n", argv[0]);
    exit(0);
  }
  emid_ga = atof(argv[2]);
  if (atof(argv[3]) > 1e-4)
    emax_ga = atof(argv[3]);
  if (emid_ga > emax_ga)
  {
    printf("emax_ga should be larger than emid_ga!\n");
    exit(0);
  }
  days = atof(argv[4]);
  if (argc >= 6)
    isroot = atoi(argv[5]);
  if (argc >= 7)
    ismon = atoi(argv[6]);

  strcpy(filename, argv[1]);
  if ((fp_out = fopen(strcat(argv[1], ".dat"), "wb+")) == NULL)
  {
    printf("cannot create output file\n");
    exit(0);
  }
  strcpy(argv[1], filename);
  if ((fp_log = fopen(strcat(argv[1], ".log"), "w+")) == NULL)
  {
    printf("cannot create log file\n");
    exit(0);
  }
  strcpy(argv[1], filename);
  if ((fp_nb = fopen(strcat(argv[1], "_tmp_nb"), "wb+")) == NULL)
  {
    printf("cannot create temp files\n");
    exit(0);
  }
  strcpy(argv[1], filename);
  fp_grid_ra_dec = fopen(strcat(argv[1], "_tmp_grid_ra_dec"), "wb+");

  fprintf(fp_log, "\n\n//////////simdata.exe//////////\n");

  eTeV = pow(10.0, emid_ga - 12.0);

  if (eTeV <= exb[0])
  {
    nTop_le = npmtb[0];
    copt = coptb[0];
  }
  else if (eTeV >= exb[nvvb - 1])
  {
    nTop_le = npmtb[0];
    copt = coptb[nvvb - 1];
  }
  else
  {
    int i = 0;
    while (eTeV + 1e-4 > exb[i])
      i++;
    nTop_le = npmtb[i - 1];
    copt = coptb[i - 1];
  }
  printf("(nTop,copt) = ( %f , %f )\n", nTop_le, copt);
  fprintf(fp_log, "(nTop,copt) = ( %f , %f )\n", nTop_le, copt);
  nTop_le = 128;
  copt = 14.4;

  //calculate the median energy of gamma and CRs.
  e_median(&emid_ga, &emid_cr);

  printf("(Emid_ga,Emid_cr) = (%g,%g)\n", emid_ga, emid_cr);
  fprintf(fp_log, "(Emid_ga,Emid_cr) = (%g,%g)\n", emid_ga, emid_cr);

  //to cancel median energy

  printf("(Emin,Emax) = (%g,%g)\n", emin_ga, emax_ga);
  fprintf(fp_log, "(Emin,Emax) = (%g,%g)\n", emin_ga, emax_ga);

  for (int tmp_i = 0, tmp_k = 0; tmp_i < NZE; tmp_i++)
  {
    NAZ0[tmp_i] = ceil(360.0 * sin((0.5 + tmp_i) * hori_sys_width * PI / 180.0) / hori_sys_width);
    NAZ[tmp_i] = tmp_k;
    tmp_k += NAZ0[tmp_i];
  }

  for (int tmp_i = 0; tmp_i < NZE; tmp_i++)
  {
    double thetamin = tmp_i * hori_sys_width;
    double thetamax = (tmp_i + 1) * hori_sys_width;
    //    area_nze[tmp_i] = PI * ( thetamax * thetamax - thetamin * thetamin )  / NAZ0[tmp_i];
    area_nze[tmp_i] = 2 * PI * (cos(thetamin * D2R) - cos(thetamax * D2R)) * R2D * R2D / NAZ0[tmp_i];
    //printf("%f %f %f %d\n",thetamin,thetamax,area_nze[tmp_i],NAZ0[tmp_i]);
  }

  if (ismon == 1)
    printf("parameters fixed\n");

  srand((unsigned)time(NULL));

  // 10^0.01 is the bin width

  //double aaa=wcda_eff_area_gamma(12.8,1.0);
  //printf("%f\n",aaa);
  N_int = (int)100 * (emax_ga - emin_ga);
  //  double  sum_bin[NZE][N_int];

  for (i = 0; i < NZE; i++)
  {
    sum_int = 0.0;
    for (n_bin = 0; n_bin < N_int; n_bin++)
    {
      bin_h = CR_intens_all(emin_cr + 0.005 + 0.01 * n_bin);
      bin_area = bin_h * 0.023293 * pow(10.0, emin_cr + 0.01 * n_bin);
      // 0.023293=10^0.01-1
      //      sum_bin[i][n_bin] = bin_area * wcda_filter_p(emid + 0.005 + 0.01*n_bin);
      sum_int += bin_area * wcda_eff_area_cr(emin_cr + 0.005 + 0.01 * n_bin, (i + 0.5) * equa_sys_width) * wcda_filter_cr(emin_cr + 0.005 + 0.01 * n_bin, (i + 0.5) * equa_sys_width);
      //printf("%d %g %g \n",n_bin,bin_h,wcda_eff_area_cr(emid + 0.005 + 0.01*n_bin, (i+0.5)*equa_sys_width));
      //printf("%d %f \n",i,bin_area* wcda_eff_area_cr(emid + 0.005 + 0.01*n_bin, (i+0.5)*equa_sys_width));
      //if(i==1.0&&n_bin==0) {printf("!!! CR_intens_p1 bin_h bin_area sum_int %e %f %f %f\n", CR_intens_p1(emid + 0.005 + 0.01*n_bin),  bin_h, bin_area, sum_int);}
    }
    ////  grid_sr is the same to area_nze[];
    //  grid_sr = 2.0 * PI * ( cos(i*hori_sys_width*D2R) - cos((i+1)*hori_sys_width*D2R) ) / NAZ0[i];

    evnum_p[i] = sum_int * area_nze[i] * D2R * D2R * dt_time * days * 86400. * 0.9354;
    //0.9354 is a parameter has no physical meaning but to calibrate to reproduce a correct bg events number according to Yao's data!!!!!!!!!

    //0.0000760631 is sr of the little grid
    //printf("evnum_p[%d]= %f \n",i,evnum_p[i]);
  }
  if (ismon == 1)
    printf("background event numbers in grid has been calculated \n");
  for (isr = 0; isr < NLST; isr++)
  {
    if (ismon == 1)
      printf("background set!   %3d%%\r", int(isr * 100 / NLST));
    for (i = 0; i < NZE; i++)
    {
      float buffer1[NAZ0[i]] = {0};
      for (j = 0; j < NAZ0[i]; j++)
      {
        if (evnum_p[i] > 50)
        {
          buffer1[j] = rand_gauss((double)evnum_p[i], (double)sqrt(evnum_p[i]));
          //            tmp_stream1 = rand_gauss((double)evnum_p[i],(double)sqrt(evnum_p[i]));
          //            if(tmp_stream1<0)  tmp_stream1 = 0.;
        }
        else
        {
          buffer1[j] = rand_possion((double)evnum_p[i]);
          //            tmp_stream1 = rand_possion((double)evnum_p[i]);
        }
        if (isnan(buffer1[j]) != 0 || isinf(buffer1[j]) != 0 || buffer1[j] < 0)
          buffer1[j] = 0.;
        //        if(isnan(tmp_stream1)!=0||isinf(tmp_stream1)!=0)  tmp_stream1 = 0.;
        //tmp_stream1=0.;
        //        fwrite(&tmp_stream1,sizeof(tmp_stream1),1,fp_out);
        cr_max += buffer1[j];
        //        cr_max+=tmp_stream1;
      }
      fwrite(buffer1, sizeof(float), NAZ0[i], fp_out);
    }
    //    printf("%d\n",isr);
    //    printf("background set!  isr = %5d in %d \r",isr,NLST);
  }
  printf("\nbg_event = %e\n", cr_max);
  fprintf(fp_log, "bg_event = %e\n", cr_max);

  if (isroot == 1)
  {
    strcpy(argv[1], filename);
    TFile hfile(strcat(argv[1], ".root"), "RECREATE");
    if (ismon == 1)
      printf("ROOT file part 1 is preparing\n");
    TH1D *hzen = new TH1D("hzen_in", "zenith dis.", 90, 0., 90.);
    TH1D *hazim = new TH1D("hazim_in", "azimuth dis.", 79, 0., 360.);
    TH1D *hra = new TH1D("hra_in", "RA dis.", 180, 0., 360.);
    TH1D *hdec = new TH1D("hdec_in", "DEC dis.", 136, -46., 90.);
    TH1D *hsig = new TH1D("hsig_in", "Signi", 101, -20, 20);
    TH2D *hradec = new TH2D("hradec_in", "RADEC dis.", 180, 0., 360., 135, -45., 90.);
    lst = 0.5 * equa_sys_width;
    for (i = 0; i < NZE; i++)
    {
      zen = (i + 0.5) * hori_sys_width;
      for (j = 0; j < NAZ0[i]; j++)
      {
        azim = (j + 0.5) * 360. / NAZ0[i];
        htoe2(zen, azim, lst, &ra, &dec);
        fwrite(&ra, sizeof(ra), 1, fp_grid_ra_dec);
        fwrite(&dec, sizeof(dec), 1, fp_grid_ra_dec);
      }
    }

    for (rewind(fp_out), isr = 0; isr < NLST; isr++)
    {
      for (rewind(fp_grid_ra_dec), i = 0; i < NZE; i++)
      {
        //        tmp_stream2=0.0;
        float buffer1[NAZ0[i]] = {0};
        fread(buffer1, sizeof(float), NAZ0[i], fp_out);
        for (j = 0; j < NAZ0[i]; j++)
        {
          //          fread(&tmp_stream1,sizeof(tmp_stream1),1,fp_out);
          tmp_stream2 += buffer1[j] / NAZ0[i];
          hzen->Fill((i + 0.5) * hori_sys_width, buffer1[j]);
          //          hazim->Fill((j-0.5+rand()/RAND_MAX)*360./NAZ0[i],tmp_stream1);
          hazim->Fill((j + 0.5) * 360. / NAZ0[i], buffer1[j]);
          fread(&ra, sizeof(ra), 1, fp_grid_ra_dec);
          fread(&dec, sizeof(dec), 1, fp_grid_ra_dec);
          ra = (ra + isr * equa_sys_width) < 360.0 ? (ra + isr * equa_sys_width) : (ra + isr * equa_sys_width - 360.0);
          hra->Fill(ra, buffer1[j]);
          hdec->Fill(dec, buffer1[j]);
          hradec->Fill(ra, dec, buffer1[j]);
        }
        fwrite(&tmp_stream2, sizeof(tmp_stream2), 1, fp_nb);
        for (j = 0; j < NAZ0[i]; j++)
        {
          //          fread(&tmp_stream1,sizeof(tmp_stream1),1,fp_out);
          hsig->Fill((buffer1[j] - tmp_stream2) / sqrt(tmp_stream2));
        }
      }
      if (ismon == 1)
        printf("writing data!   %3d%%\r", int(isr * 100 / NLST));
      //      printf("writing data! \t isr = %d in %d \r",isr,NLST);
    }

    hfile.Write();
    hfile.Close();
  }

  if (PS_INPUT == 1)
  {
    double ephi;
    double ang_reso, ang_opt, maxsigma, keepratio_ga, keepratio_cr, num_ga, num_cr;
    double rad_over_ang_reso = 7.0;
    double ga_max[300] = {0};
    double ga_max_reso0[300] = {0};
    double bg_max[300] = {0};
    double bg_max_tmp;
    double zero_num[300] = {0};
    double S_ang_opt;
    double S_ang_opt_raw;
    int i0, j0;
    int i_err;
    int imin, imax, jmin, jmax;
    double theta;
    float tmp_num_ga;
    double phi;
    int i_PS, PS_num;

    PS_num = sizeof(PS_TEVCAT) / (sizeof(double) * 7);
    double evnum_ga[PS_num][NZE];
    double ga_num[PS_num][N_int];

    //    for(n_bin=0;n_bin<N_int;n_bin++) ga_num[n_bin]=0;
    ang_opt = wcda_ang(emid_ga);
    wcda_ang_maxsigma(days, &ang_opt, &maxsigma, &keepratio_ga, &keepratio_cr, &num_ga, &num_cr);
    printf("(theta_maxsigma, maxsigma, keepratio_ga, keepratio_cr, num_ga, num_cr) = (%g, %g, %g, %g, %g, %g)\n", ang_opt, maxsigma, keepratio_ga, keepratio_cr, num_ga, num_cr);
    fprintf(fp_log, "(theta_maxsigma, maxsigma, keepratio_ga, keepratio_cr, num_ga, num_cr) = (%g, %g, %g, %g, %g, %g)\n", ang_opt, maxsigma, keepratio_ga, keepratio_cr, num_ga, num_cr);
    //!!!!!!!!!!!!!!!!!!!!!modified!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ang_opt = 0.557;
    //!!!!!!!!!!!!!!!!!!!!!modified!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ang_reso = ang_opt / 1.5852; // almost 70% within wcda_ang(emid)
    S_ang_opt = 2. * PI * (1. - cos(ang_opt * D2R)) * R2D * R2D;

    i_err = int(rad_over_ang_reso * ang_reso / hori_sys_width);
    printf("%d sources!\n", PS_num);
    fprintf(fp_log, "%d sources!\n", PS_num);

    for (i_PS = 0; i_PS < PS_num; i_PS++)
    {
      for (i = 0; i < NZE; i++)
      {
        sum_int = 0.0;
        ephi = 0.0;
        phi = 0.;
        for (n_bin = 0; n_bin < N_int; n_bin++)
        {
          bin_h = SOURCE_intens(i_PS, emin_ga + 0.005 + 0.01 * n_bin);
          //bin_h = TEST_SOURCE_intens(emin_ga + 0.005 + 0.01*n_bin);
          //bin_h = CRAB_intens_HEGRA(emin_ga + 0.005 + 0.01*n_bin);
          bin_area = bin_h * 0.023293 * pow(10.0, emin_ga + 0.01 * n_bin);
          // 0.023293=10^0.01
          phi += bin_area;
          ephi += pow(10.0, emin_ga) * bin_area;
          //        sum_bin[i][n_bin] = bin_area * wcda_eff_area_gamma(emid + 0.005 + 0.01*n_bin, (i+0.5)*equa_sys_width) * wcda_filter_ga(emid + 0.005 + 0.01*n_bin);
          ga_num[i_PS][n_bin] += bin_area * wcda_eff_area_ga(emin_ga + 0.005 + 0.01 * n_bin, (i + 0.5) * equa_sys_width) * dt_time * days * 86400.0;
          //          sum_int += bin_area * wcda_eff_area_ga(emin_ga + 0.005 + 0.01*n_bin, (i+0.5)*equa_sys_width);
          sum_int += bin_area * wcda_eff_area_ga(emin_ga + 0.005 + 0.01 * n_bin, (i + 0.5) * equa_sys_width) * wcda_filter_ga(emin_ga + 0.005 + 0.01 * n_bin, (i + 0.5) * equa_sys_width);
          //        sum_bin[i][n_bin] *= dt_time * days * 86400.0;
        }
        evnum_ga[i_PS][i] = sum_int * dt_time * days * 86400.0;
        //printf("evnum_ga[%d]= %f \n",i,evnum_ga[i]);
      }
    }
    //for(n_bin=0;n_bin<100;n_bin++)
    //{
    //  bin_area = wcda_eff_area_gamma(emid + 0.02 + 0.04*n_bin,8.);
    //  printf("%f %e\n",emid + 0.02 + 0.04*n_bin,bin_area);
    //}

    if (ismon == 1)
      printf("\npoint source event numbers in grid has been calculated.\n");
    //    printf("PHI(>Emin) = %e\n", phi);
    //    printf("E*dPHI(>Emin) = %f  (eV*cm^-2*s^-1)\n",ephi);
    printf("ang_reso (Emin) = %f (degree)\n", ang_reso);
    fprintf(fp_log, "ang_reso (Emin) = %f (degree)\n", ang_reso);
    //    fprintf(fp_log,"E*dPHI(>Emin) = %f  (eV*cm^-2*s^-1)\n",ephi);
    //    fprintf(fp_log,"angular radius (Emin) = %f (degree)\n",ang_reso);
    //double aaaa=0.;
    //equator_horizon_lst( 144., CRAB[0], CRAB[1], &ZEN, &AZI);
    //printf("======\n %f  %f\n",ZEN, AZI);

    //    for(i_PS=63;i_PS<64;i_PS++)
    //    for(i_PS=0;i_PS<PS_num;i_PS++)  {
    for (i_PS = 0; i_PS < 1; i_PS++)
    {
      for (isr = 0; isr < NLST; isr++)
      {
        lst = (isr + 0.5) * equa_sys_width;
        equator_horizon_lst(lst, PS_TEVCAT[i_PS][0], PS_TEVCAT[i_PS][1], &ZEN, &AZI);
        //equator_horizon_lst( lst, PS_TEST[0], PS_TEST[1], &ZEN, &AZI);
        //equator_horizon_lst( lst, CRAB[0], CRAB[1], &ZEN, &AZI);
        if (ZEN > zex)
          continue;
        if (AZI < 0)
          AZI += 360.;

        i0 = int(ZEN / hori_sys_width);
        if (evnum_ga[i_PS][i0] > 50)
        {
          tmp_num_ga = rand_gauss((double)evnum_ga[i_PS][i0], (double)sqrt(evnum_ga[i_PS][i0]));
        }
        else
        {
          tmp_num_ga = rand_possion((double)evnum_ga[i_PS][i0]);
        }
        //        printf("%d %f %f\n",i_PS,evnum_ga[i_PS][i0],tmp_num_ga);
        //      max+=tmp_num_ga;
        if (isnan(tmp_num_ga) != 0 || isinf(tmp_num_ga) != 0 || tmp_num_ga < 0)
          tmp_num_ga = 0.;
        ga_max_reso0[i_PS] += evnum_ga[i_PS][i0];

        bg_max_tmp = 0.;
        S_ang_opt_raw = 0.;

        i_err = int(rad_over_ang_reso * ang_reso / hori_sys_width);
        //
        imin = (i0 - i_err) > 0 ? i0 - i_err : 0;
        imax = (i0 + i_err) < NZE ? i0 + i_err : NZE;
        for (i = imin; i < imax; i++)
        {
          j0 = int(AZI * NAZ0[i] / 360.);
          jmin = (j0 - i_err) > 0 ? j0 - i_err : 0;
          jmax = (j0 + i_err) < NAZ0[i] ? j0 + i_err : NAZ0[i];
          float buffer1[jmax - jmin + 1] = {0};
          float buffer2[jmax - jmin + 1] = {0};
          fseek(fp_out, (isr * NZEAZ + NAZ[i] + jmin) * sizeof(float), SEEK_SET);
          fread(buffer1, sizeof(float), jmax - jmin + 1, fp_out);
          for (j = jmin; j < jmax; j++)
          {
            k = NAZ[i] + j;

            //            fseek(fp_out,(isr*NZEAZ+k)*sizeof(tmp_stream1),SEEK_SET);
            //            fread(&tmp_stream1,sizeof(tmp_stream1),1,fp_out);
            //printf("tmp_stream1 %f\n",tmp_stream1);
            theta = distance_horizontal((i0 + 0.5) * hori_sys_width, 360.0 * (j0 + 0.5) / NAZ0[i], (i + 0.5) * hori_sys_width, 360.0 * (j + 0.5) / NAZ0[i]);
            if (PS_TEVCAT[i_PS][6] < 0.01)
            {
              buffer2[j - jmin] = buffer1[j - jmin] + tmp_num_ga * wcda_theta_dist1(theta) * area_nze[i];
            }
            else
            {
              buffer2[j - jmin] = buffer1[j - jmin] + tmp_num_ga * bigaussian_mean(theta, sqrt(PS_TEVCAT[i_PS][6] * PS_TEVCAT[i_PS][6] + ang_reso * ang_reso)) * area_nze[i];
            }
            //printf("tmp_stream2 %f\n",tmp_stream2);
            if (buffer2[j - jmin] < 1e-12)
              zero_num[i_PS]++;
            //            fseek(fp_out,(isr*NZEAZ+k)*sizeof(tmp_stream1),SEEK_SET);
            //            fwrite(&tmp_stream2,sizeof(tmp_stream1),1,fp_out);
            if (theta < ang_opt)
            {
              bg_max_tmp += buffer1[j - jmin];
              ga_max[i_PS] += buffer2[j - jmin] - buffer1[j - jmin];
              S_ang_opt_raw += area_nze[i];
            }
          }
          fseek(fp_out, -(jmax - jmin + 1) * sizeof(float), SEEK_CUR);
          fwrite(buffer2, sizeof(float), jmax - jmin + 1, fp_out);
        }
        //  bg_max[i_PS] += bg_max_tmp * ( S_ang_opt  / S_ang_opt_raw );
        bg_max[i_PS] += bg_max_tmp;
        printf("max %f\n", bg_max_tmp);
        if (ismon == 1)
          printf("PS %4d of %4d set!   %3d%%\r", i_PS, PS_num, int(isr * 100 / NLST));
        //      printf("point source set! \t isr = %d in %d \r",isr,NLST);
      }
      //exit(0);
      printf("\n( i_PS, reso0_source_event, source_event_within_ang_opt, bg_event_within_ang_opt, raw_sigma, grids have zero event numbers ) = ( %4d, %f, %f, %f, %f, %f )\n", i_PS, ga_max_reso0[i_PS], ga_max[i_PS], bg_max[i_PS], ga_max[i_PS] / sqrt(bg_max[i_PS]), zero_num[i_PS]);
      fprintf(fp_log, "( i_PS, reso0_source_event, source_event_within_ang_opt, bg_event_within_ang_opt, raw_sigma, grids have zero event numbers ) = ( %4d, %f, %f, %f, %f, %f )\n", i_PS, ga_max_reso0[i_PS], ga_max[i_PS], bg_max[i_PS], ga_max[i_PS] / sqrt(bg_max[i_PS]), zero_num[i_PS]);
    }

    if (isroot == 1)
    {
      strcpy(argv[1], filename);
      TFile hfile(strcat(argv[1], ".root"), "UPDATE");
      if (ismon == 1)
        printf("ROOT file part2 is preparing\n");
      TH1D *hzen_ps = new TH1D("hzen_in_ps", "zenith dis.", 45, 0., 90.);
      TH1D *hazim_ps = new TH1D("hazim_in_ps", "azimuth dis.", 180, 0., 360.);
      TH1D *hra_ps = new TH1D("hra_in_ps", "RA dis.", 180, 0., 360.);
      TH1D *hdec_ps = new TH1D("hdec_in_ps", "DEC dis.", 136, -46., 90.);
      TH1D *hsig_ps = new TH1D("hsig_in_ps", "Signi", 101, -20, 20);
      TH1D *hnum_ps = new TH1D("hnum_ps", "energy.", 30, 12, 15);
      for (i = 0; i < N_int; i++)
      {
        hnum_ps->Fill(12. + (i / 100.), ga_num[i_PS][i]);
      }
      k = 0;
      lst = 0.5 * equa_sys_width;
      for (rewind(fp_grid_ra_dec), i = 0; i < NZE; i++)
      {
        zen = (i + 0.5) * hori_sys_width;
        for (j = 0; j < NAZ0[i]; j++)
        {
          azim = (j + 0.5) * 360. / NAZ0[i];
          htoe2(zen, azim, lst, &ra, &dec);
          fwrite(&ra, sizeof(ra), 1, fp_grid_ra_dec);
          fwrite(&dec, sizeof(dec), 1, fp_grid_ra_dec);
        }
      }

      for (rewind(fp_out), rewind(fp_nb), isr = 0; isr < NLST; isr++)
      {
        for (rewind(fp_grid_ra_dec), i = 0; i < NZE; i++)
        {
          fread(&tmp_stream2, sizeof(tmp_stream2), 1, fp_nb);
          float buffer1[NAZ0[i]] = {0};
          fread(buffer1, sizeof(float), NAZ0[i], fp_out);
          for (j = 0; j < NAZ0[i]; j++)
          {
            //            fread(&tmp_stream1,sizeof(tmp_stream1),1,fp_out);
            hzen_ps->Fill((i + 0.5) * hori_sys_width, buffer1[j]);
            //if(tmp_stream1>1e-6) printf("%f %f\n",(i+0.5)*hori_sys_width,tmp_stream1);
            hazim_ps->Fill((j + 0.5) * 360. / NAZ0[i], buffer1[j]);
            fread(&ra, sizeof(ra), 1, fp_grid_ra_dec);
            fread(&dec, sizeof(dec), 1, fp_grid_ra_dec);
            ra = (ra + isr * equa_sys_width) < 360.0 ? (ra + isr * equa_sys_width) : (ra + isr * equa_sys_width - 360.0);
            hra_ps->Fill(ra, buffer1[j]);
            hdec_ps->Fill(dec, buffer1[j]);
            hsig_ps->Fill((buffer1[j]) / sqrt(tmp_stream2));
          }
          //          fread(&tmp_stream2,sizeof(tmp_stream2),1,fp_nb);
          //          for(fseek(fp_out,-NAZ0[i]*sizeof(tmp_stream1),SEEK_CUR),j=0;j<NAZ0[i];j++) {
          //            fread(&tmp_stream1,sizeof(tmp_stream1),1,fp_out);
          //            hsig_ps->Fill((tmp_stream1-tmp_stream2)/sqrt(tmp_stream2));
          //          }
        }
        if (ismon == 1)
          printf("writing data!   %3d%% \r", int(isr * 100 / NLST));
        //        printf("writing data! \t isr = %d in %d \r",isr,NLST);
      }
      if (ismon == 1)
        printf("\n");
      hfile.Write();
      hfile.Close();
    }
  }

  fclose(fp_nb);
  fclose(fp_grid_ra_dec);
  fclose(fp_log);
  fclose(fp_out);

  if (ismon == 0)
  {
    strcpy(argv[1], filename);
    remove(strcat(argv[1], "_tmp_nb"));
    strcpy(argv[1], filename);
    remove(strcat(argv[1], "_tmp_grid_ra_dec"));
  }
  //  delete []n;
  //  delete []nb_tmp;
  return 0;
}

/////////////////////////////////////////////////////////////////
//////////////functions used!!!/////////////////////////////////
int rand_possion(double lambda)
{
  double x = -1, u;
  double log1, log2;
  log1 = 0;
  log2 = -lambda;
  do
  {
    u = rand() / (RAND_MAX + 1.0);
    log1 += log(u);
    x++;
  } while (log1 >= log2);
  return x;
}

double rand_gauss(double mu, double sigma)
{
  static double U, V;
  static int phase = 0;
  double Z;

  if (phase == 0)
  {
    U = rand() / (RAND_MAX + 1.0);
    V = rand() / (RAND_MAX + 1.0);
    Z = sqrt(-2.0 * log(U)) * sin(2.0 * PI * V);
  }
  else
  {
    Z = sqrt(-2.0 * log(U)) * cos(2.0 * PI * V);
  }

  phase = 1 - phase;
  return Z * sigma + mu;
}

double bigaussian_mean(double theta, double sigma)
{
  double thetamin = theta - hori_sys_width / sqrt(PI) > 0 ? theta - hori_sys_width / sqrt(PI) : 0;
  double thetamax = theta + hori_sys_width / sqrt(PI);
  // use +/- hori_sys_width/sqrt(PI) to make sure the event number in the central bin is correct
  // since the bin is a square and the integration is within a circle.
  return (exp(-thetamin * thetamin / (2 * sigma * sigma)) - exp(-thetamax * thetamax / (2 * sigma * sigma))) / (2. * PI * (cos(thetamin * D2R) - cos(thetamax * D2R)) * R2D * R2D);
}

double bigaussian(double theta, double sigma)
{
  return exp(-theta * theta / (2 * sigma * sigma)) / (2 * PI * sigma * sigma);
}

TFile *file = new TFile("effarea.root");

double wcda_eff_area_ga(double eg, double zen)
{
  double eff_area;
  int bin_num;

  eg = eg - 9.0;
  if (zen < 8.3)
    zen = 8.3;

  bin_num = ((TH1F *)file->Get("harea_p0"))->FindBin(eg, zen);
  eff_area = ((TH1F *)file->Get("harea_p0"))->GetBinContent(bin_num);
  return eff_area;
}

double wcda_eff_area_cr(double ep, double zen)
{
  double eff_area;
  int bin_num;

  ep = ep - 9.0;
  if (zen < 1.3)
    zen = 1.3;

  bin_num = ((TH1F *)file->Get("harea_p1"))->FindBin(ep, zen);
  eff_area = ((TH1F *)file->Get("harea_p1"))->GetBinContent(bin_num);
  return eff_area;
}

double wcda_filter_ga(double eg, double zen)
{
  double ratio;
  int bin_num;

  eg = eg - 9.0;
  if (zen < 8.3)
    zen = 8.3;

  bin_num = ((TH1F *)file->Get("hratio_p0"))->FindBin(eg, zen);
  ratio = ((TH1F *)file->Get("hratio_p0"))->GetBinContent(bin_num);
  return ratio;
}

double wcda_filter_cr(double ep, double zen)
{
  double ratio;
  int bin_num;

  ep = ep - 9.0;
  if (zen < 1.3)
    zen = 1.3;

  bin_num = ((TH1F *)file->Get("hratio_p1"))->FindBin(ep, zen);
  ratio = ((TH1F *)file->Get("hratio_p1"))->GetBinContent(bin_num);
  return ratio;
}

void e_median(double *emid_ga, double *emid_cr)
{
  double num_sum, num_tot, bin_width, low_edge;
  int bin_num, bin_tot;

  num_tot = ((TH1F *)file->Get("he_ga"))->Integral();
  low_edge = ((TH1F *)file->Get("he_ga"))->GetBinLowEdge(1);
  bin_width = ((TH1F *)file->Get("he_ga"))->GetBinWidth(1);
  bin_tot = ((TH1F *)file->Get("he_ga"))->FindLastBinAbove();

  num_sum = 0.;
  for (bin_num = 0; bin_num <= bin_tot; bin_num++)
  {
    num_sum += ((TH1F *)file->Get("he_ga"))->GetBinContent(bin_num);
    if (num_sum > (0.5 * num_tot))
      break;
  }
  *emid_ga = 9. + low_edge + bin_width * (bin_num);

  num_tot = ((TH1F *)file->Get("he_cr"))->Integral();
  low_edge = ((TH1F *)file->Get("he_cr"))->GetBinLowEdge(1);
  bin_width = ((TH1F *)file->Get("he_cr"))->GetBinWidth(1);
  bin_tot = ((TH1F *)file->Get("he_cr"))->FindLastBinAbove();

  num_sum = 0.;
  for (bin_num = 0; bin_num <= bin_tot; bin_num++)
  {
    num_sum += ((TH1F *)file->Get("he_cr"))->GetBinContent(bin_num);
    if (num_sum > (0.5 * num_tot))
      break;
  }
  *emid_cr = 9. + low_edge + bin_width * (bin_num);
}

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

TFile *file_ang = new TFile("angres.root");

double wcda_theta_dist(double theta)
{
  double orig_bin_cont;
  int bin_num;
  double ring_area;
  double theta_dist_density;

  bin_num = ((TH1F *)file_ang->Get("hpol"))->FindBin(theta);
  orig_bin_cont = ((TH1F *)file_ang->Get("hpol"))->GetBinContent(bin_num);

  ring_area = 2. * PI * (cos((bin_num - 1) * 0.01 * D2R) - cos(bin_num * 0.01 * D2R)) * R2D * R2D;
  theta_dist_density = orig_bin_cont / ring_area;
  if (bin_num > 700)
    theta_dist_density = 0.;
  return theta_dist_density;
}

double wcda_theta_dist1(double theta)
{
  double orig_bin_cont;
  int bin_num;
  int bin_num_le, bin_num_ue;
  double ring_area;
  double theta_dist_density;

  bin_num = ((TH1F *)file_ang->Get("hpol_rho"))->FindBin(theta);
  bin_num_le = (0 > (bin_num - hori_sys_width / (2 * 0.01))) ? 0 : (bin_num - hori_sys_width / (2 * 0.01));
  bin_num_ue = (700 < (bin_num + hori_sys_width / (2 * 0.01))) ? 700 : (bin_num + hori_sys_width / (2 * 0.01));
  //  orig_bin_cont = ((TH1F*)file_ang->Get("hpol_rho"))->GetBinContent(bin_num);
  orig_bin_cont = ((TH1F *)file_ang->Get("hpol_rho"))->Integral(bin_num_le, bin_num_ue) / (bin_num_ue - bin_num_le + 1);

  theta_dist_density = orig_bin_cont;
  if (bin_num > 700)
    theta_dist_density = 0.;
  return theta_dist_density;
}

void wcda_ang_maxsigma(double day_num, double *theta_maxsigma, double *maxsigma, double *keepratio_ga, double *keepratio_cr, double *num_ga, double *num_cr)
{
  double sigma[700];
  double n_ga, n_cr, n_cr_tot;
  double n_tmp = 0.;
  int n_maxsigma;

  n_cr_tot = ((TH1F *)file_ang->Get("hpol_cr_smth"))->Integral(1, 700);
  n_cr_tot *= day_num;
  for (int i = 1; i <= 700; i++)
  {
    n_ga = ((TH1F *)file_ang->Get("hpol_smth"))->Integral(1, i);
    n_ga *= day_num;
    //    n_cr = ((TH1F*)file_ang->Get("hpol_cr_smth"))->Integral(1,i);
    n_cr = n_cr_tot * (1. - cos(i * 0.01 * D2R)) / (1. - cos(7. * D2R));
    sigma[i - 1] = n_ga / sqrt(n_cr);
    if (sigma[i - 1] < n_tmp || n_cr < 1e-5)
      continue;
    n_tmp = sigma[i - 1];
    n_maxsigma = i;
  }

  *theta_maxsigma = n_maxsigma * 0.01;
  *maxsigma = n_tmp;
  *num_ga = ((TH1F *)file_ang->Get("hpol_smth"))->Integral(1, n_maxsigma);
  *num_cr = ((TH1F *)file_ang->Get("hpol_cr_smth"))->Integral(1, n_maxsigma);
  *num_ga *= day_num;
  *num_cr *= day_num;
  *keepratio_ga = (((TH1F *)file_ang->Get("hpol_smth"))->Integral(1, n_maxsigma)) / (((TH1F *)file_ang->Get("hpol_smth"))->Integral(1, 700));
  *keepratio_cr = (((TH1F *)file_ang->Get("hpol_cr_smth"))->Integral(1, n_maxsigma)) / (((TH1F *)file_ang->Get("hpol_cr_smth"))->Integral(1, 700));
}
