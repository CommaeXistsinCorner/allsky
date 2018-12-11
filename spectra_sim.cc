#include<stdio.h>
#include<math.h>
#include "ps_tevcat_v2.cc"

double CR_intens_p1(double ep)
{
  double k1 = 7860.; 
  double a1 = 2.66;
  double ec1 = 4e6;
  double k2 = 20.;
  double a2 = 2.4;
  double ec2 = 3e7;
  double k3 = 1.7;
  double a3 = 2.4;
  double ec3 = 2e9;
  double z = 1.;
  double crspec;

  ep = pow(10.0,ep-9.0);
  crspec = 1e-13 * (k1 * pow(ep,-a1) * exp(-ep/(z * ec1)) + k2 * pow(ep,-a2) * exp(-ep/(z * ec2)) + k3 * pow(ep,-a3) * exp(-ep/(z * ec3)));

// ev^-1 cm^-2 s^-1 sr^-1,  energy unit eg is log10(e)
  return crspec;
}

double CR_intens_p2(double ep)
{ 
  double k1 = 3550.;
  double a1 = 2.58;
  double ec1 = 4e6;
  double k2 = 20.;
  double a2 = 2.4;
  double ec2 = 3e7;
  double k3 = 1.7;
  double a3 = 2.4;
  double ec3 = 2e9;
  double z = 2.;
  double crspec;
    
  ep = pow(10.0,ep-9.0);
  crspec = 1e-13 * (k1 * pow(ep,-a1) * exp(-ep/(z * ec1)) + k2 * pow(ep,-a2) * exp(-ep/(z * ec2)) + k3 * pow(ep,-a3) * exp(-ep/(z * ec3)));
// ev^-1 cm^-2 s^-1 sr^-1
  return crspec;
}

double CR_intens_p3(double ep)
{
  double k1 = 2200.;
  double a1 = 2.63;
  double ec1 = 4e6;
  double k2 = 13.4;
  double a2 = 2.4;
  double ec2 = 3e7;
  double k3 = 1.14;
  double a3 = 2.4;
  double ec3 = 2e9;
  double z = 7.;
  double crspec;

  ep = pow(10.0,ep-9.0);
  crspec = 1e-13 * (k1 * pow(ep,-a1) * exp(-ep/(z * ec1)) + k2 * pow(ep,-a2) * exp(-ep/(z * ec2)) + k3 * pow(ep,-a3) * exp(-ep/(z * ec3)));// ev^-1 cm^-2 s^-1 sr^-1
  return crspec;
}

double CR_intens_p4(double ep)
{
  double k1 = 1430.;
  double a1 = 2.67;
  double ec1 = 4e6;
  double k2 = 13.4;
  double a2 = 2.4;
  double ec2 = 3e7;
  double k3 = 1.14;
  double a3 = 2.4;
  double ec3 = 2e9;
  double z = 13.;
  double crspec;

  ep = pow(10.0,ep-9.0);
  crspec = 1e-13 * (k1 * pow(ep,-a1) * exp(-ep/(z * ec1)) + k2 * pow(ep,-a2) * exp(-ep/(z * ec2)) +k3 * pow(ep,-a3) * exp(-ep/(z * ec3)));// ev^-1 cm^-2 s^-1 sr^-1
  return crspec;
}

double CR_intens_p5(double ep)
{
  double k1 = 2120.;
  double a1 = 2.63;
  double ec1 = 4e6;
  double k2 = 13.4;
  double a2 = 2.4;
  double ec2 = 3e7;
  double k3 = 1.14;
  double a3 = 2.4;
  double ec3 = 2e9;
  double z = 26.;
  double crspec;

  ep = pow(10.0,ep-9.0);
  crspec = 1e-13 * (k1 * pow(ep,-a1) * exp(-ep/(z * ec1)) + k2 * pow(ep,-a2) * exp(-ep/(z * ec2)) +k3 * pow(ep,-a3) * exp(-ep/(z * ec3)));// ev^-1 cm^-2 s^-1 sr^-1
  return crspec;
}

double CR_intens_all(double ep)
{
  return CR_intens_p1(ep)+CR_intens_p2(ep)+CR_intens_p3(ep)+CR_intens_p4(ep)+CR_intens_p5(ep);
}

double CR_itgr_intens_all(double ep1, double ep2)
{
  int i;
  int n = int(100*(ep2-ep1));
  double area;
  double sum = 0.;
  for(i=0;i<n;i++)
  {
    area = CR_intens_all(ep1 + 0.005 + 0.01*n) * 0.023293 * pow(10.0, ep1 + 0.01*n);
    sum += area;
  }
  return sum;
}

double CRAB[2]={83.629583,22.014472};  //HEGRA

double CRAB_intens(double eg)
{
  double gaspec;
  gaspec  = 9.0677 - 2.63 * eg - 0.43429 * pow(10.0,eg - 13.0);
  gaspec = pow(10.0,gaspec);
// ev^-1 cm^-2 s^-1 
  return gaspec;
}

double CRAB_intens_nocut(double eg)
{
  double gaspec;
  gaspec  = 9.0677 - 2.63 * eg;
  gaspec = pow(10.0,gaspec);
// ev^-1 cm^-2 s^-1 
  return gaspec;
}

double CRAB_intens_HEGRA(double eg)
{ 
  double gaspec;
  gaspec  = - 2.62 * (eg - 12.0);
  gaspec = 2.83e-23 * pow(10.0,gaspec);
// ev^-1 cm^-2 s^-1 
 return gaspec;
}

double CRAB_itgr_intens_HEGRA(double eg1, double eg2)
{
// cm^-2 s^-1 
  return 4.8114e8 * (exp(-3.73019*eg1)-exp(-3.73019*eg2));
}

double SOURCE_intens(int i, double eg)
{
  eg = pow(10.0,eg - 12.0);
  double PHI0 = PS_TEVCAT[i][2];
  double E0 = PS_TEVCAT[i][3];
  double INDEX = 0. - PS_TEVCAT[i][4];
  double ECUT = PS_TEVCAT[i][5];
  if(PS_TEVCAT[i][5]<1e-3) {return 1e-12 * PHI0 * pow((eg/E0),INDEX);}
  else {return 1e-12 * PHI0 * pow((eg/E0),INDEX) * exp(- eg/ECUT);}
}

double TEST_SOURCE_intens(double eg)
{
  eg = pow(10.0,eg - 12.0);
  double PHI0 = PS_TEST[2];
  double E0 = PS_TEST[3];
  double INDEX = 0. - PS_TEST[4];
  double ECUT = PS_TEST[5];
  if(PS_TEST[5]<1e-3) {return 1e-12 * PHI0 * pow((eg/E0),INDEX);}
  else {return 1e-12 * PHI0 * pow((eg/E0),INDEX) * exp(- eg/ECUT);}
}
