#ifndef _ANA_WID_SIM_H_
#define _ANA_WID_SIM_H_

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#ifndef YBJ_GROUP_NH
#define D2R 0.017453292519943296
#define R2D 57.295779513082322
#define PI 3.14159265358979323846
#endif

#define LHAASO_Lati 29.3586  /*  latitude  */
#define LHAASO_Lngi 100.1372 /* longtitude */

//#define  LHAASO_Lati 19.0          /*  latitude  */
//#define  LHAASO_Lngi 100.1372          /* longtitude */

//#define  LHAASO_Lati 30.102          /*  latitude  */
//#define  LHAASO_Lngi 90.522          /* longtitude */

#define YBJ_Lati 30.102                  /*  latitude  */
#define YBJ_Lngi 90.522 /* longtitude */ //ybj!! not lhaaso!!

#define PS_INPUT 1

#define equa_sys_width 0.1
#define hori_sys_width 0.1
#define zex 60
//#define TEMP_FILE "/afs/ihep.ac.cn/users/c/changxc/eos/sci_code/output"

#endif

const int NZEAZ = 1031622;

const int NZE = int(zex / hori_sys_width);
const int NRA = int(360.0 / equa_sys_width);
//const double DEC_MIN = LHAASO_Lati - zex;
//const double DEC_MAX = (LHAASO_Lati + zex)<90.0?(LHAASO_Lati + zex):90;
const int DEC_MIN = int(LHAASO_Lati - zex - 1);
const int DEC_MAX = (LHAASO_Lati + zex) < 90.0 ? int(LHAASO_Lati + zex) : 90;
const int NDEC = int((DEC_MAX - DEC_MIN) / equa_sys_width);
const int NLST = NRA;
// question, why? NLST could be smaller or even arbitrary?
const double dt_time = 1.0 / NLST;
const double beishu = NRA / 360.0;

int NAZ0[NZE];
int NAZ[NZE];

//const int NZEAZ = nzeaz();
//
//int nzeaz()
//{
//  int tmp_k=0;
//  for(int tmp_i=0; tmp_i<NZE; tmp_i++)  {
//    NAZ0[tmp_i] = ceil( 360.0 * sin((0.5+tmp_i)*hori_sys_width*PI/180.0) / hori_sys_width );
//    NAZ[tmp_i] = tmp_k;
//    tmp_k += NAZ0[tmp_i];
//  }
//  return tmp_k;
//}
