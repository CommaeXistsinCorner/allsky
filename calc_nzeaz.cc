#ifndef ANA_WID76_H
#define ANA_WID76_H 1

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#ifndef YBJ_GROUP_NH
#define D2R 0.017453292519943296
#define R2D 57.295779513082322
#define PI 3.14159265358979323846
#endif

#endif

int main(int argc, char *argv[])
{
  if (argc != 3)
  {
    printf("%s  bin_width  zen_max\n", argv[0]);
    exit(0);
  }
  double equa_sys_width = atof(argv[1]);
  double hori_sys_width = atof(argv[1]);
  double zex = atof(argv[2]);

  const int NZE = int(zex / hori_sys_width);

  int NAZ0[NZE];
  int NAZ[NZE];

  int tmp_k = 0;
  for (int tmp_i = 0; tmp_i < NZE; tmp_i++)
  {
    NAZ0[tmp_i] = ceil(360.0 * sin((0.5 + tmp_i) * hori_sys_width * PI / 180.0) / hori_sys_width);
    NAZ[tmp_i] = tmp_k;
    tmp_k += NAZ0[tmp_i];
  }

  const int NZEAZ = tmp_k;

  printf("NZEAZ= %d \n", NZEAZ);
  return 0;
}
