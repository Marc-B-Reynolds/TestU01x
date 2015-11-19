#include <stdio.h>
#include <string.h>
#include "TestU01.h"
#include "util.h"
#include "adaptive_crush.h"
/* Compile example
   $ gcc example_adaptive_crush.c adaptive_crush.c -o example_adaptive_crush -ltestu01 -lprobdist -lm */


int main (void)
{
  swrite_Basic = FALSE;
  unif01_Gen *gen;

  gen = ulcg_CreateLCG2e31 (1103515245, 12345, 3);
  bbattery_Adaptive_Crush (gen);

  return 0;
}
