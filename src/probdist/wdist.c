/*************************************************************************\
* Copyright (c) 2002 Pierre L'Ecuyer, DIRO, Université de Montréal.
* SEE: LICENSE
\*************************************************************************/

#include "wdist.h"
#include "fdist.h"

double wdist_Normal(double Junk[], double x)
{
  return fdist_Normal2 (x);
}

double wdist_ChiSquare (double W[], double x)
{
  long N = (long) W[0];
  return fdist_ChiSquare2 (N, 12, x);
}

double wdist_Unif (double Junk[], double x)
{
  return fdist_Unif(x);
}

