
#include "unif01.h"
#include "bbattery.h"

//unsigned int xorshift (void);
unsigned int dice(void);
double xorshift (void);
double MRG32k3a (void);


unif01_Gen* uweyl_CreateWeyl(double alpha, long n);
unif01_Gen* uweyl_CreateNWeyl(double alpha, long n);
unif01_Gen* uweyl_CreateSNWeyl(long m, double alpha, long n);

//#define K 0.414213562373095
#define K 0.7639320225002103035908
#define M 0x12345
double state = 0;
uint32_t n = 0;

// should fail all 15 SC
double weyl(void)
{
  state = state + K;
  state = state - (uint32_t)state;

  return state;
}

// fails 14...built-in fails 9...humm
double weyln(void)
{
  double r;

  n++;

  r  = K*n;
  r -= (uint32_t)r;
  r *= n;
  r -= (uint32_t)r;

  return r;
}

// fails SC#6...built-in passes all SC.  Humm
double weylsn()
{
  double X,Z;

  n++;
  X = n * K;
  X -= (uint32_t)X;
  X *= n;
  X -= (uint32_t)X;
  X = X * M + 0.5;
  Z = X * K;
  Z -= (uint32_t)Z;
  Z *= X;
  Z -= (uint32_t)Z;
  return Z;
}


#include <math.h>
double mod289(double x)
{
 return x - floor(x * (1.0/289.0)) * 289.0;
}

double permute(void)
{
  double x = state;

  x = mod289(((x*34.0)+1.0)*x);
  state = x;
  x *= (1.0/41.0);

  return x-(uint32_t)x;
}

int main (void) 
{
  unif01_Gen *gen;

  gen = unif01_CreateExternGen01("weyl", permute);
  //gen = uweyl_CreateWeyl(K,0);
  //gen = uweyl_CreateNWeyl(K,0);
  //gen = uweyl_CreateSNWeyl(M,K,0);
  bbattery_SmallCrush(gen);
  //bbattery_Crush(gen);
  //bbattery_BigCrush (gen);
  unif01_DeleteExternGenBits(gen);
  
  return 0;
}
