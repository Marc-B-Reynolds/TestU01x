#include "TestU01.h"
#include "util.h"
#include "bbattery.h"
#include "smultin.h"
#include "sknuth.h"
#include "smarsa.h"
#include "snpair.h"
#include "svaria.h"
#include "sstring.h"
#include "swalk.h"
#include "scomp.h"
#include "sspectral.h"
#include "swrite.h"
#include "sres.h"
#include "unif01.h"
#include "ufile.h"

#include "gofs.h"
#include "gofw.h"
#include "fdist.h"
#include "fbar.h"
#include "num.h"
#include "chrono.h"

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <limits.h>

#include "adaptive_crush.h"

#define LEN 120
#define NAMELEN 30
#define NDIM 200                  /* Dimension of extern arrays */
#define THOUSAND 1000
#define MILLION (THOUSAND * THOUSAND)
#define BILLION (THOUSAND * MILLION)

/* Parameters of Adaptive Crush tests */
#define ADAPTIVE_CRUSH_NUM 96
#define ALPHA 0.1
#define BETA 0.00000001
#define ITERATION_LIMIT 1
#define BREAK (ITERATION_LIMIT+1)

double bbattery_pVal[1 + NDIM] = { 0 };
char *bbattery_TestNames[1 + NDIM] = { 0 };
int bbattery_NTests;

void DetectIteration (double, long *, int *);


/* Gives the test number as enumerated in bbattery.tex. Some test applies
   more than one test, so the array of p-values does not correspond with
   the test number in the doc. */
static int TestNumber[1 + NDIM] = { 0 };




/*-------------------------------- Functions ------------------------------*/


static void GetName (unif01_Gen * gen, char *genName)
{
   char *p;
   int len1, len2;

   if (NULL == gen) {
      genName[0] = '\0';
      return;
   }

   /* Print only the generator name, without the parameters or seeds. */
   /* The parameters start after the first blank; name ends with ':' */
   genName[LEN] = '\0';
   len1 = strcspn (gen->name, ":");
   len1 = util_Min (LEN, len1);
   strncpy (genName, gen->name, (size_t) len1);
   genName[len1] = '\0';
   /* For Filters or Combined generators */
   p = strstr (&gen->name[1 + len1], "unif01");
   while (p != NULL) {
      len1 += 2;
      if (len1 >= LEN)
         return;
      strcat (genName, ", ");
      len2 = strcspn (p, " \0");
      len2 = util_Min (LEN - len1, len2);
      if (len2 <= 0)
         return;
      strncat (genName, p, (size_t) len2);
      len1 = strlen (genName);
      genName[len1] = '\0';
      p += len2;
      p = strstr (p, "unif01");
   }
}


/*=========================================================================*/

static void WritepVal (double p)
/*
 * Write a p-value with a nice format.
 */
{
   if (p < gofw_Suspectp) {
      gofw_Writep0 (p);

   } else if (p > 1.0 - gofw_Suspectp) {
      if (p >= 1.0 - gofw_Epsilonp1) {
         printf (" 1 - eps1");
      } else if (p >= 1.0 - 1.0e-4) {
         printf (" 1 - ");
         num_WriteD (1.0 - p, 7, 2, 2);
         /* printf (" 1 - %.2g ", 1.0 - p); */
      } else if (p >= 1.0 - 1.0e-2)
         printf ("  %.4f ", p);
      else
         printf ("   %.2f", p);
   }
}


/*=========================================================================*/

static void WriteReport (
   char *genName,                 /* Generator or file name */
   char *batName,                 /* Battery name */
   int N,                         /* Max. number of tests */
   double pVal[],                 /* p-values of the tests */
   chrono_Chrono * Timer,         /* Timer */
   lebool Flag,                  /* = TRUE for a file, FALSE for a gen */
   lebool VersionFlag,           /* = TRUE: write the version number */
   double nb                      /* Number of bits in the random file */
   )
{
   int j, co;

   printf ("\n========= Summary results of ");
   printf ("%s", batName);
   printf (" =========\n\n");
   if (VersionFlag)
     ; /* printf (" Version:          %s\n", PACKAGE_STRING); */
   /* Original code writes the version of the Package.
      The paremeter "PACKAGE_STRING" requires "config.h" in the
      same directory. We avoid to use this header file. */
   if (Flag)
      printf (" File:             ");
   else
      printf (" Generator:        ");
   printf ("%s", genName);
   if (nb > 0)
      printf ("\n Number of bits:   %.0f", nb);
   co = 0;
   /* Some of the tests have not been done: their pVal[j] < 0. */
   for (j = 0; j < N; j++) {
      if (pVal[j] >= 0.0)
         co++;
   }
   printf ("\n Number of statistics:  %1d\n", co);
   printf (" Total CPU time:   ");
   chrono_Write (Timer, chrono_hms);

   co = 0;
   for (j = 0; j < N; j++) {
      if (pVal[j] < 0.0)          /* That test was not done: pVal = -1 */
         continue;
      if ((pVal[j] < gofw_Suspectp) || (pVal[j] > 1.0 - gofw_Suspectp)) {
         co++;
         break;
      }
   }
   if (co == 0) {
      printf ("\n\n All tests were passed\n\n\n\n");
      return;
   }

   if (gofw_Suspectp >= 0.01)
      printf ("\n The following tests gave p-values outside [%.4g, %.2f]",
         gofw_Suspectp, 1.0 - gofw_Suspectp);
   else if (gofw_Suspectp >= 0.0001)
      printf ("\n The following tests gave p-values outside [%.4g, %.4f]",
         gofw_Suspectp, 1.0 - gofw_Suspectp);
   else if (gofw_Suspectp >= 0.000001)
      printf ("\n The following tests gave p-values outside [%.4g, %.6f]",
         gofw_Suspectp, 1.0 - gofw_Suspectp);
   else
      printf ("\n The following tests gave p-values outside [%.4g, %.14f]",
         gofw_Suspectp, 1.0 - gofw_Suspectp);
   printf (":\n (eps  means a value < %6.1e)", gofw_Epsilonp);
   printf (":\n (eps1 means a value < %6.1e)", gofw_Epsilonp1);
   printf (":\n\n       Test                          p-value\n");
   printf (" ----------------------------------------------\n");

   co = 0;
   for (j = 0; j < N; j++) {
      if (pVal[j] < 0.0)          /* That test was not done: pVal = -1 */
         continue;
      if ((pVal[j] >= gofw_Suspectp) && (pVal[j] <= 1.0 - gofw_Suspectp))
         continue;                /* That test was passed */
      printf (" %2d ", TestNumber[j]);
      printf (" %-30s", bbattery_TestNames[j]);
      WritepVal (pVal[j]);
      printf ("\n");
      co++;
   }

   printf (" ----------------------------------------------\n");
   if (co < N - 1) {
      printf (" All other tests were passed\n");
   }
   printf ("\n\n\n");
}


/*=========================================================================*/


/*=========================================================================*/

static void InitBat (void)
/*
 * Initializes the battery of tests: sets all p-values to -1.
 */
{
   int j;
   static int flag = 0;
   for (j = 0; j < NDIM; j++)
      bbattery_pVal[j] = -1.0;
   if (0 == flag) {
      flag++;
      for (j = 0; j < NDIM; j++)
         bbattery_TestNames[j] = util_Calloc (LEN + 1, sizeof (char));
   }
}


/*=========================================================================*/


/*=========================================================================*/

static void Adaptive_Crush (unif01_Gen * gen, int Rep[])
/*
 * A battery of stringent statistical tests for Random Number Generators
 * used in simulation.
 * Rep[i] gives the number of times that test i will be done. The default
 * values are Rep[i] = 1 for all i.
 */
{
   const int s = 30;
   const int r = 0;
   int i;
   chrono_Chrono *Timer;
   char genName[LEN + 1] = "";
   int j = -1;
   int j2 = 0;
   long size;

   Timer = chrono_Create ();
   InitBat ();
   if (swrite_Basic) {
      printf ("xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx\n"
         "                 Starting Adaptive Crush\n"
         "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx\n\n\n");
   }
   {
      sres_Basic *res;
      res = sres_CreateBasic ();
      ++j2;
      i = 0;
      size = 500 * MILLION;
      while (i < ITERATION_LIMIT){
         smarsa_SerialOver (gen, res, 1, size, 0, 4096, 2);
	 DetectIteration (res->pVal2[gofw_Mean], &size, &i);
      }
      bbattery_pVal[++j] = res->pVal2[gofw_Mean];
      TestNumber[j] = j2;
      strcpy (bbattery_TestNames[j], "SerialOver, t = 2");

      ++j2;
      i = 0;
      size = 300 * MILLION;
      while (i < ITERATION_LIMIT){
	smarsa_SerialOver (gen, res, 1, size, 0, 64, 4);
	DetectIteration (res->pVal2[gofw_Mean], &size, &i);
      }
      bbattery_pVal[++j] = res->pVal2[gofw_Mean];
      TestNumber[j] = j2;
      strcpy (bbattery_TestNames[j], "SerialOver, t = 4");

      sres_DeleteBasic (res);
   }
   {
      sres_Chi2 *res;
      res = sres_CreateChi2 ();

      ++j2;
      i = 0;
      size = 40 * MILLION;
      while (i < ITERATION_LIMIT){
	sknuth_SimpPoker (gen, res, 1, size, 0, 16, 16);
	DetectIteration (res->pVal2[gofw_Mean], &size, &i);
      }
      bbattery_pVal[++j] = res->pVal2[gofw_Mean];
      TestNumber[j] = j2;
      strcpy (bbattery_TestNames[j], "SimpPoker, d = 16");

      ++j2;
      i = 0;
      size = 40 * MILLION;
      while (i < ITERATION_LIMIT){
	sknuth_SimpPoker (gen, res, 1, size, 26, 16, 16);
	DetectIteration (res->pVal2[gofw_Mean], &size, &i);
      }
      bbattery_pVal[++j] = res->pVal2[gofw_Mean];
      TestNumber[j] = j2;
      strcpy (bbattery_TestNames[j], "SimpPoker, d = 16");

      ++j2;
      i = 0;
      size = 10 * MILLION;
      while (i < ITERATION_LIMIT){
	sknuth_SimpPoker (gen, res, 1, size, 0, 64, 64);
	DetectIteration (res->pVal2[gofw_Mean], &size, &i);
      }
      bbattery_pVal[++j] = res->pVal2[gofw_Mean];
      TestNumber[j] = j2;
      strcpy (bbattery_TestNames[j], "SimpPoker, d = 64");

      ++j2;
      i = 0;
      size = 10 * MILLION;
      while (i < ITERATION_LIMIT){
	sknuth_SimpPoker (gen, res, 1, size, 24, 64, 64);
	DetectIteration (res->pVal2[gofw_Mean], &size, &i);
      }
      bbattery_pVal[++j] = res->pVal2[gofw_Mean];
      TestNumber[j] = j2;
      strcpy (bbattery_TestNames[j], "SimpPoker, d = 64");

      ++j2;
      i = 0;
      size = 40 * MILLION;
      while (i < ITERATION_LIMIT){
	sknuth_CouponCollector (gen, res, 1, size, 0, 4);
	DetectIteration (res->pVal2[gofw_Mean], &size, &i);
      }
      bbattery_pVal[++j] = res->pVal2[gofw_Mean];
      TestNumber[j] = j2;
      strcpy (bbattery_TestNames[j], "CouponCollector, d = 4");

      ++j2;
      i = 0;
      size = 40 * MILLION;
      while (i < ITERATION_LIMIT){
	sknuth_CouponCollector (gen, res, 1, size, 28, 4);
	DetectIteration (res->pVal2[gofw_Mean], &size, &i);
      }
      bbattery_pVal[++j] = res->pVal2[gofw_Mean];
      TestNumber[j] = j2;
      strcpy (bbattery_TestNames[j], "CouponCollector, d = 4");

      ++j2;
      i = 0;
      size = 10 * MILLION;
      while (i < ITERATION_LIMIT){
	sknuth_CouponCollector (gen, res, 1, size, 0, 16);
	DetectIteration (res->pVal2[gofw_Mean], &size, &i);
      }
      bbattery_pVal[++j] = res->pVal2[gofw_Mean];
      TestNumber[j] = j2;
      strcpy (bbattery_TestNames[j], "CouponCollector, d = 16");

      ++j2;
      i = 0;
      size = 10 * MILLION;
      while (i < ITERATION_LIMIT){
	sknuth_CouponCollector (gen, res, 1, size, 26, 16);
	DetectIteration (res->pVal2[gofw_Mean], &size, &i);
      }
      bbattery_pVal[++j] = res->pVal2[gofw_Mean];
      TestNumber[j] = j2;
      strcpy (bbattery_TestNames[j], "CouponCollector, d = 16");

      ++j2;
      i = 0;
      size = 100 * MILLION;
      while (i < ITERATION_LIMIT){
	sknuth_Gap (gen, res, 1, size, 0, 0.0, 0.125);
	DetectIteration (res->pVal2[gofw_Mean], &size, &i);
      }
      bbattery_pVal[++j] = res->pVal2[gofw_Mean];
      TestNumber[j] = j2;
      strcpy (bbattery_TestNames[j], "Gap, r = 0");

      ++j2;
      i = 0;
      size = 100 * MILLION;
      while (i < ITERATION_LIMIT){
	sknuth_Gap (gen, res, 1, size, 27, 0.0, 0.125);
	DetectIteration (res->pVal2[gofw_Mean], &size, &i);
      }
      bbattery_pVal[++j] = res->pVal2[gofw_Mean];
      TestNumber[j] = j2;
      strcpy (bbattery_TestNames[j], "Gap, r = 27");

      ++j2;
      i = 0;
      size = 5 * MILLION;
      while (i < ITERATION_LIMIT){
	sknuth_Gap (gen, res, 1, size, 0, 0.0, 1.0/256.0);
	DetectIteration (res->pVal2[gofw_Mean], &size, &i);
      }
      bbattery_pVal[++j] = res->pVal2[gofw_Mean];
      TestNumber[j] = j2;
      strcpy (bbattery_TestNames[j], "Gap, r = 0");

      ++j2;
      i = 0;
      size = 5 * MILLION;
      while (i < ITERATION_LIMIT){
	sknuth_Gap (gen, res, 1, 5 * MILLION, 22, 0.0, 1.0/256.0);
	DetectIteration (res->pVal2[gofw_Mean], &size, &i);
      }
      bbattery_pVal[++j] = res->pVal2[gofw_Mean];
      TestNumber[j] = j2;
      strcpy (bbattery_TestNames[j], "Gap, r = 22");

      ++j2;
      i = 0;
      size = 500 * MILLION;
      while (i < ITERATION_LIMIT){
	sknuth_Run (gen, res, 1, size, 0, TRUE);
	DetectIteration (res->pVal2[gofw_Mean], &size, &i);
      }
      bbattery_pVal[++j] = res->pVal2[gofw_Mean];
      TestNumber[j] = j2;
      strcpy (bbattery_TestNames[j], "Run of U01, r = 0");

      ++j2;
      i = 0;
      size = 500 * MILLION;
      while (i < ITERATION_LIMIT){
	sknuth_Run (gen, res, 1, size, 15, FALSE);
	DetectIteration (res->pVal2[gofw_Mean], &size, &i);
      }
      bbattery_pVal[++j] = res->pVal2[gofw_Mean];
      TestNumber[j] = j2;
      strcpy (bbattery_TestNames[j], "Run of U01, r = 15");

      ++j2;
      i = 0;
      size = 50 * MILLION;
      while (i < ITERATION_LIMIT){
	sknuth_Permutation (gen, res, 1, size, 0, 10);
	DetectIteration (res->pVal2[gofw_Mean], &size, &i);
      }
      bbattery_pVal[++j] = res->pVal2[gofw_Mean];
      TestNumber[j] = j2;
      strcpy (bbattery_TestNames[j], "Permutation, r = 0");

      ++j2;
      i = 0;
      size = 50 * MILLION;
      while (i <ITERATION_LIMIT){
	sknuth_Permutation (gen, res, 1, size, 15, 10);
	DetectIteration (res->pVal2[gofw_Mean], &size, &i);
      }
      bbattery_pVal[++j] = res->pVal2[gofw_Mean];
      TestNumber[j] = j2;
      strcpy (bbattery_TestNames[j], "Permutation, r = 15");

      sres_DeleteChi2 (res);
   }

   {
      sknuth_Res1 *res;
      res = sknuth_CreateRes1 ();

      ++j2;
      i = 0;
      size = 10 * MILLION;
      while (i < ITERATION_LIMIT){
	sknuth_MaxOft (gen, res, 1, size, 0, MILLION / 10, 5);
	DetectIteration (res->Chi->pVal2[gofw_Mean], &size, &i);
      }
      bbattery_pVal[++j] = res->Chi->pVal2[gofw_Mean];
      TestNumber[j] = j2;
      strcpy (bbattery_TestNames[j], "MaxOft, t = 5");

      ++j2;
      i = 0;
      size = 10 * MILLION;
      while (i < ITERATION_LIMIT){
	sknuth_MaxOft (gen, res, 1, size, 0, MILLION / 10, 5);
	DetectIteration (res->Bas->pVal2[gofw_Mean], &size, &i);
      }
      bbattery_pVal[++j] = res->Bas->pVal2[gofw_Mean];
      TestNumber[j] = j2;
      strcpy (bbattery_TestNames[j], "MaxOft AD, t = 5");

      ++j2;
      i = 0;
      size = 10 * MILLION;
      while (i < ITERATION_LIMIT){
	sknuth_MaxOft (gen, res, 1, size, 0, MILLION / 10, 10);
	DetectIteration (res->Chi->pVal2[gofw_Mean], &size, &i);
      }
      bbattery_pVal[++j] = res->Chi->pVal2[gofw_Mean];
      TestNumber[j] = j2;
      strcpy (bbattery_TestNames[j], "MaxOft, t = 10");

      ++j2;
      i = 0;
      size = 10 * MILLION;
      while (i < ITERATION_LIMIT){
	sknuth_MaxOft (gen, res, 1, size, 0, MILLION / 10, 10);
	DetectIteration (res->Bas->pVal2[gofw_Mean], &size, &i);
      }
      bbattery_pVal[++j] = res->Bas->pVal2[gofw_Mean];
      TestNumber[j] = j2;
      strcpy (bbattery_TestNames[j], "MaxOft AD, t = 10");

      ++j2;
      i = 0;
      size = 10 * MILLION;
      while (i < ITERATION_LIMIT){
	sknuth_MaxOft (gen, res, 1, size, 0, MILLION / 10, 20);
	DetectIteration (res->Chi->pVal2[gofw_Mean], &size, &i);
      }
      bbattery_pVal[++j] = res->Chi->pVal2[gofw_Mean];
      TestNumber[j] = j2;
      strcpy (bbattery_TestNames[j], "MaxOft, t = 20");

      i = 0;
      size = 10 * MILLION;
      while (i < ITERATION_LIMIT){
	sknuth_MaxOft (gen, res, 1, size, 0, MILLION / 10, 20);
	DetectIteration (res->Bas->pVal2[gofw_Mean], &size, &i);
      }
      bbattery_pVal[++j] = res->Bas->pVal2[gofw_Mean];
      TestNumber[j] = j2;
      strcpy (bbattery_TestNames[j], "MaxOft AD, t = 20");

      ++j2;
      i = 0;
      size = 10 * MILLION;
      while (i < ITERATION_LIMIT){
	sknuth_MaxOft (gen, res, 1, size, 0, MILLION / 10, 30);
	DetectIteration (res->Chi->pVal2[gofw_Mean], &size, &i);
      }
      bbattery_pVal[++j] = res->Chi->pVal2[gofw_Mean];
      TestNumber[j] = j2;
      strcpy (bbattery_TestNames[j], "MaxOft, t = 30");

      i = 0;
      size = 10 * MILLION;
      while (i < ITERATION_LIMIT){
	sknuth_MaxOft (gen, res, 1, size, 0, MILLION / 10, 30);
	DetectIteration (res->Bas->pVal2[gofw_Mean], &size, &i);
      }
      bbattery_pVal[++j] = res->Bas->pVal2[gofw_Mean];
      TestNumber[j] = j2;
      strcpy (bbattery_TestNames[j], "MaxOft AD, t = 30");

      sknuth_DeleteRes1 (res);
   }
   {
      sres_Basic *res;
      res = sres_CreateBasic ();

      ++j2;
      i = 0;
      size = 10 * MILLION;
      while (i < ITERATION_LIMIT){
	svaria_SampleProd (gen, res, 1, size, 0, 10);
	DetectIteration (res->pVal2[gofw_Mean], &size, &i);
      }
      bbattery_pVal[++j] = res->pVal2[gofw_Mean];
      TestNumber[j] = j2;
      strcpy (bbattery_TestNames[j], "SampleProd, t = 10");

      ++j2;
      i = 0;
      size = 10 * MILLION;
      while (i < ITERATION_LIMIT){
	svaria_SampleProd (gen, res, 1, size, 0, 30);
	DetectIteration (res->pVal2[gofw_Mean], &size, &i);
      }
      bbattery_pVal[++j] = res->pVal2[gofw_Mean];
      TestNumber[j] = j2;
      strcpy (bbattery_TestNames[j], "SampleProd, t = 30");

      ++j2;
      i = 0;
      size = 500 * MILLION;
      while (i < ITERATION_LIMIT){
	svaria_SampleCorr (gen, res, 1, size, 0, 1);
	DetectIteration (res->pVal2[gofw_Mean], &size, &i);
      }
      bbattery_pVal[++j] = res->pVal2[gofw_Mean];
      TestNumber[j] = j2;
      strcpy (bbattery_TestNames[j], "SampleCorr");
   }
   {
      smarsa_Res2 *res2;
      sres_Chi2 *res;
      res = sres_CreateChi2 ();
      ++j2;
      i = 0;
      size = 2 * MILLION;
      while (i < ITERATION_LIMIT){
         svaria_WeightDistrib (gen, res, 1, size, 0, 256, 0.0, 0.125);
	 DetectIteration (res->pVal2[gofw_Mean], &size, &i);
      }
      bbattery_pVal[++j] = res->pVal2[gofw_Mean];
      TestNumber[j] = j2;
      strcpy (bbattery_TestNames[j], "WeightDistrib, r = 0");

      ++j2;
      i = 0;
      size = 2 * MILLION;
      while (i < ITERATION_LIMIT){
         svaria_WeightDistrib (gen, res, 1, size, 8, 256, 0.0, 0.125);
	 DetectIteration (res->pVal2[gofw_Mean], &size, &i);
      }
      bbattery_pVal[++j] = res->pVal2[gofw_Mean];
      TestNumber[j] = j2;
      strcpy (bbattery_TestNames[j], "WeightDistrib, r = 8");

      ++j2;
      i = 0;
      size = 2 * MILLION;
      while (i < ITERATION_LIMIT){
         svaria_WeightDistrib (gen, res, 1, size, 16, 256, 0.0, 0.125);
	 DetectIteration (res->pVal2[gofw_Mean], &size, &i);
      }
      bbattery_pVal[++j] = res->pVal2[gofw_Mean];
      TestNumber[j] = j2;
      strcpy (bbattery_TestNames[j], "WeightDistrib, r = 16");

      ++j2;
      i = 0;
      size = 2 * MILLION;
      while (i < ITERATION_LIMIT){
         svaria_WeightDistrib (gen, res, 1, size, 24, 256, 0.0, 0.125);
	 DetectIteration (res->pVal2[gofw_Mean], &size, &i);
      }
      bbattery_pVal[++j] = res->pVal2[gofw_Mean];
      TestNumber[j] = j2;
      strcpy (bbattery_TestNames[j], "WeightDistrib, r = 24");

      ++j2;
      i = 0;
      size = 2 * MILLION;
      while (i < ITERATION_LIMIT){
         svaria_SumCollector (gen, res, 1, size, 0, 10.0);
	 DetectIteration (res->pVal2[gofw_Mean], &size, &i);
      }
      bbattery_pVal[++j] = res->pVal2[gofw_Mean];
      TestNumber[j] = j2;
      strcpy (bbattery_TestNames[j], "SumCollector");

      ++j2;
      i = 0;
      size = MILLION;
      while (i < ITERATION_LIMIT){
         smarsa_MatrixRank (gen, res, 1, size, r, s, 2 * s, 2 * s);
	 DetectIteration (res->pVal2[gofw_Mean], &size, &i);
      }
      bbattery_pVal[++j] = res->pVal2[gofw_Mean];
      TestNumber[j] = j2;
      strcpy (bbattery_TestNames[j], "MatrixRank, 60 x 60");

      ++j2;
      i = 0;
      size = MILLION;
      while (i < ITERATION_LIMIT){
         smarsa_MatrixRank (gen, res, 1, size, 20, 10, 2 * s, 2 * s);
	 DetectIteration (res->pVal2[gofw_Mean], &size, &i);
      }
      bbattery_pVal[++j] = res->pVal2[gofw_Mean];
      TestNumber[j] = j2;
      strcpy (bbattery_TestNames[j], "MatrixRank, 60 x 60");

      ++j2;
      i = 0;
      size = 50 * THOUSAND;
      while (i < ITERATION_LIMIT){
         smarsa_MatrixRank (gen, res, 1, size, r, s, 10 * s, 10 * s);
	 DetectIteration (res->pVal2[gofw_Mean], &size, &i);
      }
      bbattery_pVal[++j] = res->pVal2[gofw_Mean];
      TestNumber[j] = j2;
      strcpy (bbattery_TestNames[j], "MatrixRank, 300 x 300");

      ++j2;
      i=0;
      size = 50 * THOUSAND;
      while (i < ITERATION_LIMIT){
         smarsa_MatrixRank (gen, res, 1, size, 20, 10, 10 * s, 10 * s);
	 DetectIteration (res->pVal2[gofw_Mean], &size, &i);
      }
      bbattery_pVal[++j] = res->pVal2[gofw_Mean];
      TestNumber[j] = j2;
      strcpy (bbattery_TestNames[j], "MatrixRank, 300 x 300");

      ++j2;
      i = 0;
      size = 2 * THOUSAND;
      while (i < ITERATION_LIMIT){
         smarsa_MatrixRank (gen, res, 1, size, r, s, 40 * s, 40 * s);
	 DetectIteration (res->pVal2[gofw_Mean], &size, &i);
      }
      bbattery_pVal[++j] = res->pVal2[gofw_Mean];
      TestNumber[j] = j2;
      strcpy (bbattery_TestNames[j], "MatrixRank, 1200 x 1200");

      ++j2;
      i = 0;
      size = 2 * THOUSAND;
      while (i < ITERATION_LIMIT){
         smarsa_MatrixRank (gen, res, 1, size, 20, 10, 40 * s, 40 * s);
	 DetectIteration (res->pVal2[gofw_Mean], &size, &i);
      }
      bbattery_pVal[++j] = res->pVal2[gofw_Mean];
      TestNumber[j] = j2;
      strcpy (bbattery_TestNames[j], "MatrixRank, 1200 x 1200");

      ++j2;
      i = 0;
      size = 20 * MILLION;
      while (i < ITERATION_LIMIT){
         smarsa_Savir2 (gen, res, 1, size, 0, 1024*1024, 30);
	 DetectIteration (res->pVal2[gofw_Mean], &size, &i);
      }
      bbattery_pVal[++j] = res->pVal2[gofw_Mean];
      TestNumber[j] = j2;
      strcpy (bbattery_TestNames[j], "Savir2");

      sres_DeleteChi2 (res);

      res2 = smarsa_CreateRes2 ();
      ++j2;
      i = 0;
      size = 100 * MILLION;
      while (i < ITERATION_LIMIT){
	 smarsa_GCD (gen, res2, 1, size, 0, 30);
	 DetectIteration (res2->GCD->pVal2[gofw_Mean], &size, &i);
      }
      bbattery_pVal[++j] = res2->GCD->pVal2[gofw_Mean];
      TestNumber[j] = j2;
      strcpy (bbattery_TestNames[j], "GCD, r = 0");

      ++j2;
      i = 0;
      size = 40 * MILLION;
      while (i < ITERATION_LIMIT){
	 smarsa_GCD (gen, res2, 1, size, 10, 20);
	 DetectIteration (res2->GCD->pVal2[gofw_Mean], &size, &i);
      }
      bbattery_pVal[++j] = res2->GCD->pVal2[gofw_Mean];
      TestNumber[j] = j2;
      strcpy (bbattery_TestNames[j], "GCD, r = 10");

      smarsa_DeleteRes2 (res2);
   }
   {
      int rr, ss;
      long ll0, ll1;
      swalk_Res *res;
      res = swalk_CreateRes ();

      ++j2;

      /* RandomWalk1 with L0 = L1 = 90. */
      rr = 0;
      ss = 30;
      ll0 = 90;
      ll1 = 90;

      i = 0;
      size = 50 * MILLION;
      while (i < ITERATION_LIMIT){
	swalk_RandomWalk1 (gen, res, 1, size, rr, ss, ll0,  ll1);
	DetectIteration (res->H[0]->pVal2[gofw_Mean], &size, &i);
      }
      bbattery_pVal[++j] = res->H[0]->pVal2[gofw_Mean];
      TestNumber[j] = j2;
      strcpy (bbattery_TestNames[j], "RandomWalk1 H, (L = 90)");

      i = 0;
      size = 50 * MILLION;
      while (i < ITERATION_LIMIT){
	swalk_RandomWalk1 (gen, res, 1, size, rr, ss, ll0,  ll1);
	DetectIteration (res->M[0]->pVal2[gofw_Mean], &size, &i);
      }
      bbattery_pVal[++j] = res->M[0]->pVal2[gofw_Mean];
      TestNumber[j] = j2;
      strcpy (bbattery_TestNames[j], "RandomWalk1 M, (L = 90)");

      i = 0;
      size = 50 * MILLION;
      while (i < ITERATION_LIMIT){
	swalk_RandomWalk1 (gen, res, 1, size, rr, ss, ll0,  ll1);
	DetectIteration (res->J[0]->pVal2[gofw_Mean], &size, &i);
      }
      bbattery_pVal[++j] = res->J[0]->pVal2[gofw_Mean];
      TestNumber[j] = j2;
      strcpy (bbattery_TestNames[j], "RandomWalk1 J, (L = 90)");

      i = 0;
      size = 50 * MILLION;
      while (i < ITERATION_LIMIT){
	swalk_RandomWalk1 (gen, res, 1, size, rr, ss, ll0,  ll1);
	DetectIteration (res->R[0]->pVal2[gofw_Mean], &size, &i);
      }
      bbattery_pVal[++j] = res->R[0]->pVal2[gofw_Mean];
      TestNumber[j] = j2;
      strcpy (bbattery_TestNames[j], "RandomWalk1 R, (L = 90)");

      i = 0;
      size = 50 * MILLION;
      while (i < ITERATION_LIMIT){
	swalk_RandomWalk1 (gen, res, 1, size, rr, ss, ll0,  ll1);
	DetectIteration (res->C[0]->pVal2[gofw_Mean], &size, &i);
      }
      bbattery_pVal[++j] = res->C[0]->pVal2[gofw_Mean];
      TestNumber[j] = j2;
      strcpy (bbattery_TestNames[j], "RandomWalk1 C, (L = 90)");

      ++j2;
      /* RandomWalk1 with L0 = L1 = 90. */
      rr = 20;
      ss = 10;
      ll0 = 90;
      ll1 = 90;

      i = 0;
      size = 10 * MILLION;
      while (i < ITERATION_LIMIT){
	swalk_RandomWalk1 (gen, res, 1, size, rr, ss, ll0,  ll1);
	DetectIteration (res->H[0]->pVal2[gofw_Mean], &size, &i);
      }
      bbattery_pVal[++j] = res->H[0]->pVal2[gofw_Mean];
      TestNumber[j] = j2;
      strcpy (bbattery_TestNames[j], "RandomWalk1 H, (L = 90)");

      i = 0;
      size = 10 * MILLION;
      while (i < ITERATION_LIMIT){
	swalk_RandomWalk1 (gen, res, 1, size, rr, ss, ll0,  ll1);
	DetectIteration (res->M[0]->pVal2[gofw_Mean], &size, &i);
      }
      bbattery_pVal[++j] = res->M[0]->pVal2[gofw_Mean];
      TestNumber[j] = j2;
      strcpy (bbattery_TestNames[j], "RandomWalk1 M, (L = 90)");

      i = 0;
      size = 10 * MILLION;
      while (i < ITERATION_LIMIT){
	swalk_RandomWalk1 (gen, res, 1, size, rr, ss, ll0,  ll1);
	DetectIteration (res->J[0]->pVal2[gofw_Mean], &size, &i);
      }
      bbattery_pVal[++j] = res->J[0]->pVal2[gofw_Mean];
      TestNumber[j] = j2;
      strcpy (bbattery_TestNames[j], "RandomWalk1 J, (L = 90)");

      i = 0;
      size = 10 * MILLION;
      while (i < ITERATION_LIMIT){
	swalk_RandomWalk1 (gen, res, 1, size, rr, ss, ll0,  ll1);
	DetectIteration (res->R[0]->pVal2[gofw_Mean], &size, &i);
      }
      bbattery_pVal[++j] = res->R[0]->pVal2[gofw_Mean];
      TestNumber[j] = j2;
      strcpy (bbattery_TestNames[j], "RandomWalk1 R, (L = 90)");

      i = 0;
      size = 10 * MILLION;
      while (i < ITERATION_LIMIT){
	swalk_RandomWalk1 (gen, res, 1, size, rr, ss, ll0,  ll1);
	DetectIteration (res->C[0]->pVal2[gofw_Mean], &size, &i);
      }
      bbattery_pVal[++j] = res->C[0]->pVal2[gofw_Mean];
      TestNumber[j] = j2;
      strcpy (bbattery_TestNames[j], "RandomWalk1 C, (L = 90)");


      ++j2;
      /* RandomWalk1 with L0 = L1 = 1000. */
      rr = 0;
      ss = 30;
      ll0 = 1000;
      ll1 = 1000;

      i = 0;
      size = 5 * MILLION;
      while (i < ITERATION_LIMIT){
	swalk_RandomWalk1 (gen, res, 1, size, rr, ss, ll0,  ll1);
	DetectIteration (res->H[0]->pVal2[gofw_Mean], &size, &i);
      }
      bbattery_pVal[++j] = res->H[0]->pVal2[gofw_Mean];
      TestNumber[j] = j2;
      strcpy (bbattery_TestNames[j], "RandomWalk1 H, (L = 1000)");

      i = 0;
      size = 5 * MILLION;
      while (i < ITERATION_LIMIT){
	swalk_RandomWalk1 (gen, res, 1, size, rr, ss, ll0,  ll1);
	DetectIteration (res->M[0]->pVal2[gofw_Mean], &size, &i);
      }
      bbattery_pVal[++j] = res->M[0]->pVal2[gofw_Mean];
      TestNumber[j] = j2;
      strcpy (bbattery_TestNames[j], "RandomWalk1 M, (L = 1000)");

      i = 0;
      size = 5 * MILLION;
      while (i < ITERATION_LIMIT){
	swalk_RandomWalk1 (gen, res, 1, size, rr, ss, ll0,  ll1);
	DetectIteration (res->J[0]->pVal2[gofw_Mean], &size, &i);
      }
      bbattery_pVal[++j] = res->J[0]->pVal2[gofw_Mean];
      TestNumber[j] = j2;
      strcpy (bbattery_TestNames[j], "RandomWalk1 J, (L = 1000)");

      i = 0;
      size = 5 * MILLION;
      while (i < ITERATION_LIMIT){
	swalk_RandomWalk1 (gen, res, 1, size, rr, ss, ll0,  ll1);
	DetectIteration (res->R[0]->pVal2[gofw_Mean], &size, &i);
      }
      bbattery_pVal[++j] = res->R[0]->pVal2[gofw_Mean];
      TestNumber[j] = j2;
      strcpy (bbattery_TestNames[j], "RandomWalk1 R, (L = 1000)");

      i = 0;
      size = 5 * MILLION;
      while (i < ITERATION_LIMIT){
	swalk_RandomWalk1 (gen, res, 1, size, rr, ss, ll0,  ll1);
	DetectIteration (res->C[0]->pVal2[gofw_Mean], &size, &i);
      }
      bbattery_pVal[++j] = res->C[0]->pVal2[gofw_Mean];
      TestNumber[j] = j2;
      strcpy (bbattery_TestNames[j], "RandomWalk1 C, (L = 1000)");

      ++j2;
      /* RandomWalk1 with L0 = L1 = 1000. */
      rr = 20;
      ss = 10;
      ll0 = 1000;
      ll1 = 1000;

      i = 0;
      size = MILLION;
      while (i < ITERATION_LIMIT){
	swalk_RandomWalk1 (gen, res, 1, size, rr, ss, ll0,  ll1);
	DetectIteration (res->H[0]->pVal2[gofw_Mean], &size, &i);
      }
      bbattery_pVal[++j] = res->H[0]->pVal2[gofw_Mean];
      TestNumber[j] = j2;
      strcpy (bbattery_TestNames[j], "RandomWalk1 H, (L = 1000)");

      i = 0;
      size = MILLION;
      while (i < ITERATION_LIMIT){
	swalk_RandomWalk1 (gen, res, 1, size, rr, ss, ll0,  ll1);
	DetectIteration (res->M[0]->pVal2[gofw_Mean], &size, &i);
      }
      bbattery_pVal[++j] = res->M[0]->pVal2[gofw_Mean];
      TestNumber[j] = j2;
      strcpy (bbattery_TestNames[j], "RandomWalk1 M, (L = 1000)");

      i = 0;
      size = MILLION;
      while (i < ITERATION_LIMIT){
	swalk_RandomWalk1 (gen, res, 1, size, rr, ss, ll0,  ll1);
	DetectIteration (res->J[0]->pVal2[gofw_Mean], &size, &i);
      }
      bbattery_pVal[++j] = res->J[0]->pVal2[gofw_Mean];
      TestNumber[j] = j2;
      strcpy (bbattery_TestNames[j], "RandomWalk1 J, (L = 1000)");

      i = 0;
      size = MILLION;
      while (i < ITERATION_LIMIT){
	swalk_RandomWalk1 (gen, res, 1, size, rr, ss, ll0,  ll1);
	DetectIteration (res->R[0]->pVal2[gofw_Mean], &size, &i);
      }
      bbattery_pVal[++j] = res->R[0]->pVal2[gofw_Mean];
      TestNumber[j] = j2;
      strcpy (bbattery_TestNames[j], "RandomWalk1 R, (L = 1000)");

      i = 0;
      size = MILLION;
      while (i < ITERATION_LIMIT){
	swalk_RandomWalk1 (gen, res, 1, size, rr, ss, ll0,  ll1);
	DetectIteration (res->C[0]->pVal2[gofw_Mean], &size, &i);
      }
      bbattery_pVal[++j] = res->C[0]->pVal2[gofw_Mean];
      TestNumber[j] = j2;
      strcpy (bbattery_TestNames[j], "RandomWalk1 C, (L = 1000)");

      ++j2;
      /* RandomWalk1 with L0 = L1 = 10000. */
      rr = 0;
      ss = 30;
      ll0 = 10000;
      ll1 = 10000;

      i = 0;
      size = MILLION / 2;
      while (i < ITERATION_LIMIT){
	swalk_RandomWalk1 (gen, res, 1, size, rr, ss, ll0,  ll1);
	DetectIteration (res->H[0]->pVal2[gofw_Mean], &size, &i);
      }
      bbattery_pVal[++j] = res->H[0]->pVal2[gofw_Mean];
      TestNumber[j] = j2;
      strcpy (bbattery_TestNames[j], "RandomWalk1 H, (L = 10000)");

      i = 0;
      size = MILLION / 2;
      while (i < ITERATION_LIMIT){
	swalk_RandomWalk1 (gen, res, 1, size, rr, ss, ll0,  ll1);
	DetectIteration (res->M[0]->pVal2[gofw_Mean], &size, &i);
      }
      bbattery_pVal[++j] = res->M[0]->pVal2[gofw_Mean];
      TestNumber[j] = j2;
      strcpy (bbattery_TestNames[j], "RandomWalk1 M, (L = 10000)");

      i = 0;
      size = MILLION / 2;
      while (i < ITERATION_LIMIT){
	swalk_RandomWalk1 (gen, res, 1, size, rr, ss, ll0,  ll1);
	DetectIteration (res->J[0]->pVal2[gofw_Mean], &size, &i);
      }
      bbattery_pVal[++j] = res->J[0]->pVal2[gofw_Mean];
      TestNumber[j] = j2;
      strcpy (bbattery_TestNames[j], "RandomWalk1 J, (L = 10000)");

      i = 0;
      size = MILLION / 2;
      while (i < ITERATION_LIMIT){
	swalk_RandomWalk1 (gen, res, 1, size, rr, ss, ll0,  ll1);
	DetectIteration (res->R[0]->pVal2[gofw_Mean], &size, &i);
      }
      bbattery_pVal[++j] = res->R[0]->pVal2[gofw_Mean];
      TestNumber[j] = j2;
      strcpy (bbattery_TestNames[j], "RandomWalk1 R, (L = 10000)");

      i = 0;
      size = MILLION / 2;
      while (i < ITERATION_LIMIT){
	swalk_RandomWalk1 (gen, res, 1, size, rr, ss, ll0,  ll1);
	DetectIteration (res->C[0]->pVal2[gofw_Mean], &size, &i);
      }
      bbattery_pVal[++j] = res->C[0]->pVal2[gofw_Mean];
      TestNumber[j] = j2;
      strcpy (bbattery_TestNames[j], "RandomWalk1 C, (L = 10000)");

      ++j2;
      /* RandomWalk1 with L0 = L1 = 10000. */
      rr = 20;
      ss = 10;
      ll0 = 10000;
      ll1 = 10000;

      i = 0;
      size = MILLION / 10;
      while (i < ITERATION_LIMIT){
	swalk_RandomWalk1 (gen, res, 1, size, rr, ss, ll0,  ll1);
	DetectIteration (res->H[0]->pVal2[gofw_Mean], &size, &i);
      }
      bbattery_pVal[++j] = res->H[0]->pVal2[gofw_Mean];
      TestNumber[j] = j2;
      strcpy (bbattery_TestNames[j], "RandomWalk1 H, (L = 10000)");

      i = 0;
      size = MILLION / 10;
      while (i < ITERATION_LIMIT){
	swalk_RandomWalk1 (gen, res, 1, size, rr, ss, ll0,  ll1);
	DetectIteration (res->M[0]->pVal2[gofw_Mean], &size, &i);
      }
      bbattery_pVal[++j] = res->M[0]->pVal2[gofw_Mean];
      TestNumber[j] = j2;
      strcpy (bbattery_TestNames[j], "RandomWalk1 M, (L = 10000)");

      i = 0;
      size = MILLION / 10;
      while (i < ITERATION_LIMIT){
	swalk_RandomWalk1 (gen, res, 1, size, rr, ss, ll0,  ll1);
	DetectIteration (res->J[0]->pVal2[gofw_Mean], &size, &i);
      }
      bbattery_pVal[++j] = res->J[0]->pVal2[gofw_Mean];
      TestNumber[j] = j2;
      strcpy (bbattery_TestNames[j], "RandomWalk1 J, (L = 10000)");

      i = 0;
      size = MILLION / 10;
      while (i < ITERATION_LIMIT){
	swalk_RandomWalk1 (gen, res, 1, size, rr, ss, ll0,  ll1);
	DetectIteration (res->R[0]->pVal2[gofw_Mean], &size, &i);
      }
      bbattery_pVal[++j] = res->R[0]->pVal2[gofw_Mean];
      TestNumber[j] = j2;
      strcpy (bbattery_TestNames[j], "RandomWalk1 R, (L = 10000)");

      i = 0;
      size = MILLION / 10;
      while (i < ITERATION_LIMIT){
	swalk_RandomWalk1 (gen, res, 1, size, rr, ss, ll0,  ll1);
	DetectIteration (res->C[0]->pVal2[gofw_Mean], &size, &i);
      }
      bbattery_pVal[++j] = res->C[0]->pVal2[gofw_Mean];
      TestNumber[j] = j2;
      strcpy (bbattery_TestNames[j], "RandomWalk1 C, (L = 10000)");

      swalk_DeleteRes (res);
   }
   {
      scomp_Res *res;
      res = scomp_CreateRes ();
      ++j2;
      i = 0;
      size = 120 * THOUSAND;
      while (i < ITERATION_LIMIT){
         scomp_LinearComp (gen, res, 1, size, 0, 1);
	 DetectIteration (res->JumpNum->pVal2[gofw_Mean], &size, &i);
      }
      bbattery_pVal[++j] = res->JumpNum->pVal2[gofw_Mean];
      TestNumber[j] = j2;
      strcpy (bbattery_TestNames[j], "LinearComp, r = 0");
      bbattery_pVal[++j] = res->JumpSize->pVal2[gofw_Mean];
      TestNumber[j] = j2;
      strcpy (bbattery_TestNames[j], "LinearComp, r = 0");

      ++j2;
      i = 0;
      size = 120 * THOUSAND;
      while (i < ITERATION_LIMIT){
	 scomp_LinearComp (gen, res, 1, size, 29, 1);
	 DetectIteration (res->JumpNum->pVal2[gofw_Mean], &size, &i);
      }
      bbattery_pVal[++j] = res->JumpNum->pVal2[gofw_Mean];
      TestNumber[j] = j2;
      strcpy (bbattery_TestNames[j], "LinearComp, r = 29");
      bbattery_pVal[++j] = res->JumpSize->pVal2[gofw_Mean];
      TestNumber[j] = j2;
      strcpy (bbattery_TestNames[j], "LinearComp, r = 29");

      scomp_DeleteRes (res);
   }
   {
      sstring_Res2 *res;
      res = sstring_CreateRes2 ();
      ++j2;
      i = 0;
      size = 1000;
      while (i < ITERATION_LIMIT){
         sstring_LongestHeadRun (gen, res, 1, size, r, s, 20 + 10 * MILLION);
	 DetectIteration (res->Chi->pVal2[gofw_Mean], &size, &i);
      }
      bbattery_pVal[++j] = res->Chi->pVal2[gofw_Mean];
      TestNumber[j] = j2;
      strcpy (bbattery_TestNames[j], "LongestHeadRun, r = 0");
      bbattery_pVal[++j] = res->Disc->pVal2;
      TestNumber[j] = j2;
      strcpy (bbattery_TestNames[j], "LongestHeadRun, r = 0");

      ++j2;
      i = 0;
      size = 300;
      while (i < ITERATION_LIMIT){
         sstring_LongestHeadRun (gen, res, 1, size, 20, 10, 20 + 10 * MILLION);
	 DetectIteration (res->Chi->pVal2[gofw_Mean], &size, &i);
      }
      bbattery_pVal[++j] = res->Chi->pVal2[gofw_Mean];
      TestNumber[j] = j2;
      strcpy (bbattery_TestNames[j], "LongestHeadRun, r = 20");
      bbattery_pVal[++j] = res->Disc->pVal2;
      TestNumber[j] = j2;
      strcpy (bbattery_TestNames[j], "LongestHeadRun, r = 20");

      sstring_DeleteRes2 (res);
   }
   {
      sres_Chi2 *res;
      res = sres_CreateChi2 ();
      ++j2;
      i = 0;
      size = 300 * MILLION;
      while (i < ITERATION_LIMIT){
	 sstring_PeriodsInStrings (gen, res, 1, size, r, s);
	 DetectIteration (res->pVal2[gofw_Mean], &size, &i);
      }
      bbattery_pVal[++j] = res->pVal2[gofw_Mean];
      TestNumber[j] = j2;
      strcpy (bbattery_TestNames[j], "PeriodsInStrings, r = 0");

      ++j2;
      i = 0;
      size = 300 * MILLION;
      while (i < ITERATION_LIMIT){
	 sstring_PeriodsInStrings (gen, res, 1, size, 15, 15);
	 DetectIteration (res->pVal2[gofw_Mean], &size, &i);
      }
      bbattery_pVal[++j] = res->pVal2[gofw_Mean];
      TestNumber[j] = j2;
      strcpy (bbattery_TestNames[j], "PeriodsInStrings, r = 15");

      sres_DeleteChi2 (res);
   }
   {
      sres_Basic *res;
      res = sres_CreateBasic ();
      ++j2;
      i = 0;
      size = 100 * MILLION;
      while (i < ITERATION_LIMIT){
         sstring_HammingWeight2 (gen, res, 1, size, r, s, MILLION);
	 DetectIteration (res->pVal2[gofw_Mean], &size, &i);
      }
      bbattery_pVal[++j] = res->pVal2[gofw_Mean];
      TestNumber[j] = j2;
      strcpy (bbattery_TestNames[j], "HammingWeight2, r = 0");

      ++j2;
      i = 0;
      size = 100 * MILLION;
      while (i < ITERATION_LIMIT){
         sstring_HammingWeight2 (gen, res, 1, size, 20, 10, MILLION);
	 DetectIteration (res->pVal2[gofw_Mean], &size, &i);
      }
      bbattery_pVal[++j] = res->pVal2[gofw_Mean];
      TestNumber[j] = j2;
      strcpy (bbattery_TestNames[j], "HammingWeight2, r = 20");

      sres_DeleteBasic (res);
   }
   {
      sstring_Res *res;
      res = sstring_CreateRes ();
      /* sstring_HammingCorr will probably be removed: less sensitive than
         svaria_HammingIndep */

      ++j2;
      i = 0;
      size = 500 * MILLION;
      while (i < ITERATION_LIMIT){
         sstring_HammingCorr (gen, res, 1, size, r, s, s);
	 DetectIteration (res->Bas->pVal2[gofw_Mean], &size, &i);
      }
      bbattery_pVal[++j] = res->Bas->pVal2[gofw_Mean];
      TestNumber[j] = j2;
      strcpy (bbattery_TestNames[j], "HammingCorr, L = 30");

      ++j2;
      i = 0;
      size = 50 * MILLION;
      while (i < ITERATION_LIMIT){
         sstring_HammingCorr (gen, res, 1, size, r, s, 10 * s);
	 DetectIteration (res->Bas->pVal2[gofw_Mean], &size, &i);
      }
      bbattery_pVal[++j] = res->Bas->pVal2[gofw_Mean];
      TestNumber[j] = j2;
      strcpy (bbattery_TestNames[j], "HammingCorr, L = 300");

      ++j2;
      i = 0;
      size = 10 * MILLION;
      while (i < ITERATION_LIMIT){
         sstring_HammingCorr (gen, res, 1, size, r, s, 40 * s);
	 DetectIteration (res->Bas->pVal2[gofw_Mean], &size, &i);
      }
      bbattery_pVal[++j] = res->Bas->pVal2[gofw_Mean];
      TestNumber[j] = j2;
      strcpy (bbattery_TestNames[j], "HammingCorr, L = 1200");

      ++j2;
      i = 0;
      size = 300 * MILLION;
      while (i < ITERATION_LIMIT){
         sstring_HammingIndep (gen, res, 1, size, r, s, s, 0);
	 DetectIteration (res->Bas->pVal2[gofw_Mean], &size, &i);
      }
      bbattery_pVal[++j] = res->Bas->pVal2[gofw_Mean];
      TestNumber[j] = j2;
      strcpy (bbattery_TestNames[j], "HammingIndep, L = 30");

      ++j2;
      i = 0;
      size = 100 * MILLION;
      while (i < ITERATION_LIMIT){
         sstring_HammingIndep (gen, res, 1, size, 20, 10, s, 0);
	 DetectIteration (res->Bas->pVal2[gofw_Mean], &size, &i);
      }
      bbattery_pVal[++j] = res->Bas->pVal2[gofw_Mean];
      TestNumber[j] = j2;
      strcpy (bbattery_TestNames[j], "HammingIndep, L = 30");

      ++j2;
      i = 0;
      size = 30 * MILLION;
      while (i < ITERATION_LIMIT){
         sstring_HammingIndep (gen, res, 1, size, r, s, 10 * s, 0);
	 DetectIteration (res->Bas->pVal2[gofw_Mean], &size, &i);
      }
      bbattery_pVal[++j] = res->Bas->pVal2[gofw_Mean];
      TestNumber[j] = j2;
      strcpy (bbattery_TestNames[j], "HammingIndep, L = 300");

      ++j2;
      i = 0;
      size = 10 * MILLION;
      while (i < ITERATION_LIMIT){
         sstring_HammingIndep (gen, res, 1, size, 20, 10, 10 * s, 0);
	 DetectIteration (res->Bas->pVal2[gofw_Mean], &size, &i);
      }
      bbattery_pVal[++j] = res->Bas->pVal2[gofw_Mean];
      TestNumber[j] = j2;
      strcpy (bbattery_TestNames[j], "HammingIndep, L = 300");

      ++j2;
      i = 0;
      size = 10 * MILLION;
      while (i < ITERATION_LIMIT){
         sstring_HammingIndep (gen, res, 1, size, r, s, 40 * s, 0);
	 DetectIteration (res->Bas->pVal2[gofw_Mean], &size, &i);
      }
      bbattery_pVal[++j] = res->Bas->pVal2[gofw_Mean];
      TestNumber[j] = j2;
      strcpy (bbattery_TestNames[j], "HammingIndep, L = 1200");

      ++j2;
      i = 0;
      size = MILLION;
      while (i < ITERATION_LIMIT){
         sstring_HammingIndep (gen, res, 1, size, 20, 10, 40 * s, 0);
	 DetectIteration (res->Bas->pVal2[gofw_Mean], &size, &i);
      }
      bbattery_pVal[++j] = res->Bas->pVal2[gofw_Mean];
      TestNumber[j] = j2;
      strcpy (bbattery_TestNames[j], "HammingIndep, L = 1200");

      sstring_DeleteRes (res);
   }
   {
      sstring_Res3 *res;
      res = sstring_CreateRes3 ();
      ++j2;
      i=0;
      size = BILLION;
      while (i < ITERATION_LIMIT){
         sstring_Run (gen, res, 1, size, r, s);
	 DetectIteration (res->NRuns->pVal2[gofw_Mean], &size, &i);
      }
      bbattery_pVal[++j] = res->NRuns->pVal2[gofw_Mean];
      TestNumber[j] = j2;
      strcpy (bbattery_TestNames[j], "Run of bits, r = 0");
      bbattery_pVal[++j] = res->NBits->pVal2[gofw_Mean];
      TestNumber[j] = j2;
      strcpy (bbattery_TestNames[j], "Run of bits, r = 0");

      ++j2;
      i = 0;
      size = BILLION;
      while (i < ITERATION_LIMIT){
         sstring_Run (gen, res, 1, size, 20, 10);
	 DetectIteration (res->NRuns->pVal2[gofw_Mean], &size, &i);
      }
      bbattery_pVal[++j] = res->NRuns->pVal2[gofw_Mean];
      TestNumber[j] = j2;
      strcpy (bbattery_TestNames[j], "Run of bits, r = 20");
      bbattery_pVal[++j] = res->NBits->pVal2[gofw_Mean];
      TestNumber[j] = j2;
      strcpy (bbattery_TestNames[j], "Run of bits, r = 20");

      sstring_DeleteRes3 (res);
   }
   {
      sres_Basic *res;
      res = sres_CreateBasic ();
      ++j2;
      i = 0;
      size = 30 + BILLION;
      while (i < ITERATION_LIMIT){
         sstring_AutoCor (gen, res, 1,  size, r, s, 1);
	 DetectIteration (res->pVal2[gofw_Mean], &size, &i);
      }
      bbattery_pVal[++j] = res->pVal2[gofw_Mean];
      TestNumber[j] = j2;
      strcpy (bbattery_TestNames[j], "AutoCor, d = 1");

      ++j2;
      i = 0;
      size = 1 + BILLION;
      while (i < ITERATION_LIMIT){
         sstring_AutoCor (gen, res, 1,  size, 20, 10, 1);
	 DetectIteration (res->pVal2[gofw_Mean], &size, &i);
      }
      bbattery_pVal[++j] = res->pVal2[gofw_Mean];
      TestNumber[j] = j2;
      strcpy (bbattery_TestNames[j], "AutoCor, d = 1");

      ++j2;
      i = 0;
      size = 31 + BILLION;
      while (i < ITERATION_LIMIT){
	 sstring_AutoCor (gen, res, 1,  size, r, s, s);
	 DetectIteration (res->pVal2[gofw_Mean], &size, &i);
      }
      bbattery_pVal[++j] = res->pVal2[gofw_Mean];
      TestNumber[j] = j2;
      strcpy (bbattery_TestNames[j], "AutoCor, d = 30");

      ++j2;
/*     util_Assert (j2 <= CRUSH_NUM, "Crush:   j2 > CRUSH_NUM");  */
      i = 0;
      size = 11 + BILLION;
      while (i < ITERATION_LIMIT){
	 sstring_AutoCor (gen, res, 1,  size, 20, 10, 10);
	 DetectIteration (res->pVal2[gofw_Mean], &size, &i);
      }
      bbattery_pVal[++j] = res->pVal2[gofw_Mean];
      TestNumber[j] = j2;
      strcpy (bbattery_TestNames[j], "AutoCor, d = 10");

      sres_DeleteBasic (res);
   }

   bbattery_NTests = ++j;
   GetName (gen, genName);
   WriteReport (genName, "Crush", bbattery_NTests,
      bbattery_pVal, Timer, FALSE, TRUE, 0.0);
   chrono_Delete (Timer);
}


/*=========================================================================*/

void bbattery_Adaptive_Crush (unif01_Gen * gen)
{
   int i;
   int Rep[NDIM + 1] = {0};
   for (i = 1; i <= ADAPTIVE_CRUSH_NUM; ++i)
      Rep[i] = 1;
   Adaptive_Crush (gen, Rep);
}


/*=========================================================================*/

void bbattery_Repeat_Adaptive_Crush (unif01_Gen * gen, int Rep[])
{
   Adaptive_Crush (gen, Rep);
}

/*=========================================================================*/
#if 0
static void WriteTime (time_t t0, time_t t1)
{
   int y1;
   double y = 0;

   y = difftime (t1, t0);
   /* printf (" Total time: %.2f sec\n\n", y); */
   printf (" Total time: ");
   y1 = y / 3600;
   printf ("%02d:", y1);
   y -= y1 * 3600.0;
   y1 = y / 60;
   printf ("%02d:", y1);
   y -= y1 * 60.0;
   printf ("%.2f\n\n", y);
}
#endif

void DetectIteration (double pvalue, long *size, int *i)
{
  if (pvalue < BETA || pvalue > 1.0-BETA)
    (*i) = BREAK;
  else if (pvalue < ALPHA || pvalue > 1.0-ALPHA){
    (*size) = 2 * (*size);
    (*i)++;
  }
  else
    (*i) = BREAK;
}
