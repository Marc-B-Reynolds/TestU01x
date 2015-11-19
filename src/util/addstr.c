/*************************************************************************\
* Copyright (c) 2002 Pierre L'Ecuyer, DIRO, Université de Montréal.
* SEE: LICENSE
\*************************************************************************/

#define _CRT_SECURE_NO_WARNINGS

#include "addstr.h"

#include <stdio.h>
#include <string.h>



#define LEN1 63



/* ========================== functions  ================================= */

void  addstr_Long (char *to, const char *add, long n)
{
   char str[LEN1 + 1];
   strcat (to, add);
   sprintf (str, "%1ld", n);
   strcat (to, str);
}

void  addstr_Uint (char *to, const char *add, unsigned int n)
{
   char str[LEN1 + 1];
   strcat (to, add);
   sprintf (str, "%1u", n);
   strcat (to, str);
}

void  addstr_Int (char *to, const char *add, int n)
{
   char str[LEN1 + 1];
   strcat (to, add);
   sprintf (str, "%1d", n);
   strcat (to, str);
}

void  addstr_Ulong (char *to, const char *add, unsigned long n)
{
   char str[LEN1 + 1];
   strcat (to, add);
   sprintf (str, "%1lu", n);
   strcat (to, str);
}

#ifdef USE_LONGLONG
void  addstr_LONG (char *to, const char *add, longlong n)
{
   char str[LEN1 + 1];
   strcat (to, add);
   sprintf (str, "%1" PRIdLEAST64, n);
   strcat (to, str);
}

void  addstr_ULONG (char *to, const char *add, ulonglong n)
{
   char str[LEN1 + 1];
   strcat (to, add);
   sprintf (str, "%1" PRIuLEAST64, n);
   strcat (to, str);
}
#endif

void  addstr_Char (char *to, const char *add, char c)
{
   char str[LEN1 + 1];
   strcat (to, add);
   sprintf (str, "%c", c);
   strcat (to, str);
}

void  addstr_Double (char *to, const char *add, double x)
{
   char str[LEN1 + 1];
   strcat (to, add);
   sprintf (str, "%.16G", x);
   strcat (to, str);
}

void  addstr_Bool (char *to, const char *add, int b)
{
   strcat (to, add);
   if (b)
      strcat (to, "TRUE");
   else
      strcat (to, "FALSE");
}

void  addstr_ArrayLong (char *to, const char *add, int high, long val[])
{
   int j;
   strcat (to, add);
   addstr_Long (to, "(", val[0]);
   for (j = 1; (j < high) && (j < 5); j++)
      addstr_Long (to, ", ", val[j]);
   if (high > 5)
      strcat (to, ", ... )");
   else
      strcat (to, ")");
}

void  addstr_ArrayUlong (char *to, const char *add, 
                         int high, unsigned long val[])
{
   int j;
   strcat (to, add);
   addstr_Ulong (to, "(", val[0]);
   for (j = 1; (j < high) && (j < 5); j++)
      addstr_Ulong (to, ", ", val[j]);
   if (high > 5)
      strcat (to, ", ... )");
   else
      strcat (to, ")");
}

void  addstr_ArrayUint (char *to, const char *add, int high, 
		        unsigned int val[])
{
   int j;
   strcat (to, add);
   addstr_Uint (to, "(", val[0]);
   for (j = 1; (j < high) && (j < 5); j++)
      addstr_Uint (to, ", ", val[j]);
   if (high > 5)
      strcat (to, ", ... )");
   else
      strcat (to, ")");
}

void  addstr_ArrayInt (char *to, const char *add, int high, int val[])
{
   int j;
   strcat (to, add);
   addstr_Int (to, "(", val[0]);
   for (j = 1; (j < high) && (j < 5); j++)
      addstr_Int (to, ", ", val[j]);
   if (high > 5)
      strcat (to, ", ... )");
   else
      strcat (to, ")");
}

void  addstr_ArrayDouble (char *to, const char *add, 
                          int high, double val[])
{
   int j;
   strcat (to, add);
   addstr_Double (to, "(", val[0]);
   for (j = 1; (j < high) && (j < 5); j++)
      addstr_Double (to, ", ", val[j]);
   if (high > 5)
      strcat (to, ", ... )");
   else
      strcat (to, ")");
}
