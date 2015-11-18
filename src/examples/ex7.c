
#include "unif01.h"
#include "bbattery.h"
#include "sres.h"

unsigned int xorshift32(void);
double xorshift64(void);
void setversion32(int x);
void setversion64(int x);

extern char* name32;
extern char* name64;

uint16_t good_32[] = {496};

int main (void) 
{
   unif01_Gen *gen;
   int i;

   //swrite_Host = FALSE;

#if 0
   for(i=0; i<2200; i++) {
     printf("========= TEST %d\n", i);
     
     setversion64(0);
     gen = unif01_CreateExternGen01 (name64, xorshift64);
     bbattery_SmallCrush(gen);
     //bbattery_Crush(gen);
     //bbattery_BigCrush (gen);
     unif01_DeleteExternGen01 (gen);
     printf("======== END OF TEST %d\n", i);
   }
#endif

#if 0
   for(i=0; i<648; i++) {
     printf("========= TEST %d\n", i);
     setversion32(i);
     gen = unif01_CreateExternGenBits(name32, xorshift32);
     bbattery_SmallCrush (gen);
     //bbattery_Crush(gen);
     //bbattery_BigCrush (gen);
     unif01_DeleteExternGenBits(gen);
     printf("======== END OF TEST %d\n", i);
   }
#endif

#if  1
   // different seed change results (all failed 4 in original brute test)
   // total failures in paran is more than eps fail
   setversion32(427);  // failing 7(3)
   setversion32(72);   // failing 5(1)
   setversion32(328);  // failing 5(1)
   setversion32(496);  // failing 6
   setversion32(389);  // failing 6
   setversion32(496);  // failing 6(3)

   setversion32(388);  // failing 4(1)

   //still failing 4 with 1 @ 1-eps
   setversion32(116);
   setversion32(282);
   setversion32(386);
   setversion32(388);
   setversion32(534);

   // still failing 4 with default seed
   setversion32(96);
   setversion32(209);
   setversion32(212);
   setversion32(272);
   setversion32(280);
   setversion32(284);
   setversion32(308);
   setversion32(404);
   setversion32(410);
   setversion32(411);
   setversion32(414);

   gen = unif01_CreateExternGenBits(name32, xorshift32);
   //bbattery_SmallCrush (gen);
   bbattery_Crush(gen);
   //bbattery_BigCrush (gen);
#endif

   return 0;
}
