
// This will run until stopped

#if 0

#include "../../Stand-alone-junk/src/TestAndSearch/PrnsTestU01.c"

#else

#include <stdio.h>
#include <stdint.h>

#if defined(_MSC_VER)
#define inline _inline
#define I2F (1.0/((1.0*(1<<22))*(1.0*(1<<30))))
_inline uint32_t __builtin_ctz(uint32_t x) { unsigned long r; _BitScanForward(&r, (unsigned long)x); return (uint32_t)r; }
#else
#include <x86intrin.h>
#define I2F 0x1p-52f
#endif

#include "util.h"
#include "unif01.h"
#include "swrite.h"

#define NAME "bs"

// if defined the lower 32-bit results are used for
// integer test, otherwise the upper. TODO: change
// to alternating
//#define USE_LOWER_BITS

// the single 64-bit state data of the generator
uint64_t state;

// POT LCG - 
#define M0 0x369dea0f31a53f85L
#define A0 0x61C8864680b583EBL

// full avalanche bit-mixer (http://zimbry.blogspot.fr/2011/09/better-bit-mixing-improving-on.html)
#if 0
#define S0 32
#define S1 29
#define S2 32
#define M1 0x4cd6944c5cc20b6dL
#define M2 0xfc12c5b19d3259e9L
#else
#define S0 31
#define S1 27
#define S2 33
#define M1 0x7fb5d329728ea185L
#define M2 0x81dadef4bc2dd44dL
#endif


#if 1
inline uint64_t next()
{
  uint64_t x = state;

  state += 0x9e3779b97f4a7c15L;

  x ^= (x >> S0); x *= M1;
  x ^= (x >> S1);	x *= M2;
  x ^= (x >> S2);

  return x;
}
#else
// dep-chain = 7
inline uint64_t next()
{
  uint64_t v = state;

  // update the state for next result
  state = state*M0;

  // perform the mix
  v ^= (v >> S0); v *= M1;
  v ^= (v >> S1); v *= M2;
  v ^= (v >> S2);

  return v;
}
#endif

#define W0 0x3504f333   // 3*2309*128413 
#define W1 0xf1bbcdcb   // 7*349*1660097 
#define M  741103597    // 13*83*686843

//--------

static inline uint32_t time_rng_hash(uint32_t t)
{
	return W1*t;
}

uint32_t rng_setup(uint32_t entity_id, uint32_t thash)
{
	return W0*entity_id ^ thash;
}

static uint32_t next_u32(void* p, void* s)
{
  uint64_t r = next();

#if defined(USE_LOWER_BITS)  
  return (uint32_t)r;
#else
  return (uint32_t)(r >> 32);
#endif
}

static double next_f64(void* p, void* s)
{
  uint64_t r = next();
  // map to [0,1)
  return (r >> 12) * I2F;
}

static void print_state(void* s)
{
  // screw big-endian
  printf ("  S = 0x%08x%08x\n", (uint32_t)(state>>32), (uint32_t)state);
}

// 
unif01_Gen* createGenerator()
{
  unif01_Gen* gen = util_Malloc(sizeof(unif01_Gen));

  gen->state   = 0;
  gen->param   = 0;
  gen->name    = NAME;
  gen->GetU01  = (void*)&next_f64;
  gen->GetBits = (void*)&next_u32;
  gen->Write   = &print_state;
  
  return gen;
}

void deleteGenerator(unif01_Gen* gen)
{
  if (gen != NULL) {
    //util_Free(gen->param);
    //util_Free(gen->state);
    util_Free(gen);
  }
}


#include "bbattery.h"

int main (void) 
{
  unif01_Gen* gen = createGenerator();
  uint64_t    s =  __rdtsc();
  uint32_t    c = -1;
  uint32_t    t;

  swrite_Basic = FALSE; // only print summary

  uint64_t z = M0;
  
  printf ("0x%08x%08x\n", (uint32_t)(z>>32), (uint32_t)z);

    do {
    state = s;
    
    // screw big-endian
    printf ("run %d -- state = 0x%08x%08x\n", -(int32_t)c,(uint32_t)(state>>32), (uint32_t)state);

    if (state != 0) {
      //bbattery_SmallCrush(gen);
      bbattery_Crush(gen);
    }
    else {
      printf ("...skipping illegal state\n");
    }
    
    // could improve the domain of the search
    t = __builtin_ctz(c)*2;
    s ^= 0x8000000000000000UL >> t;
    c -= 1;
  } while(1);
  
  deleteGenerator(gen);
  
  return 0;
}
#endif