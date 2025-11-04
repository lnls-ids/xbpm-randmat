/* Extracted from
 * https://en.wikipedia.org/wiki/Permuted_congruential_generator
 * on 2023-11-21.
 */

#include <stdint.h>

#ifndef PCG_RANDOM
#define PCG_RANDOM 1

// #ifndef RMAX
#define RMAX  9223372036854775808.0    /* 2**63 */
//#define RMAX  2147483648.0    /* 2**31 */

static uint64_t       state      = 0x4d595df4d0f33173;   // Or something seed-dependent
static uint64_t const multiplier = 6364136223846793005u;
static uint64_t const increment  = 1442695040888963407u; // Or an arbitrary odd constant

static uint32_t rotr32 (uint32_t x, unsigned r)
{
  return x >> r | x << (-r & 31);
}

uint32_t pcg32 (void)
{
  uint64_t x = state;
  unsigned count = (unsigned)(x >> 59);        /* 59 = 64 - 5      */
  
  state = x * multiplier + increment;
  x ^= x >> 18;                                /* 18 = (64 - 27)/2 */
  return rotr32((uint32_t)(x >> 27), count);   /* 27 = 32 - 5      */
}

void pcg32_init (uint64_t seed)
{
  state = seed + increment;
  (void)pcg32();
}

long double pcg_double (void)
{
  /* Assemble a 64-bit random integer. */
  uint64_t x = ((uint64_t) pcg32() << 31) | pcg32();

  //printf(" DEBUG x: %li \n", x);
  
  /* Cast to double and normalize it. */
  return x / RMAX;
  //return (double) x / RMAX;
  //return ((double)x * 100.0) / RMAX;
}

/* Generic alias for the standard PCG generator. */
//double (*pcg_rand)() = &pcg_double;

#endif
