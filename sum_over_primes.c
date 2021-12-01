#include "sum_over_primes.h"
#include <math.h>

//
// Common internal routines, needed by the implementation
//

/*                                                                              */
/* Primes Sieve                                                                 */
/*                                                                              */

// we will index primes starting with 1, as maths do
// element zero will be set to 0
static inline uint generate_prime_list_up_to(u64 N, u64** out_prime) {
    if(N < 2) return 0;

    // we do quite simple Sieve of Eratosthenes:
    // we consider only odd numbers and use bitfield for "is 2i+1 composite"
    const u64 SBS = 8*sizeof(u64);                          // sieve bucket size (in bits)
    u64* is_composite = calloc(N/(2*SBS) + 1, sizeof(u64));

#define BIT(i) (UINT64_C(1) << (i))

    is_composite[0] = 1;                                    // 1 is not prime
    // i, ii - are indices in the sieve; j,jj are the corresponding odd numbers (2*i+1 and 2*ii+1)
    for(u64 i = 1, j = 3 ; j*j <= N ; ++i, j += 2) {
        const u64 ib = i / SBS, ie = i - ib*SBS;            // bucket index and in-bucket bit index
        if(is_composite[ib] & BIT(ie)) continue;            // known composite - all multiples are marked already
        // all composite numbers below j^2 have at least one multiple less than j - they were already marked
        for(u64 jj = j*j, ii = (jj-1)/2 ; jj <= N ; ii += j, jj += 2*j) {
            const u64 iib = ii / SBS, iie = ii - iib*SBS;   // bucket index and in-bucket bit index
            is_composite[iib] |= BIT(iie);                  // mark composite
        }
    }

    uint count = 1;                                         // include 2 immediately
    for(u64 i = 0, in = (N-1)/2 ; i <= in ; ++i) {
        const u64 ib = i / SBS, ie = i - ib*SBS;            // bucket index and in-bucket bit index
        if((is_composite[ib] & BIT(ie)) == 0) ++count;      // prime
    }

    u64* prime = malloc(sizeof(u64[count+1]));
    prime[0] = 0; prime[1] = 2;                             // we sieve only odd numbers so treat 2 separately
    for(u64 i = 0, in = (N-1)/2, pi = 2 ; i <= in ; ++i) {
        const u64 ib = i / SBS, ie = i - ib*SBS;            // bucket index and in-bucket bit index
        if((is_composite[ib] & BIT(ie)) == 0) {             // prime
            assert(pi <= count); prime[pi++] = 2*i+1;
        }
    }

#undef BIT

    free(is_composite);
    *out_prime = prime;
    return count;
}

static inline void mark_composite_with_primes(u64 start, uint N, bool is_composite[static N],
    uint p_start, uint p_end, const u64 prime[static p_end]
) {
    const u64 end = start + N;
    for(uint i = p_start ; i < p_end ; ++i) { const u64 p = prime[i];
        for(u64 j = p * ((start+p-1)/p) ; j < end ; j += p)
            is_composite[j-start] = true;
    }
}

/*                                                                              */
/* Binary Search                                                                */
/*                                                                              */

/* search sorted (ascending) array for lowest index i, such that v <= a[i]  */
/* this index can be used for insertion (to keep array sorted)              */
/* comparison function acts like < (less) and is expected to be CMP(T,VT)   */
/*   i.e. if element is less than key                                       */

#define IMPL_LOWER_BOUND_VTYPE(T, VT, NAME, CMP)                                                \
    static inline uint NAME(uint n, const T elem[const static n], VT v) {                       \
        uint low = 0, high = n;                     /* start with whole array               */  \
        while(high > low) {                         /* more than one candidate              */  \
            const uint mid = low + (high - low) / 2;/* pick middle                          */  \
            if(CMP(elem[mid], v))   low = mid + 1;  /* a[mid] <  v: no need to consider mid */  \
            else                    high = mid;     /* a[mid] >= v: pick first half         */  \
        }                                                                                       \
        return low;                                                                             \
    }

#define IMPL_LOWER_BOUND(T, NAME, CMP) IMPL_LOWER_BOUND_VTYPE(T,T,NAME,CMP)

#define __LOWER_BOUND_IMPL_LESS_NUM(a,b) ((a)<(b))
IMPL_LOWER_BOUND_VTYPE(u64,u64,lower_bound_u64,__LOWER_BOUND_IMPL_LESS_NUM)

/*                                                                              */
/* prime counting function for known primes                                     */
/*                                                                              */

static inline uint primepi_known(u64 x, uint prime_count, const u64 prime[prime_count+1]) {
    const uint idx = lower_bound_u64(prime_count+1, prime, x);
    return idx <= prime_count && x == prime[idx] ? idx : idx-1;
}

/*                                                                              */
/* SQRT                                                                         */
/*                                                                              */

static inline u64 usqrt_64(u64 n) {
    u64 bit_place = UINT64_C(1) << 31;
    while(bit_place > n) bit_place >>= 1;

    u64 root = 0; while(bit_place) {
        const u64 root_cand = root | bit_place, square_cand = root_cand*root_cand;
        if(square_cand <= n) root = root_cand;
        bit_place >>= 1;
    }
    return root;
}

/*                                                                              */
/* CBRT (Cubic root)                                                            */
/*                                                                              */

static inline u64 ucbrt_64(u64 n) {
    // modified code from hacker's delight, applying 64-bit specific fixes
    u64 y = 0;
    for(int s = 63 ; s >= 0 ; s -= 3) {
        y = 2*y;
        u64 b = 3*y*(y+1) + 1;
        if((n >> s) >= b) {
            n -= (b << s); ++y;
        }
    }
    return y;
}

/*                                                                              */
/* min/max                                                                      */
/*                                                                              */

static inline u64 min_u64(u64 a, u64 b) { return a < b ? a : b; }
static inline u64 max_u64(u64 a, u64 b) { return a > b ? a : b; }

static inline uint min_u(uint a, uint b) { return a < b ? a : b; }
static inline uint max_u(uint a, uint b) { return a > b ? a : b; }

/*                                                                              */
/* alpha param                                                                  */
/*                                                                              */

static inline double get_default_alpha(u64 X) {
    const u64 X_6 = usqrt_64(ucbrt_64(X));
    const double lx = log((double)X), lx2 = lx*lx, lx3 = lx2*lx;
    const double a = 0.000681*lx3 - 0.011846*lx2 + 0.044074*lx + 0.988365;
    return a < 1 ? 1 : (a > X_6 ? X_6 : a);
}

//
// 128 bit implementation
//

#define T i128
#define TADD(x,y)       ((x)+(y))
#define TSUB(x,y)       ((x)-(y))
#define TMUL(x,y)       ((x)*(y))
#define TCONV_64(x)     ((T)(x))
#define TCONV_128(x)    (x)
#define TNEG(x)         (-(x))

#define F_POW 0
#define NAME count_primes_128
#include "sum_over_primes.inl"
#undef NAME
#undef F_POW

#define F_POW 1
#define NAME sum_primes_128
#include "sum_over_primes.inl"
#undef NAME
#undef F_POW

#define F_POW 2
#define NAME sum_primes2_128
#include "sum_over_primes.inl"
#undef NAME
#undef F_POW

#define F_POW 3
#define NAME sum_primes3_128
#include "sum_over_primes.inl"
#undef NAME
#undef F_POW

#define F_POW 4
#define NAME sum_primes4_128
#include "sum_over_primes.inl"
#undef NAME
#undef F_POW

#undef T
#undef TADD
#undef TSUB
#undef TMUL
#undef TCONV_64
#undef TCONV_128
#undef TNEG

//
// 256 bit implementation
//

#define T i256
#define TADD(x,y)       i256_add(x,y)
#define TSUB(x,y)       i256_sub(x,y)
#define TMUL(x,y)       i256_mul(x,y)
#define TCONV_64(x)     i256_from_i64(x)
#define TCONV_128(x)    i256_from_i128(x)
#define TNEG(x)         i256_neg(x)

#pragma omp declare reduction(+ : i256 : omp_out = TADD(omp_out,omp_in)) initializer(omp_priv=(i256){})

#define F_POW 0
#define NAME count_primes_256
#include "sum_over_primes.inl"
#undef NAME
#undef F_POW

#define F_POW 1
#define NAME sum_primes_256
#include "sum_over_primes.inl"
#undef NAME
#undef F_POW

#define F_POW 2
#define NAME sum_primes2_256
#include "sum_over_primes.inl"
#undef NAME
#undef F_POW

#define F_POW 3
#define NAME sum_primes3_256
#include "sum_over_primes.inl"
#undef NAME
#undef F_POW

#define F_POW 4
#define NAME sum_primes4_256
#include "sum_over_primes.inl"
#undef NAME
#undef F_POW

#undef T
#undef TADD
#undef TSUB
#undef TMUL
#undef TCONV_64
#undef TCONV_128
#undef TNEG
