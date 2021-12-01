#pragma once

#include <stdint.h>
#include <assert.h>
#include <stdbool.h>
#include <stdint.h>
#include <inttypes.h>
#include <stdlib.h>
#include <string.h>

typedef uint8_t u8; typedef int8_t i8; typedef unsigned int uint;
typedef uint64_t u64; typedef int64_t i64; typedef unsigned __int128 u128; typedef __int128 i128;

/*                                                                              */
/* 256-bit integer                                                              */
/*                                                                              */

typedef struct i256 {
    u128 l, h;
} i256;

static inline i256 i256_from_i64(i64 v) {
    return (i256){.l = v, .h = v < 0 ? -1 : 0};
}
static inline i256 i256_from_i128(i128 v) {
    return (i256){.l = v, .h = v < 0 ? -1 : 0};
}
static inline bool i256_equal(i256 a, i256 b) {
    return a.l == b.l && a.h == b.h;
}
static inline i256 i256_add(i256 a, i256 b) {
    u128 l = a.l + b.l, h = a.h + b.h;
    if(l < a.l) ++h;
    return (i256){.l = l, .h = h};
}
static inline i256 i256_sub(i256 a, i256 b) {
    u128 l = a.l - b.l, h = a.h - b.h;
    if(l > a.l) --h;
    return (i256){.l = l, .h = h};
}
static inline i256 i256_mul(i256 a, i256 b) {
    // (aB+b)*(cB+d) = (ac)B^2 + (ad+bc)B + bd
    // as we ignore the overflow, for B=2^128 we ignore the B^2 term and any overflow coming from the B term
    // but we need to be careful with the low parts multiplication, as the overflow here should be handled
    // we consider this product in a similar way, with B=2^64

    // optimize the case when lower parts of both fit into 64 bits
    // and higher parts are all zeros or all ones
    // the good chunk of multiplications would be like this
    const u128 u64_max = ~UINT64_C(0);
    if(a.l < u64_max && (a.h == 0 || ~a.h == 0) && b.l < u64_max && (b.h == 0 || ~b.h == 0)) {
        return (i256){.l = a.l*b.l, .h = a.h == b.h ? 0 : ~(u128)0};
    }

    const u128 mask = ~UINT64_C(0);
    const u128 A = a.l >> 64, B = a.l & mask, C = b.l >> 64, D = b.l & mask;
    const u128 AD = A*D, BC = B*C;
    const i256 t1 = (i256){.l = B*D, .h = a.l*b.h + a.h*b.l + A*C};
    // (ad)B and (bc)B with B=2^64 - the result wont fit into 128 bit
    const i256 t2 = (i256){.l = AD << 64, .h = AD >> 64}, t3 = (i256){.l = BC << 64, .h = BC >> 64};
    return i256_add(t1, i256_add(t2,t3));
}
static inline i256 i256_neg(i256 v) {
    u128 l = ~v.l + 1, h = ~v.h;
    if(l == 0) ++h;
    return (i256){.l = l, .h = h};
}
static inline i256 i256_udiv_u64(i256 a, u64 b, u64* rem_out) {
    assert((a.h & ((u128)1<<127)) == 0);

    // we implement simplistic "multiprecision division by digit" algorithm, assuming u64 digits
    const u128 B = (u128)1 << 64, mask = ~UINT64_C(0);
    u128 d[] = {(a.h >> 64) & mask, a.h & mask, (a.l >> 64) & mask, a.l & mask};

    u64 rem = 0;
    for(int i = 0 ; i < 4 ; ++i) {
        const u128 cur = d[i] + rem*B;
        d[i] = cur / b; rem = cur % b;
    }

    if(rem_out) *rem_out = rem;
    return (i256){.h = (d[0] << 64) | d[1], .l = (d[2] << 64) | d[3]};
}

/*                                                                              */
/* Input/Output of 128-bit integers                                             */
/*                                                                              */

// log_10(2^128-1) = 38.5, hence 39 digits max, +1 for EOL, +1 for sign
#define I128_DEC_STRLEN 41
static inline void i128_to_string(i128 n, char str[static I128_DEC_STRLEN]) {
    uint out_start = 0;
    if(n < 0) {
        n = -n; str[0] = '-'; out_start = 1;
    }

    u8 tmp[I128_DEC_STRLEN]; int dcount = 0;
    while(n) {
        tmp[dcount++] = n % 10;
        n /= 10;
    }
    assert(dcount < I128_DEC_STRLEN);
    if(dcount == 0) {
        tmp[dcount++] = 0;
    }

    for(int i = out_start, j = dcount - 1 ; j >= 0 ; ++i, --j)
        str[i] = tmp[j] + '0';
    str[dcount+out_start] = 0;
}

static inline i128 string_to_i128(uint n, const char digits[static n]) {
    i128 ret = 0;
    const bool neg = digits[0] == '-';
    for(int i = neg ? 1 : 0 ; i < n ; ++i) {
        assert('0' <= digits[i] && digits[i] <= '9');
        ret = 10*ret + (digits[i]-'0');
    }
    return neg ? -ret : ret;
}

/*                                                                              */
/* Input/Output of 256-bit integers                                             */
/*                                                                              */

// log_10(2^128-1) = 77.06, hence 78 digits max, +1 for EOL, +1 for sign
#define I256_DEC_STRLEN 80
static inline void i256_to_string(i256 n, char str[static I256_DEC_STRLEN]) {
    uint out_start = 0;
    if((n.h & ((u128)1<<127)) != 0) {
        n = i256_neg(n); str[0] = '-'; out_start = 1;
    }

    u8 tmp[I256_DEC_STRLEN]; int dcount = 0;
    while(n.h || n.l) {
        u64 rem = 0; n = i256_udiv_u64(n,10,&rem);
        tmp[dcount++] = rem;
    }
    assert(dcount < I256_DEC_STRLEN);
    if(dcount == 0) {
        tmp[dcount++] = 0;
    }

    for(int i = out_start, j = dcount - 1 ; i < dcount ; ++i, --j)
        str[i] = tmp[j] + '0';
    str[dcount] = 0;
}

static inline i256 string_to_i256(uint n, const char digits[static n]) {
    i256 ret = i256_from_i64(0);
    const bool neg = digits[0] == '-';
    for(int i = neg ? 1 : 0 ; i < n ; ++i) {
        assert('0' <= digits[i] && digits[i] <= '9');
        ret = i256_add(i256_mul(i256_from_i64(10),ret), i256_from_i64(digits[i]-'0'));
    }
    return neg ? i256_neg(ret) : ret;
}

