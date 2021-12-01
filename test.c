#include <stdio.h>

#include "sum_over_primes.h"
#include "sum_over_primes.c"

static inline u64 upow_64(u64 a, u64 b) {
    // russian peasant algorithm
    if(b == 0) return 1;                    // handle raising to 0-th power immediately
    while(b % 2 == 0) { a = a*a; b /= 2; }  // a^[2^k*b] = (a^[2^k])^b = (((a^2)..)^2)^b
    u64 ret = a; b /= 2;                    // b is odd: you can say we manually unrolled loop below
    for( ; b ; b /= 2) {
        a = a*a;                            // as above: square while pow is even
        if(b % 2 != 0) ret = a*ret;         // update result when pow is odd
    }
    return ret;
}

#define ARRAY_COUNT(arr)    (sizeof(arr)/sizeof(arr[0]))

#define I128S(s)    string_to_i128(strlen((s)), (s))
#define I256S(s)    string_to_i256(strlen((s)), (s))

#define CHECK_EQUAL_IMPL(T, NAME, T_EQ, T_STRLEN, T_TOSTR)          \
    static inline void NAME(T res, T expected) {                    \
        if(!T_EQ(res, expected)) {                                  \
            printf("FAIL!\n");                                      \
            char res_str[T_STRLEN]; T_TOSTR(res, res_str);          \
            char exp_str[T_STRLEN]; T_TOSTR(expected, exp_str);     \
            printf("expected: %s\ngot: %s\n", exp_str, res_str);    \
            exit(1);                                                \
        }                                                           \
    }

#define EQ_128(x,y) ((x) == (y))
CHECK_EQUAL_IMPL(i128,check_equal_128, EQ_128,     I128_DEC_STRLEN, i128_to_string)
CHECK_EQUAL_IMPL(i256,check_equal_256, i256_equal, I256_DEC_STRLEN, i256_to_string)
#undef EQ_128

#define TC 8

#define CALL_IMPL(T,S,ER)                                                                                                   \
    static inline T call_impl_param_ ## S(uint n, uint t, u64 X, u64 Y, u64 B, uint K) {                                    \
        switch(n) {                                                                                                         \
            case 0: return t > 1 ? count_primes_##S##_openmp_impl1(X,t,Y,B,K) : count_primes_##S##_st_impl(X,Y,B,K); break; \
            case 1: return t > 1 ? sum_primes_##S##_openmp_impl1(X,t,Y,B,K)   : sum_primes_##S##_st_impl(X,Y,B,K);   break; \
            case 2: return t > 1 ? sum_primes2_##S##_openmp_impl1(X,t,Y,B,K)  : sum_primes2_##S##_st_impl(X,Y,B,K);  break; \
            case 3: return t > 1 ? sum_primes3_##S##_openmp_impl1(X,t,Y,B,K)  : sum_primes3_##S##_st_impl(X,Y,B,K);  break; \
            case 4: return t > 1 ? sum_primes4_##S##_openmp_impl1(X,t,Y,B,K)  : sum_primes4_##S##_st_impl(X,Y,B,K);  break; \
        }                                                                                                                   \
        assert(false && "n must be in [0,4]");                                                                              \
        return ER;                                                                                                          \
    }                                                                                                                       \
    static inline T call_impl_ ## S(uint n, uint t, u64 X) {                                                                \
        switch(n) {                                                                                                         \
            case 0: return t > 1 ? count_primes_##S##_openmp(X,t) : count_primes_##S##_st(X); break;                        \
            case 1: return t > 1 ? sum_primes_##S##_openmp(X,t)   : sum_primes_##S##_st(X);   break;                        \
            case 2: return t > 1 ? sum_primes2_##S##_openmp(X,t)  : sum_primes2_##S##_st(X);  break;                        \
            case 3: return t > 1 ? sum_primes3_##S##_openmp(X,t)  : sum_primes3_##S##_st(X);  break;                        \
            case 4: return t > 1 ? sum_primes4_##S##_openmp(X,t)  : sum_primes4_##S##_st(X);  break;                        \
        }                                                                                                                   \
        assert(false && "n must be in [0,4]");                                                                              \
        return ER;                                                                                                          \
    }

CALL_IMPL(i128,128,0)
CALL_IMPL(i256,256,i256_from_i64(0))

//
// Test smaller known values, varying Y,B,K
//

static inline void test_1_impl(u64 X, uint n, i128 FX,
    u64 Y_start, u64 Y_end, u64 Y_step, u64 B_start, u64 B_end, u64 B_step
) {
    for(u64 Y = Y_start ; Y < Y_end ; Y += Y_step) for(u64 B = B_start ; B < B_end ; B += B_step) for(uint K = 0 ; K <= 7 ; ++K) {
        const i128 res_st = call_impl_param_128(n,1, X,Y,B,K), res_mp = call_impl_param_128(n,TC, X,Y,B,K);
        check_equal_128(res_st, FX); check_equal_128(res_mp, FX);
    }
}
static inline void test_1() {
    // note that in all the cases we have one value of Y less than x^(1/3) and greater than x^(1/2)
    // these should still work correctly: we should update Y to be in the expected interval
    test_1_impl(upow_64(10,6), 0, 78498, 80,1100,100, 80,1100,100);
    test_1_impl(upow_64(10,6), 1, 37550402023, 80,1100,100, 80,1100,100);
    test_1_impl(upow_64(10,6), 2, 24693298341834533, 80,1100,100, 80,1100,100);
    test_1_impl(upow_64(10,6), 3, I128S("18393235410255348281725"), 80,1100,100, 80,1100,100);
    test_1_impl(upow_64(10,6), 4, I128S("14652318678776560130517006257"), 80,1100,100, 80,1100,100);

    test_1_impl(upow_64(10,7), 0, 664579, 200,3300,150, 200,3300,150);
    test_1_impl(upow_64(10,7), 1, 3203324994356, 200,3300,150, 200,3300,150);
    test_1_impl(upow_64(10,7), 2, I128S("21113978675102768574"), 200,3300,150, 200,3300,150);
    test_1_impl(upow_64(10,7), 3, I128S("157468033886449157634628898"), 200,3300,150, 200,3300,150);
    test_1_impl(upow_64(10,7), 4, I128S("1255486194092640992925158230204098"), 200,3300,150, 200,3300,150);
}

//
// Test smaller known values, results are generated by pari/gp
//

static inline void test_2() {
    // for all the tests: as X is small, the value of Y is not that important
#define CHECK_128(n,X,FXS) do {                         \
    check_equal_128(call_impl_128(n,TC,X), I128S(FXS)); \
} while(0)
#define CHECK_256(n,X,FXS) do {                         \
    check_equal_256(call_impl_256(n,TC,X), I256S(FXS)); \
} while(0)
#define CHECK(n,X,FXS) do { CHECK_128(n,X,FXS); CHECK_256(n,X,FXS); } while(0)

    CHECK(0, upow_64(10,6), "78498");
    CHECK(1, upow_64(10,6), "37550402023");
    CHECK(2, upow_64(10,6), "24693298341834533");
    CHECK(3, upow_64(10,6), "18393235410255348281725");
    CHECK(4, upow_64(10,6), "14652318678776560130517006257");

    CHECK(0, upow_64(10,7), "664579");
    CHECK(1, upow_64(10,7), "3203324994356");
    CHECK(2, upow_64(10,7), "21113978675102768574");
    CHECK(3, upow_64(10,7), "157468033886449157634628898");
    CHECK(4, upow_64(10,7), "1255486194092640992925158230204098");

    CHECK(0, upow_64(10,8), "5761455");
    CHECK(1, upow_64(10,8), "279209790387276");
    CHECK(2, upow_64(10,8), "18433608754948081174274");
    CHECK(3, upow_64(10,8), "1375951520453237482162596242178");
    CHECK(4, upow_64(10,8), "109765422425729947150373907223888255454");

    CHECK(0, upow_64(10,9), "50847534");
    CHECK(1, upow_64(10,9), "24739512092254535");
    CHECK(2, upow_64(10,9), "16352255694497179054764665");
    CHECK(3, upow_64(10,9), "12212817300801335181685928278468661");
    CHECK_256(4, upow_64(10,9), "9745949414105865089440935358384227422742413");

    CHECK(0, upow_64(10,10), "455052511");
    CHECK(1, upow_64(10,10), "2220822432581729238");
    CHECK(2, upow_64(10,10), "14692485666215945973239505690");
    CHECK(3, upow_64(10,10), "109780001885333601058528339379120755908");
    CHECK_256(4, upow_64(10,10), "876279913324387539183894015044723229219045342750");

    CHECK(0, upow_64(10,11), "4118054813");
    CHECK(1, upow_64(10,11), "201467077743744681014");
    CHECK(0, 5*upow_64(10,11), "19308136142");
    CHECK(1, 5*upow_64(10,11), "4729818810421710631151");

    CHECK(0, upow_64(10,12), "37607912018");
    CHECK(1, upow_64(10,12), "18435588552550705911377");
    CHECK(0, 5*upow_64(10,12), "177291661649");
    CHECK(1, 5*upow_64(10,12), "435063475540860317775484");

#undef CHECK
#undef CHECK_256
#undef CHECK_128
}

//
// Test larger values, by calculating F(N) and F(N+eps), sieving [N,N+eps], and checking that we get
//   the same value for F(N+eps) if we recalculate it from known F(N) and primes from the interval
//

static inline void test_3_impl(uint fp, u64 eps, uint xc, const u64 X[static xc]) {
    // generate primes up to 10^9. This is enough to fully sieve any interval in [1,10^18]
    // we note that this takes insane amount of memory
    // thus we first calculate the contribution of primes in [N,N+eps], free the memory
    //   and only then we will actually compute F(N), F(N+eps) and compare the results
    u64* prime = 0; const uint prime_count = generate_prime_list_up_to(UINT64_C(1000000000), &prime);

    bool* sieve = malloc(sizeof(bool[eps]));
    i256 delta[xc];
    for(uint i = 0 ; i < xc ; ++i) {
        const u64 x = X[i];

        memset(sieve, 0x00, sizeof(bool[eps]));
        mark_composite_with_primes(x+1, eps, sieve, 1,prime_count+1,prime);

        delta[i] = i256_from_i64(0);
        for(uint j = 0 ; j < eps ; ++j) if(!sieve[j]) {
            const i256 p = i256_from_i64(x+1+j), p2 = i256_mul(p,p); i256 inc;
            switch(fp) {
                case 0: inc = i256_from_i64(1); break;
                case 1: inc = p;                break;
                case 2: inc = p2;               break;
                case 3: inc = i256_mul(p,p2);   break;
                case 4: inc = i256_mul(p2,p2);  break;
            }
            delta[i] = i256_add(delta[i], inc);
        }
    }
    free(prime);
    free(sieve);

    for(uint i = 0 ; i < xc ; ++i) {
        const u64 x = X[i];
        const i256 res1 = call_impl_256(fp,TC,x), res2 = call_impl_256(fp,TC,x+eps);
        const i256 res_delta = i256_sub(res2,res1);
        check_equal_256(res_delta, delta[i]);
    }
}

static inline void test_3() {
    const u64 eps = upow_64(10,6);

    {
        const u64 X[] = {
            upow_64(10,11), 5*upow_64(10,11), upow_64(10,12), 5*upow_64(10,12),
            upow_64(10,13), 5*upow_64(10,13), upow_64(10,14), 5*upow_64(10,14),
            upow_64(10,15), 5*upow_64(10,15), upow_64(10,16), 5*upow_64(10,16),
        };

        test_3_impl(1,eps,ARRAY_COUNT(X),X);
        test_3_impl(2,eps,ARRAY_COUNT(X),X);
        test_3_impl(3,eps,ARRAY_COUNT(X),X);
    }

    {
        const u64 X[] = {
            upow_64(10,11), 5*upow_64(10,11), upow_64(10,12), 5*upow_64(10,12),
            upow_64(10,13), 5*upow_64(10,13), upow_64(10,14), 5*upow_64(10,14),
        };

        test_3_impl(4,eps,ARRAY_COUNT(X),X);
    }
}


int main() {
    test_1();
    test_2();
    test_3();

    return 0;
}

