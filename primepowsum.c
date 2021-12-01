#include <stdio.h>
#include <string.h>

#include "sum_over_primes.h"
#include "sum_over_primes.c"

static void     time_init_time();
static u64      time_system_time();
static double   time_seconds_from_system(u64 stm);

#define CALL_IMPL(T,S,ER)                                                                                                   \
    static inline T call_impl_ ## S(uint n, uint t, u64 X, u64 Y, u64 B, uint K) {                                          \
        switch(n) {                                                                                                         \
            case 0: return t > 1 ? count_primes_##S##_openmp_impl1(X,t,Y,B,K) : count_primes_##S##_st_impl(X,Y,B,K); break; \
            case 1: return t > 1 ? sum_primes_##S##_openmp_impl1(X,t,Y,B,K)   : sum_primes_##S##_st_impl(X,Y,B,K);   break; \
            case 2: return t > 1 ? sum_primes2_##S##_openmp_impl1(X,t,Y,B,K)  : sum_primes2_##S##_st_impl(X,Y,B,K);  break; \
            case 3: return t > 1 ? sum_primes3_##S##_openmp_impl1(X,t,Y,B,K)  : sum_primes3_##S##_st_impl(X,Y,B,K);  break; \
            case 4: return t > 1 ? sum_primes4_##S##_openmp_impl1(X,t,Y,B,K)  : sum_primes4_##S##_st_impl(X,Y,B,K);  break; \
        }                                                                                                                   \
        assert(false && "n must be in [0,4]");                                                                              \
        return ER;                                                                                                          \
    }

CALL_IMPL(i128,128,0)
CALL_IMPL(i256,256,i256_from_i64(0))

int main(int argc, char** argv) {
    u64 X = 0; uint fp = 0; uint tc = 0; float a = 0; bool output_time = false;

    for(uint i = 1 ; i < argc ; ++i) {
        if(argv[i][0] == '-') {
            if(strstr(argv[i], "-p") == argv[i]) {
                sscanf(argv[i], "-p=%u", &fp);
            } else if(strstr(argv[i], "--power") == argv[i]) {
                sscanf(argv[i], "--power=%u", &fp);
            } else if(strstr(argv[i], "-a") == argv[i]) {
                sscanf(argv[i], "-a=%f", &a);
            } else if(strstr(argv[i], "--alpha") == argv[i]) {
                sscanf(argv[i], "--alpha=%f", &a);
            } else if(strstr(argv[i], "-t") == argv[i]) {
                sscanf(argv[i], "-t=%u", &tc);
            } else if(strstr(argv[i], "--threads") == argv[i]) {
                sscanf(argv[i], "--threads=%u", &tc);
            } else if(strcmp(argv[i], "--time") == 0) {
                output_time = true;
            } else {
                fprintf(stderr, "unknown option %s\n", argv[i]);
                return 1;
            }
        } else {
            sscanf(argv[i], " %"PRIu64" ", &X);
        }
    }

    if(fp > 4) {
        fprintf(stderr, "wrong power specified: only integers in [0,4] are supported\n");
    }

    if(X == 0 || argc == 1) {
        printf( "Usage: primepowsum [options] X\n"
                "Sums power of primes less than or equal to X\n"
                "\n"
                "Options:\n"
                "\n"
                "-p, --power=NUM    The power exponent to use\n"
                "-a, --alpha=NUM    Tuning factor: Y = X^(1/3) * alpha\n"
                "-t, --threads=NUM  Number of threads\n"
                "    --time         Also print the time elapsed in seconds\n"
        );

        return 1;
    }

    time_init_time();

    // this is a very crude estimation: we play extra safe to avoid the possibility of internal overflow
    i256 res; bool use_256 = false;
    switch(fp) {
        case 0: case 1: use_256 = false;                        break;
        case 2: use_256 = X > UINT64_C(10000000000000);         break;  // 10^13
        case 3: use_256 = X > UINT64_C(10000000000);            break;  // 10^10
        case 4: use_256 = X > UINT64_C(1000000000);             break;  // 10^9
    }

    const u64 X_3 = ucbrt_64(X), X_6 = usqrt_64(X_3);
    if(a < 1 || a >= X_6) a = get_default_alpha(X);
    const u64 B = X_3, Y = a*X_3, K = 7;

    const u64 t0 = time_system_time();
    res = use_256 ? call_impl_256(fp, tc, X,Y,B,K) : i256_from_i128(call_impl_128(fp, tc, X,Y,B,K));
    const u64 t1 = time_system_time();

    char res_str[I256_DEC_STRLEN]; i256_to_string(res, res_str);
    printf("%s\n", res_str);
    if(output_time)
        printf("%.3f\n", time_seconds_from_system(t1-t0));

    return 0;
}

static u64 _time_cycles_per_sec = 0; static double _time_sec_per_cycle = 0;
static void time_init_time() {
#if defined(_WIN32) || defined(_WIN64) || defined(__WIN32__)
    int QueryPerformanceFrequency(void* lpFrequency);
    QueryPerformanceFrequency(&_time_cycles_per_sec);
#else
    #error Not Implemented
#endif
    _time_sec_per_cycle = 1.0 / _time_cycles_per_sec;
}

static u64 time_system_time() {
#if defined(_WIN32) || defined(_WIN64) || defined(__WIN32__)
    int QueryPerformanceCounter(void* lpPerformanceCount);
    u64 ret = 0; QueryPerformanceCounter(&ret);
    return ret;
#else
    #error Not Implemented
#endif
}

static double time_seconds_from_system(u64 stm) {
    return stm * _time_sec_per_cycle;
}

