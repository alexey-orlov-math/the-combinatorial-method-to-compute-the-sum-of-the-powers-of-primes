#include <stdlib.h>
#include <stdio.h>
#include <process.h>

#include "sum_over_primes.h"
#include "sum_over_primes.c"

#define ARRAY_COUNT(arr)    (sizeof(arr)/sizeof(arr[0]))
#define I256S(s)            string_to_i256(strlen((s)), (s))

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

static inline double find_optimal_alpha(u64 X, uint fp, i256 FX, uint TC, double alpha_0, double alpha_1, double step) {
    const u64 X_2 = usqrt_64(X), X_3 = ucbrt_64(X);

    double min_alpha = -1, min_dt = 1e6;
    for(double a = alpha_0 ; a < alpha_1 ; a += step) {
        const u64 Y = a*X_3; if(Y >= X_2) break;

        char cmdline[1024] = {};
        sprintf(cmdline, "primepowsum.exe --time -p=%u -a=%f -t=8 %"PRIu64, fp, a, X);

        FILE* pipe = _popen(cmdline, "r");
        char res_str[I256_DEC_STRLEN]; float dt = 0;
        fscanf(pipe, "%s\n%f", res_str, &dt);
        _pclose(pipe);

        const i256 res = string_to_i256(strlen(res_str), res_str);
        assert(i256_equal(res,FX));

        // const u64 t0 = time_system_time();
        // const i128 res = sum_primes_128_openmp_impl(X,TC, Y,X_3,7);
        // // const i128 res = sum_primes_128_impl(X, Y,X_3,7);
        // assert(res == FX);
        // const u64 t1 = time_system_time();
        // const double dt = time_seconds_from_system(t1-t0);

        if(dt < min_dt) {
            // if(dt < 0.95*min_dt) {          // update only if the change is "significant"
            //     min_dt = dt; min_alpha = a;
            // }
            min_dt = dt; min_alpha = a;
        } else {
            if(dt > 1.5*min_dt) {           // stop if we are too far away from the optimum
                break;
            }
            if(a > min_alpha + 3*step) {    // we didnt find better alpha for a long enough interval
                break;
            }
        }
    }

    return min_alpha;
}

typedef struct AlphaInput {
    u64 X; const char* pp; i256 FX;
    double alpha_0, alpha_1; double step;
} AlphaInput;

static inline void find_alpha_1() {
    AlphaInput v[] = {
        {1*upow_64(10,9),  "10^9",    I256S("24739512092254535"), 1,5,0.5 },
        {1*upow_64(10,10), "10^10",   I256S("2220822432581729238"), 1,10,0.5 },
        {5*upow_64(10,10), "5*10^10", I256S("51814385436874673047"), 1,10,0.5 },
        {1*upow_64(10,11), "10^11",   I256S("201467077743744681014"), 4,10,0.5 },
        {5*upow_64(10,11), "5*10^11", I256S("4729818810421710631151"), 5,10,0.5 },
        {1*upow_64(10,12), "10^12",   I256S("18435588552550705911377"), 5,10,0.5 },
        {5*upow_64(10,12), "5*10^12", I256S("435063475540860317775484"), 5,10,0.5 },
        {1*upow_64(10,13), "10^13",   I256S("1699246443377779418889494"), 8,15,0.5 },
        {5*upow_64(10,13), "5*10^13", I256S("40277474360243260520915546"), 10,15,0.5 },
        {1*upow_64(10,14), "10^14",   I256S("157589260710736940541561021"), 10,15,0.5 },
        // {5*upow_64(10,14), "5*10^14", I256S("3749486993444942506049588783"), 10,20,1 },
        // {1*upow_64(10,15), "10^15",   I256S("14692398516908006398225702366"), 10,20,1 },
        // {5*upow_64(10,15), "5*10^15", I256S("350719744105802406156634876203"), 10,20,1 },
        // {1*upow_64(10,16), "10^16",   I256S("1376110854313351899159632866552"), 15,25,2 },
        // {5*upow_64(10,16), "5*10^16", I256S("32943259925529434020725057253548"), 15,25,1 },
        // {1*upow_64(10,17), "10^17",   I256S("129408626276669278966252031311350"), 15,25,2 },
        // {1*upow_64(10,18), "10^18",   I256S("12212914292949226570880576733896687"), 15,25,2 },
    };

    const uint TC = 8;
    for(uint i = 0 ; i < ARRAY_COUNT(v) ; ++i) {
       printf("%-8s: %.3f\n", v[i].pp, find_optimal_alpha(v[i].X, 1, v[i].FX, TC, v[i].alpha_0,v[i].alpha_1, v[i].step));
       fflush(stdout);
    }
}

static inline void find_alpha_2() {
    AlphaInput v[] = {
        {1*upow_64(10,10), "10^10",   I256S("14692485666215945973239505690"), 1,20,1 },
        {5*upow_64(10,10), "5*10^10", I256S("1714863031171407826702942323341"), 1,20,1 },
        {1*upow_64(10,11), "10^11",   I256S("13338380640732671147186590712800"), 1,20,1 },
        {5*upow_64(10,11), "5*10^11", I256S("1566398144419578032981266419280441"), 1,20,1 },
        {1*upow_64(10,12), "10^12",   I256S("12212907966177661747436156685876997"), 5,20,1 },
        {5*upow_64(10,12), "5*10^12", I256S("1441593988892141564900337100187358316"), 5,20,1 },
        {1*upow_64(10,13), "10^13",   I256S("11262617785640702236670513970349205634"), 5,20,1 },
        {5*upow_64(10,13), "5*10^13", I256S("1335210125295770298473184342618020082018"), 5,20,1 },
    };

    const uint TC = 8;
    for(uint i = 0 ; i < ARRAY_COUNT(v) ; ++i) {
       printf("%-8s: %.3f\n", v[i].pp, find_optimal_alpha(v[i].X, 2, v[i].FX, TC, v[i].alpha_0,v[i].alpha_1, v[i].step));
       fflush(stdout);
    }
}


int main() {
    find_alpha_1();
    // find_alpha_2();

    return 0;
}
