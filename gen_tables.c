#include <stdlib.h>
#include <stdio.h>
#include <process.h>

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

static inline void print_results(uint fp, uint vc, const u64 v[static vc], const char* l[static vc]) {
    printf("\\vspace{7mm}\n");
    printf("\\begin{table}\n");
    printf("\\caption{Values of $F_%u(x)$}\n", fp);
    printf("\\centering\n");
    printf("\\begin{tabular}{|r|r|} \\hline\n");

    char cmdline[1024] = {};
    for(uint i = 0 ; i < vc ; ++i) {
        sprintf(cmdline, "primepowsum.exe -p=%u -t=8 %"PRIu64, fp, v[i]);

        FILE* pipe = _popen(cmdline, "r");
        char res_str[I256_DEC_STRLEN]; fscanf(pipe, "%s\n", res_str);
        _pclose(pipe);

        printf("$%s$ & %s \\\\ \\hline\n", l[i], res_str);
    }

    printf("\\end{tabular}\n");
    printf("\\end{table}\n\n");
}

int main() {

    // for n = 2,3 we have enough bits for all those values
    {
        const u64 X[] = {
            upow_64(10,10), 5*upow_64(10,10), upow_64(10,11), 5*upow_64(10,11), upow_64(10,12), 5*upow_64(10,12),
            upow_64(10,13), 5*upow_64(10,13), upow_64(10,14), 5*upow_64(10,14), upow_64(10,15), 5*upow_64(10,15),
            upow_64(10,16), 5*upow_64(10,16), upow_64(10,17),
        };
        const char* L[] = {
            "10^{10}", "5 \\times 10^{10}", "10^{11}", "5 \\times 10^{11}", "10^{12}", "5 \\times 10^{12}",
            "10^{13}", "5 \\times 10^{13}", "10^{14}", "5 \\times 10^{14}", "10^{15}", "5 \\times 10^{15}",
            "10^{16}", "5 \\times 10^{16}", "10^{17}",
        };

        for(uint fp = 2 ; fp <= 3 ; ++fp)
            print_results(fp, ARRAY_COUNT(X), X, L);
    }

    // for n = 4, 256 bits are not enough past 10^15 (very crude estimation)
    {
        const u64 X[] = {
            upow_64(10,10), 5*upow_64(10,10), upow_64(10,11), 5*upow_64(10,11), upow_64(10,12), 5*upow_64(10,12),
            upow_64(10,13), 5*upow_64(10,13), upow_64(10,14), 5*upow_64(10,14), upow_64(10,15),
        };
        const char* L[] = {
            "10^{10}", "5 \\times 10^{10}", "10^{11}", "5 \\times 10^{11}", "10^{12}", "5 \\times 10^{12}",
            "10^{13}", "5 \\times 10^{13}", "10^{14}", "5 \\times 10^{14}", "10^{15}",
        };

        print_results(4, ARRAY_COUNT(X), X, L);
    }

    return 0;
}
