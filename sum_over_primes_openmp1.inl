
#define FMPIMPL1_NAME(FN) CONCAT(FN,_openmp_impl1)
#define FMPTHREAD1_NAME(FN) CONCAT(FN,_openmp_thread1)
#define FMPTHREADRESULT1_NAME(FN) CONCAT(FN,_ThreadResult1)

/*                                                                              */
/* Sum over Primes OpenMP implementation                                        */
/*   whole block processing is happening on one thread                          */
/*                                                                              */

//
// thread function for a block handling
//

typedef struct FMPTHREADRESULT1_NAME(NAME) {
    T phi_S2, phi_S3, P2;
} FMPTHREADRESULT1_NAME(NAME);

static inline FMPTHREADRESULT1_NAME(NAME) FMPTHREAD1_NAME(NAME) (
    u64 X, u64 Y, u64 x_2, uint K, uint A, uint A_x3, uint A_x_div_y_2, uint A_x4, u64 l, int B,
    uint prime_count, const u64 prime[static prime_count+1],
    const uint phi_lpf[static Y+1], const i8 phi_mu[static Y+1],
    T Fn_holder[static B+1], bool block_sieve[static B], T block_tree[static B+1], bool q_sieve[static B],
    T* out_Fn_block,
    T* out_P2_FI_block,
    T out_phi_last_block[static A_x4+1], T* out_phi_S2_block, T out_phi_S3_block[static A_x4+1]
) {
    memset(out_phi_last_block, 0x00, sizeof(T[A_x4+1]));
    memset(out_phi_S3_block, 0x00, sizeof(T[A_x4+1]));
    if(l*Y > X) {
        *out_Fn_block = *out_P2_FI_block = *out_phi_S2_block = TCONV_64(0);
        return (FMPTHREADRESULT1_NAME(NAME)) {};
    }

    T phi_S2 = TCONV_64(0), phi_S3 = TCONV_64(0), phi_S2_accum = TCONV_64(0), P2 = TCONV_64(0), P2_Fi_accum = TCONV_64(0);
    const u64 h = l + B;

    // sieve out first K primes
    memset(block_sieve, 0x00, sizeof(bool[B]));
    mark_composite_with_primes(l,B,block_sieve, 1,K+1, prime);

    //
    // phi: S3
    //

    // prepare fenwick tree for the calculation of phi for this block
    memset(block_tree, 0x00, sizeof(T[B+1]));
    for(uint i = 0 ; i < B ; ++i) if(!block_sieve[i]) FENWICK_ADD(B, block_tree, i+1, F(l+i));

    for(uint i = K+1 ; i <= A_x4 ; ++i) { const u64 p = prime[i]; const T fp = F(p);
        // we have phi calculated for b < i, thus we can find the contribution of i-th prime
        // interval for n*, intersected with [1,Y]; handle possibly empty block
        T phi_S3_block = TCONV_64(0);
        const u64 nl = max_u64(X/(h*p) + 1, 1), nh = min_u64(X/(l*p), Y);
        if(nl <= nh) {
            for(u64 n = nl ; n <= nh ; ++n) if(phi_mu[n] && phi_lpf[n] > i) {
                if(n*p <= Y) continue; // we are rounding down so we could have extended our interval to the left
                const u64 nn = X/(n*p); assert(l <= nn && nn < h);
                //
                const T v = TMUL(fp,F(n)), fs = FENWICK_SUM(B, block_tree, nn-l+1);
                phi_S3_block = phi_mu[n] > 0 ? TSUB(phi_S3_block, v) : TADD(phi_S3_block, v);
                phi_S3       = phi_mu[n] > 0 ? TSUB(phi_S3, TMUL(v,fs)) : TADD(phi_S3, TMUL(v,fs));
            }
        }

        // we can update phi value to reflect the current block
        out_phi_last_block[i] = FENWICK_SUM(B, block_tree, B);
        out_phi_S3_block[i]   = phi_S3_block;

        // "sieve" and prepare phi(n,i+1)
        for(u64 n = p * ((l+p-1)/p) ; n < h ; n += p) if(!block_sieve[n-l]) {
            block_sieve[n-l] = true;
            FENWICK_ADD(B, block_tree, n-l+1, TNEG(F(n)));
        }
    }

    // fully sieve block
    mark_composite_with_primes(l,B,block_sieve, A_x4+1,A+1, prime);

    // note that for the blocks intersecting with [1,Y], when sieving interval with primes
    //   we will strike out primes themselves.
    // while for S3 this is fine (as we go through by definition, checking that (n,m_k) = 1)
    //   for everything below we need to restore primes in the sieve
    if(l <= Y) {
        if(l == 1) block_sieve[0] = true; // we will never sieve out 1
        // unless l is prime we start at pi(l)+1 to be inside the block
        uint ps = primepi_known(l, prime_count, prime); if(prime[ps] != l) ++ps;
        for(uint i = ps ; i <= A ; ++i) {
            const u64 p = prime[i]; if(p >= h) break;
            block_sieve[p-l] = false;
        }
    }

    // update F(n) [summatory function] for the block
    T* Fn = Fn_holder+1; T Fn_accum = Fn[-1] = TCONV_64(0);
    for(uint n = 0 ; n < B ; ++n) {
        if(!block_sieve[n]) Fn_accum = TADD(Fn_accum, F(l+n));
        Fn[n] = Fn_accum;
    }

    //
    // phi: S2
    //

    if(l <= x_2) {
        for(uint i = A_x4 + 1 ; i <= A_x3 ; ++i) {
            const u64 p = prime[i]; const T fp = F(p);

            // handling the split at i = sqrt[X/Y]
            const bool case1 = i <= A_x_div_y_2;
            const u64 ql = max_u64(X/(p*h), case1 ? p+1 : p);
            const u64 qh = min_u64(X/(p*l), case1 ? Y : X/(p*p));
            if(ql > qh) continue;   // handle empty block for q

            const uint qa = primepi_known(ql, prime_count,prime)+1, qb = primepi_known(qh, prime_count,prime);
            if(qa <= qb) {          // handle empty block for q
                T sp = TCONV_64(0), spa = TCONV_64(0);
                for(uint j = qa ; j <= qb ; ++j) {
                    const u64 q = prime[j], x_pq = X/(p*q);
                    const T fq = F(q);
                    //
                    sp  = TADD(sp, TMUL(fq, Fn[x_pq-l]));
                    spa = TADD(spa, fq);
                }
                //
                phi_S2       = TADD(phi_S2, TMUL(fp,sp));
                phi_S2_accum = TADD(phi_S2_accum, TMUL(fp,spa));
            }
        }
    }

    //
    // P_2: S1
    //

    if(h > x_2) {
        // find primes q such that x/q is in our block
        const u64 ql = max_u64(X/h+1, Y+1), qh = min_u64(X/l, x_2);
        if(ql <= qh) {  // non-empty block for q
            const u64 len = qh-ql+1; assert(len < B);
            memset(q_sieve, 0x00, sizeof(bool[len]));
            mark_composite_with_primes(ql,len,q_sieve, 1,A+1, prime);

            for(uint i = 0 ; i < len ; ++i) if(!q_sieve[i]) {
                const u64 q = ql+i;
                if(q <= Y) continue; // we are rounding down so we could have extended our interval to the left
                const T fq = F(q);
                //
                P2          = TADD(P2, TMUL(fq, Fn[X/q-l]));
                P2_Fi_accum = TADD(P2_Fi_accum, fq);
            }
        }
    }

    //
    // P_2: S2
    //

    if(h > Y && l <= x_2) {
        for(int i = 0 ; i < B ; ++i) if(!block_sieve[i]) {
            const u64 p = l+i;
            if(p <= Y) continue;    // handle shorter interval (p too small: l < Y)
            if(p*p > X) break;      // handle shorter interval (p too large: h > sqrt[x])
            const T fp = F(p);
            //
            P2          = TSUB(P2, TMUL(fp,Fn[i-1]));
            P2_Fi_accum = TSUB(P2_Fi_accum, fp);
        }
    }

    *out_Fn_block = Fn_accum;
    *out_P2_FI_block = P2_Fi_accum;
    *out_phi_S2_block = phi_S2_accum;

    return (FMPTHREADRESULT1_NAME(NAME)){.phi_S2 = phi_S2, .phi_S3 = phi_S3, .P2 = P2};
}

static inline T FMPIMPL1_NAME(NAME) (u64 X, uint TC, u64 Y, u64 B, uint K) {
    //
    // compute important values, generate needed primes and tweak K if needed
    //
    const u64 x_2 = usqrt_64(X), x_3 = ucbrt_64(X);

    // check if Y is too large
    Y = Y < x_2 ? Y : x_2-1;
    // check if Y is too small
    Y = Y >= x_3 ? Y : x_3;
    u64* prime = 0; const uint prime_count = generate_prime_list_up_to(Y, &prime);

    const u64 x_4 = usqrt_64(x_2), x_div_y_2 = usqrt_64(X/Y);
    const uint A = primepi_known(Y, prime_count, prime), A_x3 = primepi_known(x_3, prime_count, prime);
    const uint A_x4 = primepi_known(x_4, prime_count, prime), A_x_div_y_2 = primepi_known(x_div_y_2, prime_count, prime);
    // check if K is too large [for the optimized special leaves processing we assume that pi(x^(1/4)) >= K
    K = K < A_x4 ? K : A_x4;

    //
    // precomputation
    //

    // least prime factor
    uint* phi_lpf = calloc(Y+1, sizeof(uint));
    for(uint i = 1 ; i <= A ; ++i) {
        const u64 p = prime[i];
        phi_lpf[p] =  i;                        // mark prime itself
        for(u64 j = p*p ; j <= Y ; j += p) {    // for n < p^2 there must be smaller prime factor
            if(!phi_lpf[j]) phi_lpf[j] = i;     // if visiting n first time: p is the least prime factor
        }
    }

    // mobius function
    i64* mu_calc = malloc(sizeof(i64[Y+1])); for(uint i = 0 ; i <= Y ; ++i) mu_calc[i] = 1;
    for(uint i = 1 ; i <= A ; ++i) {
        const u64 p = prime[i], p2 = p*p;
        for(u64 j = p ;  j <= Y ; j += p)   mu_calc[j] = -(i64)p*mu_calc[j];    // update mu (assuming squarefree)
        for(u64 j = p2 ; j <= Y ; j += p2)  mu_calc[j] = 0;                     // handle non-squarefree
    }

    i8* phi_mu = malloc(sizeof(i8[Y+1]));
    for(uint i = 0 ; i <= Y ; ++i) {
        const i64 mui = mu_calc[i];
        if(mui) {   // squarefree
            // did we process all the prime factors, or (exactly) one was omitted
            const i64 mui_abs = mui > 0 ? mui : -mui;
            phi_mu[i] = (mui_abs == i ? 1 : -1) * (mui > 0 ? 1 : -1);
        } else {    // non-squarefree
            phi_mu[i] = 0;
        }
    }
    free(mu_calc);

    //
    // intermediate values for the computation of P2 and phi
    //

    T P2 = TCONV_64(0);
    T phi_O = TCONV_64(0), phi_S1 = TCONV_64(0), phi_S2 = TCONV_64(0), phi_S3 = TCONV_64(0);
    T phi_S1_1 = TCONV_64(0), phi_S1_2 = TCONV_64(0), F_Y = TCONV_64(0), F_x_4 = TCONV_64(0);
    T Fn_accum = TCONV_64(0);

    //
    // phi: sum over ordinary nodes
    //

    if(K == 0) {
        for(u64 n = 1 ; n <= Y ; ++n) if(phi_mu[n]) {
            // phi_O += f(n)*f_sum(X/n)*mu(n)
            phi_O = TADD(phi_O, TMUL3(F(n), F_SUM(X/n), TCONV_64(phi_mu[n])));
        }
    } else {
        phi_O = PHI_ORD(X,Y,K, prime_count, prime, phi_lpf, phi_mu);
    }

    //
    // phi: S1
    // in multithreaded implementation it needs to be done separately
    //

    {
        T* Fp = malloc(sizeof(T[prime_count+1])); Fp[0] = TCONV_64(0);
        for(uint i = 1 ; i <= prime_count ; ++i) Fp[i] = TADD(Fp[i-1], F(prime[i]));

        for(uint i = A_x4+1 ; i <= A ; ++i) {
            const u64 p = prime[i]; const T fp = F(p);
            phi_S1_1 = TSUB(phi_S1_1, TMUL(fp,Fp[i]));

            if(p <= x_div_y_2) {            // sum[f(p) * F(p-1)]
                phi_S1_2 = TSUB(phi_S1_2, TMUL(fp,Fp[i-1]));
            }
            if(p > x_div_y_2 && p <= x_3) { // sum[f(p) * F(p-1) * F(x/p^2)]
                phi_S1_1 = TSUB(phi_S1_1, TMUL3(fp, Fp[i-1], Fp[primepi_known(X/(p*p), prime_count,prime)]));
            }
            if(p <= x_3) {                  // sum[f(p) * F(p) * F(p-1)]
                phi_S1_1 = TADD(phi_S1_1, TMUL3(fp,Fp[i],Fp[i-1]));
            }
        }

        F_Y = Fp[A], F_x_4 = Fp[A_x4];
        phi_S1 = TADD3(phi_S1_1, TMUL(phi_S1_2, F_Y), TMUL(F_Y, TSUB(F_Y, F_x_4)));

        free(Fp);
    }

    //
    // allocate all the things we need
    //

    // per-block per-thread summatory function values. We will access Fn[-1] if the block starts at prime, thus B+1 elements here
    T    (*Fn_holder)[TC][B+1]          = malloc(sizeof(*Fn_holder));
    // per-block per-thread  sieve: we mark *composite* numbers
    bool (*block_sieve)[TC][B]          = malloc(sizeof(*block_sieve));
    // per-block per-thread  one-based Fenwick tree
    T    (*block_tree)[TC][B+1]         = malloc(sizeof(*block_tree));
    // per-block per-thread  sieve for P_2 calculation: x/q_i; we mark *composite* numbers
    bool (*q_sieve)[TC][B]              = malloc(sizeof(*q_sieve));
    // phi values for the last elements of the previous block [0..a]
    // we will be accumulating them per-thread
    T    (*phi_last)[A_x4+1]            = calloc(A_x4+1, sizeof(T));
    T    (*phi_last_block)[TC][A_x4+1]  = malloc(sizeof(*phi_last_block));
    // per-thread accumulator of the summatory function
    T    Fn_accum_block[TC];
    // per thread block accumulator for P2
    T    p2_FI_block[TC];
    // per thread block accumulator for phi S2
    T    phi_S2_block[TC];
    // per-thread values of phi S3
    T    (*phi_S3_block)[TC][A_x4+1]    = malloc(sizeof(*phi_S3_block));

    //
    // blocks processing
    //
    for(uint bi = 0 ; ; bi += TC) {
        const u64 l = 1 + bi*B, h = l + TC*B;
        if(l*Y > X) break;

        //
        // threads processing
        //
        #pragma omp parallel for num_threads(TC) reduction(+: P2,phi_S2,phi_S3)
        for(uint t = 0 ; t < TC ; ++t) {
            FMPTHREADRESULT1_NAME(NAME) thread_res = FMPTHREAD1_NAME(NAME) (
                X,Y,x_2, K, A,A_x3,A_x_div_y_2,A_x4, l+t*B,B,
                prime_count,prime, phi_lpf,phi_mu,
                (*Fn_holder)[t], (*block_sieve)[t], (*block_tree)[t], (*q_sieve)[t],
                &Fn_accum_block[t],
                &p2_FI_block[t],
                (*phi_last_block)[t],  &phi_S2_block[t], (*phi_S3_block)[t]
            );

            P2      = TADD(P2, thread_res.P2);
            phi_S2  = TADD(phi_S2, thread_res.phi_S2);
            phi_S3  = TADD(phi_S3, thread_res.phi_S3);
        }

        // finalize per-block P2,phi values and update accumulators
        for(uint t = 0 ; t < TC ; ++t) {
            P2       = TADD(P2, TMUL(Fn_accum, p2_FI_block[t]));
            phi_S2   = TADD(phi_S2, TMUL(Fn_accum, phi_S2_block[t]));
            Fn_accum = TADD(Fn_accum, Fn_accum_block[t]);

            for(uint i = K+1 ; i <= A_x4 ; ++i) {
                phi_S3         = TADD(phi_S3, TMUL((*phi_S3_block)[t][i], (*phi_last)[i]));
                (*phi_last)[i] = TADD((*phi_last)[i], (*phi_last_block)[t][i]);
            }
        }
    }

    free(phi_S3_block); free(phi_last_block); free(phi_last);
    free(q_sieve); free(block_tree); free(block_sieve); free(Fn_holder);
    free(phi_mu); free(phi_lpf); free(prime);

    const T phi = TADD4(phi_O, phi_S1, phi_S2, phi_S3);
    const T ret = TADD(phi, TSUB(F_Y, TADD(P2, TCONV_64(1))));

    return ret;
}

#undef FMPIMPL1_NAME
#undef FMPTHREAD1_NAME
#undef FMPTHREADRESULT1_NAME
