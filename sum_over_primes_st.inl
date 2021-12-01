
#define FIMPL_NAME(FN) CONCAT(FN,_st_impl)

/*                                                                              */
/* Sum over Primes single threaded implementation                               */
/*                                                                              */

static inline T FIMPL_NAME(NAME) (u64 X, u64 Y, u64 B, uint K) {
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
    K = K < A_x4 ? K : A_x4 - 1;

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

    T P2_S1 = TCONV_64(0), P2_S2 = TCONV_64(0);
    T phi_O = TCONV_64(0), phi_S1 = TCONV_64(0), phi_S2 = TCONV_64(0), phi_S3 = TCONV_64(0);
    T phi_S1_1 = TCONV_64(0), phi_S1_2 = TCONV_64(0), F_Y = TCONV_64(0), F_x_4 = TCONV_64(0), phi_S1_Fp_accum = TCONV_64(0);
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
    // allocate all the things we need
    //

    T* phi_last = calloc(A_x4+1, sizeof(T));        // phi values for the last elements of the previous block [0..a]
    bool* block_sieve = malloc(sizeof(bool[B]));    // we proceed in blocks [l, l+B); we mark *composite* numbers
    T* block_tree = malloc(sizeof(T[B+1]));         // as above BUT we have one-based indexing

    T* Fn_holder = malloc(sizeof(T[B+1]));          // if the block starts at prime, we will access F[-1]
    T* Fn = Fn_holder+1;

    bool* q_sieve = malloc(sizeof(bool[B]));        // sieve for P_2 calculation: x/q_i; we mark *composite* numbers

    //
    // blocks processing
    //
    for(uint bi = 0 ; ; ++bi) {
        const u64 l = 1 + bi*B, h = l+B;
        if(l*Y > X) break;

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
            const u64 nl = max_u64(X/(h*p) + 1, 1), nh = min_u64(X/(l*p), Y);
            if(nl <= nh) {
                for(u64 n = nl ; n <= nh ; ++n) if(phi_mu[n] && phi_lpf[n] > i) {
                    if(n*p <= Y) continue; // we are rounding down so we could have extended our interval to the left
                    const u64 nn = X/(n*p); assert(l <= nn && nn < h);
                    // S3 += f(p*n) * mu(p*n) + phi_block(x/(n*p))
                    const T phi_block = TADD(phi_last[i], FENWICK_SUM(B, block_tree, nn-l+1));
                    const T v = TMUL3(fp,F(n),phi_block);
                    phi_S3 = phi_mu[n] > 0 ? TSUB(phi_S3, v) : TADD(phi_S3, v);
                }
            }

            // we can update phi value to reflect the current block
            phi_last[i] = TADD(phi_last[i], FENWICK_SUM(B, block_tree, B));

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
        Fn[-1] = Fn_accum;
        for(uint n = 0 ; n < B ; ++n) {
            if(!block_sieve[n]) Fn_accum = TADD(Fn_accum, F(l+n));
            Fn[n] = Fn_accum;
        }

        //
        // phi: S1
        //

        if(l <= Y && h > x_4) {
            // check if we can init the needed values: F(Y), F(x^(1/4)), and the accumulator for F(p-1)
            const u64 p_x3 = prime[A_x3];
            if(l <= p_x3 && p_x3 < h)   phi_S1_Fp_accum = Fn[(int)(p_x3-l)-1];
            if(l <= x_4 && x_4 < h)     F_x_4 = Fn[x_4-l];
            if(l <= Y && Y < h)         F_Y = Fn[Y-l];

            uint lpi = primepi_known(l, prime_count,prime); if(prime[lpi] != l) ++lpi;
            uint hpi = primepi_known(h, prime_count,prime); if(prime[hpi] == h) --hpi;
            for(uint i = max_u(lpi, A_x4+1), in = min_u(hpi,A) ; i <= in ; ++i) {
                const u64 p = prime[i]; const T fp = F(p);
                const int idx = p-l; assert(idx >= 0);
                phi_S1_1 = TSUB(phi_S1_1, TMUL(fp,Fn[idx]));                // sum[f(p) * F(p)]
                if(i <= A_x3)
                    phi_S1_1 = TADD(phi_S1_1, TMUL3(fp,Fn[idx],Fn[idx-1])); // sum[f(p) * F(p) * F(p-1)]
                if(i <= A_x_div_y_2)
                    phi_S1_2 = TSUB(phi_S1_2, TMUL(fp,Fn[idx-1]));          // sum[f(p) * F(p-1)]
            }

            // sum[f(p) * F(p-1) * F(x/p^2)]
            // we want to traverse "backwards" (from bigger primes to smaller)
            //   to be able to update values of F(p-1) without extra bookkeeping

            // we want sqrt[x/h] < p <= sqrt[x/l] and sqrt[x/Y] < p <= x^(1/3)
            const u64 ph = min_u64(usqrt_64(X/l), x_3), pl = max_u64(usqrt_64(X/h), x_div_y_2);
            const uint pa = primepi_known(pl, prime_count,prime)+1, pb = primepi_known(ph, prime_count,prime);
            for(int i = pb ; i >= pa ; --i) {
                const u64 p = prime[i]; if(p <= x_div_y_2) break;
                // sum[f(p) * F(p-1) * F(x/p^2)]
                phi_S1_1 = TSUB(phi_S1_1, TMUL3(F(p), phi_S1_Fp_accum, Fn[X/(p*p)- l]));
                // we keep F(p-1) in the accumulator, thus subtract previous prime
                phi_S1_Fp_accum = TSUB(phi_S1_Fp_accum, F(prime[i-1]));
            }

            // check if we can wrap up S1 calculation
            // F(Y)*[F(Y) - F(x^(1/4))] + S_11 + F(Y)*S_12
            if(l <= Y && Y < h)
                phi_S1 = TADD3(phi_S1_1, TMUL(phi_S1_2, F_Y), TMUL(F_Y, TSUB(F_Y, F_x_4)));
        }

        //
        // phi: S2
        //

        if(l <= x_2) {
            for(uint i = A_x4 + 1 ; i <= A_x3 ; ++i) {
                const u64 p = prime[i];

                // handling the split at i = sqrt[X/Y]
                const bool case1 = i <= A_x_div_y_2;
                const u64 ql = max_u64(X/(p*h), case1 ? p+1 : p);
                const u64 qh = min_u64(X/(p*l), case1 ? Y : X/(p*p));
                if(ql > qh) continue;   // handle empty block for q

                const uint qa = primepi_known(ql, prime_count,prime)+1, qb = primepi_known(qh, prime_count,prime);
                if(qa <= qb) {          // handle empty block for q
                    T sp = TCONV_64(0);
                    for(uint j = qa ; j <= qb ; ++j) {
                        const u64 q = prime[j], x_pq = X/(p*q);
                        // sp = sum[f(q) * F(x/pq), x/pq is in block]
                        sp = TADD(sp, TMUL(F(q), Fn[x_pq-l]));
                    }
                    // phi_S2 = sum[f(p) * sp, x/pq is in block]
                    phi_S2 = TADD(phi_S2, TMUL(F(p),sp));
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
                    // P2_S1 = sum[f(q) * F(x/q)]
                    P2_S1 = TADD(P2_S1, TMUL(F(q), Fn[X/q-l]));
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
                // P2_S2 = sum[f(p) * F(p-1)]
                P2_S2 = TADD(P2_S2, TMUL(F(p), Fn[i-1]));
            }
        }
    }

    free(q_sieve); free(Fn_holder); free(block_tree); free(block_sieve); free(phi_last);
    free(phi_mu); free(phi_lpf); free(prime);

    const T P2 = TSUB(P2_S1, P2_S2), phi = TADD4(phi_O, phi_S1, phi_S2, phi_S3);
    const T ret = TADD(phi, TSUB(F_Y, TADD(P2, TCONV_64(1))));

    return ret;
}

#undef FIMPL_NAME
