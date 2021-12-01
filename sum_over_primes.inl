
/*                                                                              */
/* Macro Magic                                                                  */
/*                                                                              */

#define TADD3(x,y,z)    TADD(x,TADD(y,z))
#define TADD4(x,y,z,w)  TADD(TADD(x,y),TADD(z,w))
#define TMUL3(x,y,z)    TMUL(x,TMUL(y,z))
#define TMUL4(x,y,z,w)  TMUL(TMUL(x,y),TMUL(z,w))

#define CONCAT(x,y) x ## y

#define FENWICK_SUM_NAME(P) CONCAT(P,_fenwick_sum)
#define FENWICK_SUM FENWICK_SUM_NAME(NAME)
#define FENWICK_ADD_NAME(P) CONCAT(P,_fenwick_add)
#define FENWICK_ADD FENWICK_ADD_NAME(NAME)

#define PHI_ORD_ITER_NAME(P) CONCAT(P,_phi_ord_iter)
#define PHI_ORD_ITER PHI_ORD_ITER_NAME(NAME)
#define PHI_ORD_NAME(P) CONCAT(P,_phi_ord)
#define PHI_ORD PHI_ORD_NAME(NAME)

#define F_POW0_NAME(P) CONCAT(P,_f_pow0)
#define F_POW0 F_POW0_NAME(NAME)
#define F_POW1_NAME(P) CONCAT(P,_f_pow1)
#define F_POW1 F_POW1_NAME(NAME)
#define F_POW2_NAME(P) CONCAT(P,_f_pow2)
#define F_POW2 F_POW2_NAME(NAME)
#define F_POW3_NAME(P) CONCAT(P,_f_pow3)
#define F_POW3 F_POW3_NAME(NAME)
#define F_POW4_NAME(P) CONCAT(P,_f_pow4)
#define F_POW4 F_POW4_NAME(NAME)

#define F_POW0_SUM_NAME(P) CONCAT(P,_f_pow0_sum)
#define F_POW0_SUM F_POW0_SUM_NAME(NAME)
#define F_POW1_SUM_NAME(P) CONCAT(P,_f_pow1_sum)
#define F_POW1_SUM F_POW1_SUM_NAME(NAME)
#define F_POW2_SUM_NAME(P) CONCAT(P,_f_pow2_sum)
#define F_POW2_SUM F_POW2_SUM_NAME(NAME)
#define F_POW3_SUM_NAME(P) CONCAT(P,_f_pow3_sum)
#define F_POW3_SUM F_POW3_SUM_NAME(NAME)
#define F_POW4_SUM_NAME(P) CONCAT(P,_f_pow4_sum)
#define F_POW4_SUM F_POW4_SUM_NAME(NAME)

/*                                                                              */
/* 1-based Fenwick Tree                                                         */
/*                                                                              */
static inline T FENWICK_SUM (uint n, T tree[static n], uint idx) {
    T ret = TCONV_64(0);
    while(idx >= 1) {
        ret  = TADD(ret, tree[idx]);
        idx -= idx & -idx;
    }
    return ret;
}
static inline void FENWICK_ADD (uint n, T tree[static n], uint idx, T x) {
    while(idx <= n) {
        tree[idx] = TADD(tree[idx], x);
        idx += idx & -idx;
    }
}

/*                                                                              */
/* common functions to sum over primes                                          */
/*                                                                              */

// For both f and summatory F, we will call them with the arguments below X, thus u64 is sufficient.
// On the other hand, the return type must be T, and all the internal calculations must be done with T
//   but we hide this complexity inside, taking u64 input.
// Also, the summatory functions are the only place where we would need to divide values of type T.
// To simplify the implementation we will check the divisibility of input and divide terms of numerator before
//   doing converting to T and doing final multiplication.

static inline T F_POW0 (u64 n) {
    return TCONV_64(1);
}
static inline T F_POW1 (u64 n) {
    return TCONV_64(n);
}
static inline T F_POW2 (u64 n) {
    const T nt = TCONV_64(n); return TMUL(nt,nt);
}
static inline T F_POW3 (u64 n) {
    const T nt = TCONV_64(n); return TMUL3(nt,nt,nt);
}
static inline T F_POW4 (u64 n) {
    const T nt = TCONV_64(n), nt2 = TMUL(nt,nt); return TMUL(nt2,nt2);
}

static inline T F_POW0_SUM (u64 n) {
    return TCONV_64(n);
}
static inline T F_POW1_SUM (u64 n_) {
    // n * (n + 1) / 2
    const u128 n = n_;
    u128 t1 = n, t2 = n+1;
    switch(n % 2) {
        case 0: t1 /= 2; break;
        case 1: t2 /= 2; break;
    }
    return TMUL(TCONV_128(t1), TCONV_128(t2));
}
static inline T F_POW2_SUM (u64 n_) {
    // n * (n+1) * (2*n+1) / 6
    const u128 n = n_;
    u128 t1 = n, t2 = n+1, t3 = 2*n+1;
    switch(n % 2) {
        case 0: t1 /= 2; break;
        case 1: t2 /= 2; break;
    }
    switch(n % 3) {
        case 0: t1 /= 3; break;
        case 1: t3 /= 3; break;
        case 2: t2 /= 3; break;
    }
    return TMUL3(TCONV_128(t1),TCONV_128(t2),TCONV_128(t3));
}
static inline T F_POW3_SUM (u64 n) {
    // 1^3 + ... + n^3 = (1 + ... + n)^2
    const T s = F_POW1_SUM(n); return TMUL(s,s);
}
static inline T F_POW4_SUM (u64 n_) {
    // n * (n+1) * (2*n+1) * (3*n^2 + 3*n - 1) / 30
    const u128 n = n_;
    u128 t1 = n, t2 = n+1, t3 = 2*n+1, t4 = 3*n*n + 3*n - 1;
    switch(n % 2) {
        case 0: t1 /= 2; break;
        case 1: t2 /= 2; break;
    }
    switch(n % 3) {
        case 0: t1 /= 3; break;
        case 1: t3 /= 3; break;
        case 2: t2 /= 3; break;
    }
    switch(n % 5) {
        case 0: t1 /= 5; break;
        case 1: t4 /= 5; break;
        case 2: t3 /= 5; break;
        case 3: t4 /= 5; break;
        case 4: t2 /= 5; break;
    }
    return TMUL4(TCONV_128(t1),TCONV_128(t2),TCONV_128(t3),TCONV_128(t4));
}

/*                                                                              */
/* determine which function we need to sum                                      */
/*                                                                              */

#ifndef F_POW
    #error F_POW must be defined
#endif

#if F_POW == 0
    #define F(n)        F_POW0(n)
    #define F_SUM(n)    F_POW0_SUM(n)
#elif F_POW == 1
    #define F(n)        F_POW1(n)
    #define F_SUM(n)    F_POW1_SUM(n)
#elif F_POW == 2
    #define F(n)        F_POW2(n)
    #define F_SUM(n)    F_POW2_SUM(n)
#elif F_POW == 3
    #define F(n)        F_POW3(n)
    #define F_SUM(n)    F_POW3_SUM(n)
#elif F_POW == 4
    #define F(n)        F_POW4(n)
    #define F_SUM(n)    F_POW4_SUM(n)
#else
    #error F_POW must be integer in [0,4]
#endif

/*                                                                              */
/* Common implementation code, independent of threading used                    */
/*                                                                              */

//
// phi: ordinary nodes
//
static inline T PHI_ORD_ITER (u64 X, u64 Y, uint K, u64 P, u64 n, const T (*phi_t)[F_POW+1][P], i8 mu_n) {
    const u64 x = X/n, q = x/P, r = x - q*P; const T qt = TCONV_64(q), Pt = TCONV_64(P);

    // k: 0 -> F_POW
#define FSV(q,FS) (q ? FS(q-1) : TCONV_64(0))
    #if F_POW == 0
        u64 binom[] = {1};          T fs[] = {qt};
    #elif F_POW == 1
        u64 binom[] = {1,1};        T fs[] = {qt, FSV(q, F_POW1_SUM)};
    #elif F_POW == 2
        u64 binom[] = {1,2,1};      T fs[] = {qt, FSV(q, F_POW1_SUM), FSV(q, F_POW2_SUM)};
    #elif F_POW == 3
        u64 binom[] = {1,3,3,1};    T fs[] = {qt, FSV(q, F_POW1_SUM), FSV(q, F_POW2_SUM), FSV(q, F_POW3_SUM)};
    #elif F_POW == 4
        u64 binom[] = {1,4,6,4,1};  T fs[] = {qt, FSV(q, F_POW1_SUM), FSV(q, F_POW2_SUM), FSV(q, F_POW3_SUM), FSV(q, F_POW4_SUM)};
    #else
        #error F_POW must be integer in [0,4]
    #endif
#undef FSV

    T qq = TCONV_64(1), PP = TCONV_64(1), phi = TCONV_64(0);
    for(uint k = 0 ; k <= F_POW ; ++k, qq = TMUL(qq,qt), PP = TMUL(PP,Pt)) {
        const T t = TADD(TMUL(fs[k],(*phi_t)[F_POW-k][P-1]), TMUL(qq,(*phi_t)[F_POW-k][r]));
        phi = TADD(phi, TMUL3(TCONV_64(binom[k]), PP, t));
    }

    // phi_O += f(n)*phi(x/n)*mu(n)
    const T ret = TMUL(F(n),phi);
    return mu_n > 0 ? ret : TNEG(ret);
}

static inline T PHI_ORD (u64 X, u64 Y, uint K, uint prime_count, const u64 prime[static prime_count+1],
    const uint phi_lpf[static Y+1], const i8 phi_mu[static Y+1]
) {
    // calculate P = m_K
    u64 phi_P = 1; for(uint i = 1 ; i <= K ; ++i) phi_P *= prime[i];

    // generate table
#define GEN_PHI_TABLE_FOR_F(TBL,F) do{                                                              \
TBL[0] = TCONV_64(0); for(uint x = 1 ; x < phi_P ; ++x) TBL[x] = TADD(TBL[x-1],F(x));               \
for(uint k = 1 ; k <= K ; ++k) {                                                                    \
    const u64 p = prime[k]; const T fp = F(p);                                                      \
    /* note that we want to calculate in place. thus we go backwards so that x/p has old value */   \
    for(uint x = phi_P-1 ; x > 0 ; --x) { TBL[x] = TSUB(TBL[x], TMUL(fp,TBL[x/p])); }               \
}                                                                                                   \
} while(0)

    T (*phi_t)[F_POW+1][phi_P] = malloc(sizeof(*phi_t));
    GEN_PHI_TABLE_FOR_F((*phi_t)[0], F_POW0);
    if(F_POW >= 1) GEN_PHI_TABLE_FOR_F((*phi_t)[1], F_POW1);
    if(F_POW >= 2) GEN_PHI_TABLE_FOR_F((*phi_t)[2], F_POW2);
    if(F_POW >= 3) GEN_PHI_TABLE_FOR_F((*phi_t)[3], F_POW3);
    if(F_POW >= 4) GEN_PHI_TABLE_FOR_F((*phi_t)[4], F_POW4);

#undef GEN_PHI_TABLE_FOR_F

    // go through ordinary nodes and calculate their contribution using phi table
    T phi_O = PHI_ORD_ITER(X,Y,K,phi_P, 1,phi_t,1);
    for(u64 n = 2 ; n <= Y ; ++n) if(phi_mu[n] && phi_lpf[n] > K) {
        const T phi_iter = PHI_ORD_ITER(X,Y,K,phi_P, n,phi_t,phi_mu[n]);
        phi_O = TADD(phi_O,phi_iter);
    }

    free(phi_t);
    return phi_O;
}

/*                                                                              */
/* Implementation                                                                */
/*                                                                              */

#include "sum_over_primes_st.inl"
#ifdef _OPENMP
#include "sum_over_primes_openmp1.inl"
#include "sum_over_primes_openmp2.inl"
#endif // _OPENMP

#define FNAIVE_NAME(FN) CONCAT(FN,_naive)
#define FST_NAME(FN) CONCAT(FN,_st)
#define FSTIMPL_NAME(FN) CONCAT(FN,_st_impl)
#define FMP_NAME(FN) CONCAT(FN,_openmp)
#define FMPIMPL_NAME(FN) CONCAT(FN,_openmp_impl1)

// used to avoid corner case handling with X too low
static inline T FNAIVE_NAME(NAME) (u64 X) {
    u64* prime = 0; const uint prime_count = generate_prime_list_up_to(X, &prime);

    T ret = TCONV_64(0);
    for(uint i = 1 ; i <= prime_count ; ++i)
        ret = TADD(ret, F(prime[i]));

    free(prime);
    return ret;
}

static inline T FST_NAME(NAME) (u64 X) {
    if(X <= 100) return FNAIVE_NAME(NAME)(X);

    const u64 X_3 = ucbrt_64(X);
    const double a = get_default_alpha(X);
    return FSTIMPL_NAME(NAME) (X, a*X_3, X_3, 7);
}
static inline T FMP_NAME(NAME) (u64 X, uint t) {
    if(X <= 100) return FNAIVE_NAME(NAME)(X);

    const u64 X_3 = ucbrt_64(X);
    const double a = get_default_alpha(X);
    return FMPIMPL_NAME(NAME) (X, t, a*X_3, X_3, 7);
}

#undef FNAIVE_NAME
#undef FST_NAME
#undef FSTIMPL_NAME
#undef FMP_NAME
#undef FMPIMPL_NAME

#undef TADD3
#undef TADD4
#undef TMUL3
#undef TMUL4

#undef CONCAT

#undef FENWICK_SUM_NAME
#undef FENWICK_SUM
#undef FENWICK_ADD_NAME
#undef FENWICK_ADD

#undef F_POW0_NAME
#undef F_POW0
#undef F_POW1_NAME
#undef F_POW1
#undef F_POW2_NAME
#undef F_POW2
#undef F_POW3_NAME
#undef F_POW3
#undef F_POW4_NAME
#undef F_POW4

#undef F_POW0_SUM_NAME
#undef F_POW0_SUM
#undef F_POW1_SUM_NAME
#undef F_POW1_SUM
#undef F_POW2_SUM_NAME
#undef F_POW2_SUM
#undef F_POW3_SUM_NAME
#undef F_POW3_SUM
#undef F_POW4_SUM_NAME
#undef F_POW4_SUM

#undef F
#undef F_SUM
