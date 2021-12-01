We have three implementations: single-threaded, openmp (1) with threads doing whole blocks processing, and openmp (2)
where phi S2, phi S3 and P2 are calculated independently with openmp threads. It turned out that the openmp (2)
is a bit slower than the openmp (1), so we will use only the first one, but we keep the code around.
We want to have "template" implementation, and usually for this macro are used, but in our case
(lots of complicated code) this is quite unwieldy, especially since we want to be able to use preprocessor inside
implementation. Thus we structure code as follows:

sum_over_primes_st.inl      | single-threaded implementation
sum_over_primes_openmp1.inl | openmp implementation that has threads processing whole block at once
sum_over_primes_openmp2.inl | openmp implementation that has phi S2, phi S3 and P2 processed independently,
                            | each one in a multithreaded fashion
                            |
sum_over_primes.inl         | defines things dependent on the "template" type T and F_POW macro;
                            | also provides fenwick tree implementation,
                            | macro for the function f to sum over primes, and its summatory function;
                            | contains some code common to all implementations and includes implementation files above
                            |
sum_over_primes.c           | defines common routines, and then include sum_over_primes.inl once for each combo
                            | of T, F_POW macro we care about. It also has macro definitions for the implementation
                            | of arithmetic with the type T

This way we only need to compile sum_over_primes.c to get all the possible implementations. It also can be directly
included into the main app c file, as we do in all the following, giving us just one file to compile to get
the executable file.

primepowsum.c               | the main application that can calculate sum of x^k over primes for x <= X and integer
                            | k >= 0. We allow k in [0,3] currently.
find_optimal_alpha.c        | utility to quickly generate optimal alpha for different inputs, allowing us to generate
                            | polynomial in log(X) for the optimal Y parameter tweaking
gen_tables.c                | generate tables with numerical results
test.c                      | various correctness tests
