gen_test_values(N) = {
    s0 = 0; s1 = 0; s2 = 0; s3 = 0; s4 = 0;
    forprime(p=1,N, s0 = s0+1; s1 = s1+p; s2 = s2+p^2; s3 = s3+p^3; s4 = s4+p^4);
    printf("0: %d\n1: %d\n2: %d\n3: %d\n4: %d\n", s0,s1,s2,s3,s4);
}

printf("10^6\n"); gen_test_values(10^6);
printf("10^7\n"); gen_test_values(10^7);
printf("10^8\n"); gen_test_values(10^8);
printf("10^9\n"); gen_test_values(10^9);
printf("10^10\n"); gen_test_values(10^10);

quit;
