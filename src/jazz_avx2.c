
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include <immintrin.h>

// 1, 4, 5, 10, 11

#define KYBER_Q 3329

#define KYBER_N 256
#define KYBER_K 3

extern
void poly_compress4_avx2_jazz(uint8_t r[128], const int16_t a[KYBER_N]);
extern
void polyvec_compress10_avx2_jazz(uint8_t r[320], const int16_t a[KYBER_N]);

extern
void poly_round_reduce_avx2(int16_t a[KYBER_N]);

int16_t reduce(int16_t a){
    a %= KYBER_Q;
    if(a > (KYBER_Q / 2)){
        a -= KYBER_Q;
    }
    if(a < -(KYBER_Q / 2)){
        a += KYBER_Q;
    }
    return a;
}

void poly_reduce(int16_t a[KYBER_N]){
    for(size_t i = 0; i < KYBER_N; i++){
        a[i] = reduce(a[i]);
    }
}

int16_t pmullw(const int16_t a, const int16_t b){
    return a * b;
}

int16_t pmulhw(const int16_t a, const int16_t b){
    return (int16_t)(((int32_t)a * b) >> 16);
}

int16_t psraw(const int16_t a, const size_t i){
    if(i == 0){
        return a;
    }
    return a >> i;
}

int16_t pmulhrsw(const int16_t a, const int16_t b){
    return (int16_t)(( (int32_t)a * b + (1 << 14)) >> 15);
}

int16_t compress_D(int16_t a, const size_t D){
    if(a < 0){
        a += KYBER_Q;
    }
    return (int16_t)(( ( ((int32_t)a) << D) + (KYBER_Q / 2) ) / KYBER_Q) & ((1 << D) - 1);
}

int16_t Barrett_quotient_1(int16_t a){
    // 19-bit suffices for D = 1
    // 315 = round(2 * 2^19 / q)
    // return ((int16_t)(((int32_t)a * 315 + (1 << 18)) >> 19)) & 0x1;
    return pmulhrsw(pmulhw(a, 315), (1 << 12)) & 0x1;
    // return psraw(pmulhw(a, 315) + (1 << 2), 3) & 0x1;
}

int16_t Barrett_quotient_4(int16_t a){
    // 16-bit suffices for D = 4
    // 315 = round(16 * 2^16 / q)
    // return ((int16_t)(((int32_t)a * 315 + (1 << 15)) >> 16)) & 0xf;
    // return ((int16_t)(((int32_t)a * 630 + (1 << 16)) >> 17)) & 0xf;
    return pmulhrsw(pmulhw(a, 630), (1 << 14)) & 0xf;
    // return psraw(pmulhw(a, 630) + (1 << 0), 1) & 0xf;
}

int16_t Barrett_quotient_5(int16_t a){
    // 15-bit suffices for D = 5
    // 315 = round(32 * 2^15 / q)
    // return ((int16_t)(((int32_t)a * 315 + (1 << 14)) >> 15)) & 0x1f;
    return pmulhrsw(a, 315) & 0x1f;
}

// 1290167 = -20553 + 20 * 2^16

int16_t Barrett_quotient_10(int16_t a){
    // this doesn't work
    // return ((int16_t)(((int32_t)a * 161271 + (1 << 18)) >> 19)) & 0x3ff;
    // 22-bit suffices for D = 10
    // 1290167 = round(1024 * 2^22 / q)
    // beware that adding prior to shifting overflows (32-bit), we must shift, add, and then shift here.
    // return ((int16_t)( ((((int32_t)a * 1290167) >> 2) + (1 << 19)) >> 20)) & 0x3ff;
    return pmulhrsw(pmulhw(a, -20553) + pmullw(a, 20), (1 << 9)) & 0x3ff;
}

int16_t Barrett_quotient_11(int16_t a){
    // 21-bit suffices for D = 11
    // 1290167 = round(2048 * 2^21 / q)
    // beware that adding prior to shifting overflows (32-bit), we must shift, add, and then shift here.
    // return ((int16_t)( ((((int32_t)a * 1290167) >> 2) + (1 << 18)) >> 19)) & 0x7ff;
    return pmulhrsw(pmulhw(a, -20553) + pmullw(a, 20), (1 << 10)) & 0x7ff;
}

void poly_compress1(uint8_t r[32], const int16_t a[KYBER_N]){

    unsigned int i,j;
    int16_t u;

    for(i=0;i<KYBER_N/8;i++) {
        r[i] = 0;
        for(j=0;j<8;j++) {
            u = a[8*i+j];

            // 19-bit precision suffices for round(2 x / q)
            // inputs are in [-q/2, ..., q/2]
            // 315 = round(2 * 2^19 / q)
            u = Barrett_quotient_1(u);

            // this is equivalent to first mapping to positive
            // standard representatives followed by
            // u = ((((uint16_t)u << 1) + KYBER_Q/2)/KYBER_Q) & 1;

            r[i] |= u << j;

        }
    }

}

void poly_compress4(uint8_t r[128], const int16_t a[KYBER_N]){

    unsigned int j,k;
    int16_t u;

    uint16_t t[8];
    for(j=0;j<KYBER_N/8;j++) {
        for(k=0;k<8;k++) {
            u  = a[8*j+k];

            // 16-bit precision suffices for round(2^4 x / q)
            // inputs are in [-q/2, ..., q/2]
            // 315 = round(16 * 2^16 / q)
            t[k] = Barrett_quotient_4(u);

            // this is equivalent to first mapping to positive
            // standard representatives followed by
            // t[j] = ((((uint16_t)u << 4) + KYBER_Q/2)/KYBER_Q) & 0xf;

        }

        r[0] = t[0] | (t[1] << 4);
        r[1] = t[2] | (t[3] << 4);
        r[2] = t[4] | (t[5] << 4);
        r[3] = t[6] | (t[7] << 4);
        r += 4;
    }

}

void poly_compress5(uint8_t r[160], const int16_t a[KYBER_N]){

    unsigned int j,k;
    int16_t u;

    uint16_t t[8];
    for(j=0;j<KYBER_N/8;j++) {
        for(k=0;k<8;k++) {
            u  = a[8*j+k];

            // 15-bit precision suffices for round(2^5 x / q)
            // inputs are in [-q/2, ..., q/2]
            // 315 = round(32 * 2^15 / q)
            t[k] = Barrett_quotient_5(u);

            // this is equivalent to first mapping to positive
            // standard representatives followed by
            // t[j] = ((((uint32_t)u << 5) + KYBER_Q/2)/KYBER_Q) & 0x1f;

        }

        r[0] = (t[0] >> 0) | (t[1] << 5);
        r[1] = (t[1] >> 3) | (t[2] << 2) | (t[3] << 7);
        r[2] = (t[3] >> 1) | (t[4] << 4);
        r[3] = (t[4] >> 4) | (t[5] << 1) | (t[6] << 6);
        r[4] = (t[6] >> 2) | (t[7] << 3);
        r += 5;
    }

}

void poly_compress10(uint8_t r[320], const int16_t a[KYBER_N]){

    unsigned int j,k;
    int16_t u;

    uint16_t t[4];
    for(j=0;j<KYBER_N/4;j++) {
        for(k=0;k<4;k++) {
            u  = a[4*j+k];

            // 22-bit suffices for round(1024 x / q)
            // inputs are in [-q/2, ..., q/2]
            // 1290167 = round(1024 * 2^22 / q)
            t[k] = Barrett_quotient_10(u);

            // this is equivalent to first mapping to positive
            // standard representatives followed by
            // t[k]  = ((((uint32_t)u << 10) + KYBER_Q/2)/ KYBER_Q) & 0x3ff;

        }

        r[0] = (t[0] >> 0);
        r[1] = (t[0] >> 8) | (t[1] << 2);
        r[2] = (t[1] >> 6) | (t[2] << 4);
        r[3] = (t[2] >> 4) | (t[3] << 6);
        r[4] = (t[3] >> 2);
        r += 5;
    }

}

void poly_compress11(uint8_t r[352], const int16_t a[KYBER_N]){

    unsigned int j,k;
    int16_t u;

    uint16_t t[8];
    for(j=0;j<KYBER_N/8;j++) {
        for(k=0;k<8;k++) {
            u  = a[8*j+k];

            // 21-bit suffices for round(2048 x / q)
            // inputs are in [-q/2, ..., q/2]
            // 1290167 = round(2048 * 2^21 / q)
            t[k] = Barrett_quotient_11(u);

            // this is equivalent to first mapping to positive
            // standard representatives followed by
            // t[k]  = ((((uint32_t)u << 11) + KYBER_Q/2)/KYBER_Q) & 0x7ff;

        }

        r[ 0] = (t[0] >>  0);
        r[ 1] = (t[0] >>  8) | (t[1] << 3);
        r[ 2] = (t[1] >>  5) | (t[2] << 6);
        r[ 3] = (t[2] >>  2);
        r[ 4] = (t[2] >> 10) | (t[3] << 1);
        r[ 5] = (t[3] >>  7) | (t[4] << 4);
        r[ 6] = (t[4] >>  4) | (t[5] << 7);
        r[ 7] = (t[5] >>  1);
        r[ 8] = (t[5] >>  9) | (t[6] << 2);
        r[ 9] = (t[6] >>  6) | (t[7] << 5);
        r[10] = (t[7] >>  3);
        r += 11;
    }

}

int main(void){

    __attribute__ ((aligned(32))) int16_t a[KYBER_K * KYBER_N];
    int16_t poly_a[KYBER_N], poly_b[KYBER_N], poly_c[KYBER_N];
    uint8_t ref[KYBER_K * 352], res[KYBER_K * 352];

    for(size_t i = 0; i < KYBER_K * KYBER_N; i++){
        a[i] = rand() % KYBER_Q;
        a[i] -= KYBER_Q / 2;
    }

    poly_compress4(ref, a);
    poly_compress4_avx2_jazz(res, a);

    assert(memcmp(ref, res, 128) == 0);

    printf("poly_compress4 is correct!\n");

    for(size_t i = 0; i < KYBER_K; i++){
        poly_compress10(ref + i * 320, a + i * KYBER_N);
    }
    polyvec_compress10_avx2_jazz(res, a);

    assert(memcmp(ref, res, KYBER_K * 320) == 0);

    printf("poly_compress10 is correct!\n");

    for(size_t i = 0; i < KYBER_N; i++){
        poly_a[i] = poly_b[i] = rand();
    }

    poly_reduce(poly_a);
    poly_round_reduce_avx2(poly_b);

    assert(memcmp(poly_a, poly_b, sizeof(int16_t) * KYBER_N) == 0);

    printf("poly_round_reduce_avx2 is correct!\n");

    printf("Test finished!\n");

}

