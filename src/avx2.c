
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include <immintrin.h>

// 1, 4, 5, 10, 11

#define KYBER_Q 3329

#define KYBER_N 256

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
    // return ((int16_t)(((int32_t)a * 1290167 + (1 << 21)) >> 22)) & 0x3ff;
    return pmulhrsw(pmulhw(a, -20553) + pmullw(a, 20), (1 << 9)) & 0x3ff;
    // return psraw(pmulhw(a, -20553) + pmullw(a, 20) + (1 << 5), 6) & 0x3ff;

}

int16_t Barrett_quotient_11(int16_t a){
    // 21-bit suffices for D = 11
    // 1290167 = round(2048 * 2^21 / q)
    // return ((int16_t)(((int32_t)a * 1290167 + (1 << 20)) >> 21)) & 0x7ff;
    return pmulhrsw(pmulhw(a, -20553) + pmullw(a, 20), (1 << 10)) & 0x7ff;
    // return psraw(pmulhw(a, -20553) + pmullw(a, 20) + (1 << 4), 5) & 0x7ff;
}


// optimize this one
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
            u = (int16_t)(((int32_t)u * 315 + (1 << 18)) >> 19) & 1;

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
            t[k] = (int16_t)(((int32_t)u * 315 + (1 << 15)) >> 16) & 0xf;

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
            t[k] = (int16_t)(((int32_t)u * 315 + (1 << 14)) >> 15) & 0x1f;

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
            t[k] = ((int16_t)(((int32_t)u * 1290167 + (1 << 21)) >> 22)) & 0x3ff;

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
            t[k] = ((int16_t)(((int32_t)u * 1290167 + (1 << 20)) >> 21)) & 0x7ff;

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

void poly_compress1_avx2(uint8_t r[32], const int16_t a[KYBER_N]){

    unsigned int i,j;
    int16_t u;

    for(i=0;i<KYBER_N/8;i++) {
        r[i] = 0;
        for(j=0;j<8;j++) {
            u = a[8*i+j];

            // 19-bit precision suffices for round(2 x / q)
            // inputs are in [-q/2, ..., q/2]
            // 315 = round(2 * 2^19 / q)
            u = (int16_t)(((int32_t)u * 315 + (1 << 18)) >> 19) & 1;

            // this is equivalent to first mapping to positive
            // standard representatives followed by
            // u = ((((uint16_t)u << 1) + KYBER_Q/2)/KYBER_Q) & 1;

            r[i] |= u << j;

        }
    }

}

void poly_compress4_avx2(uint8_t r[128], const int16_t a[KYBER_N]){

    __m256i a0, a1, a2, a3;
    const __m256i b0 = _mm256_set1_epi16(630);
    const __m256i b1 = _mm256_set1_epi16(1 << 14);
    const __m256i mask4 = _mm256_set1_epi16(0xf);
    const __m256i shift = _mm256_set1_epi16((16 << 8) + 1);
    const __m256i shuffle = _mm256_set_epi32(7, 3, 6, 2, 5, 1, 4, 0);

    for(size_t i = 0; i < KYBER_N / 64; i++){

        a0 = _mm256_loadu_si256((__m256i*)(a + i * 64 + 0 * 16));
        a1 = _mm256_loadu_si256((__m256i*)(a + i * 64 + 1 * 16));
        a2 = _mm256_loadu_si256((__m256i*)(a + i * 64 + 2 * 16));
        a3 = _mm256_loadu_si256((__m256i*)(a + i * 64 + 3 * 16));

        // ((int16_t)(((int32_t)v * 630 + (1 << 16)) >> 17)) & 0xf;

        a0 = _mm256_mulhi_epi16(a0, b0);
        a1 = _mm256_mulhi_epi16(a1, b0);
        a2 = _mm256_mulhi_epi16(a2, b0);
        a3 = _mm256_mulhi_epi16(a3, b0);

        a0 = _mm256_mulhrs_epi16(a0, b1);
        a1 = _mm256_mulhrs_epi16(a1, b1);
        a2 = _mm256_mulhrs_epi16(a2, b1);
        a3 = _mm256_mulhrs_epi16(a3, b1);

        // below are the same as before
        a0 = _mm256_and_si256(a0, mask4);
        a1 = _mm256_and_si256(a1, mask4);
        a2 = _mm256_and_si256(a2, mask4);
        a3 = _mm256_and_si256(a3, mask4);

        a0 = _mm256_packus_epi16(a0, a1);
        a2 = _mm256_packus_epi16(a2, a3);

        a0 = _mm256_maddubs_epi16(a0, shift);
        a2 = _mm256_maddubs_epi16(a2, shift);

        a0 = _mm256_packus_epi16(a0, a2);

        a0 = _mm256_permutevar8x32_epi32(a0, shuffle);

        _mm256_storeu_si256((__m256i*)(r + i * 32), a0);

    }

}

void poly_compress5_avx2(uint8_t r[128], const int16_t a[KYBER_N]){

    __m256i a0, a1;
    __m128i lo, hi;
    const __m256i b0 = _mm256_set1_epi16(315);
    const __m256i mask5 = _mm256_set1_epi16(0x1f);
    const __m256i shift1 = _mm256_set1_epi16((32 << 8) + 1);
    const __m256i shift2 = _mm256_set1_epi32((1024 << 16) + 1);
    const __m256i sllv_indx = _mm256_set1_epi64x(12);
    const __m256i shuffle = _mm256_set_epi8( 8, -1, -1, -1, -1, -1,  4,  3,
                                             2,  1,  0, -1, 12, 11, 10,  9,
                                            -1, 12, 11, 10,  9,  8, -1, -1,
                                            -1, -1, -1,  4,  3,  2,  1,  0);

    for(size_t i = 0; i < KYBER_N / 32; i++){

        a0 = _mm256_loadu_si256((__m256i*)(a + i * 32 + 0 * 16));
        a1 = _mm256_loadu_si256((__m256i*)(a + i * 32 + 1 * 16));

        // ((int16_t)(((int32_t)v * 315 + (1 << 14)) >> 15)) & 0x1f;

        a0 = _mm256_mulhrs_epi16(a0, b0);
        a1 = _mm256_mulhrs_epi16(a1, b0);

        // below are the same as before
        a0 = _mm256_and_si256(a0, mask5);
        a1 = _mm256_and_si256(a1, mask5);

        a0 = _mm256_packus_epi16(a0, a1);

        a0 = _mm256_maddubs_epi16(a0, shift1);
        a0 = _mm256_madd_epi16(a0, shift2);

        a0 = _mm256_sllv_epi32(a0, sllv_indx);
        a0 = _mm256_srli_epi64(a0, 12);

        a0 = _mm256_shuffle_epi8(a0, shuffle);

        lo = _mm256_castsi256_si128(a0);
        hi = _mm256_extracti128_si256(a0, 1);

        lo = _mm_blendv_epi8(lo, hi, _mm256_castsi256_si128(shuffle));

        _mm_storeu_si128((__m128i*)(r + i * 20 + 0), lo);
        memcpy(r + i * 20 + 16, &hi, 4);

    }

}

void poly_compress10_avx2(uint8_t r[320], const int16_t a[KYBER_N]){

    __m256i a0;
    __m256i p0, p1;
    __m128i lo, hi;

    const __m256i b0 = _mm256_set1_epi16(-20553);
    const __m256i b1 = _mm256_set1_epi16(20);
    const __m256i b2 = _mm256_set1_epi16(1 << 9);
    const __m256i mask10 = _mm256_set1_epi16(0x3ff);
    const __m256i shift = _mm256_set1_epi32((1024 << 16) + 1);
    const __m256i sllv_indx = _mm256_set1_epi64x(12);
    const __m256i shuffle = _mm256_set_epi8( 8,  4,  3,  2,  1,  0, -1, -1,
                                            -1, -1, -1, -1, 12, 11, 10,  9,
                                            -1, -1, -1, -1, -1, -1, 12, 11,
                                            10,  9,  8,  4,  3,  2,  1,  0);

    for(size_t i = 0; i < KYBER_N / 16; i++){
        a0 = _mm256_loadu_si256((__m256i*)(a + i * 16));

        // ((int16_t)(((int32_t)v * 1290167 + (1 << 21)) >> 22)) & 0x3ff;

        p0 = _mm256_mulhi_epi16(a0, b0);
        p1 = _mm256_mullo_epi16(a0, b1);
        p0 = _mm256_add_epi16(p0, p1);
        p0 = _mm256_mulhrs_epi16(p0, b2);

        // below are the same as before
        a0 = _mm256_and_si256(p0, mask10);

        a0 = _mm256_madd_epi16(a0, shift);

        a0 = _mm256_sllv_epi32(a0, sllv_indx);
        a0 = _mm256_srli_epi64(a0, 12);

        a0 = _mm256_shuffle_epi8(a0, shuffle);

        lo = _mm256_castsi256_si128(a0);
        hi = _mm256_extracti128_si256(a0, 1);

        lo = _mm_blend_epi16(lo, hi, 0xe0);

        _mm_storeu_si128((__m128i*)(r + i * 20 + 0), lo);
        memcpy(r + i * 20 + 16, &hi, 4);

    }

}

void poly_compress11_avx2(uint8_t r[352], const int16_t a[KYBER_N]){

    __m256i a0, a1;
    __m256i p0, p1;
    __m128i lo, hi, thi;

    int16_t t[16];

    const __m256i b0 = _mm256_set1_epi16(-20553);
    const __m256i b1 = _mm256_set1_epi16(20);
    const __m256i b2 = _mm256_set1_epi16(1 << 10);
    const __m256i mask11 = _mm256_set1_epi16(0x7ff);
    const __m256i shift1 = _mm256_set1_epi32((2048 << 16) + 1);
    const __m256i sllv32_indx = _mm256_set_epi32(0, 10, 0, 10, 0, 10, 0, 10);
    const __m256i srlv64_indx = _mm256_set_epi64x(30, 10, 30, 10);
    const __m256i shuffle = _mm256_set_epi8( 4,  3,  2,  1,  0,  0, -1, -1,
                                            -1, -1, 10,  9,  8,  7,  6,  5,
                                            -1, -1, -1, -1, -1, 10,  9,  8,
                                             7,  6,  5,  4,  3,  2,  1,  0);

    for(size_t i = 0; i < 16; i++){

        a0 = _mm256_loadu_si256((__m256i*)(a + i * 16));

        // ((int16_t)(((int32_t)v * 1290167 + (1 << 20)) >> 21)) & 0x7ff;

        p0 = _mm256_mulhi_epi16(a0, b0);
        p1 = _mm256_mullo_epi16(a0, b1);
        p0 = _mm256_add_epi16(p0, p1);
        p0 = _mm256_mulhrs_epi16(p0, b2);

        // below are the same as before
        a0 = _mm256_and_si256(p0, mask11);

        a0 = _mm256_madd_epi16(a0, shift1);

        a0 = _mm256_sllv_epi32(a0, sllv32_indx);
        a1 = _mm256_bsrli_epi128(a0, 8);

        a0 = _mm256_srlv_epi64(a0, srlv64_indx);
        a1 = _mm256_slli_epi64(a1, 34);

        a0 = _mm256_add_epi64(a0, a1);
        a0 = _mm256_shuffle_epi8(a0, shuffle);

        lo = _mm256_castsi256_si128(a0);
        hi = _mm256_extracti128_si256(a0, 1);
        lo = _mm_blendv_epi8(lo, hi, _mm256_castsi256_si128(shuffle));

        _mm_storeu_si128((__m128i*)(r + i * 22 +  0), lo);
        memcpy(r + i * 22 + 16, &hi, 6);

    }

}

int main(void){

    __attribute__ ((aligned(32))) int16_t a[KYBER_N];
    uint8_t ref[352], res[352];

    for(size_t i = 0; i < KYBER_N; i++){
        a[i] = rand() % KYBER_Q;
        a[i] -= KYBER_Q / 2;
    }

    // poly_compress1(ref, a);
    // poly_compress1_avx2(res, a);

    // assert(memcmp(ref, res, 32) == 0);

    poly_compress4(ref, a);
    poly_compress4_avx2(res, a);

    assert(memcmp(ref, res, 128) == 0);

    poly_compress5(ref, a);
    poly_compress5_avx2(res, a);

    assert(memcmp(ref, res, 160) == 0);

    poly_compress10(ref, a);
    poly_compress10_avx2(res, a);

    assert(memcmp(ref, res, 320) == 0);

    poly_compress11(ref, a);
    poly_compress11_avx2(res, a);

    assert(memcmp(ref, res, 352) == 0);

    printf("Test finished!\n");

}

