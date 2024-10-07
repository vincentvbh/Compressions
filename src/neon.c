
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include <arm_neon.h>

#define KYBER_Q 3329

#define KYBER_N 256

int16_t shadd(const int16_t a, const int16_t b){
    return (a + b) >> 1;
}

int16_t mul(const int16_t a, const int16_t b){
    return a * b;
}

int16_t mla(int16_t c, const int16_t a, const int16_t b){
    return c + a * b;
}

int16_t mls(int16_t c, const int16_t a, const int16_t b){
    return c - a * b;
}

int16_t sqdmulh(const int16_t a, const int16_t b){
    return (int16_t)(((int32_t)a * b * 2) >> 16);
}

int16_t sqrdmulh(const int16_t a, const int16_t b){
    return (int16_t)(((int32_t)a * b * 2 + (1 << 15)) >> 16);
}

int16_t srshr(const int16_t a, const size_t i){
    if(i == 0){
        return a;
    }
    return ((a >> (i - 1)) + 1) >> 1;
}

int16_t sshr(const int16_t a, const size_t i){
    if(i == 0){
        return a;
    }
    return a >> i;
}

int16_t srsra(int16_t c, const int16_t a, const size_t i){
    return c + srshr(a, i);
}

int16_t ssra(int16_t c, const int16_t a, const size_t i){
    return c + sshr(a, i);
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
    // return (((int32_t)a * 315 + (1 << 18)) >> 19) & 0x1;
    return srshr(sqdmulh(a, 315), 4) & 0x1;
}

int16_t Barrett_quotient_4(int16_t a){
    // 16-bit suffices for D = 4
    // 315 = round(16 * 2^16 / q)
    // return (((int32_t)a * 315 + (1 << 15)) >> 16) & 0xf;
    return shadd(sqdmulh(a, 315), 1) & 0xf;
}

int16_t Barrett_quotient_5(int16_t a){
    // 15-bit suffices for D = 5
    // 315 = round(32 * 2^15 / q)
    // return (((int32_t)a * 315 + (1 << 14)) >> 15) & 0x1f;
    return sqrdmulh(a, 315) & 0x1f;
}

int16_t Barrett_quotient_10(int16_t a){
    // 22-bit suffices for D = 10
    // 1290167 = round(1024 * 2^22 / q)
    // beware that adding prior to shifting overflows (32-bit), we must shift, add, and then shift here.
    // return ( ((((int32_t)a * 1290167) >> 1) + (1 << 20)) >> 21) & 0x3ff;
    return srshr((mla(shadd(sqdmulh(a, -20553), 0), a, 20)), 6) & 0x3ff;
}

int16_t Barrett_quotient_11(int16_t a){
    // 21-bit suffices for D = 11
    // 1290167 = round(2048 * 2^21 / q)
    // beware that adding prior to shifting overflows (32-bit), we must shift, add, and then shift here.
    // return ( ((((int32_t)a * 1290167) >> 1) + (1 << 19)) >> 20) & 0x7ff;
    return srshr((mla(shadd(sqdmulh(a, -20553), 0), a, 20)), 5) & 0x7ff;
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
            t[k] = (uint16_t)Barrett_quotient_10(u);

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

static inline
int8x16_t shift_with_holes_b(int8x16_t a, const int8x16_t mask_lo, const int8x16_t mask_hi, const size_t shift_i){
    return (int8x16_t)vorrq_s16((int16x8_t)vandq_s8(a, mask_lo),
                                 vshrq_n_s16((int16x8_t)vandq_s8(a, mask_hi), shift_i));
}

static inline
int8x16_t shift_with_holes_b_lazy(int8x16_t a, const int8x16_t mask_hi, const size_t shift_i){
    return (int8x16_t)vorrq_s16((int16x8_t)a,
                                 vshrq_n_s16((int16x8_t)vandq_s8(a, mask_hi), shift_i));
}

static inline
int8x16_t shift_with_holes_b_very_lazy(int8x16_t a, const size_t shift_i){
    return (int8x16_t)vorrq_s16((int16x8_t)a,
                                 vshrq_n_s16((int16x8_t)a, shift_i));
}

static inline
int16x8_t shift_with_holes_h(int16x8_t a, const int16x8_t mask_lo, const int16x8_t mask_hi, const size_t shift_i){
    return (int16x8_t)vorrq_s32((int32x4_t)vandq_s16(a, mask_lo),
                                 vshrq_n_s32((int32x4_t)vandq_s16(a, mask_hi), shift_i));
}

static inline
int32x4_t shift_with_holes_s(int32x4_t a, const int32x4_t mask_lo, const int32x4_t mask_hi, const size_t shift_i){
    return (int32x4_t)vorrq_s64((int64x2_t)vandq_s32(a, mask_lo),
                                 vshrq_n_s64((int64x2_t)vandq_s32(a, mask_hi), shift_i));
}

void poly_compress1_neon(uint8_t r[32], const int16_t a[KYBER_N]){

    int16x8_t tvec[16];
    int16x8_t mask1 = vdupq_n_s16(0x1);
#if __ARM_FEATURE_DOTPROD
    int8x16_t dot_v1 = {1, 2, 4, 8, 1, 2, 4, 8, 1, 2, 4, 8, 1, 2, 4, 8};
    int8x16_t dot_v4 = {1, 0, 16, 0, 1, 0, 16, 0, 1, 0, 16, 0, 1, 0, 16, 0};
    int32x4_t zero_int32x4 = {0, 0, 0, 0};
#endif

    for(size_t i = 0; i < KYBER_N / 128; i++){

        for(size_t j = 0; j < 16; j++){
            tvec[j] = vld1q_s16(a + i * 128 + 8 * j);
            tvec[j] = vqdmulhq_n_s16(tvec[j], 315);
            tvec[j] = vrshrq_n_s16(tvec[j], 4);
            tvec[j] = vandq_s16(tvec[j], mask1);
        }

#if __ARM_FEATURE_DOTPROD

        for(size_t j = 0; j < 16; j += 2){
            tvec[j] = (int16x8_t)vmovn_high_s16(vmovn_s16(tvec[j]), tvec[j + 1]);
            tvec[j] = (int16x8_t)vdotq_s32(zero_int32x4, (int8x16_t)tvec[j], dot_v1);
        }

        for(size_t j = 0; j < 16; j += 4){
            tvec[j] = vmovn_high_s32(vmovn_s32((int32x4_t)tvec[j]), (int32x4_t)tvec[j + 2]);
            tvec[j] = (int16x8_t)vdotq_s32(zero_int32x4, (int8x16_t)tvec[j], dot_v4);
        }

        tvec[0] = (int16x8_t)vmovn_high_s32(vmovn_s32((int32x4_t)tvec[0]), (int32x4_t)tvec[4]);
        tvec[8] = (int16x8_t)vmovn_high_s32(vmovn_s32((int32x4_t)tvec[8]), (int32x4_t)tvec[12]);

        tvec[0] = (int16x8_t)vmovn_high_s16(vmovn_s16(tvec[0]), tvec[8]);

#else

        for(size_t j = 0; j < 16; j += 2){
            tvec[j] = (int16x8_t)vmovn_high_s16(vmovn_s16(tvec[j]), tvec[j + 1]);
            tvec[j] = (int16x8_t)shift_with_holes_b_very_lazy((int8x16_t)tvec[j], 7);
        }

        for(size_t j = 0; j < 16; j += 4){
            tvec[j] = (int16x8_t)vmovn_high_s16(vmovn_s16(tvec[j]), tvec[j + 2]);
            tvec[j] = (int16x8_t)shift_with_holes_b_very_lazy((int8x16_t)tvec[j], 6);
        }

        tvec[0] = (int16x8_t)vmovn_high_s16(vmovn_s16(tvec[0]), tvec[4]);
        tvec[8] = (int16x8_t)vmovn_high_s16(vmovn_s16(tvec[8]), tvec[12]);
        tvec[0] = (int16x8_t)shift_with_holes_b_very_lazy((int8x16_t)tvec[0], 4);
        tvec[8] = (int16x8_t)shift_with_holes_b_very_lazy((int8x16_t)tvec[8], 4);

        tvec[0] = (int16x8_t)vmovn_high_s16(vmovn_s16(tvec[0]), tvec[8]);

#endif

        vst1q_u8(r, (uint8x16_t)tvec[0]);
        r += 16;

    }

}

void poly_compress4_neon(uint8_t r[128], const int16_t a[KYBER_N]){

    int16x8_t tvec[4];
    int16x8_t mask4 = vdupq_n_s16(0xf);
    int16x8_t one = vdupq_n_s16(1);
#if __ARM_FEATURE_DOTPROD
    int8x16_t dot_v = {1, 0, 16, 0, 1, 0, 16, 0, 1, 0, 16, 0, 1, 0, 16, 0};
    int32x4_t zero_int32x4 = {0, 0, 0, 0};
#endif

    for(size_t i = 0; i < KYBER_N / 32; i++) {

        tvec[0] = vld1q_s16(a + i * 32 + 8 * 0);
        tvec[1] = vld1q_s16(a + i * 32 + 8 * 1);
        tvec[2] = vld1q_s16(a + i * 32 + 8 * 2);
        tvec[3] = vld1q_s16(a + i * 32 + 8 * 3);

        tvec[0] = vqdmulhq_n_s16(tvec[0], 315);
        tvec[1] = vqdmulhq_n_s16(tvec[1], 315);
        tvec[2] = vqdmulhq_n_s16(tvec[2], 315);
        tvec[3] = vqdmulhq_n_s16(tvec[3], 315);

        tvec[0] = vhaddq_s16(tvec[0], one);
        tvec[1] = vhaddq_s16(tvec[1], one);
        tvec[2] = vhaddq_s16(tvec[2], one);
        tvec[3] = vhaddq_s16(tvec[3], one);

        tvec[0] = vandq_s16(tvec[0], mask4);
        tvec[1] = vandq_s16(tvec[1], mask4);
        tvec[2] = vandq_s16(tvec[2], mask4);
        tvec[3] = vandq_s16(tvec[3], mask4);

#if __ARM_FEATURE_DOTPROD

        tvec[0] = (int16x8_t)vdotq_s32(zero_int32x4, (int8x16_t)tvec[0], dot_v);
        tvec[1] = (int16x8_t)vdotq_s32(zero_int32x4, (int8x16_t)tvec[1], dot_v);
        tvec[2] = (int16x8_t)vdotq_s32(zero_int32x4, (int8x16_t)tvec[2], dot_v);
        tvec[3] = (int16x8_t)vdotq_s32(zero_int32x4, (int8x16_t)tvec[3], dot_v);

        tvec[0] = vmovn_high_s32(vmovn_s32((int32x4_t)tvec[0]), (int32x4_t)tvec[1]);
        tvec[2] = vmovn_high_s32(vmovn_s32((int32x4_t)tvec[2]), (int32x4_t)tvec[3]);

        tvec[0] = (int16x8_t)vmovn_high_s16(vmovn_s16(tvec[0]), tvec[2]);

#else

        tvec[0] = (int16x8_t)vmovn_high_s16(vmovn_s16(tvec[0]), tvec[1]);
        tvec[2] = (int16x8_t)vmovn_high_s16(vmovn_s16(tvec[2]), tvec[3]);

        tvec[0] = (int16x8_t)shift_with_holes_b_very_lazy((int8x16_t)tvec[0], 4);
        tvec[2] = (int16x8_t)shift_with_holes_b_very_lazy((int8x16_t)tvec[2], 4);

        tvec[0] = (int16x8_t)vmovn_high_s16(vmovn_s16(tvec[0]), tvec[2]);

#endif

        vst1q_s16((int16_t*)r, tvec[0]);

        r += 16;
    }

}

void poly_compress5_neon(uint8_t r[160], const int16_t a[KYBER_N]){

    int16x8_t tvec[2];
    int16x8_t mask5 = vdupq_n_s16(0x1f);
    int16x8_t mask_h_lo = {0x3ff, 0, 0x3ff, 0, 0x3ff, 0, 0x3ff, 0};
    int16x8_t mask_h_hi = {0, 0x3ff, 0, 0x3ff, 0, 0x3ff, 0, 0x3ff};
    int32x4_t mask_w_lo = {0xfffff, 0, 0xfffff, 0};
    int32x4_t mask_w_hi = {0, 0xfffff, 0, 0xfffff};
#if __ARM_FEATURE_DOTPROD
    int8x16_t dot_v = {1, 0, 32, 0, 1, 0, 32, 0, 1, 0, 32, 0, 1, 0, 32, 0};
    int32x4_t zero_int32x4 = {0, 0, 0, 0};
#else
    int8x16_t mask_b_lo = {0x1f, 0, 0x1f, 0, 0x1f, 0, 0x1f, 0, 0x1f, 0, 0x1f, 0, 0x1f, 0, 0x1f, 0};
    int8x16_t mask_b_hi = {0, 0x1f, 0, 0x1f, 0, 0x1f, 0, 0x1f, 0, 0x1f, 0, 0x1f, 0, 0x1f, 0, 0x1f};
#endif

    uint16_t t[2][8];
    for(size_t i = 0; i < KYBER_N / 16; i++) {

        tvec[0] = vld1q_s16(a + 16 * i + 8 * 0);
        tvec[1] = vld1q_s16(a + 16 * i + 8 * 1);
        tvec[0] = vqrdmulhq_n_s16(tvec[0], 315);
        tvec[1] = vqrdmulhq_n_s16(tvec[1], 315);
        tvec[0] = vandq_s16(tvec[0], mask5);
        tvec[1] = vandq_s16(tvec[1], mask5);

#if __ARM_FEATURE_DOTPROD

        tvec[0] = (int16x8_t)vdotq_s32(zero_int32x4, (int8x16_t)tvec[0], dot_v);
        tvec[1] = (int16x8_t)vdotq_s32(zero_int32x4, (int8x16_t)tvec[1], dot_v);

        tvec[0] = vmovn_high_s32(vmovn_s32((int32x4_t)tvec[0]), (int32x4_t)tvec[1]);

        tvec[0] = shift_with_holes_h(tvec[0], mask_h_lo, mask_h_hi, 6);
        tvec[0] = (int16x8_t)shift_with_holes_s((int32x4_t)tvec[0], mask_w_lo, mask_w_hi, 12);

#else

        tvec[0] = (int16x8_t)vmovn_high_s16(vmovn_s16(tvec[0]), tvec[1]);

        tvec[0] = (int16x8_t)shift_with_holes_b((int8x16_t)tvec[0], mask_b_lo, mask_b_hi, 3);
        tvec[0] = shift_with_holes_h(tvec[0], mask_h_lo, mask_h_hi, 6);
        tvec[0] = (int16x8_t)shift_with_holes_s((int32x4_t)tvec[0], mask_w_lo, mask_w_hi, 12);

#endif

        vst1q_s16((int16_t*)t[0], tvec[0]);

        r[0] = t[0][0];
        r[1] = t[0][0] >> 8;
        r[2] = t[0][1];
        r[3] = t[0][1] >> 8;
        r[4] = t[0][2];
        r[5] = t[0][4];
        r[6] = t[0][4] >> 8;
        r[7] = t[0][5];
        r[8] = t[0][5] >> 8;
        r[9] = t[0][6];

        r += 10;
    }

}

void poly_compress10_neon(uint8_t r[320], const int16_t a[KYBER_N]){

    uint16_t t[8];
    int16x8_t avec, tvec;
    int16x8_t mask10 = vdupq_n_s16(0x3ff);
    int16x8_t zero = vdupq_n_s16(0);

    int16x8_t mask_h_lo = {0x3ff, 0, 0x3ff, 0, 0x3ff, 0, 0x3ff, 0};
    int16x8_t mask_h_hi = {0, 0x3ff, 0, 0x3ff, 0, 0x3ff, 0, 0x3ff};
    int32x4_t mask_w_lo = {0xfffff, 0, 0xfffff, 0};
    int32x4_t mask_w_hi = {0, 0xfffff, 0, 0xfffff};

    for(size_t i = 0; i < KYBER_N / 8; i++){
        avec = vld1q_s16(a + 8 * i);

        tvec = vqdmulhq_n_s16(avec, -20553);
        tvec = vhaddq_s16(tvec, zero);
        tvec = vmlaq_n_s16(tvec, avec, 20);
        tvec = vrshrq_n_s16(tvec, 6);
        tvec = vandq_s16(tvec, mask10);

        tvec = shift_with_holes_h(tvec, mask_h_lo, mask_h_hi, 6);
        tvec = (int16x8_t)shift_with_holes_s((int32x4_t)tvec, mask_w_lo, mask_w_hi, 12);

        vst1q_s16((int16_t*)t, tvec);

        r[0] = t[0];
        r[1] = t[0] >> 8;
        r[2] = t[1];
        r[3] = t[1] >> 8;
        r[4] = t[2];
        r[5] = t[4];
        r[6] = t[4] >> 8;
        r[7] = t[5];
        r[8] = t[5] >> 8;
        r[9] = t[6];

        r += 10;
    }

}

void poly_compress11_neon(uint8_t r[352], const int16_t a[KYBER_N]){

    uint64_t lo64, hi64;
    __attribute__((aligned(16))) uint16_t t[8];
    int16x8_t avec, tvec;
    int16x8_t mask11 = vdupq_n_s16(0x7ff);
    int16x8_t zero = vdupq_n_s16(0);

    int16x8_t mask_h_lo = {0x7ff, 0, 0x7ff, 0, 0x7ff, 0, 0x7ff, 0};
    int16x8_t mask_h_hi = {0, 0x7ff, 0, 0x7ff, 0, 0x7ff, 0, 0x7ff};
    int32x4_t mask_w_lo = {0x3fffff, 0, 0x3fffff, 0};
    int32x4_t mask_w_hi = {0, 0x3fffff, 0, 0x3fffff};

    for(size_t i = 0; i < KYBER_N / 8; i++) {

        avec = vld1q_s16(a + 8 * i);

        tvec = vqdmulhq_n_s16(avec, -20553);
        tvec = vhaddq_s16(tvec, zero);
        tvec = vmlaq_n_s16(tvec, avec, 20);
        tvec = vrshrq_n_s16(tvec, 5);
        tvec = vandq_s16(tvec, mask11);

        tvec = shift_with_holes_h(tvec, mask_h_lo, mask_h_hi, 5);
        tvec = (int16x8_t)shift_with_holes_s((int32x4_t)tvec, mask_w_lo, mask_w_hi, 10);

        vst1q_s16((int16_t*)t, tvec);

        r[0] = t[0];
        r[1] = t[0] >> 8;
        r[2] = t[1];
        r[3] = t[1] >> 8;
        r[4] = t[2];

        memcpy(&lo64, t + 0, 8);
        memcpy(&hi64, t + 4, 8);
        hi64 = (hi64 << 4) | (lo64 >> 40);

        r[5] = hi64;
        r[6] = hi64 >> 8;
        r[7] = hi64 >> 16;
        r[8] = hi64 >> 24;
        r[9] = hi64 >> 32;
        r[10] = hi64 >> 40;

        r += 11;
    }

}

int main(void){

    int16_t a[KYBER_N];
    uint8_t ref[352], res[352];

    for(size_t i = 0; i < KYBER_N; i++){
        a[i] = rand() % KYBER_Q;
        a[i] -= KYBER_Q / 2;
    }

    poly_compress1(ref, a);
    poly_compress1_neon(res, a);

    assert(memcmp(ref, res, 32) == 0);

    poly_compress4(ref, a);
    poly_compress4_neon(res, a);

    assert(memcmp(ref, res, 128) == 0);

    poly_compress5(ref, a);
    poly_compress5_neon(res, a);

    assert(memcmp(ref, res, 160) == 0);

    poly_compress10(ref, a);
    poly_compress10_neon(res, a);

    assert(memcmp(ref, res, 320) == 0);

    poly_compress11(ref, a);
    poly_compress11_neon(res, a);

    assert(memcmp(ref, res, 352) == 0);

    printf("Test finished!\n");

}

