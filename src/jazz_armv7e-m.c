
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

// 1, 4, 5, 10, 11

#define KYBER_N 256
#define KYBER_Q 3329

int32_t mul(const int32_t a, const int32_t b){
    return a * b;
}

int32_t mla(const int32_t c, const int32_t a, const int32_t b){
    return c + a * b;
}

int32_t mls(const int32_t c, const int32_t a, const int32_t b){
    return c - a * b;
}

int32_t smmul(const int32_t a, const int32_t b){
    return ((int64_t)a * b) >> 32;
}

int32_t smmulr(const int32_t a, const int32_t b){
    return ((int64_t)a * b + (1LL << 31)) >> 32;
}

int32_t smlawx(const int32_t c, const int16_t a, const int32_t b){
    return (((int64_t)a * b) >> 16) + c;
}

int32_t ubfx(const int32_t a, size_t pos, size_t width){
    if(width < 1){
        return a;
    }
    if(pos >= 32){
        return 0;
    }
    if(pos + width > 32){
        width = 32 - pos;
    }
    if(width == 32){
        assert(pos == 0);
        return a;
    }
    return (a >> pos) & ((1 << width) - 1);
}

int16_t quotient_D_sign(int16_t a, const size_t D){
    if(a >= 0){
        // round in the non-negative case.
        return (int16_t)(( ( ((int32_t)a) << D) + (KYBER_Q / 2) ) / KYBER_Q);
    }else{
        // In C, division rounds the negative results toward zero instead of rounding-half-up.
        // Observe that for a positive real number r, round(-r) = -round(r) expect for r an half-integer for round
        // the rounding-half-up function.
        // Fortunately, since KYBER_Q is odd, a * 2^D / KYBER_Q is never a half-integer.
        // To sum up, we negate the input, round it as in the non-negative case, and return the negation of the result.
        return -(int16_t)(( ( ((int32_t)-a) << D) + (KYBER_Q / 2) ) / KYBER_Q);
    }
}

int16_t mulmod(int16_t a, const size_t D){
    return ( ((int32_t)a << D) - (int32_t)quotient_D_sign(a, D) * KYBER_Q);
}

int16_t compress_D(int16_t a, const size_t D){
    if(a < 0){
        a += KYBER_Q;
    }
    return (int16_t)(( ( ((int32_t)a) << D) + (KYBER_Q / 2) ) / KYBER_Q) & ((1 << D) - 1);
}

// ================
// Compression

int16_t Barrett_compress_1(int16_t a){
    // 32-bit suffices for D = 1
    // 2580335 = round(2 * 2^32 / q)
    // return (smlawx(1 << 15, a, 2580335) >> 16) & 0x1;
    return smmulr((int32_t)a, 2580335) & 0x1;
}

int16_t Barrett_compress_4(int16_t a){
    // 32-bit suffices for D = 4
    // 20642679 = round(16 * 2^32)
    // return (smlawx(1 << 15, a, 20642679) >> 16) & 0xf;
    return smmulr((int32_t)a, 20642679) & 0xf;

}

int16_t Barrett_compress_5(int16_t a){
    // 32-bit suffices for D = 5
    // 41285357 = round(32 * 2^32 / q)
    // return (smlawx(1 << 15, a, 41285357) >> 16) & 0x1f;
    return smmulr((int32_t)a, 41285357) & 0x1f;
}

int16_t Barrett_compress_10(int16_t a){
    // 32-bit suffices for D = 10
    // 1321131424 = round(1024 * 2^32 / q)
    // return (smlawx(1 << 15, a, 1321131424) >> 16) & 0x3f;
    return smmulr((int32_t)a, 1321131424) & 0x3ff;
}

int16_t Barrett_compress_11(int16_t a){
    // 21-bit suffices for D = 11
    // 1290167 = round(2048 * 2^21 / q)
    return (smlawx(1 << 4, a, 1290167) >> 5) & 0x7ff;
}



void poly_compress1(uint8_t r[32], const int16_t a[KYBER_N]){

    unsigned int i,j;
    int16_t u;

    for(i=0;i<KYBER_N/8;i++) {
        r[i] = 0;
        for(j=0;j<8;j++) {
            u = a[8*i+j];
            u = Barrett_compress_1(u);

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
            t[k] = Barrett_compress_4(u);

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
            t[k] = Barrett_compress_5(u);

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
            t[k] = Barrett_compress_10(u);

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
            t[k] = Barrett_compress_11(u);

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

extern void __jazz_poly_compress1(uint8_t *des, int16_t *src);
extern void __jazz_poly_compress4(uint8_t *des, int16_t *src);
extern void __jazz_poly_compress10(uint8_t *des, int16_t *src);

int main(void){

    int16_t a[KYBER_N];
    uint8_t ref[352], res[352];

    for(size_t i = 0; i < KYBER_N; i++){
        a[i] = rand() % KYBER_Q;
        a[i] -= KYBER_Q / 2;
    }

// ================
// Polynomial compression

    poly_compress1(ref, a);
    __jazz_poly_compress1(res, a);
    assert(memcmp(ref, res, 32) == 0);

    poly_compress4(ref, a);
    __jazz_poly_compress4(res, a);
    assert(memcmp(ref, res, 128) == 0);

    poly_compress10(ref, a);
    __jazz_poly_compress10(res, a);
    assert(memcmp(ref, res, 320) == 0);


    printf("Polynomial compression finished!\n");

}

