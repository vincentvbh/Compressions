
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <assert.h>

// 1, 4, 5, 10, 11

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
// Quotient

int16_t Barrett_quotient_1(int16_t a){
    // 32-bit suffices for D = 1
    // 2580335 = round(2 * 2^32 / q)
    // return (smlawx(1 << 15, a, 2580335) >> 16);
    return smmulr((int32_t)a, 2580335);
}

int16_t Barrett_quotient_4(int16_t a){
    // 32-bit suffices for D = 4
    // 20642679 = round(16 * 2^32)
    // return (smlawx(1 << 15, a, 20642679) >> 16);
    return smmulr((int32_t)a, 20642679);

}

int16_t Barrett_quotient_5(int16_t a){
    // 32-bit suffices for D = 5
    // 41285357 = round(32 * 2^32 / q)
    // return (smlawx(1 << 15, a, 41285357) >> 16);
    return smmulr((int32_t)a, 41285357);
}

int16_t Barrett_quotient_10(int16_t a){
    // 32-bit suffices for D = 10
    // 1321131424 = round(1024 * 2^32 / q)
    // return (smlawx(1 << 15, a, 1321131424) >> 16);
    return smmulr((int32_t)a, 1321131424);
}

int16_t Barrett_quotient_11(int16_t a){
    // 21-bit suffices for D = 11
    // 1290167 = round(2048 * 2^21 / q)
    return (smlawx(1 << 4, a, 1290167) >> 5);
}

// ================
// Quotient with large input

int16_t Barrett_quotient_large_11(int16_t a){
    // 24-bit suffices for D = 11
    // 10321339 = round(2048 * 2^24 / q)
    return (smlawx(1 << 7, a, 10321339) >> 8);
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

int main(void){

// ================
// Quotient for inputs in [-1664, 1664]

    for(int16_t i = -1664; i <= 1664; i++){
        assert(quotient_D_sign(i, 1) == Barrett_quotient_1(i));
    }

    for(int16_t i = -1664; i <= 1664; i++){
        assert(quotient_D_sign(i, 4) == Barrett_quotient_4(i));
    }

    for(int16_t i = -1664; i <= 1664; i++){
        assert(quotient_D_sign(i, 5) == Barrett_quotient_5(i));
    }

    for(int16_t i = -1664; i <= 1664; i++){
        assert(quotient_D_sign(i, 10) == Barrett_quotient_10(i));
    }

    for(int16_t i = -1664; i <= 1664; i++){
        assert(quotient_D_sign(i, 11) == Barrett_quotient_11(i));
    }

// ================
// Quotient for inputs in [-3328, 3328]

    for(int16_t i = -3328; i <= 3328; i++){
        assert(quotient_D_sign(i, 1) == Barrett_quotient_1(i));
    }

    for(int16_t i = -3328; i <= 3328; i++){
        assert(quotient_D_sign(i, 4) == Barrett_quotient_4(i));
    }

    for(int16_t i = -3328; i <= 3328; i++){
        assert(quotient_D_sign(i, 5) == Barrett_quotient_5(i));
    }

    for(int16_t i = -3328; i <= 3328; i++){
        assert(quotient_D_sign(i, 10) == Barrett_quotient_10(i));
    }

    for(int16_t i = -3328; i <= 3328; i++){
        assert(quotient_D_sign(i, 11) == Barrett_quotient_large_11(i));
    }

// ================
// Compression

    for(int16_t i = -1664; i <= 1664; i++){
        assert(compress_D(i, 1) == Barrett_compress_1(i));
    }

    for(int16_t i = -1664; i <= 1664; i++){
        assert(compress_D(i, 4) == Barrett_compress_4(i));
    }

    for(int16_t i = -1664; i <= 1664; i++){
        assert(compress_D(i, 5) == Barrett_compress_5(i));
    }

    for(int16_t i = -1664; i <= 1664; i++){
        assert(compress_D(i, 10) == Barrett_compress_10(i));
    }

    for(int16_t i = -1664; i <= 1664; i++){
        assert(compress_D(i, 11) == Barrett_compress_11(i));
    }

    printf("Test finished!\n");

}

