
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <assert.h>

// we can use ubfx here

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
    // 2580335 = round(2 * 2^32 / q)
    return (smlawx(1 << 15, a, 2580335) >> 16) & 0x1;
    // return smmulr((int32_t)a, 2580335) & 0x1;
}

int16_t Barrett_quotient_4(int16_t a){
    // 16-bit suffices for D = 4
    // 315 = round(16 * 2^16 / q)
    // return ((int16_t)(((int32_t)a * 315 + (1 << 15)) >> 16)) & 0xf;
    // 20642679 = round(16 * 2^32)
    return (smlawx(1 << 15, a, 20642679) >> 16) & 0xf;
    // return smmulr((int32_t)a, 20642679) & 0xf;

}

int16_t Barrett_quotient_5(int16_t a){
    // 15-bit suffices for D = 5
    // 315 = round(32 * 2^15 / q)
    // return ((int16_t)(((int32_t)a * 315 + (1 << 14)) >> 15)) & 0x1f;
    // 41285357 = round(32 * 2^32 / q)
    return (smlawx(1 << 15, a, 41285357) >> 16) & 0x1f;
    // return smmulr((int32_t)a, 41285357) & 0x1f;
}

int16_t Barrett_quotient_10(int16_t a){
    // this doesn't work
    // return ((int16_t)(((int32_t)a * 161271 + (1 << 18)) >> 19)) & 0x3ff;
    // 22-bit suffices for D = 10
    // 1321131424 = round(1024 * 2^32 / q)
    return (smlawx(1 << 15, a, 1321131424) >> 16) & 0x3ff;
    // return smmulr((int32_t)a, 1321131424) & 0x3ff;
}

int16_t Barrett_quotient_11(int16_t a){
    // 21-bit suffices for D = 11
    // 1290167 = round(2048 * 2^21 / q)
    // beware that adding prior to shifting overflows (32-bit), we must shift, add, and then shift here.
    // return ((int16_t)(((int32_t)a * 1290167 + (1 << 20)) >> 21)) & 0x7ff;
    return ubfx(smlawx(1 << 4, a, 1290167), 5, 11);
    // return (smlawx(1 << 4, a, 1290167) >> 5) & 0x7ff;
}

int main(void){

    for(int16_t i = -1664; i <= 1664; i++){
        assert(compress_D(i, 1) == Barrett_quotient_1(i));
    }

    for(int16_t i = -1664; i <= 1664; i++){
        assert(compress_D(i, 4) == Barrett_quotient_4(i));
    }

    for(int16_t i = -1664; i <= 1664; i++){
        assert(compress_D(i, 5) == Barrett_quotient_5(i));
    }

    for(int16_t i = -1664; i <= 1664; i++){
        assert(compress_D(i, 10) == Barrett_quotient_10(i));
    }

    for(int16_t i = -1664; i <= 1664; i++){
        assert(compress_D(i, 11) == Barrett_quotient_11(i));
    }

    printf("Test finished!\n");

}

