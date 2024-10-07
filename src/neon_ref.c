
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <assert.h>

// 1, 4, 5, 10, 11

#define KYBER_Q 3329

int16_t shadd(const int16_t a, const int16_t b){
    return (a + b) >> 1;
}

int16_t mul(const int16_t a, const int16_t b){
    return a * b;
}

int16_t mla(const int16_t c, const int16_t a, const int16_t b){
    return c + a * b;
}

int16_t mls(const int16_t c, const int16_t a, const int16_t b){
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

