
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <assert.h>

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
    // return (((int32_t)a * 315 + (1 << 18)) >> 19) & 0x1;
    return pmulhrsw(pmulhw(a, 315), (1 << 12)) & 0x1;
    // return psraw(pmulhw(a, 315) + (1 << 2), 3) & 0x1;
}

int16_t Barrett_quotient_4(int16_t a){
    // 16-bit suffices for D = 4
    // 315 = round(16 * 2^16 / q)
    // return (((int32_t)a * 315 + (1 << 15)) >> 16) & 0xf;
    // return (((int32_t)a * 630 + (1 << 16)) >> 17) & 0xf;
    return pmulhrsw(pmulhw(a, 630), (1 << 14)) & 0xf;
    // return psraw(pmulhw(a, 630) + (1 << 0), 1) & 0xf;
}

int16_t Barrett_quotient_5(int16_t a){
    // 15-bit suffices for D = 5
    // 315 = round(32 * 2^15 / q)
    // return (((int32_t)a * 315 + (1 << 14)) >> 15) & 0x1f;
    return pmulhrsw(a, 315) & 0x1f;
}

// 1290167 = -20553 + 20 * 2^16

int16_t Barrett_quotient_10(int16_t a){
    // 22-bit suffices for D = 10
    // 1290167 = round(1024 * 2^22 / q)
    // beware that adding prior to shifting overflows (32-bit), we must shift, add, and then shift here.
    // return ( ((((int32_t)a * 1290167) >> 1) + (1 << 20)) >> 21) & 0x3ff;
    return pmulhrsw(pmulhw(a, -20553) + pmullw(a, 20), (1 << 9)) & 0x3ff;
}

int16_t Barrett_quotient_11(int16_t a){
    // 21-bit suffices for D = 11
    // 1290167 = round(2048 * 2^21 / q)
    // beware that adding prior to shifting overflows (32-bit), we must shift, add, and then shift here.
    // return ( ((((int32_t)a * 1290167) >> 1) + (1 << 19)) >> 20) & 0x7ff;
    return pmulhrsw(pmulhw(a, -20553) + pmullw(a, 20), (1 << 10)) & 0x7ff;
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

