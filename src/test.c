
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <assert.h>

// 1, 4, 5, 10, 11

#define KYBER_Q 3329

int16_t Barrett_floor_reduce(int16_t a){
    return a - (int16_t)(((int32_t)a * 20159) >> 26) * KYBER_Q;
}

int16_t Barrett_round_reduce(int16_t a){
    return a - (int16_t)(((int32_t)a * 20159 + (1 << 25)) >> 26) * KYBER_Q;
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
    return ((int16_t)(((int32_t)a * 315 + (1 << 18)) >> 19)) & 0x1;
}

int16_t Barrett_quotient_2(int16_t a){
    // 18-bit suffices for D = 2
    // 315 = round(4 * 2^18 / q)
    return ((int16_t)(((int32_t)a * 315 + (1 << 17)) >> 18)) & 0x3;
}

int16_t Barrett_quotient_3(int16_t a){
    // 17-bit suffices for D = 3
    // 315 = round(8 * 2^17 / q)
    return ((int16_t)(((int32_t)a * 315 + (1 << 16)) >> 17)) & 0x7;
}

int16_t Barrett_quotient_4(int16_t a){
    // 16-bit suffices for D = 4
    // 315 = round(16 * 2^16 / q)
    return ((int16_t)(((int32_t)a * 315 + (1 << 15)) >> 16)) & 0xf;
}

int16_t Barrett_quotient_5(int16_t a){
    // 15-bit suffices for D = 5
    // 315 = round(32 * 2^15 / q)
    return ((int16_t)(((int32_t)a * 315 + (1 << 14)) >> 15)) & 0x1f;
}

int16_t Barrett_quotient_6(int16_t a){
    // 14-bit suffices for D = 6
    // 315 = round(64 * 2^14 / q)
    return ((int16_t)(((int32_t)a * 315 + (1 << 13)) >> 14)) & 0x3f;
}

int16_t Barrett_quotient_7(int16_t a){
    // 13-bit suffices for D = 7
    // 315 = round(128 * 2^13 / q)
    return ((int16_t)(((int32_t)a * 315 + (1 << 12)) >> 13)) & 0x7f;
}

int16_t Barrett_quotient_8(int16_t a){
    // 12-bit suffices for D = 8
    // 315 = round(256 * 2^12 / q)
    return ((int16_t)(((int32_t)a * 161271 + (1 << 20)) >> 21)) & 0xff;
}

int16_t Barrett_quotient_9(int16_t a){
    // this doesn't work
    // return ((int16_t)(((int32_t)a * 315 + (1 << 10)) >> 11)) & 0x1ff;
    // 20-bit suffices for D = 9
    // 161271 = round(512 * 2^20 / q)
    return ((int16_t)(((int32_t)a * 1290167 + (1 << 22)) >> 23)) & 0x1ff;
}

int16_t Barrett_quotient_10(int16_t a){
    // this doesn't work
    // return ((int16_t)(((int32_t)a * 161271 + (1 << 18)) >> 19)) & 0x3ff;
    // 22-bit suffices for D = 10
    // 1290167 = round(1024 * 2^22 / q)
    return ((int16_t)(((int32_t)a * 1290167 + (1 << 21)) >> 22)) & 0x3ff;
}

int16_t Barrett_quotient_11(int16_t a){
    // 21-bit suffices for D = 11
    // 1290167 = round(2048 * 2^21 / q)
    return ((int16_t)(((int32_t)a * 1290167 + (1 << 20)) >> 21)) & 0x7ff;
}

int main(void){

    int16_t ubound, t;

    ubound = 0;
    for(int16_t i = -3328; i <= 3229; i++){
        t = Barrett_floor_reduce(i);
        if(t > ubound){
            ubound = t;
        }
        if((t < 0) && (t < -ubound)){
            ubound = -t;
        }
    }
    printf("ubound of floor reduction: %4d\n", ubound);

    ubound = 0;
    for(int16_t i = -32767; i < 32767; i++){
        t = Barrett_round_reduce(i);
        if(t > ubound){
            ubound = t;
        }
        if((t < 0) && (t < -ubound)){
            ubound = -t;
        }
    }
    printf("ubound of round reduction: %4d\n", ubound);

    for(int16_t i = -1664; i <= 3329; i++){
        assert(compress_D(i, 1) == Barrett_quotient_1(i));
    }

    for(int16_t i = -1664; i <= 3329; i++){
        assert(compress_D(i, 2) == Barrett_quotient_2(i));
    }

    for(int16_t i = -1664; i <= 3329; i++){
        assert(compress_D(i, 3) == Barrett_quotient_3(i));
    }

    for(int16_t i = -1664; i <= 3329; i++){
        assert(compress_D(i, 4) == Barrett_quotient_4(i));
    }

    for(int16_t i = -1664; i <= 3329; i++){
        assert(compress_D(i, 5) == Barrett_quotient_5(i));
    }

    for(int16_t i = -1664; i <= 3329; i++){
        assert(compress_D(i, 6) == Barrett_quotient_6(i));
    }

    for(int16_t i = -1664; i <= 3329; i++){
        assert(compress_D(i, 7) == Barrett_quotient_7(i));
    }

    for(int16_t i = -1664; i <= 3329; i++){
        assert(compress_D(i, 8) == Barrett_quotient_8(i));
    }

    for(int16_t i = -1664; i <= 3329; i++){
        assert(compress_D(i, 9) == Barrett_quotient_9(i));
    }

    for(int16_t i = -1664; i <= 1664; i++){
        assert(compress_D(i, 10) == Barrett_quotient_10(i));
    }

    for(int16_t i = -1664; i <= 1664; i++){
        assert(compress_D(i, 11) == Barrett_quotient_11(i));
    }

    printf("Test finished!\n");

}

