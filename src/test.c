
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

int16_t Barrett_compress_1(int16_t a){
    // 19-bit suffices for D = 1
    // 315 = round(2 * 2^19 / q)
    return (((int32_t)a * 315 + (1 << 18)) >> 19) & 0x1;
}

int16_t Barrett_compress_2(int16_t a){
    // 18-bit suffices for D = 2
    // 315 = round(4 * 2^18 / q)
    return (((int32_t)a * 315 + (1 << 17)) >> 18) & 0x3;
}

int16_t Barrett_compress_3(int16_t a){
    // 17-bit suffices for D = 3
    // 315 = round(8 * 2^17 / q)
    return (((int32_t)a * 315 + (1 << 16)) >> 17) & 0x7;
}

int16_t Barrett_compress_4(int16_t a){
    // 16-bit suffices for D = 4
    // 315 = round(16 * 2^16 / q)
    return (((int32_t)a * 315 + (1 << 15)) >> 16) & 0xf;
}

int16_t Barrett_compress_5(int16_t a){
    // 15-bit suffices for D = 5
    // 315 = round(32 * 2^15 / q)
    return (((int32_t)a * 315 + (1 << 14)) >> 15) & 0x1f;
}

int16_t Barrett_compress_6(int16_t a){
    // 14-bit suffices for D = 6
    // 315 = round(64 * 2^14 / q)
    return (((int32_t)a * 315 + (1 << 13)) >> 14) & 0x3f;
}

int16_t Barrett_compress_7(int16_t a){
    // 13-bit suffices for D = 7
    // 315 = round(128 * 2^13 / q)
    return (((int32_t)a * 315 + (1 << 12)) >> 13) & 0x7f;
}

int16_t Barrett_compress_8(int16_t a){
    // 12-bit suffices for D = 8
    // 315 = round(256 * 2^12 / q)
    return (((int32_t)a * 161271 + (1 << 20)) >> 21) & 0xff;
}

int16_t Barrett_compress_9(int16_t a){
    // 20-bit suffices for D = 9
    // 161271 = round(512 * 2^20 / q)
    return (((int32_t)a * 161271 + (1 << 19)) >> 20) & 0x1ff;
}

int16_t Barrett_compress_10(int16_t a){
    // 22-bit suffices for D = 10
    // 1290167 = round(1024 * 2^22 / q)
    // beware that adding prior to shifting overflows (32-bit), we must shift, add, and then shift here.
    return ( ((((int32_t)a * 1290167) >> 1) + \
                    (1 << 20)) >> 21) & 0x3ff;
}

int16_t Barrett_compress_11(int16_t a){
    // 21-bit suffices for D = 11
    // 1290167 = round(2048 * 2^21 / q)
    // beware that adding prior to shifting overflows (32-bit), we must shift, add, and then shift here.
    return ( ((((int32_t)a * 1290167) >> 1) + \
                    (1 << 19)) >> 20) & 0x7ff;
}

int16_t Barrett_quotient_1(int16_t a){
    // 19-bit suffices for D = 1
    // 315 = round(2 * 2^19 / q)
    return (((int32_t)a * 315 + (1 << 18)) >> 19);
}

int16_t Barrett_quotient_2(int16_t a){
    // 18-bit suffices for D = 2
    // 315 = round(4 * 2^18 / q)
    return (((int32_t)a * 315 + (1 << 17)) >> 18);
}

int16_t Barrett_quotient_3(int16_t a){
    // 17-bit suffices for D = 3
    // 315 = round(8 * 2^17 / q)
    return (((int32_t)a * 315 + (1 << 16)) >> 17);
}

int16_t Barrett_quotient_4(int16_t a){
    // 16-bit suffices for D = 4
    // 315 = round(16 * 2^16 / q)
    return (((int32_t)a * 315 + (1 << 15)) >> 16);
}

int16_t Barrett_quotient_5(int16_t a){
    // 15-bit suffices for D = 5
    // 315 = round(32 * 2^15 / q)
    return (((int32_t)a * 315 + (1 << 14)) >> 15);
}

int16_t Barrett_quotient_6(int16_t a){
    // 14-bit suffices for D = 6
    // 315 = round(64 * 2^14 / q)
    return (((int32_t)a * 315 + (1 << 13)) >> 14);
}

int16_t Barrett_quotient_7(int16_t a){
    // 13-bit suffices for D = 7
    // 315 = round(128 * 2^13 / q)
    return (((int32_t)a * 315 + (1 << 12)) >> 13);
}

int16_t Barrett_quotient_8(int16_t a){
    // 12-bit suffices for D = 8
    // 315 = round(256 * 2^12 / q)
    return (((int32_t)a * 161271 + (1 << 20)) >> 21);
}

int16_t Barrett_quotient_9(int16_t a){
    // 20-bit suffices for D = 9
    // 161271 = round(512 * 2^20 / q)
    return (((int32_t)a * 161271 + (1 << 19)) >> 20);
}

int16_t Barrett_quotient_10(int16_t a){
    // 22-bit suffices for D = 10
    // 1290167 = round(1024 * 2^22 / q)
    // beware that adding prior to shifting overflows (32-bit), we must shift, add, and then shift here.
    return ( ((((int32_t)a * 1290167) >> 1) + \
                    (1 << 20)) >> 21);
}

int16_t Barrett_quotient_11(int16_t a){
    // 21-bit suffices for D = 11
    // 1290167 = round(2048 * 2^21 / q)
    // beware that adding prior to shifting overflows (32-bit), we must shift, add, and then shift here.
    return ( ((((int32_t)a * 1290167) >> 1) + \
                    (1 << 19)) >> 20);
}

int main(void){

    int16_t ubound, t;

    ubound = 0;
    for(int32_t i = -32768; i < 32768; i++){
        t = Barrett_floor_reduce(i);
        assert(t >= 0);
        if(t > ubound){
            ubound = t;
        }
        if((t < 0) && (t < -ubound)){
            ubound = -t;
        }
    }
    printf("ubound of floor reduction: %4d\n", ubound);

    ubound = 0;
    for(int32_t i = -32768; i < 32768; i++){
        t = Barrett_round_reduce(i);
        if(t > ubound){
            ubound = t;
        }
        if((t < 0) && (t < -ubound)){
            ubound = -t;
        }
    }
    printf("ubound of round reduction: %4d\n", ubound);

// ================
// Modular multiplication

    for(size_t j = 1; j <= 11; j++){
        for(int16_t i = -1664; i <= 1664; i++){
            t = mulmod(i, j);
            assert((-KYBER_Q / 2 <= t) && (t <= KYBER_Q / 2));
        }
    }

// ================
// Quotient

    for(int16_t i = -1664; i <= 1664; i++){
        assert(quotient_D_sign(i, 1) == Barrett_quotient_1(i));
    }

    for(int16_t i = -1664; i <= 1664; i++){
        assert(quotient_D_sign(i, 2) == Barrett_quotient_2(i));
    }

    for(int16_t i = -1664; i <= 1664; i++){
        assert(quotient_D_sign(i, 3) == Barrett_quotient_3(i));
    }

    for(int16_t i = -1664; i <= 1664; i++){
        assert(quotient_D_sign(i, 4) == Barrett_quotient_4(i));
    }

    for(int16_t i = -1664; i <= 1664; i++){
        assert(quotient_D_sign(i, 5) == Barrett_quotient_5(i));
    }

    for(int16_t i = -1664; i <= 1664; i++){
        assert(quotient_D_sign(i, 6) == Barrett_quotient_6(i));
    }

    for(int16_t i = -1664; i <= 1664; i++){
        assert(quotient_D_sign(i, 7) == Barrett_quotient_7(i));
    }

    for(int16_t i = -1664; i <= 1664; i++){
        assert(quotient_D_sign(i, 8) == Barrett_quotient_8(i));
    }

    for(int16_t i = -1664; i <= 1664; i++){
        assert(quotient_D_sign(i, 9) == Barrett_quotient_9(i));
    }

    for(int16_t i = -1664; i <= 1664; i++){
        assert(quotient_D_sign(i, 10) == Barrett_quotient_10(i));
    }

    for(int16_t i = -1664; i <= 1664; i++){
        assert(quotient_D_sign(i, 11) == Barrett_quotient_11(i));
    }

// ================
// Compression

    for(int16_t i = -1664; i <= 1664; i++){
        assert(compress_D(i, 1) == Barrett_compress_1(i));
    }

    for(int16_t i = -1664; i <= 1664; i++){
        assert(compress_D(i, 2) == Barrett_compress_2(i));
    }

    for(int16_t i = -1664; i <= 1664; i++){
        assert(compress_D(i, 3) == Barrett_compress_3(i));
    }

    for(int16_t i = -1664; i <= 1664; i++){
        assert(compress_D(i, 4) == Barrett_compress_4(i));
    }

    for(int16_t i = -1664; i <= 1664; i++){
        assert(compress_D(i, 5) == Barrett_compress_5(i));
    }

    for(int16_t i = -1664; i <= 1664; i++){
        assert(compress_D(i, 6) == Barrett_compress_6(i));
    }

    for(int16_t i = -1664; i <= 1664; i++){
        assert(compress_D(i, 7) == Barrett_compress_7(i));
    }

    for(int16_t i = -1664; i <= 1664; i++){
        assert(compress_D(i, 8) == Barrett_compress_8(i));
    }

    for(int16_t i = -1664; i <= 1664; i++){
        assert(compress_D(i, 9) == Barrett_compress_9(i));
    }

    for(int16_t i = -1664; i <= 1664; i++){
        assert(compress_D(i, 10) == Barrett_compress_10(i));
    }

    for(int16_t i = -1664; i <= 1664; i++){
        assert(compress_D(i, 11) == Barrett_compress_11(i));
    }



    printf("Test finished!\n");

}

