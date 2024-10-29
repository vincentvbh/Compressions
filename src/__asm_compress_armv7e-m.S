
.syntax unified
.thumb

.align 4
.global __asm_compress_1
.type __asm_compress_1, %function
__asm_compress_1:
    push.w {lr}

    movw r3, #:lower16:2580335
    movt r3, #:upper16:2580335

    ldrsh.w r2, [r1]

    smmulr r2, r2, r3
    ubfx.w r2, r2, #0, #1

    strh.w r2, [r0]

    pop.w {pc}

.align 4
.global __asm_compress_4
.type __asm_compress_4, %function
__asm_compress_4:
    push.w {lr}

    movw r3, #:lower16:20642679
    movt r3, #:upper16:20642679

    ldrsh.w r2, [r1]

    smmulr r2, r2, r3
    ubfx.w r2, r2, #0, #4

    strh.w r2, [r0]

    pop.w {pc}

.align 4
.global __asm_compress_10
.type __asm_compress_10, %function
__asm_compress_10:
    push.w {lr}

    movw r3, #:lower16:1321131424
    movt r3, #:upper16:1321131424

    ldrsh.w r2, [r1]

    smmulr r2, r2, r3
    ubfx.w r2, r2, #0, #10

    strh.w r2, [r0]

    pop.w {pc}

