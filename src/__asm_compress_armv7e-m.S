
.syntax unified
.thumb

.align 4
.global __asm_compress1
.type __asm_compress1, %function
__asm_compress1:
    push.w {lr}

    movw r3, #:lower16:2580335
    movt r3, #:upper16:2580335

    ldrsh.w r2, [r1]

    smmulr r2, r2, r3
    ubfx.w r2, r2, #0, #1

    strh.w r2, [r0]

    pop.w {pc}

.align 4
.global __asm_compress4
.type __asm_compress4, %function
__asm_compress4:
    push.w {lr}

    movw r3, #:lower16:20642679
    movt r3, #:upper16:20642679

    ldrsh.w r2, [r1]

    smmulr r2, r2, r3
    ubfx.w r2, r2, #0, #4

    strh.w r2, [r0]

    pop.w {pc}

.align 4
.global __asm_compress10
.type __asm_compress10, %function
__asm_compress10:
    push.w {lr}

    movw r3, #:lower16:1321131424
    movt r3, #:upper16:1321131424

    ldrsh.w r2, [r1]

    smmulr r2, r2, r3
    ubfx.w r2, r2, #0, #10

    strh.w r2, [r0]

    pop.w {pc}

.align 4
.global __asm_poly_compress1
.type __asm_poly_compress1, %function
__asm_poly_compress1:
    push.w {r4-r12, lr}

    movw r2, #:lower16:2580335
    movt r2, #:upper16:2580335

.rept 32

    ldrsh.w r5, [r1, #2]
    ldrsh.w r6, [r1, #4]
    ldrsh.w r7, [r1, #6]
    ldrsh.w r8, [r1, #8]
    ldrsh.w r9, [r1, #10]
    ldrsh.w r10, [r1, #12]
    ldrsh.w r11, [r1, #14]
    ldrsh.w r4, [r1], #16

    smmulr r4, r4, r2
    smmulr r5, r5, r2
    smmulr r6, r6, r2
    smmulr r7, r7, r2
    smmulr r8, r8, r2
    smmulr r9, r9, r2
    smmulr r10, r10, r2
    smmulr r11, r11, r2

    ubfx.w r4, r4, #0, #1
    ubfx.w r5, r5, #0, #1
    ubfx.w r6, r6, #0, #1
    ubfx.w r7, r7, #0, #1
    ubfx.w r8, r8, #0, #1
    ubfx.w r9, r9, #0, #1
    ubfx.w r10, r10, #0, #1
    ubfx.w r11, r11, #0, #1

    orr.w r4, r4, r5, lsl #1
    orr.w r4, r4, r6, lsl #2
    orr.w r4, r4, r7, lsl #3
    orr.w r4, r4, r8, lsl #4
    orr.w r4, r4, r9, lsl #5
    orr.w r4, r4, r10, lsl #6
    orr.w r4, r4, r11, lsl #7

    strb.w r4, [r0], #1

.endr

    pop.w {r4-r12, pc}

.align 4
.global __asm_poly_compress4
.type __asm_poly_compress4, %function
__asm_poly_compress4:
    push.w {r4-r12, lr}

    movw r2, #:lower16:20642679
    movt r2, #:upper16:20642679

.rept 64

    ldrsh.w r5, [r1, #2]
    ldrsh.w r6, [r1, #4]
    ldrsh.w r7, [r1, #6]
    ldrsh.w r4, [r1], #8

    smmulr r4, r4, r2
    smmulr r5, r5, r2
    smmulr r6, r6, r2
    smmulr r7, r7, r2

    ubfx.w r4, r4, #0, #4
    ubfx.w r5, r5, #0, #4
    ubfx.w r6, r6, #0, #4
    ubfx.w r7, r7, #0, #4

    orr.w r4, r4, r5, lsl #4
    orr.w r4, r4, r6, lsl #8
    orr.w r4, r4, r7, lsl #12

    strh.w r4, [r0], #2

.endr

    pop.w {r4-r12, pc}

.align 4
.global __asm_poly_compress10
.type __asm_poly_compress10, %function
__asm_poly_compress10:
    push.w {r4-r12, lr}

    movw r2, #:lower16:1321131424
    movt r2, #:upper16:1321131424

.rept 64

    ldrsh.w r5, [r1, #2]
    ldrsh.w r6, [r1, #4]
    ldrsh.w r7, [r1, #6]
    ldrsh.w r4, [r1], #8

    smmulr r4, r4, r2
    smmulr r5, r5, r2
    smmulr r6, r6, r2
    smmulr r7, r7, r2

    ubfx.w r4, r4, #0, #10
    ubfx.w r5, r5, #0, #10
    ubfx.w r6, r6, #0, #10
    ubfx.w r7, r7, #0, #10

    orr.w r4, r4, r5, lsl #10
    orr.w r4, r4, r6, lsl #20
    orr.w r4, r4, r7, lsl #30
    lsr.w r7, r7, #2

    strb.w r7, [r0, #4]
    str.w r4, [r0], #5

.endr

    pop.w {r4-r12, pc}

