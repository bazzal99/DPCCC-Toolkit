/**
 * @file encoder.c
 * @brief DPCCC RSC encoder and circular state initialisation.
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "config.h"
#include "encoder.h"
#include "interleaver.h"
#include "utils.h"

/* =========================================================================
 * Global encoder tables
 * ========================================================================= */

int G1[] = G1_INIT;
int G2[] = G2_INIT;
int CIRCULAR_STATES[10][32];

static int DATA_PUNCT[]    = DATA_PUNCT_INIT;
static int PARITY_K[]      = PARITY_K_INIT;
static int PARITY_K_Q[]    = PARITY_K_Q_INIT;


/* =========================================================================
 * Circular state initialisation
 * ========================================================================= */

void Initialize_States(int circularstates[], int framesize)
{
    int i, j, BoolCircular;
    int *d         = malloc(framesize * sizeof(int));
    int *myparity  = malloc(framesize * sizeof(int));
    int prevstate[DELAYS], nextstate[DELAYS];
    int nextstateint;

    circularstates[0] = 0;
    for (i = 1; i < (int)pow(2, DELAYS); i++)
        circularstates[i] = 9999;

    do {
        BoolCircular = 1;
        generateRandomBinaryArray(d, framesize);

        int_to_binary(0, DELAYS, prevstate);
        Conv_Encode_number_of_delays(prevstate, d, framesize,
                                     G1, G2, myparity, nextstate, DELAYS);
        nextstateint = binary_to_int(nextstate, DELAYS);

        Conv_Encode_number_of_delays(nextstate, d, framesize,
                                     G1, G2, myparity, nextstate, DELAYS);
        circularstates[binary_to_int(nextstate, DELAYS)] = nextstateint;

        for (j = 1; j < (int)pow(2, DELAYS); j++)
            if (circularstates[j] == 9999) { BoolCircular = 0; break; }

    } while (BoolCircular == 0);

    free(d);
    free(myparity);
}


/* =========================================================================
 * Encoding
 * ========================================================================= */

int Encode_1_bit_number_of_delays(int prevstate[], int input,
                                   int g1[], int g2[],
                                   int nextstate[], int num_delays)
{
    int i, parity, x = 0;
    for (i = num_delays; i > 0; i--) x ^= g2[i] & prevstate[i - 1];
    x ^= input;
    parity      = g1[0] ? x : 0;
    nextstate[0] = x;
    for (i = 1; i < num_delays; i++)      nextstate[i] = prevstate[i - 1];
    for (i = 1; i < num_delays + 1; i++)  if (g1[i]) parity ^= prevstate[i - 1];
    return parity;
}

void Conv_Encode_number_of_delays(int prevstate[], int input[], int size,
                                   int g1[], int g2[], int parity[],
                                   int nextstate[], int delays)
{
    int i;
    int theprevstate[delays];
    copy_array(prevstate, theprevstate, delays);
    for (i = 0; i < size; i++) {
        parity[i] = Encode_1_bit_number_of_delays(theprevstate, input[i],
                                                   g1, g2, nextstate, delays);
        copy_array(nextstate, theprevstate, delays);
    }
}

/**
 * Encode the full DPCCC frame.
 * Each segment is tail-biting encoded independently using its own
 * circular starting state.
 */
void Encode(int input[], int parity[],
            int circularstates[MAX_SEGMENTS][32],
            int number_of_repetition[])
{
    int i, j, xx, startstate, nextstateint;
    int prevstate[DELAYS], nextstate[DELAYS];

    if (!CIRCULAR_ENCODING) {
        int zeros[DELAYS];
        for (i = 0; i < DELAYS; i++) zeros[i] = 0;
        Conv_Encode_number_of_delays(zeros, input, number_of_repetition[0],
                                     G1, G2, parity, nextstate, DELAYS);
        return;
    }

    for (i = 0; i < MAX_SEGMENTS; i++) {
        if (number_of_repetition[i] == 0) break;
        xx = 0;
        for (j = 0; j < i; j++) xx += number_of_repetition[j];

        /* Pass 1: determine ending state */
        int_to_binary(0, DELAYS, prevstate);
        Conv_Encode_number_of_delays(prevstate, input + xx,
                                     number_of_repetition[i],
                                     G1, G2, parity + xx, nextstate, DELAYS);
        nextstateint = binary_to_int(nextstate, DELAYS);

        /* Pass 2: encode from circular starting state */
        startstate = circularstates[i][nextstateint];
        int_to_binary(startstate, DELAYS, prevstate);
        Conv_Encode_number_of_delays(prevstate, input + xx,
                                     number_of_repetition[i],
                                     G1, G2, parity + xx, nextstate, DELAYS);
    }
}

void Encodes(int input[], int parity[],
             int circularstates[], int framesize)
{
    int prevstate[DELAYS], nextstate[DELAYS];
    int startstate, nextstateint;

    int_to_binary(0, DELAYS, prevstate);
    Conv_Encode_number_of_delays(prevstate, input, framesize,
                                 G1, G2, parity, nextstate, DELAYS);
    nextstateint = binary_to_int(nextstate, DELAYS);

    startstate = circularstates[nextstateint];
    int_to_binary(startstate, DELAYS, prevstate);
    Conv_Encode_number_of_delays(prevstate, input, framesize,
                                 G1, G2, parity, nextstate, DELAYS);
}

void DrawTrellis(int g1[], int g2[], int delays, int trellis[][4])
{
    int i, x[8], nextstate[8];
    for (i = 0; i < (int)pow(2, delays); i++) {
        int_to_binary(i, delays, x);
        trellis[i][0] = Encode_1_bit_number_of_delays(x, 0, g1, g2, nextstate, delays);
        trellis[i][1] = binary_to_int(nextstate, delays);
        trellis[i][2] = Encode_1_bit_number_of_delays(x, 1, g1, g2, nextstate, delays);
        trellis[i][3] = binary_to_int(nextstate, delays);
    }
}


/* =========================================================================
 * Puncturing
 * ========================================================================= */

void puncmask_all_frame_FrameSize(int x[], int framesize, int type)
{
    for (int i = 0; i < framesize; i++) {
        if      (type == 0) x[i] = DATA_PUNCT[i % PUNCT_MASK];
        else if (type == 1) x[i] = PARITY_K[i % PUNCT_MASK];
        else                x[i] = PARITY_K_Q[i % PUNCT_MASK];
    }
}


/* =========================================================================
 * Distance evaluation
 * ========================================================================= */

int parity_distance_circular_With_Puncturing(int x[], int size,
                                              int type, int frame_size)
{
    int *frame = calloc(frame_size, sizeof(int));
    int *parity = malloc(frame_size * sizeof(int));
    int *mask   = malloc(frame_size * sizeof(int));
    int *cs     = malloc(32 * sizeof(int));
    int d = 0;

    for (int i = 0; i < size; i++) frame[x[i]] = 1;
    Initialize_States(cs, frame_size);
    Encodes(frame, parity, cs, frame_size);
    puncmask_all_frame_FrameSize(mask, frame_size, type);
    Array_AND(parity, mask, parity, frame_size);
    d = Number_Of_Ones(parity, frame_size);

    free(frame); free(parity); free(mask); free(cs);
    return d;
}

int Turbo_distance_circular_With_Puncturing(int x[], int size,
                                             int **Original_address)
{
    int *interleaver  = calloc(FRAME_SIZE, sizeof(int));
    int *impulse_resp = malloc(FRAME_SIZE * sizeof(int));
    int *mask         = malloc(FRAME_SIZE * sizeof(int));
    int *cs           = malloc(32 * sizeof(int));
    int frame_size = SIZE / CUTS;
    int d = size;

    for (int i = 0; i < size; i++) {
        int layer = x[i] % CUTS;
        int pos   = x[i] / CUTS;
        int mapped = Original_address[layer][pos];
        if (mapped == -1 || mapped> SIZE || mapped<0) { free(interleaver); free(impulse_resp); free(mask); free(cs); return 9999; }
        interleaver[mapped] = 1;
    }
    Initialize_States(cs, frame_size);
    for (int l = 0; l < CUTS; l++) {
        Encodes(interleaver + l * frame_size, impulse_resp + l * frame_size, cs, frame_size);
        puncmask_all_frame_FrameSize(mask + l * frame_size, frame_size, 2);
    }
    Array_AND(impulse_resp, mask, impulse_resp, CUTS * frame_size);
    d += Number_Of_Ones(impulse_resp, CUTS * frame_size);

    free(interleaver); free(impulse_resp); free(mask); free(cs);
    return d;
}

int Parity_calc(int output[], int Places[][MAX_SEGMENTS],
                int OriginalAddress[], int number_of_repetition[MAX_SEGMENTS],
                int circular_states[MAX_SEGMENTS][32],
                int DataMask[], int Parity_Mask[])
{
    int input[SIZE], d[FRAME_SIZE], d_interleaved[FRAME_SIZE], cp[FRAME_SIZE];
    int p = 0;

    copy_array(output, input, SIZE);
    Non_uniform_repeate(input, d, Places);
    Array_AND(input, DataMask, input, SIZE);
    p += Number_Of_Ones(input, SIZE);

    ARP_Interleave(d, d_interleaved, OriginalAddress);
    Encode(d_interleaved, cp, circular_states, number_of_repetition);
    Array_AND(cp, Parity_Mask, cp, FRAME_SIZE);
    p += Number_Of_Ones(cp, FRAME_SIZE);

    return p;
}
