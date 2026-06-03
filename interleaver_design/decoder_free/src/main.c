/**
 * @file main.c
 * @brief DPCCC interleaver design - decoder-free RTZ-based method.
 *
 * Designs CUTS independent ARP sub-interleavers for a DPCCC code using
 * a decoder-free mHD estimation technique. The algorithm builds each
 * sub-interleaver layer by layer, rejecting any shift that produces a
 * codeword below the distance threshold.
 *
 * Configuration:
 *   ../../shared/config.h   -- encoder, frame structure, puncturing
 *   include/design_config.h -- threshold, search depth
 *
 * Build: make    Run: ./dpccc_decoder_free
 *
 * Output:
 *   result.txt  -- valid interleavers in C-initialiser format
 *   temps.txt   -- progress log
 *
 * Reference:
 *   M. Bazzal et al., "Distance-centric joint interleaver and structural
 *   code design for concatenated convolutional codes,"
 *   IEEE Open Journal of the Communications Society, 2025.
 *
 * Author:  Mohammad Bazzal, IMT Atlantique, Lab-STICC, Brest, France
 * License: MIT
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "../../../shared/config.h"
#include "../../../shared/encoder.h"
#include "../../../shared/interleaver.h"
#include "../../../shared/utils.h"
#include "../include/design_config.h"
#include "../include/designer.h"

int main(void)
{
    srand(time(NULL));

    int small_fs = SIZE / CUTS;
    Initialize_States(CIRCULAR_STATES[0], small_fs);
    Initialize_States(CIRCULAR_STATES[1], SIZE);

    int circular_states[MAX_SEGMENTS][32];
    int Places[SIZE][MAX_SEGMENTS];
    int number_of_repetition[MAX_SEGMENTS];
    int i, j;

    for (i = 0; i < SIZE; i++)
        for (j = 0; j < MAX_SEGMENTS; j++)
            Places[i][j] = -1;
    for (i = 0; i < MAX_SEGMENTS; i++) number_of_repetition[i] = 0;

    /* Special-coupling DPCCC structure */
    number_of_repetition[0] = SIZE;
    for (i = 0; i < CUTS; i++) number_of_repetition[i + 1] = SIZE / CUTS;

    for (i = 0; i < SIZE; i++) Places[i][0] = i;
    int beginning = SIZE;
    for (i = 0; i < CUTS; i++)
        for (j = 0; j < SIZE / CUTS; j++)
            Places[i + j * CUTS][1] = beginning++;

    for (i = 0; i < MAX_SEGMENTS; i++) {
        if (number_of_repetition[i] == 0) break;
        Initialize_States(circular_states[i], number_of_repetition[i]);
    }

    int data_ones    = 0;
    int pk_ones      = 0;
    int pkq_ones     = 0;
    int tmp_d[]  = DATA_PUNCT_INIT;
    int tmp_k[]  = PARITY_K_INIT;
    int tmp_kq[] = PARITY_K_Q_INIT;
    for (i = 0; i < PUNCT_MASK; i++) { data_ones += tmp_d[i]; pk_ones += tmp_k[i]; pkq_ones += tmp_kq[i]; }
    int sbits_data    = data_ones * (SIZE / PUNCT_MASK);
    int sbits_parity1 = pk_ones  * (SIZE / PUNCT_MASK);
    int sbits_parity2 = pkq_ones * (SIZE / (PUNCT_MASK * CUTS));
    double R = (double)SIZE / (sbits_data + sbits_parity1 + sbits_parity2 * CUTS);
    printf("Code rate R = %.4f\n", R);

    FILE *f = fopen("result.txt", "w"); fclose(f);

    int threshold = DESIGN_THRESHOLD;
    for (int t = 0; t < 100; t++) {
        Build_Multiple_Interleaver_Iteratively(Places, number_of_repetition,
                                               circular_states, R, threshold);
        if (!WITH_MULTIPLICITY) threshold++;
    }
    return 0;
}
