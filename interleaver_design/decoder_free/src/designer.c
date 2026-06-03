/**
 * @file designer.c
 * @brief Decoder-free DPCCC interleaver design algorithm.
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <pthread.h>
#include <unistd.h>
#include "../../../shared/config.h"
#include "../../../shared/encoder.h"
#include "../../../shared/interleaver.h"
#include "../../../shared/utils.h"
#include "../include/design_config.h"
#include "../include/rtz_search.h"
#include "../include/designer.h"

int StopRecursion = 0;

void Build_Multiple_Interleaver_Iteratively(
    int Places[][MAX_SEGMENTS],
    int number_of_repetition[MAX_SEGMENTS],
    int circular_states[MAX_SEGMENTS][32],
    double R, int threshold)
{
    int i, j, l;
    int small_fs = SIZE / CUTS;

    /* Allocate per-cut interleaver tables */
    int **OriginalAddress   = malloc(CUTS * sizeof(int *));
    int **InterleavedAddress = malloc(CUTS * sizeof(int *));
    for (i = 0; i < CUTS; i++) {
        OriginalAddress[i]    = malloc(small_fs * sizeof(int));
        InterleavedAddress[i] = malloc(small_fs * sizeof(int));
    }

    int **Total_S  = malloc(CUTS * sizeof(int *));
    for (i = 0; i < CUTS; i++) Total_S[i] = malloc(ARP_Q * sizeof(int));
    int *Total_P   = malloc(CUTS * sizeof(int));
    int *S         = malloc(ARP_Q * sizeof(int));

    int results[small_fs], count = 0;
    Choose_Period(results, &count, small_fs);

    FILE *temps;
    long long iterations_done;

rep_begin:
    iterations_done = 0;

    for (l = 0; l < CUTS; l++) {
    repeate:
        iterations_done++;
        if (iterations_done > 30000) goto rep_begin;

        int P;
        for (i = 0; i < 100; i++) {
            if (SHUFFLE) shuffleArray(results, count);
            P = results[0];
            StopRecursion = 0;
            Generate_ARP_Interleaver_Iteratively(l, 0, OriginalAddress,
                                                 InterleavedAddress,
                                                 threshold, S, P);
            if (i == 5) goto rep_begin;
            if (StopRecursion == 1) {
                ARPInterleaverInitialization(small_fs, P, ARP_Q, S,
                                             InterleavedAddress[l],
                                             OriginalAddress[l]);
                StopRecursion = 0;
                break;
            }
        }

        if (VERBOSE) {
            printf("Cut %d:  S = [", l);
            for (i = 0; i < ARP_Q; i++) printf("%d%s", S[i], i < ARP_Q-1 ? ", " : "");
            printf("]  P = %d\n", P);
        }

        temps = fopen("temps.txt", "a");
        fprintf(temps, "Cut %d:  S = [", l);
        for (i = 0; i < ARP_Q; i++) fprintf(temps, "%d%s", S[i], i < ARP_Q-1 ? ", " : "");
        fprintf(temps, "]  P = %d\n", P);
        fclose(temps);

        copy_array(S, Total_S[l], ARP_Q);
        Total_P[l] = P;
    }

    /* Write result in C-initialiser format */
    FILE *file = fopen("result.txt", "a");
    fprintf(file, "int total_S[%d][%d] = { ", CUTS, ARP_Q);
    for (i = 0; i < CUTS; i++) {
        fprintf(file, "{");
        for (j = 0; j < ARP_Q; j++)
            fprintf(file, (j == ARP_Q-1) ? "%d}" : "%d, ", Total_S[i][j]);
        fprintf(file, (i == CUTS-1) ? " };" : " ,\n ");
    }
    fprintf(file, "\nint total_P[%d] = {", CUTS);
    for (i = 0; i < CUTS; i++)
        fprintf(file, (i == CUTS-1) ? "%d};\n" : "%d, ", Total_P[i]);
    fprintf(file, "distance = %d\n\n", threshold);
    fclose(file);

    for (i = 0; i < CUTS; i++) { free(OriginalAddress[i]); free(InterleavedAddress[i]); free(Total_S[i]); }
    free(OriginalAddress); free(InterleavedAddress);
    free(Total_S); free(Total_P); free(S);
}


void Generate_ARP_Interleaver_Iteratively(int the_layer, int current_layer,
                                           int **OriginalAddress,
                                           int **InterleavedAddress,
                                           int threshold, int *S, int P)
{
    int small_fs = SIZE / CUTS;

    /* All Q layers filled: mark success */
    if (current_layer == ARP_Q) {
        StopRecursion = 1;
        return;
    }

    int array_s[small_fs];
    int size = 0;
    Add_ARP_LAYER(the_layer, current_layer, OriginalAddress, InterleavedAddress,
                  threshold, S, P, array_s, &size, SHIFTS_PER_LAYER);

    for (int i = 0; i < size && StopRecursion == 0; i++) {
        S[current_layer] = array_s[i];
        ARPInterleaverInitialization_Per_Layer(small_fs, current_layer, P,
                                               ARP_Q, S,
                                               InterleavedAddress[the_layer],
                                               OriginalAddress[the_layer]);
        Generate_ARP_Interleaver_Iteratively(the_layer, current_layer + 1,
                                              OriginalAddress, InterleavedAddress,
                                              threshold, S, P);
        if (StopRecursion == 0)
            Arp_DE_Interleave_Layer(small_fs, current_layer, P, ARP_Q, S,
                                    InterleavedAddress[the_layer],
                                    OriginalAddress[the_layer]);
    }
}


void Add_ARP_LAYER(int the_layer, int current_shift_layer,
                   int **OriginalAddress, int **InterleavedAddress,
                   int threshold, int *S, int P,
                   int *array_s, int *size, int max_size_layer)
{
    int small_fs = SIZE / CUTS;
    int lcm      = (ARP_Q * PUNCT_MASK) / gcd(ARP_Q, PUNCT_MASK);
    int all_sl[small_fs];
    int reach = small_fs;
    int mult  = 0;

    for (int i = 0; i < reach; i++) all_sl[i] = i;
    if (SHUFFLE) shuffleArray(all_sl, reach);

    int g2[] = IW2_GIRTH_BOUNDS;
    int p2[] = IW2_PERIOD_BUDGETS;
    int gM[] = IWM_GIRTH_BOUNDS;
    int pM[] = IWM_STEP_BUDGETS;

    for (int pp = 0; pp < reach && *size < max_size_layer; pp++) {
        int sl = all_sl[pp];
        /* ARP regularity check */
        S[current_shift_layer] = sl;
        int ok = 1;
        for (int j = 0; j < current_shift_layer; j++)
            if ((P * current_shift_layer + sl) % ARP_Q ==
                (P * j + S[j]) % ARP_Q) { ok = 0; break; }
        if (!ok) continue;

        ARPInterleaverInitialization_Per_Layer(small_fs, current_shift_layer,
                                               P, ARP_Q, S,
                                               InterleavedAddress[the_layer],
                                               OriginalAddress[the_layer]);

        /* RTZ search over all active positions */
        int rejected = 0;
        for (int i = 0; i < SIZE && !rejected; i++) {
            if (i % CUTS > the_layer) continue;  /* skip inactive layers */
            int layer = i % CUTS;
            int pos   = i / CUTS;

            /* Weight-2 periodic kernel */
            for (int jj = 0; jj < IW2_MAX_WEIGHT && !rejected; jj++) {
                if (Min_Girth_One_Layer_JUMP_1_7_Iteratively_ARP_INTERLEAVER(
                        current_shift_layer, pos, layer,
                        OriginalAddress, InterleavedAddress,
                        threshold, &mult,
                        p2[jj], p2[jj], g2[jj], 2, ENCODER_PERIOD) == 0)
                    rejected = 1;
            }
            /* General weight-3+ kernel */
            for (int jj = 3; jj <= IWM_MAX_WEIGHT && !rejected; jj++) {
                if (Min_Girth_One_Layer_JUMP_1_7_Iteratively_ARP_INTERLEAVER(
                        current_shift_layer, pos, layer,
                        OriginalAddress, InterleavedAddress,
                        threshold, &mult,
                        pM[jj], pM[jj], gM[jj], 0, 1) == 0)
                    rejected = 1;
            }
        }

        Arp_DE_Interleave_Layer(small_fs, current_shift_layer, P, ARP_Q, S,
                                 InterleavedAddress[the_layer],
                                 OriginalAddress[the_layer]);

        if (!rejected) {
            array_s[*size] = sl;
            (*size)++;
        }
    }
}
