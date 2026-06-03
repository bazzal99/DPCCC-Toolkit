/**
 * @file main.c
 * @brief DPCCC Monte Carlo BER/FER simulator.
 *
 * Simulates the Bit Error Rate (BER) and Frame Error Rate (FER) of a
 * DPCCC code over AWGN using a turbo (log-MAP) decoder.
 *
 * The interleaver parameters are loaded via ARP_SEMI_Initialization()
 * in shared/interleaver.c -- update that function with your designed S, P.
 *
 * Configuration:
 *   ../../shared/config.h   -- encoder, frame size, puncturing, SNR range
 *   include/sim_config.h    -- termination mode, RNG seed, target errors
 *
 * Build: make    Run: ./dpccc_sim
 *
 * Output: result.txt -- BER and FER arrays per SNR point
 *
 * Note: This simulator is single-threaded and will be slow for large K
 * or high SNR. Consider parallelising Build_SNR() with pthreads for
 * production use.
 *
 * Author:  Mohammad Bazzal, IMT Atlantique, Lab-STICC, Brest, France
 * License: MIT
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "../../shared/config.h"
#include "../../shared/encoder.h"
#include "../../shared/interleaver.h"
#include "../../shared/decoder.h"
#include "../../shared/utils.h"
#include "../include/sim_config.h"

/* =========================================================================
 * Forward declarations
 * ========================================================================= */
void Build_SNR(int Places[][MAX_SEGMENTS],
               int OriginalAddress[], int InterleavedAddress[],
               int number_of_repetition[MAX_SEGMENTS],
               int circular_states[MAX_SEGMENTS][32],
               double R);


/* =========================================================================
 * main
 * ========================================================================= */
int main(void)
{
    unsigned int seed = (RNG_SEED == 0) ? (unsigned)time(NULL) : RNG_SEED;
    srand(seed);

    int circular_states[MAX_SEGMENTS][32];
    int Places[SIZE][MAX_SEGMENTS];
    int number_of_repetition[MAX_SEGMENTS];
    int OriginalAddress[FRAME_SIZE], InterleavedAddress[FRAME_SIZE];
    int i, j;

    for (i = 0; i < SIZE; i++)
        for (j = 0; j < MAX_SEGMENTS; j++)
            Places[i][j] = -1;
    for (i = 0; i < MAX_SEGMENTS; i++) number_of_repetition[i] = 0;

#if SEMI_CIRCULAR_TERMINATION == 2
    /* Special coupling: data block + CUTS parity sub-blocks */
    number_of_repetition[0] = SIZE;
    for (i = 0; i < CUTS; i++) number_of_repetition[i+1] = SIZE / CUTS;

    for (i = 0; i < SIZE; i++) Places[i][0] = i;
    int beginning = SIZE;
    for (i = 0; i < CUTS; i++)
        for (j = 0; j < SIZE / CUTS; j++)
            Places[i + j*CUTS][1] = beginning++;
#else
    /* Legacy: single full-frame encoder */
    number_of_repetition[0] = FRAME_SIZE;
    for (i = 0; i < SIZE; i++)
        for (j = 0; j < 2; j++)
            Places[i][j] = i + j * SIZE;
#endif

    if (CIRCULAR_ENCODING) {
        for (i = 0; i < MAX_SEGMENTS; i++) {
            if (number_of_repetition[i] == 0) break;
            Initialize_States(circular_states[i], number_of_repetition[i]);
        }
    }

    /* Load interleaver (hardcoded paper parameters — edit in interleaver.c) */
    ARP_SEMI_Initialization(OriginalAddress, InterleavedAddress,
                            number_of_repetition);

    /* Compute code rate */
    int tmp_d[]  = DATA_PUNCT_INIT;
    int tmp_k[]  = PARITY_K_INIT;
    int tmp_kq[] = PARITY_K_Q_INIT;
    int d_ones = 0, k_ones = 0, kq_ones = 0;
    for (i = 0; i < PUNCT_MASK; i++) {
        d_ones += tmp_d[i]; k_ones += tmp_k[i]; kq_ones += tmp_kq[i];
    }
    double R = (double)SIZE / (d_ones*(SIZE/PUNCT_MASK) +
                               k_ones*(SIZE/PUNCT_MASK) +
                               kq_ones*(SIZE/(PUNCT_MASK*CUTS))*CUTS);
    printf("Rate = %.4f\n", R);

    Build_SNR(Places, OriginalAddress, InterleavedAddress,
              number_of_repetition, circular_states, R);
    return 0;
}


/* =========================================================================
 * Build_SNR — main simulation loop
 * ========================================================================= */
void Build_SNR(int Places[][MAX_SEGMENTS],
               int OriginalAddress[], int InterleavedAddress[],
               int number_of_repetition[MAX_SEGMENTS],
               int circular_states[MAX_SEGMENTS][32],
               double R)
{
    int num_states = (int)pow(2, DELAYS);
    int trellis[8][4];
    DrawTrellis(G1, G2, DELAYS, trellis);

    int input[SIZE], d[FRAME_SIZE], d_int[FRAME_SIZE], cp[FRAME_SIZE], dec[SIZE];
    double R_X[SIZE], R_Y[FRAME_SIZE];
    double L_X[SIZE], L_X_perm[FRAME_SIZE], L_Y[FRAME_SIZE];
    double La_D[FRAME_SIZE], La_D_perm[FRAME_SIZE], L_D[FRAME_SIZE], L_tmp[FRAME_SIZE];
    double log_alpha[8][FRAME_SIZE+1], log_beta[8][FRAME_SIZE+1], log_gamma[8][2][FRAME_SIZE];
    double Alpha_init[MAX_SEGMENTS][8], Beta_init[MAX_SEGMENTS][8];

    int DataMask[SIZE], ParityMask[FRAME_SIZE];
    puncmask_all_frame_FrameSize(DataMask,   SIZE,        0);
    puncmask_all_frame_FrameSize(ParityMask,  SIZE,        1);
    for (int i = 0; i < CUTS; i++)
        puncmask_all_frame_FrameSize(ParityMask + SIZE + i*(SIZE/CUTS), SIZE/CUTS, 2);

    int    target_errors = TARGET_FRAME_ERRORS;
    double BER[40] = {0}, FER[40] = {0};
    int    snr_idx = 0;

    FILE *file;
    char filename[] = "result.txt";
    file = fopen(filename, "w"); fclose(file);

    for (double snr = SNR_MIN; snr <= SNR_MAX; snr += SNR_STEP) {
        double EbNo    = pow(10.0, snr / 10.0);
        double sigma_sq = 1.0 / (2.0 * R * EbNo);
        long long nerr = 0, total_frames = 0;
        int frame_err = 0;

        while (frame_err < target_errors) {
            generateRandomBinaryArray(input, SIZE);

            Non_uniform_repeate(input, d, Places);
            ARP_Interleave(d, d_int, OriginalAddress);
            Encode(d_int, cp, circular_states, number_of_repetition);

            for (int i = 0; i < SIZE;       i++) R_X[i] = 2.0*input[i] - 1.0;
            for (int i = 0; i < FRAME_SIZE; i++) R_Y[i] = 2.0*cp[i]    - 1.0;

            addAWGN(R_X, R_X, SIZE,       sigma_sq);
            addAWGN(R_Y, R_Y, FRAME_SIZE, sigma_sq);

            for (int i = 0; i < SIZE; i++)
                L_X[i] = DataMask[i] ? (2.0/sigma_sq)*R_X[i] : 0.0;
            for (int i = 0; i < FRAME_SIZE; i++) {
                L_Y[i]  = ParityMask[i] ? (2.0/sigma_sq)*R_Y[i] : 0.0;
                La_D[i] = 0.0;
            }

            Non_uniform_repeate_double(L_X, L_X_perm, Places);
            ARP_Interleave_double(L_X_perm, L_X_perm, OriginalAddress);
            ARP_Interleave_double(La_D,     La_D_perm, OriginalAddress);

            for (int j = 0; j < MAX_SEGMENTS && number_of_repetition[j]; j++)
                for (int i = 0; i < num_states; i++)
                    Alpha_init[j][i] = Beta_init[j][i] = log(1.0/num_states);

            for (int iter = 0; iter < MAX_ITERATIONS; iter++) {
                int offset = 0;
                for (int j = 0; j < MAX_SEGMENTS && number_of_repetition[j]; j++) {
                    decode(La_D_perm+offset, L_X_perm+offset, L_Y+offset,
                           L_D+offset, num_states, trellis,
                           log_alpha, log_beta, log_gamma,
                           Alpha_init[j], Beta_init[j],
                           number_of_repetition[j]);
                    offset += number_of_repetition[j];
                }
                ARP_DE_Interleave_double(L_D, L_D, InterleavedAddress);
                if (iter == MAX_ITERATIONS-1) break;
                for (int i = 0; i < SIZE; i++)
                    for (int j = 0; j < MAX_SEGMENTS && Places[i][j] != -1; j++)
                        L_tmp[Places[i][j]] = Calculate_Extrinsic(Places,i,j,L_D,L_X,La_D);
                for (int i = 0; i < FRAME_SIZE; i++) La_D[i] = L_tmp[i];
                ARP_Interleave_double(La_D, La_D_perm, OriginalAddress);
            }

            int errors = 0;
            for (int i = 0; i < SIZE; i++) {
                dec[i] = (L_D[i] > 0) ? 1 : 0;
                errors += (dec[i] != input[i]);
            }
            nerr += errors;
            total_frames++;
            if (errors) frame_err++;
        }

        BER[snr_idx] = (double)nerr  / ((double)SIZE * total_frames);
        FER[snr_idx] = (double)frame_err / total_frames;

        printf("SNR = %.2f dB  BER = %.4e  FER = %.4e  frames = %lld\n",
               snr, BER[snr_idx], FER[snr_idx], total_frames);

        /* Adaptive error target */
        if      (BER[snr_idx] > 1e-2) target_errors = TARGET_FRAME_ERRORS;
        else if (BER[snr_idx] > 1e-4) target_errors = TARGET_FRAME_ERRORS;
        else if (BER[snr_idx] > 1e-5) target_errors = 200;
        else if (BER[snr_idx] > 1e-6) target_errors = 150;
        else                           target_errors = 100;

        /* Write to file */
        file = fopen(filename, "a");
        fprintf(file, "SNR = %.2f\n", snr);
        fprintf(file, "BER = [");
        for (int i = 0; i <= snr_idx; i++)
            fprintf(file, (i == snr_idx) ? "%.7e]\n" : "%.7e, ", BER[i]);
        fprintf(file, "FER = [");
        for (int i = 0; i <= snr_idx; i++)
            fprintf(file, (i == snr_idx) ? "%.7e]\n\n" : "%.7e, ", FER[i]);
        fclose(file);

        snr_idx++;
    }
}
