/**
 * @file designer.c
 * @brief N-impulse DPCCC interleaver design implementation.
 *
 * Uses a parallel threaded decoder to evaluate each ARP candidate.
 * Each thread tests a different pair (or triple) of pinned bit positions
 * and reports if any decoding produces a codeword below the threshold.
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <pthread.h>
#include <unistd.h>
#include "../../../shared/config.h"
#include "../../../shared/encoder.h"
#include "../../../shared/interleaver.h"
#include "../../../shared/decoder.h"
#include "../../../shared/utils.h"
#include "../include/design_config.h"
#include "../include/designer.h"

/* =========================================================================
 * Globals
 * ========================================================================= */

int leave          = 0;
int TOTAL_DISTANCE = 0;

static pthread_mutex_t lock;

/* =========================================================================
 * Thread data structure
 * ========================================================================= */

typedef struct {
    int  thread_id;
    int  (*Places)[MAX_SEGMENTS];
    int  *OriginalAddress;
    int  *InterleavedAddress;
    int  *number_of_repetition;
    int  (*circular_states)[32];
    double R;
    double snr;
    int  layer;
} ThreadData;


/* =========================================================================
 * Thread worker
 * ========================================================================= */

static void *impulse_thread(void *arg)
{
    ThreadData *data  = (ThreadData *)arg;
    int thread_id     = data->thread_id;
    double R          = data->R;
    double snr        = data->snr;
    int  layer        = data->layer;
    int  (*Places)[MAX_SEGMENTS]  = data->Places;
    int  *OriginalAddress         = data->OriginalAddress;
    int  *InterleavedAddress      = data->InterleavedAddress;
    int  *number_of_repetition    = data->number_of_repetition;
    int  (*circular_states)[32]   = data->circular_states;

    int num_states  = (int)pow(2, DELAYS);
    double EbNo     = pow(10.0, snr / 10.0);
    double sigma_sq = 1.0 / (2.0 * R * EbNo);

    int input[SIZE], d[FRAME_SIZE], d_int[FRAME_SIZE], cp[FRAME_SIZE], dec[FRAME_SIZE];
    double R_X[SIZE], R_Y[FRAME_SIZE];
    double L_X[SIZE], L_X_perm[FRAME_SIZE], L_Y[FRAME_SIZE];
    double La_D[FRAME_SIZE], La_D_perm[FRAME_SIZE], L_D[FRAME_SIZE], L_tmp[FRAME_SIZE];
    double log_alpha[8][FRAME_SIZE+1], log_beta[8][FRAME_SIZE+1], log_gamma[8][2][FRAME_SIZE];
    double Alpha_init[MAX_SEGMENTS][8], Beta_init[MAX_SEGMENTS][8];

    int DataMask[SIZE], ParityMask[FRAME_SIZE];
    puncmask_all_frame_FrameSize(DataMask,  SIZE,        0);
    puncmask_all_frame_FrameSize(ParityMask, SIZE,        1);
    for (int i = 0; i < CUTS; i++)
        puncmask_all_frame_FrameSize(ParityMask + SIZE + i*(SIZE/CUTS), SIZE/CUTS, 2);

    int trellis[8][4];
    DrawTrellis(G1, G2, DELAYS, trellis);

    int pos1 = thread_id;

    for (int pos2 = 0; pos2 < SIZE; pos2++) {
        if (leave) pthread_exit(NULL);

        /* All-zero input */
        for (int i = 0; i < SIZE; i++) input[i] = 0;
        Non_uniform_repeate(input, d, Places);
        ARP_Interleave(d, d_int, OriginalAddress);
        Encode(d_int, cp, circular_states, number_of_repetition);

        /* Modulate */
        for (int i = 0; i < SIZE;       i++) R_X[i] = 2.0*input[i] - 1.0;
        for (int i = 0; i < FRAME_SIZE; i++) R_Y[i] = 2.0*cp[i]    - 1.0;

        addAWGN(R_X, R_X, SIZE,       sigma_sq);
        addAWGN(R_Y, R_Y, FRAME_SIZE, sigma_sq);

        /* LLRs */
        for (int i = 0; i < SIZE; i++)
            L_X[i] = DataMask[i] ? (2.0/sigma_sq)*R_X[i] : 0.0;
        for (int i = 0; i < FRAME_SIZE; i++) {
            L_Y[i]  = ParityMask[i] ? (2.0/sigma_sq)*R_Y[i] : 0.0;
            La_D[i] = 0.0;
        }

        /* Pin impulse positions */
        L_X[pos1] = 1e10;
#if IMPULSE_N >= 2
        L_X[(pos1 + pos2) % SIZE] = 1e10;
#endif
#if IMPULSE_N >= 3
        L_X[(pos1 + 2*pos2) % SIZE] = 1e10;
#endif

        /* Mask out future layers */
        for (int i = 0; i < SIZE; i++)
            if (i % CUTS > layer) L_X[i] = -1e10;
        for (int i = layer+1; i < CUTS; i++)
            for (int j = 0; j < SIZE/CUTS; j++)
                L_X[SIZE + i*(SIZE/CUTS) + j] = -1e10;

        Non_uniform_repeate_double(L_X, L_X_perm, Places);
        ARP_Interleave_double(L_X_perm, L_X_perm, OriginalAddress);
        ARP_Interleave_double(La_D, La_D_perm, OriginalAddress);

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

        for (int i = 0; i < SIZE; i++) dec[i] = (L_D[i] > 0) ? 1 : 0;

        /* Reject if decoded bits outside designed layers */
        int bad = 0;
        for (int i = 0; i < SIZE; i++)
            if (dec[i] && i % CUTS > layer) { bad = 1; break; }
        if (bad) continue;

        int dist = Parity_calc(dec, Places, OriginalAddress,
                               number_of_repetition, circular_states,
                               DataMask, ParityMask);
        if (dist > 0 && dist < TOTAL_DISTANCE) {
            pthread_mutex_lock(&lock);
            leave = 1;
            printf("  Rejected: found d=%d < threshold=%d\n", dist, TOTAL_DISTANCE);
            pthread_mutex_unlock(&lock);
            pthread_exit(NULL);
        }
    }
    return NULL;
}


/* =========================================================================
 * Public API
 * ========================================================================= */

int get_max_threads(void)
{
    long n = sysconf(_SC_NPROCESSORS_ONLN);
    return (n > 0) ? (int)n : 4;
}

int Impulse_Method_parallel(int dmin_threshold,
                             int Places[][MAX_SEGMENTS],
                             int OriginalAddress[],
                             int InterleavedAddress[],
                             int number_of_repetition[MAX_SEGMENTS],
                             int circular_states[MAX_SEGMENTS][32],
                             double snr, double R, int layer)
{
    int num_threads = get_max_threads();
    pthread_t      threads[num_threads];
    ThreadData     tdata[num_threads];
    pthread_attr_t attr;

    pthread_mutex_init(&lock, NULL);
    pthread_attr_init(&attr);
    pthread_attr_setstacksize(&attr, 64 * 1024 * 1024);  /* 64 MB stack */

    TOTAL_DISTANCE = dmin_threshold;
    leave = 0;

    for (int rep = 0; rep < IMPULSE_REPEATE && !leave; rep++) {
        for (int i = 0; i < num_threads; i++) {
            tdata[i] = (ThreadData){
                .thread_id          = i,
                .Places             = Places,
                .OriginalAddress    = OriginalAddress,
                .InterleavedAddress = InterleavedAddress,
                .number_of_repetition = number_of_repetition,
                .circular_states    = circular_states,
                .R   = R,
                .snr = snr,
                .layer = layer
            };
            if (pthread_create(&threads[i], &attr, impulse_thread, &tdata[i]) != 0)
                fprintf(stderr, "Failed to create thread %d\n", i);
        }
        for (int i = 0; i < num_threads; i++)
            pthread_join(threads[i], NULL);
    }

    pthread_attr_destroy(&attr);
    pthread_mutex_destroy(&lock);
    return leave ? -1 : 1;
}


void Build_Multiple_Interleaver_Iteratively(
    int Places[][MAX_SEGMENTS],
    int number_of_repetition[MAX_SEGMENTS],
    int circular_states[MAX_SEGMENTS][32],
    double R, int threshold)
{
    int small_fs = SIZE / CUTS;
    int i, j, l;

    int **OriginalAddress    = malloc(CUTS * sizeof(int *));
    int **InterleavedAddress = malloc(CUTS * sizeof(int *));
    for (i = 0; i < CUTS; i++) {
        OriginalAddress[i]    = malloc(small_fs * sizeof(int));
        InterleavedAddress[i] = malloc(small_fs * sizeof(int));
    }

    /* Flat full-frame maps (needed by Impulse_Method_parallel) */
    int *full_orig = malloc(FRAME_SIZE * sizeof(int));
    int *full_int  = malloc(FRAME_SIZE * sizeof(int));

    int **Total_S = malloc(CUTS * sizeof(int *));
    for (i = 0; i < CUTS; i++) Total_S[i] = malloc(ARP_Q * sizeof(int));
    int *Total_P  = malloc(CUTS * sizeof(int));
    int *S        = malloc(ARP_Q * sizeof(int));

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
        if (SHUFFLE) shuffleArray(results, count);
        P = results[rand() % count];

        /* Random ARP candidate */
        generateRandomArpInterleaver(S, &P, ARP_Q, small_fs);
        ARPInterleaverInitialization(small_fs, P, ARP_Q, S,
                                     InterleavedAddress[l], OriginalAddress[l]);

        /* Build flat full-frame map for this partial set of cuts */
        for (i = 0; i < SIZE; i++) full_orig[i] = i;
        for (i = 0; i < CUTS; i++) {
            if (i <= l) {
                for (j = 0; j < small_fs; j++)
                    full_orig[SIZE + i*small_fs + j] = OriginalAddress[i][j] + SIZE + i*small_fs;
            } else {
                for (j = 0; j < small_fs; j++)
                    full_orig[SIZE + i*small_fs + j] = SIZE + i*small_fs + j;
            }
        }
        for (i = 0; i < FRAME_SIZE; i++)
            full_int[full_orig[i]] = i;

        int sol = Impulse_Method_parallel(threshold, Places, full_orig, full_int,
                                          number_of_repetition, circular_states,
                                          IMPULSE_SNR, R, l);
        if (sol == -1) goto repeate;

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

    /* Write result */
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
    free(full_orig); free(full_int);
    free(Total_S); free(Total_P); free(S);
}
