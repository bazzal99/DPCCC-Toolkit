/**
 * @file decoder.h
 * @brief Log-MAP (BCJR) decoder for one RSC trellis segment.
 *
 * Implements the forward-backward algorithm in the log domain using
 * the max-star (Jacobi logarithm) approximation.
 * Set FAST_MAX_STAR in config.h for the two-piece linear approximation.
 */
#ifndef DECODER_H
#define DECODER_H
#include "config.h"

/**
 * Full log-MAP decode of one segment.
 * Computes log_gamma, log_alpha, log_beta in one call and produces
 * the a-posteriori LLR L_D[].
 *
 * Alpha_Initiale_value[] and Beta_Initiale_value[] are updated in place
 * to serve as initial conditions for the next iteration (tail-biting decoding).
 */
void decode(double La_D[], double L_X[], double L_Y[], double L_D[],
            int num_states, int trellis[][4],
            double log_alpha[][FRAME_SIZE+1],
            double log_beta[][FRAME_SIZE+1],
            double log_gamma[][2][FRAME_SIZE],
            double Alpha_Initiale_value[], double Beta_Initiale_value[],
            int frame_size);

void compute_log_gamma(double La_D[], double L_X[], double L_Y[],
                        int num_states, int trellis[][4],
                        double log_gamma[][2][FRAME_SIZE], int frame_size);

void compute_log_alpha(int numstates, int trellis[][4],
                        double log_alpha[][FRAME_SIZE+1],
                        double log_gamma[][2][FRAME_SIZE],
                        double Alpha_Initiale_value[], int frame_size);

void compute_log_beta(int numstates, int trellis[][4],
                       double log_beta[][FRAME_SIZE+1],
                       double log_gamma[][2][FRAME_SIZE],
                       double Beta_Initiale_value[], int frame_size);

/** max*(x,y) = max(x,y) + log(1 + exp(-|x-y|)) */
double max_star(double x, double y);

#endif /* DECODER_H */
