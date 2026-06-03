/**
 * @file decoder.c
 * @brief Log-MAP (BCJR) decoder implementation.
 */
#include <math.h>
#include "config.h"
#include "decoder.h"

double max_star(double x, double y)
{
#ifdef FAST_MAX_STAR
    double xmy = fabs(x - y);
    double z   = (xmy > 2.5) ? 0.0 : 0.6 - 0.24 * xmy;
    return (x > y) ? x + z : y + z;
#else
    double max_xy = (x > y) ? x : y;
    return max_xy + log(1.0 + exp(-fabs(x - y)));
#endif
}

void compute_log_gamma(double La_D[], double L_X[], double L_Y[],
                        int num_states, int trellis[][4],
                        double log_gamma[][2][FRAME_SIZE], int frame_size)
{
    for (int i = 0; i < frame_size; i++)
        for (int j = 0; j < num_states; j++) {
            double x0 = 2 * trellis[j][0] - 1;
            double x1 = 2 * trellis[j][2] - 1;
            log_gamma[j][0][i] = (-La_D[i] - L_X[i] + x0 * L_Y[i]) / 2.0;
            log_gamma[j][1][i] = ( La_D[i] + L_X[i] + x1 * L_Y[i]) / 2.0;
        }
}

void compute_log_alpha(int numstates, int trellis[][4],
                        double log_alpha[][FRAME_SIZE+1],
                        double log_gamma[][2][FRAME_SIZE],
                        double Alpha_Initiale_value[], int frame_size)
{
    for (int i = 0; i < numstates; i++)
        log_alpha[i][0] = Alpha_Initiale_value[i];

    for (int i = 1; i < frame_size + 1; i++) {
        for (int j = 0; j < numstates; j++) {
            int s0 = -1, s1 = -1;
            for (int k = 0; k < numstates; k++) {
                if (trellis[k][1] == j) s0 = k;
                if (trellis[k][3] == j) s1 = k;
            }
            log_alpha[j][i] = max_star(
                log_alpha[s0][i-1] + log_gamma[s0][0][i-1],
                log_alpha[s1][i-1] + log_gamma[s1][1][i-1]);
        }
        double norm = log_alpha[0][i];
        for (int j = 1; j < numstates; j++) norm = max_star(log_alpha[j][i], norm);
        for (int j = 0; j < numstates; j++) log_alpha[j][i] -= norm;
    }
}

void compute_log_beta(int numstates, int trellis[][4],
                       double log_beta[][FRAME_SIZE+1],
                       double log_gamma[][2][FRAME_SIZE],
                       double Beta_Initiale_value[], int frame_size)
{
    for (int i = 0; i < numstates; i++)
        log_beta[i][frame_size] = Beta_Initiale_value[i];

    for (int i = frame_size - 1; i >= 0; i--) {
        for (int j = 0; j < numstates; j++)
            log_beta[j][i] = max_star(
                log_gamma[j][0][i] + log_beta[trellis[j][1]][i+1],
                log_gamma[j][1][i] + log_beta[trellis[j][3]][i+1]);
        double norm = log_beta[0][i];
        for (int j = 1; j < numstates; j++) norm = max_star(log_beta[j][i], norm);
        for (int j = 0; j < numstates; j++) log_beta[j][i] -= norm;
    }
}

void decode(double La_D[], double L_X[], double L_Y[], double L_D[],
            int num_states, int trellis[][4],
            double log_alpha[][FRAME_SIZE+1],
            double log_beta[][FRAME_SIZE+1],
            double log_gamma[][2][FRAME_SIZE],
            double Alpha_Initiale_value[], double Beta_Initiale_value[],
            int frame_size)
{
    compute_log_gamma(La_D, L_X, L_Y, num_states, trellis, log_gamma, frame_size);
    compute_log_alpha(num_states, trellis, log_alpha, log_gamma, Alpha_Initiale_value, frame_size);
    compute_log_beta(num_states, trellis, log_beta, log_gamma, Beta_Initiale_value, frame_size);

    for (int i = 0; i < frame_size; i++) {
        double llr1 = log_alpha[0][i] + log_gamma[0][1][i] + log_beta[trellis[0][3]][i+1];
        for (int j = 1; j < num_states; j++)
            llr1 = max_star(log_alpha[j][i] + log_gamma[j][1][i] + log_beta[trellis[j][3]][i+1], llr1);
        double llr0 = log_alpha[0][i] + log_gamma[0][0][i] + log_beta[trellis[0][1]][i+1];
        for (int j = 1; j < num_states; j++)
            llr0 = max_star(log_alpha[j][i] + log_gamma[j][0][i] + log_beta[trellis[j][1]][i+1], llr0);
        L_D[i] = llr1 - llr0;
    }
    for (int i = 0; i < num_states; i++) {
        Alpha_Initiale_value[i] = log_alpha[i][frame_size];
        Beta_Initiale_value[i]  = log_beta[i][0];
    }
}
