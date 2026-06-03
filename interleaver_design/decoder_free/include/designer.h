/**
 * @file designer.h
 * @brief Decoder-free DPCCC interleaver design algorithm.
 *
 * Designs a set of CUTS independent ARP sub-interleavers, one per
 * DPCCC parity sub-block, such that the minimum Hamming distance of
 * the resulting DPCCC code meets or exceeds DESIGN_THRESHOLD.
 *
 * The algorithm proceeds cut by cut (l = 0 … CUTS-1):
 *   For each cut l, run the ARP layer-by-layer search using the
 *   DPCCC-aware RTZ kernel (rtz_search.h), which respects the
 *   partial-encoder structure (only cuts 0..l are active).
 *
 * The main entry point is Build_Multiple_Interleaver_Iteratively().
 *
 * Output file: result.txt
 *   Format per solution:
 *     int total_S[CUTS][Q] = { {s0_0,..,s0_{Q-1}}, ..., {sL_0,..,sL_{Q-1}} };
 *     int total_P[CUTS]    = {P_0, ..., P_{L-1}};
 *     distance = <d>
 */

#ifndef DESIGNER_H
#define DESIGNER_H

#include "../../../shared/config.h"

/**
 * Run CUTS rounds of the ARP design search and write the result.
 *
 * @param Places               DPCCC repetition table [SIZE][MAX_SEGMENTS].
 * @param number_of_repetition Segment lengths; zero-terminated.
 * @param circular_states      Per-segment circular states [MAX_SEGMENTS][32].
 * @param R                    Code rate (for impulse-method evaluation).
 * @param threshold            Minimum mHD target.
 */
void Build_Multiple_Interleaver_Iteratively(
    int Places[][MAX_SEGMENTS],
    int number_of_repetition[MAX_SEGMENTS],
    int circular_states[MAX_SEGMENTS][32],
    double R, int threshold);

/**
 * Recursive layer builder: find a valid shift S[current_layer] for one
 * ARP sub-layer of cut `the_layer`.
 */
void Generate_ARP_Interleaver_Iteratively(int the_layer, int current_layer,
                                           int **OriginalAddress,
                                           int **InterleavedAddress,
                                           int threshold, int *S, int P);

/**
 * Find all valid shift values for ARP sub-layer `current_shift_layer`
 * of cut `the_layer`, storing up to `max_size_layer` candidates in
 * array_s[].
 */
void Add_ARP_LAYER(int the_layer, int current_shift_layer,
                   int **OriginalAddress, int **InterleavedAddress,
                   int threshold, int *S, int P,
                   int *array_s, int *size, int max_size_layer);

/** Global flag: set to 1 when a complete valid interleaver is found. */
extern int StopRecursion;

#endif /* DESIGNER_H */
