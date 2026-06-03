/**
 * @file designer.h
 * @brief N-impulse DPCCC interleaver design using the decoder.
 *
 * The N-impulse method evaluates an ARP interleaver candidate by:
 *   1. Encoding the all-zero codeword.
 *   2. Adding AWGN noise.
 *   3. Pinning IMPULSE_N decoder inputs to a very high LLR value
 *      (forcing the decoder to "believe" those bits are 1).
 *   4. Running the turbo decoder for MAX_ITERATIONS iterations.
 *   5. Checking if the decoded word is a valid codeword with
 *      Hamming weight < DESIGN_THRESHOLD.
 *
 * If such a word is found, the interleaver is rejected.
 * The method is parallelised across CPU cores using pthreads.
 *
 * The outer loop runs over all cuts l = 0 … CUTS-1 independently.
 * For each cut, the ARP parameters are chosen randomly from the set
 * of valid periods, and the sub-layer shift vector S[] is built using
 * the layer-by-layer search with the impulse evaluation as the oracle.
 *
 * Output: result.txt (same format as decoder_free tool)
 */

#ifndef DESIGNER_H
#define DESIGNER_H

#include "../../../shared/config.h"

/**
 * Run CUTS rounds of N-impulse design and write results.
 *
 * @param Places               DPCCC repetition table.
 * @param number_of_repetition Segment lengths.
 * @param circular_states      Per-segment circular states.
 * @param R                    Code rate.
 * @param threshold            Minimum mHD target.
 */
void Build_Multiple_Interleaver_Iteratively(
    int Places[][MAX_SEGMENTS],
    int number_of_repetition[MAX_SEGMENTS],
    int circular_states[MAX_SEGMENTS][32],
    double R, int threshold);

/**
 * Multi-threaded N-impulse evaluation of one interleaver candidate.
 *
 * Spawns get_max_threads() threads, each testing a different pair of
 * impulse positions.
 *
 * @param dmin_threshold       Distance threshold for rejection.
 * @param Places               DPCCC repetition table.
 * @param OriginalAddress      Full-frame ARP forward map.
 * @param InterleavedAddress   Full-frame ARP inverse map.
 * @param number_of_repetition Segment lengths.
 * @param circular_states      Per-segment circular states.
 * @param snr                  Eb/N0 (dB) for the AWGN channel.
 * @param R                    Code rate.
 * @param layer                Current cut index (controls which LLRs are masked).
 * @return                     -1 if a codeword below threshold was found, 1 otherwise.
 */
int Impulse_Method_parallel(int dmin_threshold,
                             int Places[][MAX_SEGMENTS],
                             int OriginalAddress[],
                             int InterleavedAddress[],
                             int number_of_repetition[MAX_SEGMENTS],
                             int circular_states[MAX_SEGMENTS][32],
                             double snr, double R, int layer);

/** Returns the number of online CPU cores. */
int get_max_threads(void);

/** Global flag: 1 when a low-weight codeword has been found. */
extern int leave;

/** Global distance used as the rejection threshold across threads. */
extern int TOTAL_DISTANCE;

#endif /* DESIGNER_H */
