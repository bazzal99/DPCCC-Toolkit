/**
 * @file interleaver.h
 * @brief DPCCC ARP interleaver and repetition mapping.
 *
 * The DPCCC uses a two-level interleaver structure:
 *
 *  1. Non-uniform repetition (Non_uniform_repeate):
 *     Maps the K-bit information block to the FRAME_SIZE expanded frame,
 *     placing each bit according to the Places[][] table.
 *
 *  2. ARP permutation (ARP_Interleave / ARP_DE_Interleave):
 *     Permutes the expanded frame using the ARP forward map OriginalAddress[].
 *
 * The ARP map is π(i) = (P * i + S[i mod Q]) mod segment_size,
 * applied independently to each of the CUTS parity sub-block segments.
 * OriginalAddress[][] is therefore a 2D array: [CUTS][SIZE/CUTS].
 */

#ifndef INTERLEAVER_H
#define INTERLEAVER_H

#include "config.h"

/* =========================================================================
 * Non-uniform repetition mapping
 * ========================================================================= */

/**
 * Expand the K-bit information block `input` to the full frame `new_input`
 * using the Places[][] repetition table.
 *
 * Places[i][j] = position in the expanded frame where the j-th copy of
 * bit i is placed. -1 terminates the list for bit i.
 *
 * @param input      Information block, length SIZE.
 * @param new_input  Expanded frame, length FRAME_SIZE.
 * @param Places     Repetition table, size [SIZE][MAX_SEGMENTS].
 */
void Non_uniform_repeate(int input[], int new_input[],
                          int Places[][MAX_SEGMENTS]);

/** Double-precision version of Non_uniform_repeate (for LLR processing). */
void Non_uniform_repeate_double(double input[], double new_input[],
                                 int Places[][MAX_SEGMENTS]);


/* =========================================================================
 * ARP interleaver (full-frame, integer and double variants)
 * ========================================================================= */

/**
 * Apply the ARP permutation to `data[]`, writing result to `output[]`.
 * output[OriginalAddress[i]] = data[i] for all i in [0, FRAME_SIZE).
 */
void ARP_Interleave(int data[], int output[], int OriginalAddress[]);
void ARP_Interleave_double(double data[], double output[], int OriginalAddress[]);

/**
 * Apply the inverse ARP permutation.
 * output[InterleavedAddress[i]] = data[i] for all i in [0, FRAME_SIZE).
 */
void ARP_DE_Interleave(int data[], int output[], int InterleavedAddress[]);
void ARP_DE_Interleave_double(double data[], double output[], int InterleavedAddress[]);


/* =========================================================================
 * Per-cut ARP initialisation (for decoder_free design tool)
 * ========================================================================= */

/**
 * Fully initialise one ARP sub-interleaver for segment `cut`.
 *
 * OriginalAddress[j] = (P * j + S[j mod Q]) mod segment_size
 *
 * @param Size       Segment length = SIZE/CUTS.
 * @param Period     ARP period P.
 * @param PeriodARP  Sub-layer period Q.
 * @param Shift      Shift vector S[], length Q.
 * @param InterleavedAddress  Output: inverse map.
 * @param OriginalAddress     Output: forward map.
 */
void ARPInterleaverInitialization(int Size, int Period, int PeriodARP,
                                  int Shift[],
                                  int InterleavedAddress[],
                                  int OriginalAddress[]);

/** Initialise sub-layers 0 … layer only (partial initialisation). */
void ARPInterleaverInitialization_Till_Layer(int Size, int layer,
                                             int Period, int Q,
                                             int Shift[],
                                             int InterleavedAddress[],
                                             int OriginalAddress[]);

/** Initialise only sub-layer `layer`. */
void ARPInterleaverInitialization_Per_Layer(int Size, int layer,
                                            int Period, int Q,
                                            int Shift[],
                                            int InterleavedAddress[],
                                            int OriginalAddress[]);

/** Undo (de-initialise) sub-layer `layer`. */
void Arp_DE_Interleave_Layer(int Size, int layer, int Period, int Q,
                             int Shift[],
                             int InterleavedAddress[],
                             int OriginalAddress[]);


/* =========================================================================
 * Full-frame ARP initialisation (for montecarlo)
 * ========================================================================= */

/**
 * Initialise the complete DPCCC interleaver from hardcoded interleaver
 * parameters found in the paper.
 *
 * Fills OriginalAddress[FRAME_SIZE] and InterleavedAddress[FRAME_SIZE]
 * and sets number_of_repetition[].
 */
void ARP_SEMI_Initialization(int OriginalAddress[], int InterleavedAddress[],
                              int number_of_repetition[MAX_SEGMENTS]);

/**
 * S-random interleaver initialisation (alternative to ARP for testing).
 */
void S_Random_Interleaver_Initialization(int s, int OriginalAddress[],
                                          int InterleavedAddress[],
                                          int number_of_repetition[MAX_SEGMENTS]);

/** Generate an S-random permutation of length `size` with spread parameter `s`. */
void S_Random(int s, int output[], int size);


/* =========================================================================
 * Period selection
 * ========================================================================= */

/** Fill `periods[]` with valid ARP periods for frame size `framesize`. */
void Choose_Period(int periods[], int *size, int framesize);

/** Compute the minimum spread of linear permutation P*i mod K. */
int Min_Spread(int p);

/** GCD via Euclid's algorithm. */
int gcd(int a, int b);

/** Coprimality test. */
int isRelativePrime(int a, int b);

/** Fill `result[]` with all integers coprime with `number` in [1, number]. */
void saveRelativePrimes(int number, int *result, int *count);

/** Evaluate ARP at position i. */
int inter(int i, int p, int s[]);

/** Check if array arr[0..Q-1] has all distinct values. */
int isArrayUnique(int arr[], int Q);

/** Generate a random valid ARP interleaver. */
void generateRandomArpInterleaver(int S[], int *P, int Q, int frame_size);


/* =========================================================================
 * Extrinsic information (used by decoder)
 * ========================================================================= */

/**
 * Compute the extrinsic LLR for bit i at repetition index j.
 * Used during the turbo decoding iteration to pass information
 * between the sub-decoders.
 */
double Calculate_Extrinsic(int Places[][MAX_SEGMENTS], int row, int column,
                            double L_D[], double L_X[], double La_D[]);

#endif /* INTERLEAVER_H */
