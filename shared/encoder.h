/**
 * @file encoder.h
 * @brief DPCCC RSC encoder and circular state initialisation.
 *
 * The DPCCC encoder encodes a frame of length FRAME_SIZE as a sequence
 * of independently tail-biting (circular) RSC segments:
 *   - Segment 0: data block, length SIZE
 *   - Segments 1..CUTS: parity sub-blocks, each length SIZE/CUTS
 *
 * All segments share the same RSC generator polynomials G1/G2.
 * Each segment has its own set of circular starting states, stored in
 * circular_states[segment_index][state].
 *
 * Global state:
 *   G1[], G2[]         — generator polynomials (set from config.h macros)
 *   CIRCULAR_STATES[]  — circular states for the global small-frame encoder
 */

#ifndef ENCODER_H
#define ENCODER_H

#include "config.h"

/* =========================================================================
 * Global encoder tables
 * ========================================================================= */

extern int G1[4];
extern int G2[4];

/**
 * CIRCULAR_STATES[CUTS+1][32]: circular starting states for the global
 * encoder instance. Populated by Initialize_States().
 */
extern int CIRCULAR_STATES[10][32];


/* =========================================================================
 * Circular state initialisation
 * ========================================================================= */

/**
 * Find and store the circular (tail-biting) starting states for a
 * single RSC encoder of frame length `framesize`.
 *
 * @param circularstates  Output: circularstates[end_state] = start_state.
 *                        Must have at least 2^DELAYS entries.
 * @param framesize       Length of the segment to be tail-biting encoded.
 */
void Initialize_States(int circularstates[], int framesize);


/* =========================================================================
 * Encoding
 * ========================================================================= */

/**
 * Encode the full DPCCC frame.
 *
 * The input frame of length FRAME_SIZE is split into segments according to
 * number_of_repetition[]. Each segment is tail-biting encoded independently.
 *
 * @param input              Binary input frame, length FRAME_SIZE.
 * @param parity             Output parity bits, length FRAME_SIZE.
 * @param circularstates     Per-segment circular states,
 *                           circularstates[seg][2^DELAYS].
 * @param number_of_repetition  Segment lengths; zero-terminated.
 */
void Encode(int input[], int parity[],
            int circularstates[MAX_SEGMENTS][32],
            int number_of_repetition[]);

/**
 * Encode a single segment of length `framesize` using tail-biting.
 * Used for sub-frame initialisation.
 *
 * @param input         Binary input, length `framesize`.
 * @param parity        Output parity bits, length `framesize`.
 * @param circularstates  Circular states for this segment.
 * @param framesize     Segment length.
 */
void Encodes(int input[], int parity[],
             int circularstates[], int framesize);

/**
 * Encode `size` bits starting from prevstate[], writing to parity[].
 *
 * @param prevstate  Initial encoder state (length `delays`).
 * @param input      Input bit array, length `size`.
 * @param size       Number of bits to encode.
 * @param g1         Feedforward polynomial.
 * @param g2         Feedback polynomial.
 * @param parity     Output parity bits.
 * @param nextstate  Final encoder state.
 * @param delays     Number of memory elements.
 */
void Conv_Encode_number_of_delays(int prevstate[], int input[], int size,
                                   int g1[], int g2[], int parity[],
                                   int nextstate[], int delays);

/**
 * Encode one bit and return the parity output.
 *
 * @param prevstate   Current encoder state.
 * @param input       Input bit (0 or 1).
 * @param g1          Feedforward polynomial.
 * @param g2          Feedback polynomial.
 * @param nextstate   Next encoder state.
 * @param num_delays  Number of memory elements.
 * @return            Parity bit.
 */
int Encode_1_bit_number_of_delays(int prevstate[], int input,
                                   int g1[], int g2[],
                                   int nextstate[], int num_delays);

/**
 * Build the 4-column trellis table for the RSC encoder.
 *
 * trellis[state] = {parity_out_0, next_state_0, parity_out_1, next_state_1}
 *
 * @param g1      Feedforward polynomial.
 * @param g2      Feedback polynomial.
 * @param delays  Number of memory elements.
 * @param trellis Output table, size [2^delays][4].
 */
void DrawTrellis(int g1[], int g2[], int delays, int trellis[][4]);


/* =========================================================================
 * Puncturing masks
 * ========================================================================= */

/**
 * Fill `x[0..framesize-1]` with the puncturing mask of type `type`.
 *
 * @param x          Output mask array.
 * @param framesize  Length of mask to fill.
 * @param type       0 = data, 1 = parity-K, 2 = parity-K/Q.
 */
void puncmask_all_frame_FrameSize(int x[], int framesize, int type);


/* =========================================================================
 * Distance evaluation
 * ========================================================================= */

/**
 * Compute the Hamming weight of the DPCCC codeword for input pattern x[].
 *
 * The input pattern is first expanded via the DPCCC repetition mapping
 * (Non_uniform_repeate), then interleaved and encoded.
 *
 * @param x                Positions of 1s in the K-bit information block.
 * @param size             Input weight.
 * @param Original_address Full-frame ARP interleaver forward map.
 * @param Places           DPCCC repetition mapping table.
 * @param circular_states  Per-segment circular states.
 * @param number_of_repetition  Segment lengths.
 * @return                 Total codeword Hamming weight.
 */
int Turbo_distance_circular_With_Puncturing(int x[], int size,
                                             int **Original_address);

/**
 * Compute the punctured parity weight for input x[] on one RSC segment.
 *
 * @param x          Positions of 1s in the segment.
 * @param size       Input weight.
 * @param type       Puncturing mask type (0, 1, or 2).
 * @param frame_size Segment length.
 * @return           Punctured parity weight.
 */
int parity_distance_circular_With_Puncturing(int x[], int size,
                                              int type, int frame_size);

/**
 * Compute the full codeword weight of a decoded hard-decision output.
 * Used by the impulse method to check if a decoded word is a valid
 * low-weight codeword.
 */
int Parity_calc(int output[],
                int Places[][MAX_SEGMENTS],
                int OriginalAddress[],
                int number_of_repetition[MAX_SEGMENTS],
                int circular_states[MAX_SEGMENTS][32],
                int DataMask[], int Parity_Mask[]);

#endif /* ENCODER_H */
