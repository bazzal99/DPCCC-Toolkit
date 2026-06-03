/**
 * @file config.h
 * @brief Shared DPCCC encoder and frame configuration.
 *
 * Edit this file before compiling any of the three tools.
 * All tools share this single configuration file.
 *
 * Paper reference:
 *   M. Bazzal, J. Nadal, S. Weithoffer, C. A. Nour, C. Douillard,
 *   "A Novel Parallel Concatenated Convolutional Code Structure Based
 *   on Frame Decomposition," ISTC 2025.
 *   DOI: 10.1109/ISTC65386.2025.11154568
 */

#ifndef CONFIG_H
#define CONFIG_H

#include <math.h>

/* =========================================================================
 * RSC component encoder
 * =========================================================================
 *
 * Two encoder configurations from the paper:
 *
 *   G1 = G(1, 15/13)₈  — single parity, no puncturing
 *        G1_INIT = {1, 1, 0, 1}  (feedforward: 1 + D + D³)
 *        G2_INIT = {1, 0, 1, 1}  (feedback:    1 + D² + D³)
 *
 *   G2 = G(1, 15/13, 17/13)₈  — double parity, with systematic puncturing
 *        Same G1/G2 polynomials, but add a second feedforward polynomial
 *        and apply puncturing masks below.
 *
 * Default: G1 configuration (matches most results in Table I).
 */

/** Number of memory elements. Trellis has 2^DELAYS states. */
#define DELAYS  3

/** RSC feedforward polynomial G1 (binary coefficients, g[0]=input tap). */
#define G1_INIT  {1, 1, 0, 1}

/** RSC feedback polynomial G2 (denominator). */
#define G2_INIT  {1, 0, 1, 1}

/**
 * Encoder period p.
 * For G(1, 15/13)₈: period = 7.
 * IMPORTANT: choose L (= CUTS) such that L mod ENCODER_PERIOD != 0.
 * (See Section IV-B of the paper — L mod p = 0 degrades performance.)
 */
#define ENCODER_PERIOD  7

/**
 * Encoding mode:
 *   1 = Circular / tail-biting (TB)  — used for all DPCCC sub-blocks
 *   0 = Zero-termination (ZT)
 */
#define CIRCULAR_ENCODING  1


/* =========================================================================
 * DPCCC frame structure (Section III of the paper)
 * =========================================================================
 *
 *  FRAME_SIZE = SIZE + CUTS * (SIZE / CUTS)
 *
 *  The frame consists of:
 *    - Data block:          SIZE bits
 *    - Parity sub-blocks:   CUTS * (SIZE/CUTS) = SIZE bits
 *
 *  Overall code rate (no puncturing) = SIZE / (SIZE + SIZE) = 1/2
 *  To get rate 1/3: add the full-length first parity stream (see encoder.c).
 *
 *  Recommended L (= CUTS) values from Table I of the paper:
 *    K = 128:   L = 2  (mHD = 28 for G1)
 *    K = 1024:  L = 2 or L = 4  (mHD = 46 for G1)
 *    K = 6144:  L = 8  (mHD = 57 for G1, 72 for G2)
 *
 *  Always verify: CUTS mod ENCODER_PERIOD != 0
 *    L=2: 2 mod 7 = 2 ✓
 *    L=4: 4 mod 7 = 4 ✓
 *    L=8: 8 mod 7 = 1 ✓
 *    L=7: 7 mod 7 = 0 ✗  (never use)
 *    L=14: 14 mod 7 = 0 ✗ (never use)
 */

/** Information block length K in bits. */
#define SIZE  1024

/**
 * Number of DPCCC parity sub-blocks L (decomposition factor).
 * Each sub-block has length SIZE/CUTS bits.
 */
#define CUTS  2

/** Total encoded frame length = SIZE (data) + SIZE (all parities). */
#define FRAME_SIZE  (SIZE + CUTS * (SIZE / CUTS))

/** Maximum number of independently-encoded segments. */
#define MAX_SEGMENTS  30


/* =========================================================================
 * Puncturing patterns
 * =========================================================================
 *
 * For G1 (no puncturing): all masks = {1} (keep all bits).
 * For G2 (from paper, Table I footnote):
 *   data puncturing:    [10100000]
 *   first parity:       [11101111]
 *   sub-block parity:   [01011010]
 *
 * Uncomment the G2 masks to switch to the G2 configuration.
 */

/** Puncturing mask period (bits per puncturing cycle). */
#define PUNCT_MASK  1

/* --- G1: no puncturing (default) --- */
#define DATA_PUNCT_INIT    {1}
#define PARITY_K_INIT      {1}
#define PARITY_K_Q_INIT    {1}

/* --- G2: with systematic puncturing (uncomment for G2) ---
#undef  PUNCT_MASK
#define PUNCT_MASK  8
#define DATA_PUNCT_INIT    {1, 0, 1, 0, 0, 0, 0, 0}
#define PARITY_K_INIT      {1, 1, 1, 0, 1, 1, 1, 1}
#define PARITY_K_Q_INIT    {0, 1, 0, 1, 1, 0, 1, 0}
*/


/* =========================================================================
 * ARP interleaver
 * ========================================================================= */

/**
 * Number of ARP sub-layers Q.
 * Must divide SIZE/CUTS (the sub-block length).
 */
#define ARP_Q  4

/**
 * Minimum spread offset for period selection.
 * A period P is accepted if Min_Spread(P) >= sqrt(2*K/CUTS) - SPREAD_OFFSET.
 */
#define SPREAD_OFFSET  2


/* =========================================================================
 * Monte Carlo simulation (montecarlo tool)
 * ========================================================================= */

/** Minimum Eb/N0 (dB) for the SNR sweep. */
#define SNR_MIN   0.0

/** Maximum Eb/N0 (dB) for the SNR sweep. */
#define SNR_MAX   5.0

/** Eb/N0 step size (dB). */
#define SNR_STEP  0.25

/**
 * Maximum number of log-MAP decoding iterations.
 * Paper uses 8 iterations for all DPCCC and LTE simulations.
 */
#define MAX_ITERATIONS  8

/**
 * max* approximation in the log-MAP decoder.
 * Uncomment to use the two-piece linear approximation (max-log-MAP).
 * Default: exact log-MAP.
 */
/* #define FAST_MAX_STAR */

#endif /* CONFIG_H */
