/**
 * @file config.h
 * @brief Shared DPCCC encoder and puncturing configuration.
 *
 * This file is included by all three tools in the DPCCC-Toolkit:
 *   - interleaver_design/decoder_free
 *   - interleaver_design/double_impulse
 *   - montecarlo
 *
 * Edit the values here before compiling any tool.
 */

#ifndef CONFIG_H
#define CONFIG_H

#include <math.h>

/* =========================================================================
 * RSC component encoder
 * ========================================================================= */

/** Number of memory elements (delays). Trellis has 2^DELAYS states. */
#define DELAYS  3

/**
 * RSC generator polynomials (binary, coefficient order: g[0]=direct input,
 * g[1..DELAYS]=tap from shift register).
 * Default: G(1, 13/15) in octal — same as LTE.
 * G1 = feedforward (numerator), G2 = feedback (denominator).
 */
#define G1_INIT  {1, 1, 0, 1}
#define G2_INIT  {1, 0, 1, 1}

/**
 * Period of the RSC encoder's free impulse response.
 * For G(1, 13/15): period = 7.
 */
#define ENCODER_PERIOD  7

/**
 * Encoding mode:
 *   1 = Circular / tail-biting (TB) — used for DPCCC sub-blocks
 *   0 = Zero-termination (ZT)
 */
#define CIRCULAR_ENCODING  1

/* =========================================================================
 * DPCCC frame structure
 * =========================================================================
 *
 * The DPCCC frame of length FRAME_SIZE = K + K/CUTS * CUTS consists of:
 *   - One data block of length K = SIZE bits
 *   - CUTS parity sub-blocks, each of length K/CUTS bits
 *
 * MAX_SEGMENTS is the maximum number of separately-encoded segments
 * (data block + parity sub-blocks).
 */

/** Information block length K in bits. */
#define SIZE  1024

/** Number of DPCCC parity sub-blocks (decomposition factor L). */
#define CUTS  2

/** Total frame length = SIZE + CUTS * (SIZE/CUTS) = 2*SIZE. */
#define FRAME_SIZE  (SIZE + CUTS * (SIZE / CUTS))

/** Maximum number of encoded segments (data + parity sub-blocks). */
#define MAX_SEGMENTS  30

/* =========================================================================
 * Puncturing patterns
 * =========================================================================
 *
 * Three independent patterns, each of period PUNCT_MASK:
 *   type 0: data bits
 *   type 1: parity bits of the first RSC (data block)
 *   type 2: parity bits of the sub-block RSCs
 *
 * 1 = transmit, 0 = puncture.
 */

/** Puncturing mask period. */
#define PUNCT_MASK  1

#define DATA_PUNCT_INIT    {1}
#define PARITY_K_INIT      {1}
#define PARITY_K_Q_INIT    {1}

/* =========================================================================
 * ARP interleaver
 * ========================================================================= */

/** Number of ARP sub-layers Q. Must divide SIZE. */
#define ARP_Q  4

/* =========================================================================
 * SNR sweep (used by montecarlo)
 * ========================================================================= */

#define SNR_MIN   0.0
#define SNR_MAX   8.0
#define SNR_STEP  0.25

/** Maximum number of decoding iterations. */
#define MAX_ITERATIONS  20

/* =========================================================================
 * max_star approximation
 *
 * Define FAST_MAX_STAR to use a two-piece linear approximation of
 * log(1 + exp(-|x-y|)).  Undefine (default) for exact log-MAP.
 * ========================================================================= */
/* #define FAST_MAX_STAR */

#endif /* CONFIG_H */
