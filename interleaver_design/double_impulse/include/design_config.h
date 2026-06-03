/**
 * @file design_config.h
 * @brief Configuration specific to the N-impulse DPCCC interleaver design tool.
 *
 * The N-impulse method estimates the minimum Hamming distance of a
 * DPCCC code by injecting N simultaneous impulses into the decoder
 * (using very high LLR values) and checking whether the decoder
 * converges to a low-weight codeword.
 *
 * Setting IMPULSE_N = 1 gives the single-impulse method.
 * Setting IMPULSE_N = 2 gives the double-impulse method.
 * Setting IMPULSE_N = 3 gives the triple-impulse method.
 *
 * Edit the values below before compiling interleaver_design/double_impulse.
 */

#ifndef DESIGN_CONFIG_H
#define DESIGN_CONFIG_H

/** Minimum acceptable mHD. */
#define DESIGN_THRESHOLD  40

/** Number of simultaneous impulses to inject (1, 2, or 3). */
#define IMPULSE_N  2

/** Number of impulse-method repetitions per thread per candidate. */
#define IMPULSE_REPEATE  1

/** Number of shift candidates per ARP sub-layer (search breadth). */
#define SHIFTS_PER_LAYER  1

/** Randomise the search order (recommended). */
#define SHUFFLE  1

/**
 * Track multiplicity of the minimum distance.
 *   0 = distance only (faster)
 *   1 = count multiplicities
 */
#define WITH_MULTIPLICITY  0

/** Print progress to stdout. */
#define VERBOSE  1

/** SNR used for the impulse-method decoder. */
#define IMPULSE_SNR  5.75

#endif /* DESIGN_CONFIG_H */
