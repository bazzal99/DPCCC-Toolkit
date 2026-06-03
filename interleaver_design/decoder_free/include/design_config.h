/**
 * @file design_config.h
 * @brief Configuration specific to the decoder-free DPCCC interleaver design tool.
 *
 * Edit the values below before compiling interleaver_design/decoder_free.
 */

#ifndef DESIGN_CONFIG_H
#define DESIGN_CONFIG_H

/* =========================================================================
 * Design search parameters
 * ========================================================================= */

/** Minimum acceptable mHD. The search rejects interleavers below this. */
#define DESIGN_THRESHOLD  46

/** Number of shift candidates to keep per ARP sub-layer (search breadth). */
#define SHIFTS_PER_LAYER  2

/** Randomise the search order for periods and shift values (recommended). */
#define SHUFFLE  1

/**
 * Track and report the multiplicity of the minimum distance.
 *   0 = report distance only (faster)
 *   1 = count multiplicities (slightly slower, uses repeates.txt)
 */
#define WITH_MULTIPLICITY  0

/** Number of impulse-method parallel repetitions per thread. */
#define IMPULSE_REPEATE  1

/** Print search progress to stdout. */
#define VERBOSE  1

/* =========================================================================
 * RTZ search hyperparameters (circular mode)
 * ========================================================================= */

/** Girth bounds for the periodic (weight-2, jump=ENCODER_PERIOD) kernel. */
#define IW2_GIRTH_BOUNDS   {2, 4, 6, 8}

/** Period budgets for the periodic (weight-2) kernel. */
#define IW2_PERIOD_BUDGETS {10, 10, 7, 4}

/** Maximum weight tested by the periodic kernel. */
#define IW2_MAX_WEIGHT  4

/** Girth bounds for the general (weight-3+, jump=1) kernel. */
#define IWM_GIRTH_BOUNDS   {0, 1, 2, 3, 4, 5, 6, 7, 8, 9}

/** Step budgets for the general kernel. */
#define IWM_STEP_BUDGETS   {-1, -1, -1, 50, 28, 15, 10, 7, 5, 4}

/** Maximum weight tested by the general kernel. */
#define IWM_MAX_WEIGHT  9

#endif /* DESIGN_CONFIG_H */
