/**
 * @file sim_config.h
 * @brief Configuration specific to the DPCCC Monte Carlo BER/FER simulator.
 *
 * Edit the values below before compiling montecarlo/.
 * The interleaver parameters (S, P) used by ARP_SEMI_Initialization()
 * are hardcoded in interleaver.c — update them there to test your
 * designed interleavers.
 */

#ifndef SIM_CONFIG_H
#define SIM_CONFIG_H

/** Termination mode for the DPCCC decoder:
 *   0 = full frame, single encoder (legacy)
 *   1 = K + K + K/Q + ... (non-uniform repetition, non-special coupling)
 *   2 = K + K/CUTS repeated CUTS times (special coupling — used in paper)
 */
#define SEMI_CIRCULAR_TERMINATION  2

/** Fixed RNG seed for reproducible runs. Set to 0 to use time(NULL). */
#define RNG_SEED  12345

/** Target number of frame errors at each SNR point before moving on. */
#define TARGET_FRAME_ERRORS  250

#endif /* SIM_CONFIG_H */
