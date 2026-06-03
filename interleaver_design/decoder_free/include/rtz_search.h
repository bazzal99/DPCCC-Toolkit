/**
 * @file rtz_search.h
 * @brief DPCCC-aware RTZ graph-girth search kernel for decoder-free design.
 *
 * This kernel differs from the standard ARP interleaver search in that it
 * operates on the DPCCC correlation graph, where:
 *   - Each position i in the K-bit block has a CUTS-dimensional label:
 *       layer = i % CUTS,  position = floor(i / CUTS)
 *   - The search respects the DPCCC structure: it only traverses positions
 *     whose layer index is <= the currently designed cut (the_layer).
 *   - Per-cut ARP interleavers are stored in OriginalAddress[CUTS][SIZE/CUTS]
 *     and InterleavedAddress[CUTS][SIZE/CUTS].
 *
 * Two jump modes are supported, controlled by the `jump_type` parameter:
 *   jump_type 0 or 1: unit steps (general RTZ, iw=3+) going forward
 *   jump_type 2 or 3: steps of `jump` = ENCODER_PERIOD (periodic RTZ, iw=2)
 *
 * These are combined in the unified kernel
 * Min_Girth_One_Layer_JUMP_1_7_Iteratively_ARP_INTERLEAVER().
 */

#ifndef RTZ_SEARCH_H
#define RTZ_SEARCH_H

/**
 * Top-level entry: search for RTZ cycles starting from `position` in the
 * DPCCC interleaver, for the cut being designed (`current_shift_layer`).
 *
 * @param current_shift_layer  Index of the cut whose shift is being chosen.
 * @param position             Starting position in [0, SIZE).
 * @param layer                Layer of the starting position (position % CUTS).
 * @param OriginalAddress      Per-cut forward maps, size [CUTS][SIZE/CUTS].
 * @param InterleavedAddress   Per-cut inverse maps, size [CUTS][SIZE/CUTS].
 * @param threshold            Current distance threshold.
 * @param mult                 Multiplicity counter (updated in place).
 * @param max_RTZ_size_start_0 Step budget for the original-domain traversal.
 * @param max_RTZ_size_start_1 Step budget for the interleaved-domain traversal.
 * @param mingirth             Maximum cycle length to explore.
 * @param jump_type            0/1 = unit steps, 2/3 = periodic (ENCODER_PERIOD) steps.
 * @param jump                 Step size (1 or ENCODER_PERIOD).
 * @return                     0 if a distance < threshold was found, 1 otherwise.
 */
int Min_Girth_One_Layer_JUMP_1_7_Iteratively_ARP_INTERLEAVER(
    int current_shift_layer, int position, int layer,
    int **OriginalAddress, int **InterleavedAddress,
    int threshold, int *mult,
    int max_RTZ_size_start_0, int max_RTZ_size_start_1,
    int mingirth, int jump_type, int jump);

/** Recursive implementation of the above. */
int Min_Girth_One_Ref_JUMP_1_7_Iteratively_ARP_INTERLEAVER(
    int current_shift_layer, int the_layer, int my_layer,
    int **OriginalAddress, int **InterleavedAddress,
    int girth, int ref, int the_ref, int starting,
    int *path, int pathLen,
    int threshold, int *mult,
    int max_RTZ_size_start_0, int max_RTZ_size_start_1,
    int mingirth, int jump_type, int jump);

#endif /* RTZ_SEARCH_H */
