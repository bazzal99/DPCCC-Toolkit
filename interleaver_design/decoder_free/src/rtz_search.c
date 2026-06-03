/**
 * @file rtz_search.c
 * @brief DPCCC-aware RTZ graph-girth search kernel.
 *
 * Each position i in [0, SIZE) has a DPCCC label:
 *   layer    = i % CUTS          (which parity sub-block it belongs to)
 *   position = floor(i / CUTS)   (its index within that sub-block)
 *
 * The search explores cycles in the per-cut ARP correlation graph,
 * respecting the partial-encoder constraint: only positions whose
 * layer <= the_layer are active during the design of cut the_layer.
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../../../shared/config.h"
#include "../../../shared/encoder.h"
#include "../../../shared/utils.h"
#include "../include/rtz_search.h"

int Min_Girth_One_Layer_JUMP_1_7_Iteratively_ARP_INTERLEAVER(
    int current_shift_layer, int position, int layer,
    int **OriginalAddress, int **InterleavedAddress,
    int threshold, int *mult,
    int max_RTZ_size_start_0, int max_RTZ_size_start_1,
    int mingirth, int jump_type, int jump)
{
    int path[40];
    int ref     = position * CUTS + layer;    /* global position index */
    int the_ref = ref;
    return Min_Girth_One_Ref_JUMP_1_7_Iteratively_ARP_INTERLEAVER(
        current_shift_layer, layer, layer,
        OriginalAddress, InterleavedAddress,
        0, ref, the_ref, 0,
        path, 0,
        threshold, mult,
        max_RTZ_size_start_0, max_RTZ_size_start_1,
        mingirth, jump_type, jump);
}

int Min_Girth_One_Ref_JUMP_1_7_Iteratively_ARP_INTERLEAVER(
    int current_shift_layer, int the_layer, int my_layer,
    int **OriginalAddress, int **InterleavedAddress,
    int girth, int ref, int the_ref, int starting,
    int *path, int pathLen,
    int threshold, int *mult,
    int max_RTZ_size_start_0, int max_RTZ_size_start_1,
    int mingirth, int jump_type, int jump)
{
    int i, d, x = 1;
    int next_layer, next_pos;
    int small_fs = SIZE / CUTS;

    if (ref == -1 || ref>SIZE || ref<0) return 1;
    /* Symmetry pruning */
    if (starting == 0 && ref < the_ref)               return 1;
    /* Layer constraint: only positions in layers 0..the_layer */
    if (starting == 0 && ref % CUTS > the_layer)       return 1;

    /* Bounds check */
    if (starting == 0) {
        next_layer = ref % CUTS;
        next_pos   = ref / CUTS;
        if (OriginalAddress[next_layer][next_pos] == -1) return 1;
    }
    if (starting == 1 && InterleavedAddress[my_layer][ref] == -1) return 1;

    /* Record path node */
    if (starting == 1)
        path[pathLen] = InterleavedAddress[my_layer][ref] * CUTS + my_layer;
    else
        path[pathLen] = ref;
    pathLen++;

    if (girth > mingirth) return 1;

    /* Cycle detected */
    if (search(path[pathLen-1], path, pathLen-1) != -1 && girth != 0) {
        d = Turbo_distance_circular_With_Puncturing(path, pathLen-1,
                                                    OriginalAddress);
        if (d < threshold) {
            printf("RTZ cycle: d=%d\n", d);
            return 0;
        }
        return 1;
    }

    /* Expand neighbours */
    if (starting == 0) {
        /* Original domain -> Interleaved domain */
        int budget = (jump_type == 0 || jump_type == 1)
                     ? max_RTZ_size_start_0
                     : max_RTZ_size_start_0;
        int step   = (jump_type == 2 || jump_type == 3) ? jump : 1;

        for (i = step; i <= budget * step; i += step) {
            if (x == 0) return 0;
            int cm_minus = circular_modulus(ref - i, SIZE);
            int cm_plus  = circular_modulus(ref + i, SIZE);

            next_layer = cm_minus % CUTS;
            next_pos   = cm_minus / CUTS;
            if (next_layer <= the_layer)
                x &= Min_Girth_One_Ref_JUMP_1_7_Iteratively_ARP_INTERLEAVER(
                         current_shift_layer, the_layer, next_layer,
                         OriginalAddress, InterleavedAddress,
                         girth+1, OriginalAddress[next_layer][next_pos],
                         the_ref, 1, path, pathLen,
                         threshold, mult, max_RTZ_size_start_0,
                         max_RTZ_size_start_1, mingirth, jump_type, jump);

            if (x == 0) return 0;
            next_layer = cm_plus % CUTS;
            next_pos   = cm_plus / CUTS;
            if (next_layer <= the_layer)
                x &= Min_Girth_One_Ref_JUMP_1_7_Iteratively_ARP_INTERLEAVER(
                         current_shift_layer, the_layer, next_layer,
                         OriginalAddress, InterleavedAddress,
                         girth+1, OriginalAddress[next_layer][next_pos],
                         the_ref, 1, path, pathLen,
                         threshold, mult, max_RTZ_size_start_0,
                         max_RTZ_size_start_1, mingirth, jump_type, jump);
        }
    } else {
        /* Interleaved domain -> Original domain */
        int budget = max_RTZ_size_start_1;
        int step   = (jump_type == 1 || jump_type == 3) ? jump : 1;

        for (i = step; i <= budget * step; i += step) {
            if (x == 0) return 0;
            int cm_minus = circular_modulus(ref - i, small_fs);
            int cm_plus  = circular_modulus(ref + i, small_fs);

            x &= Min_Girth_One_Ref_JUMP_1_7_Iteratively_ARP_INTERLEAVER(
                     current_shift_layer, the_layer, my_layer,
                     OriginalAddress, InterleavedAddress,
                     girth+1, cm_minus * CUTS + my_layer,
                     the_ref, 0, path, pathLen,
                     threshold, mult, max_RTZ_size_start_0,
                     max_RTZ_size_start_1, mingirth, jump_type, jump);

            if (x == 0) return 0;
            x &= Min_Girth_One_Ref_JUMP_1_7_Iteratively_ARP_INTERLEAVER(
                     current_shift_layer, the_layer, my_layer,
                     OriginalAddress, InterleavedAddress,
                     girth+1, cm_plus * CUTS + my_layer,
                     the_ref, 0, path, pathLen,
                     threshold, mult, max_RTZ_size_start_0,
                     max_RTZ_size_start_1, mingirth, jump_type, jump);
        }
    }
    return x;
}
