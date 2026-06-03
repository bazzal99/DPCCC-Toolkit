/**
 * @file utils.h
 * @brief General-purpose utility functions shared across the DPCCC toolkit.
 */

#ifndef UTILS_H
#define UTILS_H

#include <stdio.h>

/* =========================================================================
 * Integer / binary conversion
 * ========================================================================= */

void int_to_binary(int num, int size, int binary_arr[]);
int  binary_to_int(int binary_arr[], int size);
int  circular_modulus(int value, int range);


/* =========================================================================
 * Array utilities
 * ========================================================================= */

int  search(int target, int arr[], int n);
int  Number_Of_Ones(int a[], int size);
void Array_XOR(int a[], int b[], int c[], int size);
void Array_AND(int a[], int b[], int c[], int size);
void copy_array(int source[], int destination[], int size);


/* =========================================================================
 * Sorting and shuffling
 * ========================================================================= */

void bubbleSort(int arr[], int size);
void swap(int *a, int *b);
void shuffleArray(int arr[], int size);
void generateRandomBinaryArray(int *array, int length);
void generateRandom_num_ones(int *array, int length, int num_ones);


/* =========================================================================
 * Noise and signal processing
 * ========================================================================= */

/** Add AWGN noise with variance `sigmasquare` to `binary_signal`. */
void addAWGN(double *binary_signal, double *noisy_signal,
             int length, double sigmasquare);

/** Generate Box-Muller Gaussian noise samples. */
void generateNoise(double *noise, int length);

double MIN(double a, double b);
int    sign(double value);


/* =========================================================================
 * File utilities
 * ========================================================================= */

/** Remove duplicate lines from `filename` in-place. Returns duplicate count. */
int  removeDuplicateLines(const char *filename);

/** Append a sorted copy of `path[0..pathLen-1]` as one line to `filename`. */
void print_in_file_sorted(const char *filename, int path[], int pathLen);

#endif /* UTILS_H */
