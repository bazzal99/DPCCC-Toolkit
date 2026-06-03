/**
 * @file utils.c
 * @brief General-purpose utility function implementations.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "config.h"
#include "utils.h"

#define HASH_TABLE_SIZE  1009
#define MAX_LINE_LEN     60

/* =========================================================================
 * Integer / binary conversion
 * ========================================================================= */

void int_to_binary(int num, int size, int binary_arr[])
{
    for (int i = size - 1; i >= 0; i--) { binary_arr[i] = num % 2; num /= 2; }
}

int binary_to_int(int binary_arr[], int size)
{
    int dec = 0;
    for (int i = 0; i < size; i++)
        dec += binary_arr[i] * (int)pow(2, size - 1 - i);
    return dec;
}

int circular_modulus(int value, int range)
{
    return (value % range + range) % range;
}

/* =========================================================================
 * Array utilities
 * ========================================================================= */

int search(int target, int arr[], int n)
{
    for (int i = 0; i < n; i++) if (arr[i] == target) return i;
    return -1;
}

int Number_Of_Ones(int a[], int size)
{
    int s = 0;
    for (int i = 0; i < size; i++) s += a[i];
    return s;
}

void Array_XOR(int a[], int b[], int c[], int size)
{
    for (int i = 0; i < size; i++) c[i] = a[i] ^ b[i];
}

void Array_AND(int a[], int b[], int c[], int size)
{
    for (int i = 0; i < size; i++) c[i] = a[i] & b[i];
}

void copy_array(int source[], int destination[], int size)
{
    for (int i = 0; i < size; i++) destination[i] = source[i];
}

/* =========================================================================
 * Sorting / shuffling
 * ========================================================================= */

void bubbleSort(int arr[], int size)
{
    for (int i = 0; i < size - 1; i++)
        for (int j = 0; j < size - i - 1; j++)
            if (arr[j] > arr[j + 1]) { int t = arr[j]; arr[j] = arr[j+1]; arr[j+1] = t; }
}

void swap(int *a, int *b) { int t = *a; *a = *b; *b = t; }

void shuffleArray(int arr[], int size)
{
    for (int i = size - 1; i > 0; i--)
        swap(&arr[i], &arr[rand() % (i + 1)]);
}

void generateRandomBinaryArray(int *array, int length)
{
    for (int i = 0; i < length; i++) array[i] = rand() % 2;
}

void generateRandom_num_ones(int *array, int length, int num_ones)
{
    for (int i = 0; i < length; i++) array[i] = 0;
    while (num_ones > 0) {
        int pos = rand() % length;
        if (array[pos] == 0) { array[pos] = 1; num_ones--; }
    }
}

/* =========================================================================
 * Noise and signal
 * ========================================================================= */

#define TWO_PI  6.28318530717959

void generateNoise(double *noise, int length)
{
    for (int i = 0; i < length; i++) {
        double u1 = ((double)rand() + 1.) / ((double)RAND_MAX + 1.);
        double u2 = ((double)rand() + 1.) / ((double)RAND_MAX + 1.);
        noise[i]  = sqrt(-2.0 * log(u2)) * cos(TWO_PI * u1);
    }
}

void addAWGN(double *binary_signal, double *noisy_signal,
             int length, double sigmasquare)
{
    double noise[length];
    generateNoise(noise, length);
    for (int i = 0; i < length; i++)
        noisy_signal[i] = binary_signal[i] + sqrt(sigmasquare) * noise[i];
}

double MIN(double a, double b) { return (a < b) ? a : b; }
int    sign(double value) { return (value > 0) ? 1 : -1; }

/* =========================================================================
 * File utilities
 * ========================================================================= */

static unsigned long hash_line(const char *str)
{
    unsigned long h = 5381;
    int c;
    while ((c = *str++)) h = ((h << 5) + h) + c;
    return h % HASH_TABLE_SIZE;
}

struct Node { char line[MAX_LINE_LEN]; struct Node *next; };

static void ht_insert(struct Node *table[], const char *line)
{
    unsigned long idx = hash_line(line);
    struct Node *node = malloc(sizeof(struct Node));
    if (!node) { perror("malloc"); exit(1); }
    strncpy(node->line, line, MAX_LINE_LEN - 1);
    node->line[MAX_LINE_LEN - 1] = '\0';
    node->next = table[idx];
    table[idx] = node;
}

static int ht_find(struct Node *table[], const char *line)
{
    unsigned long idx = hash_line(line);
    struct Node *cur  = table[idx];
    while (cur) { if (strcmp(cur->line, line) == 0) return 1; cur = cur->next; }
    return 0;
}

int removeDuplicateLines(const char *filename)
{
    FILE *f = fopen(filename, "r");
    if (!f) { perror("fopen"); return -1; }

    struct Node *table[HASH_TABLE_SIZE];
    for (int i = 0; i < HASH_TABLE_SIZE; i++) table[i] = NULL;

    FILE *tmp = fopen("_dedup_tmp.txt", "w");
    if (!tmp) { perror("fopen tmp"); fclose(f); return -1; }

    char line[MAX_LINE_LEN];
    int dups = 0;
    while (fgets(line, sizeof(line), f)) {
        if (ht_find(table, line)) { dups++; }
        else { ht_insert(table, line); fputs(line, tmp); }
    }
    fclose(f); fclose(tmp);
    if (remove(filename) != 0)               { perror("remove"); return -1; }
    if (rename("_dedup_tmp.txt", filename) != 0) perror("rename");
    return dups;
}

void print_in_file_sorted(const char *filename, int path[], int pathLen)
{
    int newpath[pathLen];
    copy_array(path, newpath, pathLen);
    bubbleSort(newpath, pathLen);
    FILE *f = fopen(filename, "a");
    for (int i = 0; i < pathLen; i++)
        fprintf(f, (i == pathLen - 1) ? "%d" : "%d,", newpath[i]);
    fprintf(f, "\n");
    fclose(f);
}
