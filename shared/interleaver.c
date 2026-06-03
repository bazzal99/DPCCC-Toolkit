/**
 * @file interleaver.c
 * @brief DPCCC ARP interleaver and repetition mapping implementations.
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "config.h"
#include "interleaver.h"
#include "encoder.h"
#include "utils.h"

/* =========================================================================
 * Non-uniform repetition
 * ========================================================================= */

void Non_uniform_repeate(int input[], int new_input[], int Places[][MAX_SEGMENTS])
{
    for (int i = 0; i < SIZE; i++)
        for (int j = 0; j < MAX_SEGMENTS; j++) {
            if (Places[i][j] == -1) break;
            new_input[Places[i][j]] = input[i];
        }
}

void Non_uniform_repeate_double(double input[], double new_input[], int Places[][MAX_SEGMENTS])
{
    for (int i = 0; i < SIZE; i++)
        for (int j = 0; j < MAX_SEGMENTS; j++) {
            if (Places[i][j] == -1) break;
            new_input[Places[i][j]] = input[i];
        }
}

/* =========================================================================
 * Full-frame ARP permutation
 * ========================================================================= */

void ARP_Interleave(int data[], int output[], int OriginalAddress[])
{
    int mydata[FRAME_SIZE];
    copy_array(data, mydata, FRAME_SIZE);
    for (int i = 0; i < FRAME_SIZE; i++)
        output[OriginalAddress[i]] = mydata[i];
}

void ARP_DE_Interleave(int data[], int output[], int InterleavedAddress[])
{
    int mydata[FRAME_SIZE];
    copy_array(data, mydata, FRAME_SIZE);
    for (int i = 0; i < FRAME_SIZE; i++)
        output[InterleavedAddress[i]] = mydata[i];
}

void ARP_Interleave_double(double data[], double output[], int OriginalAddress[])
{
    double mydata[FRAME_SIZE];
    for (int i = 0; i < FRAME_SIZE; i++) mydata[i] = data[i];
    for (int i = 0; i < FRAME_SIZE; i++) output[OriginalAddress[i]] = mydata[i];
}

void ARP_DE_Interleave_double(double data[], double output[], int InterleavedAddress[])
{
    double mydata[FRAME_SIZE];
    for (int i = 0; i < FRAME_SIZE; i++) mydata[i] = data[i];
    for (int i = 0; i < FRAME_SIZE; i++) output[InterleavedAddress[i]] = mydata[i];
}

/* =========================================================================
 * Per-segment ARP initialisation
 * ========================================================================= */

void ARPInterleaverInitialization(int Size, int Period, int PeriodARP,
                                  int Shift[], int InterleavedAddress[], int OriginalAddress[])
{
    for (int i = 0; i < Size; i++) {
        OriginalAddress[i] = (int)((i * Period + Shift[i % PeriodARP]) % Size);
        InterleavedAddress[OriginalAddress[i]] = i;
    }
}

void ARPInterleaverInitialization_Per_Layer(int Size, int layer, int Period,
                                            int Q, int Shift[],
                                            int InterleavedAddress[], int OriginalAddress[])
{
    for (int i = layer; i < Size; i += Q) {
        OriginalAddress[i] = (int)((i * Period + Shift[i % Q]) % Size);
        InterleavedAddress[OriginalAddress[i]] = i;
    }
}

void ARPInterleaverInitialization_Till_Layer(int Size, int layer, int Period,
                                             int Q, int Shift[],
                                             int InterleavedAddress[], int OriginalAddress[])
{
    for (int j = 0; j <= layer; j++)
        for (int i = j; i < Size; i += Q) {
            OriginalAddress[i] = (int)((i * Period + Shift[i % Q]) % Size);
            InterleavedAddress[OriginalAddress[i]] = i;
        }
}

void Arp_DE_Interleave_Layer(int Size, int layer, int Period, int Q,
                             int Shift[], int InterleavedAddress[], int OriginalAddress[])
{
    for (int i = layer; i < Size; i += Q) {
        OriginalAddress[i] = (int)((i * Period + Shift[i % Q]) % Size);
        InterleavedAddress[OriginalAddress[i]] = -1;
        OriginalAddress[i] = -1;
    }
}

/* =========================================================================
 * Period selection
 * ========================================================================= */

int gcd(int a, int b) { return (b == 0) ? a : gcd(b, a % b); }
int isRelativePrime(int a, int b) { return gcd(a, b) == 1; }

void saveRelativePrimes(int number, int *result, int *count)
{
    for (int i = 1; i <= number; i++)
        if (isRelativePrime(number, i)) result[(*count)++] = i;
}

int Min_Spread(int p)
{
    int upper = (int)sqrt(2.0 * SIZE);
    int interleaved[SIZE];
    int smin = SIZE;
    for (int i = 0; i < SIZE; i++) interleaved[(p * i) % SIZE] = i;
    for (int i = 0; i < SIZE; i++)
        for (int j = 1; j < upper; j++)
            if (i != j) {
                int x = abs(j) + abs(interleaved[i] - interleaved[(i + j) % SIZE]);
                if (x < smin) smin = x;
            }
    return smin;
}

void Choose_Period(int periods[], int *size, int framesize)
{
    int result[framesize], count = 0;
    saveRelativePrimes(framesize, result, &count);
    for (int p = 0; p < count; p++)
        if (Min_Spread(result[p]) >= ((int)sqrt(2.0 * framesize) - 2))
            periods[(*size)++] = result[p];
}

int inter(int i, int p, int s[]) { return (p * i + s[i % ARP_Q]) % SIZE; }

int isArrayUnique(int arr[], int Q)
{
    for (int i = 0; i < Q - 1; i++)
        for (int j = i + 1; j < Q; j++)
            if (arr[i] == arr[j]) return 0;
    return 1;
}

void generateRandomArpInterleaver(int S[], int *P, int Q, int frame_size)
{
    int primes[frame_size], count = 0;
    saveRelativePrimes(frame_size, primes, &count);
    *P = primes[rand() % count];
    do {
        for (int i = 0; i < Q; i++) S[i] = rand() % frame_size;
    } while (!isArrayUnique(S, Q));
}

/* =========================================================================
 * Full-frame interleaver (hardcoded paper interleavers)
 * ========================================================================= */

void ARP_SEMI_Initialization(int OriginalAddress[], int InterleavedAddress[],
                              int number_of_repetition[MAX_SEGMENTS])
{
    int i, j, beginning;

    /* --- Interleaver parameters from the paper (K=6144, rate 4/5) ---
     * Replace these arrays with your designed S[], P values.
     * total_S[l] and total_P[l] are the ARP parameters for cut l. */
    int Q0 = 4;
    int total_S[8][4] = {
        {1, 417,   9, 313},
        {1, 582, 446, 371},
        {3, 161, 151, 685},
        {2, 426, 406, 766},
        {0, 332, 376, 224},
        {3, 407, 764, 382},
        {0, 247, 161,  81},
        {2, 311, 306, 601}
    };
    int total_P[8] = {11, 27, 5, 21, 15, 19, 3, 7};

    /* Systematic part: identity */
    for (i = 0; i < SIZE; i++) OriginalAddress[i] = i;

    /* Parity sub-blocks */
    beginning = SIZE;
    for (i = 0; i < CUTS; i++) {
        int seg_size = SIZE / CUTS;
        for (j = 0; j < seg_size; j++) {
            OriginalAddress[beginning + j] =
                (int)((j * total_P[i] + total_S[i][j % Q0]) % seg_size)
                + beginning;
        }
        beginning += seg_size;
    }

    for (i = 0; i < FRAME_SIZE; i++)
        InterleavedAddress[OriginalAddress[i]] = i;
}

void S_Random(int s, int output[], int size)
{
    int i, j, x, flag;
    for (i = 0; i < size; i++) {
        do {
            x = rand() % size;
            flag = 0;
            for (j = 0; j < s && (i - j - 1) >= 0; j++)
                if (abs(x - output[i - j - 1]) < s) { flag = 1; break; }
            for (j = 0; j < i; j++)
                if (x == output[j]) { flag = 1; break; }
        } while (flag);
        output[i] = x;
    }
}

void S_Random_Interleaver_Initialization(int s, int OriginalAddress[],
                                          int InterleavedAddress[],
                                          int number_of_repetition[MAX_SEGMENTS])
{
    int temp[FRAME_SIZE];
    int beginning = 0, mysize;

    for (int i = 0; i < SIZE; i++) OriginalAddress[i] = i;
    mysize = FRAME_SIZE - SIZE;
    beginning = SIZE;
    S_Random(s, temp, mysize);
    for (int i = 0; i < mysize; i++)
        OriginalAddress[beginning + i] = temp[i] + beginning;
    for (int i = 0; i < FRAME_SIZE; i++)
        InterleavedAddress[OriginalAddress[i]] = i;
}

/* =========================================================================
 * Extrinsic LLR computation (turbo decoder message passing)
 * ========================================================================= */

double Calculate_Extrinsic(int Places[][MAX_SEGMENTS], int row, int column,
                            double L_D[], double L_X[], double La_D[])
{
    double extrinsic = 0.0;
    for (int i = 0; i < MAX_SEGMENTS; i++) {
        int e = Places[row][i];
        if (e == -1) break;
        if (e != Places[row][column])
            extrinsic += L_D[e] - L_X[Places[row][0]] - La_D[e];
    }
    return extrinsic;
}
