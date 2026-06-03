# DPCCC-Toolkit

A complete toolkit for the design and simulation of **Doubly Punctured Concatenated Convolutional Codes (DPCCC)** — a novel turbo code structure introduced in:

> M. Bazzal, J. Nadal, S. Weithoffer, C. Abdel Nour, C. Douillard,
> *"Distance-centric joint interleaver and structural code design for concatenated convolutional codes,"*
> IEEE Open Journal of the Communications Society, 2025.

**Code:** https://github.com/bazzal99/DPCCC-Toolkit

---

## What is a DPCCC?

A DPCCC encodes a K-bit information block through:
1. One RSC encoder (tail-biting) of length K — producing the systematic and first parity stream
2. L independent RSC encoders (tail-biting), each of length K/L — producing L parity sub-streams

The interleaver connecting blocks 1 and 2 is a set of L independent ARP (Almost Regular Permutation) sub-interleavers, one per parity sub-block. This structure allows flexible rate control and efficient iterative decoding.

---

## Repository structure

```
DPCCC-Toolkit/
├── shared/                        ← Shared code used by all three tools
│   ├── config.h                   ← All tunable parameters (edit first)
│   ├── encoder.h / encoder.c      ← DPCCC RSC encoder, circular states, distance
│   ├── interleaver.h / interleaver.c  ← ARP interleaver, repetition mapping
│   ├── decoder.h / decoder.c      ← Log-MAP (BCJR) decoder
│   └── utils.h / utils.c          ← Binary conversion, noise, file I/O
│
├── interleaver_design/
│   ├── decoder_free/              ← Tool 1: decoder-free mHD-based design
│   │   ├── include/
│   │   │   ├── design_config.h    ← Threshold, search depth
│   │   │   ├── rtz_search.h/c    ← DPCCC-aware RTZ graph search
│   │   │   └── designer.h/c      ← Layer-by-layer ARP design
│   │   ├── src/main.c
│   │   └── Makefile
│   │
│   └── double_impulse/            ← Tool 2: N-impulse decoder-based design
│       ├── include/
│       │   ├── design_config.h    ← IMPULSE_N (1/2/3), threshold, SNR
│       │   └── designer.h/c      ← Threaded impulse method
│       ├── src/main.c
│       └── Makefile
│
└── montecarlo/                    ← Tool 3: BER/FER Monte Carlo simulator
    ├── include/sim_config.h       ← Termination mode, RNG seed
    ├── src/main.c
    └── Makefile
```

---

## Configuration

All shared parameters live in **`shared/config.h`**:

```c
#define SIZE      1024   /* Information block length K                  */
#define CUTS      2      /* Number of parity sub-blocks L               */
#define ARP_Q     4      /* ARP sub-layer period Q                      */
#define DELAYS    3      /* RSC memory elements (2^DELAYS trellis states)*/
#define G1_INIT   {1,1,0,1}  /* RSC feedforward polynomial G1          */
#define G2_INIT   {1,0,1,1}  /* RSC feedback polynomial G2             */
#define PUNCT_MASK  1         /* Puncturing period                      */
/* ... SNR range, max iterations, etc. */
```

Tool-specific parameters are in each tool's `include/design_config.h` or `include/sim_config.h`.

---

## Tool 1: Decoder-free interleaver design

Uses a decoder-free RTZ (Return-To-Zero) graph-girth search to find ARP sub-interleavers that maximise the minimum Hamming distance (mHD) without running any decoder.

```bash
cd interleaver_design/decoder_free
# Edit include/design_config.h: set DESIGN_THRESHOLD, SHIFTS_PER_LAYER, etc.
make
./dpccc_decoder_free
```

**Output:** `result.txt` — valid interleavers in C-initialiser format:
```c
int total_S[2][4] = { {1, 3, 7, 2}, {0, 5, 4, 1} };
int total_P[2]    = {11, 7};
distance = 46
```

---

## Tool 2: N-impulse interleaver design

Uses the DPCCC decoder directly as the evaluation oracle. Set `IMPULSE_N` to 1, 2, or 3 for single, double, or triple impulse injection.

```bash
cd interleaver_design/double_impulse
# Edit include/design_config.h: set IMPULSE_N, DESIGN_THRESHOLD, IMPULSE_SNR
make
./dpccc_impulse
```

The search is parallelised across all available CPU cores using pthreads.

---

## Tool 3: Monte Carlo BER/FER simulator

Simulates the BER and FER of a DPCCC code over AWGN using the log-MAP decoder.

```bash
cd montecarlo
# Edit include/sim_config.h: set RNG_SEED, TARGET_FRAME_ERRORS
# Edit shared/interleaver.c: update ARP_SEMI_Initialization() with your S, P
make
./dpccc_sim
```

**Output:** `result.txt` — BER and FER arrays per SNR point:
```
SNR = 4.00
BER = [1.2340000e-02, 3.4560000e-03, ...]
FER = [4.5670000e-02, 1.2340000e-02, ...]
```

> **Note:** The simulator is intentionally single-threaded. For large K or deep
> waterfall regions, consider parallelising `Build_SNR()` — the structure is
> already designed for it.

---

## Building all tools at once

```bash
cd interleaver_design/decoder_free && make && cd ../..
cd interleaver_design/double_impulse && make && cd ../..
cd montecarlo && make && cd ..
```

Each tool requires GCC with C99 support. The impulse tool additionally requires pthreads (`-lpthread`).

---

## Example results (from the paper, K=6144, rate 4/5, CUTS=8)

| Tool | Design time | mHD achieved |
|------|------------|--------------|
| Decoder-free (Tool 1) | ~minutes | d = 10–13 |
| Double-impulse (Tool 2) | ~hours | d = 10–13 |
| Monte Carlo (Tool 3) | ~days (K=6144) | — |

---

## License

MIT License. See [LICENSE](LICENSE).

---

## Citation

```bibtex
@article{bazzal2025dpccc,
  author  = {Bazzal, Mohammad and Nadal, Jeremy and Weithoffer, Stefan
             and {Abdel Nour}, Charbel and Douillard, Catherine},
  title   = {Distance-centric joint interleaver and structural code design
             for concatenated convolutional codes},
  journal = {IEEE Open Journal of the Communications Society},
  year    = {2025}
}
```
