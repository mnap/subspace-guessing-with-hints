# Code Supplement to Paper
Paper: **Subspace Guessing and Rank-Metric Solvers with Hints** by *Anmoal Porwal, Harrison Banda, Jan Brinkmann, Anna Baumeister, Juliane Krämer, Antonia Wachter-Zeh* (Full version: https://eprint.iacr.org/2026/132).

This repository can be used to
- estimate the complexity of attacking the MinRank and RSD problems (for Mirath and RYDE parameter sets) under various fractions of hints on the secret decomposition matrices `S` and `C'`. It uses the modified kernel search and (row-space) GRS algorithms described in the paper.
- compute, for each fixed `S`-hint fraction, the minimum `C'`-hint fraction ("threshold") at which these attacks become polynomial-time.

## Installation
This project requires Python 3.13 or newer. No external runtime dependencies are required.
You should therefore be able to run:
```bash
python main.py --help
# alternative with uv:
# uv run python main.py --help
```

## Usage
If you are using `uv`, replace `python` with `uv run python` for all commands below.

Show usage and options:
```bash
python main.py --help
```

Options:
```
  -h, --help            show this help message and exit
  --family {all,mirath,ryde}
                        Restrict the default analysis to one parameter family.
  --name NAME           Analyze only the named parameter set(s). Can be passed multiple times.
  --mode {all,entry,bit}
                        Choose whether to run entry-level, bit-level, or both models. Bit-level analysis is only run for q = 2^nu with nu > 1; other parameter sets are
                        skipped.
  --trials TRIALS       Monte Carlo trials used to estimate t_c.
  --s-grid-steps S_GRID_STEPS
                        Number of intervals for the S-hint fraction grid between 0 and 1; use 0 to run only fraction_S = 0.
  --c-grid-steps C_GRID_STEPS
                        Number of intervals for the C-hint fraction grid between 0 and 1 in the average-complexity tables only; use 0 to run only fraction_C = 0. This option
                        does not affect threshold tables.
  --seed SEED           Base RNG seed for repeatable Monte Carlo estimates.
  --report {thresholds,average,both}
                        Choose which summary to print: threshold table, average-complexity table, or both.
  --averaging {harmonic,mean}
                        How to average sampled complexities: harmonic means 1 / average(1/E), mean means average(E).
  --distribution {random,balanced}
                        How to sample C' hints: random placement with fixed total weight, or as balanced across columns as possible.
```

Run for all parameter sets and reproduce the values in the full version of the paper:
```bash
python main.py # might take a few minutes
```

The above default command `python main.py` is equivalent to the following:
```bash
# Run for all parameter sets, compute both Fq-entry and bit-level formulas ("mode" option),
# use 2000 Monte Carlo trials, only fraction_S = 0 and a 10-step grid for C' in the
# average-complexity tables, RNG seed 0, print both summaries, average with 1 / average(1/E),
# and use random C' hint placement.
python main.py \
  --family all \
  --mode all \
  --trials 2000 \
  --s-grid-steps 0 \
  --c-grid-steps 10 \
  --seed 0 \
  --report both \
  --averaging harmonic \
  --distribution random
```

For just one parameter set, you can run, e.g., `python main.py --name Mirath-1a` (see `MIRATH_PARAMS` and `RYDE_PARAMS` in `main.py` for all possible names).

## How the Estimates Are Computed

For each parameter set, the code evaluates the exponential factor of the complexity formulas from the paper:
- the kernel-search algorithm with hints for Mirath / MinRank parameter sets
- the row-space GRS algorithm with hints for RYDE / RSD parameter sets

The polynomial prefactor is ignored. In other words, the program reports only the exponential part of the attack cost.

### Hint Models
Given a fixed total number of hints `h_c` in `C'`, `--distribution random` (default) places the hints uniformly at random in the matrix. The option `--distribution balanced` distributes them as evenly across columns as possible. The balanced distribution is the worst case for fixed `h_c`, i.e., it yields the highest complexity for the kernel-search-with-hints and GRS-with-hints attacks described in the paper.

Entry-level analysis treats each known field element as one hint. Bit-level analysis treats each known bit in the binary representation of a field element as one hint. It is only run when the field size is of the form `q = 2^nu` with `nu > 1`; parameter sets with `q = 2` are therefore skipped in bit mode.

### Averaging Rules

Let `E_i` be the exponential complexity computed in trial `i`. With `--averaging mean`, the reported value is

```text
E_final = (sum_i E_i) / TRIALS.
```

With `--averaging harmonic` (default), it is

```text
E_final = 1 / ((sum_i 1 / E_i) / TRIALS).
```

The harmonic rule is equivalent to averaging the corresponding success probabilities over all trials and it aligns better with the usual definition of workfactor in cryptography since that is defined in terms of average success probability. Note when `--distribution balanced` is used, all `E_i` are the same, so both averaging rules give the same result.

### Average-Complexity Table

Consider the command

```bash
python main.py \
  --name Mirath-1a \
  --mode entry \
  --trials 2000 \
  --s-grid-steps 2 \
  --c-grid-steps 10 \
  --report average \
  --averaging harmonic \
  --distribution random
```

For `Mirath-1a`, the parameters are `m = 16`, `n = 16`, and `r = 4`. Hence:
- `S` has dimensions `m x r = 16 x 4`, so it contains `64` entry-level coordinates
- `C'` has dimensions `r x (n - r) = 4 x 12`, so it contains `48` entry-level coordinates

The options `--s-grid-steps 2` and `--c-grid-steps 10` define the hint-fraction grid for the average-complexity table:
- the `S` fractions are `0`, `0.5`, and `1.0`
- the `C'` fractions are `0`, `0.1`, `0.2`, ..., `1.0`

So the code evaluates the average complexity on all pairs
`(fraction_S, fraction_C)` in that grid.

For example, for the pair `(0, 0.2)`:
- `fraction_S = 0` means `h_s = 0` entry hints in `S`
- `fraction_C = 0.2` means `h_c = round(0.2 * 48) = 10` entry hints in `C'`

Since `--distribution random` and `--trials 2000` are used, the code performs `2000` trials. In each trial, it places these `10` hints uniformly at random in `C'`, computes the corresponding exponential complexity `E_i` using the entry-level kernel-search-with-hints formula, and then averages the resulting values according to the selected averaging rule.

With `--averaging harmonic`, the reported complexity is

```text
E_final = 1 / ((sum_i 1 / E_i) / 2000).
```

This value is computed for every grid point. The tables report `log2(E_final)` and are grouped as follows:
- `ENTRY | MINRANK`
- `ENTRY | RSD`
- `BIT | MINRANK`
- `BIT | RSD`

Within each table, the column header is printed once and all relevant parameter sets are listed underneath it.

### Threshold Table

The threshold table is derived from the same underlying complexity computation, but summarized differently.

For each `S`-hint fraction in the grid, the code finds the minimum `C'`-hint fraction for which `E_final <= 2`. This is done by binary search over the `C'`-hint count.

Only `--s-grid-steps` affects the threshold table. The option `--c-grid-steps` is not used there.

With `--s-grid-steps 2`, the threshold table therefore contains one row for each of the three `S` fractions: `0`, `0.5`, `1.0`.

For each of these rows, the code reports the smallest `C'` hint fraction that makes the attack polynomial-time under the `E_final <= 2` convention, using the selected distribution and averaging mode.

The threshold output is also grouped into the same four tables:
- `ENTRY | MINRANK`
- `ENTRY | RSD`
- `BIT | MINRANK`
- `BIT | RSD`

Each row corresponds to one parameter set and one `S` fraction. It contains:
- `min_fraction_C`, the minimum `C'` hint fraction
- `log2_avg_complexity_at_threshold`, the base-2 logarithm of `E_final` at that fraction

The threshold tables use `E_final <= 2` as the convention for calling an attack polynomial-time. Under the harmonic averaging rule, `E_final <= 2` is equivalent to requiring that the average success probability is at least `0.5`. Because each `E_i` is at least `1`, the alternative condition `E_final <= 1` would require every sampled complexity to equal `1` and would therefore almost never be satisfied.

## Acknowledgments
Parts of this code were generated or revised with the help of OpenAI Codex v0.125.0 using model gpt-5.4.
