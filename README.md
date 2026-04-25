# Code Supplement to Paper
Paper: **Subspace Guessing and Rank-Metric Solvers with Hints** by *Anmoal Porwal, Harrison Banda, Jan Brinkmann, Anna Baumeister, Juliane Krämer, Antonia Wachter-Zeh* (Full version: https://eprint.iacr.org/2026/132).

This repository can be used to
- estimate the complexity of attacking the MinRank and RSD problems (for Mirath and RYDE parameter sets) under various fractions of hints on the secret decomposition matrices `S` and `C'`. It uses the modified kernel search and (row-space) GRS algorithms described in the paper.
- compute the minimum hint fraction at which these attacks become polynomial-time.

## Installation
This project requires Python 3.13 or newer. No external runtime dependencies are required.
You should therefore be able to run:
```bash
python main.py --help
# alternative with uv:
# uv run python main.py --help
```

## Usage
If you are using `uv`, replace `python` with `uv run python` for all commands.

Show all options:
```bash
python main.py --help
```

Run for all parameter sets (this prints both the threshold table and the average-complexity table):
```bash
python main.py # might take a few minutes
```

Run only for one parameter set:
```bash
python main.py --name Mirath-1a
python main.py --name RYDE-1
```

Set number of Monte Carlo trials and the hint-fraction grids:

```bash
python main.py --trials 5000 --s-grid-steps 20 --c-grid-steps 20
```

Choose which summary to compute and print:

```bash
python main.py --report thresholds
python main.py --report average
python main.py --report both
```

Choose the `C'` hint-placement model:

```bash
python main.py --distribution random
python main.py --distribution balanced
```
Given a fixed total number of hints `h_c` in the matrix `C'`, the option `distribution=random` places these hints uniformly at random in `C'`. The option `distribution=balanced` distributes them as evenly across columns as possible. Note the balanced distribution is the worst case: for fixed `h_c`, it yields the highest complexity for the attacks described in the paper (kernel-search-with-hints and GRS-with-hints).

```bash
python main.py --averaging mean
python main.py --averaging harmonic
```
With `mean`, the average complexity is computed as `(sum_i E_i) / TRIALS`, where `E_i` is the complexity in the `i`th trial. With `harmonic`, it is computed as `1 / ((sum_i 1 / E_i) / TRIALS)`.
The `harmonic` strategy is equivalent to averaging the corresponding success probabilities over all trials. When `distribution=balanced`, all `E_i` are the same, so both averaging strategies give the same result.

The default command `python main.py` is equivalent to the following:
```bash
# Run for all parameter sets, compute both Fq-entry and bit-level formulas ("mode" option),
# use 2000 Monte Carlo trials, a 2-step grid for S and a 10-step grid for C', RNG seed 0,
# print both summaries, average with 1 / average(1/E), and use random C' hint placement.
python main.py \
  --family all \
  --mode all \
  --trials 2000 \
  --s-grid-steps 2 \
  --c-grid-steps 10 \
  --seed 0 \
  --report both \
  --averaging harmonic \
  --distribution random
```

To reproduce the values in the full-version paper, run the following:
```bash
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

## How The Estimates Are Computed

For each parameter set, the code evaluates the exponential factor of the complexity formulas from the paper:
- the kernel-search algorithm with hints for Mirath / MinRank parameter sets
- the row-space GRS algorithm with hints for RYDE / RSD parameter sets

The polynomial prefactor is ignored. In other words, the program reports only the exponential part of the attack cost.

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

The options `--s-grid-steps 2` and `--c-grid-steps 10` define the hint-fraction grid:
- the `S` fractions are `0`, `0.5`, and `1.0`
- the `C'` fractions are `0`, `0.1`, `0.2`, ..., `1.0`

So the code evaluates the average complexity on all pairs
`(fraction_S, fraction_C)` in that grid.

For example, for the pair `(0, 0.2)`:
- `fraction_S = 0` means `h_s = 0` entry hints in `S`
- `fraction_C = 0.2` means `h_c = round(0.2 * 48) = 10` entry hints in `C'`

Since `--distribution random` and `--trials 2000` are used, the code performs `2000` trials. In each trial, it places these `10` hints uniformly at random in `C'`, computes the corresponding exponential complexity `E_i` using the entry-level kernel-search-with-hints formula, and then averages the resulting values according to the selected averaging rule.

With `--averaging harmonic`, the reported value is

```text
E_final = 1 / ((sum_i 1 / E_i) / 2000).
```

This value is computed for every grid point and printed as the average-complexity table.

### Threshold Table

The threshold table is derived from the same underlying complexity computation, but summarized differently.

For each `S`-hint fraction in the grid, the code finds the minimum `C'`-hint fraction for which the sampled average complexity is at most `2`. This is done by binary search over the `C'`-hint count.

With `--s-grid-steps 2`, the threshold table therefore contains one row for each of the three `S` fractions:
- `0`
- `0.5`
- `1.0`

For each of these rows, the code reports the smallest `C'` hint fraction that makes the attack polynomial-time under the convention `average(E) <= 2`, using the selected distribution and averaging mode.

Note: the threshold table uses the condition `average(E) <= 2` under the selected averaging mode. Under the `harmonic` averaging rule, this is equivalent to requiring that the average success probability is at least `0.5`. (Further, note that `E` can never be below `1`, so the alternative condition `average(E) <= 1` would only succeed when all sampled complexities `E_i` are equal to `1`, and hence the condition would basically never be satisfied.)

## Footnote
This code was written with the help of OpenAI Codex v0.125.0 (model gpt-5.4).