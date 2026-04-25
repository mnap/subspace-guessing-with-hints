import argparse
from types import SimpleNamespace

from simulation import AnalysisConfig, analyze_parameter_sets


DEFAULT_TRIALS = 2000
DEFAULT_S_GRID_STEPS = 2
DEFAULT_C_GRID_STEPS = 10
FRACTION_PRECISION = 2
LOG2_COMPLEXITY_PRECISION = 0

MIRATH_PARAMS = [
    # Mirath Specification Document, Table 4
    # For MinRank analysis we interpret `k` as the theorem's `K`.
    SimpleNamespace(name="Mirath-1a", family="mirath", q=16, m=16, n=16, k=143, r=4, bitsec=158),
    SimpleNamespace(name="Mirath-3a", family="mirath", q=16, m=19, n=19, k=195, r=5, bitsec=225),
    SimpleNamespace(name="Mirath-5a", family="mirath", q=16, m=22, n=22, k=255, r=6, bitsec=301),
    SimpleNamespace(name="Mirath-1b", family="mirath", q=2, m=42, n=42, k=1443, r=4, bitsec=158),
    SimpleNamespace(name="Mirath-3b", family="mirath", q=2, m=50, n=50, k=2024, r=5, bitsec=224),
    SimpleNamespace(name="Mirath-5b", family="mirath", q=2, m=56, n=56, k=2499, r=6, bitsec=289),
]

RYDE_PARAMS = [
    # RYDE Specification Document, Table 4
    SimpleNamespace(name="RYDE-1", family="ryde", q=2, m=53, n=53, k=45, r=4),
    SimpleNamespace(name="RYDE-3", family="ryde", q=2, m=61, n=61, k=51, r=5),
    SimpleNamespace(name="RYDE-5", family="ryde", q=2, m=61, n=57, k=44, r=7),
]

ALL_PARAMS = MIRATH_PARAMS + RYDE_PARAMS


def parse_args():
    """Parse CLI options for selecting parameter sets and analysis settings."""
    parser = argparse.ArgumentParser(
        description=(
            "Estimate the fraction of hints needed for polynomial-time kernel search "
            "(MinRank) and GRS row-space (RSD) under configurable hint placement."
        )
    )
    parser.add_argument(
        "--family",
        choices=("all", "mirath", "ryde"),
        default="all",
        help="Restrict the default analysis to one parameter family.",
    )
    parser.add_argument(
        "--name",
        action="append",
        default=[],
        help="Analyze only the named parameter set(s). Can be passed multiple times.",
    )
    parser.add_argument(
        "--mode",
        choices=("all", "entry", "bit"),
        default="all",
        help="Choose whether to run entry-level, bit-level, or both models. Bit-level analysis is only run for q = 2^nu with nu > 1; other parameter sets are skipped.",
    )
    parser.add_argument(
        "--trials",
        type=int,
        default=DEFAULT_TRIALS,
        help="Monte Carlo trials used to estimate t_c.",
    )
    parser.add_argument(
        "--s-grid-steps",
        type=int,
        default=DEFAULT_S_GRID_STEPS,
        help="Number of intervals for the S-hint fraction grid between 0 and 1; use 0 to run only fraction_S = 0.",
    )
    parser.add_argument(
        "--c-grid-steps",
        type=int,
        default=DEFAULT_C_GRID_STEPS,
        help="Number of intervals for the C-hint fraction grid between 0 and 1 in the average-complexity tables only; use 0 to run only fraction_C = 0. This option does not affect threshold tables.",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=0,
        help="Base RNG seed for repeatable Monte Carlo estimates.",
    )
    parser.add_argument(
        "--report",
        choices=("thresholds", "average", "both"),
        default="both",
        help="Choose which summary to print: threshold table, average-complexity table, or both.",
    )
    parser.add_argument(
        "--averaging",
        choices=("harmonic", "mean"),
        default="harmonic",
        help="How to average sampled complexities: harmonic means 1 / average(1/E), mean means average(E).",
    )
    parser.add_argument(
        "--distribution",
        choices=("random", "balanced"),
        default="random",
        help="How to sample C' hints: random placement with fixed total weight, or as balanced across columns as possible.",
    )
    return parser.parse_args()


def select_params(args):
    """Filter the parameter sets according to CLI family and name selectors."""
    params = ALL_PARAMS
    if args.family != "all":
        params = [param for param in params if param.family == args.family]
    if args.name:
        wanted = set(args.name)
        params = [param for param in params if param.name in wanted]
    return params


def format_fraction(value):
    """Pretty-print a threshold fraction, preserving infeasible points."""
    if value is None:
        return "infeasible"
    return f"{value:.{FRACTION_PRECISION}f}"


def grouped_results(analysis_results):
    """Group analysis results by mode and attack in display order."""
    display_order = (("entry", "minrank"), ("entry", "rsd"), ("bit", "minrank"), ("bit", "rsd"))
    grouped = {key: [] for key in display_order}
    for result in analysis_results:
        key = (result["mode"], result["attack"])
        if key in grouped:
            grouped[key].append(result)
    return [(mode, attack, grouped[(mode, attack)]) for mode, attack in display_order if grouped[(mode, attack)]]


def print_aligned_table(rows, left_align_columns):
    """Print a table with per-column width alignment."""
    if not rows:
        return

    widths = [max(len(row[index]) for row in rows) for index in range(len(rows[0]))]
    for row in rows:
        formatted = []
        for index, cell in enumerate(row):
            if index in left_align_columns:
                formatted.append(cell.ljust(widths[index]))
            else:
                formatted.append(cell.rjust(widths[index]))
        print("  ".join(formatted))


def print_threshold_summary(analysis_results, args):
    """Print grouped minimum C'-hint threshold tables."""
    print(f"Monte Carlo trials = {args.trials}")
    print(f"S-grid steps       = {args.s_grid_steps}")
    print(f"Base seed          = {args.seed}")
    print(f"Averaging          = {args.averaging}")
    print(f"Distribution       = {args.distribution}")

    for mode, attack, results in grouped_results(analysis_results):
        fraction_s_values = [point["fraction_s"] for point in results[0]["points"]]
        rows = [["parameter_set", "fraction_S", "min_fraction_C", "log2_avg_complexity_at_threshold"]]

        print()
        print(f"THRESHOLDS | {mode.upper()} | {attack.upper()}")

        for fraction_s in fraction_s_values:
            for result in results:
                params = result["params"]
                point = next(point for point in result["points"] if point["fraction_s"] == fraction_s)
                row = [params.name, f"{point['fraction_s']:.{FRACTION_PRECISION}f}", format_fraction(point["threshold"])]
                value = point["log2_average_complexity_at_threshold"]
                if value is None:
                    row.append("n/a")
                elif abs(value) < 5e-4:
                    row.append(f"{0:.{LOG2_COMPLEXITY_PRECISION}f}")
                else:
                    row.append(f"{value:.{LOG2_COMPLEXITY_PRECISION}f}")
                rows.append(row)

        print_aligned_table(rows, left_align_columns={0})


def print_average_complexity_summary(analysis_results, args):
    """Print grouped average-complexity tables."""
    if not analysis_results:
        return

    print()
    print(f"Average-complexity C-grid steps = {args.c_grid_steps}")
    print(f"Average-complexity metric       = log2(avg complexity) with averaging={args.averaging}")

    for mode, attack, results in grouped_results(analysis_results):
        fraction_c_values = results[0]["fraction_c_values"]
        rows = [["parameter_set", "fraction_S"] + [f"c={fraction_c:.{FRACTION_PRECISION}f}" for fraction_c in fraction_c_values]]

        print()
        print(f"AVG COMPLEXITY | {mode.upper()} | {attack.upper()}")

        for result in results:
            params = result["params"]
            for row in result["rows"]:
                values = [params.name, f"{row['fraction_s']:.{FRACTION_PRECISION}f}"]
                for cell in row["values"]:
                    value = cell["log2_average_complexity"]
                    if value == float("inf"):
                        values.append("inf")
                    elif abs(value) < 5e-4:
                        values.append(f"{0:.{LOG2_COMPLEXITY_PRECISION}f}")
                    else:
                        values.append(f"{value:.{LOG2_COMPLEXITY_PRECISION}f}")
                rows.append(values)

        print_aligned_table(rows, left_align_columns={0})


def main():
    """Run the CLI entry point."""
    args = parse_args()
    params = select_params(args)
    if not params:
        raise SystemExit("No parameter sets matched the requested filters.")

    config = AnalysisConfig(
        trials=args.trials,
        s_grid_steps=args.s_grid_steps,
        c_grid_steps=args.c_grid_steps,
        base_seed=args.seed,
        mode=args.mode,
        report=args.report,
        averaging=args.averaging,
        distribution=args.distribution,
    )
    results = analyze_parameter_sets(params, config)
    if args.report in ("thresholds", "both"):
        print_threshold_summary(results["thresholds"], args)
    if args.report in ("average", "both"):
        print_average_complexity_summary(results["average_complexity"], args)


if __name__ == "__main__":
    main()
