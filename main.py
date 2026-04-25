import argparse
from types import SimpleNamespace

from simulation import AnalysisConfig, scan_parameter_sets


DEFAULT_TRIALS = 2000
DEFAULT_S_GRID_STEPS = 10
DEFAULT_C_GRID_STEPS = 10

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
    parser = argparse.ArgumentParser(
        description=(
            "Estimate the fraction of hints needed for polynomial-time kernel search "
            "(MinRank) and GRS row-space (RSD) under random hint placement."
        )
    )
    parser.add_argument(
        "--family",
        choices=("all", "mirath", "ryde"),
        default="all",
        help="Restrict the default scan to one preset family.",
    )
    parser.add_argument(
        "--name",
        action="append",
        default=[],
        help="Analyze only the named preset(s). Can be passed multiple times.",
    )
    parser.add_argument(
        "--attack",
        choices=("all", "minrank", "rsd"),
        default="all",
        help="Restrict which attack formulas are evaluated.",
    )
    parser.add_argument(
        "--mode",
        choices=("all", "entry", "bit"),
        default="all",
        help="Choose whether to run entry-level, bit-level, or both models.",
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
        help="Number of intervals for the C-hint fraction grid between 0 and 1 in the average-complexity report; use 0 to run only fraction_C = 0.",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=0,
        help="Base RNG seed for repeatable Monte Carlo estimates.",
    )
    parser.add_argument(
        "--show-average-complexity",
        action="store_true",
        help="Also print a separate grid of log2 average complexity under the selected averaging strategy.",
    )
    parser.add_argument(
        "--averaging",
        choices=("reciprocal", "mean_e"),
        default="reciprocal",
        help="How to average sampled complexities: reciprocal means 1 / average(1/E), mean_e means average(E).",
    )
    parser.add_argument(
        "--distribution",
        choices=("fixed", "balanced"),
        default="fixed",
        help="How to sample C' hints: fixed total without replacement, or as balanced across columns as possible.",
    )
    return parser.parse_args()

def select_params(args):
    params = ALL_PARAMS
    if args.family != "all":
        params = [param for param in params if param.family == args.family]
    if args.name:
        wanted = set(args.name)
        params = [param for param in params if param.name in wanted]
    return params


def format_fraction(value):
    if value is None:
        return "infeasible"
    return f"{value:.3f}"


def print_threshold_summary(analysis_results, args):
    print(f"Monte Carlo trials = {args.trials}")
    print(f"S-grid steps       = {args.s_grid_steps}")
    print(f"Base seed          = {args.seed}")
    print(f"Averaging          = {args.averaging}")
    print(f"Distribution       = {args.distribution}")

    for result in analysis_results:
        params = result["params"]
        print()
        print(
            f"{result['attack'].upper()} | {result['mode'].upper()} | {params.name} "
            f"(q={params.q}, m={params.m}, n={params.n}, k={params.k}, r={params.r})"
        )
        print("fraction_S  min_fraction_C  h_c  log2_avg_complexity  p_avg_at_threshold")
        for point in result["points"]:
            threshold = point["threshold"]
            threshold_str = format_fraction(threshold)
            h_c_str = "n/a" if threshold is None else str(point["h_c"])
            avg_str = "n/a" if threshold is None else f"{point['log2_average_complexity']:.3f}"
            probability_str = "n/a" if threshold is None else f"{point['p_avg']:.6f}"
            print(f"{point['fraction_s']:.3f}       {threshold_str:>12}  {h_c_str:>3}  {avg_str:>19}  {probability_str}")


def print_average_complexity_summary(analysis_results, args):
    if not analysis_results:
        return

    print()
    print(f"Average-complexity C-grid steps = {args.c_grid_steps}")
    print(f"Average-complexity metric       = log2(avg complexity) with averaging={args.averaging}")

    for result in analysis_results:
        params = result["params"]
        print()
        print(
            f"AVG COMPLEXITY | {result['attack'].upper()} | {result['mode'].upper()} | {params.name} "
            f"(q={params.q}, m={params.m}, n={params.n}, k={params.k}, r={params.r})"
        )
        header = ["fraction_S"] + [f"c={fraction_c:.3f}" for fraction_c in result["fraction_c_values"]]
        print("  ".join(header))
        for row in result["rows"]:
            values = [f"{row['fraction_s']:.3f}"]
            for cell in row["values"]:
                value = cell["log2_average_complexity"]
                if value == float("inf"):
                    values.append("inf")
                elif abs(value) < 5e-4:
                    values.append("0.000")
                else:
                    values.append(f"{value:.3f}")
            print("  ".join(values))


def main():
    args = parse_args()
    params = select_params(args)
    if not params:
        raise SystemExit("No parameter sets matched the requested filters.")

    config = AnalysisConfig(
        trials=args.trials,
        s_grid_steps=args.s_grid_steps,
        c_grid_steps=args.c_grid_steps,
        base_seed=args.seed,
        attack=args.attack,
        mode=args.mode,
        show_average_complexity=args.show_average_complexity,
        averaging=args.averaging,
        distribution=args.distribution,
    )
    results = scan_parameter_sets(params, config)
    print_threshold_summary(results["thresholds"], args)
    if args.show_average_complexity:
        print_average_complexity_summary(results["average_complexity"], args)


if __name__ == "__main__":
    main()
