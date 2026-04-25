import math
import random
from dataclasses import dataclass


@dataclass(frozen=True)
class AnalysisConfig:
    """Configuration for threshold and average-complexity analyses."""
    trials: int = 2000
    s_grid_steps: int = 2
    c_grid_steps: int = 10
    base_seed: int = 0
    mode: str = "all"
    report: str = "both"
    averaging: str = "harmonic"
    distribution: str = "random"


def q_to_nu(q):
    """Return nu = log2(q) for q = 2^nu, which is required in bit-level mode."""
    if q < 2 or q & (q - 1):
        raise ValueError(f"bit-level analysis requires q to be a power of two, got q={q}")
    return int(math.log2(q))


def supports_bit_level(params):
    """Return whether bit-level analysis applies to this parameter set."""
    return params.q > 2 and not (params.q & (params.q - 1))


def fraction_grid(grid_steps):
    """Build an inclusive uniform grid on [0, 1] with the requested number of intervals."""
    if grid_steps <= 0:
        return [0.0]
    return [index / grid_steps for index in range(grid_steps + 1)]


def capacity_s(params, mode):
    """Total number of hintable coordinates in S for the selected hint model."""
    if mode == "entry":
        return params.m * params.r
    return q_to_nu(params.q) * params.m * params.r


def capacity_c(params, mode):
    """Total number of hintable coordinates in C' for the selected hint model."""
    if mode == "entry":
        return params.r * (params.n - params.r)
    return q_to_nu(params.q) * params.r * (params.n - params.r)


def kernel_search_complexity(params, h_s, h_c, t_c):
    """Return the kernel-search base, exponent, and ell from the entry-level formula."""
    reduced_k = max(0, params.k - h_s)
    ell = max(0, min(params.n - params.r, params.n - params.r - math.ceil(reduced_k / params.m)))
    exponent = params.r * math.ceil(reduced_k / params.m) - h_c + t_c
    return params.q, exponent, ell


def kernel_search_bit_complexity(params, h_s, h_c, t_c):
    """Return the kernel-search base, exponent, and ell from the bit-level formula."""
    nu = q_to_nu(params.q)
    reduced_k = max(0.0, params.k - (h_s / nu))
    ell = max(0, min(params.n - params.r, params.n - params.r - math.ceil(reduced_k / params.m)))
    exponent = nu * params.r * math.ceil(reduced_k / params.m) - h_c + t_c
    return 2, exponent, ell


def row_space_grs_complexity(params, h_s, h_c, t_c):
    """Return the row-space GRS base, exponent, and ell from the entry-level formula."""
    ell = max(0, min(params.n - params.r, params.n - params.k - params.r + math.floor(h_s / params.m)))
    exponent = params.r * (params.k - math.floor(h_s / params.m)) - h_c + t_c
    return params.q, exponent, ell


def row_space_grs_bit_complexity(params, h_s, h_c, t_c):
    """Return the row-space GRS base, exponent, and ell from the bit-level formula."""
    nu = q_to_nu(params.q)
    ell = max(0, min(params.n - params.r, params.n - params.k - params.r + math.floor(h_s / (nu * params.m))))
    exponent = nu * params.r * (params.k - math.floor(h_s / (nu * params.m))) - h_c + t_c
    return 2, exponent, ell


def attack_for_params(params):
    """Return the attack associated with a parameter family."""
    if params.family == "mirath":
        return "minrank"
    if params.family == "ryde":
        return "rsd"
    raise ValueError(f"Unsupported parameter family: {params.family}")


def sample_t_c_random(total_capacity, num_columns, ell, h_c, rng):
    """Sample t_c when exactly h_c C' hints are placed uniformly without replacement."""
    if num_columns <= 0 or ell <= 0 or h_c <= 0:
        return 0
    column_capacity = total_capacity // num_columns
    h_c = min(h_c, total_capacity)
    column_counts = [0] * num_columns
    for position in rng.sample(range(total_capacity), h_c):
        column_counts[position // column_capacity] += 1
    return sum(sorted(column_counts)[: min(ell, num_columns)])


def sample_t_c_balanced(total_capacity, num_columns, ell, h_c):
    """Compute t_c for the most balanced placement of h_c hints across C' columns."""
    if num_columns <= 0 or ell <= 0 or h_c <= 0:
        return 0
    h_c = min(h_c, total_capacity)
    base = h_c // num_columns
    remainder = h_c % num_columns
    small_columns = max(0, num_columns - remainder)
    smallest_sum = min(ell, small_columns) * base
    if ell > small_columns:
        smallest_sum += (ell - small_columns) * (base + 1)
    return smallest_sum


def sample_trial(params, attack, mode, target_h_s, target_h_c, distribution, rng):
    """Sample one hint pattern and return the resulting exponential complexity term."""
    if attack == "minrank" and mode == "entry":
        complexity_function = kernel_search_complexity
        total_capacity = capacity_c(params, mode)
    elif attack == "minrank" and mode == "bit":
        complexity_function = kernel_search_bit_complexity
        total_capacity = capacity_c(params, mode)
    elif attack == "rsd" and mode == "entry":
        complexity_function = row_space_grs_complexity
        total_capacity = capacity_c(params, mode)
    elif attack == "rsd" and mode == "bit":
        complexity_function = row_space_grs_bit_complexity
        total_capacity = capacity_c(params, mode)
    else:
        raise ValueError(f"Unsupported attack/mode combination: {attack}/{mode}")

    num_columns = params.n - params.r
    base, exponent_without_t_c, ell = complexity_function(params, target_h_s, target_h_c, 0)

    if distribution == "random":
        actual_h_c = target_h_c
        t_c = sample_t_c_random(total_capacity, num_columns, ell, actual_h_c, rng)
    elif distribution == "balanced":
        actual_h_c = target_h_c
        t_c = sample_t_c_balanced(total_capacity, num_columns, ell, actual_h_c)
    else:
        raise ValueError(f"Unsupported distribution: {distribution}")

    exponent = exponent_without_t_c
    exponent += t_c
    return {
        "base": base,
        "exponent": exponent,
        "ell": ell,
        "h_c": actual_h_c,
        "t_c": t_c,
    }


def average_complexity(params, attack, mode, h_s, h_c, trials, seed, averaging, distribution):
    """Estimate average exponential complexity over random C' hint placements."""
    if trials <= 0:
        raise ValueError("trials must be positive")

    rng = random.Random(seed)
    total_inverse = 0.0
    total_raw = 0.0
    last_ell = None

    for _ in range(trials):
        trial = sample_trial(params, attack, mode, h_s, h_c, distribution, rng)
        exponent = trial["exponent"]
        if exponent < 0:
            raise ValueError(
                "Encountered E < 1 for a sampled hint pattern; "
                f"distribution={distribution}, exponent={exponent}, h_c={trial['h_c']}, t_c={trial['t_c']}"
            )
        e_value = trial["base"] ** exponent
        total_raw += e_value
        total_inverse += 1.0 / e_value
        last_ell = trial["ell"]

    mean_raw = total_raw / trials
    mean_inverse = total_inverse / trials

    if averaging == "harmonic":
        avg_complexity = math.inf if mean_inverse == 0.0 else 1.0 / mean_inverse
    elif averaging == "mean":
        avg_complexity = mean_raw
    else:
        raise ValueError(f"Unsupported averaging strategy: {averaging}")

    return {
        "h_s": h_s,
        "h_c": h_c,
        "ell": last_ell,
        "mean": mean_raw,
        "average_complexity": avg_complexity,
        "log2_average_complexity": math.inf if avg_complexity == math.inf else math.log2(avg_complexity),
    }


def threshold_for_fraction_s(params, attack, mode, fraction_s, config):
    """Binary-search the smallest h_c whose sampled average complexity is at most 2."""
    max_h_s = capacity_s(params, mode)
    max_h_c = capacity_c(params, mode)
    h_s = min(max_h_s, int(round(fraction_s * max_h_s)))
    seed_base = (
        config.base_seed
        + 17 * (params.m + params.n + params.k + params.r)
        + 1009 * h_s
        + (0 if attack == "minrank" else 1)
        + (0 if mode == "entry" else 1000003)
    )

    left = 0
    right = max_h_c
    best_result = None
    while left <= right:
        mid = (left + right) // 2
        result = average_complexity(
            params=params,
            attack=attack,
            mode=mode,
            h_s=h_s,
            h_c=mid,
            trials=config.trials,
            seed=seed_base + 9173 * mid,
            averaging=config.averaging,
            distribution=config.distribution,
        )
        if result["average_complexity"] <= 2.0:
            best_result = result
            right = mid - 1
        else:
            left = mid + 1

    if best_result is None:
        return {
            "fraction_s": fraction_s,
            "h_s": h_s,
            "ell": None,
            "threshold": None,
            "h_c": None,
            "mean": None,
            "average_complexity": None,
            "log2_average_complexity_at_threshold": None,
        }

    return {
        "fraction_s": fraction_s,
        "h_s": h_s,
        "ell": best_result["ell"],
        "threshold": best_result["h_c"] / max_h_c if max_h_c else 0.0,
        "h_c": best_result["h_c"],
        "mean": best_result["mean"],
        "average_complexity": best_result["average_complexity"],
        "log2_average_complexity_at_threshold": best_result["log2_average_complexity"],
    }


def average_complexity_grid(params, attack, mode, config):
    """Evaluate average complexity on the full (fraction_S, fraction_C) grid."""
    max_h_s = capacity_s(params, mode)
    max_h_c = capacity_c(params, mode)
    fraction_s_values = fraction_grid(config.s_grid_steps)
    fraction_c_values = fraction_grid(config.c_grid_steps)
    rows = []

    for fraction_s in fraction_s_values:
        h_s = min(max_h_s, int(round(fraction_s * max_h_s)))
        values = []
        for fraction_c in fraction_c_values:
            h_c = min(max_h_c, int(round(fraction_c * max_h_c)))
            seed = (
                config.base_seed
                + 17 * (params.m + params.n + params.k + params.r)
                + 1009 * h_s
                + 9173 * h_c
                + (0 if attack == "minrank" else 1)
                + (0 if mode == "entry" else 1000003)
            )
            result = average_complexity(
                params=params,
                attack=attack,
                mode=mode,
                h_s=h_s,
                h_c=h_c,
                trials=config.trials,
                seed=seed,
                averaging=config.averaging,
                distribution=config.distribution,
            )
            values.append(
                {
                    "fraction_c": fraction_c,
                    "h_c": h_c,
                    "mean": result["mean"],
                    "average_complexity": result["average_complexity"],
                    "log2_average_complexity": result["log2_average_complexity"],
                }
            )
        rows.append({"fraction_s": fraction_s, "h_s": h_s, "values": values})

    return {
        "params": params,
        "attack": attack,
        "mode": mode,
        "fraction_c_values": fraction_c_values,
        "rows": rows,
    }


def analyze_parameter_sets(params_list, config):
    """Run the requested analyses for each parameter set."""
    threshold_results = []
    average_results = []

    for params in params_list:
        attack = attack_for_params(params)
        if config.mode == "all":
            modes = ["entry"]
            if supports_bit_level(params):
                modes.append("bit")
        elif config.mode == "bit":
            if not supports_bit_level(params):
                continue
            modes = ["bit"]
        else:
            modes = [config.mode]

        for mode in modes:
            threshold_results.append(
                {
                    "params": params,
                    "attack": attack,
                    "mode": mode,
                    "points": [
                        threshold_for_fraction_s(params, attack, mode, fraction_s, config)
                        for fraction_s in fraction_grid(config.s_grid_steps)
                    ],
                }
            )
            if config.report in ("average", "both"):
                average_results.append(average_complexity_grid(params, attack, mode, config))

    return {
        "thresholds": threshold_results,
        "average_complexity": average_results,
    }
