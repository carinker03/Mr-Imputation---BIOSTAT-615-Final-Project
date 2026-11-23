"""A small, inspectable haploid reference-copying HMM.

The implementation uses Python lists rather than an API or global state so the
recurrences can be read directly and tested independently of FastAPI.
"""

from dataclasses import dataclass
from math import isfinite
from typing import List, Optional, Sequence


@dataclass
class HMMResult:
    """Scaled forward-backward arrays and posterior allele probabilities."""

    alpha: List[List[float]]
    beta: List[List[float]]
    gamma: List[List[float]]
    scaling_factors: List[float]
    allele_probabilities: List[List[float]]


def _validate_inputs(
    target: Sequence[Optional[int]],
    panel: Sequence[Sequence[int]],
    epsilon: float,
    redraw_probability: float,
) -> None:
    """Validate the bounded binary fixture inputs before running recurrences."""
    if not target:
        raise ValueError("target must contain at least one site")
    if not panel:
        raise ValueError("panel must contain at least one haplotype")
    site_count = len(target)
    if any(len(haplotype) != site_count for haplotype in panel):
        raise ValueError("every panel haplotype must match target length")
    if any(allele not in (0, 1, None) for allele in target):
        raise ValueError("target alleles must be 0, 1, or None for missing")
    if any(allele not in (0, 1) for haplotype in panel for allele in haplotype):
        raise ValueError("panel alleles must be binary")
    if not 0 <= epsilon <= 1:
        raise ValueError("epsilon must be between 0 and 1")
    if not 0 <= redraw_probability <= 1:
        raise ValueError("redraw_probability must be between 0 and 1")


def _emission(
    target_allele: Optional[int],
    panel_allele: int,
    epsilon: float,
) -> float:
    """Return the observation likelihood for one hidden reference state."""
    if target_allele is None:
        return 1.0
    return 1.0 - epsilon if target_allele == panel_allele else epsilon


def call_alleles(allele_probabilities: Sequence[Sequence[float]]) -> List[int]:
    """Choose the most probable allele, resolving exact ties to allele 0."""
    return [
        1 if probabilities[1] > probabilities[0] else 0
        for probabilities in allele_probabilities
    ]


def run_forward_backward(
    target: Sequence[Optional[int]],
    panel: Sequence[Sequence[int]],
    epsilon: float = 0.01,
    redraw_probability: float = 0.03,
) -> HMMResult:
    """Compute scaled forward-backward posteriors for a binary haploid target.

    ``None`` denotes a missing target observation. The uniform-redraw
    transition permits an O(K) update by separating the same-row term from
    the sum over all previous states.
    """
    _validate_inputs(target, panel, epsilon, redraw_probability)

    site_count = len(target)
    panel_size = len(panel)
    initial_probability = 1.0 / panel_size
    alpha = [[0.0] * panel_size for _ in range(site_count)]
    beta = [[0.0] * panel_size for _ in range(site_count)]
    scaling_factors = [0.0] * site_count

    for state in range(panel_size):
        alpha[0][state] = initial_probability * _emission(
            target[0], panel[state][0], epsilon
        )
    scaling_factors[0] = sum(alpha[0])
    if scaling_factors[0] == 0:
        raise ValueError("initial evidence has zero probability")
    alpha[0] = [value / scaling_factors[0] for value in alpha[0]]

    for site in range(1, site_count):
        previous_total = sum(alpha[site - 1])
        for state in range(panel_size):
            transition_sum = (
                (1.0 - redraw_probability) * alpha[site - 1][state]
                + redraw_probability * previous_total / panel_size
            )
            alpha[site][state] = transition_sum * _emission(
                target[site], panel[state][site], epsilon
            )
        scaling_factors[site] = sum(alpha[site])
        if scaling_factors[site] == 0:
            raise ValueError("site evidence has zero probability")
        alpha[site] = [
            value / scaling_factors[site] for value in alpha[site]
        ]

    beta[-1] = [1.0] * panel_size
    for site in range(site_count - 2, -1, -1):
        next_total = 0.0
        for next_state in range(panel_size):
            next_total += _emission(
                target[site + 1], panel[next_state][site + 1], epsilon
            ) * beta[site + 1][next_state]
        for state in range(panel_size):
            same_state_term = (
                (1.0 - redraw_probability)
                * _emission(target[site + 1], panel[state][site + 1], epsilon)
                * beta[site + 1][state]
            )
            redraw_term = redraw_probability * next_total / panel_size
            beta[site][state] = (
                same_state_term + redraw_term
            ) / scaling_factors[site + 1]

    gamma = []
    for site in range(site_count):
        unnormalized = [
            alpha[site][state] * beta[site][state]
            for state in range(panel_size)
        ]
        gamma_total = sum(unnormalized)
        gamma.append([value / gamma_total for value in unnormalized])

    allele_probabilities = []
    for site in range(site_count):
        probability_one = sum(
            gamma[site][state]
            for state in range(panel_size)
            if panel[state][site] == 1
        )
        allele_probabilities.append([1.0 - probability_one, probability_one])

    all_values = alpha + beta + gamma + [scaling_factors] + allele_probabilities
    if not all(isfinite(value) for row in all_values for value in row):
        raise ValueError("inference produced a non-finite value")

    return HMMResult(
        alpha=alpha,
        beta=beta,
        gamma=gamma,
        scaling_factors=scaling_factors,
        allele_probabilities=allele_probabilities,
    )
