"""Deterministic synthetic haplotype generation for Imputation Explorer.

This is a small educational simulator, not a population-genetics package.
Haplotypes copy founder patterns within blocks, while population-specific
founder weights and low mutation rates create shared but distinguishable
synthetic populations.
"""

from dataclasses import dataclass
from typing import Dict, List, Mapping, Optional, Sequence, Tuple

import numpy as np


POPULATIONS = ("A", "B", "C")
DEFAULT_SEED = 615
DEFAULT_SITE_COUNT = 400
DEFAULT_PANEL_SIZE = 40
DEFAULT_CANDIDATES_PER_POPULATION = 100
FOUNDER_COUNT = 6
DEFAULT_BLOCK_SIZE = 20
DEFAULT_MUTATION_RATE = 0.02

# These weights deliberately overlap: population labels influence sampling,
# but no population receives a unique algorithmic advantage.
POPULATION_FOUNDER_WEIGHTS = {
    "A": np.array([0.35, 0.25, 0.15, 0.10, 0.10, 0.05]),
    "B": np.array([0.10, 0.20, 0.30, 0.20, 0.15, 0.05]),
    "C": np.array([0.15, 0.10, 0.10, 0.20, 0.30, 0.15]),
}


@dataclass(frozen=True)
class SyntheticExperiment:
    """One reproducible synthetic target and its reference candidate pools."""

    seed: int
    target_population: str
    site_count: int
    block_size: int
    mutation_rate: float
    target_haplotype: np.ndarray
    reference_haplotypes: np.ndarray
    reference_population_labels: Tuple[str, ...]
    reference_ids: Tuple[str, ...]
    candidate_haplotypes: np.ndarray
    candidate_population_labels: Tuple[str, ...]
    candidate_ids: Tuple[str, ...]
    site_alleles: Tuple[Tuple[str, str], ...]
    population_counts: Mapping[str, int]


def _validate_configuration(
    seed: int,
    site_count: int,
    panel_size: int,
    target_population: str,
    target_representation: float,
    candidates_per_population: int,
    block_size: int,
    mutation_rate: float,
) -> None:
    if not isinstance(seed, (int, np.integer)):
        raise ValueError("seed must be an integer")
    if site_count <= 0:
        raise ValueError("site_count must be positive")
    if panel_size <= 0:
        raise ValueError("panel_size must be positive")
    if target_population not in POPULATIONS:
        raise ValueError("target_population must be one of A, B, or C")
    if not 0.0 <= target_representation <= 1.0:
        raise ValueError("target_representation must be between 0 and 1")
    if candidates_per_population <= 0:
        raise ValueError("candidates_per_population must be positive")
    if panel_size > candidates_per_population * len(POPULATIONS):
        raise ValueError("panel_size exceeds the available candidate pool")
    quotas = _panel_counts(
        panel_size, target_population, target_representation
    )
    if max(quotas.values()) > candidates_per_population:
        raise ValueError(
            "a population quota exceeds its available candidate haplotypes"
        )
    if block_size <= 0:
        raise ValueError("block_size must be positive")
    if not 0.0 <= mutation_rate <= 1.0:
        raise ValueError("mutation_rate must be between 0 and 1")


def _panel_counts(
    panel_size: int,
    target_population: str,
    target_representation: float,
) -> Dict[str, int]:
    """Return deterministic integer quotas using a largest-remainder rule."""
    target_count = min(
        panel_size,
        int(np.floor(panel_size * target_representation + 0.5)),
    )
    remaining = panel_size - target_count
    other_populations = [
        population
        for population in POPULATIONS
        if population != target_population
    ]
    first_other_count = remaining // 2 + remaining % 2
    counts = {population: 0 for population in POPULATIONS}
    counts[target_population] = target_count
    counts[other_populations[0]] = first_other_count
    counts[other_populations[1]] = remaining // 2
    return counts


def _make_founder_blocks(
    rng: np.random.Generator,
    site_count: int,
    block_size: int,
    mutation_rate: float,
) -> np.ndarray:
    """Create shared founder patterns with strong within-block linkage."""
    block_count = int(np.ceil(site_count / block_size))
    founders = np.empty((FOUNDER_COUNT, site_count), dtype=np.int8)
    for founder in range(FOUNDER_COUNT):
        for block in range(block_count):
            start = block * block_size
            stop = min(start + block_size, site_count)
            block_base = rng.integers(0, 2)
            founders[founder, start:stop] = block_base
            mutations = rng.random(stop - start) < mutation_rate
            founders[founder, start:stop][mutations] = 1 - block_base
    return founders


def _sample_haplotype(
    rng: np.random.Generator,
    founders: np.ndarray,
    population: str,
    block_size: int,
) -> np.ndarray:
    """Sample one chromosome by choosing a founder independently per block."""
    site_count = founders.shape[1]
    haplotype = np.empty(site_count, dtype=np.int8)
    weights = POPULATION_FOUNDER_WEIGHTS[population]
    for start in range(0, site_count, block_size):
        stop = min(start + block_size, site_count)
        founder = rng.choice(FOUNDER_COUNT, p=weights)
        haplotype[start:stop] = founders[founder, start:stop]
    return haplotype


def _site_alleles(site_count: int) -> Tuple[Tuple[str, str], ...]:
    """Assign varied fictional nucleotide pairs to synthetic site indices."""
    pairs = (("A", "G"), ("C", "T"), ("A", "C"), ("G", "T"))
    return tuple(pairs[site % len(pairs)] for site in range(site_count))


def simulate_experiment(
    seed: int = DEFAULT_SEED,
    site_count: int = DEFAULT_SITE_COUNT,
    panel_size: int = DEFAULT_PANEL_SIZE,
    target_population: str = "A",
    target_representation: float = 0.5,
    candidates_per_population: int = DEFAULT_CANDIDATES_PER_POPULATION,
    block_size: int = DEFAULT_BLOCK_SIZE,
    mutation_rate: float = DEFAULT_MUTATION_RATE,
) -> SyntheticExperiment:
    """Generate a reproducible synthetic target and selected reference panel.

    The returned arrays use shape ``(number_of_haplotypes, site_count)``.
    Therefore ``reference_haplotypes[j, t]`` is the allele at site ``t`` for
    reference row ``j``, matching ``run_forward_backward``.
    """
    _validate_configuration(
        seed,
        site_count,
        panel_size,
        target_population,
        target_representation,
        candidates_per_population,
        block_size,
        mutation_rate,
    )
    stream_sequence = np.random.SeedSequence(int(seed))
    founder_sequence, target_sequence, reference_sequence, ordering_sequence = (
        stream_sequence.spawn(4)
    )
    founder_rng = np.random.default_rng(founder_sequence)
    target_rng = np.random.default_rng(target_sequence)
    reference_rng = np.random.default_rng(reference_sequence)
    ordering_rng = np.random.default_rng(ordering_sequence)

    founders = _make_founder_blocks(
        founder_rng, site_count, block_size, mutation_rate
    )
    target_haplotype = _sample_haplotype(
        target_rng, founders, target_population, block_size
    )

    candidate_haplotypes = []
    candidate_labels = []
    candidate_ids = []
    ordering = {}
    for population in POPULATIONS:
        ordering[population] = ordering_rng.permutation(candidates_per_population)
        for candidate_number in range(candidates_per_population):
            candidate_haplotypes.append(
                _sample_haplotype(
                    reference_rng, founders, population, block_size
                )
            )
            candidate_labels.append(population)
            candidate_ids.append(f"{population}-{candidate_number:03d}")

    candidate_haplotypes_array = np.asarray(candidate_haplotypes, dtype=np.int8)
    counts = _panel_counts(
        panel_size, target_population, target_representation
    )
    selected_indices = []
    for population in POPULATIONS:
        population_offset = POPULATIONS.index(population) * candidates_per_population
        selected_indices.extend(
            population_offset + ordering[population][: counts[population]]
        )
    selected_indices_array = np.asarray(selected_indices, dtype=int)

    return SyntheticExperiment(
        seed=int(seed),
        target_population=target_population,
        site_count=site_count,
        block_size=block_size,
        mutation_rate=mutation_rate,
        target_haplotype=target_haplotype,
        reference_haplotypes=candidate_haplotypes_array[selected_indices_array],
        reference_population_labels=tuple(
            candidate_labels[index] for index in selected_indices
        ),
        reference_ids=tuple(candidate_ids[index] for index in selected_indices),
        candidate_haplotypes=candidate_haplotypes_array,
        candidate_population_labels=tuple(candidate_labels),
        candidate_ids=tuple(candidate_ids),
        site_alleles=_site_alleles(site_count),
        population_counts=counts,
    )
