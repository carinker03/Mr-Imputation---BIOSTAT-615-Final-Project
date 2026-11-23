"""Deterministic observation masks for synthetic target haplotypes."""

from dataclasses import dataclass
from math import isfinite
from typing import Iterable, List, Optional, Sequence, Tuple

import numpy as np


DEFAULT_MASK_SEED = 615
DEFAULT_OBSERVED_FRACTION = 0.40


@dataclass(frozen=True)
class MaskedTarget:
    """A target truth paired with an immutable input mask and HMM view.

    ``observed_mask[t]`` is true when the target allele is available to the
    inference function. ``hmm_target`` contains the allele at observed sites
    and ``None`` at missing sites.
    """

    truth: Tuple[int, ...]
    observed_mask: Tuple[bool, ...]
    hmm_target: Tuple[Optional[int], ...]
    mode: str
    seed: Optional[int]
    requested_observed_fraction: Optional[float]
    observed_count: int
    missing_count: int


def _validate_target(target: Sequence[int]) -> Tuple[int, ...]:
    """Validate and copy a binary target so callers' data cannot be changed."""
    values = tuple(target)
    if not values:
        raise ValueError("target must contain at least one SNP")
    if any(
        isinstance(allele, bool)
        or not isinstance(allele, (int, np.integer))
        or allele not in (0, 1)
        for allele in values
    ):
        raise ValueError("target must contain only binary alleles 0 and 1")
    return tuple(int(allele) for allele in values)


def _validate_seed(seed: int) -> int:
    if not isinstance(seed, (int, np.integer)):
        raise ValueError("seed must be an integer")
    return int(seed)


def _validate_positions(
    positions: Iterable[int],
    site_count: int,
    name: str,
) -> Tuple[int, ...]:
    values = tuple(positions)
    if any(
        isinstance(position, bool)
        or not isinstance(position, (int, np.integer))
        or position < 0
        or position >= site_count
        for position in values
    ):
        raise ValueError(
            f"{name} must contain unique positions between 0 and {site_count - 1}"
        )
    if len(set(values)) != len(values):
        raise ValueError(f"{name} must not contain duplicate positions")
    return tuple(int(position) for position in values)


def _build_result(
    truth: Tuple[int, ...],
    observed_mask: Sequence[bool],
    mode: str,
    seed: Optional[int],
    requested_observed_fraction: Optional[float],
) -> MaskedTarget:
    mask = tuple(bool(value) for value in observed_mask)
    hmm_target = tuple(
        allele if observed else None
        for allele, observed in zip(truth, mask)
    )
    observed_count = sum(mask)
    return MaskedTarget(
        truth=truth,
        observed_mask=mask,
        hmm_target=hmm_target,
        mode=mode,
        seed=seed,
        requested_observed_fraction=requested_observed_fraction,
        observed_count=observed_count,
        missing_count=len(truth) - observed_count,
    )


def create_automatic_mask(
    target: Sequence[int],
    observed_fraction: float = DEFAULT_OBSERVED_FRACTION,
    seed: int = DEFAULT_MASK_SEED,
) -> MaskedTarget:
    """Reveal the highest fixed random-priority sites at the requested density."""
    truth = _validate_target(target)
    if not isinstance(observed_fraction, (int, float, np.number)):
        raise ValueError("observed_fraction must be numeric")
    observed_fraction = float(observed_fraction)
    if not isfinite(observed_fraction) or not 0.0 <= observed_fraction <= 1.0:
        raise ValueError("observed_fraction must be between 0 and 1")
    seed = _validate_seed(seed)

    requested_count = int(np.floor(len(truth) * observed_fraction + 0.5))
    priorities = np.random.default_rng(seed).random(len(truth))
    observed_indices = np.argsort(-priorities, kind="stable")[:requested_count]
    observed_mask = [False] * len(truth)
    for index in observed_indices:
        observed_mask[int(index)] = True
    return _build_result(
        truth,
        observed_mask,
        mode="automatic",
        seed=seed,
        requested_observed_fraction=observed_fraction,
    )


def create_manual_mask(
    target: Sequence[int],
    observed_positions: Iterable[int] = (),
    missing_positions: Iterable[int] = (),
) -> MaskedTarget:
    """Create an exact manual mask, initially observing every target site."""
    truth = _validate_target(target)
    observed_positions = _validate_positions(
        observed_positions, len(truth), "observed_positions"
    )
    missing_positions = _validate_positions(
        missing_positions, len(truth), "missing_positions"
    )
    if set(observed_positions).intersection(missing_positions):
        raise ValueError(
            "a position cannot be both observed and missing"
        )
    observed_mask = [True] * len(truth)
    for position in missing_positions:
        observed_mask[position] = False
    for position in observed_positions:
        observed_mask[position] = True
    return _build_result(
        truth,
        observed_mask,
        mode="manual",
        seed=None,
        requested_observed_fraction=None,
    )


def update_manual_mask(
    masked_target: MaskedTarget,
    observed_positions: Iterable[int] = (),
    missing_positions: Iterable[int] = (),
) -> MaskedTarget:
    """Mask or unmask individual sites without changing the target truth."""
    truth = _validate_target(masked_target.truth)
    observed_positions = _validate_positions(
        observed_positions, len(truth), "observed_positions"
    )
    missing_positions = _validate_positions(
        missing_positions, len(truth), "missing_positions"
    )
    if set(observed_positions).intersection(missing_positions):
        raise ValueError(
            "a position cannot be both observed and missing"
        )
    observed_mask = list(masked_target.observed_mask)
    for position in missing_positions:
        observed_mask[position] = False
    for position in observed_positions:
        observed_mask[position] = True
    return _build_result(
        truth,
        observed_mask,
        mode="manual",
        seed=masked_target.seed,
        requested_observed_fraction=None,
    )
