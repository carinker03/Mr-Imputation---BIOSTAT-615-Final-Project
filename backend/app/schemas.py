"""Validated request and response models for the local Mr. Imputation API."""

from typing import List, Literal, Optional

from pydantic import BaseModel, Field


class ExperimentRequest(BaseModel):
    """Settings used to generate one deterministic synthetic experiment."""

    seed: int = Field(default=615)
    site_count: int = Field(default=40, ge=1, le=400)
    panel_size: int = Field(default=40, ge=1, le=300)
    target_population: Literal["A", "B", "C"] = "A"
    target_representation: float = Field(default=0.5, ge=0.0, le=1.0)
    observed_fraction: float = Field(default=0.4, ge=0.0, le=1.0)
    candidates_per_population: int = Field(default=100, ge=1, le=1000)
    mask_seed: Optional[int] = None


class ExperimentResponse(BaseModel):
    """JSON-safe synthetic data and its exact generated settings."""

    experiment_id: str
    input_id: str
    seed: int
    mask_seed: int
    target_population: str
    site_count: int
    panel_size: int
    target_representation: float
    observed_fraction: float
    target_truth: List[int]
    observed_target: List[Optional[int]]
    observed_mask: List[bool]
    reference_haplotypes: List[List[int]]
    reference_population_labels: List[str]
    reference_ids: List[str]
    population_counts: dict
    site_alleles: List[List[str]]


class ImputeRequest(BaseModel):
    """Observed target and selected reference panel sent to the HMM."""

    observed_target: List[Optional[int]]
    reference_haplotypes: List[List[int]]
    epsilon: float = Field(default=0.01, ge=0.0, le=1.0)
    redraw_probability: float = Field(default=0.03, ge=0.0, le=1.0)


class ImputeResponse(BaseModel):
    """HMM posterior arrays and deterministic allele calls."""

    alpha: List[List[float]]
    beta: List[List[float]]
    gamma: List[List[float]]
    scaling_factors: List[float]
    allele_probabilities: List[List[float]]
    allele_calls: List[int]
