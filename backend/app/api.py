"""FastAPI routes for the first Mr. Imputation vertical slice."""

import hashlib
import json
import os
from typing import Any, Optional

import numpy as np
from fastapi import FastAPI, HTTPException
from fastapi.middleware.cors import CORSMiddleware

from .hmm import call_alleles, run_forward_backward
from .masking import create_automatic_mask
from .schemas import (
    ExperimentRequest,
    ExperimentResponse,
    ImputeRequest,
    ImputeResponse,
)
from .simulation import simulate_experiment


app = FastAPI(title="Mr. Imputation API", version="0.1.0")
app.add_middleware(
    CORSMiddleware,
    allow_origins=[
        os.getenv("FRONTEND_ORIGIN", "http://localhost:5173"),
    ],
    allow_credentials=False,
    allow_methods=["GET", "POST"],
    allow_headers=["*"],
)


def _stable_id(prefix: str, values: Any) -> str:
    serialized = json.dumps(values, sort_keys=True, separators=(",", ":"))
    digest = hashlib.sha256(serialized.encode("utf-8")).hexdigest()[:12]
    return f"{prefix}-{digest}"


def _derive_mask_seed(seed: int, explicit_mask_seed: Optional[int]) -> int:
    if explicit_mask_seed is not None:
        return explicit_mask_seed
    sequence = np.random.SeedSequence(seed).spawn(2)
    return int(sequence[1].generate_state(1, dtype=np.uint32)[0])


@app.get("/health")
def health() -> dict:
    """Return a simple readiness response for local development."""
    return {"status": "ok", "service": "mr-imputation-api"}


@app.post("/experiments", response_model=ExperimentResponse)
def generate_experiment(request: ExperimentRequest) -> ExperimentResponse:
    """Generate synthetic data and apply an automatic observation mask."""
    mask_seed = _derive_mask_seed(request.seed, request.mask_seed)
    try:
        experiment = simulate_experiment(
            seed=request.seed,
            site_count=request.site_count,
            panel_size=request.panel_size,
            target_population=request.target_population,
            target_representation=request.target_representation,
            candidates_per_population=request.candidates_per_population,
        )
        masked_target = create_automatic_mask(
            experiment.target_haplotype.tolist(),
            observed_fraction=request.observed_fraction,
            seed=mask_seed,
        )
    except ValueError as error:
        raise HTTPException(status_code=422, detail=str(error)) from error

    identity = {
        "seed": request.seed,
        "site_count": request.site_count,
        "panel_size": request.panel_size,
        "target_population": request.target_population,
        "target_representation": request.target_representation,
        "observed_fraction": request.observed_fraction,
        "mask_seed": mask_seed,
    }
    return ExperimentResponse(
        experiment_id=_stable_id("experiment", {"seed": request.seed}),
        input_id=_stable_id("input", identity),
        seed=request.seed,
        mask_seed=mask_seed,
        target_population=request.target_population,
        site_count=experiment.site_count,
        panel_size=experiment.reference_haplotypes.shape[0],
        target_representation=request.target_representation,
        observed_fraction=request.observed_fraction,
        target_truth=list(masked_target.truth),
        observed_target=list(masked_target.hmm_target),
        observed_mask=list(masked_target.observed_mask),
        reference_haplotypes=experiment.reference_haplotypes.tolist(),
        reference_population_labels=list(experiment.reference_population_labels),
        reference_ids=list(experiment.reference_ids),
        population_counts=dict(experiment.population_counts),
        site_alleles=[list(pair) for pair in experiment.site_alleles],
    )


@app.post("/impute", response_model=ImputeResponse)
def impute(request: ImputeRequest) -> ImputeResponse:
    """Run HMM inference using only observed target data and references."""
    try:
        result = run_forward_backward(
            request.observed_target,
            request.reference_haplotypes,
            epsilon=request.epsilon,
            redraw_probability=request.redraw_probability,
        )
    except ValueError as error:
        raise HTTPException(status_code=422, detail=str(error)) from error

    return ImputeResponse(
        alpha=result.alpha,
        beta=result.beta,
        gamma=result.gamma,
        scaling_factors=result.scaling_factors,
        allele_probabilities=result.allele_probabilities,
        allele_calls=call_alleles(result.allele_probabilities),
    )
