# Imputation Explorer: implementation specification

Addressed to the coding agent building Clara Rinker's BIOSTAT 615 final project.
Prepared October 26, 2026, from the completed project interview. Working title; a final product name is not a prerequisite.

## 1. Purpose and audience

Build a public, free educational web application that makes genotype imputation visible and explains how reference-panel representation, panel size, and observed input information affect predictions. The primary audience is students with basic probability and programming knowledge who are beginning statistical genetics. Visual explanations should also welcome beginners; optional mathematical detail should support graduate statistical computing instruction.

The central experience is a controlled synthetic experiment: inspect known truth, mask selected variants, impute using different panels, and compare predictions and uncertainty. Users must be able to inspect the actual dynamic programming calculations that generated their results.

Use a playful science-explainer introduction inspired by the owner's interest in the Mr. DNA sequence, with an original identity and a precise, uncluttered workspace. Do not reuse film characters, audio, or artwork. This is an educational simulator, not a clinical tool or production imputation service.

## 2. Commitments and scope

The core implementation is haploid: one chromosome copy at selected biallelic SNP positions. The three aligned tracks are ground truth, masked input, and imputed output. A separate reference haplotype grid makes the panel visible. A locus is an existing aligned position, not a missing piece of physical DNA; imputation predicts an allele rather than inserting sequence. Columns are selected variant sites, not every DNA base.

Call collections of samples reference panels, not reference genomes. Explain that two homologous chromosome copies inherited from parents are different from the complementary strands of one DNA double helix. Diploid targets are an extension, not required.

Use entirely synthetic populations for the complete deliverable. Real public data, admixture, repeated-target simulations, and adjustable HMM parameters are extensions. Never use restricted research data or real-person genotype uploads.

### Checkpoint: November 23, 2026

Deliver a coherent, independently usable public application with:

- Short, skippable, replayable text-guided introduction and a basic About/research context section.
- Three fictional populations; one small target region by default.
- Python simulation, HMM forward-backward inference, frequency baseline, and evaluation.
- Three aligned tracks and a reference haplotype grid.
- Panel representation presets and slider, advanced panel-size control, and automatic masking presets and slider.
- Forward and backward step-through with selectable probability cells, worked arithmetic, and posterior allele probabilities.
- Accuracy, uncertainty, counts, and factual explanations.
- Clear experiment/input/run identity in the data model, even if persistent history is added afterward.
- Tested public deployment, plus a locally runnable presentation backup.

### Post-checkpoint priorities, in agreed order

1. Browser-persistent experiment history, controlled comparisons, JSON export/import, and CSV export.
2. Manual masking, including individual sites and contiguous dragged ranges.
3. Longer-region mode with summaries and a detailed inspection window.
4. Animated, narrated version of the introduction, retaining text accessibility.

These are planned final-version features, implemented sequentially within the time budget. Optional prediction-confidence filtering follows the core results and history work; it must not displace checkpoint correctness. Preserve a complete text-guided experience if animation cannot be finished. Document any deferred feature honestly.

Final submission deadline: Thursday, December 17, 2026, 3:30 PM America/Detroit (ET). The owner will produce a 5–10 minute video explaining the application, computational course connection, and usefulness. Target a December 10 feature freeze.

## 3. User journey and interface

Open with a prepared reproducible example and a short guided story: known variants, concealed variants, reference chromosomes, probabilistic prediction, evaluation against held-out truth. Explain that truth is available only because this is a simulation. Skip and replay must always work.

The workspace contains:

1. Experiment identity and target population, plus a deliberate New experiment action.
2. Input subexperiment controls: region, observed-site density, and later automatic/manual masking.
3. Reference controls: representation and size, with counts by population.
4. Run imputation action, clearly indicating changed settings and stale results.
5. Aligned truth/input/output tracks with shared position selection.
6. Reference grid and prediction probabilities.
7. Results and comparison history.
8. Expandable HMM inspector and accessible About/glossary/research content.

No setting change silently generates a new underlying target. Replaying animation never reruns the statistical model. Draft controls do not overwrite completed runs. Use an explicit Run action so slider movement does not flood the backend or history. Retain results with their exact settings; show a changed-settings indicator when needed.

### Defaults and proposed numerical limits

These are implementation defaults, not empirical claims about chip coverage. Adjust for measured performance while retaining the intended distinctions:

| Control | Default | Proposed range/presets |
| --- | --- | --- |
| Teaching region | 40 SNP sites | Fixed 40-site prefix initially |
| Longer region | 200 SNP sites | 100, 200, 400; cap 400 |
| Reference panel | 40 haplotypes | 10–100, integer steps |
| Target-population representation | 50% | 0–100%; presets 10%, 50%, 90% |
| Observed sites | 40% | 0–100%; presets 15%, 40%, 80% |
| Target population | Fictional population A | A, B, C |
| Emission mismatch probability | 0.01 | Fixed in core implementation |
| HMM redraw probability per adjacent site | 0.03 | Fixed in core implementation |

Show realized counts and percentages when rounding to whole haplotypes or sites. Use neutral labels such as low, moderate, and high representation; never suggest EUR or EAS are universal population-generic panels. The non-target fraction is split deterministically between the other two populations, with the rounding rule recorded.

Observed-site presets mean sparse, moderate, and dense observations in this toy model. Do not label the middle preset an empirically typical chip percentage. Missingness is absence of a measured allele, not deletion of genomic material.

### Visual encoding

Use consistent colors for A/C/G/T with letters always visible. Missing alleles use a gray dashed tile and an accessible missing label. Correctness uses a separate outline or icon, never a reassignment of allele colors. Confidence uses a separate probability bar or annotation. Include a legend.

Use aligned rectangular tracks rather than a twisting helix where positional comparisons would become harder. A decorative helix may appear in the introduction. Keyboard selection must provide the same site details as hover/click. Respect reduced-motion preferences and ensure readable contrast.

The reference grid shows all rows for a small panel or a scrollable/virtualized selection for larger panels. Always show full inference panel size and displayed row count; the visible subset must not silently become the statistical panel. Label fictional population membership. Highlight evidence and posterior reference support, not merely visual matches.

## 4. Experiment hierarchy and reproducibility

Use three explicit levels:

- **Experiment:** an immutable master seed, simulation version, fictional target population, complete synthetic target region, and reference candidate pools. New experiment starts another independent target; a deliberate target-population change also creates a new experiment. Retain the previous one in history.
- **Input subexperiment:** the same underlying target with a particular region length and exact observation mask. Mask or region changes create a distinct input subexperiment. Clearly label these as different inputs.
- **Run:** a selected reference composition, panel size, fixed or advanced model parameters, and computed outputs for one input subexperiment.

Default comparisons are between runs within the same input subexperiment. Across-subexperiment or across-experiment comparisons require explicit selection and list all changed inputs. Do not attribute differences to representation alone when size, mask, region, model, or target also changed.

Use separate deterministic random streams for population founders, target generation, reference generation, candidate ordering, and masking. Changing one control must not consume another stream and regenerate unrelated data. Stable candidate IDs and per-population orderings support nested sampling. Increasing panel size at unchanged composition reuses eligible reference chromosomes as far as integer quotas permit. Record actual membership; composition changes necessarily replace some members.

Generate the complete capped target region and full candidate pools once per experiment. Smaller regions are prefixes of this region. Automatic masking assigns fixed random priorities to sites; at a selected region length, reveal the highest-priority requested count. Thus density changes are nested for a fixed region. On region expansion, preserve every existing site's observation state and assign new sites using the stored rule; the realized overall density may differ slightly. Display actual density and record this behavior. Manual masks are exact stored boolean arrays; returning to automatic mode creates another input subexperiment rather than silently editing a prior one.

There is no globally unchanged set of reference samples when composition changes; explain this rather than promising to isolate composition perfectly. Stable pools reduce incidental resampling, but comparisons remain examples from the chosen synthetic data.

## 5. Synthetic data model

Build a lightweight simulator, not a full population-genetics package. Use a shared collection of ancestral block haplotypes, population-specific block-frequency weights, occasional population-enriched alleles, and recombination between blocks. Each chromosome is sampled independently from this generative model. Adjacent sites must share meaningful linkage patterns; independently sampling each site's allele is insufficient for the intended HMM lesson.

Suggested construction: a 400-site region partitioned into variable blocks; shared binary founder patterns per block; population-specific founder weights with overlapping support; low-rate independent mutations. Generate independent reference candidate chromosomes and target chromosomes using distinct streams. Freeze founder structure and parameters under a versioned configuration. Targets must not be sampled by copying a complete realized reference individual. Sharing simulated ancestry is intentional; leaking held-out target alleles is not.

Assign a fixed pair of distinct nucleotide alleles to each site and encode computations as 0/1. All samples share the same site catalog and allele orientation. Include varied nucleotide pairs across sites for visual diversity. Site numbers are synthetic indices, not real genomic coordinates or rsIDs. Sites may become monomorphic in a selected panel; handle this explicitly.

Simulation code produces ground truth before masking and does not invoke imputation. The inference function accepts only observed target alleles, the mask, panel alleles, and HMM parameters. Evaluation receives truth separately. Do not select or alter examples dynamically to force matched panels to win. Use a documented default seed and allow explicit new experiments. Validate that informative haplotype structure exists and report variability honestly; a small experiment can favor any panel.

A selected panel may lack a true target allele entirely. The model below can then be confidently wrong. Show this limitation and do not recover the hidden allele from truth.

## 6. Statistical method

Implement a simplified haploid reference-copying HMM. This is an educational model, not a reimplementation of a production server. The course connection is dynamic programming over aligned loci, not minimum edit distance, sequence alignment, or insertion/deletion costs.

Let L be the number of sites and K the selected reference haplotypes. Reference allele h[j,t] is binary. Hidden state Z[t] identifies the reference row matching the target locally. Initial state probability is uniform, pi[j] = 1/K. With redraw probability rho[t], define:

```text
T_t(i,j) = (1 - rho[t]) * I(i = j) + rho[t] / K
```

A redraw can select the same row, so rho is not literally the probability of changing row. In the core model, rho is constant per adjacent synthetic SNP. Explain that it is a simplified model parameter, not an estimated biological recombination rate; do not imply a real physical region length.

For mismatch parameter epsilon:

```text
e_t(j) = 1                         if target site t is missing
         1 - epsilon               if observed allele equals h[j,t]
         epsilon                   otherwise
```

Missing observations contribute no allele preference. Epsilon provides a mismatch tolerance; simulated genotyping error is not required in the first version.

### Scaled forward-backward algorithm

Use float64 and consistent per-site scaling to prevent numerical underflow. The following is the canonical recurrence, using sites 0 through L-1:

```text
u_0(j) = pi[j] * e_0(j)
c_0 = sum_j u_0(j)
alpha_0(j) = u_0(j) / c_0

u_t(j) = e_t(j) * sum_i [alpha_(t-1)(i) * T_t(i,j)]
c_t = sum_j u_t(j)
alpha_t(j) = u_t(j) / c_t

beta_(L-1)(j) = 1
beta_t(i) = sum_j [T_(t+1)(i,j) * e_(t+1)(j) * beta_(t+1)(j)] / c_(t+1)

gamma_t(j) = alpha_t(j) * beta_t(j)
gamma_t = gamma_t / sum_j gamma_t(j)

P_t(allele = a) = sum_j [gamma_t(j) * I(h[j,t] = a)]
```

For the uniform-redraw transition, implement the forward and backward sums in O(K) per site rather than materializing an L by K by K tensor. Preserve a way to reconstruct an individual transition's contribution for explanation. Total core inference is O(LK) time and O(LK) memory when retaining instructional arrays. A naive O(LK²) reference implementation is useful only in small verification fixtures.

The allele posterior is a weighted combination of panel alleles; no mutation is added at missing sites in this first model. Predictions are marginal maximum-probability alleles, not a Viterbi path. Resolve exact ties deterministically to allele 0 and document the tie. Preserve observed input alleles in the output track; imputation fills masked sites only.

An optional Viterbi overlay may be added later. Label it the single most likely reference-state path, which differs from the marginal probabilities used for imputation. Do not make Viterbi required for the checkpoint.

### Required HMM inspector

Provide forward/backward navigation, play/pause, phase labels, and final posterior view. Forward steps show normalized filtering probabilities based on evidence through that site. Backward steps show likelihood contributions from later evidence; beta values are not a normalized state-probability distribution. Final gamma combines both directions.

Selecting a cell shows its definition, symbolic recurrence, actual substituted values, normalization factor, and result. Expand predecessor contributions on demand rather than rendering all arithmetic at once. Explain emission 1 at missing sites. Display full-precision-derived values with enough readable digits, clearly separating rounded text from the actual computation. Include initial and terminal boundary cases.

The same site selection must connect the input tile, reference column, probability cell, and resulting allele probabilities. In longer mode, a detail window shows calculations from inference over the entire selected region, including evidence outside the visible window. An independently run short-region experiment may yield different probabilities; label that distinction.

## 7. Baseline, metrics, and interpretation

At each masked site the baseline uses only the allele frequencies in the selected panel:

```text
baseline P_t(1) = count of panel haplotypes carrying allele 1 / K
baseline call = argmax over alleles, with the same tie rule
```

It ignores both target observations and linkage context. Evaluate the baseline and HMM on exactly the same masked sites and panel. The HMM can do worse in an individual example.

For the masked-site set M, report:

- Number of evaluated sites, correct calls, and incorrect calls.
- Accuracy = correct calls / |M|.
- Per-site predicted allele probability and confidence = max(P(0), P(1)).
- Mean confidence, explicitly distinguished from empirical accuracy.
- Mean binary Brier score = mean over M of (P_t(1) - truth_t)^2; lower is better.
- Side-by-side HMM and baseline results, with accuracy differences in percentage points.

Keep advanced metrics collapsed by default. Do not call haploid accuracy genotype concordance without a qualification. Do not present correlation across a short target region as standard per-variant dosage r² across individuals. Per-call posterior probability, empirical accuracy, and server-reported variant quality are distinct concepts. Discuss the distinction in help text, not as a claim that a toy threshold reproduces an industry standard.

If confidence filtering is implemented, show unresolved sites using an explicit threshold. Retain original calls/probabilities. Report coverage = resolved masked sites / all masked sites and accuracy among resolved sites alongside unfiltered accuracy; zero resolved sites produces N/A, not perfect accuracy. Filtering is a display/decision policy and need not rerun the HMM.

Edge cases: zero masked sites gives N/A accuracy/Brier and a clear no-sites-to-impute message; no observations gives valid prior/panel-informed predictions with an explanation; monomorphic panel columns and exact ties must not fail; all probabilities must remain finite.

All takeaway text is deterministic and rule-based. No LLM calls, paid APIs, or external interpretation service. Summaries name what changed, state actual outcomes, note input differences and small sample counts, and avoid causal/general claims from a single example. If confidence rises while accuracy falls, say so. A panel with better target representation is not guaranteed to win. Population labels are not algorithmic bonuses; they determine sampled panel membership only.

## 8. History and exports

Persist experiments, input subexperiments, completed runs, and selection state locally in the browser, preferably IndexedDB. No accounts or server-side personal history. Explain that storage is device/browser-specific and may be cleared. Provide clear-history and per-experiment delete controls with appropriate confirmation for removal of stored work. New experiment does not delete old work.

Comparison tables default to the current input subexperiment and focus on panel composition/size changes. Permit explicit cross-input comparison with prominent difference labels. Include timestamps, IDs, actual counts, metrics, and model version. Store immutable snapshots rather than pointers to mutable current settings.

JSON exports must be versioned, sufficient to restore the experiment, and contain seed/PRNG identity, generator and model versions, parameters, exact masks, site catalog, target truth, reference candidate data or exact reproducible generation artifacts, selected panel IDs, numerical results, and grouping metadata. Storing realized synthetic arrays avoids dependence on a future library's random-number implementation. Include no credentials. Validate import schemas and bounds; reject incompatible versions with an explanatory message. Restoration must preserve original results and identify explicit recomputation separately.

CSV export offers a run-summary table and a site-level results table, downloadable separately or selected from the same export menu. Summary rows include experiment/input/run IDs, settings, denominators, both methods' metrics, and threshold if used. Site rows include IDs, synthetic position, nucleotide alleles, truth, observed/missing status, panel frequency, HMM probability/call/confidence/correctness, baseline equivalents, and unresolved status. Blank correctness for observed sites; do not count them in accuracy. Include version identifiers and explain columns in the UI. Use ordinary UTF-8 CSV readable by R and Python.

## 9. Architecture and service contract

Use React with TypeScript and Vite for the frontend. Use Python with FastAPI/Pydantic for the backend and NumPy for numeric arrays; pytest for numerical validation. Pin compatible dependencies when building. SVG/CSS and an existing lightweight animation approach are sufficient; avoid an elaborate graphics engine.

Suggested repository layout:

```text
frontend/                 React application and interaction tests
backend/app/api.py        FastAPI routes and validation
backend/app/simulation.py Pure synthetic-data generation
backend/app/masking.py    Mask operations and identity rules
backend/app/hmm.py        Pure inference and trace construction
backend/app/baseline.py   Independent-site comparator
backend/app/metrics.py    Evaluation with explicit truth input
backend/app/schemas.py    Versioned request/response contracts
backend/tests/            Numerical and integration tests
fixtures/                 Small independently checked examples
docs/                     Glossary, sources, deployment instructions
SPEC.md                   This implementation contract
README.md                 Setup, run, test, architecture, limitations
```

Statistics functions must run independently of FastAPI so the owner can inspect them in Python. Keep truth out of inference arguments and traces. Server orchestration may generate and evaluate truth, but never use it to select reference rows or fill masked sites.

Proposed API: GET /health; POST /experiments to generate bounded synthetic data; POST /impute for observations, mask, selected panel, and model settings; POST /evaluate for a separate evaluation operation. These may be orchestrated through one frontend action, but preserve the functional boundary. Responses identify experiment/input/run IDs, request ID, generator/model/schema versions, selected settings, probabilities, scaling factors, and bounded trace arrays. Separate evaluation can also be done in a small pure frontend routine if matched to Python fixtures; prefer Python as the canonical implementation.

Use stateless requests or self-contained experiment payloads; do not rely on free-host process memory for persistence. Limit L and K server-side. Avoid returning a transition tensor; return parameters and arrays sufficient to reconstruct cell contributions. Reject stale responses when users start another run. Offer clear loading, cold-start, retry, invalid-request, and offline states. Cache identical immutable requests only if the key covers all input arrays, parameters, and versions.

Configure allowed frontend origins, environment-specific API URL, sensible request timeouts, and bounded request sizes. No database, authentication, genetic-data uploads, or paid API is necessary.

## 10. Deployment and costs

The finished application must be publicly accessible without installation, a visitor account, or payment. Development and hosting must use free options; do not start a paid plan or require billing without a new decision from the owner.

Select a free static frontend host and a free Python-capable backend host during the first implementation week. Provider names and availability are deliberately not asserted here: verify current official terms, inactivity behavior, build/runtime limits, account/card requirements, and public-use eligibility before choosing. Record the actual provider decision and deploy commands in README. Free tiers may sleep or change; disclose loading delays in the app and maintain a local launch path. If no suitable zero-cost Python service is available, report the constraint early rather than silently changing the agreed architecture or adding charges.

Prove deployment with a minimal React-to-FastAPI round trip before building the full visual experience. Test from a browser without a developer login. A local backup and exported deterministic examples support the presentation but do not replace the public-site goal.

## 11. Scientific context and content

Include concise blurbs linked to observations users can make: representation, panel size, linkage, missingness, uncertainty, and limitations of one example. All fictional populations must be explicitly labeled synthetic; do not attach real population names to invented frequency models or imply quantitative applicability to them.

Primary real-world motivation supplied by the owner:

Carlson, J. C., Krishnan, M., Liu, S., et al. (2026). The Soifua Manuia reference panel with 2,570 Samoan haplotypes improves genotype imputation quality among Samoans. Communications Biology 9, 1191. https://doi.org/10.1038/s42003-026-10788-7

The supplied paper describes a panel combining 1,285 Samoan individuals with 1000 Genomes participants and comparisons with other panels. Explain this as evidence from that study, not a promise that specificity always improves every prediction. Cite the paper directly for any numerical claim; distinguish its variant quality and evaluation measures from this app's toy metrics.

Describe the owner's Mexico City reference-panel project as ongoing research asking whether improved regional representation improves imputation; no results are established here. Discuss admixture and Indigenous ancestry respectfully and precisely: identity, geography, and genetic ancestry are distinct, and useful local haplotype matches can vary across a chromosome. Do not assume a city, nationality, or community defines a homogeneous genetic population. Use further primary citations only after verifying them; never invent empirical findings.

The course materials on dynamic programming and HMMs are pedagogical context. The project proposal assignment supplies the checkpoint and submission format; the interview supplies the final deadline. The owner will separately handle the Google Doc, repository URL, Canvas submission, and video. Writing this spec does not create a repository or submit anything.

## 12. Correctness and acceptance checks

Numerical correctness is essential because the app teaches the arithmetic. Implement meaningful tests rather than relying on visual plausibility:

1. On a tiny K=2, L=3 fixture, independently enumerate all hidden paths and compare likelihood, state posteriors, and missing-allele probabilities against scaled forward-backward results within 1e-10.
2. Compare optimized O(LK) recurrences with a naive transition-matrix implementation on several small seeded inputs.
3. Verify transition rows and posterior distributions sum to one, missing emissions equal one, terminal beta equals one, and all arrays stay finite at allowed bounds.
4. Test K=1, all-observed, all-missing, ties, monomorphic columns, and a true allele absent from the panel.
5. Verify changing hidden truth while keeping inference inputs identical cannot change inference results. Changing representation cannot regenerate target or mask.
6. Check deterministic seeds, nested reference selection where possible, automatic mask nesting, region-prefix preservation, and versioned JSON round trips.
7. Independently calculate metric counts, Brier score, baseline frequencies, and threshold coverage on small fixtures, including empty denominators.
8. Assert displayed cell contributions and normalization reproduce the backend cell from unrounded values. Test forward, backward, and posterior boundary cases.
9. End-to-end: intro skip/replay, controls/run, changed-settings indicator, site inspection, keyboard navigation, stale-response handling, and correct grouping after mask/panel changes.
10. For final expansions: history survives reload, exports retain exact settings, imported records restore correctly, manual drag masks the intended range, and cross-input comparisons show differences.

Benchmark the maximum proposed region/panel before shipping. Target warm compute under two seconds for the bounded workload and smooth local inspector stepping; report actual measurements rather than asserting untested performance. Cold starts are measured separately. Virtualize large grids and avoid serializing unnecessary trace expansions.

Checkpoint acceptance: a classmate can open the public link, understand the introduction, change panel representation with input fixed, run both methods, identify predictions and errors, and inspect a worked forward/backward calculation. The owner can explain why dynamic programming is used and demonstrate an independent numerical check.

Final acceptance additionally covers each completed expansion, reproducible downloads, clear limitations, a reliable local launch, and no required paid service. Do not claim deferred features are implemented.

## 13. Schedule and risk management

Budget approximately 5–7 owner hours/week: roughly 20–28 hours to the checkpoint and 17–24 more afterward. These are planning limits, not guaranteed feature estimates. Favor one complete vertical slice over disconnected polished screens.

| Period | Main outcome |
| --- | --- |
| Oct 26–Nov 1 | Repository/setup, free-host feasibility, deployed API round trip, contracts, tiny numerical fixture |
| Nov 2–8 | Independent simulator, HMM, baseline, metrics, correctness checks |
| Nov 9–15 | Three tracks, reference grid, controls, coherent default example |
| Nov 16–23 | Worked HMM inspector, text introduction/About, checkpoint deployment and testing |
| Nov 24–30 | Persistent grouped history, JSON restoration/export, CSV export |
| Dec 1–7 | Manual masking, then longer mode as capacity allows |
| Dec 8–10 | Narrated animation if feasible; otherwise improve existing explanations; feature freeze |
| Dec 11–17 | Fixes, deployment rehearsal, owner video preparation and final submission |

If core calculations or worked-cell inspection slip, defer cosmetic animation and extensions first. If local storage grows, cap retained runs with an explicit notice and export option; never silently discard experiments. If simulation effects are weak, inspect LD and sampling design rather than adding an accuracy bonus by population label. If effects are unrealistically deterministic, increase overlap and variation and document the model's limits.

Longer regions can add informative evidence but can also introduce harder sites; neither accuracy nor confidence must increase monotonically. Default summaries should explain observed outcomes without imposing an expected ordering. Panel size, representation, and missingness comparisons must remain independently controllable.

## 14. Extensions, explicitly outside the required core

- Multiple independent targets and aggregate evaluation across repetitions, with uncertainty and valid per-variant metrics where defined.
- Admixed targets with local ancestry segments and appropriately composed references.
- Two homologous chromosome copies; phase/genotype modeling needs a separately specified statistical extension, not two decorative complementary strands.
- Verified public-data example, including provenance and licensing.
- Adjustable HMM mismatch/redraw parameters with explanations and recorded settings.
- Optional Viterbi overlay, clearly distinguished from marginal posterior inference.

Do not implement these before the agreed priorities without revisiting scope. The central deliverable remains an understandable, reproducible educational experiment with inspectable Python computation.
