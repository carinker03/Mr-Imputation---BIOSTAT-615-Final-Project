import { useState } from "react";

type Population = "A" | "B" | "C";
type OptionalAllele = 0 | 1 | null;

interface ExperimentSettings {
  seed: number;
  site_count: number;
  panel_size: number;
  target_population: Population;
  target_representation: number;
  observed_fraction: number;
}

interface ExperimentResponse {
  experiment_id: string;
  input_id: string;
  seed: number;
  mask_seed: number;
  target_population: Population;
  site_count: number;
  panel_size: number;
  target_representation: number;
  observed_fraction: number;
  target_truth: number[];
  observed_target: OptionalAllele[];
  observed_mask: boolean[];
  reference_haplotypes: number[][];
  reference_population_labels: Population[];
  reference_ids: string[];
  population_counts: Record<Population, number>;
  site_alleles: [string, string][];
}

interface ImputeResponse {
  allele_probabilities: [number, number][];
  allele_calls: number[];
}

const initialSettings: ExperimentSettings = {
  seed: 615,
  site_count: 40,
  panel_size: 40,
  target_population: "A",
  target_representation: 0.5,
  observed_fraction: 0.4,
};

const API_BASE_URL = (import.meta.env.VITE_API_BASE_URL ?? "/api").replace(
  /\/$/,
  "",
);

function alleleLabel(allele: number | null): string {
  return allele === null ? "missing" : String(allele);
}

function baseLabel(
  allele: number | null,
  site: number,
  siteAlleles: [string, string][],
): string {
  return allele === null ? "?" : siteAlleles[site][allele];
}

function App() {
  const [settings, setSettings] = useState(initialSettings);
  const [experiment, setExperiment] = useState<ExperimentResponse | null>(null);
  const [imputation, setImputation] = useState<ImputeResponse | null>(null);
  const [loading, setLoading] = useState<"generate" | "impute" | null>(null);
  const [error, setError] = useState<string | null>(null);

  function updateSetting<K extends keyof ExperimentSettings>(
    key: K,
    value: ExperimentSettings[K],
  ) {
    setSettings((current) => ({ ...current, [key]: value }));
  }

  async function generateExperiment() {
    setLoading("generate");
    setError(null);
    setImputation(null);
    try {
      const response = await fetch(`${API_BASE_URL}/experiments`, {
        method: "POST",
        headers: { "Content-Type": "application/json" },
        body: JSON.stringify(settings),
      });
      if (!response.ok) {
        throw new Error(await response.text());
      }
      setExperiment((await response.json()) as ExperimentResponse);
    } catch (requestError) {
      setError(
        requestError instanceof Error
          ? requestError.message
          : "Experiment generation failed.",
      );
      setExperiment(null);
    } finally {
      setLoading(null);
    }
  }

  async function runImputation() {
    if (!experiment) return;
    setLoading("impute");
    setError(null);
    try {
      const response = await fetch(`${API_BASE_URL}/impute`, {
        method: "POST",
        headers: { "Content-Type": "application/json" },
        body: JSON.stringify({
          observed_target: experiment.observed_target,
          reference_haplotypes: experiment.reference_haplotypes,
        }),
      });
      if (!response.ok) {
        throw new Error(await response.text());
      }
      setImputation((await response.json()) as ImputeResponse);
    } catch (requestError) {
      setError(
        requestError instanceof Error
          ? requestError.message
          : "Imputation failed.",
      );
    } finally {
      setLoading(null);
    }
  }

  return (
    <main className="app-shell">
      <header className="app-header">
        <p className="eyebrow">Synthetic statistical genetics lab</p>
        <h1>Mr. Imputation</h1>
        <p>
          Explore how synthetic reference panels and observed SNPs influence
          haploid imputation.
        </p>
        <p className="info-copy">
          This educational simulator creates a fictional target chromosome,
          hides some of its alleles, and uses a reference panel to make
          probabilistic predictions. Because this is synthetic data, the
          hidden truth is available for visual comparison and learning.
        </p>
      </header>

      <section className="settings-panel" aria-labelledby="settings-heading">
        <div>
          <p className="eyebrow">1. Create an experiment</p>
          <h2 id="settings-heading">Experiment settings</h2>
          <p className="info-copy">
            Choose the target population, region size, panel composition, and
            observation density. These controls prepare an experiment; the
            model runs only when you press the button, so slider changes do not
            silently create new results.
          </p>
        </div>
        <div className="settings-grid">
          <label>
            Target population
            <select
              value={settings.target_population}
              onChange={(event) =>
                updateSetting("target_population", event.target.value as Population)
              }
            >
              <option value="A">Synthetic population A</option>
              <option value="B">Synthetic population B</option>
              <option value="C">Synthetic population C</option>
            </select>
          </label>
          <label>
            SNP region length
            <input
              type="number"
              min="1"
              max="400"
              value={settings.site_count}
              onChange={(event) =>
                updateSetting("site_count", Number(event.target.value))
              }
            />
          </label>
          <label>
            Reference panel size
            <input
              type="number"
              min="1"
              max="300"
              value={settings.panel_size}
              onChange={(event) =>
                updateSetting("panel_size", Number(event.target.value))
              }
            />
          </label>
          <label>
            Target-population representation
            <input
              type="range"
              min="0"
              max="1"
              step="0.1"
              value={settings.target_representation}
              onChange={(event) =>
                updateSetting("target_representation", Number(event.target.value))
              }
            />
            <span>{Math.round(settings.target_representation * 100)}%</span>
          </label>
          <label>
            Observation density
            <input
              type="range"
              min="0"
              max="1"
              step="0.1"
              value={settings.observed_fraction}
              onChange={(event) =>
                updateSetting("observed_fraction", Number(event.target.value))
              }
            />
            <span>{Math.round(settings.observed_fraction * 100)}%</span>
          </label>
          <label>
            Random seed
            <input
              type="number"
              value={settings.seed}
              onChange={(event) => updateSetting("seed", Number(event.target.value))}
            />
          </label>
        </div>
        <button onClick={generateExperiment} disabled={loading !== null}>
          {loading === "generate" ? "Generating…" : "Generate Experiment"}
        </button>
        <p className="helper-text">
          Changing settings prepares a new request only. Nothing runs until you
          click the button.
        </p>
      </section>

      {error && <p className="error-banner" role="alert">{error}</p>}

      {!experiment && !loading && (
        <section className="empty-state">
          <h2>No experiment generated yet</h2>
          <p>Choose settings above to create a reproducible synthetic region.</p>
        </section>
      )}

      {experiment && (
        <>
          <section className="metadata-card">
            <div>
              <span className="eyebrow">Generated experiment</span>
              <strong>{experiment.experiment_id}</strong>
            </div>
            <div>
              <span className="eyebrow">Target</span>
              <strong>Population {experiment.target_population}</strong>
            </div>
            <div>
              <span className="eyebrow">Seeds</span>
              <strong>{experiment.seed} / mask {experiment.mask_seed}</strong>
            </div>
          </section>

          <section className="visual-panel">
            <div className="section-heading">
              <div>
                <p className="eyebrow">2. Inspect the input</p>
                <h2>Aligned synthetic SNP tracks</h2>
                <p className="info-copy">
                  Truth is the complete simulated chromosome. Observed is the
                  information given to the model. A missing tile means that no
                  allele was measured at that SNP; it does not mean DNA was
                  deleted.
                </p>
              </div>
              <span>{experiment.site_count} sites</span>
            </div>
            <Track
              label="Truth"
              values={experiment.target_truth}
              siteAlleles={experiment.site_alleles}
            />
            <Track
              label="Observed"
              values={experiment.observed_target}
              siteAlleles={experiment.site_alleles}
            />
            <div className="legend">
              <span><i className="allele-zero">0</i> first synthetic base / allele 0</span>
              <span><i className="allele-one">1</i> second synthetic base / allele 1</span>
              <span><i className="missing">?</i> missing</span>
            </div>
          </section>

          <section className="visual-panel">
            <div className="section-heading">
              <div>
                <p className="eyebrow">Reference evidence</p>
                <h2>Reference panel</h2>
                <p className="info-copy">
                  Each row is one synthetic haplotype the HMM may copy locally.
                  Population labels describe how rows were sampled; they are
                  not bonuses applied by the inference algorithm.
                </p>
              </div>
              <span>{experiment.panel_size} haplotypes</span>
            </div>
            <div className="reference-grid">
              {experiment.reference_haplotypes.map((haplotype, index) => (
                <div className="reference-row" key={experiment.reference_ids[index]}>
                  <span className="reference-label">
                    {experiment.reference_population_labels[index]} ·{" "}
                    {experiment.reference_ids[index]}
                  </span>
                  <div className="snp-row">
                    {haplotype.map((allele, site) => (
                      <span className={`snp-tile allele-${allele}`} key={site}>
                        <strong>{experiment.site_alleles[site][allele]}</strong>
                        <small>{allele}</small>
                      </span>
                    ))}
                  </div>
                </div>
              ))}
            </div>
          </section>

          <section className="visual-panel">
            <div className="section-heading">
              <div>
                <p className="eyebrow">3. Run the model</p>
                <h2>Imputation</h2>
                <p className="info-copy">
                  The forward-backward HMM combines observed alleles with
                  linkage patterns in the reference panel. It returns
                  probabilities for each allele, then chooses the more
                  probable allele at each site; exact ties resolve to allele 0.
                </p>
              </div>
              <button onClick={runImputation} disabled={loading !== null}>
                {loading === "impute" ? "Running…" : "Run Imputation"}
              </button>
            </div>
            {!imputation && <p className="empty-inline">Run the HMM to see posterior allele probabilities.</p>}
            {imputation && (
              <>
                <p className="result-note">
                  The imputed row preserves observed alleles and fills only
                  sites that were missing from the input.
                </p>
                <p className="interpretation-copy">
                  Expected dosage is the posterior probability of allele 1 in
                  this haploid model. Confidence is the larger of P(0) and
                  P(1), so it describes certainty in this prediction, not
                  empirical accuracy. A confident prediction can still be
                  wrong when the selected panel lacks useful support.
                </p>
                <Track
                  label="Truth"
                  values={experiment.target_truth}
                  siteAlleles={experiment.site_alleles}
                />
                <Track
                  label="Imputed"
                  values={experiment.observed_target.map(
                    (allele, site) => allele ?? imputation.allele_calls[site],
                  )}
                  siteAlleles={experiment.site_alleles}
                />
                <div className="probability-grid">
                  {imputation.allele_probabilities.map(([probabilityZero, probabilityOne], site) => {
                    const confidence = Math.max(probabilityZero, probabilityOne);
                    const wasObserved = experiment.observed_mask[site];
                    return (
                      <div className="probability-cell" key={site}>
                        <span className="site-label">SNP {site}</span>
                        <strong>
                          {wasObserved ? "Observed" : "Imputed"} · Call:{" "}
                          {experiment.site_alleles[site][imputation.allele_calls[site]]}{" "}
                          <small>({imputation.allele_calls[site]})</small>
                        </strong>
                        <span>Expected dosage: {probabilityOne.toFixed(3)}</span>
                        <span>Confidence: {confidence.toFixed(3)}</span>
                        <div className="probability-bar">
                          <span style={{ width: `${confidence * 100}%` }} />
                        </div>
                        <small>
                          P(0): {probabilityZero.toFixed(3)} · P(1):{" "}
                          {probabilityOne.toFixed(3)}
                        </small>
                        <small title="Standard population-level imputation r² requires multiple targets or repeated draws.">
                          r²: N/A for one synthetic target
                        </small>
                      </div>
                    );
                  })}
                </div>
              </>
            )}
          </section>
        </>
      )}
    </main>
  );
}

function Track({
  label,
  values,
  siteAlleles,
}: {
  label: string;
  values: Array<number | null>;
  siteAlleles: [string, string][];
}) {
  return (
    <div className="track">
      <strong className="track-label">{label}</strong>
      <div className="track-content">
        <div className="snp-row">
          {values.map((allele, site) => (
            <span
              className={`snp-tile ${allele === null ? "missing" : `allele-${allele}`}`}
              key={site}
              title={`SNP ${site}: ${alleleLabel(allele)}`}
            >
              <strong>{baseLabel(allele, site, siteAlleles)}</strong>
              <small>{allele === null ? "—" : allele}</small>
            </span>
          ))}
        </div>
        <div className="position-row" aria-label={`${label} SNP positions`}>
          {values.map((_, site) => (
            <span key={site}>{site}</span>
          ))}
        </div>
      </div>
    </div>
  );
}

export default App;
