import unittest
from unittest.mock import Mock

import numpy as np

from backend.app.hmm import run_forward_backward
from backend.app.masking import create_automatic_mask
from backend.app.simulation import simulate_experiment


class StatisticalPipelineIntegrationTests(unittest.TestCase):
    def _run_pipeline(self, experiment, masked_target):
        return run_forward_backward(
            masked_target.hmm_target,
            experiment.reference_haplotypes.tolist(),
            epsilon=0.01,
            redraw_probability=0.03,
        )

    def test_simulation_masking_and_hmm_outputs_are_compatible(self):
        experiment = simulate_experiment(
            seed=701,
            site_count=24,
            panel_size=12,
            candidates_per_population=20,
        )
        masked_target = create_automatic_mask(
            experiment.target_haplotype.tolist(),
            observed_fraction=0.5,
            seed=702,
        )
        result = self._run_pipeline(experiment, masked_target)

        site_count = 24
        panel_size = 12
        self.assertEqual(len(result.alpha), site_count)
        self.assertEqual(len(result.beta), site_count)
        self.assertEqual(len(result.gamma), site_count)
        self.assertEqual(len(result.allele_probabilities), site_count)
        for site in range(site_count):
            self.assertEqual(len(result.alpha[site]), panel_size)
            self.assertEqual(len(result.beta[site]), panel_size)
            self.assertEqual(len(result.gamma[site]), panel_size)
            self.assertEqual(len(result.allele_probabilities[site]), 2)

        for distribution in result.alpha + result.gamma:
            self.assertTrue(all(0.0 <= value <= 1.0 for value in distribution))
            self.assertAlmostEqual(sum(distribution), 1.0, delta=1e-10)
        for distribution in result.allele_probabilities:
            self.assertTrue(all(0.0 <= value <= 1.0 for value in distribution))
            self.assertAlmostEqual(sum(distribution), 1.0, delta=1e-10)
        self.assertTrue(
            all(
                np.isfinite(value)
                for row in result.beta
                for value in row
            )
        )

    def test_masking_and_inference_preserve_truth_and_reference_panel(self):
        experiment = simulate_experiment(
            seed=703,
            site_count=24,
            panel_size=12,
            candidates_per_population=20,
        )
        original_target = experiment.target_haplotype.copy()
        original_panel = experiment.reference_haplotypes.copy()
        masked_target = create_automatic_mask(
            experiment.target_haplotype.tolist(),
            observed_fraction=0.5,
            seed=704,
        )
        self._run_pipeline(experiment, masked_target)

        np.testing.assert_array_equal(experiment.target_haplotype, original_target)
        np.testing.assert_array_equal(
            experiment.reference_haplotypes, original_panel
        )
        self.assertEqual(masked_target.truth, tuple(original_target.tolist()))

    def test_repeating_pipeline_with_identical_seeds_is_identical(self):
        first_experiment = simulate_experiment(
            seed=705,
            site_count=24,
            panel_size=12,
            candidates_per_population=20,
        )
        second_experiment = simulate_experiment(
            seed=705,
            site_count=24,
            panel_size=12,
            candidates_per_population=20,
        )
        first_mask = create_automatic_mask(
            first_experiment.target_haplotype.tolist(),
            observed_fraction=0.5,
            seed=706,
        )
        second_mask = create_automatic_mask(
            second_experiment.target_haplotype.tolist(),
            observed_fraction=0.5,
            seed=706,
        )
        first_result = self._run_pipeline(first_experiment, first_mask)
        second_result = self._run_pipeline(second_experiment, second_mask)

        self.assertEqual(first_mask, second_mask)
        self.assertEqual(first_result, second_result)

    def test_panel_composition_change_does_not_change_target_or_mask(self):
        target_representation_heavy = simulate_experiment(
            seed=707,
            site_count=24,
            panel_size=12,
            target_representation=0.8,
            candidates_per_population=20,
        )
        other_representation_heavy = simulate_experiment(
            seed=707,
            site_count=24,
            panel_size=12,
            target_representation=0.2,
            candidates_per_population=20,
        )
        first_mask = create_automatic_mask(
            target_representation_heavy.target_haplotype.tolist(),
            observed_fraction=0.5,
            seed=708,
        )
        second_mask = create_automatic_mask(
            other_representation_heavy.target_haplotype.tolist(),
            observed_fraction=0.5,
            seed=708,
        )

        np.testing.assert_array_equal(
            target_representation_heavy.target_haplotype,
            other_representation_heavy.target_haplotype,
        )
        self.assertEqual(first_mask.observed_mask, second_mask.observed_mask)
        self.assertNotEqual(
            target_representation_heavy.reference_population_labels,
            other_representation_heavy.reference_population_labels,
        )

    def test_hmm_receives_observed_target_not_hidden_truth(self):
        experiment = simulate_experiment(
            seed=709,
            site_count=24,
            panel_size=12,
            candidates_per_population=20,
        )
        masked_target = create_automatic_mask(
            experiment.target_haplotype.tolist(),
            observed_fraction=0.5,
            seed=710,
        )
        observed_only_hmm = Mock(wraps=run_forward_backward)
        observed_only_hmm(
            masked_target.hmm_target,
            experiment.reference_haplotypes.tolist(),
            epsilon=0.01,
            redraw_probability=0.03,
        )

        self.assertEqual(
            observed_only_hmm.call_args.args[0],
            masked_target.hmm_target,
        )
        self.assertEqual(len(observed_only_hmm.call_args.args), 2)
        self.assertNotIn("truth", observed_only_hmm.call_args.kwargs)
        self.assertIn(None, observed_only_hmm.call_args.args[0])


if __name__ == "__main__":
    unittest.main()
