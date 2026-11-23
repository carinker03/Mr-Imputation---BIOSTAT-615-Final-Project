import unittest

import numpy as np

from backend.app.simulation import (
    POPULATIONS,
    simulate_experiment,
)


class SimulationTests(unittest.TestCase):
    def test_identical_seed_and_parameters_are_reproducible(self):
        first = simulate_experiment(seed=11, site_count=40, panel_size=10)
        second = simulate_experiment(seed=11, site_count=40, panel_size=10)
        np.testing.assert_array_equal(
            first.target_haplotype, second.target_haplotype
        )
        np.testing.assert_array_equal(
            first.reference_haplotypes, second.reference_haplotypes
        )
        self.assertEqual(first.reference_ids, second.reference_ids)

    def test_different_seeds_can_produce_different_data(self):
        first = simulate_experiment(seed=11, site_count=40, panel_size=10)
        second = simulate_experiment(seed=12, site_count=40, panel_size=10)
        self.assertFalse(np.array_equal(first.target_haplotype, second.target_haplotype))

    def test_haplotypes_are_binary_and_have_requested_lengths(self):
        result = simulate_experiment(seed=13, site_count=37, panel_size=17)
        self.assertEqual(result.target_haplotype.shape, (37,))
        self.assertEqual(result.reference_haplotypes.shape, (17, 37))
        self.assertTrue(np.all(np.isin(result.target_haplotype, [0, 1])))
        self.assertTrue(np.all(np.isin(result.reference_haplotypes, [0, 1])))

    def test_reference_labels_match_rows_and_use_only_fictional_populations(self):
        result = simulate_experiment(seed=14, site_count=30, panel_size=20)
        self.assertEqual(len(result.reference_population_labels), 20)
        self.assertEqual(len(result.reference_ids), 20)
        self.assertTrue(
            set(result.reference_population_labels).issubset(set(POPULATIONS))
        )
        self.assertEqual(
            result.reference_haplotypes.shape[0],
            len(result.reference_population_labels),
        )

    def test_panel_size_and_composition_controls_are_respected(self):
        result = simulate_experiment(
            seed=15,
            site_count=30,
            panel_size=20,
            target_population="B",
            target_representation=0.5,
        )
        self.assertEqual(result.population_counts, {"A": 5, "B": 10, "C": 5})
        self.assertEqual(
            {
                population: result.reference_population_labels.count(population)
                for population in POPULATIONS
            },
            result.population_counts,
        )

    def test_target_is_stable_when_panel_settings_change(self):
        smaller = simulate_experiment(
            seed=16, site_count=40, panel_size=10, target_representation=0.1
        )
        larger = simulate_experiment(
            seed=16, site_count=40, panel_size=40, target_representation=0.9
        )
        np.testing.assert_array_equal(
            smaller.target_haplotype, larger.target_haplotype
        )
        np.testing.assert_array_equal(
            smaller.candidate_haplotypes, larger.candidate_haplotypes
        )

    def test_linkage_is_stronger_within_blocks_than_across_boundaries(self):
        result = simulate_experiment(
            seed=17,
            site_count=200,
            panel_size=40,
            block_size=20,
            mutation_rate=0.01,
        )
        within_block = []
        across_boundary = []
        for haplotype in result.candidate_haplotypes:
            for site in range(result.site_count - 1):
                agreement = haplotype[site] == haplotype[site + 1]
                if (site + 1) % result.block_size == 0:
                    across_boundary.append(agreement)
                else:
                    within_block.append(agreement)
        self.assertGreater(np.mean(within_block), np.mean(across_boundary))
        self.assertGreater(np.mean(within_block), 0.8)

    def test_population_weights_produce_deterministic_frequency_differences(self):
        result = simulate_experiment(
            seed=18,
            site_count=200,
            panel_size=60,
            candidates_per_population=100,
        )
        frequencies = {}
        for population in POPULATIONS:
            rows = [
                index
                for index, label in enumerate(result.candidate_population_labels)
                if label == population
            ]
            frequencies[population] = result.candidate_haplotypes[rows].mean()
        self.assertGreater(max(frequencies.values()) - min(frequencies.values()), 0.05)

    def test_invalid_configuration_values_raise_clear_errors(self):
        invalid_cases = (
            {"site_count": 0},
            {"panel_size": 0},
            {"target_population": "D"},
            {"target_representation": 1.1},
            {"candidates_per_population": 0},
            {"block_size": 0},
            {"mutation_rate": -0.1},
            {"panel_size": 301, "candidates_per_population": 100},
            {
                "panel_size": 300,
                "target_representation": 1.0,
                "candidates_per_population": 100,
            },
        )
        for parameters in invalid_cases:
            with self.subTest(parameters=parameters):
                with self.assertRaises(ValueError):
                    simulate_experiment(**parameters)

    def test_site_catalog_has_requested_length_and_varied_allele_pairs(self):
        result = simulate_experiment(seed=19, site_count=25, panel_size=10)
        self.assertEqual(len(result.site_alleles), 25)
        self.assertGreater(len(set(result.site_alleles)), 1)
        for allele_zero, allele_one in result.site_alleles:
            self.assertNotEqual(allele_zero, allele_one)


if __name__ == "__main__":
    unittest.main()
