import itertools
import math
import unittest

from backend.app.hmm import call_alleles, run_forward_backward


def naive_forward_backward(target, panel, epsilon, redraw_probability):
    """Run the same HMM with an explicit transition matrix.

    This intentionally uses O(LK^2) nested summations and exists only as an
    independent numerical reference for the optimized implementation.
    """
    panel_size = len(panel)
    site_count = len(target)
    transition = [
        [
            (1.0 - redraw_probability if previous == current else 0.0)
            + redraw_probability / panel_size
            for current in range(panel_size)
        ]
        for previous in range(panel_size)
    ]

    def emission(site, state):
        if target[site] is None:
            return 1.0
        return (
            1.0 - epsilon
            if target[site] == panel[state][site]
            else epsilon
        )

    alpha = [[0.0] * panel_size for _ in range(site_count)]
    scaling_factors = [0.0] * site_count
    for state in range(panel_size):
        alpha[0][state] = emission(0, state) / panel_size
    scaling_factors[0] = sum(alpha[0])
    alpha[0] = [value / scaling_factors[0] for value in alpha[0]]

    for site in range(1, site_count):
        for current in range(panel_size):
            alpha[site][current] = emission(site, current) * sum(
                alpha[site - 1][previous] * transition[previous][current]
                for previous in range(panel_size)
            )
        scaling_factors[site] = sum(alpha[site])
        alpha[site] = [
            value / scaling_factors[site] for value in alpha[site]
        ]

    beta = [[0.0] * panel_size for _ in range(site_count)]
    beta[-1] = [1.0] * panel_size
    for site in range(site_count - 2, -1, -1):
        for previous in range(panel_size):
            beta[site][previous] = sum(
                transition[previous][current]
                * emission(site + 1, current)
                * beta[site + 1][current]
                for current in range(panel_size)
            ) / scaling_factors[site + 1]

    gamma = []
    for site in range(site_count):
        unnormalized = [
            alpha[site][state] * beta[site][state]
            for state in range(panel_size)
        ]
        total = sum(unnormalized)
        gamma.append([value / total for value in unnormalized])

    allele_probabilities = []
    for site in range(site_count):
        probability_one = sum(
            gamma[site][state]
            for state in range(panel_size)
            if panel[state][site] == 1
        )
        allele_probabilities.append([1.0 - probability_one, probability_one])

    return {
        "alpha": alpha,
        "beta": beta,
        "gamma": gamma,
        "scaling_factors": scaling_factors,
        "allele_probabilities": allele_probabilities,
    }


def enumerate_posteriors(target, panel, epsilon, redraw_probability):
    """Calculate the fixture result by enumerating every hidden state path."""
    panel_size = len(panel)
    site_count = len(target)
    path_weights = []
    for path in itertools.product(range(panel_size), repeat=site_count):
        weight = 1.0 / panel_size
        for site in range(site_count):
            if site:
                weight *= (
                    (1.0 - redraw_probability)
                    if path[site - 1] == path[site]
                    else 0.0
                ) + redraw_probability / panel_size
            if target[site] is None:
                emission = 1.0
            elif target[site] == panel[path[site]][site]:
                emission = 1.0 - epsilon
            else:
                emission = epsilon
            weight *= emission
        path_weights.append((path, weight))

    likelihood = sum(weight for _, weight in path_weights)
    posterior = [[0.0] * panel_size for _ in range(site_count)]
    for path, weight in path_weights:
        for site, state in enumerate(path):
            posterior[site][state] += weight / likelihood
    return posterior


class ForwardBackwardTests(unittest.TestCase):
    def setUp(self):
        self.target = [0, None, 1]
        self.panel = [[0, 0, 1], [1, 1, 0]]
        self.epsilon = 0.1
        self.redraw_probability = 0.2

    def test_matches_independent_path_enumeration(self):
        result = run_forward_backward(
            self.target, self.panel, self.epsilon, self.redraw_probability
        )
        expected = enumerate_posteriors(
            self.target, self.panel, self.epsilon, self.redraw_probability
        )
        for actual_site, expected_site in zip(result.gamma, expected):
            for actual, expected_value in zip(actual_site, expected_site):
                self.assertAlmostEqual(actual, expected_value, places=10)

    def test_optimized_result_matches_naive_matrix_reference(self):
        fixtures = [
            (
                [0, None, 1],
                [[0, 0, 1], [1, 1, 0], [0, 1, 1]],
                0.1,
                0.2,
            ),
            ([None, 1, None], [[0, 1, 0]], 0.25, 0.0),
            (
                [1, None, 0, 1],
                [[0, 0, 1, 1], [1, 1, 0, 0]],
                0.001,
                0.75,
            ),
        ]
        attributes = (
            "alpha",
            "beta",
            "gamma",
            "scaling_factors",
            "allele_probabilities",
        )
        for target, panel, epsilon, redraw_probability in fixtures:
            with self.subTest(
                target=target,
                panel=panel,
                epsilon=epsilon,
                redraw_probability=redraw_probability,
            ):
                optimized = run_forward_backward(
                    target, panel, epsilon, redraw_probability
                )
                naive = naive_forward_backward(
                    target, panel, epsilon, redraw_probability
                )
                for attribute in attributes:
                    optimized_values = getattr(optimized, attribute)
                    naive_values = naive[attribute]
                    if attribute == "scaling_factors":
                        optimized_values = [optimized_values]
                        naive_values = [naive_values]
                    self.assertEqual(len(optimized_values), len(naive_values))
                    for optimized_row, naive_row in zip(
                        optimized_values, naive_values
                    ):
                        self.assertEqual(len(optimized_row), len(naive_row))
                        for optimized_value, naive_value in zip(
                            optimized_row, naive_row
                        ):
                            self.assertAlmostEqual(
                                optimized_value,
                                naive_value,
                                delta=1e-10,
                            )

    def test_distributions_and_scaling_are_valid(self):
        result = run_forward_backward(self.target, self.panel)
        for distribution in result.alpha + result.gamma:
            self.assertAlmostEqual(sum(distribution), 1.0)
        for distribution in result.allele_probabilities:
            self.assertAlmostEqual(sum(distribution), 1.0)
        self.assertEqual(result.beta[-1], [1.0, 1.0])
        self.assertTrue(all(math.isfinite(value) for value in result.scaling_factors))

    def test_transition_rows_sum_to_one_at_boundary_redraw_values(self):
        for panel_size in (1, 2, 5):
            for redraw_probability in (0.0, 0.25, 1.0):
                transition = [
                    [
                        (1.0 - redraw_probability if previous == current else 0.0)
                        + redraw_probability / panel_size
                        for current in range(panel_size)
                    ]
                    for previous in range(panel_size)
                ]
                with self.subTest(
                    panel_size=panel_size,
                    redraw_probability=redraw_probability,
                ):
                    for row in transition:
                        self.assertAlmostEqual(sum(row), 1.0, delta=1e-10)

    def test_posterior_state_and_allele_distributions_sum_to_one(self):
        result = run_forward_backward(
            [0, None, 1, None],
            [[0, 0, 1, 1], [1, 1, 0, 0], [0, 1, 1, 0]],
            epsilon=0.2,
            redraw_probability=0.6,
        )
        for site in range(len(result.gamma)):
            self.assertAlmostEqual(sum(result.gamma[site]), 1.0, delta=1e-10)
            self.assertAlmostEqual(
                sum(result.allele_probabilities[site]),
                1.0,
                delta=1e-10,
            )

    def test_all_observed_target_is_supported(self):
        target = [0, 1, 1, 0]
        result = run_forward_backward(
            target,
            [[0, 1, 1, 0], [1, 0, 0, 1]],
            epsilon=0.01,
            redraw_probability=0.03,
        )
        self.assertEqual(len(result.gamma), len(target))
        for site in range(len(target)):
            self.assertAlmostEqual(sum(result.gamma[site]), 1.0, delta=1e-10)
            self.assertAlmostEqual(
                sum(result.allele_probabilities[site]),
                1.0,
                delta=1e-10,
            )
            observed_probability = result.allele_probabilities[site][target[site]]
            other_probability = result.allele_probabilities[site][1 - target[site]]
            self.assertGreater(observed_probability, other_probability)

    def test_all_missing_target_is_supported(self):
        result = run_forward_backward(
            [None, None, None],
            [[0, 0, 1], [1, 1, 0], [1, 0, 1]],
            epsilon=0.01,
            redraw_probability=0.03,
        )
        for site in range(3):
            self.assertAlmostEqual(sum(result.gamma[site]), 1.0)
            self.assertAlmostEqual(
                sum(result.allele_probabilities[site]),
                1.0,
            )
        self.assertEqual(result.gamma, [[1.0 / 3.0] * 3] * 3)

    def test_exact_posterior_tie_calls_allele_zero(self):
        result = run_forward_backward(
            [None],
            [[0], [1]],
            epsilon=0.01,
            redraw_probability=0.03,
        )
        self.assertEqual(result.allele_probabilities, [[0.5, 0.5]])
        self.assertEqual(call_alleles(result.allele_probabilities), [0])

    def test_missing_site_has_no_emission_preference(self):
        result = run_forward_backward(
            [None], [[0], [1]], epsilon=0.01, redraw_probability=0.03
        )
        self.assertEqual(result.alpha[0], [0.5, 0.5])
        self.assertEqual(result.allele_probabilities[0], [0.5, 0.5])

    def test_k_one_and_monomorphic_panel_are_supported(self):
        result = run_forward_backward([0, None, 1], [[0, 0, 1]])
        self.assertEqual(result.gamma, [[1.0], [1.0], [1.0]])
        self.assertEqual(
            result.allele_probabilities,
            [[1.0, 0.0], [1.0, 0.0], [0.0, 1.0]],
        )

    def test_rejects_mismatched_panel_lengths(self):
        with self.assertRaises(ValueError):
            run_forward_backward([0, 1], [[0]])


if __name__ == "__main__":
    unittest.main()
