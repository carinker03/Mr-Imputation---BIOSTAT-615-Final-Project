import unittest

import numpy as np

from backend.app.hmm import run_forward_backward
from backend.app.masking import (
    create_automatic_mask,
    create_manual_mask,
    update_manual_mask,
)


class MaskingTests(unittest.TestCase):
    def setUp(self):
        self.target = [0, 1, 0, 1, 1, 0, 0, 1, 0, 1]

    def test_automatic_masks_are_reproducible(self):
        first = create_automatic_mask(self.target, 0.4, seed=21)
        second = create_automatic_mask(self.target, 0.4, seed=21)
        self.assertEqual(first, second)

    def test_automatic_mask_count_uses_requested_density(self):
        result = create_automatic_mask(self.target, observed_fraction=0.4, seed=22)
        self.assertEqual(result.observed_count, 4)
        self.assertEqual(result.missing_count, 6)
        self.assertEqual(sum(result.observed_mask), 4)
        self.assertEqual(
            sum(value is None for value in result.hmm_target),
            6,
        )

    def test_density_boundaries_produce_all_missing_or_all_observed(self):
        all_missing = create_automatic_mask(self.target, 0.0, seed=23)
        all_observed = create_automatic_mask(self.target, 1.0, seed=23)
        self.assertEqual(all_missing.observed_mask, (False,) * len(self.target))
        self.assertEqual(all_missing.hmm_target, (None,) * len(self.target))
        self.assertEqual(all_observed.observed_mask, (True,) * len(self.target))
        self.assertEqual(all_observed.hmm_target, tuple(self.target))

    def test_manual_masking_and_unmasking_positions(self):
        initial = create_manual_mask(
            self.target,
            missing_positions=(1, 5),
        )
        updated = update_manual_mask(
            initial,
            observed_positions=(1,),
            missing_positions=(8,),
        )
        self.assertTrue(updated.observed_mask[1])
        self.assertFalse(updated.observed_mask[5])
        self.assertFalse(updated.observed_mask[8])
        self.assertEqual(updated.truth, tuple(self.target))
        self.assertEqual(updated.hmm_target[1], self.target[1])
        self.assertIsNone(updated.hmm_target[8])

    def test_truth_is_not_changed(self):
        original = list(self.target)
        result = create_automatic_mask(self.target, 0.5, seed=24)
        self.assertEqual(self.target, original)
        self.assertEqual(result.truth, tuple(original))

    def test_masked_target_is_compatible_with_hmm(self):
        result = create_automatic_mask(self.target, 0.5, seed=25)
        panel = [
            [0, 1, 0, 1, 1, 0, 0, 1, 0, 1],
            [1, 0, 1, 0, 0, 1, 1, 0, 1, 0],
        ]
        inference = run_forward_backward(result.hmm_target, panel)
        self.assertEqual(len(inference.gamma), len(self.target))
        self.assertTrue(
            all(
                allele is None or allele in (0, 1)
                for allele in result.hmm_target
            )
        )

    def test_automatic_density_changes_are_nested_for_same_seed(self):
        sparse = create_automatic_mask(self.target, 0.3, seed=26)
        dense = create_automatic_mask(self.target, 0.7, seed=26)
        for sparse_observed, dense_observed in zip(
            sparse.observed_mask, dense.observed_mask
        ):
            if sparse_observed:
                self.assertTrue(dense_observed)

    def test_invalid_targets_and_positions_raise_errors(self):
        invalid_targets = ([], [0, 2], [True, 0], ["0", 1])
        for target in invalid_targets:
            with self.subTest(target=target):
                with self.assertRaises(ValueError):
                    create_automatic_mask(target)
        with self.assertRaises(ValueError):
            create_automatic_mask(self.target, observed_fraction=-0.1)
        with self.assertRaises(ValueError):
            create_automatic_mask(self.target, observed_fraction=1.1)
        with self.assertRaises(ValueError):
            create_manual_mask(self.target, missing_positions=(len(self.target),))
        with self.assertRaises(ValueError):
            create_manual_mask(
                self.target,
                observed_positions=(2,),
                missing_positions=(2,),
            )

    def test_manual_all_missing_and_all_observed_masks(self):
        all_missing = create_manual_mask(
            self.target,
            missing_positions=range(len(self.target)),
        )
        all_observed = create_manual_mask(self.target)
        self.assertEqual(all_missing.missing_count, len(self.target))
        self.assertEqual(all_observed.observed_count, len(self.target))

    def test_automatic_mask_does_not_depend_on_numpy_global_state(self):
        np.random.seed(123)
        first = create_automatic_mask(self.target, 0.5, seed=27)
        np.random.seed(999)
        second = create_automatic_mask(self.target, 0.5, seed=27)
        self.assertEqual(first.observed_mask, second.observed_mask)


if __name__ == "__main__":
    unittest.main()
