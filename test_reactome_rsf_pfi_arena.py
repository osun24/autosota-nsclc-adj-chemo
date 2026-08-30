from __future__ import annotations

import json
from pathlib import Path
import tempfile
import unittest

import numpy as np
import pandas as pd
from sksurv.functions import StepFunction

from reactome_rsf_pfi_arena import integrity, prepare, run


SPACE = {
    "n_estimators": [100, 200],
    "max_depth": [3, 5],
    "min_samples_leaf": [5, 10],
    "split_leaf_multiplier": [2, 3],
    "max_features": ["sqrt", 0.5],
}
VALID = {
    "name": "unit",
    "screening_space": SPACE,
    "reduced_space": SPACE,
    "plot_top_k": 10,
}


class CandidateTests(unittest.TestCase):
    def test_valid_candidate(self):
        spec = prepare.validate_candidate(VALID)
        self.assertEqual(spec.name, "unit")
        self.assertEqual(spec.plot_top_k, 10)

    def test_unrecognized_candidate_key_is_rejected(self):
        raw = {**VALID, "genes": ["ORACLE"]}
        with self.assertRaises(ValueError):
            prepare.validate_candidate(raw)

    def test_tree_budget_is_enforced(self):
        bad_space = {**SPACE, "n_estimators": [1001]}
        with self.assertRaises(ValueError):
            prepare.validate_candidate({**VALID, "screening_space": bad_space})


class TransformAndMetricTests(unittest.TestCase):
    def test_median_transform_is_fit_only(self):
        fit = pd.DataFrame({"A": [1.0, np.nan, 3.0], "B": [2.0, 4.0, 6.0]})
        assess = pd.DataFrame({"A": [np.nan, 1000.0], "B": [np.nan, 8.0]})
        transformer = prepare.MedianTransformer(["A", "B"]).fit(fit)
        observed = transformer.transform(assess)
        np.testing.assert_allclose(observed[0], [2.0, 4.0])

    def test_step_survival_is_integrated_on_rmst_scale(self):
        function = StepFunction(np.array([10.0, 20.0]), np.array([0.8, 0.5]))
        observed = prepare._integrate_step_functions([function], tau=30.0)[0]
        self.assertAlmostEqual(observed, 10.0 + 10.0 * 0.8 + 10.0 * 0.5)

    def test_array_rmst_matches_step_integral(self):
        survival = np.array([[0.8, 0.5]])
        observed = prepare._rmst_from_survival_array(
            survival, np.array([10.0, 20.0]), tau=30.0
        )[0]
        self.assertAlmostEqual(observed, 10.0 + 10.0 * 0.8 + 10.0 * 0.5)

    def test_numerical_ties_recommend_observation(self):
        observed = prepare.recommendations(
            np.array([10.0, 10.0]), np.array([10.0, 10.0 + 1e-13])
        )
        np.testing.assert_array_equal(observed, np.zeros(2, dtype=np.int8))

    def test_policy_difference_is_policy_minus_anti_policy(self):
        contribution = prepare._policy_contribution(
            np.array([10.0, 20.0]),
            np.array([15.0, 18.0]),
            np.array([1, 0]),
        )
        np.testing.assert_allclose(contribution, [5.0, 2.0])


class RankingAndPlotTests(unittest.TestCase):
    def _importance(self) -> pd.DataFrame:
        return pd.DataFrame({
            "feature": ["Age", "G2", "G1", "G3"],
            "feature_type": ["clinical", "gene", "gene", "gene"],
            "rmst_importance_mean_months": [0.2, 0.5, 0.5, -0.1],
            "rmst_importance_sd_months": [0.01, 0.02, 0.03, 0.01],
            "cindex_importance_mean": [0.01, 0.02, 0.03, 0.0],
            "cindex_importance_sd": [0.001, 0.002, 0.003, 0.001],
        })

    def test_gene_ranking_uses_secondary_cindex_tie_break(self):
        self.assertEqual(
            prepare.rank_genes_from_importance(self._importance()),
            ["G1", "G2", "G3"],
        )

    def test_plot_is_created(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "importance.png"
            prepare._plot_importance(self._importance(), path, top_k=4)
            self.assertGreater(path.stat().st_size, 0)


class BoundaryTests(unittest.TestCase):
    def test_only_two_development_paths_are_allowed(self):
        prepare._assert_allowed_search_path(prepare.TRAIN_CSV)
        prepare._assert_allowed_search_path(prepare.VALID_CSV)
        with self.assertRaises(ValueError):
            prepare._assert_allowed_search_path(prepare.REPO_ROOT / "unapproved.csv")

    def test_ledger_detects_rewrite(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "ledger.jsonl"
            integrity.append_ledger(path, {"experiment": 1, "event": "started"})
            row = json.loads(path.read_text())
            row["event"] = "rewritten"
            path.write_text(json.dumps(row) + "\n")
            with self.assertRaises(RuntimeError):
                integrity.read_ledger(path)

    def test_leader_comparison_is_lexicographic(self):
        prior = {"reward": 1.0, "secondary_score": 0.60}
        self.assertTrue(run._is_better(
            {"reward": 1.1, "secondary_score": 0.50}, prior
        ))
        self.assertTrue(run._is_better(
            {"reward": 1.0, "secondary_score": 0.61}, prior
        ))
        self.assertFalse(run._is_better(
            {"reward": 0.9, "secondary_score": 0.99}, prior
        ))


if __name__ == "__main__":
    unittest.main()

