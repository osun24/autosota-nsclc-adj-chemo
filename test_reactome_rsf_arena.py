from __future__ import annotations

import json
from pathlib import Path
import tempfile
import unittest

import numpy as np
from sksurv.functions import StepFunction

from reactome_rsf_arena import integrity, prepare


VALID = {
    "name": "unit",
    "selector": "dr_pathway",
    "n_genes": 12,
    "benefit_threshold_months": 0.25,
    "rsf": {
        "n_estimators": 300,
        "max_depth": 6,
        "min_samples_leaf": 16,
        "min_samples_split": 32,
        "max_features": 0.5,
    },
}


class CandidateTests(unittest.TestCase):
    def test_valid_candidate(self):
        spec = prepare.validate_candidate(VALID)
        self.assertEqual(spec.n_genes, 12)

    def test_explicit_gene_interaction_key_is_rejected(self):
        raw = {**VALID, "n_gene_interactions": 0}
        with self.assertRaises(ValueError):
            prepare.validate_candidate(raw)

    def test_genomic_features_are_main_effects_only(self):
        clinical, genomic = prepare.feature_names(["G1", "G2"])
        self.assertEqual(genomic, clinical + ["G1", "G2"])
        self.assertFalse(any("*ACT" in name for name in genomic))

    def test_feature_budget_is_enforced(self):
        raw = {**VALID, "n_genes": 33}
        with self.assertRaises(ValueError):
            prepare.validate_candidate(raw)

    def test_tree_budget_allows_1000_but_not_1001(self):
        accepted = {**VALID, "rsf": {**VALID["rsf"], "n_estimators": 1000}}
        self.assertEqual(prepare.validate_candidate(accepted).n_estimators, 1000)
        rejected = {**VALID, "rsf": {**VALID["rsf"], "n_estimators": 1001}}
        with self.assertRaises(ValueError):
            prepare.validate_candidate(rejected)

    def test_coordinated_split_leaf_is_enforced(self):
        raw = {**VALID, "rsf": {**VALID["rsf"], "min_samples_split": 20}}
        with self.assertRaises(ValueError):
            prepare.validate_candidate(raw)


class MetricTests(unittest.TestCase):
    def test_step_survival_is_integrated_on_rmst_scale(self):
        function = StepFunction(np.array([10.0, 20.0]), np.array([0.8, 0.5]))
        observed = prepare._integrate_step_functions([function], tau=30.0)[0]
        self.assertAlmostEqual(observed, 10.0 + 10.0 * 0.8 + 10.0 * 0.5)

    def test_aipw_arm_scores_reduce_to_outcome_regression_off_arm(self):
        phi0, phi1 = prepare.aipw_arm_scores(
            treatment=np.array([0, 1]),
            propensity=np.array([0.25, 0.25]),
            ipcw_time=np.array([10.0, 20.0]),
            mu0=np.array([12.0, 12.0]),
            mu1=np.array([15.0, 15.0]),
        )
        self.assertEqual(phi1[0], 15.0)
        self.assertEqual(phi0[1], 12.0)

    def test_identical_policies_have_zero_increment(self):
        phi0 = np.array([10.0, 20.0, 30.0, 40.0])
        phi1 = np.array([15.0, 18.0, 35.0, 38.0])
        rec = np.array([1, 0, 1, 0])
        result = prepare._bootstrap_increment(
            phi0, phi1, rec, phi0, phi1, rec, draws=50
        )
        self.assertEqual(result["mean"], 0.0)
        self.assertEqual(result["selection_lcb"], 0.0)

    def test_numerical_ties_never_recommend_act(self):
        rec = prepare.recommendations(np.array([0.0, 1e-13, -1e-13]), threshold=0.0)
        np.testing.assert_array_equal(rec, np.zeros(3, dtype=np.int8))

    def test_act_usage_summary_preserves_patient_distribution(self):
        summary = prepare._summarize_act_usage([
            {
                "tree_contains_act_split_fraction": 0.5,
                "patient_tree_path_act_fraction": np.array([0.0, 0.5]),
                "patient_terminal_difference_tree_fraction": np.array([0.0, 0.25]),
            },
            {
                "tree_contains_act_split_fraction": 1.0,
                "patient_tree_path_act_fraction": np.array([0.5, 1.0]),
                "patient_terminal_difference_tree_fraction": np.array([0.5, 1.0]),
            },
        ])
        self.assertEqual(summary["tree_contains_act_split_fraction"], 0.75)
        self.assertEqual(summary["patient_tree_paths_traversing_act_fraction"], 0.5)
        terminal = summary["per_patient_different_terminal_tree_fraction"]
        self.assertEqual(terminal["mean"], 0.4375)
        self.assertEqual(terminal["median"], 0.375)
        self.assertEqual(terminal["nonzero_patient_fraction"], 0.75)


class BoundaryTests(unittest.TestCase):
    def test_train_and_validation_paths_are_allowed(self):
        prepare._assert_allowed_search_path(prepare.TRAIN_CSV)
        prepare._assert_allowed_search_path(prepare.VALID_CSV)

    def test_test_and_arbitrary_paths_are_rejected(self):
        with self.assertRaises(ValueError):
            prepare._assert_allowed_search_path(prepare.REPO_ROOT / "affyfRMATest.csv")
        with self.assertRaises(ValueError):
            prepare._assert_allowed_search_path(prepare.REPO_ROOT / "other.csv")

    def test_ledger_detects_rewrite(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "ledger.jsonl"
            integrity.append_ledger(path, {"experiment": 1, "event": "started"})
            row = json.loads(path.read_text())
            row["event"] = "rewritten"
            path.write_text(json.dumps(row) + "\n")
            with self.assertRaises(RuntimeError):
                integrity.read_ledger(path)


if __name__ == "__main__":
    unittest.main()
