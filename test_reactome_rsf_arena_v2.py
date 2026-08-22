from __future__ import annotations

import json
from pathlib import Path
import tempfile
import unittest

import numpy as np
import pandas as pd

from reactome_rsf_arena_v2 import integrity, prepare


VALID = {
    "name": "unit_module",
    "selector": "dr_gene",
    "n_genes": 8,
    "representation": "module",
    "module_count": 2,
    "benefit_threshold_months": 0.0,
    "rsf": {
        "n_estimators": 300,
        "max_depth": 6,
        "min_samples_leaf": 16,
        "min_samples_split": 32,
        "max_features": 0.5,
    },
}


class CandidateTests(unittest.TestCase):
    def test_valid_module_candidate(self):
        spec = prepare.validate_candidate(VALID)
        self.assertEqual(spec.module_count, 2)
        self.assertEqual(spec.representation, "module")

    def test_module_count_is_bounded(self):
        raw = {**VALID, "module_count": 5}
        with self.assertRaises(ValueError):
            prepare.validate_candidate(raw)

    def test_raw_representation_has_one_feature_per_gene(self):
        raw = {**VALID, "representation": "raw", "module_count": 8}
        self.assertEqual(prepare.validate_candidate(raw).module_count, 8)
        with self.assertRaises(ValueError):
            prepare.validate_candidate({**raw, "module_count": 2})

    def test_tree_budget_allows_1000_but_not_1001(self):
        accepted = {**VALID, "rsf": {**VALID["rsf"], "n_estimators": 1000}}
        self.assertEqual(prepare.validate_candidate(accepted).n_estimators, 1000)
        rejected = {**VALID, "rsf": {**VALID["rsf"], "n_estimators": 1001}}
        with self.assertRaises(ValueError):
            prepare.validate_candidate(rejected)


class ModuleTests(unittest.TestCase):
    def _frame(self) -> pd.DataFrame:
        n = 6
        values = {
            column: np.linspace(0.0, 1.0, n)
            for column in prepare.CLINICAL_COLUMNS
        }
        values[prepare.TREATMENT] = np.array([0, 1, 0, 1, 0, 1])
        values.update({
            "G1": np.arange(n, dtype=float),
            "G2": np.arange(n, dtype=float) + 2,
            "G3": np.arange(n, dtype=float) * 2,
            "G4": np.arange(n, dtype=float) * 3,
        })
        return pd.DataFrame(values)

    def test_module_transform_reduces_four_genes_to_two_features(self):
        frame = self._frame()
        transformer = prepare.FeatureTransformer(
            ["G1", "G2", "G3", "G4"], "module", 2
        ).fit(frame.iloc[:4])
        matrix = transformer.transform(frame.iloc[4:])
        self.assertEqual(matrix.shape, (2, len(prepare.CLINICAL_COLUMNS) + 2))
        self.assertEqual(transformer.names[-2:], ["REACTOME_MODULE_1", "REACTOME_MODULE_2"])

    def test_clinical_transform_contains_no_gene_features(self):
        frame = self._frame()
        transformer = prepare.FeatureTransformer([], "raw", 0).fit(frame)
        self.assertEqual(transformer.transform(frame).shape[1], len(prepare.CLINICAL_COLUMNS))


class MetricTests(unittest.TestCase):
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
        observed = prepare.recommendations(np.array([0.0, 1e-13]), threshold=0.0)
        np.testing.assert_array_equal(observed, np.zeros(2, dtype=np.int8))


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


if __name__ == "__main__":
    unittest.main()
