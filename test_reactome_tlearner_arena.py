from __future__ import annotations

import json
from pathlib import Path
import tempfile
import unittest

import numpy as np
import pandas as pd

from reactome_tlearner_arena import integrity, prepare


ARM = {
    "n_estimators": 300,
    "max_depth": 6,
    "min_samples_leaf": 8,
    "min_samples_split": 16,
    "max_features": 0.5,
}
VALID = {
    "name": "unit_tlearner",
    "selector": "dr_gene",
    "n_genes": 8,
    "representation": "module",
    "module_count": 2,
    "benefit_threshold_months": 0.0,
    "tlearner": {"observation": ARM, "act": {**ARM, "min_samples_leaf": 10, "min_samples_split": 20}},
}


class CandidateTests(unittest.TestCase):
    def test_valid_candidate_has_separate_arm_specs(self):
        spec = prepare.validate_candidate(VALID)
        self.assertEqual(spec.observation.min_samples_leaf, 8)
        self.assertEqual(spec.act.min_samples_leaf, 10)

    def test_both_arm_specs_are_required(self):
        with self.assertRaises(ValueError):
            prepare.validate_candidate({**VALID, "tlearner": {"act": ARM}})

    def test_tree_budget_rejects_1001(self):
        bad = {**VALID, "tlearner": {"observation": {**ARM, "n_estimators": 1001}, "act": ARM}}
        with self.assertRaises(ValueError):
            prepare.validate_candidate(bad)


class FeatureTests(unittest.TestCase):
    def _frame(self) -> pd.DataFrame:
        n = 8
        values = {column: np.linspace(0.0, 1.0, n) for column in prepare.CLINICAL_COLUMNS}
        values[prepare.TREATMENT] = np.array([0, 1] * 4)
        values.update({f"G{i}": np.arange(n, dtype=float) * i for i in range(1, 5)})
        return pd.DataFrame(values)

    def test_treatment_is_absent_from_features(self):
        transformer = prepare.FeatureTransformer([], "raw", 0).fit(self._frame())
        self.assertNotIn(prepare.TREATMENT, transformer.names)
        self.assertEqual(len(transformer.names), len(prepare.PRETREATMENT_COLUMNS))

    def test_module_transform_is_fit_only_and_compressed(self):
        frame = self._frame()
        transformer = prepare.FeatureTransformer(["G1", "G2", "G3", "G4"], "module", 2).fit(frame.iloc[:6])
        self.assertEqual(transformer.transform(frame.iloc[6:]).shape, (2, len(prepare.PRETREATMENT_COLUMNS) + 2))


class BoundaryTests(unittest.TestCase):
    def test_dependency_and_local_manifest_verify(self):
        integrity.verify_lock()

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

