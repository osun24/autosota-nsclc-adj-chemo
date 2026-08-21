from __future__ import annotations

import tempfile
from pathlib import Path
import unittest

import numpy as np

from reactome_rsf.msigdb import parse_gmt
from reactome_rsf.run import (
    _alignment_rmst_difference,
    _assert_training_or_validation_path,
    select_available_reactome_genes,
)


class ReactomeRsfTests(unittest.TestCase):
    def test_parse_gmt_deduplicates_members_preserving_order(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            path = Path(temp_dir) / "tiny.gmt"
            path.write_text("SET_A\tdesc\tTP53\tEGFR\tTP53\nSET_B\tdesc\tKRAS\n")
            parsed = parse_gmt(path)
        self.assertEqual(parsed["SET_A"], ("TP53", "EGFR"))
        self.assertEqual(parsed["SET_B"], ("KRAS",))

    def test_gene_selection_requires_both_expression_splits(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            path = Path(temp_dir) / "tiny.gmt"
            path.write_text("SET_A\tdesc\tTP53\tEGFR\tKRAS\n")
            genes, total = select_available_reactome_genes(
                path,
                ["TP53", "EGFR", "train_only"],
                ["TP53", "KRAS", "valid_only"],
            )
        self.assertEqual(total, 3)
        self.assertEqual(genes, ["TP53"])

    def test_sealed_test_path_is_rejected(self):
        with self.assertRaisesRegex(ValueError, "Sealed test"):
            _assert_training_or_validation_path("affyfRMATest.csv")

    def test_delta_rmst_is_aligned_minus_not_aligned(self):
        result = _alignment_rmst_difference(
            time_values=np.array([50.0, 60.0, 10.0, 20.0]),
            event=np.array([1, 1, 1, 1], dtype=bool),
            aligned=np.array([True, True, False, False]),
            tau=60.0,
        )
        self.assertEqual(result["aligned_n"], 2)
        self.assertEqual(result["not_aligned_n"], 2)
        self.assertGreater(result["delta_rmst_months"], 0)


if __name__ == "__main__":
    unittest.main()
