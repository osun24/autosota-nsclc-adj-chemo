from __future__ import annotations

import json
import os
from pathlib import Path
import subprocess
import tempfile
import unittest
from unittest.mock import patch

import numpy as np
import pandas as pd

from rsf_tlearner_permimp_arena import integrity, prepare, run
from rsf_tlearner_permimp_arena.train import CANDIDATE


def _table() -> pd.DataFrame:
    return pd.DataFrame({
        "gene": ["OBS_ONLY", "ACT_ONLY", "NEGATIVE", "TIE_B", "TIE_A"],
        "observation_selection_score": [.30, 0, -.10, .05, .05],
        "act_selection_score": [0, .40, -.20, .02, .02],
        "observation_policy_pfi_mean": [.40, 0, -.05, .06, .06],
        "act_policy_pfi_mean": [0, .50, -.10, .03, .03],
        "joint_policy_pfi_mean": [.35, .45, -.08, .04, .04],
        "observation_cindex_pfi_mean": [.02, 0, 0, .01, .01],
        "act_cindex_pfi_mean": [0, .03, 0, .01, .01],
        "observation_positive_fraction": [1, 0, 0, .5, .5],
        "act_positive_fraction": [0, 1, 0, .5, .5],
    })


class CandidateContractTests(unittest.TestCase):
    def test_schema_and_independent_caps(self):
        spec = prepare.validate_candidate(CANDIDATE)
        self.assertEqual(spec.max_panel_genes, {"observation": 16, "act": 16})
        self.assertIsNone(spec.parent_run_id)

    def test_unknown_keys_and_cap_33_rejected(self):
        with self.assertRaises(ValueError):
            prepare.validate_candidate({**CANDIDATE, "genes": ["FORBIDDEN"]})
        bad = {**CANDIDATE, "max_panel_genes": {"observation": 33, "act": 0}}
        with self.assertRaises(ValueError):
            prepare.validate_candidate(bad)

    def test_empty_one_arm_is_valid(self):
        raw = {**CANDIDATE, "max_panel_genes": {"observation": 0, "act": 12}}
        self.assertEqual(prepare.validate_candidate(raw).max_panel_genes["observation"], 0)

    def test_act_is_absent_from_both_feature_matrices(self):
        n = 8
        frame = pd.DataFrame({column: np.arange(n, dtype=float) for column in prepare.PRETREATMENT_COLUMNS})
        frame[prepare.TREATMENT] = [0, 1] * 4
        frame["G"] = np.linspace(0, 1, n)
        for genes in ([], ["G"]):
            transformer = prepare.locked_t.FeatureTransformer(genes, "raw", len(genes)).fit(frame)
            self.assertNotIn(prepare.TREATMENT, transformer.names)


class ImportanceTests(unittest.TestCase):
    def test_obs_only_and_act_only_are_not_forced_shared(self):
        spec = prepare.validate_candidate({**CANDIDATE, "max_panel_genes": {"observation": 1, "act": 1}})
        panels, _ = prepare.select_panels(_table(), spec)
        self.assertEqual(panels["observation"], ["OBS_ONLY"])
        self.assertEqual(panels["act"], ["ACT_ONLY"])

    def test_unused_arm_is_exact_zero(self):
        self.assertEqual(prepare.component_drop(used=False, baseline=9.0, permuted=-4.0), 0.0)
        self.assertEqual(prepare.component_drop(used=True, baseline=9.0, permuted=7.5), 1.5)

    def test_negative_pfi_excluded_and_caps_apply(self):
        spec = prepare.validate_candidate({**CANDIDATE, "max_panel_genes": {"observation": 2, "act": 2}})
        panels, _ = prepare.select_panels(_table(), spec)
        self.assertNotIn("NEGATIVE", panels["observation"] + panels["act"])
        self.assertLessEqual(len(panels["observation"]), 2)
        self.assertLessEqual(len(panels["act"]), 2)

    def test_diverged_parent_pools_cannot_cross_migrate(self):
        table = _table().copy()
        table["observation_in_pool"] = [True, False, True, True, True]
        table["act_in_pool"] = [False, True, True, True, True]
        spec = prepare.validate_candidate({**CANDIDATE, "max_panel_genes": {"observation": 5, "act": 5}})
        panels, pools = prepare.select_panels(table, spec)
        self.assertNotIn("ACT_ONLY", panels["observation"] + pools["observation"])
        self.assertNotIn("OBS_ONLY", panels["act"] + pools["act"])

    def test_ranking_is_deterministic_with_symbol_tiebreak(self):
        first = prepare.rank_arm_importance(_table(), "observation")["gene"].tolist()
        second = prepare.rank_arm_importance(_table().sample(frac=1, random_state=7), "observation")["gene"].tolist()
        self.assertEqual(first, second)
        self.assertLess(first.index("TIE_A"), first.index("TIE_B"))

    def test_paired_permutations_are_identical_and_reproducible(self):
        a = prepare.paired_permutation_indices(50, context=2, gene_index=9, repeat=1)
        b = prepare.paired_permutation_indices(50, context=2, gene_index=9, repeat=1)
        np.testing.assert_array_equal(a, b)
        self.assertFalse(np.array_equal(a, prepare.paired_permutation_indices(50, context=2, gene_index=10, repeat=1)))

    def test_outer_assessment_outcomes_are_not_an_input_to_screening(self):
        fit = pd.DataFrame({"fit": [1, 2, 3]})
        outer_a = pd.DataFrame({"OS_STATUS": [0, 0]})
        outer_b = pd.DataFrame({"OS_STATUS": [1, 1]})
        with patch.object(prepare, "_importance_table", return_value=_table()) as mocked:
            left, _ = prepare.select_panels(prepare._importance_table(fit, {}, None, context=1, smoke=True), prepare.validate_candidate(CANDIDATE))
            right, _ = prepare.select_panels(prepare._importance_table(fit, {}, None, context=1, smoke=True), prepare.validate_candidate(CANDIDATE))
        self.assertEqual(left, right)
        self.assertFalse(outer_a.equals(outer_b))
        self.assertTrue(all(call.args[0] is fit for call in mocked.call_args_list))


class BoundaryTests(unittest.TestCase):
    def test_environment_and_manifests_verify(self):
        integrity.verify_lock()

    def test_ledger_integrity_and_cumulative_budget(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "ledger.jsonl"
            integrity.append_ledger(path, {"event": "finished", "experiment": 1, "wall_clock_seconds": 100})
            integrity.append_ledger(path, {"event": "finished", "experiment": 2, "status": "failed", "wall_clock_seconds": 200})
            rows = integrity.read_ledger(path)
            self.assertEqual(run.cumulative_elapsed(rows), 300)
            self.assertEqual(run.next_timeout(rows), 3600)
            altered = json.loads(path.read_text().splitlines()[0]); altered["experiment"] = 99
            lines = path.read_text().splitlines(); lines[0] = json.dumps(altered); path.write_text("\n".join(lines) + "\n")
            with self.assertRaises(RuntimeError):
                integrity.read_ledger(path)

    def test_unique_root_and_completed_parent(self):
        root = prepare.validate_candidate(CANDIDATE)
        run.validate_round_contract(1, root, None)
        with self.assertRaises(RuntimeError):
            run.validate_round_contract(2, root, None)
        child = prepare.validate_candidate({**CANDIDATE, "parent_run_id": "run_001_parent"})
        with tempfile.TemporaryDirectory() as directory:
            parent = Path(directory) / "run_001_parent"; parent.mkdir()
            with self.assertRaises(RuntimeError):
                run.validate_round_contract(2, child, parent)
            (parent / "result.json").write_text("{}")
            run.validate_round_contract(2, child, parent)

    def test_parent_pool_hashes(self):
        genes = ["G1", "G2"]
        with tempfile.TemporaryDirectory() as directory:
            parent = Path(directory) / "run_001_parent"; parent.mkdir()
            context = {"fit_index_sha256": "fit", "next_pools": {"observation": genes, "act": ["G2"]}, "next_pool_sha256": {"observation": prepare._pool_hash(genes), "act": prepare._pool_hash(["G2"])}}
            lineage = {"run_id": parent.name, "eligible_gene_universe_sha256": prepare._pool_hash(genes), "contexts": {"full_development": context}}
            path = parent / "pool_lineage.json"; path.write_text(json.dumps(lineage))
            (parent / "result.json").write_text(json.dumps({"artifacts": {"pool_lineage_sha256": integrity.sha256_file(path)}}))
            prepare.load_parent_pools(parent, parent.name, genes)
            lineage["contexts"]["full_development"]["next_pools"]["act"] = ["G1"]
            path.write_text(json.dumps(lineage))
            (parent / "result.json").write_text(json.dumps({"artifacts": {"pool_lineage_sha256": integrity.sha256_file(path)}}))
            with self.assertRaises(RuntimeError):
                prepare.load_parent_pools(parent, parent.name, genes)

    def test_nomination_is_strict_and_gate_dependent(self):
        threshold = float(prepare.BUDGET["run20_reward_to_beat"])
        self.assertFalse(run.qualifies_for_nomination({"eligible": True, "reward": threshold}))
        self.assertFalse(run.qualifies_for_nomination({"eligible": False, "reward": threshold + 1}))
        self.assertTrue(run.qualifies_for_nomination({"eligible": True, "reward": threshold + 1e-6}))


@unittest.skipUnless(os.environ.get("RUN_ARENA_SMOKE") == "1", "set RUN_ARENA_SMOKE=1 for the synthetic smoke integration")
class SmokeIntegrationTests(unittest.TestCase):
    def test_smoke_command_and_required_artifacts(self):
        completed = subprocess.run([os.sys.executable, "-m", "rsf_tlearner_permimp_arena.run", "--smoke"], cwd=Path(__file__).resolve().parent, check=False)
        self.assertEqual(completed.returncode, 0)


if __name__ == "__main__":
    unittest.main()
