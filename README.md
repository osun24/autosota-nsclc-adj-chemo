# autosota-nsclc-adj-chemo
Testing AutoSOTA against an NSCLC adjuvant chemotherapy treatment recommendation problem

## AutoSOTA Arena Scaffold

This repository now has isolated arenas for future autonomous research sessions:

- `xgb_arena/prepare.py` is frozen utility code for train/validation loading, IPTW, labels, and validation metrics. It has no test-set path constant and guards against loading any path containing `test`.
- `xgb_arena/train.py` is the agent-editable XGBoost-Cox pipeline with a single `run() -> dict` entry point. The default loop budget is 10 Optuna trials x 2 bootstraps.
- `rsf_arena/` repeats the same structure for Random Survival Forest.
- `deepsurv_arena/` repeats the same structure for DeepSurv. It is GPU-friendly, but CPU-safe by default.
- `reactome_rsf_pfi_arena/` screens every development-available Reactome gene
  in an RSF, computes held-out permutation importance for five-year causal-RMST
  policy difference (with C-index secondary), and optimizes a final top-0..32
  gene panel plus forest hyperparameters.
- Each arena has `program.md`, `log.md`, and `runs/`.
- `red_lines.md` is the shared sealed-test and no-leakage rule set.
- `finalize.py` is HUMAN ONLY for XGBoost.
- `finalize-rsf.py`, `finalize-deepsurv.py`, and `finalize-all.py` are HUMAN ONLY finalization scripts for the other arenas and combined comparison.

To smoke-test the scaffold without touching the sealed test set:

```bash
python xgb_arena/train.py
python rsf_arena/train.py
python deepsurv_arena/train.py
```

To launch a later autonomous loop, agree on a fixed iteration budget first, then instruct the agent to follow one arena's `program.md`. The loop must compare candidates only on validation metrics and must not run any `finalize*.py` script.

After the autonomous phase is complete and `affyfRMATest.csv` has been added locally, perform a single-shot sealed evaluation:

```bash
python finalize.py --xgb-run-dir xgb_arena/runs/<chosen_run_dir>
python finalize-rsf.py --rsf-run-dir rsf_arena/runs/<chosen_run_dir>
python finalize-deepsurv.py --deepsurv-run-dir deepsurv_arena/runs/<chosen_run_dir>
python finalize-all.py \
  --xgb-run-dir xgb_arena/runs/<chosen_xgb_run> \
  --rsf-run-dir rsf_arena/runs/<chosen_rsf_run> \
  --deepsurv-run-dir deepsurv_arena/runs/<chosen_deepsurv_run>
```

By default, finalization refits the chosen configuration on Train+Validation before evaluating the sealed test set. Use `--no-refit-train-valid` to evaluate saved validation model artifacts directly.

## Colab GPU Bridge for DeepSurv

If Colab is the only available GPU source, keep the local autoresearch loop in
control and use Colab only as a manually started GPU worker. The bridge lives in
`colab_bridge/`.

Local side:

```bash
python colab_bridge/local_colab_queue.py submit \
  --queue-root "/path/to/Google Drive/autosota_colab_queue" \
  --job-id deepsurv_iter_001 \
  --require-pushed \
  --env DEEPSURV_ARENA_N_TRIALS=20 \
  --env DEEPSURV_ARENA_BOOTSTRAPS=2 \
  --env DEEPSURV_ARENA_EPOCHS=200
python colab_bridge/local_colab_queue.py watch \
  --queue-root "/path/to/Google Drive/autosota_colab_queue" \
  --job-id deepsurv_iter_001
python colab_bridge/local_colab_queue.py collect \
  --queue-root "/path/to/Google Drive/autosota_colab_queue" \
  --job-id deepsurv_iter_001
```

Colab side:

```python
from google.colab import drive
drive.mount("/content/drive")

!python /content/drive/MyDrive/path/to/colab_gpu_worker.py \
  --queue-root /content/drive/MyDrive/autosota_colab_queue \
  --data-dir /content/drive/MyDrive/nsclc_train_valid_only \
  --repo-url https://github.com/osun24/autosota-nsclc-adj-chemo.git
```

The Colab data directory must contain only `affyfRMATrain.csv`,
`affyfRMAValidation.csv`, and `LOOCV_Genes2.csv`. Do not put
`affyfRMATest.csv` there. Jobs are leased with heartbeat files; stale running
jobs are requeued after the lease expires, so a Colab disconnect can be picked
up by restarting the worker.
