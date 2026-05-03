# autosota-nsclc-adj-chemo
Testing AutoSOTA against an NSCLC adjuvant chemotherapy treatment recommendation problem

## AutoSOTA Arena Scaffold

This repository now has an isolated XGBoost arena for future autonomous research sessions:

- `xgb_arena/prepare.py` is frozen utility code for train/validation loading, IPTW, labels, and validation metrics. It has no test-set path constant and guards against loading any path containing `test`.
- `xgb_arena/train.py` is the agent-editable XGBoost-Cox pipeline with a single `run() -> dict` entry point. The default loop budget is 10 Optuna trials x 2 bootstraps.
- `xgb_arena/program.md` is the human-editable instruction file for the later autonomous loop.
- `xgb_arena/log.md` is append-only during autonomous iterations.
- `xgb_arena/runs/` stores model artifacts and metadata for completed runs.
- `red_lines.md` is the shared sealed-test and no-leakage rule set.
- `finalize.py` is HUMAN ONLY. Run it manually after adding `affyfRMATest.csv` locally.

To smoke-test the scaffold without touching the sealed test set:

```bash
python xgb_arena/train.py
```

To launch a later autonomous loop, agree on a fixed iteration budget first, then instruct the agent to follow `xgb_arena/program.md`. The loop must compare candidates only on validation metrics and must not run `finalize.py`.

After the autonomous phase is complete and `affyfRMATest.csv` has been added locally, perform a single-shot sealed evaluation:

```bash
python finalize.py --xgb-run-dir xgb_arena/runs/<chosen_run_dir>
```

By default, `finalize.py` refits the chosen XGB configuration on Train+Validation before evaluating the sealed test set. Use `--no-refit-train-valid` to evaluate the saved validation model artifact directly.
