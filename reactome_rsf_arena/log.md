# Reactome causal-RMST RSF experiment log

Fixed budget: 20 full experiments. Baseline is defined in `train.py` but has not
yet consumed a full experiment. Smoke-test metrics are not scientific feedback.

Append each prespecified hypothesis before editing `train.py`, then complete the
same entry after the launcher returns.

Every completed entry must include the mandatory `train_act_usage_cg` and
`validation_act_usage_cg` lines from `program.md`, even if the run is
ineligible. Each line reports the ACT split-tree fraction, counterfactual
patient-tree path-traversal fraction, and the mean/median/p10/p90/nonzero
summary of each patient's fraction of trees reaching different ACT-versus-OBS
terminal nodes.
