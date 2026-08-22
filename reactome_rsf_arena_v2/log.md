# Reactome causal-RMST RSF v2 experiment log

Fixed budget: 20 full experiments. The starting candidate carries forward the
legitimate pathway-level stability result from v1; it has not consumed a v2
experiment. Smoke metrics are infrastructure-only.

Read `v1_findings.md` before proposing the first hypothesis. Every adaptive
probe, including train-only structural or covariate probes, must consume a full
experiment unless it is part of the locked smoke check.

Append each prespecified hypothesis before editing `train.py`, then complete
the same entry after the launcher returns using the schema in `program.md`.
