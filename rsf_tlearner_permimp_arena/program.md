# Program

1. The root run screens all development-eligible Reactome genes independently in both treatment arms.
2. Every outer fitting partition uses two inner folds. Shallow arm forests generate held-out OBS-component, ACT-component, joint-policy, and arm C-index permutation drops. Outer assessment rows do not determine panels.
3. An arm selection score is mean policy PFI minus one standard error. Only positive scores enter that arm's final raw-gene panel, capped at 32. The full deterministic ranking supplies its next pool.
4. Later candidates name a completed parent. Each outer fold and the full-development fit inherit only the corresponding parent pool after hashes and fitting-index lineage are verified.
5. The locked repeated two-by-four outer evaluator compares the genomic T-learner to the locked clinical T-learner with three seeds, IPCW-AIPW value, bootstrap selection LCB, repeat penalty, overlap, stability, use, and discrimination gates.
6. The human launcher records a hash-chained attempt, environment, elapsed time, code snapshot, log, results, PFI artifacts, panels, and lineage. Only eligible runs can lead; only a leader strictly better than frozen run 20 can become the test nominee.
