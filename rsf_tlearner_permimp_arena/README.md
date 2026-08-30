# RSF T-Learner Permutation-Importance Arena

This isolated arena combines the locked Reactome RSF T-learner evaluation with progressive, fit-only permutation-importance selection. OBS and ACT receive independent raw-gene panels; a gene unused by one arm has exactly zero importance for that arm.

Run the reduced contract check with:

```bash
python -m rsf_tlearner_permimp_arena.run --smoke
```

Run one scientific attempt with:

```bash
python -m rsf_tlearner_permimp_arena.run
```

Only `train.py` is candidate-editable. Run 1 must be the unique root (`parent_run_id: null`). Every later run must name a completed parent and inherits that parent's same-fold next pools. The launcher limits the search to 20 attempts and 21,600 cumulative seconds. Failures and timeouts consume an attempt and elapsed time.

An eligible result must also have reward strictly above the frozen T-learner run-20 reward (`-3.5746550301982576`) before `test_nominee.txt` is populated. The sealed test set is never loaded here.
