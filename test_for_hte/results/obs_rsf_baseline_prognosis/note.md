The headline number is real but must not be read as printed
Tuning went fine (full-OBS C-index 0.702; held-out RMST stable at 45.0–47.6 months across the 5 outer folds). The table in baseline_prognosis_table.csv:

Aligned	Discordant
N	672	362
Mean predicted OBS RMST	48.9 mo	40.6 mo
SMD	+1.05	
An SMD of 1.05 is a massive prognostic imbalance — so the confounding check fires hard. But pooling the two treatment arms hides a sign reversal:

within actual-OBS patients (n=588/294): SMD +1.70 — aligned patients have better prognosis
within actual-ACT patients (n=84/68): SMD −1.04 — aligned patients have worse prognosis
The pooled +1.05 is a Simpson's-paradox composite of two opposite gradients weighted by a heavily unbalanced cohort (882 OBS vs 152 ACT). Since 85% of patients actually received OBS, "aligned" is mostly recommended-OBS, and the pooled column is close to a restatement of "recommendation correlates with prognosis."

Why that correlation exists
The T-learner's benefit estimate is strongly anti-correlated with baseline OBS prognosis: corr(benefit, rmst0) = −0.82. That is partly mechanical — benefit = rmst1 − rmst0 has rmst0 in it with a negative sign. But it isn't only mechanical: the independent nested baseline model, which never saw the policy model, reproduces it at corr = −0.76. So "recommend ACT" genuinely tracks "poor prognosis under observation" (mean rmst0 38.5 vs 50.6 months).

The practical consequence: aligned patients are healthier at baseline largely regardless of treatment, and the observed outcomes agree (within actual-OBS, mean OS 56.2 mo aligned vs 46.8 mo discordant). A naive aligned-vs-discordant survival comparison would be confounded to the point of being uninterpretable — which is exactly what this script was built to reveal.