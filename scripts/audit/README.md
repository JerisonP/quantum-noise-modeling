# Independent NumPy checks

These scripts re-implement the old generators and estimators in NumPy, line
for line, to confirm the findings in `docs/AUDIT.md` without trusting the
Julia code. Requirements: `numpy`, `scipy`.

| Script | Confirms |
|--------|----------|
| `audit_acf_ou_gamma_wk.py` | A1 (1/f ACF theory mismatch), A2 (OU demeaning bias), B1 (Γ overflow → NaN), B2 (Wiener–Khinchin offset) |
| `audit_alpha_sweep.py` | A3 (α-sweep targets: exact finite-N expected periodogram) |
| `calibrate_statistical_tests.py` | the 5-SE statistical tests in `test/` are calibrated (z ≈ N(0,1) over 200 seeds) |

Run: `python3 scripts/audit/<name>.py` (the first two take about 1 minute and use about 1 GB of RAM).
