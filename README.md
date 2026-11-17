# Type 2 diabetes core defects: gene ranking, Mendelian randomisation and disease model

Code for: S. Maharana, "Candidate genes for the three core defects of type 2 diabetes: tissue expression, Mendelian randomisation and an uncertainty-aware model of restoring function" (submitted to *Computers in Biology and Medicine*).

## Scripts (`scripts/`)

| Script | What it does | Paper section |
|---|---|---|
| `additive_prior_sensitivity.py` | Prior sensitivity of the additive model; binomial test of the 3-of-11 agreement | 2.7, 3.2, 3.3 |
| `additive_drug_check.py` | External check of the additive model against a thiazolidinedione and metformin | 2.7, 3.3 |
| `t2d_disposition_model.py` | Multiplicative (disposition-index) model fitted to baseline HbA1c and three drugs by MCMC; leave-one-drug-out check; restoration of each defect | 3.3 (revision) |
| `mr_coloc_checks.py` | MR sensitivity (F-statistic, IVW, MR-Egger, weighted median, Steiger) and colocalisation (coloc ABF) | 2.5, 3.2 |
| `workflow_figure.py` | Fig. 1 / graphical abstract | – |

## Run

```bash
pip install -r requirements.txt
python3 scripts/additive_prior_sensitivity.py
python3 scripts/additive_drug_check.py
python3 scripts/t2d_disposition_model.py          # writes drug_check.tsv, restoration.tsv, restoration.png
python3 scripts/mr_coloc_checks.py eqtl.tsv t2d.tsv instruments.tsv --cases 80154 --controls 853816 --genes genes.txt
python3 scripts/workflow_figure.py
```

Random seed: 20261018.

## Input formats for `mr_coloc_checks.py` (tab-separated)

- `eqtl.tsv`: gene, snp, beta_eqtl, se_eqtl, eaf, n_eqtl (InsPIRE nominal cis-eQTLs, effect allele aligned to T2D)
- `t2d.tsv`: snp, beta_t2d, se_t2d (DIAMANTE European T2D)
- `instruments.tsv`: gene, snp (independent lead eQTLs used in the main MR)
- `genes.txt` (optional): one gene per line (e.g. the 11 overlap genes and the Table 4 genes)

Check the DIAMANTE case and control counts for the release you use before running.

## Data sources (all public)

- NCBI GEO: GSE38642, GSE50397, GSE25724, GSE20966, GSE29221, GSE18732, GSE23343, GSE15653, GSE64998
- InsPIRE islet eQTLs: Zenodo record 3408356
- DIAMANTE T2D summary statistics: diagram-consortium.org

## Results (`results/`) and figures (`figures/`)

Outputs of `t2d_disposition_model.py` and the manuscript figures.

## Licence

MIT (see LICENSE).
