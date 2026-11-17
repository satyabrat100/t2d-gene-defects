"""
MR sensitivity and colocalisation checks for the islet MR (Section 2.5 / 3.2).

Inputs (tab-separated, one row per SNP):
  eqtl.tsv : gene, snp, beta_eqtl, se_eqtl, eaf, n_eqtl          (InsPIRE nominal cis-eQTL; same allele as T2D)
  t2d.tsv  : snp, beta_t2d, se_t2d                               (DIAMANTE, aligned to the eQTL effect allele)
  instruments.tsv : gene, snp                                    (independent lead eQTLs used in the main MR)

Outputs:
  mr_sensitivity.tsv : gene, n_iv, F_min, IVW, IVW_p, Egger, Egger_p, Egger_intercept_p,
                       WMedian, WMedian_p, Steiger_ok
  coloc.tsv          : gene, nsnps, PP.H0..PP.H4  (coloc ABF, Giambartolomei et al. 2014 default priors)

Usage:
  python3 mr_coloc_checks.py eqtl.tsv t2d.tsv instruments.tsv --cases 80154 --controls 853816 --genes genes.txt
  (genes.txt: optional list, e.g. the 11 overlap genes and the 6 genes of Table 4)
"""
import argparse
import numpy as np
import pandas as pd
from scipy import stats

def ivw(bx, by, sy):
    r = by / bx; se = sy / np.abs(bx); w = 1 / se**2
    b = np.sum(w * r) / np.sum(w); s = np.sqrt(1 / np.sum(w))
    if len(r) > 1:  # multiplicative random effects
        phi = max(1.0, np.sum(w * (r - b) ** 2) / (len(r) - 1)); s *= np.sqrt(phi)
    return b, s

def egger(bx, by, sy):
    s = np.sign(bx); bx, by = bx * s, by * s
    X = np.column_stack([np.ones_like(bx), bx]); W = np.diag(1 / sy**2)
    XtWX = X.T @ W @ X; coef = np.linalg.solve(XtWX, X.T @ W @ by)
    res = by - X @ coef; dof = len(by) - 2
    sigma2 = max(1.0, (res**2 / sy**2).sum() / dof)
    cov = np.linalg.inv(XtWX) * sigma2; se = np.sqrt(np.diag(cov))
    p = 2 * stats.t.sf(np.abs(coef / se), dof)
    return coef[1], se[1], p[1], p[0]

def weighted_median(bx, by, sy, nboot=1000, seed=1):
    def wm(r, w):
        o = np.argsort(r); r, w = r[o], w[o] / w.sum()
        c = np.cumsum(w) - 0.5 * w; return np.interp(0.5, c, r)
    r = by / bx; w = (bx / sy) ** 2; est = wm(r, w)
    rng = np.random.default_rng(seed)  # resamples T2D effects; eQTL error ignored (second order)
    boot = [wm(rng.normal(by, sy) / bx, w) for _ in range(nboot)]
    return est, np.std(boot)

def r2_from_z(beta, se, n):
    z2 = (beta / se) ** 2; return z2 / (z2 + n - 2)

def steiger(beta_e, se_e, n_e, beta_t, se_t, cases, controls):
    r2_exp = r2_from_z(beta_e, se_e, n_e).sum()
    # liability-scale approximation for a binary outcome (Lee et al. 2012 not applied; observed scale)
    r2_out = r2_from_z(beta_t, se_t, cases + controls).sum()
    return r2_exp > r2_out

def logbf(beta, se, prior_sd):
    v = se**2; r = prior_sd**2 / (prior_sd**2 + v); z = beta / se
    return 0.5 * (np.log(1 - r) + r * z**2)

def logsum(x):
    m = np.max(x); return m + np.log(np.sum(np.exp(x - m)))

def coloc_abf(d, sd_eqtl=0.15, sd_t2d=0.20, p1=1e-4, p2=1e-4, p12=1e-5):
    l1 = logbf(d.beta_eqtl.values, d.se_eqtl.values, sd_eqtl)
    l2 = logbf(d.beta_t2d.values, d.se_t2d.values, sd_t2d)
    s1, s2, s12 = logsum(l1), logsum(l2), logsum(l1 + l2)
    h = np.array([0.0, np.log(p1) + s1, np.log(p2) + s2,
                  0.0,
                  np.log(p12) + s12])
    # H3 term: p1 p2 (sum_i sum_j!=i) = p1 p2 (S1 S2 - S12)
    a, b = s1 + s2, s12
    h[3] = np.log(p1) + np.log(p2) + a + np.log1p(-np.exp(b - a)) if b < a else -np.inf
    pp = np.exp(h - logsum(h)); return pp

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("eqtl"); ap.add_argument("t2d"); ap.add_argument("instruments")
    ap.add_argument("--cases", type=int, required=True); ap.add_argument("--controls", type=int, required=True)
    ap.add_argument("--genes"); ap.add_argument("--out", default=".")
    a = ap.parse_args()
    e = pd.read_csv(a.eqtl, sep="\t"); t = pd.read_csv(a.t2d, sep="\t"); iv = pd.read_csv(a.instruments, sep="\t")
    m = e.merge(t, on="snp")
    genes = open(a.genes).read().split() if a.genes else sorted(iv.gene.unique())
    rows, crows = [], []
    for g in genes:
        x = m[m.gene == g]
        v = x[x.snp.isin(iv[iv.gene == g].snp)]
        if len(v):
            bx, by, sy = v.beta_eqtl.values, v.beta_t2d.values, v.se_t2d.values
            F = (bx / v.se_eqtl.values) ** 2
            b, s = ivw(bx, by, sy)
            r = dict(gene=g, n_iv=len(v), F_min=F.min(), IVW=b, IVW_p=2 * stats.norm.sf(abs(b / s)),
                     Steiger_ok=bool(steiger(bx, v.se_eqtl.values, v.n_eqtl.values[0], by, sy, a.cases, a.controls)))
            if len(v) >= 3:
                eb, es, ep, ip = egger(bx, by, sy); wb, ws = weighted_median(bx, by, sy)
                r.update(Egger=eb, Egger_p=ep, Egger_intercept_p=ip, WMedian=wb, WMedian_p=2 * stats.norm.sf(abs(wb / ws)))
            rows.append(r)
        if len(x) >= 10:
            pp = coloc_abf(x); crows.append(dict(gene=g, nsnps=len(x), **{f"PP.H{i}": pp[i] for i in range(5)}))
    pd.DataFrame(rows).to_csv(f"{a.out}/mr_sensitivity.tsv", sep="\t", index=False)
    pd.DataFrame(crows).to_csv(f"{a.out}/coloc.tsv", sep="\t", index=False)
    print(pd.DataFrame(rows).to_string(index=False)); print(pd.DataFrame(crows).round(3).to_string(index=False))

if __name__ == "__main__":
    main()
