"""
Disposition-index (multiplicative) T2D disease model, fitted to three drugs.

Replaces the additive equilibrium  Q = q0 - phi_S*S - phi_B*B + phi_H*H
with a glucose balance in which insulin action is the PRODUCT of secretion and sensitivity:

    hepatic production          = H                       (H: hepatic glucose output, [0,1])
    glucose disposal            = Q * (g0 + k * S * B * Q)  (insulin ~ B*Q, action ~ S*insulin)
    steady state  ->  k*S*B*Q^2 + g0*Q - H = 0   ->   Q = (-g0 + sqrt(g0^2 + 4 k S B H)) / (2 k S B)
    HbA1c (%)     = (Q + 2.59) / 1.59            (ADAG relation)

Scaling: a healthy profile (B = S = 0.9, H = 0.3) gives Q = 5.3 mmol/L; u is the share of
disposal that is insulin-independent in health, so g0 = u*0.3/5.3 and k = (1-u)*0.3/(0.81*5.3^2).

Unknowns (u, cohort mean B0, cohort mean S0) get weak priors and are fitted by random-walk
Metropolis MCMC (common random numbers for the virtual cohort) against four targets:
    baseline HbA1c 7.6%  |  sulfonylurea  -1.51% (Hirst 2013)  [maximal beta-cell drive: B -> 1]
    thiazolidinedione  -1.30% (Aronoff 2000, -1.00 to -1.60)    [glucose disposal +54%: S x 1.54]
    metformin          -1.12% (Hirst 2012, 0.92-1.32)           [EGP -19%: H x 0.81; disposal +13%: S x 1.13]
Leave-one-drug-out: fit on baseline + two drugs, predict the third.

Outputs: drug_check.tsv, restoration.tsv, restoration.png
Usage:   python3 t2d_disposition_model.py [--draws 6000] [--people 200] [--seed 20261018]
"""
import argparse
import numpy as np
import pandas as pd

HB = lambda q: (q + 2.59) / 1.59
TARGETS = {"baseline": (7.60, 0.20), "SU": (-1.51, 0.14), "TZD": (-1.30, 0.15), "MET": (-1.12, 0.10)}

def glucose(B, S, H, u):
    g0 = u * 0.3 / 5.3
    k = (1 - u) * 0.3 / (0.81 * 5.3 ** 2)
    a = k * np.clip(S, 1e-3, 1) * np.clip(B, 1e-3, 1)
    return (-g0 + np.sqrt(g0 ** 2 + 4 * a * H)) / (2 * a)

def cohort(rng, people):
    return rng.normal(0, 0.10, (3, people))          # fixed person-level deviations (B, S, H)

def outcomes(th, z):
    u, B0, S0 = th
    B = np.clip(B0 + z[0], 0.05, 1); S = np.clip(S0 + z[1], 0.05, 1); H = np.clip(0.65 + z[2], 0.3, 1)
    base = HB(glucose(B, S, H, u))
    return {"baseline": base.mean(),
            "SU": (HB(glucose(np.ones_like(B), S, H, u)) - base).mean(),
            "TZD": (HB(glucose(B, np.minimum(S * 1.54, 1), H, u)) - base).mean(),
            "MET": (HB(glucose(B, np.minimum(S * 1.13, 1), H * 0.81, u)) - base).mean()}, (B, S, H, base)

def logpost(th, z, keys):
    u, B0, S0 = th
    if not (0 <= u <= 0.6 and 0.1 <= B0 <= 0.85 and 0.1 <= S0 <= 0.85):
        return -np.inf
    lp = -0.5 * (((B0 - 0.5) / 0.2) ** 2 + ((S0 - 0.5) / 0.2) ** 2)   # weak priors
    o, _ = outcomes(th, z)
    return lp - 0.5 * sum(((o[k] - TARGETS[k][0]) / TARGETS[k][1]) ** 2 for k in keys)

def mcmc(z, keys, rng, n=6000, burn=2000):
    th = np.array([0.2, 0.5, 0.5]); lp = logpost(th, z, keys); step = np.array([0.05, 0.03, 0.03]); ch = []
    for i in range(n):
        prop = th + rng.normal(0, step); lq = logpost(prop, z, keys)
        if np.log(rng.uniform()) < lq - lp:
            th, lp = prop, lq
        if i >= burn:
            ch.append(th.copy())
    return np.array(ch)[::4]

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--draws", type=int, default=6000); ap.add_argument("--people", type=int, default=200)
    ap.add_argument("--seed", type=int, default=20261018)
    a = ap.parse_args(); rng = np.random.default_rng(a.seed)

    z = cohort(rng, a.people)
    rows = []
    for held in ["SU", "TZD", "MET", None]:
        keys = ["baseline"] + [k for k in ("SU", "TZD", "MET") if k != held]
        ch = mcmc(z, keys, rng, n=a.draws)
        pr = {k: np.array([outcomes(t, z)[0][k] for t in ch]) for k in TARGETS}
        for k in TARGETS:
            m, lo, hi = np.percentile(pr[k], [50, 5, 95])
            rows.append(dict(fit_on="+".join(keys), predicted=k, out_of_sample=(k == held),
                             median=round(m, 2), q05=round(lo, 2), q95=round(hi, 2), trial=TARGETS[k][0]))
    dc = pd.DataFrame(rows); dc.to_csv("drug_check.tsv", sep="\t", index=False)
    print(dc[dc.out_of_sample | (dc.fit_on.str.count("\\+") == 3)].to_string(index=False))

    for i, nm in enumerate(["u", "B0", "S0"]):
        print(nm, np.round(np.percentile(ch[:, i], [50, 5, 95]), 2))
    res, full = [], {"X": [], "Y": [], "Z": []}
    curves = {d: [] for d in ("X beta-cell", "Y insulin sensitivity", "Z liver output")}
    fs = np.linspace(0, 1, 11)
    for t in ch:
        _, (Bp, Sp, Hp, bp) = outcomes(t, z); u = t[0]
        for d in curves:
            row = []
            for f in fs:
                Bn, Sn, Hn = {"X beta-cell": (Bp + f * (1 - Bp), Sp, Hp), "Y insulin sensitivity": (Bp, Sp + f * (1 - Sp), Hp),
                              "Z liver output": (Bp, Sp, Hp - f * (Hp - 0.3))}[d]
                new = HB(glucose(Bn, Sn, Hn, u)); row.append(((new - bp).mean(), (new < 6.5).mean()))
            curves[d].append(row)
    for d in curves:
        arr = np.array(curves[d])                                # draws x shares x 2
        for f in (0.2, 0.5, 1.0):
            j = int(round(f * 10)); m, lo, hi = np.percentile(arr[:, j, 0], [50, 5, 95])
            res.append(dict(defect=d, share=f, median=round(m, 2), q05=round(lo, 2), q95=round(hi, 2),
                            below_6_5=round(float(arr[:, j, 1].mean()), 2) if f == 1.0 else None))
    rt = pd.DataFrame(res); rt.to_csv("restoration.tsv", sep="\t", index=False); print(rt.to_string(index=False))
    X, Y, Zl = (np.array(curves[d])[:, -1, 0] for d in curves)
    print("P(X>Y) =", round(float((X < Y).mean()), 2), " P(Y>Z) =", round(float((Y < Zl).mean()), 2),
          " P(Z>X) =", round(float((Zl < X).mean()), 2))

    try:
        import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
        fig, ax = plt.subplots(figsize=(7, 4), dpi=200)
        for d, c in zip(curves, ["#2a6db5", "#d07a1f", "#3c9a5f"]):
            arr = np.array(curves[d])[:, :, 0]; q = np.percentile(arr, [50, 5, 95], axis=0)
            ax.plot(fs * 100, q[0], color=c, lw=2, label=d); ax.fill_between(fs * 100, q[1], q[2], color=c, alpha=0.15)
        ax.set_xlabel("Share of lost function restored (%)"); ax.set_ylabel("Change in HbA1c at steady state (%)")
        ax.set_title("Disposition-index model fitted to three drugs (median, 5-95%)"); ax.legend(frameon=False)
        plt.tight_layout(); plt.savefig("restoration.png")
    except ImportError:
        pass

if __name__ == "__main__":
    main()
