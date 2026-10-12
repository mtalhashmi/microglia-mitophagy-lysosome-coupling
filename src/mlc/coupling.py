"""The Mitophagy-Lysosome Coupling Index (MLCI) and its statistics."""
import numpy as np
import pandas as pd
from scipy import stats, sparse
from .stats import ols, design


def partial_spearman(x, y, covars=None):
    """Spearman correlation of x and y after regressing covariates out of the ranks."""
    rx = stats.rankdata(x).astype(float); ry = stats.rankdata(y).astype(float)
    if covars is not None and covars.shape[1] > 0:
        Z = np.column_stack([np.ones(len(rx)), covars])
        rx = rx - Z @ np.linalg.lstsq(Z, rx, rcond=None)[0]
        ry = ry - Z @ np.linalg.lstsq(Z, ry, rcond=None)[0]
    return float(np.corrcoef(rx, ry)[0, 1])


def coupling_per_sample(obs, a, b, sample_col="sample_id", covars=(), min_cells=100, n_boot=100, seed=0):
    """One MLCI per biological sample: partial Spearman r, Fisher z, bootstrap 95% CI."""
    rng = np.random.default_rng(seed); rows = []
    for s, d in obs.groupby(sample_col, observed=True):
        if len(d) < min_cells: continue
        xa, xb = d[a].to_numpy(float), d[b].to_numpy(float)
        C = d[list(covars)].to_numpy(float) if len(covars) else None
        r = partial_spearman(xa, xb, C); boots = []
        for _ in range(n_boot):
            i = rng.integers(0, len(d), len(d))
            boots.append(partial_spearman(xa[i], xb[i], None if C is None else C[i]))
        lo, hi = np.percentile(boots, [2.5, 97.5]) if n_boot else (np.nan, np.nan)
        rows.append(dict(sample_id=s, n_cells=len(d), r=r, z=np.arctanh(np.clip(r, -0.999, 0.999)),
                         ci_lo=lo, ci_hi=hi))
    return pd.DataFrame(rows)


def _rowmean(X, idx):
    return np.asarray(X[:, idx].mean(axis=1)).ravel()


def fast_score(X, gene_idx, ctrl_idx):
    """Gene-set score = mean(set genes) - mean(expression-matched control genes), per cell."""
    return _rowmean(X, gene_idx) - _rowmean(X, ctrl_idx)


def _bins(gene_means, n_bins):
    order = np.argsort(gene_means); b = np.empty(len(gene_means), int)
    b[order] = np.arange(len(gene_means)) * n_bins // len(gene_means)
    return b


def matched_controls(gene_means, set_idx, n_ctrl=5, n_bins=25, seed=0):
    rng = np.random.default_rng(seed); b = _bins(gene_means, n_bins)
    pool = np.setdiff1d(np.arange(len(gene_means)), set_idx); ctrl = []
    for g in set_idx:
        cand = pool[b[pool] == b[g]]
        ctrl.extend(rng.choice(cand, size=min(n_ctrl, len(cand)), replace=False))
    return np.unique(ctrl)


def matched_random_set(gene_means, set_idx, n_bins=25, rng=None):
    """Random gene set with the same size and the same expression-level profile as set_idx."""
    rng = rng or np.random.default_rng(); b = _bins(gene_means, n_bins)
    pool = np.setdiff1d(np.arange(len(gene_means)), set_idx); out = []
    for g in set_idx:
        cand = np.setdiff1d(pool[b[pool] == b[g]], out)
        if len(cand) == 0: cand = np.setdiff1d(pool, out)
        out.append(rng.choice(cand))
    return np.array(out)


def age_slope(z, age, cohort):
    """Pooled within-cohort slope of z on age (each cohort centred separately)."""
    df = pd.DataFrame(dict(z=np.asarray(z, float), age=np.asarray(age, float), c=np.asarray(cohort)))
    df["z"] -= df.groupby("c")["z"].transform("mean"); df["age"] -= df.groupby("c")["age"].transform("mean")
    den = (df["age"] ** 2).sum()
    return float((df["age"] * df["z"]).sum() / den) if den > 0 else np.nan


def dersimonian_laird(beta, se):
    """Random-effects meta-analysis of per-cohort slopes."""
    beta, se = np.asarray(beta, float), np.asarray(se, float)
    w = 1 / se ** 2; fe = np.sum(w * beta) / np.sum(w); Q = np.sum(w * (beta - fe) ** 2); k = len(beta)
    C = np.sum(w) - np.sum(w ** 2) / np.sum(w)
    tau2 = max(0.0, (Q - (k - 1)) / C) if C > 0 else 0.0
    ws = 1 / (se ** 2 + tau2); re = np.sum(ws * beta) / np.sum(ws); se_re = np.sqrt(1 / np.sum(ws))
    return dict(beta=re, se=se_re, ci_lo=re - 1.96 * se_re, ci_hi=re + 1.96 * se_re, z=re / se_re,
                p=2 * stats.norm.sf(abs(re / se_re)), tau2=tau2,
                I2=max(0.0, (Q - (k - 1)) / Q) if Q > 0 else 0.0, k=k)


def cohort_slopes(R, y="z", age="log2_age", by=("level", "cohort")):
    """OLS slope of y on age inside every group (cohort x adjustment level)."""
    rows = []
    for key, d in R.groupby(list(by)):
        if d[age].nunique() < 2 or len(d) < 4: continue
        X, names = design(d, [age]); f = ols(d[y], X, names).loc[age]
        keys = key if isinstance(key, tuple) else (key,)
        rows.append({**dict(zip(by, keys)), "n": len(d), "beta": f.beta, "se": f.se, "p": f.p})
    return pd.DataFrame(rows)


def perm_null_age_slope(Xs, gene_means, setA, setB, age, cohort, n_perm=200, seed=0):
    """Null distribution of the age slope for random expression-matched gene-set pairs.
    Xs: dict sample -> cells x genes (log-normalised, dense or sparse); age/cohort follow the order of Xs."""
    rng = np.random.default_rng(seed); keys = list(Xs); null = []
    for _ in range(n_perm):
        ra = matched_random_set(gene_means, setA, rng=rng); rb = matched_random_set(gene_means, setB, rng=rng)
        ca = matched_controls(gene_means, ra, seed=int(rng.integers(1e9)))
        cb = matched_controls(gene_means, rb, seed=int(rng.integers(1e9)))
        z = [np.arctanh(np.clip(partial_spearman(fast_score(Xs[k], ra, ca), fast_score(Xs[k], rb, cb)),
                                -.999, .999)) for k in keys]
        null.append(age_slope(np.array(z), age, cohort))
    return np.array(null)
