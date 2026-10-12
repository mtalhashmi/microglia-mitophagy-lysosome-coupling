"""Gene rewiring, replicated aging DE, TF activity, human bulk test, candidate ranking."""
import numpy as np
import pandas as pd
from scipy import stats
from .stats import ols, design, ols_multi, bh_fdr
from . import genesets as G


# ------------------------------------------------------------ rewiring (notebook 09)
def partial_corr_matrix(X, covars):
    R = stats.rankdata(X, axis=0).astype(float); Z = np.column_stack([np.ones(len(R)), covars])
    R = R - Z @ np.linalg.lstsq(Z, R, rcond=None)[0]
    return np.corrcoef(R, rowvar=False)


def cross_arm_z(ds, genes_a, genes_b, covar_cols, min_cells=100):
    """Per sample: Fisher-z matrix of partial correlations between arm A genes and arm B genes."""
    A = [g for g in genes_a if g in ds.var_names]; B = [g for g in genes_b if g in ds.var_names]
    idx = ds.gene_idx(A + B); out = {}
    for s in ds.obs["sample_id"].unique():
        m = (ds.obs["sample_id"] == s).values
        if m.sum() < min_cells: continue
        X = ds.L[m][:, idx].toarray(); C = partial_corr_matrix(X, ds.obs.loc[m, covar_cols].values)[: len(A), len(A):]
        out[s] = pd.DataFrame(np.arctanh(np.clip(C, -.999, .999)), index=A, columns=B)
    return out


def rewiring_test(Z, age, cohort, n_perm=2000, seed=0):
    """Per gene: how much weaker is its link to the other arm in old vs young samples?
    Permutation shuffles old/young labels within cohort."""
    ids = list(Z)
    A = sorted(set.intersection(*[set(Z[s].index) for s in ids]))
    B = sorted(set.intersection(*[set(Z[s].columns) for s in ids]))
    M = np.stack([Z[s].loc[A, B].values for s in ids]); age = np.asarray(age, float); cohort = np.asarray(cohort)
    old = np.zeros(len(ids), bool)
    for c in np.unique(cohort):
        i = np.where(cohort == c)[0]; old[i] = age[i] >= np.median(age[i])
    d = M[~old].mean(0) - M[old].mean(0)
    oa, ob = d.mean(1), d.mean(0)
    rng = np.random.default_rng(seed)
    na, nb = [], []
    for _ in range(n_perm):
        perm = old.copy()
        for c in np.unique(cohort):
            i = np.where(cohort == c)[0]; perm[i] = rng.permutation(old[i])
        dd = M[~perm].mean(0) - M[perm].mean(0); na.append(dd.mean(1)); nb.append(dd.mean(0))
    pa = (1 + (np.abs(np.array(na)) >= np.abs(oa)).sum(0)) / (n_perm + 1)
    pb = (1 + (np.abs(np.array(nb)) >= np.abs(ob)).sum(0)) / (n_perm + 1)
    t = pd.concat([pd.DataFrame(dict(gene=A, arm="mitophagy", delta=oa, p=pa)),
                   pd.DataFrame(dict(gene=B, arm="lysosome", delta=ob, p=pb))])
    t["fdr"] = bh_fdr(t.p.values)
    return t.sort_values("p").reset_index(drop=True), M[~old].mean(0), M[old].mean(0), A, B


# ------------------------------------------------------------ replicated aging DE (notebook 10)
def pseudobulk(ds, min_cells=100):
    cols = {}
    for s in ds.obs["sample_id"].unique():
        m = (ds.obs["sample_id"] == s).values
        if m.sum() >= min_cells: cols[s] = np.asarray(ds.X[m].sum(axis=0)).ravel()
    return pd.DataFrame(cols, index=ds.var_names)


def age_de(counts, age_log2, min_cpm=1.0):
    """Linear model on log2-CPM per gene (cohort-wise). Replicate-aware; a lightweight alternative to DESeq2."""
    cpm = counts / counts.sum() * 1e6; keep = (cpm > min_cpm).mean(axis=1) >= 0.5
    Y = np.log2(cpm[keep] + 1).T.values; X = np.column_stack([np.ones(len(age_log2)), age_log2])
    B, SE, P = ols_multi(Y, X)
    z = np.sign(B[1]) * stats.norm.isf(np.clip(P[1], 1e-300, 1) / 2)
    return pd.DataFrame({"beta": B[1], "se": SE[1], "p": P[1], "z": z}, index=counts.index[keep])


def stouffer(zs):
    Z = pd.concat(zs, axis=1).dropna(); k = Z.shape[1]
    out = pd.DataFrame({"meta_z": Z.sum(axis=1) / np.sqrt(k)})
    out["p"] = 2 * stats.norm.sf(out.meta_z.abs()); out["fdr"] = bh_fdr(out.p.values)
    out["concordant"] = (np.sign(Z).nunique(axis=1) == 1); return out.join(Z.add_prefix("z_"))


# ------------------------------------------------------------ TF activity (notebook 11)
def get_net(organism="mouse"):
    """TF-target network. Uses decoupler's CollecTRI when installed (pip install decoupler==1.8.0);
    otherwise a built-in mini network (TFEB -> CLEAR targets)."""
    try:
        import decoupler as dc
        net = dc.get_collectri(organism=organism, split_complexes=False)[["source", "target", "weight"]]
        print("using CollecTRI network:", net.source.nunique(), "TFs"); return net
    except Exception as e:
        print("decoupler not available -> built-in TFEB-family mini network (", type(e).__name__, ")")
        return pd.DataFrame([("Tfeb", g, 1.0) for g in G.CLEAR_TARGETS], columns=["source", "target", "weight"])


def tf_activity(ds, net, min_targets=5):
    """Per-cell activity = weighted mean of z-scored target expression (excluding the TF itself)."""
    out = {}
    for tf, sub in net.groupby("source"):
        sub = sub[(sub.target != tf)]; idx = ds.gene_idx(sub.target)
        if len(idx) < min_targets: continue
        w = sub.set_index("target").weight.reindex(ds.var_names[idx]).fillna(1).values
        X = ds.L[:, idx].toarray(); X = (X - X.mean(0)) / (X.std(0) + 1e-6)
        out[tf] = X @ w / np.sum(np.abs(w))
    return pd.DataFrame(out, index=ds.obs_names)


# ------------------------------------------------------------ human bulk (notebook 12)
def bulk_scores(counts, set_a, set_b):
    cpm = counts / counts.sum() * 1e6; l = np.log2(cpm + 1); l = l[l.mean(axis=1) > 1]
    z = l.sub(l.mean(axis=1), axis=0).div(l.std(axis=1) + 1e-9, axis=0)
    ia, ib = [g for g in set_a if g in z.index], [g for g in set_b if g in z.index]
    return pd.DataFrame({"mito": z.loc[ia].mean(), "lyso": z.loc[ib].mean()}), (len(ia), len(ib))


def interaction_test(d, covars=(), n_perm=5000, seed=0):
    """lyso ~ mito * age (+covariates). Permutation p for the mito:age interaction."""
    d = d.copy()
    d["age_s"] = (d.age - d.age.mean()) / d.age.std()
    d["mito_s"] = (d.mito - d.mito.mean()) / d.mito.std()
    d["inter"] = d.age_s * d.mito_s
    X = np.column_stack([np.ones(len(d)), d.mito_s, d.age_s, d.inter] + [d[c].to_numpy(float) for c in covars])
    f = ols(d.lyso, X, ["const", "mito", "age", "mito:age"] + list(covars)); obs = f.loc["mito:age", "t"]
    rng = np.random.default_rng(seed); null = []
    for _ in range(n_perm):
        a = rng.permutation(d.age_s.values); Xp = X.copy(); Xp[:, 2] = a; Xp[:, 3] = a * d.mito_s.values
        null.append(ols(d.lyso, Xp).iloc[3].t)
    p = (1 + (np.abs(null) >= abs(obs)).sum()) / (n_perm + 1); return f, p, np.array(null)


# ------------------------------------------------------------ candidate ranking (notebook 15)
def rank_product(E, cols, n_perm=5000, seed=0):
    R = E[cols].fillna(0).rank(ascending=False, method="average"); rp = np.exp(np.log(R).mean(axis=1))
    rng = np.random.default_rng(seed); nmin = np.empty(n_perm)
    for i in range(n_perm):
        sh = np.column_stack([rng.permutation(R[c].values) for c in cols])
        nmin[i] = np.exp(np.log(sh).mean(axis=1)).min()
    E = E.copy(); E["rank_product"] = rp; E["p_familywise"] = [(1 + (nmin <= v).sum()) / (n_perm + 1) for v in rp]
    return E.sort_values("rank_product")
