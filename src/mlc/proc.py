"""Normalisation, gene-set scoring, clustering, state labels, display embedding."""
import numpy as np
import pandas as pd
from scipy import sparse
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans
from sklearn.manifold import TSNE
from . import genesets as G
from .coupling import fast_score, matched_controls


def normalize(ds, target=1e4):
    X = ds.X.tocsr().astype("float32"); tot = np.asarray(X.sum(axis=1)).ravel(); tot[tot == 0] = 1
    L = sparse.diags(target / tot) @ X; L.data = np.log1p(L.data); ds.L = sparse.csr_matrix(L); return ds


def gene_means(ds): return np.asarray(ds.L.mean(axis=0)).ravel()


def score_genes(ds, genes, seed=0):
    """Per-cell score for a gene set (mean minus expression-matched controls). Returns (scores, n_found)."""
    idx = ds.gene_idx(genes)
    if len(idx) < 3: return np.zeros(ds.n_obs), len(idx)
    gm = gene_means(ds); ctrl = matched_controls(gm, idx, seed=seed)
    return fast_score(ds.L, idx, ctrl), len(idx)


def prefix_genes(ds, prefixes): return [g for g in ds.var_names if g.startswith(tuple(prefixes))]


def hvg(ds, n_top=2000):
    m = gene_means(ds); sq = np.asarray(ds.L.multiply(ds.L).mean(axis=0)).ravel(); v = np.maximum(sq - m ** 2, 0)
    disp = np.log((v + 1e-9) / (m + 1e-9)); b = pd.qcut(m, 20, labels=False, duplicates="drop")
    z = pd.Series(disp).groupby(b).transform(lambda x: (x - x.mean()) / (x.std() + 1e-9)).values
    z = np.array(z, dtype=float); z[m < 0.0125] = -np.inf
    return np.argsort(-z)[:n_top]


def pca(ds, n_comps=30, n_top=2000, seed=0, center_by=None):
    idx = hvg(ds, n_top); X = ds.L[:, idx].toarray(); X = (X - X.mean(0)) / (X.std(0) + 1e-6); X = np.clip(X, -10, 10)
    if center_by:                                   # display only: remove per-sample offsets
        X = X - pd.DataFrame(X).groupby(ds.obs[center_by].values).transform("mean").values
    model = PCA(n_components=min(n_comps, X.shape[1] - 1), svd_solver="randomized", random_state=seed)
    ds.obsm["pca"] = model.fit_transform(X)
    return ds


def cluster(ds, k=8, seed=0):
    km = KMeans(n_clusters=k, n_init=10, random_state=seed)
    ds.obs["cluster"] = km.fit_predict(ds.obsm["pca"][:, :20]).astype(str)
    return ds


def add_scores(ds, panels, prefix="", seed=0):
    out = {}
    for name, genes in panels.items():
        s, n = score_genes(ds, genes, seed); ds.obs[prefix + name] = s; out[name] = n
    return out


def cluster_table(ds):
    t = ds.obs.groupby("cluster")[["micro", "nonmicro"]].mean(); t["n_cells"] = ds.obs.groupby("cluster").size()
    t["contaminant"] = t["nonmicro"] > t["micro"]; return t.round(3)


def label_states(ds, thr=0.5):
    """Cell-level state call from marker-panel scores: the most enriched ACTIVATED panel if its score exceeds thr,
    otherwise 'homeostatic'. Returns the mean panel score of every state (for the heat map)."""
    act = [k for k in G.STATE_PANELS if k != "homeostatic"]
    S = ds.obs[["st_" + k for k in act]].copy(); S.columns = act
    best = S.idxmax(axis=1); ds.obs["state"] = np.where(S.max(axis=1) > thr, best, "homeostatic")
    P = ds.obs[["st_" + k for k in G.STATE_PANELS]].copy(); P.columns = list(G.STATE_PANELS)
    return P.groupby(ds.obs["state"].values).mean().round(2)


def display_embedding(ds, max_cells=3000, seed=0):
    """2-D t-SNE for figures only (never used for statistics)."""
    rng = np.random.default_rng(seed); idx = rng.choice(ds.n_obs, min(max_cells, ds.n_obs), replace=False)
    P = ds.obsm["pca"][idx, :20]
    emb = TSNE(n_components=2, init="pca", perplexity=30, random_state=seed).fit_transform(P)
    return idx, emb
