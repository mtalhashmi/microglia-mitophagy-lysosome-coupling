"""Quality control: per-sample MAD filters and a marker-based doublet/contamination filter."""
import numpy as np
import pandas as pd
from scipy import stats


def add_qc(ds, species="mouse"):
    prefix = "mt-" if species == "mouse" else "MT-"
    mt = ds.var_names.str.startswith(prefix)
    tot = np.asarray(ds.X.sum(axis=1)).ravel(); ng = np.asarray((ds.X > 0).sum(axis=1)).ravel()
    mtc = np.asarray(ds.X[:, np.where(mt)[0]].sum(axis=1)).ravel() if mt.any() else np.zeros(ds.n_obs)
    ds.obs["total_counts"] = tot; ds.obs["n_genes"] = ng
    ds.obs["pct_mt"] = 100 * mtc / np.maximum(tot, 1); ds.obs["depth"] = np.log1p(tot)
    return ds


def _mad_out(x, n):
    x = np.asarray(x, float); med = np.median(x); mad = stats.median_abs_deviation(x)
    return (x < med - n * mad) | (x > med + n * mad)


def mad_filter(ds, nmads=5, mt_cap=15.0, min_genes=200):
    """Remove outlier cells separately inside every sample. Returns (filtered dataset, log table)."""
    keep = np.ones(ds.n_obs, bool); log = []
    for s in ds.obs["sample_id"].unique():
        m = np.where((ds.obs["sample_id"] == s).values)[0]; o = ds.obs.iloc[m]
        bad = (_mad_out(np.log1p(o["total_counts"]), nmads) | _mad_out(np.log1p(o["n_genes"]), nmads)
               | (o["pct_mt"].values > mt_cap) | (o["n_genes"].values < min_genes))
        keep[m[bad]] = False
        log.append(dict(sample_id=s, n_before=len(m), n_removed=int(bad.sum())))
    return ds.subset_cells(keep), pd.DataFrame(log)
