"""Small, dependency-free statistics (numpy/scipy only)."""
import numpy as np
import pandas as pd
from scipy import stats


def ols(y, X, names=None):
    """Ordinary least squares. X must include the intercept column. Returns a tidy table."""
    y = np.asarray(y, float); X = np.asarray(X, float); n, k = X.shape
    XtXi = np.linalg.pinv(X.T @ X); beta = XtXi @ X.T @ y; res = y - X @ beta
    dof = max(n - k, 1); s2 = res @ res / dof; se = np.sqrt(np.diag(XtXi) * s2)
    t = beta / np.where(se > 0, se, np.nan); p = 2 * stats.t.sf(np.abs(t), dof)
    return pd.DataFrame({"beta": beta, "se": se, "t": t, "p": p}, index=names or [f"x{i}" for i in range(k)])


def design(df, cols, cohort=None):
    """Intercept + numeric columns (+ cohort dummies when cohort is given)."""
    X = [np.ones(len(df))] + [df[c].to_numpy(float) for c in cols]; names = ["const"] + list(cols)
    if cohort is not None:
        for c in sorted(df[cohort].unique())[1:]:
            X.append((df[cohort] == c).to_numpy(float)); names.append(f"{cohort}_{c}")
    return np.column_stack(X), names


def ols_multi(Y, X):
    """Fit many responses at once (columns of Y). Returns beta, se, p for every column and predictor."""
    Y = np.asarray(Y, float); X = np.asarray(X, float); n, k = X.shape
    XtXi = np.linalg.pinv(X.T @ X); B = XtXi @ X.T @ Y; R = Y - X @ B
    dof = max(n - k, 1); s2 = (R ** 2).sum(0) / dof
    SE = np.sqrt(np.outer(np.diag(XtXi), s2)); T = B / np.where(SE > 0, SE, np.nan)
    return B, SE, 2 * stats.t.sf(np.abs(T), dof)


def bh_fdr(p):
    p = np.asarray(p, float); n = len(p); order = np.argsort(p); ranked = p[order] * n / (np.arange(n) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]; out = np.empty(n); out[order] = np.minimum(ranked, 1)
    return out


def zscore_rows(df):
    return df.sub(df.mean(axis=1), axis=0).div(df.std(axis=1).replace(0, np.nan), axis=0)
