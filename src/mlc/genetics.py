"""Human genetics helpers: orthology, gene-set enrichment vs matched random sets, Bayesian colocalisation."""
import numpy as np
import pandas as pd
from scipy.special import logsumexp


def load_orthologs(path):
    """MGI HOM_MouseHumanSequence.rpt -> table of mouse/human symbol pairs."""
    t = pd.read_csv(path, sep="\t")
    h = t[t["Common Organism Name"] == "human"][["DB Class Key", "Symbol", "EntrezGene ID"]]
    m = t[t["Common Organism Name"] == "mouse, laboratory"][["DB Class Key", "Symbol"]]
    return m.merge(h, on="DB Class Key", suffixes=("_mouse", "_human"))


def to_human(genes, orth):
    return sorted(set(orth[orth.Symbol_mouse.isin(genes)].Symbol_human))


def set_enrichment(gene_table, hset, expressed, n_perm=10000, seed=0, score="ZSTAT", exclude=None):
    """Mean gene-level GWAS Z of hset vs random same-size sets drawn from microglia-EXPRESSED genes only."""
    g = gene_table if exclude is None else gene_table[~exclude]
    obs_genes = g[g.SYMBOL.isin(hset)]; pool = g[g.SYMBOL.isin(expressed)][score].values
    rng = np.random.default_rng(seed); n = len(obs_genes); obs = obs_genes[score].mean()
    null = np.array([rng.choice(pool, n, replace=False).mean() for _ in range(n_perm)])
    return dict(n_genes=n, mean_z=obs, null_mean=null.mean(), p=(1 + (null >= obs).sum()) / (n_perm + 1)), null


def _labf(beta, se, sd):
    V = se ** 2; z = beta / se; r = sd ** 2 / (sd ** 2 + V)
    return 0.5 * (np.log(1 - r) + r * z ** 2)


def coloc_abf(b1, se1, b2, se2, type1="quant", type2="quant", p1=1e-4, p2=1e-4, p12=1e-5):
    """Approximate-Bayes-factor colocalisation (Giambartolomei 2014; Wakefield ABFs). Returns PP.H0..H4."""
    sd = {"quant": 0.15, "cc": 0.2}
    l1 = _labf(np.asarray(b1), np.asarray(se1), sd[type1]); l2 = _labf(np.asarray(b2), np.asarray(se2), sd[type2])
    s1, s2, s12 = logsumexp(l1), logsumexp(l2), logsumexp(l1 + l2)
    h3 = np.log(p1) + np.log(p2) + (s1 + s2) + np.log1p(-np.exp(np.clip(s12 - (s1 + s2), -700, -1e-12)))
    lH = np.array([0.0, np.log(p1) + s1, np.log(p2) + s2, h3, np.log(p12) + s12])
    pp = np.exp(lH - logsumexp(lH)); return dict(zip(["PP.H0", "PP.H1", "PP.H2", "PP.H3", "PP.H4"], pp))
