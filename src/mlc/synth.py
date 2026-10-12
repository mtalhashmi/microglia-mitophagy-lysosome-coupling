"""Simulated sample data. Purpose: run and understand the whole pipeline in minutes, with a KNOWN answer.
Built in: MLCI falls with age, specific to the mitophagy-lysosome pair; DAM cells accumulate with age;
a handful of 'interface' genes carry the coupling; TFEB-regulon activity tracks coupling."""
import json
import numpy as np
import pandas as pd
from scipy import sparse
from . import genesets as G
from .config import DATA, RAW, PROC, SEED
from .dataset import Dataset

INTERFACE_MITO = ["Pink1", "Prkn", "Bnip3", "Bnip3l", "Optn", "Tbk1"]
INTERFACE_LYSO = ["Ctsd", "Ctsb", "Lamp1", "Lamp2", "Atp6v0a1", "Gba", "Mcoln1", "Tpp1"]


def gene_panel(n_total=1500):
    ribo = [f"Rpl{i}" for i in range(1, 31)] + [f"Rps{i}" for i in range(1, 21)]
    ox = [f"Ndufa{i}" for i in range(1, 14)] + [f"Ndufb{i}" for i in range(1, 12)] + ["Sdha", "Sdhb", "Sdhc", "Sdhd",
          "Uqcrc1", "Uqcrc2", "Uqcrb", "Uqcrq", "Cox4i1", "Cox5a", "Cox5b", "Cox6a1", "Cox6b1", "Cox6c", "Cox7a2", "Cox7c",
          "Atp5f1a", "Atp5f1b", "Atp5f1c", "Atp5pb"]
    mt = ["mt-Co1", "mt-Co2", "mt-Co3", "mt-Nd1", "mt-Nd2", "mt-Nd4", "mt-Cytb", "mt-Atp6"]
    core = (G.MITOPHAGY + G.LYSOSOME + G.AUTOPHAGY_CORE + G.PROTEASOME + G.MICROGLIA + G.NON_MICROGLIA + G.DISSOCIATION
            + sum(G.STATE_PANELS.values(), []) + ribo + ox + mt + ["Tfeb", "Tfe3", "Mitf", "Apoe", "Cd68"] + G.CLEAR_TARGETS)
    core = list(dict.fromkeys(core)); fill = [f"Gene{i:04d}" for i in range(n_total - len(core))]
    return core + fill


def _sample_counts(rng, genes, gi, age, n, truth):
    G_ = len(genes); mu = truth["mu"]; beta = truth["beta"]; idx = lambda L: np.array([gi[g] for g in L if g in gi])
    dam_p = np.clip(0.10 + 0.25 * (age - 3) / 21, 0.02, 0.45)
    st = rng.choice(4, n, p=[1 - dam_p - 0.08, dam_p, 0.05, 0.03])            # homeo, DAM, IFN, prolif
    contam, dbl, lowq = rng.random(n) < .04, rng.random(n) < .02, rng.random(n) < .03
    u, v, t, r_ox, e_rib, e_pro, e_core = rng.normal(size=(7, n)); tau = rng.normal(0, 0.10)
    rho = np.clip(0.78 - 0.02 * (age - 3) + tau, 0.1, 0.92)
    F_l = rho * u + np.sqrt(1 - rho ** 2) * v + 0.3 * t; F_m = u + 0.3 * r_ox
    eta = np.tile(mu, (n, 1)) + beta[None, :] * ((age - 12) / 12)
    Lm = np.full(len(G.MITOPHAGY), 0.10); Ll = np.full(len(G.LYSOSOME), 0.10)
    for i, g in enumerate(G.MITOPHAGY): Lm[i] = 0.45 if g in INTERFACE_MITO else 0.10
    for i, g in enumerate(G.LYSOSOME): Ll[i] = 0.50 if g in INTERFACE_LYSO else 0.10
    eta[:, idx(G.MITOPHAGY)] += Lm * F_m[:, None]; eta[:, idx(G.LYSOSOME)] += Ll * F_l[:, None] + 2.0 * tau
    eta[:, idx(G.CLEAR_TARGETS)] += 0.15 * F_l[:, None] + 4.0 * tau
    eta[:, [gi["Tfeb"]]] += 0.4 * F_l[:, None] + 1.0 * tau
    ox = idx([g for g in genes if g.startswith(G.OXPHOS_PREFIX)]); eta[:, ox] += 0.25 * r_ox[:, None]
    rb = idx([g for g in genes if g.startswith(G.RIBOSOME_PREFIX)]); eta[:, rb] += 0.25 * e_rib[:, None]
    eta[:, idx(G.PROTEASOME)] += 0.25 * (0.5 * t + 0.87 * e_pro)[:, None]
    eta[:, idx(G.AUTOPHAGY_CORE)] += 0.25 * (0.6 * t + 0.8 * e_core)[:, None]
    hs, dm, ifn, pr = (idx(G.STATE_PANELS[k]) for k in ("homeostatic", "DAM", "interferon", "proliferating"))
    for k, (genes_k, amp) in enumerate([(hs, .8), (dm, 1.3), (ifn, 1.6), (pr, 1.6)]):
        eta[np.ix_(st == k, genes_k)] += amp
    eta[np.ix_(st == 1, idx(G.LYSOSOME))] += 0.4                                  # DAM raises lysosome genes (confound)
    stress = rng.exponential(0.6, n); eta[:, idx(G.DISSOCIATION)] += 0.9 * stress[:, None]
    mic, non = idx(G.MICROGLIA), idx(G.NON_MICROGLIA)
    eta[np.ix_(contam, mic)] -= 2.5; eta[np.ix_(contam, non)] += 2.8; eta[np.ix_(dbl, non)] += 1.8
    mt = idx([g for g in genes if g.startswith("mt-")]); eta[np.ix_(lowq, mt)] += 2.8
    depth = np.exp(rng.normal(7.4, 0.4, n)) * np.where(lowq, 0.4, 1.0)
    P = np.exp(eta - eta.max(1, keepdims=True)); P /= P.sum(1, keepdims=True)
    X = rng.poisson(depth[:, None] * P).astype("float32")
    truth_cells = pd.DataFrame({"true_state": np.array(["homeostatic", "DAM", "interferon", "proliferating"])[st],
                                "true_contam": contam, "true_doublet": dbl, "true_rho": rho})
    return X, truth_cells


def make_sample_data(seed=SEED, cells_per_sample=450, n_genes=1500):
    rng = np.random.default_rng(seed); genes = gene_panel(n_genes); gi = {g: i for i, g in enumerate(genes)}
    mu = rng.normal(-1.8, 1.1, len(genes))
    for g in G.MITOPHAGY + G.LYSOSOME + G.AUTOPHAGY_CORE + G.PROTEASOME: mu[gi[g]] = rng.normal(-1.0, 0.4)
    for g in G.MICROGLIA + sum(G.STATE_PANELS.values(), []): mu[gi[g]] = rng.normal(-0.6, 0.5)
    for g in G.NON_MICROGLIA: mu[gi[g]] = -4.0
    for g in genes:
        if g.startswith("mt-"): mu[gi[g]] = -1.0
        if g.startswith(G.RIBOSOME_PREFIX): mu[gi[g]] = -0.3
    beta = rng.normal(0, 0.04, len(genes))
    for g in G.MITOPHAGY + G.LYSOSOME: beta[gi[g]] = -0.08
    for g in G.CLEAR_TARGETS: beta[gi[g]] = -0.12
    for g in INTERFACE_MITO: beta[gi[g]] = -0.20
    for g in INTERFACE_LYSO: beta[gi[g]] = -0.15
    for g in ["Apoe", "Cd68", "Lpl", "Cst7", "Itgax"]: beta[gi[g]] = 0.25
    truth = {"mu": mu, "beta": beta}
    cohorts = {"SIM_A": [3] * 4 + [12] * 3 + [24] * 4, "SIM_B": [4] * 4 + [22] * 4}; sheet = []
    PROC.mkdir(parents=True, exist_ok=True)
    for c, ages in cohorts.items():
        Xs, obs = [], []
        for j, a in enumerate(ages):
            sid = f"{c}_m{j+1:02d}"; X, tc = _sample_counts(rng, genes, gi, a, cells_per_sample, truth)
            o = tc.copy(); o.index = [f"{sid}_cell{i:04d}" for i in range(len(o))]; o["sample_id"] = sid; o["cohort"] = c; o["age_months"] = a
            Xs.append(sparse.csr_matrix(X)); obs.append(o)
            sheet.append(dict(cohort=c, sample_id=sid, gsm=f"GSMSIM{len(sheet)+1:03d}", species="mouse", tissue="microglia",
                              region="hippocampus", age_months=a, sex="male" if j % 2 else "female", condition="control",
                              replicate_unit="mouse", eligible_primary=True, files="SIM"))
        ds = Dataset(sparse.vstack(Xs), pd.concat(obs), pd.DataFrame(index=genes)); ds.save(PROC / f"{c}_raw")
    pd.DataFrame(sheet).to_csv(DATA / "sample_sheet_sample.csv", index=False)
    (DATA / "cohorts_sample.json").write_text(json.dumps({"mouse": list(cohorts), "human_sn": []}))
    # ---------------- human / genetics / perturbation side files
    sim = RAW / "sim"; sim.mkdir(parents=True, exist_ok=True)
    hm = [g.upper() for g in genes]
    rows = [["DB Class Key", "Common Organism Name", "NCBI Taxon ID", "Symbol", "EntrezGene ID"]]
    for i, g in enumerate(genes): rows += [[i, "mouse, laboratory", 10090, g, 10000 + i], [i, "human", 9606, g.upper(), 20000 + i]]
    pd.DataFrame(rows[1:], columns=rows[0]).to_csv(sim / "HOM_MouseHumanSequence_sim.rpt", sep="\t", index=False)
    n_s = 80; age = rng.uniform(30, 95, n_s); rho = np.clip(1.15 - 0.02 * (age - 30), 0.02, 0.95)
    m = rng.normal(size=n_s); l = rho * m + np.sqrt(1 - rho ** 2) * rng.normal(size=n_s)
    hg = [g.upper() for g in G.MITOPHAGY + G.LYSOSOME] + [f"GENE{i:04d}" for i in range(300)]
    eta = rng.normal(5, 1, len(hg))[:, None] + np.zeros((1, n_s))
    for i, g in enumerate(hg):
        if g.title() in G.MITOPHAGY or g in [x.upper() for x in G.MITOPHAGY]: eta[i] += 0.6 * m
        if g in [x.upper() for x in G.LYSOSOME]: eta[i] += 0.6 * l
    for i, g in enumerate(hg):
        if g in [x.upper() for x in G.MITOPHAGY + G.LYSOSOME]: eta[i] += -0.012 * (age - 60)
    lib = np.exp(rng.normal(0, 0.25, n_s)); cnt = rng.poisson(np.exp(eta) * lib[None, :] * 3).astype(int)
    pd.DataFrame(cnt, index=hg, columns=[f"H{i+1:02d}" for i in range(n_s)]).to_csv(sim / "human_bulk_counts.csv.gz")
    pd.DataFrame({"sample": [f"H{i+1:02d}" for i in range(n_s)], "age": age.round(1),
                  "sex": rng.choice(["M", "F"], n_s)}).to_csv(sim / "human_bulk_meta.csv", index=False)
    z = rng.normal(0, 1, len(genes)); truth_ifc = [g.upper() for g in INTERFACE_MITO + INTERFACE_LYSO]
    for i, g in enumerate(genes):
        if g.upper() in truth_ifc: z[i] += 2.2
    gl = pd.DataFrame({"GENE": 20000 + np.arange(len(genes)), "CHR": rng.integers(1, 23, len(genes)),
                       "START": rng.integers(1e6, 2e8, len(genes)), "ZSTAT": z, "SYMBOL": hm})
    gl["P"] = 1 - __import__("scipy.stats").stats.norm.cdf(gl.ZSTAT); gl.to_csv(sim / "ad_gene_level_sim.csv", index=False)
    z2 = rng.normal(0, 1, len(genes)); g2 = gl.copy(); g2["ZSTAT"] = z2; g2["P"] = 1 - __import__("scipy.stats").stats.norm.cdf(z2)
    g2.to_csv(sim / "height_gene_level_sim.csv", index=False)
    pd.DataFrame({"gene": [h for h in hm if not h.startswith("GENE")] + [h for h in hm if h.startswith("GENE")][:700]}
                 ).to_csv(sim / "human_microglia_expressed_sim.csv", index=False)
    lf = rng.normal(0, 0.15, len(genes))
    for g in G.MITOPHAGY + G.LYSOSOME: lf[gi[g]] = 2.5 * beta[gi[g]] + rng.normal(0, 0.06)
    pd.DataFrame({"Gene": genes, "logFC": lf, "FDR": np.clip(np.abs(rng.normal(0.05, 0.1, len(genes))), 1e-4, 1)}
                 ).to_csv(sim / "cccp_sim.csv.gz", index=False)
    for shared in (True, False): make_coloc_locus(rng, shared).to_csv(sim / f"coloc_locus_{'shared' if shared else 'distinct'}.csv", index=False)
    json.dump({"interface_mito": INTERFACE_MITO, "interface_lyso": INTERFACE_LYSO}, open(sim / "truth.json", "w"))
    return sheet


def make_coloc_locus(rng, shared=True, n_snp=300, n_ind=3000):
    """Two association studies (GWAS-like and eQTL-like) at one locus with correlated SNPs."""
    maf = rng.uniform(0.1, 0.5, n_snp); g = np.zeros((n_ind, n_snp)); z = rng.normal(size=(n_ind, n_snp)); lat = np.zeros_like(z)
    lat[:, 0] = z[:, 0]
    for j in range(1, n_snp): lat[:, j] = 0.92 * lat[:, j - 1] + np.sqrt(1 - 0.92 ** 2) * z[:, j]
    thr = np.array([np.quantile(lat[:, j], 1 - maf[j]) for j in range(n_snp)]); g = (lat > thr).astype(float) + (lat > thr + 0.8)
    c1 = 120; c2 = 120 if shared else 230
    def study(c):
        y = 0.35 * g[:, c] + rng.normal(size=n_ind); rows = []
        for j in range(n_snp):
            x = g[:, j] - g[:, j].mean(); b = (x @ (y - y.mean())) / (x @ x); res = y - y.mean() - b * x
            rows.append((b, np.sqrt(res @ res / (n_ind - 2) / (x @ x))))
        return np.array(rows)
    a, b = study(c1), study(c2)
    return pd.DataFrame({"snp": [f"rs{i}" for i in range(n_snp)], "maf": maf / 2, "beta_gwas": a[:, 0], "se_gwas": a[:, 1],
                         "beta_eqtl": b[:, 0], "se_eqtl": b[:, 1], "n": n_ind})
