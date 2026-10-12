"""All figure functions. Each returns the matplotlib figure and saves a PNG into figures/<profile>/."""
import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
from .config import FIG

NAVY, TEAL, ORANGE, RED, GREY, PURPLE = "#1F3A5F", "#1B7F79", "#E07A1F", "#B83232", "#8896A8", "#6B4BA1"
STATE_COL = {"homeostatic": TEAL, "DAM": ORANGE, "interferon": PURPLE, "proliferating": "#C9A227", "chemokine_inflam": RED, "MHC_II": GREY}
plt.rcParams.update({"font.size": 8, "axes.spines.top": False, "axes.spines.right": False, "figure.dpi": 110,
                     "axes.titlesize": 8.5, "axes.titleweight": "bold", "axes.titlelocation": "left"})


def _save(fig, name):
    fig.tight_layout(); fig.savefig(FIG / name, dpi=170); return fig


def fig_qc(obs_before, obs_after, name="fig_qc.png"):
    fig, ax = plt.subplots(1, 3, figsize=(9.5, 2.8))
    for a, col, t in zip(ax, ["n_genes", "pct_mt", "total_counts"], ["genes per cell", "% mitochondrial reads", "UMI per cell"]):
        ss = sorted(obs_before.sample_id.unique())
        a.boxplot([obs_before[obs_before.sample_id == s][col] for s in ss], positions=np.arange(len(ss)) - .18, widths=.3,
                  patch_artist=True, boxprops=dict(facecolor="#D8DEE8"), medianprops=dict(color="k"), showfliers=False)
        a.boxplot([obs_after[obs_after.sample_id == s][col] for s in ss], positions=np.arange(len(ss)) + .18, widths=.3,
                  patch_artist=True, boxprops=dict(facecolor=TEAL, alpha=.8), medianprops=dict(color="k"), showfliers=False)
        a.set_xticks(range(len(ss))); a.set_xticklabels([s.split("_")[-1] for s in ss], rotation=90, fontsize=6); a.set_title(t)
    ax[0].text(0.02, 0.97, "grey = before QC, teal = after", transform=ax[0].transAxes, fontsize=6.5, va="top")
    return _save(fig, name)


def fig_clusters(Z, tbl, name="fig_states.png"):
    fig, ax = plt.subplots(1, 2, figsize=(9, 3.0), gridspec_kw={"width_ratios": [1, 1.2]})
    t = tbl[["micro", "nonmicro"]]; im = ax[0].imshow(t.values, cmap="RdBu_r", vmin=-max(abs(t.values.max()), 1), vmax=max(abs(t.values.max()), 1), aspect="auto")
    ax[0].set_yticks(range(len(t))); ax[0].set_yticklabels([f"cl {i} (n={n})" for i, n in zip(t.index, tbl.n_cells)], fontsize=6.5)
    ax[0].set_xticks([0, 1]); ax[0].set_xticklabels(["microglia\nmarkers", "non-microglia\nmarkers"]); ax[0].set_title("a  identity score by cluster")
    for i in range(len(t)):
        for j in range(2): ax[0].text(j, i, f"{t.values[i, j]:.2f}", ha="center", va="center", fontsize=6.5)
    ax[1].imshow(Z.values, cmap="viridis", aspect="auto"); ax[1].set_xticks(range(Z.shape[1])); ax[1].set_xticklabels(Z.columns, rotation=40, ha="right", fontsize=6.5)
    ax[1].set_yticks(range(len(Z))); ax[1].set_yticklabels(Z.index, fontsize=6.5); ax[1].set_title("b  marker-panel score in each called state")
    return _save(fig, name)


def fig_overview(emb, obs, name="fig_overview.png"):
    fig, ax = plt.subplots(1, 3, figsize=(10, 3.1))
    for a, col, t in zip(ax, ["cohort", "state", "age_months"], ["cohort", "state", "age (months)"]):
        v = obs[col].values
        if col == "age_months":
            sc_ = a.scatter(emb[:, 0], emb[:, 1], c=v.astype(float), s=2, cmap="viridis"); plt.colorbar(sc_, ax=a, shrink=.7)
        else:
            for k, lv in enumerate(sorted(set(v))):
                m = v == lv; a.scatter(emb[m, 0], emb[m, 1], s=2, color=STATE_COL.get(lv, [TEAL, ORANGE, PURPLE][k % 3]), label=lv)
            a.legend(markerscale=4, frameon=False, fontsize=6.5, loc="best")
        a.set_title(t); a.set_xticks([]); a.set_yticks([])
    return _save(fig, name)


def fig_mlci(R, meta, null, obs_slope, spec, name="fig2_mlci.png"):
    fig, ax = plt.subplots(1, 4, figsize=(13, 3.1)); cols = [TEAL, ORANGE, PURPLE, NAVY]
    d = R[R.level == "M2"]
    for k, (c, g) in enumerate(d.groupby("cohort")):
        ax[0].errorbar(g.age_months, g.r, yerr=[g.r - g.ci_lo, g.ci_hi - g.r], fmt="o", ms=3.5, color=cols[k], lw=.7, label=c)
        b = np.polyfit(g.log2_age, g.z, 1); xs = np.linspace(g.log2_age.min(), g.log2_age.max(), 20)
        ax[0].plot(2 ** xs, np.tanh(np.polyval(b, xs)), color=cols[k], lw=1)
    ax[0].set_xlabel("age (months)"); ax[0].set_ylabel("MLCI (partial Spearman r)"); ax[0].legend(frameon=False, fontsize=6.5); ax[0].set_title("a  one point = one mouse")
    s = meta["slopes"]; s = s[s.level == "M2"]
    for i, (_, r) in enumerate(s.iterrows()): ax[1].errorbar(r.beta, i, xerr=1.96 * r.se, fmt="s", color=GREY, ms=4)
    m = meta["pooled"]; ax[1].errorbar(m["beta"], len(s), xerr=1.96 * m["se"], fmt="D", color=RED, ms=6)
    ax[1].set_yticks(range(len(s) + 1)); ax[1].set_yticklabels(list(s.cohort) + ["pooled"]); ax[1].axvline(0, color="k", lw=.6)
    ax[1].set_xlabel("slope of z per doubling of age"); ax[1].set_title("b  forest plot (M2)")
    ax[2].hist(null, bins=30, color="#C9D1DC"); ax[2].axvline(obs_slope, color=RED); ax[2].set_xlabel("slope: random gene-set pairs"); ax[2].set_ylabel("count")
    ax[2].set_title("c  observed vs matched null")
    sp = spec.groupby("pair").beta.mean().sort_values(); ax[3].barh(sp.index, sp.values, color=[RED if "mitophagy-lysosome" == i else GREY for i in sp.index])
    ax[3].axvline(0, color="k", lw=.6); ax[3].set_xlabel("mean slope of z"); ax[3].set_title("d  specificity controls")
    return _save(fig, name)


def fig_state_res(T, comp, name="fig3_states.png"):
    fig, ax = plt.subplots(1, 2, figsize=(9, 3.1))
    for k, (st, g) in enumerate(T.groupby("state")):
        a = g.groupby("age_months").z.mean(); ax[0].plot(a.index, a.values, "o-", ms=3.5, color=STATE_COL.get(st, GREY), label=st)
    ax[0].set_xlabel("age (months)"); ax[0].set_ylabel("mean MLCI (z) within state"); ax[0].legend(frameon=False, fontsize=6.5); ax[0].set_title("a  decoupling inside each state")
    comp.groupby("age_months").mean().plot.bar(stacked=True, ax=ax[1], color=[STATE_COL.get(c, GREY) for c in comp.columns if c != "age_months"], width=.8, legend=False)
    ax[1].set_ylabel("fraction of cells"); ax[1].set_title("b  state composition changes too"); ax[1].legend(frameon=False, fontsize=6.5, bbox_to_anchor=(1, 1))
    return _save(fig, name)


def fig_rewire(Mj, Mo, A, B, tbl, name="fig4_rewiring.png"):
    fig, ax = plt.subplots(1, 3, figsize=(11, 3.4)); v = max(np.abs(Mj).max(), np.abs(Mo).max())
    for a, M, t in zip(ax[:2], [Mj, Mo], ["a  young: cross-arm link (z)", "b  old"]):
        im = a.imshow(M, cmap="RdBu_r", vmin=-v, vmax=v, aspect="auto"); a.set_title(t)
        a.set_xticks(range(len(B))); a.set_xticklabels(B, rotation=90, fontsize=4.5); a.set_yticks(range(len(A))); a.set_yticklabels(A, fontsize=5)
    top = tbl.head(12).iloc[::-1]; ax[2].barh(top.gene, top.delta, color=[ORANGE if x == "mitophagy" else TEAL for x in top.arm])
    ax[2].set_xlabel("weaker link in old (delta z)"); ax[2].set_title("c  top interface genes")
    return _save(fig, name)


def fig_de(sig, highlight, name="fig_age_de.png"):
    fig, ax = plt.subplots(1, 2, figsize=(8.5, 3.1))
    zc = [c for c in sig.columns if c.startswith("z_")]
    ax[0].scatter(sig[zc[0]], sig[zc[1]], s=3, color=GREY, alpha=.5)
    h = sig.loc[sig.index.intersection(highlight)]; ax[0].scatter(h[zc[0]], h[zc[1]], s=9, color=RED)
    ax[0].axhline(0, lw=.5, color="k"); ax[0].axvline(0, lw=.5, color="k"); ax[0].set_xlabel(zc[0]); ax[0].set_ylabel(zc[1]); ax[0].set_title("a  aging z per cohort (red = coupling genes)")
    ax[1].scatter(sig.meta_z, -np.log10(sig.p), s=3, color=GREY); ax[1].scatter(h.meta_z, -np.log10(h.p), s=9, color=RED)
    ax[1].set_xlabel("meta z (old vs young)"); ax[1].set_ylabel("-log10 p"); ax[1].set_title("b  replicated aging signature")
    return _save(fig, name)


def fig_tf(TF, tfs, name="fig_tf.png"):
    fig, ax = plt.subplots(1, 2, figsize=(8.5, 3.1))
    t = tfs.sort_values("p_mlci").head(8).iloc[::-1]; ax[0].barh(t.tf, t.beta_mlci, color=TEAL); ax[0].set_xlabel("effect on MLCI per unit TF activity"); ax[0].set_title("a  TF activity vs MLCI (age-adjusted)")
    tf = tfs.sort_values("p_mlci").tf.iloc[0]; g = TF
    for c, d in g.groupby("cohort"): ax[1].scatter(d[tf], d.mlci, s=16, label=c)
    ax[1].set_xlabel(f"{tf} regulon activity (sample mean)"); ax[1].set_ylabel("MLCI (z)"); ax[1].legend(frameon=False, fontsize=6.5); ax[1].set_title("b  best candidate regulator")
    return _save(fig, name)


def fig_human(d, f, p, name="fig5_human.png"):
    fig, ax = plt.subplots(1, 2, figsize=(8.5, 3.1))
    sc_ = ax[0].scatter(d.mito, d.lyso, c=d.age, cmap="viridis", s=18); plt.colorbar(sc_, ax=ax[0], label="age (y)")
    ax[0].set_xlabel("mitophagy score"); ax[0].set_ylabel("lysosome score"); ax[0].set_title("a  human microglia (bulk)")
    lo, hi = d.age <= d.age.median(), d.age > d.age.median()
    for m, c, lab in [(lo, TEAL, "younger half"), (hi, ORANGE, "older half")]:
        b = np.polyfit(d.mito[m], d.lyso[m], 1); xs = np.linspace(d.mito.min(), d.mito.max(), 10); ax[0].plot(xs, np.polyval(b, xs), color=c, label=f"{lab} slope {b[0]:.2f}")
    ax[0].legend(frameon=False, fontsize=6.5)
    ax[1].bar(["mito", "age", "mito:age"], f.loc[["mito", "age", "mito:age"], "t"], color=[GREY, GREY, RED]); ax[1].axhline(0, color="k", lw=.6)
    ax[1].set_ylabel("t statistic"); ax[1].set_title(f"b  model terms (interaction perm p = {p:.4f})")
    return _save(fig, name)


def fig_genetics(null, stat, coloc, name="fig6_genetics.png"):
    fig, ax = plt.subplots(1, 2, figsize=(8.5, 3.0))
    ax[0].hist(null, bins=40, color="#C9D1DC"); ax[0].axvline(stat["mean_z"], color=RED); ax[0].set_xlabel("mean gene-level AD Z of random expressed sets")
    ax[0].set_title(f"a  interface genes (n={stat['n_genes']}) p = {stat['p']:.4f}")
    c = coloc.set_index("locus")[["PP.H3", "PP.H4"]]; c.plot.bar(ax=ax[1], color=[GREY, RED], width=.7); ax[1].set_ylabel("posterior probability")
    ax[1].set_title("b  colocalisation with microglia eQTL"); ax[1].set_xticklabels(c.index, rotation=0); ax[1].legend(frameon=False, fontsize=6.5)
    return _save(fig, name)


def fig_perturb(J, rho, p, null, name="fig_perturbation.png"):
    fig, ax = plt.subplots(1, 2, figsize=(8.5, 3.0))
    ax[0].scatter(J.age_z, J.pert_lfc, s=14, color=NAVY); ax[0].axhline(0, lw=.5, color="k"); ax[0].axvline(0, lw=.5, color="k")
    ax[0].set_xlabel("aging signature (meta z)"); ax[0].set_ylabel("CCCP log fold-change"); ax[0].set_title(f"a  coupling genes: rho = {rho:.2f}")
    ax[1].hist(null, bins=30, color="#C9D1DC"); ax[1].axvline(rho, color=RED); ax[1].set_xlabel("rho in random same-size gene sets"); ax[1].set_title(f"b  empirical p = {p:.4f}")
    return _save(fig, name)


def fig_candidates(E, name="fig_candidates.png"):
    t = E.head(12).iloc[::-1]; fig, ax = plt.subplots(1, 2, figsize=(9, 3.4), gridspec_kw={"width_ratios": [1, 1.4]})
    ax[0].barh(t.index, t.rank_product, color=TEAL); ax[0].set_xlabel("rank product (lower = stronger)"); ax[0].set_title("a  top candidate genes")
    cols = [c for c in ["abs_delta", "abs_meta_z", "tf_target", "magma_z", "pp_h4", "pert_concord", "human_concord"] if c in t]
    M = t[cols].rank(pct=True).values; ax[1].imshow(M, cmap="Blues", aspect="auto"); ax[1].set_xticks(range(len(cols))); ax[1].set_xticklabels(cols, rotation=40, ha="right", fontsize=6.5)
    ax[1].set_yticks(range(len(t))); ax[1].set_yticklabels(t.index, fontsize=6.5); ax[1].set_title("b  evidence layers (darker = stronger)")
    return _save(fig, name)


def fig_robust(rob, name="fig_robustness.png"):
    fig, ax = plt.subplots(figsize=(5.5, 2.8))
    for i, (_, r) in enumerate(rob.iterrows()): ax.errorbar(r.beta, i, xerr=[[r.beta - r.ci_lo], [r.ci_hi - r.beta]], fmt="o", color=NAVY, ms=4)
    ax.set_yticks(range(len(rob))); ax.set_yticklabels(rob.variant, fontsize=7); ax.axvline(0, color="k", lw=.6)
    ax.set_xlabel("pooled slope of MLCI (z per doubling of age)"); ax.set_title("Robustness: every variant keeps the same sign")
    return _save(fig, name)
