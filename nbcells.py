"""Notebook definitions: every cell = code + plain-language explanation. Used to build .ipynb files, run the sample data and write the PDF."""
NBS = []
HEADER = '''import sys; sys.path.insert(0, "../src")
import numpy as np, pandas as pd, matplotlib.pyplot as plt
from mlc.config import *
from mlc import genesets as G
pd.set_option("display.width", 150); pd.set_option("display.max_columns", 20)
print("profile:", PROFILE, "| cohorts:", cohorts())'''
HEADER_WHAT = ("Loads the tools and the shared settings. <b>sys.path.insert</b> tells Python where the project code lives; "
               "<b>mlc.config</b> gives every notebook the same folders (PROC = processed data, RES = results tables, FIG = figures) and the profile "
               "(<i>sample</i> = built-in simulated data, <i>real</i> = your downloads). It also tells pandas to print wide tables on one screen. This cell is identical in every notebook.")
HEADER_FAIL = "ModuleNotFoundError: mlc  ->  the notebook is not inside the notebooks/ folder, or src/mlc is missing. Open the notebook from inside notebooks/."


def nb(id, title, goal, real_only=False, intro=""):
    n = dict(id=id, title=title, goal=goal, real_only=real_only, intro=intro, cells=[])
    NBS.append(n); return n


def cell(n, code, what, expect="", fail="", fig=None, run=True, header=False):
    n["cells"].append(dict(code=code.strip("\n"), what=what, expect=expect, fail=fail, fig=fig or [], run=run and not n["real_only"], header=header))


def start(n): cell(n, HEADER, HEADER_WHAT, fail=HEADER_FAIL, header=True)


# ================================================================== 00
n = nb("00", "Set up, simulated sample data, freeze the gene sets",
       "Create a realistic <b>simulated</b> dataset (2 cohorts, 19 mice, about 8,500 cells) with a known answer, so you can run the entire pipeline before touching real data; then freeze the gene sets.",
       intro="Why simulate? Real downloads are large and their file formats vary. The simulated data have the same structure as real data, so if every notebook runs here, "
             "any later problem is about your files, not your installation. The simulation contains: (1) coupling that falls with age, (2) more disease-associated (DAM) cells with age, "
             "(3) a few contaminating non-microglia cells, low-quality cells and doublets, (4) six 'interface' genes that carry the coupling.")
start(n)
cell(n, '''from mlc import synth
synth.make_sample_data()          # about 10 seconds
set_profile("sample")''',
     "<b>synth.make_sample_data()</b> writes two cohorts (SIM_A: 11 mice at 3, 12, 24 months; SIM_B: 8 mice at 4 and 22 months), a sample sheet, human bulk data, a gene-level GWAS table, "
     "a CCCP table and two colocalisation loci into data/. <b>set_profile</b> records that you are in sample mode.",
     expect="profile set to sample - restart this notebook's kernel before continuing",
     fail="PermissionError -> the project folder is read-only (for example inside a zip). Copy the folder to Documents first.")
cell(n, '''sheet = pd.read_csv(DATA / "sample_sheet_sample.csv")
print(sheet.head(6).to_string(index=False))
print(sheet.groupby(["cohort", "age_months"]).size())''',
     "The <b>sample sheet</b> has one row per biological sample (mouse). This is the most important table in the project: every statistic later uses the mouse, not the cell, as the unit.",
     expect="Prints 6 rows of the sheet and a count of mice per cohort and age (4+3+4 for SIM_A, 4+4 for SIM_B).")
cell(n, '''from mlc.dataset import Dataset
ds = Dataset.load(PROC / "SIM_A_raw")
print(ds)
print("median UMI per cell:", int(np.median(np.asarray(ds.X.sum(axis=1)).ravel())))
print("cells per sample:"); print(ds.obs.sample_id.value_counts().sort_index().head(4))''',
     "Loads one simulated cohort. A <b>Dataset</b> holds a sparse count matrix (cells x genes), a table about cells (obs) and a table about genes (var). "
     "<b>Note:</b> columns that start with <i>true_</i> exist only in simulated data (they store the answer); real data never has them and the analysis never uses them.",
     expect="4950 cells x 1500 genes; median UMI per cell about 1,600; 450 cells per sample.")
cell(n, '''import hashlib, json, datetime
sets = {k: getattr(G, k) for k in dir(G) if k.isupper()}
blob = json.dumps(sets, sort_keys=True).encode()
digest = hashlib.sha256(blob).hexdigest()
(RES / "genesets_frozen.json").write_text(json.dumps(sets, indent=1))
(RES / "genesets_frozen.sha256").write_text(f"{digest}  {datetime.date.today()}\\n")
print(digest)''',
     "<b>Freezing.</b> Writes every gene list into one file and computes its SHA-256 fingerprint (any change to any gene changes the fingerprint). "
     "Paste this fingerprint into your OSF pre-registration. It proves the gene sets were fixed <i>before</i> you looked at results.",
     expect="A 64-character hexadecimal string.")

# ================================================================== 01 (real only)
n = nb("01", "Real data: download files and build the sample sheet",
       "Download the public GEO files and create <b>data/sample_sheet.csv</b>. <b>Skip this notebook in the dry run</b> (it needs internet and your choices).", real_only=True,
       intro="This notebook is run once per cohort. Everything downloaded goes to data/raw/ and is never uploaded to GitHub.")
start(n)
cell(n, '''set_profile("real")      # then: Kernel > Restart, and re-run the first cell''',
     "Switches the project from simulated to real data. All real outputs go to separate folders (data/processed/real, results/real), so simulated results never mix with real ones.", run=False)
cell(n, '''from mlc.dataset import geo_list, geo_download, geo_sample_table
url, files = geo_list("GSE121654", "suppl")      # also try GSE129788
for f in files: print(f)''',
     "Asks NCBI GEO which supplementary files the series offers. This removes the need to hunt for download buttons.",
     expect="A list of file names such as GSE121654_RAW.tar.", fail="URLError -> no internet or GEO is busy; retry later. The same files can be downloaded by hand from the GEO page.", run=False)
cell(n, '''name = "GSE121654_RAW.tar"                     # EDIT: choose from the list above
geo_download(url + name, RAW / "GSE121654" / name)''',
     "Downloads one file to data/raw/GSE121654/. If it already exists it is skipped, so the cell is safe to re-run. Large files can take several minutes.", run=False)
cell(n, '''import tarfile
tarfile.open(RAW / "GSE121654" / name).extractall(RAW / "GSE121654")
from mlc.dataset import inventory
inventory(RAW / "GSE121654")''',
     "Unpacks the archive and prints every file with its size and what kind of file it looks like (10x matrix, genes, barcodes, dense table, h5ad). "
     "Use <b>inventory</b> whenever a loader complains: it shows what you actually have.", run=False)
cell(n, '''for acc in ["GSE121654", "GSE129788"]:
    t = geo_sample_table(acc, RAW / "geo_meta")
    t.to_csv(RAW / f"{acc}_geo_samples.csv", index=False)
    print(acc, t.shape); print(t.head(3).T)''',
     "Downloads each series' sample table (title, age, sex, tissue as written by the authors) and saves it as CSV. Open the CSV in Excel; you will copy its information into the sample sheet.", run=False)
cell(n, '''template = pd.DataFrame(columns=["cohort", "sample_id", "gsm", "species", "tissue", "region", "age_months", "sex",
                                 "condition", "replicate_unit", "eligible_primary", "files", "prefix"])
template.to_csv(DATA / "sample_sheet.csv", index=False)
print("Now open data/sample_sheet.csv in Excel and fill ONE ROW PER MOUSE/DONOR (see the table in the guide).")''',
     "Creates an empty sample sheet with the right column names. Fill one row per <b>independent biological sample</b>. <b>files</b> = folder (or file) holding that sample's counts; "
     "<b>prefix</b> = the file-name start when several samples share one folder (for example GSM5416205_6mth). <b>eligible_primary</b> = TRUE only for adult, control, independent microglia samples.", run=False)
cell(n, '''import json
from mlc.dataset import validate_sheet
ok = validate_sheet(pd.read_csv(SHEET))
(DATA / "cohorts_real.json").write_text(json.dumps({"mouse": sorted(ok.cohort.unique()), "human_sn": []}))''',
     "Checks the inclusion rule (at least 3 young and 3 old samples per cohort) and prints a plain-language verdict per cohort. Cohorts that fail are not used for the main tests. "
     "It also records which cohorts the next notebooks should process.", expect="e.g. 'GSE121654: 12 eligible samples; young=4 mid=0 old=8 -> OK'", run=False)

# ================================================================== 02 (real only)
n = nb("02", "Real data: load every sample into one object per cohort",
       "Read the downloaded count files (any common format) and save one raw dataset per cohort. <b>Skip in the dry run.</b>", real_only=True)
start(n)
cell(n, '''from mlc.dataset import read_any, concat, inventory
sheet = pd.read_csv(SHEET)
sheet = sheet[sheet.eligible_primary.astype(str).str.lower() == "true"]
inventory(sheet.files.iloc[0])           # sanity: what is in the first sample's folder?''',
     "Reads the sample sheet (eligible rows only) and shows the files of the first sample so you can confirm the paths are right before loading everything.", run=False)
cell(n, '''for c, sub in sheet.groupby("cohort"):
    parts = []
    for _, r in sub.iterrows():
        pre = r["prefix"] if "prefix" in sub.columns and isinstance(r["prefix"], str) else None
        d = read_any(r["files"], prefix=pre)
        d.obs.index = [f"{r.sample_id}_{i}" for i in d.obs.index]
        d.obs["sample_id"], d.obs["cohort"], d.obs["age_months"] = r.sample_id, c, float(r.age_months)
        parts.append(d); print(r.sample_id, d.n_obs, "cells")
    ds = concat(parts, join="outer"); ds.save(PROC / f"{c}_raw"); print(c, ds)''',
     "<b>read_any</b> detects the format (10x triplet, dense table, h5ad) automatically. Cell names are prefixed with the sample ID so two samples can never share a cell name. "
     "One combined dataset per cohort is saved as <i>{cohort}_raw</i>, exactly the name the simulated data used, so notebooks 03 onwards run unchanged.",
     fail="'No complete 10x triplet' -> run inventory(folder) and check file names; if a file is a single dense table pass its path in the 'files' column. "
          "'cell names must be unique' -> two rows of the sheet point to the same file.", run=False)

# ================================================================== 03
n = nb("03", "Quality control, sample by sample",
       "Remove low-quality cells using thresholds computed <b>inside each sample</b> (median absolute deviation), so one bad mouse cannot distort the others.")
start(n)
cell(n, '''from mlc.dataset import Dataset
from mlc import qc, plots
dsets = {}
for c in cohorts()["mouse"]:
    dsets[c] = qc.add_qc(Dataset.load(PROC / f"{c}_raw"), species="mouse")
    print(c, dsets[c].n_obs, "cells in", dsets[c].obs.sample_id.nunique(), "samples")
print(dsets[cohorts()["mouse"][0]].obs[["total_counts", "n_genes", "pct_mt"]].describe().round(1))''',
     "Loads each raw cohort and adds three quality numbers to every cell: <b>total_counts</b> (UMIs), <b>n_genes</b> (genes detected), <b>pct_mt</b> (percent of reads from mitochondrial genes; a high value means a dying cell). "
     "Use <i>species='human'</i> for human data (mitochondrial genes start with MT-).",
     expect="Counts of cells and a describe() table: median UMI about 1,600, median genes about 740.")
cell(n, '''qcd, logs = {}, []
for c, ds in dsets.items():
    f, log = qc.mad_filter(ds, nmads=5, mt_cap=15, min_genes=200)
    log["cohort"] = c; logs.append(log); qcd[c] = f
    print(c, "kept", f.n_obs, "of", ds.n_obs, f"({100 * f.n_obs / ds.n_obs:.1f}%)")
logs = pd.concat(logs); logs.to_csv(RES / "qc_log.csv", index=False); print(logs.head(5))''',
     "<b>mad_filter</b> removes, within each sample, cells whose UMI count or gene count is more than 5 median-absolute-deviations from that sample's median, cells with more than 15% mitochondrial reads, and cells with fewer than 200 genes. "
     "The log (cells removed per sample) is saved: reviewers like to see it. For single-nucleus data use <i>mt_cap=5</i>.",
     expect="About 98-99% kept per cohort (simulated data are cleaner than real data; real cohorts often keep 85-95%).")
cell(n, '''c0 = cohorts()["mouse"][0]
plots.fig_qc(dsets[c0].obs, qcd[c0].obs); plt.show()''',
     "Plots the three quality measures for every sample before (grey) and after (teal) filtering. After filtering the long low-quality tails should be gone while the bulk of cells is untouched. "
     "If a whole sample looks very different from the others, note it in the paper's supplement.",
     fig=["fig_qc.png"])
cell(n, '''for c, f in qcd.items():
    f.save(PROC / f"{c}_qc")
print("saved:", [p.name for p in PROC.glob("*_qc")])''',
     "Saves the filtered cohorts so the next notebook can continue from here. Every notebook saves a checkpoint, so you can stop and restart any time without repeating earlier steps.")

# ================================================================== 04
n = nb("04", "Keep only microglia and label their states",
       "Cluster cells, remove contaminating cell types and doublets using marker genes, then label each cluster as homeostatic, DAM, interferon, proliferating, and so on.",
       intro="States matter: DAM cells express more lysosomal genes, so a change in the <i>mix</i> of states could look like a change in coupling. "
             "Labelling states lets later notebooks separate the two.")
start(n)
cell(n, '''from mlc.dataset import Dataset
from mlc import proc, plots
mg = {}
for c in cohorts()["mouse"]:
    ds = proc.normalize(Dataset.load(PROC / f"{c}_qc"))
    proc.pca(ds, n_comps=30); proc.cluster(ds, k=8, seed=SEED)
    mg[c] = ds; print(c, ds.n_obs, "cells,", ds.obs.cluster.nunique(), "clusters")''',
     "<b>normalize</b>: divides each cell's counts by its total and multiplies by 10,000, then takes log(1+x), so deeply and shallowly sequenced cells become comparable. "
     "<b>pca</b>: picks the 2,000 most variable genes and compresses them into 30 axes. <b>cluster</b>: groups cells into 8 clusters (k-means on the first 20 axes).",
     expect="Each cohort: cells, 8 clusters.", fail="MemoryError -> too many cells for your RAM; set n_top=1000 in proc.pca or process cohorts one at a time.")
cell(n, '''for c, ds in mg.items():
    proc.add_scores(ds, {"micro": G.MICROGLIA, "nonmicro": G.NON_MICROGLIA})
    print(c); print(proc.cluster_table(ds))''',
     "Scores every cell for <b>microglia identity</b> (P2ry12, Tmem119, Cx3cr1 ...) and <b>non-microglia identity</b> (neurons, oligodendrocytes, astrocytes, endothelium, lymphocytes). "
     "A score is the average of the gene set minus the average of expression-matched control genes, so 0 means 'nothing special'. The table lists each cluster's mean scores; "
     "<b>contaminant = True</b> when the non-microglia score beats the microglia score.",
     expect="Mostly False; one or two small clusters are True (contaminants; scores about -1.9 and +1.3).")
cell(n, '''for c, ds in mg.items():
    t = proc.cluster_table(ds)
    drop = list(t.index[t.contaminant])                     # EDIT to override, e.g. drop = ["3"]
    nm = ds.obs.nonmicro; thr = nm.median() + 5 * (nm - nm.median()).abs().median() * 1.4826
    flagged = (nm > thr).values                              # cells with non-microglia genes switched on
    keep = ~ds.obs.cluster.isin(drop).values & ~flagged
    print(c, "dropped clusters:", drop, "| flagged cells:", int(flagged.sum()), "| kept", int(keep.sum()), "of", ds.n_obs)
    mg[c] = ds.subset_cells(keep)''',
     "Removes (1) whole contaminant clusters and (2) individual cells that express non-microglia genes unusually strongly "
     "(more than 5 robust standard deviations above the median), which also catches most doublets. "
     "This is a simple marker-based doublet filter; if you also want Scrublet, run it outside this guide and drop its flagged cells here.",
     expect="One or two small clusters dropped plus about 150-200 flagged cells per cohort (about 8% of cells removed in total).")
cell(n, '''for c, ds in mg.items():
    proc.normalize(ds); proc.pca(ds, n_comps=30); proc.cluster(ds, k=6, seed=SEED)
    ds.obs["dissoc"], _ = proc.score_genes(ds, G.DISSOCIATION)
    proc.add_scores(ds, G.STATE_PANELS, prefix="st_")
    Z = proc.label_states(ds, thr=1.0)
    print(c, ds.obs.state.value_counts().to_dict())''',
     "With contaminants gone, the cells are re-clustered (the clusters are used only for display). Each cell also receives a <b>dissociation-stress score</b> (Fos, Jun, Hspa1a ...; genes switched on by tissue dissociation, not by biology) "
     "and a score for every state's marker panel. A cell is called <b>DAM, interferon, proliferating, chemokine or MHC-II</b> if its best activated panel scores above <b>thr = 1.0</b> (about one log unit above matched controls); otherwise it is <b>homeostatic</b>. "
     "Raise thr for stricter calls; real data are noisier, so try 0.6 to 1.2 and check the heat map.",
     expect="e.g. {'homeostatic': 3258, 'DAM': 820, 'interferon': 250, 'proliferating': 231, 'MHC_II': 92, 'chemokine_inflam': 26}. Rare states have few cells and a few are false calls (noisy 4-gene panels), which is why later notebooks analyse only the large states.")
cell(n, '''c0 = cohorts()["mouse"][0]
plots.fig_clusters(Z, proc.cluster_table(mg[c0])); plt.show()
if "true_state" in mg[c0].obs:       # simulated data only: how well did we recover the truth?
    print(pd.crosstab(mg[c0].obs.state, mg[c0].obs.true_state))''',
     "Left: identity scores per cluster (all clusters should now be microglia). Right: mean score of every marker panel (columns) in every called state (rows): the diagonal should be bright. "
     "In the dry run the table compares our calls with the simulated truth; nearly every simulated DAM cell should be called DAM.",
     fig=["fig_states.png"])
cell(n, '''for c, ds in mg.items():
    ds.save(PROC / f"{c}_microglia"); print(c, ds.obs.state.value_counts().to_dict())''',
     "Saves the clean microglia with their state labels.")

# ================================================================== 05
n = nb("05", "Overview figure (Figure 1)",
       "Draw both cohorts in one map coloured by cohort, state and age. <b>For display only</b>: no statistic anywhere uses this map.")
start(n)
cell(n, '''from mlc.dataset import Dataset, concat
from mlc import proc, plots
parts = [Dataset.load(PROC / f"{c}_microglia") for c in cohorts()["mouse"]]
A = concat(parts, join="inner")
proc.normalize(A); proc.pca(A, n_comps=30, center_by="sample_id")
idx, emb = proc.display_embedding(A, max_cells=3000)
plots.fig_overview(emb, A.obs.iloc[idx]); plt.show()''',
     "Joins the cohorts (shared genes only), repeats normalisation and PCA with each sample's average removed (so batch offsets do not dominate the picture), and draws a 2-D t-SNE map of 3,000 random cells. "
     "<b>Important:</b> removing offsets can change correlations, which is why this map is never used for coupling statistics.",
     expect="Three panels. Cohorts should be intermixed; states form groups; age shows a gradient along the DAM direction.",
     fig=["fig_overview.png"], fail="t-SNE takes a minute for 3,000 cells; lower max_cells to 2,000 on slow computers.")

# ================================================================== 06
n = nb("06", "Check the gene sets and score every cell",
       "Confirm that the frozen gene sets are detected in your data, then compute one score per cell for each module (mitophagy, lysosome and four negative-control modules).")
start(n)
cell(n, '''from mlc.dataset import Dataset
from mlc import proc
modules = lambda ds: {"mitophagy": G.MITOPHAGY, "lysosome": G.LYSOSOME, "autophagy_core": G.AUTOPHAGY_CORE,
                      "proteasome": G.PROTEASOME, "ribosome": proc.prefix_genes(ds, G.RIBOSOME_PREFIX),
                      "oxphos": proc.prefix_genes(ds, G.OXPHOS_PREFIX)}
print("genes shared by mitophagy and lysosome sets:", set(G.MITOPHAGY) & set(G.LYSOSOME))
rows = []
for c in cohorts()["mouse"]:
    ds = Dataset.load(PROC / f"{c}_microglia")
    for name, genes in modules(ds).items():
        rows.append(dict(cohort=c, module=name, found=len(ds.gene_idx(genes)), total=len(genes)))
cov = pd.DataFrame(rows); cov["fraction"] = (cov.found / cov.total).round(2)
print(cov.pivot(index="module", columns="cohort", values="fraction"))''',
     "<b>Coverage check.</b> The two primary sets must not share genes (they do not). For each cohort the table shows what fraction of every gene set exists in the data. "
     "Below about 0.7 the module is unreliable in that cohort: say so in the Methods. Control modules are built from gene-name prefixes (Rpl/Rps = ribosome; Ndu/Sdh/Uqcr/Cox/Atp5 = oxidative phosphorylation).",
     expect="All 1.0 in the simulation (real data are usually 0.85-1.0).", fail="All fractions 0 -> gene symbols are in a different style (Ensembl IDs, or capital letters). Check ds.var_names[:10].")
cell(n, '''for c in cohorts()["mouse"]:
    ds = Dataset.load(PROC / f"{c}_microglia")
    proc.add_scores(ds, modules(ds), prefix="m_", seed=SEED)
    d = pd.get_dummies(ds.obs.state, prefix="state", drop_first=True).astype(float)
    for col in d: ds.obs[col] = d[col].values
    ds.save(PROC / f"{c}_scored")
    print(c, [x for x in ds.obs.columns if x.startswith("m_")])''',
     "Computes <b>m_mitophagy, m_lysosome, m_ribosome ...</b> for every cell (set average minus matched-control average) and adds 0/1 columns for the cell states (used later to adjust for state). "
     "One column per module, one row per cell, saved as <i>{cohort}_scored</i>.",
     expect="Six m_ columns per cohort.")
cell(n, '''ds = Dataset.load(PROC / f"{cohorts()['mouse'][0]}_scored")
o = ds.obs[ds.obs.sample_id == ds.obs.sample_id.iloc[0]]
print(o[[c for c in o.columns if c.startswith("m_")]].corr(method="spearman").round(2))''',
     "A first look at the raw relationships in one mouse: the correlation of every module with every other. The mitophagy-lysosome correlation here is <b>not yet adjusted</b> for sequencing depth or cell state; "
     "notebook 07 does that. Do not interpret these numbers.",
     expect="A 6 x 6 correlation table with ones on the diagonal.")

# ================================================================== 07
n = nb("07", "The core result: the Mitophagy-Lysosome Coupling Index (MLCI) across age",
       "One MLCI per mouse, a meta-analysis across cohorts, a random-gene-set null and four negative-control pairs. This notebook tests hypotheses H1 and H2.",
       intro="MLCI = partial Spearman correlation between the mitophagy score and the lysosome score <i>across the cells of one mouse</i>, after removing the effect of nuisance variables. "
             "Three levels: <b>M0</b> raw; <b>M1</b> adjusted for depth and mitochondrial %; <b>M2</b> (primary) also adjusted for stress score and cell state. "
             "Report all three: they show how much each confounder matters.")
start(n)
cell(n, '''from mlc.dataset import Dataset
from mlc.coupling import *
from mlc import plots
sheet = pd.read_csv(SHEET).set_index("sample_id")
obsd = {c: Dataset.load(PROC / f"{c}_scored").obs for c in cohorts()["mouse"]}
levels = {"M0": [], "M1": ["depth", "pct_mt"], "M2": ["depth", "pct_mt", "dissoc", "STATE"]}
rows = []
for c, o in obsd.items():
    state_cols = [x for x in o.columns if x.startswith("state_")]
    for lv, cov in levels.items():
        cv = [x for x in cov if x != "STATE"] + (state_cols if "STATE" in cov else [])
        r = coupling_per_sample(o, "m_mitophagy", "m_lysosome", covars=cv, min_cells=MIN_CELLS, n_boot=100, seed=SEED)
        r["level"], r["cohort"] = lv, c; rows.append(r)
R = pd.concat(rows, ignore_index=True)
R["age_months"] = R.sample_id.map(sheet.age_months).astype(float); R["log2_age"] = np.log2(R.age_months)
R.to_csv(RES / "mlci_per_sample.csv", index=False)
print(R[R.level == "M2"].groupby(["cohort", "age_months"]).r.mean().round(2))''',
     "<b>The central computation.</b> For every mouse: take its cells, regress the nuisance variables out of the ranked scores, correlate the residuals (r), convert to Fisher z (a scale suited to averaging), "
     "and bootstrap the cells 100 times for a 95% interval. <b>depth</b> = log total UMIs; <b>pct_mt</b> = mitochondrial %; <b>dissoc</b> = stress score; <b>state_*</b> = cell-state indicators.",
     expect="Average MLCI (r) falling with age: about 0.39, 0.32, 0.19 at 3, 12, 24 months (SIM_A) and 0.28 to 0.13 (SIM_B).",
     fail="'KeyError: dissoc' -> notebook 04 was not run for that cohort. 'Empty table' -> every sample has fewer than MIN_CELLS cells; lower MIN_CELLS in config.py.")
cell(n, '''S = cohort_slopes(R)
pooled = {lv: dersimonian_laird(g.beta, g.se) for lv, g in S.groupby("level")}
P_ = pd.DataFrame(pooled).T
print(S.round(3)); print(P_[["beta", "ci_lo", "ci_hi", "p", "I2", "k"]].round(4))
S.to_csv(RES / "mlci_cohort_slopes.csv", index=False); P_.to_csv(RES / "mlci_meta.csv")''',
     "<b>Age effect.</b> Within each cohort, a straight-line fit of MLCI (z) against log2(age) gives a slope (change in z each time age doubles) with its standard error. "
     "Slopes are then pooled across cohorts by a <b>random-effects meta-analysis</b>. A negative pooled slope with a confidence interval below 0 supports H1. "
     "<b>I2</b> = percentage of variation between cohorts that is real heterogeneity; with only 2-3 cohorts treat it as rough.",
     expect="Slope about -0.07 in each cohort; pooled M2 slope -0.067 (95% interval about -0.088 to -0.047), p far below 0.05. M0, M1 and M2 almost identical here; in real data they usually differ more.")
cell(n, '''from mlc import genesets as G
n_perm = 200                      # use 1000+ for real analyses (takes longer)
nulls, obs_slopes = {}, {}
for c in cohorts()["mouse"]:
    ds = Dataset.load(PROC / f"{c}_scored")
    keep = np.asarray((ds.L > 0).mean(axis=0)).ravel() > 0.05             # genes seen in >5% of cells
    ds = ds.subset_genes(keep); means = np.asarray(ds.L.mean(axis=0)).ravel()
    iA, iB = ds.gene_idx(G.MITOPHAGY), ds.gene_idx(G.LYSOSOME)
    ids = R[(R.cohort == c) & (R.level == "M0")].sample_id.tolist()
    Xs = {s: ds.L[(ds.obs.sample_id == s).values].toarray() for s in ids}
    ages = np.log2(sheet.loc[ids, "age_months"].astype(float).values); coh = np.array([c] * len(ids))
    ca, cb = matched_controls(means, iA), matched_controls(means, iB)
    z = [np.arctanh(partial_spearman(fast_score(Xs[s], iA, ca), fast_score(Xs[s], iB, cb))) for s in ids]
    obs_slopes[c] = age_slope(z, ages, coh)
    nulls[c] = perm_null_age_slope(Xs, means, iA, iB, ages, coh, n_perm=n_perm, seed=SEED)
    p = (1 + (np.abs(nulls[c]) >= abs(obs_slopes[c])).sum()) / (n_perm + 1)
    print(c, "observed slope", round(obs_slopes[c], 3), "| null sd", round(nulls[c].std(), 3),
          "| permutation p", round(p, 4))
pooled_null = np.mean(list(nulls.values()), axis=0); pooled_obs = np.mean(list(obs_slopes.values()))
p_pool = (1 + (np.abs(pooled_null) >= abs(pooled_obs)).sum()) / (n_perm + 1)
print("pooled observed", round(pooled_obs, 3), "permutation p", round(p_pool, 4))
pd.DataFrame({"null": pooled_null}).to_csv(RES / "perm_null_pooled.csv", index=False)''',
     "<b>The matched random null.</b> Is a decline this steep special to <i>these</i> gene sets, or would any two gene sets of the same size and expression level decline like this? "
     "Each permutation draws a random set for each arm (same size, same expression-level profile as your real sets), recomputes the MLCI in every mouse and the age slope. "
     "The permutation p is the fraction of random pairs with a slope at least as extreme as yours. The smallest possible p is 1/(n_perm+1) = 0.005 here.",
     expect="Observed slope about -0.07; null spread (sd) about 0.011-0.014; permutation p = 0.005 (the smallest value possible with 200 permutations).",
     fail="Slow? Lower n_perm to 100. 'MemoryError' -> a cohort has too many cells: keep only the first 1,500 cells per sample when building Xs.")
cell(n, '''pairs = [("m_mitophagy", "m_lysosome"), ("m_ribosome", "m_lysosome"), ("m_proteasome", "m_lysosome"),
         ("m_oxphos", "m_mitophagy"), ("m_autophagy_core", "m_lysosome")]
rows = []
for c, o in obsd.items():
    cv = ["depth", "pct_mt", "dissoc"] + [x for x in o.columns if x.startswith("state_")]
    for a, b in pairs:
        r = coupling_per_sample(o, a, b, covars=cv, min_cells=MIN_CELLS, n_boot=0)
        r["pair"], r["cohort"] = f"{a[2:]}-{b[2:]}", c; rows.append(r)
SP = pd.concat(rows); SP["log2_age"] = np.log2(SP.sample_id.map(sheet.age_months).astype(float))
spec = cohort_slopes(SP, by=("pair", "cohort")); spec.to_csv(RES / "specificity_slopes.csv", index=False)
print(spec.pivot(index="pair", columns="cohort", values="beta").round(3))''',
     "<b>Specificity (H2).</b> Repeats the age-slope for four negative-control pairs (ribosome-lysosome, proteasome-lysosome, oxidative-phosphorylation-mitophagy, core autophagy-lysosome). "
     "If the mitophagy-lysosome slope is clearly more negative than all of these, the effect is specific; if the control pairs decline just as much, the finding is a general loss of co-expression, not coupling.",
     expect="Mitophagy-lysosome row about -0.07; control rows between -0.02 and +0.02.")
cell(n, '''fig = plots.fig_mlci(R, {"slopes": S, "pooled": pooled["M2"]}, pooled_null, pooled_obs, spec)
plt.show()''',
     "<b>Figure 2.</b> (a) One point per mouse with its bootstrap interval and a fitted line per cohort. (b) Forest plot: cohort slopes and the pooled estimate (red) with 95% intervals. "
     "(c) Your observed slope (red line) against the matched-random-set null. (d) The mitophagy-lysosome pair (red) against the negative controls.",
     fig=["fig2_mlci.png"])

# ================================================================== 08
n = nb("08", "Is the decoupling inside cell states, or just a change in the state mix?",
       "Compute MLCI within each state separately and show how state proportions change with age (hypothesis H3).")
start(n)
cell(n, '''from mlc.dataset import Dataset
from mlc.coupling import *
from mlc import plots
sheet = pd.read_csv(SHEET).set_index("sample_id")
obsd = {c: Dataset.load(PROC / f"{c}_scored").obs for c in cohorts()["mouse"]}
rows = []
for c, o in obsd.items():
    for st, d in o.groupby("state"):
        r = coupling_per_sample(d, "m_mitophagy", "m_lysosome", covars=["depth", "pct_mt", "dissoc"],
                                min_cells=30, n_boot=0)
        r["state"], r["cohort"] = st, c; rows.append(r)
T = pd.concat(rows, ignore_index=True)
T["age_months"] = T.sample_id.map(sheet.age_months).astype(float); T["log2_age"] = np.log2(T.age_months)
print(T.groupby(["state", "cohort"]).size().unstack())''',
     "Splits each mouse's cells by state and computes the MLCI inside each state (state is now the grouping, so it is not adjusted for). "
     "States with fewer than 30 cells in a mouse are skipped (<b>min_cells=30</b>), so rare states may have too few mice to analyse.",
     expect="A table: number of mice per state and cohort (homeostatic and DAM have most).")
cell(n, '''S_state = cohort_slopes(T, by=("state", "cohort"))
res = {st: dersimonian_laird(g.beta, g.se) for st, g in S_state.groupby("state") if len(g) >= 2}
out = pd.DataFrame(res).T[["beta", "ci_lo", "ci_hi", "p", "k"]].round(3)
print(out); out.to_csv(RES / "mlci_state_meta.csv")''',
     "For each state with data from at least two cohorts, pools the age slope across cohorts. "
     "<b>Reading it:</b> a negative slope with an interval below 0 inside the homeostatic state means decoupling happens within cells of the same type, not only because more cells became DAM.",
     expect="Homeostatic: slope about -0.06 with an interval clearly below 0. DAM: negative but with a wide interval that includes 0 (few young DAM cells): report that honestly.")
cell(n, '''comp = (pd.concat([o[["sample_id", "state"]] for o in obsd.values()])
        .groupby("sample_id").state.value_counts(normalize=True).unstack(fill_value=0))
comp["age_months"] = comp.index.map(sheet.age_months).astype(float)
comp.to_csv(RES / "state_proportions.csv")
print(comp.groupby("age_months").mean().round(2))''',
     "State composition per mouse (fraction of its cells in each state), averaged by age. "
     "This shows what <i>does</i> change by composition (DAM rises with age), so you can state honestly that both effects exist.",
     expect="DAM fraction rising from about 5% to 25-30% with age.")
cell(n, '''plots.fig_state_res(T, comp); plt.show()''',
     "<b>Figure 3.</b> (a) Mean MLCI inside each state by age: lines decline in every large state. (b) Composition by age: DAM expands. Together: decoupling is a within-state property, in addition to a composition shift.",
     fig=["fig3_states.png"])

# ================================================================== 09
n = nb("09", "Which genes carry the decoupling? (gene-level rewiring)",
       "Find the mitophagy and lysosome genes whose link to the other arm weakens most with age: the <b>interface genes</b>. These become your candidate list.")
start(n)
cell(n, '''from mlc.dataset import Dataset
from mlc.analysis import cross_arm_z, rewiring_test
from mlc import plots
sheet = pd.read_csv(SHEET).set_index("sample_id")
Z, cohs = {}, {}
for c in cohorts()["mouse"]:
    ds = Dataset.load(PROC / f"{c}_scored")
    cv = ["depth", "pct_mt", "dissoc"] + [x for x in ds.obs.columns if x.startswith("state_")]
    z = cross_arm_z(ds, G.MITOPHAGY, G.LYSOSOME, cv, min_cells=MIN_CELLS)
    Z.update(z); cohs.update({s: c for s in z}); print(c, len(z), "samples")
ages = [float(sheet.loc[s, "age_months"]) for s in Z]; coh = [cohs[s] for s in Z]''',
     "For every mouse, computes the adjusted correlation between <i>each</i> mitophagy gene and <i>each</i> lysosome gene (19 x 44 pairs; Fisher-z scale). "
     "The result is one small matrix per mouse. Nothing is averaged over cells that belong to a different mouse.",
     expect="'SIM_A 11 samples' and 'SIM_B 8 samples'.")
cell(n, '''tbl, Mj, Mo, A, B = rewiring_test(Z, ages, coh, n_perm=1000, seed=SEED)
tbl.to_csv(RES / "rewiring_genes.csv", index=False)
print(tbl.head(12).round(4))''',
     "<b>delta</b> for a gene = its average link to the other arm in young mice minus in old mice (positive = weaker with age). "
     "Young/old are defined within each cohort (below or above the cohort's median age). "
     "The p-value comes from shuffling the young/old labels 1,000 times <i>within cohort</i> (so the replicated design is respected), then Benjamini-Hochberg <b>fdr</b> controls false discoveries across all genes. "
     "Use n_perm=5000 or more for real analyses.",
     expect="Top genes include Ctsd, Lamp2, Mcoln1, Atp6v0a1, Lamp1 (lysosome) and Pink1, Bnip3l (mitophagy); FDR about 0.02-0.03.",
     fail="'min() arg is an empty sequence' -> no mouse passed MIN_CELLS; check notebook 07 ran.")
cell(n, '''plots.fig_rewire(Mj, Mo, A, B, tbl); plt.show()''',
     "<b>Figure 4a-c.</b> Heat maps of the average cross-arm link in young and old mice (rows mitophagy genes, columns lysosome genes; red = positive link): the bright block fades with age. "
     "Panel c ranks genes by loss of coupling. Genes with a small fdr are your candidate list for notebooks 13 to 15.",
     fig=["fig4_rewiring.png"])
cell(n, '''import json
truth = json.load(open(RAW / "sim" / "truth.json"))
true_set = set(truth["interface_mito"] + truth["interface_lyso"])
top = set(tbl.head(14).gene)
print("simulated interface genes found in the top 14:", len(top & true_set), "of", len(true_set), sorted(top & true_set))''',
     "<b>Simulated data only:</b> checks the method against the known answer. The top 14 genes should contain most of the 14 planted interface genes. "
     "On real data there is no such check, which is why the permutation p-values and replication matter.",
     expect="11 of 14 planted genes recovered (about 10-13 is typical).")

# ================================================================== 10
n = nb("10", "A replicated aging signature (differential expression)",
       "Which genes change with age, using the mouse as the unit of replication: sum counts per mouse, fit a model per gene, combine cohorts.",
       intro="This replaces the descriptive fold-change ranking of the first guide. Cells from one mouse are added together first ('pseudobulk'), so the statistics see 19 mice, not 8,500 cells.")
start(n)
cell(n, '''from mlc.dataset import Dataset
from mlc.analysis import pseudobulk, age_de, stouffer
from mlc import plots
sheet = pd.read_csv(SHEET).set_index("sample_id")
zs = {}
for c in cohorts()["mouse"]:
    ds = Dataset.load(PROC / f"{c}_scored")
    pb = pseudobulk(ds, MIN_CELLS)
    age = np.log2(sheet.loc[pb.columns, "age_months"].astype(float).values)
    de = age_de(pb, age); de.to_csv(RES / f"age_de_{c}.csv"); zs[c] = de["z"]
    print(c, "mice:", pb.shape[1], "| genes tested:", len(de), "| p<0.05:", int((de.p < 0.05).sum()))''',
     "<b>pseudobulk</b> sums the counts of all cells of one mouse. <b>age_de</b> converts to log2 counts-per-million and fits, for every gene, a straight line against log2(age); the sign and size of the slope become a z-score. "
     "It keeps genes with at least 1 CPM in half the mice. This is a light-weight alternative to DESeq2; for the final paper you can repeat it with pyDESeq2 (optional) and check the top genes agree.",
     expect="Each cohort: 11 and 8 mice; all 1,500 genes tested; about 360 genes with p<0.05 (real data: 15,000+ genes tested).")
cell(n, '''sig = stouffer(zs)
sig.to_csv(RES / "age_signature_meta.csv")
print(sig[sig.concordant].sort_values("p").head(10).round(4))''',
     "<b>Stouffer's method</b> combines the cohorts' z-scores (sum divided by square root of the number of cohorts); "
     "<b>concordant</b> means every cohort moves in the same direction. Genes that are significant <i>and</i> concordant form the replicated aging signature.",
     expect="DAM genes up with age (Itgax, Lpl, Cst7, Cd68, Trem2, Apoe; meta z about 6-7) and planted interface genes down (for example Bnip3).")
cell(n, '''plots.fig_de(sig, G.MITOPHAGY + G.LYSOSOME); plt.show()''',
     "(a) The z-score of every gene in cohort 1 against cohort 2: genes along the diagonal replicate. Red = mitophagy and lysosome genes. (b) Volcano-style view of the combined signature.",
     fig=["fig_age_de.png"])

# ================================================================== 11
n = nb("11", "Which transcription factors track the coupling?",
       "Estimate transcription-factor activity from the expression of each factor's target genes, then ask which factor explains differences in MLCI between mice (H5, regulatory arm).")
start(n)
cell(n, '''from mlc.dataset import Dataset
from mlc.analysis import get_net, tf_activity
from mlc.stats import ols, design, bh_fdr
from mlc import plots
sheet = pd.read_csv(SHEET).set_index("sample_id")
net = get_net("mouse")
print(net.groupby("source").size().sort_values(ascending=False).head(5))''',
     "<b>get_net</b> loads a table of 'transcription factor -> target gene' relations. If you installed <i>decoupler</i> (pip install decoupler==1.8.0) it downloads CollecTRI, a curated network of hundreds of factors; "
     "otherwise it uses a built-in mini-network of the TFEB 'CLEAR' lysosomal-biogenesis targets so the notebook still runs.",
     expect="Either many TFs (CollecTRI) or a single TFEB entry with 21 targets (built-in).",
     fail="A warning that decoupler is not available is normal in the dry run, not an error.")
cell(n, '''rows = []
for c in cohorts()["mouse"]:
    ds = Dataset.load(PROC / f"{c}_scored")
    act = tf_activity(ds, net); act["sample_id"] = ds.obs.sample_id.values
    per = act.groupby("sample_id").mean(); per["cohort"] = c; rows.append(per)
TF = pd.concat(rows)
TF["log2_age"] = np.log2(sheet.loc[TF.index, "age_months"].astype(float).values)
R = pd.read_csv(RES / "mlci_per_sample.csv")
TF["mlci"] = R[R.level == "M2"].set_index("sample_id").z.reindex(TF.index)
TF = TF.dropna(subset=["mlci"]); print(TF.shape, "mice")''',
     "Per cell, the activity of a factor = the average (z-scored) expression of its target genes. Then averaged per mouse, so each mouse has one activity per factor next to its MLCI. "
     "<b>Caution:</b> activity is inferred from target expression, not measured; TFEB is controlled mostly by moving into the nucleus, which RNA data cannot see.",
     expect="(19, ...) = 19 mice.")
cell(n, '''tfs = [c for c in TF.columns if c not in ("cohort", "log2_age", "mlci")]; out = []
for t in tfs:
    d = TF.assign(tfv=TF[t]); X, nm = design(d, ["tfv", "log2_age"], cohort="cohort"); f = ols(d.mlci, X, nm)
    Xa, na = design(d, ["log2_age"], cohort="cohort"); fa = ols(d[t], Xa, na)
    out.append(dict(tf=t, beta_mlci=f.loc["tfv", "beta"], p_mlci=f.loc["tfv", "p"],
                    beta_age=fa.loc["log2_age", "beta"], p_age=fa.loc["log2_age", "p"]))
tfs_df = pd.DataFrame(out); tfs_df["fdr_mlci"] = bh_fdr(tfs_df.p_mlci.values)
tfs_df.sort_values("p_mlci").to_csv(RES / "tf_vs_mlci.csv", index=False)
print(tfs_df.sort_values("p_mlci").head(8).round(4))''',
     "For each factor, two small regressions across mice: (1) <b>MLCI ~ factor activity + age + cohort</b> (does activity explain coupling beyond age?) and (2) <b>activity ~ age + cohort</b> (does activity itself change with age?). "
     "<b>beta_mlci</b> positive = higher activity goes with stronger coupling. With hundreds of factors, use <b>fdr_mlci</b>, not the raw p.",
     expect="TFEB: beta_mlci about +0.09, p about 0.005; its age effect is negative but not significant (p about 0.2).")
cell(n, '''plots.fig_tf(TF, tfs_df); plt.show()''',
     "(a) Strength of association with MLCI for the top factors. (b) The best candidate, one point per mouse, coloured by cohort. If the points rise from left to right within each cohort, mice with higher regulon activity keep tighter coupling.",
     fig=["fig_tf.png"])

# ================================================================== 12
n = nb("12", "Human validation: bulk microglia, and what to do with single-nucleus data",
       "Test whether the mitophagy-lysosome relationship also weakens with age in people (H4). The dry run uses simulated human data in the same format as GSE99074.")
start(n)
cell(n, '''from mlc.genetics import load_orthologs, to_human
from mlc.analysis import bulk_scores, interaction_test
from mlc import plots
SIM = RAW / "sim"          # REAL DATA: replace with your downloaded files (see the guide)
orth = load_orthologs(SIM / "HOM_MouseHumanSequence_sim.rpt")      # real: RAW / "HOM_MouseHumanSequence.rpt"
H_MITO, H_LYSO = to_human(G.MITOPHAGY, orth), to_human(G.LYSOSOME, orth)
print(len(H_MITO), "human mitophagy genes;", len(H_LYSO), "human lysosome genes")
counts = pd.read_csv(SIM / "human_bulk_counts.csv.gz", index_col=0)
meta = pd.read_csv(SIM / "human_bulk_meta.csv").set_index("sample")
print(counts.shape, "genes x samples; ages", meta.age.min(), "to", meta.age.max())''',
     "Maps mouse gene symbols to human ones with the MGI orthology table (some genes map to several, some to none), and loads a human count table (genes x samples) with each person's age and sex. "
     "In real data (GSE99074) the count file name and the metadata column names differ: print <i>counts.columns[:5]</i> and the GEO sample table and adjust the two paths.",
     expect="19 and 44 human genes; (363, 80) genes x samples.", fail="Zero genes mapped -> the orthology file or symbol style is wrong (print orth.head()).")
cell(n, '''d, (na, nb_) = bulk_scores(counts, H_MITO, H_LYSO)
d["age"] = meta.age.reindex(d.index); d["sex_M"] = (meta.sex.reindex(d.index) == "M").astype(float)
f, p, null = interaction_test(d, covars=["sex_M"], n_perm=2000, seed=SEED)
print(f.round(4)); print("permutation p for mito:age interaction =", round(p, 4))''',
     "Each person gets a <b>mitophagy score</b> and a <b>lysosome score</b> (mean of z-scored log2-CPM of the genes). The model is <i>lyso ~ mito x age + sex</i>. "
     "The row <b>mito:age</b> asks whether the mitophagy-lysosome relationship changes with age; a <b>negative</b> t-statistic means it weakens, as in mouse. "
     "The permutation p shuffles ages 2,000 times so the p-value does not rely on normality.",
     expect="mito:age t about -2.8; permutation p about 0.03 (the main mito effect is strongly positive).")
cell(n, '''plots.fig_human(d, f, p); plt.show()''',
     "<b>Figure 5.</b> (a) Each dot is a person; lines show the mitophagy-lysosome slope in the younger and older half: the older half is flatter. (b) The model's t-statistics; the red bar is the interaction.",
     fig=["fig5_human.png"])
cell(n, '''from mlc.stats import ols_multi
lg = np.log2(counts / counts.sum() * 1e6 + 1); a_ = meta.age.reindex(lg.columns).values
X = np.column_stack([np.ones(len(a_)), (a_ - a_.mean()) / a_.std()])
B, SE, Pv = ols_multi(lg.T.values, X)
pd.DataFrame({"gene": lg.index, "beta_age": B[1], "p": Pv[1]}).to_csv(RES / "human_gene_age.csv", index=False)
print("saved human gene-level age associations:", lg.shape[0], "genes")''',
     "Also records, for every gene, how its expression changes with age in the human samples. Notebook 15 uses this as one evidence layer (does the mouse aging direction agree with human?).")
cell(n, '''print("Human single-nucleus cohorts (GSE148822, GSE157827): repeat notebooks 02 to 07 with species='human'")''',
     "<b>Human single-nucleus data.</b> Use the same notebooks: in 03 call <i>qc.add_qc(..., species='human')</i> with <i>mt_cap=5</i>; in 04 use human gene symbols (capital letters); in 06 use the human gene lists "
     "<b>H_MITO</b> and <b>H_LYSO</b> from above (and capitalised prefixes RPL/RPS, NDU/SDH/UQCR/COX/ATP5); each person is one sample and has few microglia, so lower MIN_CELLS to 50.", run=True)

# ================================================================== 13
n = nb("13", "Human genetics: are the interface genes enriched for Alzheimer's risk?",
       "Use gene-level Alzheimer's GWAS results to test the interface genes against random sets of microglia-expressed genes, and test colocalisation with a microglia eQTL (H5).",
       intro="Dry run: a simulated gene-level table in the same shape as MAGMA output. Real data: run MAGMA first (instructions after the cells) and point the path to <i>ad_gene.genes.out</i> joined with gene symbols.")
start(n)
cell(n, '''from mlc.genetics import load_orthologs, to_human, set_enrichment, coloc_abf
from mlc import plots
SIM = RAW / "sim"
orth = load_orthologs(SIM / "HOM_MouseHumanSequence_sim.rpt")
gt = pd.read_csv(SIM / "ad_gene_level_sim.csv")             # real: ad_gene.genes.out merged with symbols (see below)
expressed = pd.read_csv(SIM / "human_microglia_expressed_sim.csv").gene
rew = pd.read_csv(RES / "rewiring_genes.csv")
iface = rew[rew.fdr < 0.1].gene.tolist()
if len(iface) < 5: iface = rew.head(10).gene.tolist(); print("few genes at FDR<0.1; using top 10 by p-value")
hset = to_human(iface, orth); print(len(iface), "interface genes ->", len(hset), "human genes:", hset[:8])''',
     "Takes the interface genes from notebook 09 (FDR below 0.1), converts them to human symbols, and loads the gene-level Alzheimer's results (one Z-score per gene; larger = stronger association). "
     "<b>expressed</b> = genes expressed in human microglia: the comparison pool must be microglia-expressed genes, otherwise any microglia gene set looks 'enriched' just because microglia genes are enriched in Alzheimer's risk.",
     expect="21 interface genes mapped to 21 human symbols.")
cell(n, '''excl = ((gt.CHR == 19) & gt.START.between(44e6, 46.5e6)) | ((gt.CHR == 6) & gt.START.between(25e6, 34e6))   # APOE, MHC
stat, null = set_enrichment(gt, hset, expressed, n_perm=10000, seed=SEED, exclude=excl)
print("Alzheimer's:", {k: round(v, 4) for k, v in stat.items()})
ctrl = pd.read_csv(SIM / "height_gene_level_sim.csv")       # an unrelated trait as a specificity control
stat_c, null_c = set_enrichment(ctrl, hset, expressed, n_perm=10000, seed=SEED, exclude=excl)
print("control trait:", {k: round(v, 4) for k, v in stat_c.items()})''',
     "Compares the average gene-level Z of the interface genes with 10,000 random sets of the same size drawn from microglia-expressed genes (<b>empirical p</b>). "
     "The APOE and MHC regions are excluded because they dominate Alzheimer's signal. The same test on an <b>unrelated trait</b> (height) is the specificity control: it should be <i>not</i> significant.",
     expect="Alzheimer's: mean Z about 1.6 against a null of about 0.07, p = 0.0001; control trait p about 0.65.")
cell(n, '''rows = []
for locus, gene in [("shared", "Ctsd"), ("distinct", "Pink1")]:      # real: one row per candidate gene and its eQTL file
    d = pd.read_csv(SIM / f"coloc_locus_{locus}.csv")
    r = coloc_abf(d.beta_gwas, d.se_gwas, d.beta_eqtl, d.se_eqtl, "quant", "quant")
    rows.append(dict(locus=gene, gene=gene, **{k: round(v, 3) for k, v in r.items()}))
coloc = pd.DataFrame(rows); coloc.to_csv(RES / "coloc_results.csv", index=False); print(coloc)''',
     "<b>Colocalisation</b> asks: do the disease association and the microglia eQTL at a locus come from the <i>same</i> causal variant (PP.H4 near 1) or from two nearby different ones (PP.H3 near 1)? "
     "The function is a transparent implementation of the standard Bayesian method (approximate Bayes factors; default priors p1=p2=1e-4, p12=1e-5). "
     "Each row of the real analysis needs, for one locus, the same SNPs with beta and standard error from both studies.",
     expect="Ctsd locus: PP.H4 about 1.0; Pink1 locus: PP.H3 about 1.0.", fail="'Different SNP counts' -> merge the two tables on rsID first and keep common SNPs; check both use the same genome build.")
cell(n, '''plots.fig_genetics(null, stat, coloc); plt.show()''',
     "<b>Figure 6a-b.</b> (a) Real mean Z (red) against the random expressed-gene null. (b) Posterior probabilities for the two loci.",
     fig=["fig6_genetics.png"])

# ================================================================== 14
n = nb("14", "Do experimental perturbations resemble the aging signature?",
       "Compare the replicated aging signature with a mitochondrial-stress perturbation (CCCP), using only the coupling genes and a proper random null.")
start(n)
cell(n, '''from scipy.stats import spearmanr
from mlc import plots
cccp = pd.read_csv(RAW / "sim" / "cccp_sim.csv.gz")        # real: GSE179611_DMSOLPS_vs_CCCPLPS-ExpDiff.csv.gz
print(cccp.columns.tolist())
GENE, LFC = "Gene", "logFC"                                 # EDIT for real data after reading the printed column names
pert = cccp.set_index(GENE)[LFC]
age = pd.read_csv(RES / "age_signature_meta.csv", index_col=0)
set_genes = [g for g in set(G.MITOPHAGY + G.LYSOSOME) if g in age.index and g in pert.index]
J = pd.DataFrame({"age_z": age.loc[set_genes, "meta_z"], "pert_lfc": pert.loc[set_genes]})
rho, _ = spearmanr(J.age_z, J.pert_lfc); print(len(set_genes), "coupling genes; Spearman rho =", round(rho, 3))''',
     "Loads a perturbation table (CCCP versus control, log fold-change per gene). <b>GENE</b> and <b>LFC</b> are the two column names; in the real GEO file they differ (commonly 'gene', 'log2FoldChange' or 'logFC'), so read the printed names and edit. "
     "Then correlates, over the coupling genes only, the aging z-score with the perturbation log fold-change.",
     expect="63 coupling genes; rho about 0.26.", fail="KeyError 'Gene' -> edit GENE/LFC to match the printed column names exactly (capital letters count).")
cell(n, '''rng = np.random.default_rng(SEED); pool = [g for g in age.index if g in pert.index]; null = []
for _ in range(5000):
    r = rng.choice(pool, len(set_genes), replace=False)
    null.append(spearmanr(age.loc[r, "meta_z"], pert.loc[r])[0])
null = np.array(null); p = (1 + (np.abs(null) >= abs(rho)).sum()) / 5001
print("empirical p =", round(p, 4)); plots.fig_perturb(J, rho, p, null); plt.show()''',
     "<b>The honest comparison (here borderline: p about 0.06; report it as such, do not round it down).</b> A transcriptome-wide correlation is easy to get and hard to interpret. Instead, 5,000 random gene sets of the same size give the correlation expected by chance; "
     "the empirical p says whether the coupling genes agree with the perturbation more than random genes do. Repeat for every extra perturbation dataset you add.",
     fig=["fig_perturbation.png"])
cell(n, '''conc = pd.DataFrame({"pert_concord": (np.sign(age.meta_z.reindex(pert.index)) == np.sign(pert)).astype(float)}).dropna()
conc.to_csv(RES / "perturb_concordance.csv"); print(conc.pert_concord.mean().round(3), "of genes agree in direction")''',
     "Saves, for every gene, whether aging and the perturbation move it in the same direction (1) or not (0); notebook 15 uses it as one evidence layer.")

# ================================================================== 15
n = nb("15", "Rank the candidate genes transparently",
       "Combine seven independent evidence layers into one ranking (rank product) with a family-wise permutation p-value, replacing the arbitrary 'bridge score' of the first guide.")
start(n)
cell(n, '''from mlc.genetics import load_orthologs
from mlc.analysis import get_net, rank_product
from mlc import plots
SIM = RAW / "sim"
orth = load_orthologs(SIM / "HOM_MouseHumanSequence_sim.rpt")
h = orth.drop_duplicates("Symbol_mouse").set_index("Symbol_mouse").Symbol_human
rew = pd.read_csv(RES / "rewiring_genes.csv").set_index("gene")
sig = pd.read_csv(RES / "age_signature_meta.csv", index_col=0)
gz = pd.read_csv(SIM / "ad_gene_level_sim.csv").set_index("SYMBOL").ZSTAT
col = pd.read_csv(RES / "coloc_results.csv").set_index("gene")["PP.H4"]
conc = pd.read_csv(RES / "perturb_concordance.csv", index_col=0).pert_concord
hga = pd.read_csv(RES / "human_gene_age.csv").set_index("gene").beta_age
E = pd.DataFrame(index=rew.index)
E["abs_delta"] = rew.delta.abs()
E["abs_meta_z"] = sig.meta_z.abs().reindex(E.index) * sig.concordant.reindex(E.index).astype(float)
E["tf_target"] = E.index.isin(get_net("mouse").target).astype(float)
E["magma_z"] = h.reindex(E.index).map(gz)
E["pp_h4"] = col.reindex(E.index)
E["pert_concord"] = conc.reindex(E.index)
hs = np.sign(h.reindex(E.index).map(hga)); E["human_concord"] = (hs == np.sign(sig.meta_z.reindex(E.index))).astype(float)
E = E.fillna(0); E.to_csv(RES / "evidence_table.csv"); print(E.describe().round(2).loc[["count", "mean", "max"]])''',
     "Builds one table with a row per candidate gene and seven columns, each from an earlier notebook: <b>abs_delta</b> (loss of coupling), <b>abs_meta_z</b> (replicated aging change, only if all cohorts agree), "
     "<b>tf_target</b> (target of a regulator in the network), <b>magma_z</b> (human Alzheimer's gene-level Z), <b>pp_h4</b> (colocalisation), <b>pert_concord</b> (agrees with CCCP), <b>human_concord</b> (agrees with human aging). "
     "Missing evidence is filled with 0 (no support), never with a guess.",
     expect="Count, mean and max for 7 columns over 63 genes (some columns are constant in the simulation, which is harmless: ties share ranks).", fail="FileNotFoundError -> run notebooks 07 to 14 first; each writes one of the input files.")
cell(n, '''cols = ["abs_delta", "abs_meta_z", "tf_target", "magma_z", "pp_h4", "pert_concord", "human_concord"]
ranked = rank_product(E, cols, n_perm=5000, seed=SEED)
ranked.to_csv(RES / "candidate_ranking.csv"); print(ranked.head(10).round(3)[["rank_product", "p_familywise"] + cols])''',
     "<b>Rank product.</b> Within each layer genes are ranked (1 = strongest). A gene's score is the geometric mean of its seven ranks, so a gene must be good in <i>several</i> layers, not brilliant in one. "
     "The permutation p shuffles ranks independently within each layer 5,000 times and asks how often the best gene by chance would score as well as yours (<b>family-wise</b>, so it is already corrected for testing all genes).",
     expect="Ctsd first (family-wise p about 0.003), then Bnip3 (about 0.08), Ctsb, Gba, Pink1...: only the gene with support from many layers (including colocalisation) is significant after correcting for all genes. This is the honest result.")
cell(n, '''plots.fig_candidates(ranked); plt.show()''',
     "<b>Figure 6c-d.</b> (a) The top genes by rank product. (b) Which evidence layers support each gene (darker = stronger rank). "
     "Pick 3-5 genes for experiments <i>before</i> reading their literature, then check novelty.",
     fig=["fig_candidates.png"])

# ================================================================== 16
n = nb("16", "Robustness table and figure inventory",
       "Show that the main finding survives reasonable alternative choices (the 'multiverse'), and list every figure and table to deposit.")
start(n)
cell(n, '''from mlc.dataset import Dataset
from mlc.coupling import *
from mlc import plots
sheet = pd.read_csv(SHEET).set_index("sample_id")
obsd = {c: Dataset.load(PROC / f"{c}_scored").obs for c in cohorts()["mouse"]}
variants = {"primary (M2)": dict(cov=["depth", "pct_mt", "dissoc", "STATE"], min_cells=100, drop=None),
            "no stress covariate": dict(cov=["depth", "pct_mt", "STATE"], min_cells=100, drop=None),
            "min 200 cells": dict(cov=["depth", "pct_mt", "dissoc", "STATE"], min_cells=200, drop=None),
            "without DAM cells": dict(cov=["depth", "pct_mt", "dissoc"], min_cells=60, drop="DAM"),
            "no covariates (M0)": dict(cov=[], min_cells=100, drop=None)}
rows = []
for name, v in variants.items():
    parts = []
    for c, o in obsd.items():
        o = o if v["drop"] is None else o[o.state != v["drop"]]
        sc_ = [x for x in o.columns if x.startswith("state_")] if "STATE" in v["cov"] else []
        r = coupling_per_sample(o, "m_mitophagy", "m_lysosome", covars=[x for x in v["cov"] if x != "STATE"] + sc_,
                                min_cells=v["min_cells"], n_boot=0); r["cohort"], r["level"] = c, "x"; parts.append(r)
    Rv = pd.concat(parts); Rv["log2_age"] = np.log2(Rv.sample_id.map(sheet.age_months).astype(float))
    m = dersimonian_laird(*cohort_slopes(Rv)[["beta", "se"]].values.T); rows.append(dict(variant=name, **m))
rob = pd.DataFrame(rows); rob.to_csv(RES / "robustness.csv", index=False)
print(rob[["variant", "beta", "ci_lo", "ci_hi", "p"]].round(3))''',
     "Repeats the pooled age-slope under five analysis variants (no stress adjustment, stricter minimum cell count, DAM cells removed, no covariates...). "
     "<b>What reviewers want to see:</b> the same sign and a similar size in every row. If one variant flips the answer, report it and explain.",
     expect="Five rows, all about -0.07 with intervals below 0.")
cell(n, '''plots.fig_robust(rob); plt.show()''', "Forest plot of the five variants.", fig=["fig_robustness.png"])
cell(n, '''inv = pd.DataFrame({"file": sorted(p.name for p in FIG.glob("*.png"))})
inv.to_csv(RES / "figure_inventory.csv", index=False)
print("figures:"); print(inv.to_string(index=False))
print("\\nresult tables:", sorted(p.name for p in RES.glob("*.csv"))[:30])''',
     "Lists every figure and results table. Open the PNGs in <b>Inkscape</b> to assemble multi-panel journal figures (export as SVG/PDF; Arial or Helvetica at 5-7 pt). Deposit the tables with the paper.")

# ================================================================== 17
n = nb("17", "Release the code and data (GitHub + Zenodo)",
       "Write the README and citation file, then publish with GitHub Desktop and mint a DOI. <b>Run after your real-data analysis is final.</b>", real_only=True)
start(n)
cell(n, '''readme = """# mlcoupling-microglia
Mitophagy-Lysosome Coupling Index (MLCI) in aging microglia: code for the analysis in <PAPER TITLE>.

## Run
1. `conda env create -f environment.yml` then `conda activate mlc`
2. `jupyter lab`, open notebooks/ and run 00 to 16 in order.
3. Data accessions: GSE121654, GSE129788, GSE99074, GSE148822, GSE157827; GWAS Catalog GCST90027158.
Pre-registration: <OSF LINK>. Frozen gene-set hash: results/real/genesets_frozen.sha256
"""
(ROOT / "README.md").write_text(readme)
cff = """cff-version: 1.2.0
message: "If you use this software, please cite it."
title: "mlcoupling-microglia"
authors:
  - family-names: "Hashmi"
    given-names: "FIRSTNAME"
version: 1.0.0
license: MIT
repository-code: "https://github.com/YOURUSER/mlcoupling-microglia"
"""
(ROOT / "CITATION.cff").write_text(cff)
(ROOT / ".gitignore").write_text("data/raw/\\ndata/processed/\\n*.h5ad\\n.ipynb_checkpoints/\\n__pycache__/\\n")
print("written: README.md, CITATION.cff, .gitignore")''',
     "Writes three small files that every public repository should have. <b>Edit the &lt;placeholders&gt;</b> (paper title, OSF link, your first name, your GitHub user name) before publishing.", run=False)
