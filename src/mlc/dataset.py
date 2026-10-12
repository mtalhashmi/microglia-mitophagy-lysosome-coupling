"""Light single-cell container + file readers (numpy/pandas/scipy only; no scanpy needed)."""
import gzip, re, json, urllib.request
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.io import mmread


class Dataset:
    """X = raw counts (cells x genes, sparse); L = log-normalised (set by proc.normalize);
    obs = per-cell table; var = per-gene table; obsm = dict of per-cell arrays (e.g. PCA)."""

    def __init__(self, X, obs, var, L=None, obsm=None, uns=None):
        self.X = sparse.csr_matrix(X)
        self.obs = obs.copy(); self.var = var.copy()
        self.obs.index = self.obs.index.astype(str); self.var.index = self.var.index.astype(str)
        self.L = None if L is None else sparse.csr_matrix(L)
        self.obsm = dict(obsm or {}); self.uns = dict(uns or {})

    n_obs = property(lambda s: s.X.shape[0])
    n_vars = property(lambda s: s.X.shape[1])
    var_names = property(lambda s: s.var.index)
    obs_names = property(lambda s: s.obs.index)

    def __repr__(self):
        return f"Dataset: {self.n_obs} cells x {self.n_vars} genes | obs: {list(self.obs.columns)[:8]}..."

    def subset_cells(self, mask):
        mask = np.asarray(mask)
        idx = np.where(mask)[0] if mask.dtype == bool else mask
        return Dataset(self.X[idx], self.obs.iloc[idx], self.var,
                       None if self.L is None else self.L[idx],
                       {k: v[idx] for k, v in self.obsm.items()}, self.uns)

    def subset_genes(self, mask):
        mask = np.asarray(mask)
        idx = np.where(mask)[0] if mask.dtype == bool else mask
        return Dataset(self.X[:, idx], self.obs, self.var.iloc[idx],
                       None if self.L is None else self.L[:, idx], self.obsm, self.uns)

    def gene_idx(self, genes):
        """Positions of the genes that exist in this dataset (missing genes are skipped)."""
        pos = self.var_names.get_indexer(list(genes))
        return pos[pos >= 0]

    def save(self, folder):
        folder = Path(folder); folder.mkdir(parents=True, exist_ok=True)
        sparse.save_npz(folder / "counts.npz", self.X)
        if self.L is not None: sparse.save_npz(folder / "lognorm.npz", self.L)
        self.obs.to_csv(folder / "obs.csv.gz"); self.var.to_csv(folder / "var.csv.gz")
        if self.obsm: np.savez_compressed(folder / "obsm.npz", **self.obsm)
        (folder / "uns.json").write_text(json.dumps(self.uns, default=str))

    @staticmethod
    def load(folder):
        folder = Path(folder)
        X = sparse.load_npz(folder / "counts.npz")
        L = sparse.load_npz(folder / "lognorm.npz") if (folder / "lognorm.npz").exists() else None
        obs = pd.read_csv(folder / "obs.csv.gz", index_col=0, low_memory=False)
        var = pd.read_csv(folder / "var.csv.gz", index_col=0)
        obsm = dict(np.load(folder / "obsm.npz")) if (folder / "obsm.npz").exists() else {}
        uns = json.loads((folder / "uns.json").read_text()) if (folder / "uns.json").exists() else {}
        return Dataset(X, obs, var, L, obsm, uns)


def concat(datasets, join="outer"):
    names = [set(d.var_names) for d in datasets]
    genes = sorted(set.union(*names) if join == "outer" else set.intersection(*names))
    genes = pd.Index(genes); mats = []
    for d in datasets:
        idx = genes.get_indexer(d.var_names); coo = d.X.tocoo(); keep = idx[coo.col] >= 0
        mats.append(sparse.csr_matrix((coo.data[keep], (coo.row[keep], idx[coo.col[keep]])),
                                      shape=(d.n_obs, len(genes))))
    obs = pd.concat([d.obs for d in datasets])
    assert obs.index.is_unique, "cell names must be unique across samples"
    return Dataset(sparse.vstack(mats).tocsr(), obs, pd.DataFrame(index=genes))


# ------------------------------------------------------------------ readers
def _open(p):
    p = str(p)
    return gzip.open(p, "rb") if p.endswith(".gz") else open(p, "rb")


def _find(folder, pattern, prefix=None):
    hits = sorted(Path(folder).glob(f"{prefix or ''}*{pattern}*"))
    return hits[0] if hits else None


def read_10x(folder, prefix=None):
    """Read a 10x-style triplet (matrix.mtx, genes.tsv/features.tsv, barcodes.tsv), gz or plain."""
    mtx = _find(folder, "matrix.mtx", prefix)
    gen = _find(folder, "features.tsv", prefix) or _find(folder, "genes.tsv", prefix)
    bar = _find(folder, "barcodes.tsv", prefix)
    if not (mtx and gen and bar):
        raise FileNotFoundError(f"No complete 10x triplet in {folder} (prefix={prefix}). Run inventory(folder).")
    with _open(mtx) as f: M = mmread(f).tocsr()
    g = pd.read_csv(gen, header=None, sep="\t"); b = pd.read_csv(bar, header=None, sep="\t")
    sym = g[1] if g.shape[1] > 1 else g[0]
    if M.shape[0] == len(g): M = M.T                      # 10x stores genes x cells
    var = pd.DataFrame(index=pd.Index(_dedupe(sym.astype(str).tolist())))
    return Dataset(M, pd.DataFrame(index=b[0].astype(str)), var)


def _dedupe(names):
    seen = {}; out = []
    for n in names:
        if n in seen: seen[n] += 1; out.append(f"{n}-{seen[n]}")
        else: seen[n] = 0; out.append(n)
    return out


_BARCODE = re.compile(r"^[ACGTN]{10,}([-_.].*)?$|^[A-Za-z0-9_.-]*[ACGT]{12,}[A-Za-z0-9_.-]*$")


def read_dense(path, genes_in_rows=None):
    """Read a dense count table (.txt/.tsv/.csv, optionally gz). Orientation is auto-detected."""
    p = str(path); sep = "," if re.search(r"\.csv(\.gz)?$", p) else "\t"
    df = pd.read_csv(p, sep=sep, index_col=0)
    if genes_in_rows is None:
        looks_barcode = np.mean([bool(_BARCODE.match(str(i))) for i in df.index[:30]])
        genes_in_rows = looks_barcode < 0.5
    if genes_in_rows: df = df.T
    return Dataset(sparse.csr_matrix(df.values.astype("float32")),
                   pd.DataFrame(index=df.index.astype(str)),
                   pd.DataFrame(index=pd.Index(_dedupe(df.columns.astype(str).tolist()))))


def read_h5ad(path):
    try:
        import anndata
    except ImportError as e:
        raise ImportError("Reading .h5ad needs anndata: run  pip install anndata") from e
    a = anndata.read_h5ad(path)
    return Dataset(a.X, a.obs, a.var)


def read_any(path, prefix=None):
    """Auto-detect: folder with 10x triplet, dense table, or .h5ad."""
    p = Path(path)
    if p.is_dir():
        if _find(p, "matrix.mtx", prefix): return read_10x(p, prefix)
        raise ValueError(f"{p}: no matrix.mtx found; point to a single table file instead")
    if p.suffix == ".h5ad": return read_h5ad(p)
    return read_dense(p)


def _kind(name):
    n = name.lower()
    if "matrix.mtx" in n: return "10x matrix"
    if "genes.tsv" in n or "features.tsv" in n: return "10x genes/features"
    if "barcodes.tsv" in n: return "10x barcodes"
    if re.search(r"\\.(txt|tsv|csv)(\\.gz)?$", n): return "dense table"
    if n.endswith(".h5ad"): return "h5ad"
    if n.endswith((".tar", ".zip")): return "archive (extract it)"
    return "other"


def inventory(folder):
    """Print every file with size and a guess of what it is. Use this when a loader complains."""
    for f in sorted(Path(folder).rglob("*")):
        if f.is_file():
            print(f"{f.stat().st_size/1e6:9.1f} MB  {_kind(f.name):22s} {f.relative_to(folder)}")


# ------------------------------------------------------------------ GEO helpers (need internet)
def geo_stub(acc): return acc[:-3] + "nnn"


def geo_list(acc, kind="suppl"):
    """List files GEO offers for a series (kind: 'suppl' or 'matrix')."""
    url = f"https://ftp.ncbi.nlm.nih.gov/geo/series/{geo_stub(acc)}/{acc}/{kind}/"
    html = urllib.request.urlopen(url, timeout=60).read().decode()
    return url, [h for h in re.findall(r'href="([^"?/][^"]*)"', html) if not h.startswith("http")]


def geo_download(url, dest):
    dest = Path(dest); dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists(): print("already downloaded:", dest.name); return dest
    print("downloading", url)
    urllib.request.urlretrieve(url, dest)
    print("saved", dest, f"{dest.stat().st_size/1e6:.1f} MB")
    return dest


def parse_series_matrix(text):
    """Turn a GEO *_series_matrix.txt content into a tidy per-sample table."""
    rows = {}; n_char = 0
    for ln in text.splitlines():
        if not ln.startswith("!Sample_"): continue
        parts = ln.rstrip("\n").split("\t")
        key = parts[0][len("!Sample_"):]
        vals = [v.strip().strip('"') for v in parts[1:]]
        if key.startswith("characteristics"):
            for i, v in enumerate(vals):
                if ":" in v:
                    k, val = v.split(":", 1); rows.setdefault(k.strip().lower(), [""] * len(vals))[i] = val.strip()
        elif key in ("title", "geo_accession", "source_name_ch1", "organism_ch1"):
            rows[key] = vals
    return pd.DataFrame(rows)


def geo_sample_table(acc, dest):
    url, files = geo_list(acc, "matrix"); out = []
    for f in [x for x in files if x.endswith("series_matrix.txt.gz")]:
        p = geo_download(url + f, Path(dest) / f)
        with gzip.open(p, "rt", errors="ignore") as fh: out.append(parse_series_matrix(fh.read()))
    df = pd.concat(out, ignore_index=True); df.insert(0, "cohort", acc); return df


def validate_sheet(sheet, min_per_age_group=3):
    """Check the sample sheet against the inclusion rule and print a plain-language report."""
    need = ["cohort", "sample_id", "species", "age_months", "eligible_primary", "replicate_unit"]
    miss = [c for c in need if c not in sheet.columns]
    if miss: raise ValueError(f"sample sheet is missing columns: {miss}")
    ok = sheet[sheet.eligible_primary.astype(str).str.lower() == "true"].copy()
    ok["age_months"] = ok.age_months.astype(float)
    for c, d in ok.groupby("cohort"):
        counts = d.groupby("age_months").size()
        grp = pd.cut(d.age_months, [0, 9, 15, 100], labels=["young(<9)", "mid(9-15)", "old(>15)"]).value_counts()
        y_, o_ = grp.get("young(<9)", 0), grp.get("old(>15)", 0)
        verdict = "OK" if (y_ >= min_per_age_group and o_ >= min_per_age_group) else "TOO FEW replicates"
        mid = grp.get("mid(9-15)", 0)
        print(f"{c}: {len(d)} eligible samples; young={y_} mid={mid} old={o_} -> {verdict}")
    return ok
