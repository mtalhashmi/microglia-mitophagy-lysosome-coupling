"""Paths and settings shared by every notebook. 'profile' = sample (built-in simulated data) or real."""
import json, os
from pathlib import Path

ROOT = Path(os.environ.get("MLC_ROOT", Path(__file__).resolve().parents[2]))
DATA = ROOT / "data"; RAW = DATA / "raw"
DATA.mkdir(parents=True, exist_ok=True); RAW.mkdir(parents=True, exist_ok=True)
_pf = DATA / "profile.txt"
PROFILE = _pf.read_text().strip() if _pf.exists() else "sample"
PROC = DATA / "processed" / PROFILE
RES = ROOT / "results" / PROFILE
FIG = ROOT / "figures" / PROFILE
for _p in (PROC, RES, FIG): _p.mkdir(parents=True, exist_ok=True)
SEED = 0
MIN_CELLS = 100            # minimum microglia per sample for an MLCI estimate
SHEET = DATA / ("sample_sheet_sample.csv" if PROFILE == "sample" else "sample_sheet.csv")


def set_profile(name):
    assert name in ("sample", "real")
    _pf.write_text(name); print("profile set to", name, "- restart this notebook's kernel before continuing")


def cohorts():
    """Cohort lists, written by notebook 00 (sample) or notebook 01 (real)."""
    f = DATA / f"cohorts_{PROFILE}.json"
    return json.loads(f.read_text()) if f.exists() else {"mouse": [], "human_sn": []}
