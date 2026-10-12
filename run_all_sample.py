import sys, io, json, time, traceback, os, contextlib
import matplotlib; matplotlib.use("Agg")
sys.path.insert(0, "."); os.chdir("notebooks")
import nbcells
only = sys.argv[1:]  # optional list of notebook ids
out = {}
for n in nbcells.NBS:
    if n["real_only"] or (only and n["id"] not in only): continue
    ns = {"__name__": "__main__"}; res = []; t0 = time.time()
    for i, c in enumerate(n["cells"]):
        if not c["run"]: res.append(""); continue
        buf = io.StringIO()
        try:
            with contextlib.redirect_stdout(buf): exec(compile(c["code"], f"nb{n['id']}_cell{i}", "exec"), ns)
            res.append(buf.getvalue())
        except Exception:
            res.append(buf.getvalue() + "\nERROR\n" + traceback.format_exc()); print(f"NB{n['id']} cell {i} FAILED\n", traceback.format_exc()[-1500:]); break
    print(f"notebook {n['id']}: OK ({time.time()-t0:.1f}s)" if not any("ERROR" in r for r in res) else f"notebook {n['id']}: FAILED"); out[n["id"]] = res
    json.dump(out, open("../outputs.json", "w"))
