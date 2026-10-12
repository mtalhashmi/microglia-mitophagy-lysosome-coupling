import json, re, sys
sys.path.insert(0, "."); import nbcells
def md(t):
    t = re.sub(r"<b>(.*?)</b>", r"**\1**", t); t = re.sub(r"<i>(.*?)</i>", r"*\1*", t)
    return t.replace("&lt;", "<").replace("&gt;", ">")
for n in nbcells.NBS:
    cells = [{"cell_type": "markdown", "metadata": {}, "source": md(f"# Notebook {n['id']}: {n['title']}\n\n**Goal.** {n['goal']}" + (f"\n\n{n['intro']}" if n["intro"] else "") + ("\n\n> REAL DATA ONLY - skip this notebook in the dry run." if n["real_only"] else "")).splitlines(True)}]
    for i, c in enumerate(n["cells"], 1):
        cells.append({"cell_type": "markdown", "metadata": {}, "source": md(f"**Cell {i}.** {c['what']}" + (f"\n\n*You should see:* {c['expect']}" if c["expect"] else "")).splitlines(True)})
        cells.append({"cell_type": "code", "execution_count": None, "metadata": {}, "outputs": [], "source": c["code"].splitlines(True)})
    nb = {"cells": cells, "metadata": {"kernelspec": {"display_name": "Python 3", "language": "python", "name": "python3"}, "language_info": {"name": "python"}}, "nbformat": 4, "nbformat_minor": 4}
    fn = f"notebooks/{n['id']}_{re.sub('[^a-z0-9]+','_',n['title'].lower()).strip('_')[:40]}.ipynb"
    json.dump(nb, open(fn, "w"), indent=1)
print("notebooks written")
