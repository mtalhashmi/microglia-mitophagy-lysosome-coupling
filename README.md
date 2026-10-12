# mlcoupling-microglia

Mitophagy-Lysosome Coupling Index (MLCI) in aging microglia: a replicated, pre-registered, specificity-controlled pipeline.
Needs only numpy, pandas, scipy, scikit-learn, matplotlib (no scanpy).

## Quick start (dry run on built-in simulated data, about 2 minutes)
    conda env create -f environment.yml
    conda activate mlc
    python run_all_sample.py          # runs notebooks 00, 03-16 headless and reports OK / FAILED
    jupyter lab                       # then open notebooks/ and run 00, 03 ... 16 in order

## Real data
Run notebooks 01 and 02 (download + load), then 03-16. Simulated and real outputs are kept apart
(data/processed/{sample,real}, results/{sample,real}).

## Layout
src/mlc/   dataset.py (readers), qc.py, proc.py, coupling.py (MLCI + statistics), analysis.py, genetics.py, plots.py, synth.py
notebooks/ 00 ... 17, one step each; nbcells.py = the source of every notebook cell
