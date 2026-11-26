# Manuscript figures

This folder redraws the seven figures used in the revised DSA manuscript.
It reads saved numerical outputs; it does not rerun simulations or solve a PDE.
The PDF filenames match the current Overleaf manuscript.

## Rebuild

Use Python 3.9--3.12 with the supplied pinned requirements. From this folder:

```sh
python3 -m pip install -r requirements.txt
python3 make_manuscript_figures.py
```

The `figures/` folder will contain vector PDFs, PNG previews, and an input-hash
manifest. Optional `--data PATH` and `--output PATH` arguments select other folders.

## Appendix B.1 inputs

The continuous-splitting spatial histogram uses the new alpha=1 Expanse run:
`data/dsa_forward_path_dens_pp_split_0_alpha_1.00_expanse_20261123.dat`.
It combines 12,800,000 initial seed histories from four batches.
The two comparison histograms are Juan's saved terminal-weighted direct and
splitting-matched importance-sampling outputs.
