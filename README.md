# ETV4 context representations

Code and data supporting *Sequence and family context shape ETV4 representations in protein language models*. This exploratory study examines mutation/context responses in ESM-2 and NID/prompt responses in ProFam. Regulatory correspondence remains unresolved.

**Private review package.** Public release, licences and archival identifiers are pending. The [repository audit](REPOSITORY-AUDIT.md) separates current capabilities from remaining release work.

## Reproduce figures and tables without model inference

Use Python 3.11 and Libertinus fonts (available with TeX Live):

```sh
python3 -m venv .venv
.venv/bin/python -m pip install -r requirements-reproduction.txt
.venv/bin/python repro/prepare_figures.py
.venv/bin/python verify_reproduction.py
```

The generator reads seven frozen inputs under `repro/inputs/` and writes four figures and eleven tables under `repro/figures/` and `repro/data/`. It does not load protein models or download checkpoints. Set `LIBERTINUS_FONT_DIR` if the Libertinus Sans OTF directory cannot be found automatically. Verification requires exact table equality, source hashes, numerical bounds and descriptive control counts. Reviewed figure PDFs remain in `figures/`; regenerated PDF metadata/font rendering can vary across library versions.

## Supplementary information

The [supplement PDF](docs/etv4-supplement.pdf) is a convenience copy; the submission version should also accompany the preprint. Its [LaTeX source](supplement/main.tex) mirrors the author-edited `supplement.tex` in Overleaf. It has no displayed author/date and no separate AI-assistance section. Overleaf is authoritative during editing; the old automatic prose generator was removed to prevent overwriting author edits.

Compile the synchronized source with:

```sh
cd supplement
sh compile.sh
```

The PDF is `supplement/build/main.pdf`. The manuscript itself remains in Overleaf and is not duplicated as an editable paper here.

## Contents

| Location | Role |
|---|---|
| `data/` | Eleven full-precision result tables; column definitions in `DATA-DICTIONARY.md` |
| `repro/` | Verified presentation reconstruction and saved summaries |
| `historical-execution/original/` | 650M/ProFam scientific adapters, analysis, inputs, configuration and checks |
| `historical-execution/3b/` | 3B scientific adapters, analysis, inputs, configuration and checks |
| `historical-execution/3b-analysis/` | Cross-model comparison code |
| `provenance/` | Source hashes, array manifests, verification and release records |
| `supplement/`, `docs/`, `figures/` | Supplemental source, convenience PDF and reviewed figures |

Historical execution code is retained for inspection; it is not advertised as a portable end-to-end runner. Its workspace guards and required external assets are described in `REPRODUCIBILITY.md`. The 711 biological/control activation files are inventoried but not included. No inference was run during packaging or cleanup.

Picasso smoke/hook-debug scripts, unrelated workflow utilities and obsolete document generators were removed from the current tree. Scientific controls and numerical qualification records that support the results remain. Prior commits and the frozen study workspace preserve provenance.
