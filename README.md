# ETV4 context representations

Private preparation of the supplement and computational package for *Sequence and family context shape ETV4 representations in protein language models*. This exploratory case study examines paired mutation/context responses in ESM-2 and NID/prompt responses in ProFam. Regulatory correspondence remains unresolved. No public release or preprint DOI exists yet.

## Reproduce the presentation without inference

Use Python 3.11 and Libertinus fonts (available with TeX Live). From the repository root:

```sh
python3 -m venv .venv
.venv/bin/python -m pip install -r requirements-reproduction.txt
.venv/bin/python repro/prepare_figures.py
.venv/bin/python verify_reproduction.py
python3 supplement/build_source.py
cd supplement
sh compile.sh
```

The first script reads seven saved inputs under `repro/inputs/`; it regenerates four figures and eleven tables under `repro/figures/` and `repro/data/`. It neither imports a protein-model library nor downloads weights. If necessary, set `LIBERTINUS_FONT_DIR` to the directory containing Libertinus Sans OTF files. Verification checks exact table equality, frozen hashes, numerical bounds and descriptive control counts. Reviewed figure PDFs remain in `figures/`; regenerated graphics can differ in metadata/font rendering across Matplotlib versions. The supplemental PDF is `supplement/build/main.pdf`.

The checked seven-page review PDF is [docs/etv4-supplement-v1.pdf](docs/etv4-supplement-v1.pdf). Verification details are in [provenance/package-review.md](provenance/package-review.md).

## Contents and reproducibility scope

- `data/`: eleven full-precision source tables; see `DATA-DICTIONARY.md`.
- `repro/`: portable presentation code and saved scientific summaries.
- `supplement/`: standalone Libertinus source and deterministic table builder.
- `historical-execution/original/`: original 650M/ProFam adapters, analysis, tests, inputs, configuration and verification receipts.
- `historical-execution/3b/`: 3B adapters, analysis, postprocessing and execution configuration.
- `historical-execution/3b-analysis/`: original local orchestration record, including historical machine paths.
- `provenance/`: source identities and checks.

The historical scripts are an inspectable execution record, not a newly certified one-command inference workflow. They retain stage guards, historical workspace paths and checkpoint requirements. `REPRODUCIBILITY.md` explains remaining dependencies. No model experiments were rerun while packaging.

The manuscript remains authoritative in Overleaf and is not duplicated here as an editable paper. Supplement integration is pending author review. Public release, code/data licences, complete author metadata and raw-activation archival arrangements remain decisions for the release checkpoint.
