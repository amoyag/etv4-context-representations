# ETV4 context representations

**Code and data supporting the manuscript _Sequence and family context shape ETV4 representations in protein language models_.**

How does a protein language model represent a folded domain when its surrounding sequence or its protein-family context changes? We examine this question in ETV4, a transcription factor whose ETS domain binds DNA and whose flanking regions contribute to regulation. Published experiments on ETV4 constructs and mutations provide the biological reference for this exploratory computational case study.

We compare two kinds of context. With ESM-2 (650M and 3B), we measure how the representation response to a mutation differs between full-length ETV4 and an ETS + CID construct. With ProFam, we measure how adding the N-terminal inhibitory region (NID) changes the representation of a shared ETS-containing region under different family prompts. CID denotes the C-terminal inhibitory region. The models, comparisons and readouts differ; their numerical scales are not directly interchangeable.

This repository makes the study's comparisons, controls, positive observations and negative findings inspectable. It contains processed results, figure-generation code, historical model-execution and analysis code, input definitions, provenance records and supplementary information.

**Publication status:** manuscript under author review; no preprint DOI has been assigned. This repository is currently a private review package. Public access, licences and an archival release remain pending. A private repository link does not give proposal reviewers access.

## Start here

| To understand or inspect… | Open… |
|---|---|
| The protein regions and experimental design | [Figure 1: study design](figures/fig1-design.pdf) |
| Mutation responses in ESM-2 | [Figure 2: mutation responses](figures/fig2-mutation-responses.pdf) |
| Native-order and shuffled-sequence controls | [Figure 3: order controls](figures/fig3-order-controls.pdf) |
| Family prompting and NID context in ProFam | [Figure 4: family conditioning](figures/fig4-profam-conditioning.pdf) |
| Detailed methods and supplemental tables | [Supplementary information PDF](docs/etv4-supplement.pdf) |
| Numerical results and column meanings | [Result tables](data/) and [data dictionary](DATA-DICTIONARY.md) |
| What can currently be reproduced | [Reproducibility guide](REPRODUCIBILITY.md) |

## What the study finds

In ESM-2, the Y401N/Y403G double mutant gives the largest mutation response in ETS + CID and the largest difference in that response between ETS + CID and full-length ETV4 at the primary blocks. It does **not** give the largest ratio of those quantities. The fixed ratio criterion fails at every tested block in both model sizes. Native-order sensitivity meets the specified criterion at the three tested 3B blocks, but not at the three tested 650M blocks.

In ProFam, family prompting changes the representation response to adding NID. The direction of the matched norm comparison reverses with model depth: homologous prompts give smaller responses than all ten non-homologous settings at blocks 8 and 11, but larger responses in eight of ten settings at block 16. Prompt settings are not independent biological replicates. Because ProFam is autoregressive, ETS positions cannot attend to the downstream CID; this comparison does not test CID-to-ETS regulation.

These observations establish context dependence under explicit, controlled comparisons and show why its interpretation depends on the model, depth and readout. They do not establish that the measured changes encode regulation, explain an allosteric mechanism or predict DNA-binding affinity. The contribution is a worked case study for investigating domain–context representations, including controls that expose the limits of a proposed regulatory interpretation.

## Reproduce the figures and tables

No GPU or model checkpoint is needed for this step. Use Python 3.11 and Libertinus fonts (available with TeX Live):

```sh
git clone https://github.com/amoyag/etv4-context-representations.git
cd etv4-context-representations
python3 -m venv .venv
.venv/bin/python -m pip install -r requirements-reproduction.txt
.venv/bin/python repro/prepare_figures.py
.venv/bin/python verify_reproduction.py
```

The generator reads seven frozen inputs under `repro/inputs/` and writes four figures and eleven tables under `repro/figures/` and `repro/data/`. It does not run protein models. Set `LIBERTINUS_FONT_DIR` if the Libertinus Sans OTF directory cannot be found automatically. Verification checks exact table equality, source hashes, numerical bounds and descriptive control counts. Reviewed figure PDFs are in `figures/`; regenerated PDF metadata and font rendering can vary across environments.

This reproduces the presentation and checks the saved numerical results. The historical code is available for inspection, but is not yet a portable end-to-end inference or activation-analysis runner. See [the reproducibility guide](REPRODUCIBILITY.md) for external assets, workspace dependencies and the distinction between these levels of reproduction.

## Data and repository contents

| Location | Contents |
|---|---|
| [data/](data/) | Eleven full-precision result tables |
| [repro/](repro/) | Figure/table generation and frozen summary inputs |
| [historical-execution/original/](historical-execution/original/) | ESM-2 650M and ProFam adapters, analysis, scientific inputs, configuration and checks |
| [historical-execution/3b/](historical-execution/3b/) | ESM-2 3B adapters, analysis, scientific inputs, configuration and checks |
| [historical-execution/3b-analysis/](historical-execution/3b-analysis/) | Cross-model comparison code |
| [provenance/](provenance/) | Source hashes, activation manifests, verification and release records |
| [supplement/](supplement/), [docs/](docs/), [figures/](figures/) | Supplemental LaTeX, review PDF and figures |

The 711 biological/control activation files are saved model representations: 217 ESM-2 650M outputs, 277 ProFam outputs and 217 ESM-2 3B outputs. Their identities and checksums are recorded in [the activation manifests](provenance/activation-manifests/). Their NPZ payloads are held outside Git and are **not yet available for download from this repository or a public archive**. Model weights are not redistributed. Public archival links will be added when the corresponding records are released.

Scientific controls and numerical qualification records are retained. Machine-setup smoke tests and unrelated workflow utilities are excluded from the current tree. No new model inference was performed to prepare this repository.

## Supplementary information

The [PDF](docs/etv4-supplement.pdf) is a convenient reading copy; the final supplement should also accompany the preprint. Its [LaTeX source](supplement/main.tex) is synchronized from the author-edited Overleaf project. Compile it with:

```sh
cd supplement
sh compile.sh
```

The resulting PDF is `supplement/build/main.pdf`. The main manuscript is maintained separately in Overleaf; its public preprint link will be added after posting.

## Citation and reuse

Citation metadata and persistent identifiers will be added to the archived release after the manuscript metadata are finalized. Until then, identify the repository URL and exact commit when referring to this review package.

Licences for original code and derived materials have not yet been adopted. Third-party inputs and model assets retain their own terms. See [rights and release status](RIGHTS-AND-RELEASE.md) and [the repository audit](REPOSITORY-AUDIT.md) for the remaining release work.
