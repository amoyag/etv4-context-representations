# Reproducibility levels

## Verified from processed results

The default command sequence in the README uses saved summaries to rebuild every displayed source table and four figure panels. The supplement tables were prepared from these same results, original run plan and prompt identities; its current LaTeX source is edited in Overleaf and synchronized deliberately. This checks the presentation and arithmetic, not how hidden states were produced.

## Historical analysis from stored activations

Original `analyse_results.py` imports `run_experiments.py` helpers, reads activation manifests and verifies file/configuration hashes. The 3B path also needs its `postprocess/` verifier and the scripts in `historical-execution/3b/scripts`. Analysis of activations needs the NPZ files named in the run manifests; these large arrays are not in Git. Previously completed arithmetic-verification receipts are included, but this packaging session has not rerun that level of analysis. A public data archive and its licence must be resolved before claiming complete downloadable reproduction.

## Historical inference

The adapters require the original model checkpoints, regression weights/configuration where applicable, recorded dependency versions and backend qualification. Original assets are identified by the included model manifests; the 3B checkpoint manifest identifies its assets. The original adapters use paths relative to their historical workspace; The machine-specific 3B orchestration wrapper has been removed. Stage guards intentionally remain unchanged. Historical flags such as `E2_authorized: false` record pre-inference preparation and are not overwritten to make packaging scripts runnable. Retained coverage and arithmetic-verification records document subsequent completed runs, but the release does not include the complete internal activation/authorization workspace.

Original inference dependencies are preserved in `historical-execution/original/environment/requirements-macmini.lock.txt`. The 3B execution configuration records PyTorch/CUDA/hardware and numerical floors. The small presentation environment in `requirements-reproduction.txt` is separate. Model packages/weights and third-party papers are not vendored. No activation regeneration, fresh qualification, new model run, changed threshold or head scan is authorized by the reproduction commands above.

## Public-release checks still open

Choose licences for original code and derived data; review attribution and redistribution of UniProt-derived sequences; decide activation archival scope; replace historical host paths in a separately tested runner only if full inference portability is required; confirm final authors and citation metadata. These open release choices do not prevent the private saved-summary package from being inspected and tested now.

The missing end-to-end dependencies and release actions are listed in `REPOSITORY-AUDIT.md`. The three array manifests under `provenance/activation-manifests/` identify 711 biological/control outputs; their NPZ payloads remain outside Git.
