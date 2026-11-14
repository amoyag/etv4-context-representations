# Reproducibility levels

## Verified from processed results

The default command sequence in the README uses saved summaries to rebuild every displayed source table and four figure panels. The supplement tables were prepared from these same results, original run plan and prompt identities; its current LaTeX source is edited in Overleaf and synchronized deliberately. This checks the presentation and arithmetic, not how hidden states were produced.

## Verified analysis from stored activations

`repro/analyse_activations.py` is a portable NumPy-only adaptation of the two historical analysis scripts. It changes file access and command-line handling, preserving the scientific calculations and leaving the historical scripts unchanged. It never imports the model adapters, loads weights or runs inference.

Extract all three activation ZIP files into one directory (with `650m/`, `profam/`, `3b/` subdirectories). Run the command in the README with a new output directory. Each manifest must match the complete fixed plan; each array checksum, configuration hash, native origin and shape is checked. Every output field is compared with the frozen summaries, including all 107,102 numbers and all categorical results. Float comparisons use rtol=1e-12 and atol=1e-14, solely for float64/library reduction differences; these tolerances do not modify protocol floors or classifications.

This path was tested against all locally held arrays on 14 November 2026. The arrays are openly available at https://doi.org/10.5281/zenodo.24220637. The archive also includes the 51 original and 20 3B numerical-qualification arrays with their receipts. Their presence supports inspection of numerical floors; the portable biological/control analysis consumes the frozen floors rather than recalculating them.

## Historical inference

The adapters require the original model checkpoints, regression weights/configuration where applicable, recorded dependency versions and backend qualification. Original assets are identified by the included model manifests; the 3B checkpoint manifest identifies its assets. The original adapters use paths relative to their historical workspace; The machine-specific 3B orchestration wrapper has been removed. Stage guards intentionally remain unchanged. Historical flags such as `E2_authorized: false` record pre-inference preparation and are not overwritten to make packaging scripts runnable. Retained coverage and arithmetic-verification records document subsequent completed runs, but the release does not include the complete internal activation/authorization workspace.

Original inference dependencies are preserved in `historical-execution/original/environment/requirements-macmini.lock.txt`. The 3B execution configuration records PyTorch/CUDA/hardware and numerical floors. The small presentation environment in `requirements-reproduction.txt` is separate. Model packages/weights and third-party papers are not vendored. No activation regeneration, fresh qualification, new model run, changed threshold or head scan is authorized by the reproduction commands above.

## Remaining release steps

The author approved MIT for original code and CC BY 4.0 for derived materials and authorized public GitHub access. The Zenodo deposit is published and linked. Finalize the manuscript citation after preprint posting. Full model-inference portability is outside the current package; the tested activation-analysis path avoids needing inference to check the published results.
