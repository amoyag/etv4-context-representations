# Repository and sharing audit — 14 November 2026

## Conclusion

The private repository is a tested package for reproducing the manuscript presentation from saved results. It is not yet a complete public study archive. All tracked files were inventoried and assigned a disposition in `provenance/repository-file-audit.tsv`; source identities, script imports/entry points, reproduction checks, supplemental source/PDF and release documentation were examined. This is not a new model run or independent scientific review.

## Where each item belongs

| Material | Placement | Status and reason |
|---|---|---|
| Analysis/metric/plotting code | GitHub and versioned archive | Keep: necessary to inspect and reproduce reported calculations |
| Exact input sequences, prompts, masks, seeds and model/block configuration | GitHub and archive | Keep: define the computational comparisons; source attribution retained |
| Full-precision processed results and data dictionary | GitHub and archive | Keep: eleven tables and canonical summaries support figures and negative findings |
| Relevant controls, numerical floors, verification receipts and source hashes | GitHub and archive | Keep: they support interpretation; they are distinct from setup-only smoke tests |
| Software/checkpoint versions and upstream download identities | GitHub | Keep: document assets without redistributing third-party model weights |
| Supplement PDF | Preprint submission; optional GitHub convenience copy | Kept here, but it should not exist only in a mutable repository |
| Supplement LaTeX | GitHub snapshot; Overleaf editing source | Keep synchronized. No automatic generator may overwrite author edits |
| Raw biological/control activations and qualification data needed to audit numerical floors | Suitable data archive, linked from GitHub/manuscript | Missing public deposit. 711 biological/control files are now inventoried; qualification payload archival scope remains to be finalized |
| Model weights and installed third-party source | Upstream provider | Do not duplicate; retain exact identifiers/hashes and access instructions |
| Picasso smoke/hook-debug code | Internal infrastructure history | Removed from current repository; not evidence for a manuscript result |
| Local monitoring, staging, workflow tests and old manuscript generators | Internal project history | Removed where unrelated; the guard helper imported by the archived scientific adapter remains to preserve that source |
| Credentials, private chats/notebooks, caches, virtual environments, scheduler/account setup and failed infrastructure logs | Outside release | Not needed for publication or included in the selected release package |

## Findings addressed

1. The manuscript previously lacked the actual GitHub URL. Its availability paragraph now links the repository and accurately states private status and omitted raw arrays/weights.
2. The repository supplement was stale: it still showed an author/date and an AI-assistance section. It now mirrors the current Overleaf supplement and has been recompiled.
3. The supplement generator could recreate superseded text and overwrite author edits. It was removed; the LaTeX source is edited in Overleaf and copied deliberately.
4. Twelve setup/workflow/obsolete files were removed, including `smoke_3b.py` and `test_hook_layout.py`. Scientific model controls, arithmetic tests and numerical qualification records were retained.
5. Array identities were implicit in the external workspace. Three activation manifests now identify all 711 biological/control files. They are not substitutes for the array payloads.

## Material gaps before public release

- **Public availability:** a private GitHub link is not public data/code access. Keep it private until the author authorizes release; make the availability statement match actual access at submission.
- **Licensing:** original code and derived data do not yet have release licences. Suggested choices for author review are MIT for original code and CC BY 4.0 for derived data/documentation where applicable; third-party sequence/model terms must be preserved. These are suggestions, not licences already applied or universal server requirements.
- **Persistent version:** create a versioned code/data release and archive it (for example in Zenodo) after approval. A code DOI is distinct from a preprint DOI. This is a reproducibility recommendation, not a claim that either candidate server universally requires a GitHub DOI.
- **Raw data:** archive the stored activations necessary to recompute the summaries, with manifests/checksums and the numerical-qualification evidence needed for the frozen floors. Until then, claim saved-summary reproduction only. Do not state that every study output is publicly available.
- **Executable scope:** original analysis imports historical stage guards and needs WORK_STATE/receipts absent from the release. The 3B runner expects a package manifest and external checkpoints that are also not supplied as a runnable installation. Full execution portability remains unverified. Preserve the archival code and implement a separately documented portable entry point if end-to-end rerunning is a release objective; do not silently bypass old safeguards or rerun models.
- **Citation/authorship:** finalize contributors, version and citation metadata (such as `CITATION.cff`) after author confirmation. Do not invent an archive or manuscript DOI.
- **Final consistency:** freeze matching manuscript, supplement, data and code versions; verify their links and identifiers; apply the selected server's current fields and policies.

## Requirements versus recommendations

No server or journal has yet been selected, so there is no single completed compliance checklist. Open access to the paper, access to underlying data/code, and reuse licences are separate matters. Neither candidate requires every artifact to be stored specifically on GitHub.

Preprints.org's current instructions require associated data to be available absent legal/confidentiality constraints; they encourage a recognized data repository and permit supplemental files with the preprint. LaTeX submissions must include compilation inputs. Its AI-reporting requirements apply to the manuscript; they do not require duplicating an extensive AI declaration in this supplement. [Author instructions](https://www.preprints.org/instructions-for-authors), checked 14 November 2026.

GitHub documents Zenodo archiving and release DOIs and recommends including a licence so readers know how work may be reused. [GitHub documentation](https://docs.github.com/en/repositories/archiving-a-github-repository/referencing-and-citing-content), checked 14 November 2026. Zenodo's GitHub integration requires a public repository; no activation/publication was performed.

The current bioRxiv submission page returned HTTP 403 during this check. Its exact current requirements remain to be verified at venue selection; medRxiv policies and other journals' policies were not substituted for bioRxiv requirements.

## Verification scope

The checks cover preserved source identities, eleven regenerated TSVs, all fixed ESM metric rows/bounds, masks, ProFam ordering counts, supplemental source synchronization, LaTeX compilation/rendering and relevant arithmetic/plan tests. No new model experiment or full raw-activation rerun was performed. See the dated verification record.
