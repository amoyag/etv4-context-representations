# E0 source audit — 10 November 2026

Actor: Codex `/root`, standalone E0 coordinator. No model was loaded or executed. Source identities and checksums are recorded in `manifest.json`; extracted tables, sequences and comparisons are in the two registries under `data/derived/`.

## Identity and native sequences

Currie et al., *Structured and disordered regions cooperatively mediate DNA-binding autoinhibition of ETS factors ETV1, ETV4, and ETV5*, Nucleic Acids Research 45 (2017), 2223–2241, DOI [10.1093/nar/gkx068](https://doi.org/10.1093/nar/gkx068). The original main and supplementary PDFs match the setup SHA-256 values. The main Methods, pp. 2224–2226, specify human cDNA, pET28 expression in E. coli, point mutagenesis, purification and EMSA. Supplementary S1 p. 5 contains the full ETV4 sequence. Its text layer contains repeated alignment fragments; the 12 distinct, ordered ETV4 fragments reconstruct 484 residues. They match UniProt P43268, canonical sequence version 3, exactly. The rendered S1 page was inspected, including the termini and mutated segment. The source-defined ETS domain is 339–420; the assay construct 337–430 includes adjacent residues and must not be renamed the structural domain.

All seven substituted residue identities match this sequence: Y401, Y403, E404, E423, E425, F428 and F432. Seven WT boundaries from S3 plus the five double/single-mutant genotypes in FL and 337–436 yield 17 distinct sequences. The two WT sequences used in S10 are already among the seven. S10 does not create additional WT constructs. The 17 rows identify reported native constructs; they do not mean 17 independently tabulated numeric affinities.

The His tag was cleaved during EMSA protein preparation (main pp. 2224–2225); the residual vector sequence was not specified. Native-only inputs are exact with respect to S1 and the named substitutions, not a claim to reproduce every amino acid of the purified reagent. No invented scar is added. The NMR ligation construct has separately described cloning additions including an added cysteine and TEV cleavage sequence preceding the ETS construct (pp. 2225–2226); it is not the EMSA WT 165–436 input. Structural/NMR observations serve as interpretation context only, not additional labels for these 17 inputs. This omission limits extrapolation to the exact expressed reagent; E0 does not resolve an unreported scar.

## Numeric WT endpoints

Supplementary Table S3, PDF p. 4, was read as text and visually checked. KD units are **10^-11 M**, errors are **SEM**, and the final column is the reported replicate count. Values below are transcribed directly, not digitised:

| Native boundaries | KD in 10^-11 M | SEM | n |
|---|---:|---:|---:|
| 337–430 | 6.1 | 0.6 | 25 |
| 337–436 | 12 | 2 | 23 |
| 337–484 | 11 | 2 | 4 |
| 165–430 | 12 | 1 | 3 |
| 165–436 | 60 | 10 | 11 |
| 165–484 | 66 | 9 | 18 |
| 1–484 | 83 | 8 | 35 |

S2 repeats the 337–430 and FL entries; these are not independent additional observations. The registry also preserves S3's reported fold-inhibition values and SEM. Do not replace reported rounded ratios by recomputing means. From the four reported KD means, the factorial ratio is 60 × 6.1 / (12 × 12) = 2.5416667 and its natural log is 0.9328200. This is a derived descriptive contrast, not a new experiment or an uncertainty estimate. Raw titrations, pairing and cross-construct covariance are unavailable; no bootstrap, error propagation assuming independence, or new significance test is performed.

## Mutant endpoints and replicate provenance

Figure S10, supplementary PDF p. 14, was inspected in the original rendered page. S10B compares KD(FL)/KD(337–436) for WT and five mutant groups. The paired constructs differ in residues 1–336 **and** 437–484. The source interprets this ratio as NID-associated inhibition, but the computational comparison cannot attribute every mutant effect solely to the N flank.

Y401N/Y403G weakens DNA binding in the CID construct and relieves the FL-versus-CID inhibition. E404K strengthens DNA binding in both contexts while the source reports preserved autoinhibition. E423K, E425K and F428A/F432A are reported not to alter DNA binding/autoinhibition. The latter absolute-affinity assertions are in the caption with “data not shown”; they are qualitative evidence, not hidden numeric data. Preserved here means the source's reported interpretation, not a newly established equivalence.

S10B n values are WT 23; Y401N/Y403G 5; E404K 4; E423K 4; E425K 4; F428A/F432A 6. S10A displays example isotherms with n = 7 (CID WT), 14 (FL WT), 3 (CID Y401N/Y403G), 4 (FL Y401N/Y403G). Main Figure 7D gives FL mutant n = 6 and CID mutant n = 5. These locators report different subsets/endpoints. We keep S10A and S10B records separate and do not resolve them by averaging counts or borrowing one n for another endpoint. No mutant KD or ratio is digitised, fitted or invented. The primary genotype analysis uses source-supported categories.

The shared assay is EMSA with the 27-bp ETS-site duplex specified in Methods, approximately 10^-12 M labeled DNA, 45 min at room temperature in 25 mM Tris pH 7.9, 0.1 mM EDTA, 60 mM KCl, 6 mM MgCl2, 200 micrograms/ml BSA, 10 mM DTT, 2.5 ng/microlitre poly(dIdC), 10% glycerol. KD was obtained by nonlinear fitting; the paper also assesses an alternative Hill fit. No cell-level regulation or organism-level phenotype is supplied by this assay.

## Conditioning sequences and source limits

UniProt JSON snapshots were downloaded on 10 November 2026 using the public REST API. The exact requests and snapshot hashes are in `manifest.json`. The homologue query returned five reviewed canonical ETV1/5 sequences. Their best 82-residue ETS match is 95.12% identical to ETV4; full-length identity to focal ETV4 is approximately 51–55%. They are close domains in different full-length sequences, not near-duplicate full-length focal inputs. `prompt-selection.json` records every eligibility decision, matching covariate and accession. Regulatory phenotypes of these prompts are not labels for this analysis.

The control query returned a bounded first page of 500 reviewed human nuclear proteins. Exact length matching was possible, but amino-acid-composition and taxonomy matching was imperfect (per-member L1 frequency differences approximately 0.31–0.51). These controls cannot isolate a pure effect of homology. Composition-shuffled homologues preserve exact length and composition as an additional model-only control. Overlapping prompt sets and repeated controls are computational sensitivity settings, not biological replication. Sequence membership in pretraining has not been ruled out. The current P43268 sequence predates training, and ProFam's sources include UniRef and AFDB; no claim of held-out evolutionary generalisation is allowed.

Optional ETV1/5 measured biological contrasts are excluded from this frozen MVP. In particular, the S2 ETV1 full-length boundary is 479, whereas current canonical P50549 is 477 residues. Reconciling that isoform/sequence would be separate work. It does not affect the fully identified ETV1/5 sequences used only as prompts.

## Source-to-design disposition

The native ETV4 primary sequence identities, boundaries and qualitative mutant contrasts are resolved. Exact EMSA cloning scars and raw replicate data remain unavailable and limit the interpretation as stated above. No essential primary numeric endpoint depends on reconstructing either. All comparisons use the actual shared native coordinates; the causal ProFam readout is restricted to information available at those positions. Source-based regulatory claims remain retrospective and single-protein. See the separate novelty record for the narrower possible computational contribution.
