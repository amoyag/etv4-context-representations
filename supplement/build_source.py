"""Generate supplemental tables from frozen saved results; no model inference."""
from pathlib import Path
import csv,json
ROOT=Path(__file__).resolve().parents[1]
def read(name):
 with (ROOT/'data'/name).open() as f:return list(csv.DictReader(f,delimiter='\t'))
def esc(s):
 s=str(s)
 for a,b in [('&',r'\&'),('%',r'\%'),('_',r'\_'),('#',r'\#')]:s=s.replace(a,b)
 return s
def table(caption,headers,rows,widths=None):
 spec=widths or '@{}'+'l'*len(headers)+'@{}'
 head=' & '.join(headers)+r'\\'
 return '\n'.join([(r'\Needspace{'+str(len(rows)+7)+r'\baselineskip}') if len(rows)<=10 else '',r'\begingroup\small',r'\begin{longtable}{'+spec+'}',r'\caption{'+caption+r'}\\',r'\toprule '+head+r'\midrule\endfirsthead',r'\toprule '+head+r'\midrule\endhead',*[' & '.join(str(x) for x in row)+r'\\' for row in rows],r'\bottomrule\end{longtable}',r'\endgroup'])
run=json.loads((ROOT/'historical-execution/original/design/run-plan.json').read_text())
original=json.loads((ROOT/'historical-execution/original/design/execution-config.json').read_text())
three=json.loads((ROOT/'provenance/execution-config-3b.json').read_text())
prompts=json.loads((ROOT/'repro/inputs/exports/results-v1/data/derived/prompt-panel.json').read_text())
parts=[r'''\documentclass[a4paper,11pt]{scrartcl}
\usepackage[margin=2.35cm]{geometry}
\usepackage{fontspec,unicode-math}
\setmainfont{Libertinus Serif}\setsansfont{Libertinus Sans}\setmonofont{Libertinus Mono}\setmathfont{Libertinus Math}
\usepackage{microtype,amsmath,booktabs,longtable,array,graphicx,needspace}
\usepackage[hidelinks]{hyperref}
\setkomafont{disposition}{\normalfont\bfseries}
\setlength{\emergencystretch}{2em}
\renewcommand{\thesection}{S\arabic{section}}
\renewcommand{\thetable}{S\arabic{table}}
\title{Supplementary material\\\large Sequence and family context shape ETV4 representations in protein language models}
\author{Aurelio Moya Garc\'ia}\date{13 November 2026; review version 1}
\begin{document}\maketitle
This supplement documents the existing computational comparisons, controls and provenance. It does not add model experiments. The complete manuscript remains under author review. Tables retain all tested blocks and unfavourable outcomes; the accompanying TSV files retain full numerical precision.
\section{Biological regions, constructs and observations}
Human ETV4 uses native UniProt P43268 sequence version 3. Region names follow the manuscript. NID is predominantly disordered; CID includes the inhibitory H4 helix and is not wholly disordered. ETS + CID (337--436) also includes the two-residue junction. Shortened CID denotes 421--430. The published observations are from Currie et al., \emph{Nucleic Acids Research} 45, 2223--2241 (2017), \href{https://doi.org/10.1093/nar/gkx068}{doi:10.1093/nar/gkx068}, including Supplementary Tables S1, S3 and S10. They are biological reference observations, not outputs of the protein models.
''']
parts.append(table('Region definitions used throughout the study.',['Region','Native residues'],[['N-terminal extension','1--164'],['NID','165--336'],['Junction','337--338'],['Folded ETS domain','339--420'],['CID','421--436'],['C-terminal tail','437--484']]))
parts.append(r'''Seven WT constructs and five mutant genotypes in each of two contexts give 17 distinct native input sequences. Y401N/Y403G and E404K affect ETS; E423K, E425K and F428A/F432A affect CID. The model-input registry includes the entire input sequence even when a substituted residue is outside the readout mask. Y401N/Y403G weakens DNA binding in ETS + CID and relieves NID-associated inhibition in the published comparison. E404K increases affinity while inhibition is reported as preserved. The CID substitutions are controls with binding and inhibition reported as preserved; some source statements lack displayed quantitative affinity data. We do not infer numerical equivalence from those statements.

Full length versus ETS + CID changes the N-terminal extension, NID and the C-terminal tail. It does not isolate NID. The ProFam WT comparison instead adds NID to ETS + CID, keeping the same prompt before each target.
''')
regions=['ETS + shortened CID','ETS + CID','ETS + CID + tail','NID + ETS + shortened CID','NID + ETS + CID','NID + ETS + CID + tail','Full length']
rows=read('published-WT-affinities.tsv')
parts.append(table('Published WT affinities: means and reported SEM, converted to nM. Lower $K_D$ means tighter DNA binding. Replicate counts remain in the source TSV.',['Construct','$K_D$ (nM)','SEM (nM)'],[[regions[i],f"{float(r['KD_nM']):.3f}",f"{float(r['SEM_nM']):.3f}"] for i,r in enumerate(rows)]))
parts.append(r'''\section{Exact readout masks and model measurements}
The structural ETS domain is 339--420 (82 residues). The shared readout window 337--430 has 94 residues: two junction, 82 ETS and ten CID positions. The ESM mutation readout excludes the union of substituted positions inside this window: 401, 403, 404, 423, 425 and 428. Its 88 positions comprise two junction, 79 ETS and seven CID positions. F432 remains in the input but is outside the window. The WT factorial and ProFam readouts use all 94 positions. The complete masks are also machine-readable in the original run plan.

\noindent\textbf{ESM mutation mask (88 positions):} ''')
parts.append(', '.join(map(str,run['masks']['genotype_primary']))+r'.\par')
parts.append(r'''At each position, the post-block vector is divided by its Euclidean norm. For a vector difference $w$, RMS denotes $\sqrt{|S|^{-1}\sum_{i\in S}\|w_i\|_2^2}$. Let $u$ be mutant minus WT in ETS + CID and $v$ mutant minus WT in full length. Then $D=\mathrm{RMS}(u)$ and $M=\mathrm{RMS}(v-u)$. $M$ is not a difference of scalar magnitudes. $D$ is not a binding energy. $M/D$ is a prespecified diagnostic without an established biological normalization.

The WT contrast is $J=z(NC)-z(N)-z(C)+z(0)$, with $0$ = ETS + shortened CID, $N$ = NID + ETS + shortened CID, $C$ = ETS + CID, and $NC$ = NID + ETS + CID. ProFam measures $q_P=z(NC\mid P)-z(C\mid P)$ for prompt $P$. Its left-to-right architecture prevents ETS residues from seeing downstream CID and the whole readout from seeing 431--436. The resulting zero ProFam $J$ is an architectural check.
\section{Numerical allowances and every fixed ESM block}
Numerical qualification used artificial inputs and froze $\eta$ before biological inference. Denominators must exceed $4\eta$. Ratio bounds are $L=\max(0,M-8\eta)/(D+4\eta)$ and $U=(M+8\eta)/(D-4\eta)$. They are sensitivity bounds under the measured implementation discrepancy, not confidence intervals or a guarantee covering all numerical error.
''')
eta=[]
for name,config in [('650M',original['models']['esm2_650m']),('3B',three['models']['esm2_3b']),('ProFam',original['models']['profam1'])]:
 for block,value in config['eta'].items():eta.append([name,block,f'{1000*value:.3f}',f'{4000*value:.3f}',f'{8000*value:.3f}'])
parts.append(table('Frozen numerical allowances, all displayed in units of $10^{-3}$. Exact values remain in the execution configurations.',['Model','Block',r'$\eta$',r'$4\eta$',r'$8\eta$'],eta))
rows=read('esm-all-fixed-blocks.tsv')
parts.append(table('All fixed-block ESM mutation measurements. Bounds use the model/block-specific allowance.',['Model/block','Genotype','$M$','$D$','$M/D$','$L$','$U$'],[[r['model']+'/'+r['block'],r['genotype'].replace('_','/'),*[f'{float(r[k]):.3f}' for k in ['M','D','M_over_D','L','U']]] for r in rows]))
parts.append(r'''The original directional rule required resolved denominators, raw-$M$ separation, ratio-bound separation and native-order separation for Y401N/Y403G. The overall specificity classification remained unmet in both models at every fixed block. The 3B native-order condition passed without rescuing ratio-bound separation. These failures do not establish absence of regulatory information.
\section{Sequence-order controls and seeds}
Twenty full-length controls shuffle residues 1--336 (N-terminal extension + NID), preserving composition, length and endpoints, with each shuffled prefix shared by WT and every mutant. Seeds are 13001--13020. WT factorial controls shuffle NID with seeds 11001--11020 and the six-residue CID completion with seeds 12001--12020, preserving endpoints. Only four CID-completion positions are permutable; repeated orders and one native-order draw are retained. Exceedance counts are descriptive, not $p$-values.
''')
parts.append(table('Y401N/Y403G fixed native-order test. Margin and required allowance are in units of $10^{-3}$. The unrounded comparison determines the outcome.',['Model/block','Margin',r'Required $8\eta$','Pass'],[[r['model']+'/'+r['block'],f"{1000*float(r['margin']):.3f}",f"{1000*float(r['required_8eta']):.3f}",r['passes']] for r in read('esm-native-order-margins.tsv')]))
parts.append(r'''At the primary block of each ESM model, ten of twenty WT factorial controls equal or exceed the native RMS($J$). This is a different contrast from mutation-dependent native-order sensitivity. All individual WT and mutation-control values are retained in the data files.
\section{Protein prompts and all ProFam response magnitudes}
Five pairs selected from five reviewed ETV1/ETV5 entries (seed 20261110) were used in both orders. Controls are exactly length-matched nonhomologous prompts and composition-shuffled homologous prompts; the empty prompt supplies a baseline. Prompt-shuffle seeds are $14000+10s+m$, with set $s=1,\ldots,5$ and member $m=0,1$. There are 31 named settings. The ten nonhomologous settings contain only four distinct ordered prompts. They are not independent biological replicates. The exact sequences, ordering, hashes and token-prefix lengths are in \texttt{prompt-panel.json}.
''')
pr=[]
for key,p in prompts.items():pr.append([esc(key),esc(', '.join(p['accessions']) or 'none'),'/'.join(str(len(s)) for s in p['sequences']) or '0',p['token_prefix_length']])
parts.append(table('Prompt identity and lengths. Accessions for shuffled homologues identify the source sequences; sequences were shuffled before use. Prefix tokens include the model document/separator tokens.',['Setting','Accessions','AA lengths','Prefix tokens'],pr))
responses={}
for r in read('profam-prompt-responses.tsv'):responses.setdefault(r['prompt'],{})[r['block']]=r['rms_q']
parts.append(table('RMS response to adding NID under every prompt at all fixed ProFam blocks.',['Setting','Block 8','Block 11','Block 16'],[[esc(key),*[f'{float(responses[key][b]):.3f}' for b in ['8','11','16']]] for key in prompts]))
parts.append(r'''Homologous responses are smaller than matched nonhomologous responses in all ten settings at blocks 8 and 11. At block 16 they are larger in eight of ten settings. These counts and paired differences of magnitudes were prepared descriptively for the manuscript from existing outputs. They are not the same as vector differences between responses. Prompt-order and NID-shuffle values are included in the accompanying TSVs. This ordering reversal does not establish a reversal of biological inhibition.
\section{Execution and reproducibility records}
Original ESM-2 650M and ProFam inference used float32 on Apple MPS with CPU qualification; their primary blocks were 22 and 11. The original environment recorded Python 3.11.15, PyTorch 2.15.0, fair-esm 2.0.0, ProFam 0.1.8, Transformers 4.57.6 and NumPy 2.4.6. The 3B follow-up used an A100-SXM4-40GB, PyTorch 2.15.0+cu132, CUDA 13.2, TF32 disabled, batch size one and primary block 24. Model weights were fixed; summaries used float64. Inference seed was 19001.

There were 217 biological/control forwards for each ESM model and 277 for ProFam. The 3B artificial qualification was a separate 20-forward job. These counts are computational workloads, not independent samples. The follow-up was motivated by the 650M result and is not an independent confirmation or an isolated parameter-count comparison.

Separate arithmetic implementations previously checked 78,789 original and 28,311 follow-up values. Those checks were same-agent verification. Their receipts are included; packaging did not repeat model inference. The default repository workflow reconstructs figures and tables from saved summaries. Historical adapters, inputs, masks, seeds, configuration and checkpoint identities are retained separately. Raw activations and model weights are outside the Git repository; a full activation-based reproduction requires those assets. Public archival arrangements remain pending.
\section{AI assistance and human scientific direction}
OpenAI Codex assisted source reconciliation, design and code proposals, authorized computational execution, numerical/source checks, interpretation proposals, figure preparation and manuscript/supplement drafting. Earlier interpretability literature synthesis used Claude Code. Exact deployed assistant versions are incomplete and are not inferred from recommended settings. ESM-2 and ProFam are the scientific models studied, not manuscript authors.

Aurelio questioned the biological justification of M/D, clarified regions and constructs, requested comparable plot scales, adopted the combined target-context/family-context argument and deferred attention-head experiments. Those exchanges changed the interpretation and writing. Final author approval is still pending. Longhand reflection was offered; its use was not reported. A separate Codex coherence check did not constitute independent validation of the full study. The collaboration was not compared with a human-only workflow.
\begin{longtable}{@{}p{3.6cm}p{11cm}@{}}
\caption{Task-specific research assistance. Available product identities are reported without reconstructing missing deployed model versions.}\\
\toprule Task & Assistance, direction and verification\\\midrule\endfirsthead
\toprule Task & Assistance, direction and verification\\\midrule\endhead
Biological source/design & Codex reconciled sources, constructs and comparison proposals. Aurelio directed the scientific question and later revised the interpretation.\\
Implementation/execution & Codex developed code and carried out earlier authorized model execution. Frozen configurations and completed-run receipts identify those operations.\\
Numerical checks & Codex performed separate arithmetic checks and source/output verification. These were not independent scientific review.\\
Literature synthesis & Codex and earlier Claude Code work supplied reading notes and syntheses. Source reading depth varies and is recorded.\\
Interpretation/writing & Codex proposed reasoning and prose; Aurelio questioned M/D, clarified the system, chose scope and revised the argument in exchanges.\\
Figures/supplement & Codex prepared plots, captions, source tables and this supplement from existing results. Fixed blocks and negative outcomes were retained.\\
Final scientific approval & Remains with the human author. Neither tool use nor a completed PDF constitutes final approval or a measured collaborative advantage.\\
\bottomrule\end{longtable}
\section{Machine-readable data and release status}
The eleven TSV files cover all fixed ESM metrics, native-order margins, primary shuffled responses, WT factorial controls, all ProFam prompt responses, NID shuffles, vector contrasts, prompt-order contrasts, matched magnitude differences and published WT affinities. \texttt{DATA-DICTIONARY.md} defines their columns. Frozen source hashes and a reconstruction report accompany the repository.

The repository is private during preparation. Code/data licences, complete author metadata, activation archival scope and public release remain subject to author decisions. The present supplement is a review artifact, not a submitted or published version.
\end{document}
''')
(ROOT/'supplement/main.tex').write_text('\n'.join(parts))
print('Supplement source generated from saved data')
