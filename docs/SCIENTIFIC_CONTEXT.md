# Scientific context and provenance

## Question and interpretation

The candidate sequence is the input; its model-predicted structure is the query; an experimentally determined structure is the reference. The endpoint is predicted fold preservation relative to that reference. A candidate with high structural similarity has computational support for a similar global fold, not evidence of measured folding, stability, activity, binding, or biological function.

[TM-align](https://doi.org/10.1093/nar/gki524) finds a structural residue alignment and superposition. TM-score downweights distant pairs relative to squared-distance measures such as RMSD. The pipeline uses `TMalign prediction.pdb reference.pdb`, parses each score's explicit `Chain_1`/`Chain_2` normalization label, and retains the reference-normalized score as `tm_score`. Both lengths, aligned length, RMSD (Å), and aligned identical-residue fraction are retained. No sequence alignment or canonical sequence length is substituted for the structural normalization length.

Reference normalization anchors comparisons to the supplied observed structure. For 1ND4 chain A, that is 255 Cα positions, not the 264-residue deposited sequence. Examine `aligned_length / wt_length` and `aligned_length / predicted_length` alongside both TM-scores. A short query matching a small domain may have a much higher query-normalized score than reference-normalized score. Missing coordinates reduce the structure being assessed; a high score against an incomplete reference cannot establish preservation of omitted regions.

Rank or compare distributions within the same target, reference, chain, and prediction model. A global score can hide local active-site rearrangements, altered side chains, interface defects, or domain orientation changes. Ligand-bound references can differ from unliganded predictions. Single-chain scoring does not assess assemblies, ligand recognition, dynamics, or oligomerization. There is no universal TM-score threshold that validates a design.

The deposited-sequence baseline asks how a prediction of the reference construct compares with its measured conformation. It contextualizes a variant score but does not fully isolate prediction error from mutation effects: prediction uncertainty, model biases, state/ligand differences, and construct differences remain coupled. The pipeline does not analyze confidence tracks, and similarity to a known reference is not an independent benchmark of model generalization or absence of training-set overlap.

## Bundled experimental references

The original `.ent` files are retained unchanged. Example sequences are extracted from chain-A `SEQRES`; observed Cα records come from `ATOM`. No missing residues are filled from coordinates, and no unknown SEQRES residues are silently dropped.

| Detail | TEM (`TEM_pdb.ent`) | APH (`APH_pdb.ent`) |
| --- | --- | --- |
| Accession | [1AXB](https://www.rcsb.org/structure/1AXB) | [1ND4](https://www.rcsb.org/structure/1ND4) |
| Protein and organism | TEM-1 β-lactamase, *Escherichia coli* | Aminoglycoside 3′-phosphotransferase IIa, *Klebsiella pneumoniae* |
| Selected chain | A (only protein chain) | A (B is present but excluded) |
| Method / resolution | X-ray diffraction / 2.00 Å | X-ray diffraction / 2.10 Å |
| Deposited / observed Cα length | 263 / 263 | 264 / 255 |
| Construct context | DBREF maps to P62593 residues 24–286, excluding the precursor's N-terminal signal sequence; PDB numbering starts at 26 and has numbering gaps | DBREF maps to P00552 residues 1–264; `REMARK 465` lists residues 1–9 unresolved in both A and B |
| Bound components excluded from scoring | Phosphonate acylation-transition-state analog (FOS), covalently associated with Ser70 | Kanamycin, ions, acetate, and waters |

The local TEM `COMPND` marks both `ENGINEERED: YES` and `MUTATION: YES`; current RCSB lists zero mutations. No specific sequence substitution is identified by the local header. This discrepancy is unresolved here; the example uses the deposited sequence and describes the construct rather than asserting exact canonical WT identity. APH is also marked engineered in the local header; RCSB lists zero mutations, but its ligand-bound state and missing N-terminal coordinates still matter. “WT reference” names the intended comparison role, not a guarantee that every experimental construct is an unmodified full-length apo WT.

Primary experimental reports are [Maveyraud, Pratt & Samama (1998), DOI 10.1021/bi972501b](https://doi.org/10.1021/bi972501b) for 1AXB and [Nurizzo et al. (2003), DOI 10.1016/S0022-2836(03)00121-9](https://doi.org/10.1016/S0022-2836(03)00121-9) for 1ND4. Header metadata and observed lengths were checked locally; experimental identities, methods, resolutions, and ligand context were cross-checked against those RCSB entries on 2026-10-18.

## Folding service and coordinate contract

Checked against available official documentation on 2026-10-18; no live inference was performed:

| Implementation | Official source / verification |
| --- | --- |
| Bearer-authenticated `POST https://biohub.ai/api/v1/fold` with JSON `sequence` and `model` | [Biohub fold API](https://docs.biohub.ai/api/protein/models/fold) documents endpoint, authentication, and JSON coordinates |
| Default `esmfold2-fast-2026-05` | [Biohub model page](https://biohub.ai/models/esmfold2) uses this identifier; the API reference names it as the default |
| JSON `coordinates` interpreted as `[L,37,3]` with Cα at atom index 1 | [Official SDK fold response](https://github.com/Biohub/esm/blob/main/esm/sdk/forge.py), [ProteinChain atom37 implementation](https://github.com/Biohub/esm/blob/main/esm/utils/structure/protein_chain.py), and [atom ordering constants](https://github.com/Biohub/esm/blob/main/esm/utils/residue_constants.py); local installed SDK also inspected |

The API reference describes nested coordinate arrays but does not spell out every dimension or atom order. The atom37 interpretation follows the official SDK's protein representation and constants (`N`, `CA`, `C`, …). The converter requires exactly 37 coordinate triples per sequence position and exact sequence/coordinate outer-length agreement; an alternate layout is rejected rather than guessed. This contract has not been confirmed with a new live response in this cleanup.

Only Cα coordinates are written to predicted PDBs, even though the service model supports fuller atomic output. An all-null Cα triple is omitted with a warning and a PDB remark listing missing sequence positions; sequence-based residue numbers are preserved. Partial-null, nonfinite, malformed, or PDB-unrepresentable Cα coordinates fail conversion. Fewer than three usable Cα atoms fail. Non-Cα slots are not used for scoring. Coordinate-derived B factors are zero placeholders, not pLDDT. Confidence output is not a reported endpoint.

A raw `pdb` response remains supported for compatibility, but it is not promised by the documented JSON schema. It must have one model/chain, standard residues, finite Cα coordinates, sequence-matching residue identities, and positions in 1..L without insertion codes. It is reduced to one Cα per residue before caching. Prepared references follow the same observed-coordinate policy described in [usage](USAGE.md).

The endpoint and model are retained, with no silent fallback to another provider or original ESMFold. Only `sequence` and `model` are explicitly requested; other inference settings use server defaults, which may change. Response model identity is checked when supplied, and recorded in cache metadata; absence of that field remains an uncertainty. Exact remote weights/defaults and deterministic reproducibility cannot be established by a local cache signature.

## Model references

The official [Biohub ESM repository](https://github.com/Biohub/esm) links ESMFold2 to **Language Modeling Materializes a World Model of Protein Biology**, [2026 preprint](https://www.biorxiv.org/content/10.64898/2026.06.03.729735). The preprint page itself was access-blocked during verification, so the model page, API documentation, and SDK are the checked sources for service behavior. We do not infer undocumented details of the Fast deployment from the preprint.

[Lin et al., *Science* 2023, 379:1123–1130, DOI 10.1126/science.ade2574](https://doi.org/10.1126/science.ade2574) describes the original ESMFold and provides historical context only. It does not document the deployed ESMFold2-Fast model. The previous EasyNano citation was removed because a downstream application preprint is not authoritative documentation of this service or model implementation.
