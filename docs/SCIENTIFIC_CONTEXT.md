# Scientific Context

## Overarching Goal

The point of this repository is to measure whether a designed or mutated amino-acid sequence still looks like the wild-type protein in three-dimensional space.

That is a different question from ordinary sequence similarity. Two sequences can be different while preserving the same overall fold, and two sequences can remain fairly similar while introducing a structural failure that matters biologically. If the working hypothesis is that a sequence should retain the WT architecture, then structural comparison is the right level of analysis.

This repository turns that idea into a repeatable batch workflow:

1. fold each candidate sequence with the ESM API
2. compare the predicted structure to the correct WT/reference structure with TM-align
3. store the resulting similarity score in a form that can be analyzed across large batches

## Why Use TM-align

TM-align is a structure alignment method built around the TM-score rather than around RMSD alone.

That matters because RMSD has a well-known weakness: a small number of badly placed residues can inflate the score even when the global fold is still largely correct. For variant screening and batch comparison, that is often the wrong emphasis. In many applications you care more about whether the overall topology is retained than about a few local coordinate deviations.

The TM-align paper addressed exactly that problem. Its scoring function weights close residue pairs more heavily and downweights large local deviations, making the score more sensitive to global structural agreement. In practice, that is why TM-score is widely used as a fold-similarity metric rather than just a raw fit error.

Reference:

- Zhang Y, Skolnick J. TM-align: a protein structure alignment algorithm based on the TM-score. *Nucleic Acids Research* 33(7):2302-2309. https://pmc.ncbi.nlm.nih.gov/articles/PMC1084323/

## What TM-score Means In This Repository

TM-align reports two TM-scores for a pairwise comparison. The alignment is the same, but the normalization length differs:

- one score is normalized by the query or predicted structure length
- the other is normalized by the reference structure length

For this repository, the primary field is:

- `tm_score`: TM-score normalized by the WT/reference structure length

The secondary field is:

- `tm_score_query`: TM-score normalized by the predicted-query length

That choice follows the biological question being asked.

The WT/reference structure is the anchor. It is the structure whose preservation matters. When a candidate sequence is longer, shorter, slightly truncated, or otherwise perturbed, normalizing by the WT length keeps the score tied to the reference object rather than to the idiosyncrasies of the candidate.

In plain terms: this repository is built to answer "how well does this candidate reproduce the WT structure?" not "how internally coherent is the candidate on its own preferred normalization scale?"

## Why Not Use RMSD As The Main Output

RMSD is still useful, and the pipeline records it. But it is not the best headline metric here.

Reasons:

- RMSD is sensitive to a few large local deviations
- RMSD depends strongly on the aligned subset and the overall length scale
- RMSD alone is a poor single-number summary of fold preservation across large and diverse batches

TM-score is not magic, but it is usually a better first-pass structural similarity metric for this exact use case.

## Why Use The ESM API

The repository is intentionally API-based. It does not try to switch among local folding models or maintain separate execution paths.

That design choice keeps the operational story simple:

- one folding backend
- one request format
- one response interpretation layer
- one checkpointed pipeline

The current implementation uses:

- `POST https://biohub.ai/api/v1/fold`
- default model: `esmfold2-fast-2026-05`

The live endpoint was checked during development to confirm the response shape used here. In the observed response, the API returns coordinate arrays, which the pipeline converts into a PDB when raw PDB text is not included.

## What This Pipeline Is Good For

This repository is a good fit when:

- you already know which WT/reference structure each candidate should be compared against
- you need a consistent structural similarity score across many candidate sequences
- you need resumability because folding large batches may take time or fail intermittently
- you want outputs that can be read directly into analysis notebooks or downstream scoring pipelines

Typical examples:

- variant libraries tied to a known parent protein
- design campaigns where candidates are supposed to preserve a known fold
- structure-preservation checks after mutation, recombination, or optimization

## What This Pipeline Does Not Do

It does not try to solve every structural comparison problem.

In particular:

- it does not discover the correct reference automatically
- it does not infer biological function from structure similarity alone
- it does not replace local structural analysis of active sites or interfaces
- it does not treat ESM confidence outputs as the main endpoint

Those are different questions. This repository stays narrow: sequence in, fold, align to WT, report a structural similarity score that is easy to batch.

## Practical Interpretation

The most defensible way to use the output is comparatively:

- rank candidates within the same target
- compare distributions across candidate sets
- inspect low-scoring outliers directly

The score is most meaningful when the reference structure is fixed and the candidate set is coherent. That is the setting this repository was designed for.
