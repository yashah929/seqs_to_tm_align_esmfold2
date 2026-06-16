# Scientific Context

## Overarching Goal

The purpose of this repository is to measure whether a designed or mutated amino-acid sequence remains close to the wild-type structure in three-dimensional space.

That is a different question from ordinary sequence similarity. Two sequences can be different while preserving the same overall fold, and two sequences can remain fairly similar while introducing a structural failure that matters biologically. If the working hypothesis is that a sequence should retain the WT architecture, then structural comparison is the right level of analysis.

The workflow is simple:

1. fold each candidate sequence with the ESM API
2. compare the predicted structure to the correct WT/reference structure with TM-align
3. store the resulting similarity score in a form that can be analyzed across large batches

## Why Use TM-align

TM-align is a structure alignment method built around the TM-score rather than around RMSD alone.

That matters because RMSD has a well-known weakness: a small number of badly placed residues can inflate the score even when the global fold is still largely correct. For variant screening and batch comparison, that is often the wrong emphasis. In many applications you care more about whether the overall topology is retained than about a few local coordinate deviations.

The TM-align paper addressed this directly. Its scoring function weights close residue pairs more heavily and downweights large local deviations, making the score more responsive to global structural agreement. That is the reason TM-score remains useful in large structural screens where fold preservation, rather than local coordinate perfection, is the main endpoint.

Primary reference:

- Zhang Y, Skolnick J. TM-align: a protein structure alignment algorithm based on the TM-score. *Nucleic Acids Research* 33(7):2302-2309. https://pmc.ncbi.nlm.nih.gov/articles/PMC1084323/

## What TM-score Means In This Repository

TM-align reports two TM-scores for a pairwise comparison. The alignment is the same, but the normalization length differs:

- one score is normalized by the query or predicted structure length
- the other is normalized by the reference structure length

For this repository, the primary field is:

- `tm_score`: TM-score normalized by the WT/reference structure length

The secondary field is:

- `tm_score_query`: TM-score normalized by the predicted-query length

That choice follows from the biological question.

The WT/reference structure is the anchor. It is the structure whose preservation matters. When a candidate sequence is longer, shorter, slightly truncated, or otherwise perturbed, normalizing by the WT length keeps the score tied to the reference object rather than to the idiosyncrasies of the candidate.

Put plainly, the repository is meant to answer "how well does this candidate reproduce the WT structure?" rather than "how well does the candidate score after being normalized to itself?"

## Why Not Use RMSD As The Main Output

RMSD is still useful, and the pipeline records it. But it is not the best headline metric here.

Reasons:

- RMSD is sensitive to a few large local deviations
- RMSD depends strongly on the aligned subset and the overall length scale
- RMSD alone is a poor single-number summary of fold preservation across large and diverse batches

TM-score is not magic, but it is usually a better first-pass structural similarity metric for this exact use case.

## Why Use The ESM API

The repository is intentionally API-based. It does not try to maintain separate local-model and remote-model code paths.

That design choice keeps the operational story simple:

- one folding backend
- one request format
- one response interpretation layer
- one checkpointed pipeline

The current implementation uses:

- `POST https://biohub.ai/api/v1/fold`
- default model: `esmfold2-fast-2026-05`

The endpoint was checked during development to confirm the response shape used here. In the observed response, the API returns coordinate arrays, which the pipeline converts into a PDB when raw PDB text is absent.

For background on the model family, the relevant reference is the original ESMFold paper:

- Lin Z, Akin H, Rao R, Hie B, Zhu Z, et al. Evolutionary-scale prediction of atomic-level protein structure with a language model. *Science*. 2023;379(6637):1123-1130. https://www.science.org/doi/10.1126/science.ade2574

For the newer `ESMFold2-Fast` / `ESMFold2` naming used by the current API model, I have not located a standalone methods paper. The most direct citable source I found is a recent preprint that explicitly uses those model names:

- Hu Y, Cheng W, Wang J, Liu Y. EasyNano: rapid epitope-targeted nanobody CDR design via differentiable distogram optimization with ESMFold2. *arXiv*. 2026. https://arxiv.org/abs/2606.12772

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

The most defensible use of the output is comparative:

- rank candidates within the same target
- compare distributions across candidate sets
- inspect low-scoring outliers directly

The score is most informative when the reference structure is fixed and the candidate set is coherent. That is the setting this repository was built for.

## Additional References

- Zhang Y, Skolnick J. TM-align: a protein structure alignment algorithm based on the TM-score. *Nucleic Acids Research*. 2005;33(7):2302-2309. https://pmc.ncbi.nlm.nih.gov/articles/PMC1084323/
- Lin Z, Akin H, Rao R, Hie B, Zhu Z, et al. Evolutionary-scale prediction of atomic-level protein structure with a language model. *Science*. 2023;379(6637):1123-1130. https://www.science.org/doi/10.1126/science.ade2574
- Hu Y, Cheng W, Wang J, Liu Y. EasyNano: rapid epitope-targeted nanobody CDR design via differentiable distogram optimization with ESMFold2. *arXiv*. 2026. https://arxiv.org/abs/2606.12772
