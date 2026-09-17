# ClonoMappeR

An R package to map BCR sequences against public, antigen-specific B-cell clone databases.

Given a set of query BCR contigs (e.g. from single-cell VDJ sequencing), ClonoMappeR finds the closest known antigen-specific clone for each cell by matching V/J gene usage and measuring CDR3 amino acid sequence similarity. You can then apply a distance threshold to call cells as antigen-specific.

Documentation: <https://hungms.github.io/ClonoMappeR/>

## Installation

```r
install.packages("devtools")
devtools::install_github("hungms/ClonoMappeR")
```

## Overview

The package has two main functions.

`get_reference()` retrieves a curated reference database of antigen-specific BCR contigs for a given context (e.g. Sars-CoV-2), optionally filtered by organism, publication and epitope.

`find_publicBCR()` matches a query dataset against a reference. It selects query contigs whose heavy/light chain V and J genes match a reference contig, then computes how similar each query CDR3 amino acid sequence is to the reference CDR3 sequences, and reports the closest match per cell.

## Quick start

```r
library(ClonoMappeR)
library(dplyr)

# Example query: plasma cell BCRs after Sars-CoV-2 mRNA-1273 vaccination
path  <- system.file("extdata", "GSE219098_Plasmablast_v1D14.csv", package = "ClonoMappeR")
query <- read.csv(path, row.names = 2)   # row names are the cell barcodes

# Reference: known Sars-CoV-2 S2-specific BCRs
reference <- get_reference(context = "SarsCoV2", org = "human", epitope = "S2")

output <- find_publicBCR(
  query       = query,
  reference   = reference,
  CDRH3       = "CDRH3",   # required: query column holding the CDRH3 AA sequence
  VH          = "VH",      # optional: query column holding the heavy chain V gene
  JH          = "JH",      # optional: query column holding the heavy chain J gene
  dist_method = "levenshtein"
)

# Call antigen-specific cells
ag_specific <- output %>% filter(CDRH3_dist < 0.1)
```

## Input format

### Query

The query is a `data.frame` or `data.table` where **each row is one cell, or one unique heavy/light chain pair**. Row names are used as cell identifiers and are returned in the output as the `barcodes` column; if the query has no meaningful row names, rows are numbered `1:N` instead.

Query columns can be named anything you like — you tell `find_publicBCR()` which column is which via its arguments:

| Argument | Required | Contents |
| -------- | -------- | -------- |
| `CDRH3`  | Yes      | Heavy chain CDR3 amino acid sequence |
| `VH`     | No       | Heavy chain V gene |
| `JH`     | No       | Heavy chain J gene |
| `CDRL3`  | No       | Light chain CDR3 amino acid sequence |
| `VL`     | No       | Light chain V gene |
| `JL`     | No       | Light chain J gene |

Any argument left as `NA` (the default) is excluded from matching. Columns not referenced by an argument are dropped before matching and do not appear in the output, so keep your own metadata in a separate table and re-join it on `barcodes` afterwards.

### Reference

The reference is also a `data.frame` or `data.table`, one contig per row. Unlike the query, **its columns must use the canonical names below** — `find_publicBCR()` internally prefixes every reference column with `ref_` and looks for `ref_CDRH3`, `ref_VH`, and so on. Do not add the `ref_` prefix yourself.

Sequence columns (these must exist for whichever chains you are matching on):

| Column  | Required | Contents |
| ------- | -------- | -------- |
| `CDRH3` | Yes      | Heavy chain CDR3 amino acid sequence |
| `VH`    | Only if `VH` is set in the call | Heavy chain V gene |
| `JH`    | Only if `JH` is set in the call | Heavy chain J gene |
| `CDRL3` | Only if `CDRL3` is set in the call | Light chain CDR3 amino acid sequence |
| `VL`    | Only if `VL` is set in the call | Light chain V gene |
| `JL`    | Only if `JL` is set in the call | Light chain J gene |

Metadata columns are all optional. They are carried through to the output with a `ref_` prefix so you can trace which reference contig a query cell matched. The bundled databases use:

| Column | Contents |
| ------ | -------- |
| `context` | Antigen context, e.g. `SarsCoV2` |
| `strain` | Strain or variant the contig was tested against, e.g. `SARS-CoV2_WT`. Multiple values are separated by `;` |
| `epitope` | Epitope or antigen sub-domain, e.g. `S2`, `RBD` |
| `org` | `human` or `mouse` |
| `publication` | Source accession or PMID, e.g. `GSE219098` |
| `name` | Contig or clone identifier from the source study |
| `binding` | Logical; whether the contig binds the antigen |

`binding` is the one metadata column that changes behaviour: if the reference contains both `TRUE` and `FALSE` values, non-binders are dropped automatically so that only true binders are used as the reference. A reference with no `binding` column, or with only one value, is used as-is.

A minimal custom reference therefore looks like:

```r
reference <- data.frame(
  VH    = c("IGHV4-59", "IGHV3-66"),
  JH    = c("IGHJ5",    "IGHJ3"),
  CDRH3 = c("CAKGIYSSSSYWFGPW", "CARDHSGHALDIW"),
  name  = c("clone_1", "clone_2")
)
```

### Validation rules

Both query and reference pass through the same checks before matching, and rows failing any of them are silently dropped (the counts removed are reported in the console):

* **CDR3 sequences** must be uppercase amino acid letters only and at least 5 residues long. Sequences containing `*`, `_`, lowercase letters, digits or gaps are removed, as are `NA`, `""`, `"NA"` and `"None"`. The bundled databases use IMGT junction sequences including the leading `C` and trailing `W`/`F`, e.g. `CARDHSGHALDIW`.
* **V and J genes** must start with `IGH`, `IGK` or `IGL`, e.g. `IGHV4-59`, `IGHJ5`, `IGKV1-9`, `IGLJ3`.
* Gene matching is an **exact string comparison**, so query and reference must use the same nomenclature. In particular, strip allele suffixes (`IGHV4-59*01` will not match `IGHV4-59`) and make sure ambiguous calls are formatted consistently.

## How matching works

`find_publicBCR()` runs the following steps in order:

1. Drop invalid V/J genes and CDR3 sequences (< 5 AA) from both the query and the reference.
2. Keep only rows whose V/J gene combination is present in *both* the query and the reference (skipped if no gene arguments are supplied).
3. Compute the CDR3 amino acid distance for every remaining query–reference pair. Hamming distance only compares sequences of equal length and is normalised by the reference CDR3 length; Levenshtein distance handles unequal lengths and is normalised by the longer of the two sequences. Both therefore range from 0 to 1.
4. For each query cell, keep the single reference contig with the **minimum** distance — the CDRH3 distance when matching on heavy chain only, or the mean of the CDRH3 and CDRL3 distances when matching on both.

Setting `dist_method = c("hamming", "levenshtein")` runs both and returns one row per cell per method. Use `ncores` to parallelise, and `output_dir` plus `output_name` to write the result to CSV.

## Output format

The result is a `data.table` with one row per matched query cell per distance method:

| Column | Description |
| ------ | ----------- |
| `barcodes` | Query cell identifier, taken from the query row names |
| `VH`, `JH`, `CDRH3`, ... | The matched query columns, renamed to their canonical names |
| `CDRH3_length`, `CDRL3_length` | CDR3 lengths in amino acids |
| `ref_*` | All columns of the matched reference contig, including its metadata |
| `CDRH3_dist_method`, `CDRL3_dist_method` | `hamming` or `levenshtein` |
| `CDRH3_dist`, `CDRL3_dist` | Normalised CDR3 distance from 0 to 1. `0` is an exact match to the reference sequence |
| `mean_dist` | Mean of the per-chain distances; equals `CDRH3_dist` when only the heavy chain is used |
| `match_method` | Which fields were used for matching, e.g. `VHJHCDRH3` or `VHJHCDRH3_VLJLCDRL3` |

Cells with no reference match at all are absent from the output. The output also carries two internal helper columns, `gene_key` and `rn`, which can safely be ignored or dropped.

## Setting a distance threshold

There is no universal cut-off — how strict to be depends on your antigen, depth of sequencing and tolerance for false positives. As a starting point we recommend a minimum Levenshtein distance of **less than 0.1**. See the [Getting Started vignette](https://hungms.github.io/ClonoMappeR/articles/ClonoMappeR.html), which benchmarks thresholds by comparing known antigen-positive against antigen-negative sorted cells.

```r
ag_specific_bcr <- output %>%
  filter(CDRH3_dist_method == "levenshtein") %>%
  filter(CDRH3_dist < 0.1)

nrow(ag_specific_bcr)
```

## Reference database collection

We are continuously expanding the database, which currently holds 103,928 contigs across seven antigen contexts. `get_reference(context = ...)` accepts `SarsCoV2`, `Vaccinia`, `Tetanus`, `Measles`, `Mumps`, `Hpylori` and `NP`, and can be narrowed further with the `org`, `publication` and `epitope` arguments.

| Context | Organism | Publication | Epitopes | Contigs | DOI |
| ------- | -------- | ----------- | -------- | ------: | --- |
| SarsCoV2 | human | PMID32805021 | N, NTD, RBD, S, S1, S2 | 8,227 | 10.1093/bioinformatics/btaa739 |
| SarsCoV2 | mouse | PMID32805021 | M, N, NTD, RBD, S, S1, S2 | 270 | 10.1093/bioinformatics/btaa739 |
| SarsCoV2 | human | GSE219098 | S2 | 83,802 | 10.1016/j.celrep.2023.112780 |
| SarsCoV2 | human | GSE252959 | S | 524 | 10.1038/s41467-024-48570-0 |
| SarsCoV2 | human | GSE253857 | S | 602 | 10.1038/s41467-024-48570-0 |
| Vaccinia | human | PMID36130603 | B5 | 150 | 10.1016/j.immuni.2022.08.019 |
| Vaccinia | human | PMID40865529 | A35 | 8 | 10.1016/j.cell.2025.08.004 |
| Vaccinia | human | PMID40432083 | D8 | 3 | 10.3390/vaccines13050471 |
| Tetanus | human | GSE252959 | TT | 286 | 10.1038/s41467-024-48570-0 |
| Tetanus | human | GSE253857 | TT | 429 | 10.1038/s41467-024-48570-0 |
| Tetanus | human | PMID35855325 | TT | 37 | 10.1093/nargab/lqac049 |
| Measles | human | PMID42102820 | FE-1a to FE-5, HE-1a to HE-4 | 19 | 10.1016/j.chom.2026.04.010 |
| Measles | human | PMID26187412 | - | 1 | 10.1016/j.immuni.2015.06.016 |
| Mumps | human | PMID26187412 | - | 1 | 10.1016/j.immuni.2015.06.016 |
| Hpylori | human | PMID12117924 | antiurease, unspecified | 4 | 10.1128/IAI.70.8.4158-4164.2002 |
| NP | mouse | GSE154634 | NP-KLH | 3,386 | 10.1038/s41590-021-00936-y |
| NP | mouse | GSE240813 | NP-KLH | 6,179 | 10.1038/s41590-024-01831-y |

Two entries include non-binding (`binding == FALSE`) contigs alongside their binders: SarsCoV2 / GSE219098 (56,147 of its 83,802 contigs) and Vaccinia / PMID40432083 (2 of 3). `find_publicBCR()` drops these automatically, but they can be retained deliberately as a negative control when benchmarking thresholds — see the Getting Started vignette. Every other entry contains binders only.

## Citation

Code was inspired by the following publication in Cell Reports:

```bibtex
@article{LOPESDEASSIS2023112780,
title = {Tracking B cell responses to the SARS-CoV-2 mRNA-1273 vaccine},
journal = {Cell Reports},
volume = {42},
number = {7},
pages = {112780},
year = {2023},
issn = {2211-1247},
doi = {https://doi.org/10.1016/j.celrep.2023.112780},
url = {https://www.sciencedirect.com/science/article/pii/S221112472300791X},
author = {Felipe {Lopes de Assis} and Kenneth B. Hoehn and Xiaozhen Zhang and Lela Kardava and Connor D. Smith and Omar {El Merhebi} and Clarisa M. Buckner and Krittin Trihemasava and Wei Wang and Catherine A. Seamon and Vicky Chen and Paul Schaughency and Foo Cheung and Andrew J. Martins and Chi-I Chiang and Yuxing Li and John S. Tsang and Tae-Wook Chun and Steven H. Kleinstein and Susan Moir},
keywords = {SARS-CoV-2, mRNA vaccine, B cells, immunological memory, single-cell profiling, BCR repertoire}
}
```
