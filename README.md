# PerturbAtlas

![Version](https://img.shields.io/badge/version-v1.0.0-blue)
![Python](https://img.shields.io/badge/python-3.9%2B-blue)
![License](https://img.shields.io/badge/license-MIT-green)
![Status](https://img.shields.io/badge/status-stable-success)
> **A lightweight Python framework for perturbation signature retrieval
> and biologically interpretable program inference from transcriptomic
> gene sets.**

**Current release:** `PerturbAtlas-K562-v1.0`
Current status: **Version 1.0.0 (Initial Public Release)**

PerturbAtlas is designed to behave like a **biologically informed search
engine** rather than a black-box classifier. Given a transcriptomic gene
set, it retrieves biologically similar perturbation signatures from a
curated reference atlas, predicts the dominant biological program, and
quantifies the strength of evidence supporting each prediction.

------------------------------------------------------------------------

## Overview

Traditional enrichment analyses identify pathways or biological
processes associated with a gene set, but they do not directly answer
**which perturbation is most likely to have generated the observed
transcriptional program**.

PerturbAtlas addresses this problem by treating perturbation signatures
as searchable biological programs. Instead of returning only enriched
pathways or functional annotations, it retrieves the most similar
perturbation signatures, predicts the dominant biological program
represented by the query, and provides an interpretable evidence score
describing the reliability of the prediction.

The first public release, **PerturbAtlas-K562-v1.0**, was constructed
from **105 curated K562 cell line Perturb-seq perturbation signatures** organized
into three biologically meaningful transcriptional programs.

------------------------------------------------------------------------

## Who is PerturbAtlas for?

PerturbAtlas is intended for researchers working with:

-   Bulk RNA-seq differential expression analyses
-   CRISPR perturbation experiments
-   Perturb-seq datasets
-   Gene regulatory network studies
-   Functional genomics
-   Transcriptomic signature interpretation

------------------------------------------------------------------------

## Design Philosophy

PerturbAtlas was developed around four principles:

-   Biological interpretability
-   Transparency
-   Simplicity
-   Extensibility

Every prediction can be traced back to retrieved perturbation
signatures. The package favors transparent biological reasoning over
opaque machine-learning predictions.

------------------------------------------------------------------------

## Key Features

  -----------------------------------------------------------------------
  Feature                             Description
  ----------------------------------- -----------------------------------
  Biological program classification   Predicts the dominant
                                      transcriptional program represented
                                      by a query.

  Perturbation retrieval              Retrieves the most similar
                                      perturbation signatures.

  Evidence scoring                    Combines multiple
                                      validation-derived metrics into a
                                      single evidence score.

  Biological interpretation           Produces concise human-readable
                                      summaries.

  Master regulator prioritization     Highlights regulatory perturbations
                                      with strong influence.

  Fully interpretable                 No black-box inference.
  -----------------------------------------------------------------------

------------------------------------------------------------------------

## Current Atlas

  Property                  Value
  ------------------------- ------------------------
  Atlas                     PerturbAtlas-K562-v1.0
  Organism                  *Homo sapiens*
  Cell line                 K562
  Dataset                   K562 Perturb-seq
  Reference perturbations   105
  Biological modules        3
  Reference genes           1540

Modules:

-   Myeloid differentiation
-   Erythroid differentiation
-   Developmental / cell-state transition

------------------------------------------------------------------------

## Installation

### Install from GitHub

```bash
git clone https://github.com/UCSOLMAZ/PerturbAtlas.git
cd PerturbAtlas
python -m pip install -e .
```

### Install from PyPI *(coming soon)*

```bash
pip install perturbatlas
```

*PyPI distribution is planned for a future release.*

------------------------------------------------------------------------

## Quick Start
The following example loads the released K562 atlas, queries a small gene set, and prints a biological interpretation.
``` python
from perturbatlas import K562Atlas

atlas = K562Atlas.load("data/PerturbAtlas-K562-v1.0.pkl")

result = atlas.query([
    "SPI1","CEBPA","CEBPE","CSF3R","MPO","ELANE"
])

result.summary()
result.interpret()
```

------------------------------------------------------------------------

## Example Output

```text
==================================================
PerturbAtlas-K562-v1.0 Query Result
==================================================

Predicted module:
  myeloid

Module scores:
  myeloid         0.0860
  developmental   0.0278
  erythroid       0.0076

Evidence score:
   Very High (1.00)

Top perturbations:
  1. CEBPA      | recall=1.00 | precision=0.13 | module=myeloid


```
## Query Workflow

```text
Input gene set
      |
      v
Module scoring
      |
      v
Perturbation retrieval
      |
      v
Evidence scoring
      |
      v
Biological interpretation
```

------------------------------------------------------------------------

## Evidence Score
```
  Component                      Weight
  ---------------------------- --------
  Module separation                 40%
  Top perturbation agreement        40%
  Mean retrieval recall             20%

     Score Interpretation
  -------- ----------------
     >=%0.85 Very High
     >=%0.70 High
     >=%0.50 Moderate
     <0.50 Low

Low-evidence predictions are intentionally reported conservatively.

------------------------------------------------------------------------
```
## Validation

### Leave-One-Perturbation-Out (LOPO)

**Purpose:** Recover known perturbation programs.

**Result:** 90.7% module accuracy, 98.1% Top-3 retrieval.

**Conclusion:** Excellent recovery of known perturbational programs.

### Random Gene Controls

**Purpose:** Evaluate specificity.

**Result:** Mean evidence score of approximately 0.06 across 400 random 
gene-set queries with no high-confidence predictions.

### Query-size Robustness

Recommendation: use 20 or more informative genes for optimal
performance.

### Noise Robustness

Twenty informative genes mixed with eighty random genes retained
96.3% module prediction accuracy.

### Cross-module Specificity

Evidence decreases for mixed biological states and increases again as
one program becomes dominant.

------------------------------------------------------------------------

## API

### K562Atlas

-   `query()`
-   `info()`
-   `save()`
-   `load()`
-   `modules()`
-   `perturbations()`

### AtlasQueryResult

-   `summary()`
-   `interpret()`
-   `top_perturbations()`
-   `to_csv()`
-   `evidence()`
-   `evidence_score()`
-   `confidence_score()`

------------------------------------------------------------------------

## Repository Structure

```text
PerturbAtlas/
|
|-- perturbatlas/
|   |-- __init__.py
|   |-- atlas.py
|   |-- result.py
|
|-- data/
|   |-- PerturbAtlas-K562-v1.0.pkl
|
|-- examples/
|   |-- QuickStart.ipynb
|   |-- ExampleQueries.ipynb
|   |-- Benchmark.ipynb
|
|-- README.md
|-- LICENSE
|-- pyproject.toml
|-- requirements.txt
`-- .gitignore
```

------------------------------------------------------------------------

## Data Source

PerturbAtlas-K562-v1.0 was constructed using publicly available K562 Perturb-seq data. The current atlas is based on curated perturbation signatures derived from the following study:

**Norman TM, Horlbeck MA, Replogle JM, Ge AY, Xu A, Jost M, Gilbert LA, Weissman JS.**

*Exploring genetic interaction manifolds constructed from rich single-cell phenotypes.*

**Science.** 2019;365(6455):786-793.

DOI: https://doi.org/10.1126/science.aax4438

The processed dataset used during atlas construction was obtained from:

https://www.kaggle.com/datasets/alexandervc/scrnaseq-crisprperturbseq-normanweissman

We gratefully acknowledge the authors of the original Perturb-seq study for making these data publicly available.


------------------------------------------------------------------------
## Citation

If you use PerturbAtlas in your research, please cite both the software and the reference dataset used to construct the released atlas.

**PerturbAtlas**

Ufuk Solmaz.

*PerturbAtlas: A lightweight Python framework for perturbation signature retrieval and biologically interpretable program inference from transcriptomic gene sets.*

GitHub repository (current release)

*https://github.com/UCSOLMAZ/PerturbAtlas*

**Reference dataset**

Norman TM, Horlbeck MA, Replogle JM, Ge AY, Xu A, Jost M, Gilbert LA, Weissman JS.

*Exploring genetic interaction manifolds constructed from rich single-cell phenotypes.*

Science. 2019;365(6455):786-793.

DOI: https://doi.org/10.1126/science.aax4438

------------------------------------------------------------------------

## Roadmap

The current release (**PerturbAtlas-K562-v1.0**) represents the first public version of the framework. Planned future developments include:

### Near-term

- Additional cell-type perturbation atlases
- Expanded benchmarking and validation
- Drug perturbation atlas support
- Improved documentation and tutorials
- PyPI distribution (`pip install perturbatlas`)

### Long-term

- Single-cell perturbation atlases
- Cross-cell consensus perturbation atlases
- Automatic atlas download and version management
- Interactive visualization utilities
- Integration with additional transcriptomic resources


------------------------------------------------------------------------
## License

PerturbAtlas is released under the **MIT License**.

You are free to use, modify, and distribute the software in accordance with the terms of the license.

See the `LICENSE` file for complete license information.

------------------------------------------------------------------------
## Author

**Ufuk C. Solmaz**

Developer and maintainer of **PerturbAtlas**

GitHub: https://github.com/UCSOLMAZ/

For questions, suggestions, bug reports, or collaborations, please open an issue or discussion on the GitHub repository.

------------------------------------------------------------------------

## Contributing

Contributions, feature requests, bug reports, and suggestions are welcome.

If you encounter a bug or have ideas for improving PerturbAtlas, please open an issue or submit a pull request.

Community contributions are greatly appreciated.
