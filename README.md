# CanonicalTockySeq: Canonical Tocky Analysis for Temporal Gradient on a Transcriptomic Manifold

<a href="https://github.com/MonoTockyLab/CanonicalTockySeq">
<img src="man/figures/CanonicalTockySeq.jpg" align="center" width=100%>
</a>


**Author:** Dr Masahiro Ono

**Date:** 13 March 2026

## Introduction

Experimentally resolving temporal transcriptional dynamics at single-cell resolution *in vivo* is a major challenge, as standard scRNA-seq only provides cross-sectional "snapshots". **CanonicalTockySeq** bridges this gap by integrating a molecular clock of T-cell receptor (TCR) signalling—based on the **Nr4a3-Tocky** Fluorescent Timer—with scRNA-seq to establish an experimentally anchored temporal reference.

By using landmark Tocky fractions (New, Persistent, and Arrested) as biological ground truths, the package constructs a transcriptomic manifold in canonical space with unique **conical geometry**. This allows for the decoupling of **temporal progression** (geodesic angle) from **signalling strength** (radial intensity).


## Key Features

The **CanonicalTockySeq** framework provides a robust methodology for capturing continuous, endogenous temporal resolution anchored to biological history:

#### **CanonicalTockySeq (RDA Architecture)**
* Implements a supervised **Canonical Redundancy Analysis (RDA)** designed for single-cell Tocky data.


* Constrains high-dimensional transcriptomic space using reference Tocky signatures to derive a model-driven path of differentiation.


* Projects individual single cells onto this structure to visualize their state while preserving biological heterogeneity.


#### **GradientTockySeq**

* Utilizes **Piecewise Spherical Linear Interpolation (SLERP)** to connect landmark biplot vectors (New → Persistent → Arrested).


* Generates a curved, **geodesic path** along the transcriptomic hypersphere, preventing artifacts inherent to standard polynomial trajectories.


* Maps cells to a 0–90° **Tocky Time** coordinate (angle) and a normalized **Tocky Intensity** (magnitude).


#### **Temporal Cascade Analysis**
* Orders genes by their peak activation timing to reveal sequential transcriptional programs.

* Facilitates the identification of immediate-early committers, intermediate effectors, and terminal exhaustion markers.


## Availability

* **CanonicalTockySeq** is available at GitHub: [CanonicalTockySeq](https://github.com/MonoTockyLab/CanonicalTockySeq).

## Package Documentation

The **CanonicalTockySeq** package documentation is available online:

* **Website**: [https://MonoTockyLab.github.io/CanonicalTockySeq/](https://MonoTockyLab.github.io/CanonicalTockySeq)

## Citation

If you use `CanonicalTockySeq` in your research, please cite the relevant publication and the R package.

### bioRxiv preprint

```bibtex
@article{Ono2026CanonicalTocky,
  author = {Ono, Masahiro and others},
  title = {Temporal Mechanisms of T-Cell Fate Decisions under Immune Checkpoint Blockade Resolved by CanonicalTockySeq},
  elocation-id = {2026.03.10.710825},
  year = {2026},
  doi = {10.64898/2026.03.10.710825},
  publisher = {Cold Spring Harbor Laboratory},
  journal = {bioRxiv},
  url = {https://www.biorxiv.org/content/10.64898/2026.03.10.710825v1}
}
```

### R package

```bibtex
@Manual{OnoCanonicalTockySeq,
  title = {CanonicalTockySeq: Canonical Tocky Analysis for Temporal Gradient on a Transcriptomic Manifold},
  author = {Masahiro Ono},
  year = {2026},
  note = {R package version 0.1.0.9000},
  url = {https://github.com/MonoTockyLab/CanonicalTockySeq}
}
```

## License

`CanonicalTockySeq` is licensed under the Apache License 2.0. See the `LICENSE` file for details.

## Copyright and intellectual property

Copyright © Masahiro Ono.

Original graphical content in this repository, including characters and logos, is protected by copyright unless otherwise stated.

A patent application relating to aspects of the methodology has been filed.

## Contact

For questions about the package, please contact: [m.ono@imperial.ac.uk](mailto:m.ono@imperial.ac.uk)


## The Ono Lab (MonoTockyLab)

<img src="man/figures/MonoLab.jpg" alt="MonoTockyLab" align="center" width="40%">

**The Masahiro Ono Lab (MonoTockyLab)** develops experimental and computational approaches to study immune cell dynamics, with a particular focus on the temporal regulation of gene expression in T cells.

The lab is known for the development of **Tocky** (*Timer of cell kinetics and activity*), a platform that uses Fluorescent Timer proteins to analyse transcriptional and signalling dynamics *in vivo* at single-cell resolution. Our research integrates mouse genetics, immunology, flow cytometry, single-cell omics, and computational modelling.

Current research directions include:

- cancer immunology and immunotherapy
- temporal mechanisms of T cell activation, differentiation, and tolerance
- **Foxp3 transcriptional dynamics** and their regulation in vivo
- computational methods for time-resolved single-cell analysis, including **CanonicalTockySeq**

**Principal Investigator**: Dr Masahiro Ono, Reader in Immunology at Imperial College London.

Dr Ono is the creator of **Tocky**, spanning both its transgenic reporter systems and associated analytical frameworks.



## Contact and More

**Email**:
<a href="mailto:m.ono@imperial.ac.uk">
<img src="https://upload.wikimedia.org/wikipedia/commons/e/ec/Circle-icons-mail.svg" alt="Email" width="10%">
</a>

**Personal Homepage**:
<a href="http://monotockylab.github.io">
<img src="man/figures/MonoLab.jpg" alt="MonoTockyLab Homepage" align="center" width="30%"/>
</a>

**GitHub**:
<a href="https://github.com/MonoTockyLab">
<img src="https://github.githubassets.com/images/modules/logos_page/GitHub-Mark.png" alt="GitHub" align="center" width="70" height="70"/>
</a>

**Twitter**:
<a href="https://twitter.com/MonoTockyLab">
<img src="https://upload.wikimedia.org/wikipedia/commons/6/6f/Logo_of_Twitter.svg" alt="Twitter" align="center" width="50" height="50"/>
</a>

**Professional Homepage**: [Imperial College London - Masahiro Ono](https://www.imperial.ac.uk/people/m.ono)
