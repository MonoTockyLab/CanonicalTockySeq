# CanonicalTockySeq: Canonical Tocky Analysis for Temporal Gradient on a Transcriptomic Manifold

<a href="[https://monotockylab.github.io/CanonicalTockySeq/](https://www.google.com/search?q=https://monotockylab.github.io/CanonicalTockySeq/)">
<img src="man/figures/CanonicalTockySeq.jpg" align="left" width=100%>
</a>


**Author:** Dr Masahiro Ono

**Date:** 21 February 2026

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

* **CanonicalTockySeq** is available at GitHub: [CanonicalTockySeq](https://www.google.com/search?q=https://github.com/MonoTockyLab/CanonicalTockySeq).

## Package Documentation

The **CanonicalTockySeq** package documentation is available online:

* **Website**: [https://MonoTockyLab.github.io/CanonicalTockySeq/](https://www.google.com/search?q=https://MonoTockyLab.github.io/CanonicalTockySeq/)

## Copyright and Citation Guidelines

### Copyright

All code and original graphical content within the CanonicalTockySeq package, including anime-like characters and logos, are copyrighted by [Masahiro Ono](https://monotockylab.github.io/). A patent application related to the methodologies employed within this package has been filed and is pending. The intellectual property is held under Imperial College London and Masahiro Ono.

### Usage

The CanonicalTockySeq code is available on GitHub and is intended for public viewing, verification, and use related to the associated academic publication. For permissions or inquiries beyond the scope of the license, please contact: <a href="mailto:m.ono@imperial.ac.uk">m.ono@imperial.ac.uk</a>.

### Citing CanonicalTockySeq

If you use the **CanonicalTockySeq** package or any of its components in a scientific publication or in any other public work, please cite the GitHub repository for now. (This section will be updated with specific paper details once published).

**BibTeX Entry:**

```bibtex
@Manual{OnoCanonicalTockySeq,
  title = {CanonicalTockySeq: Canonical Redundancy Analysis for Tocky Differentiation Trajectories},
  author = {Masahiro Ono},
  year = {2026},
  note = {R package version 0.0.1},
  url = {[https://github.com/MonoTockyLab/CanonicalTockySeq](https://github.com/MonoTockyLab/CanonicalTockySeq)}
}
```

### License

**Apache License 2.0**: The CanonicalTockySeq package is licensed under the Apache License, Version 2.0 (the "License"); you may not use this file except in compliance with the License. You may obtain a copy of the License at:

[http://www.apache.org/licenses/LICENSE-2.0](http://www.apache.org/licenses/LICENSE-2.0)

Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the License for the specific language governing permissions and limitations under the License.

### Warranty

This software is distributed WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

#### Why Citation Is Important

Citing software you've used is crucial for acknowledging contributions and ensuring reproducibility, which are critical for scientific progress.

* Giving credit to the developers and researchers who have contributed to the tools you utilize respects and acknowledges their intellectual contributions.
* Proper citations allow other researchers to access the same tools and versions, thus replicating and verifying your scientific results.

Citations are integral to the scientific ecosystem; they help trace the evolution of ideas and enable others to build upon existing research.

We kindly remind our users that **citing software is as important as citing academic articles in maintaining the integrity of the scientific record.**

#### Further Resources

For additional guidance on citation practices and maintaining research integrity, we recommend visiting the [Committee on Publication Ethics (COPE)](https://publicationethics.org/), which offers valuable resources and support for adhering to ethical practices in scholarly publishing.

## The Ono Lab (MonoTockyLab)

<img src="man/figures/MonoLab.jpg" alt="MonoTockyLab" align="center" width="40%">

**The Masahiro Ono Lab (MonoTockyLab)** offers innovative approaches to analyzing omics and flow cytometric data. The lab is particularly well-known for their development of Timer-of-cell-kinetics-and-Activity (**Tocky**) and integrated analysis of immunological data using both experiments and computational analysis.

**Principal Investigator**: Dr. Masahiro Ono, Reader in Immunology at Imperial College London.

Dr. Ono is **the creator and developer of Tocky**. He innovated the transgenic and computational technologies that constitute Tocky.

In 2008, Dr. Ono initiated his pioneering transition from molecular immunology to becoming an **Integrated Experimental and Computational Immunologist**, demonstrating his visionary leadership and pioneering spirit in the development and application of multidimensional analysis and computational methods to address experimental and immunological problems. Tocky represents one of the fusion technologies that Dr. Ono has both created and developed.

Tocky employs the Fluorescent Timer protein to analyze the temporal dynamics of cell activities and development *in vivo*. His lab integrates molecular biology, immunology, and computational analysis to develop novel research tools, thereby enhancing the understanding of immune cell biology.

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
