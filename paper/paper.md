---
title: 'crystract: A Crystallography Package in R for .cif Data Processing'
authors:
  - name: Don Ngo
    orcid: 0009-0001-2779-2146
    affiliation: 1
  - name: Julia M. Huebner
    orcid: 0000-0003-2048-6629
    affiliation: 2
      - name: Marc Spitzner
    affiliation: 3
  - name: Shaunna M. Morrison
    orcid: 0000-0002-1712-8057
    affiliation: 4
    - name: Anirudh Prabhu
    orcid: 0000-0002-9921-6084
    affiliation: 1
affiliations:
 - name: Earth and Planets Laboratory, Carnegie Institution for Science, Washington, DC, USA
   index: 1
 - name: Technische Universität Dresden, Dresden, Germany
   index: 2
 - name: Elbe Flugzeugwerke GmbH, Dresden, Germany
   index: 3
   -name: Department of Earth and Planetary Sciences, Rutgers University, Piscataway, NJ, USA
   index: 4
date: 15 September 2025
bibliography: paper.bib
---

# Summary

The Crystallographic Information File (CIF) is the standard format for disseminating crystal structure data, yet parsing and analyzing these files for large-scale computational and statistical analysis is a significant bottleneck for research in the chemical and material sciences [@Hall1991]. `crystract` is an R package designed to provide an efficient, open-source solution for the batch processing and statistical analysis of CIF files. The package streamlines the extraction of metadata, unit cell parameters, atomic coordinates, and symmetry operations for a singular or hundreds of CIF files at the same time. From the information stored in the CIF file, our package can determine the first-neighbor bonding environment around each symmetrically independent atom in a unit cell, with all interatomic distances and bond angles, while propagating experimental uncertainties. Furthermore, a comprehensive workflow for efficient extraction and processing of these data is provided, centered around the calculation of average interatomic distances, as identified as a key parameter for a streamlined comparison of structural features of different compounds. One of the most important features within this workflow is the package's ability to process positional or occupational disorder. To handle positional disorder, a filtering function to eliminate non-physical distances is provided. Occupational disorder is taken into account via the possibility of calculating an occupancy-weighted average interatomic distance. Additionally, filtering functions to calculate average distances, only including user-specified elements or atomic positions, are available. This paper outlines the architecture and core functionalities of `crystract`, demonstrating its utility with a practical example.

# Statement of need

Despite the standardization provided by the CIF format, significant practical barriers remain for researchers aiming to perform high-throughput computational analysis, particularly in batch. The first barrier is technical: CIF files, while standardized, often exhibit syntactic variations or errors depending on their originating software or laboratory, which can cause simplistic parsers to fail. The second barrier is conceptual: a CIF file does not typically contain an explicit list of all atoms in the unit cell. Instead, it reports the unique atoms in the asymmetric unit and a set of symmetry operations. A complete structural analysis, therefore, requires the correct application of these symmetry operations to generate the full unit cell and its atomic contents—a non-trivial, error-prone computational task that must be handled by specialized software. This complexity often forces researchers into a fragmented and inefficient workflow, piecing together disparate tools for data validation, structure generation, geometric analysis, and final statistical modeling.

A number of excellent software tools have been developed to address some of these challenges, yet a comprehensive review reveals a specific and critical gap. The most mature ecosystem for computational materials science currently resides in Python, where `pymatgen` stands as a powerful and widely adopted standard [@Ong2013]. This extensive library offers robust CIF parsing and a host of advanced analysis functions, including the calculation of interatomic distances and the prediction of bonding environments using sophisticated, data-mining approaches like the CrystalNN algorithm [@Pan2020]. While `pymatgen` offers robust CIF parsing and advanced analysis, it has limitations for a complete, statistically rigorous workflow: it does not natively compute bond angles, nor does it programmatically propagate the experimental uncertainties reported in CIFs for its derived geometric quantities. Furthermore, batch processing requires the user to write custom scripts to loop over files. Filtering presents another challenge: although `pymatgen` includes a function to filter for minimum and maximum distances, they have to be specified one element at a time, one file at a time, by the user. This makes it cumbersome to use `pymatgen` for large scale datasets and studies.

Other specialized Python tools like `cifkit`—a Python based command-line tool—and its companion software, `cif-bond-analyzer` (CBA)—designed for the singular task of exporting bond lists—provide fast, lightweight, and native batch processing for geometric analysis. However, their scope is intentionally limited to interatomic distances and bond analysis; they also do not provide native functions for calculating bond angles or propagating experimental uncertainties.

Other specialized tools, such as the CCDC's Mercury [@Groom2016] or the IUCr's enCIFer, are indispensable for interactive 3D visualization and formal CIF syntax validation, respectively. However, their primary design as graphical user interface (GUI) applications makes them ill-suited for the automated, scriptable, and reproducible workflows required for modern data science without the possibility of handling a high throughput.

Within the R ecosystem, the landscape is sparse. The `cry` package provides basic statistics for crystallography computation and falls short of large-scale research and ML based applications [@Roveda2017]. Its design, centered on a custom S3 object system, is not optimized for the high-throughput, data-frame-centric workflows required for large-scale statistical analysis. `cry` is limited to the analysis of crystallographic parameters and diffraction data from individual files. It is not equipped for the geometric analysis of atomic structures, as it is unable to apply symmetry operations to generate a full unit cell from asymmetric coordinates; therefore, unable to calculate the resulting interatomic distances and angles, or handle the structural disorder commonly found in real materials.

To our knowledge, no package or software—whether in Python, as GUIs, or within R—provides a single, integrated solution capable of combining the automated batch processing of large collections of CIF files and the systematic propagation of experimental uncertainties, while providing additional features indispensable for the structural analysis of large datasets. Such a research software landscape forces researchers into creating fragmented and inefficient workflows, piecing together disparate tools for structure generation, geometric analysis, and finally statistical modeling.

| Task | CIFkit | CBA (CIF Bond analyzer) | pymatgen | cry | crystract |
|------------|------------|------------|------------|------------|------------|
| Read/parse singular CIF file | Yes | No, uses cifcit for this | Partially, via importers like CifParser | Yes | Yes |
| Batch / high-throughput processing of many CIFs | Yes | Yes | Yes, but certain tasks require manual scripting | No | Yes |
| Supercell / unit cell generation / lattice operations | Yes, can generate unit cell and supercell via +/- 1 shifts | Yes, when computing minimum bond lengths for site | Yes, structure operations, transformation, supercell, etc. | No | Yes, can generate unit cell and supercell via +/- 1 shifts |
| Coordination number determination | Yes | Yes | Yes | No | Yes |
| Calculation of interatomic distances | Yes | Yes | Yes | No | Yes |
| Calculation of bond angles | No | No | No | No | Yes |
| Error propagation | No | No | No | No | Yes |
| Handling of occupational disorder | Yes | No, uses cifkit | Partial, can read occupancies, but does not use them | No | Yes |
| Handling of structural disorder | No | No | No | No | Yes, via filtering function |
| Filtering for specific atoms or crystallographic sites | Yes | Yes | No | No | Yes |
| Calculation of weighted average accounting for disorder. | No | No | No | No | Yes |
| Output to multiple formats | No outputs of extracted data, only figures such as histograms or visualizations | No outputs of extracted data, only figures such as histograms or visualizations | No | No | Yes |

: Overview of functionalities provided by existing packages in Python and R.

# Crystract

To address these needs, we have developed "`crystract`", an open-source R package designed to provide a seamless, robust, and statistically-minded workflow for crystallographic analysis. `crystract` provides an end-to-end toolkit that operates entirely within the R environment. Its primary contributions are fourfold.

First, it provides a robust and efficient engine for parsing and processing large batches of CIF files. It is designed from the ground up around R's data-centric paradigm, directly presenting all extracted and calculated data in tidy data frames [@Wickham2014] ready for immediate manipulation and analysis with the wider R ecosystem.

Second, it offers comprehensive geometric analysis. This output includes not only the CIF file's core metadata but also a rich set of derived attributes essential for crystallographic research: a complete list of atomic coordinates after the application of symmetry operations to generate the full cell, all interatomic distances based on predicted bonded pairs for direct neighbors within the coordination sphere of an individual atom using the Crystal NN algorithm [@Pan2020], and bond angles—a feature not natively available in any other command-line batch-processing tools.

Third, and most uniquely, `crystract` introduces a capability largely absent in other programmatic tools: the rigorous propagation of experimental uncertainties from the CIF through all derived geometric quantities. This feature, grounded in standard error propagation theory [@Ku1966] facilitates a more sophisticated and honest statistical treatment of structural data, allowing researchers to quantify the confidence in their calculated results.

Fourth, as the most valuable feature for the user community, it provides a suite of integrated filtering functions that are essential for high-throughput analysis and handling complex structures, thus exceeding the capabilities offered by currently existing software tools. Our workflow is centered around the calculation of average interatomic distances, as a key parameter in the structural comparison of different compounds. The average atomic distance can be calculated for all symmetrically independent atoms or a subset of these by prior application of the filtering function based on the user-specified element or atomic site. Furthermore, occupational or positional disorder can be handled via a filtering function to exclude non-physical distances in the process of calculating the average distance. This filtering is based on the automatic recognition of the elements in a structure and the calculation of an expected interatomic distance from their covalent radii [@Emsley1998; @Pyykko2009] or a user-specified list of atomic radii. Partial site occupation occurs in many real compounds and is indispensable for deriving structure-property trends or predicting new functional materials [@Jakob2025]. If one simply calculates interatomic distances from the coordinates given in the CIF file without accounting for partial occupation, one obtains non-physical distances between atoms that cannot coexist in the same local configuration. Such artifacts would distort any statistical measure by including interactions that never occur in reality. While taking these features native to real crystal structures into account, a weighted average distance can be calculated, either for all individual atoms or a user-defined subset by the application of the available filtering functions for user-specified atoms or crystallographic sites.

# Implications

The availability of packages like "`crystract`", "`pymatgen`", "`cifkit`", and others bring the field of crystallography, mineralogy, and materials science a step closer to realizing the transformative potential of data-driven science. `crystract` is one part of a larger effort made to develop AI methods to perform comparative studies across hundreds or even thousands of structures thus providing researchers with the foundations to derive overarching structure-property relationships of minerals and materials to ultimately employ known compounds for new applications or predict new materials.

`crystract` creates a launchpad for integrating crystallographic data into machine learning pipelines for a variety of application areas. Our package's capability to handle positional or occupational disorder, propagate uncertainties, and generate comprehensive data outputs provide an excellent feature set for predictive modeling efforts.

Finally, as an open source package, "`crystract`" will be a community-driven, transparent resource that invites extensions and improvements from other researchers who want to use crystallographic data in their own scientific explorations.

# Acknowledgements

The authors would like to thank Michael Baitinger, Robert T. Downs, Jolyon Ralph, and Xiaogang Ma for their discussions on crystallography, and cyberinfrastructure development. D.N. has been supported by the Earth and Planetary Science Interdisciplinary Internship at Carnegie Science (a National Science Foundation REU). Additionally, funding and support for this project was provided by the Carnegie Institution for Science and a private foundation.

# References
