# CLiNC (Cell Lineage from Normalized Covariance)


CLiNC (Cell Lineage from Normalized Covariance) is a method to reconstruct developmental hierarchies from clonal barcoding data. The method is described in REF. Briefly, the model underlying CLiNC assumes that all barcodes are deposited as a synchronous moment in differentiation and that differentiation events are not directly coupled to cell division (as in asymmetric division). 

#### Algorithm overview
The input data is a matrix of barcode counts in across cell types. In principle these counts should represent numbers of cells (as opposed to numbers of sequencing reads). The output is an inferred cell type hierarchy and a list of putative tree violations. The only parameter is the false-discovery rate for detection of conformal symmetry violations (default 5%). The CLiNCs pipeline includes the following steps:

1. Calculate normalized covariance between each pair of cell types
2. Use neighbor-joining to iterative form a cell type hierarchy
3. Identify statistically significant deviations from conformal symmetry
4. Use symmetry violations to infer putative differentiation pathways that violate the hierarchy
5. Detect possible tree distortions and apply a post-hoc correction


## Usage

CLiNC is available as a [python package](https://github.com/AllonKleinLab/CLiNC/tree/master/clinc_python) and an [R package](https://github.com/AllonKleinLab/CLiNC/tree/master/clinc_R). Detection of cross-hierarchy differentiation pathways is only available in python.
