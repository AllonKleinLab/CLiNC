# Cell-Lineage-from-Normalized-Covariance (R)

Cell-Lineage-from-Normalized-Covariance (CLiNC) is a method to reconstruct developmental hierarchies from clonal barcoding data. The method is described in REF. Briefly, the model underlying CLiNC assumes that all barcodes are deposited as a synchronous moment in differentiation and that differentiation events are not directly coupled to cell division (as in asymmetric division). 

#### Algorithm overview
The input data is a matrix of barcode counts in across cell types. In principle these counts should represent numbers of cells (as opposed to numbers of sequencing reads). The output is an inferred cell type hierarchy and a list of putative tree violations. The only parameter is the false-discovery rate for detection of conformal symmetry violations (default 5%). The CLiNCs pipeline includes the following steps:

1. Calculate normalized covariance between each pair of cell types
2. Use neighbor-joining to iterative form a cell type hierarchy
3. Identify statistically significant deviations from conformal symmetry
4. Use symmetry violations to infer putative differentiation pathways that violate the hierarchy

Note that step 4 is only available in the python package. 

## Installation ##

Install from github

```
devtools::install_github("AllonKleinLab/Cell-Lineage-from-Normalized-Covariance/clinc_R/clinc")
```

## Usage ##

#### Setup ####

Download or clone this github repository by running in the terminal
```
https://github.com/AllonKleinLab/Cell-Lineage-from-Normalized-Covariance.git
```

Open an R console and load the CLiNC library ```library('clinc')```

Find the full path to the R example directory ```Cell-Lineage-from-Normalized-Covariance/clinc_R/example``` and set this as the current path in R by running ```setwd(<path to example directory>)```

#### Load data and set parameters ####

```
input_data_path <- 'barcode_counts.tsv'
output_directory <- 'example_output'
symmetry_violation_FDR <- 0.05

make_output_dir(output_directory)
barcode_counts <- load_data(input_data_path)
celltype_names <- names(barcode_counts)
```

#### Calculate and plot normalized covariance ####

```
X <- get_normalized_covariance(barcode_counts)
plot_normalized_covariance(output_directory, X)
```
The plot will be saved as a pdf in the output directory


#### Build and plot the cell type hierarchy ####

```
hierarchy_data <- build_hierarchy(barcode_counts)
parent_map <- hierarchy_data[[1]]
node_groups <- hierarchy_data[[2]]
plot_hierarchy(parent_map, celltype_names)
```

#### Calculate and plot symmetry violations ####

```
violations_data <- detect_symmetry_violations(barcode_counts, parent_map, symmetry_violation_FDR)
violations <- violations_data[[1]]
plot_data <- violations_data[[2]]
plot_violations(output_directory, plot_data)
```
