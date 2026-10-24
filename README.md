# CellProjectedPhenotypes


# Diffusion Cell-Projected Phenotypes (Diffusion CPP)

Diffusion CPP extends the [Cell-Projected Phenotypes (CPP)](https://github.com/KellisLab/CellProjectedPhenotypes) framework by incorporating diffusion geometry to capture nonlinear transcriptional relationships in single-cell RNA-sequencing data.

Rather than constructing cellular neighborhoods directly in PCA space, Diffusion CPP uses diffusion-map coordinates derived from a cosine-based nearest-neighbor graph. These coordinates are integrated into the CPP framework to identify phenotype-associated cellular neighborhoods and calculate continuous cell-level phenotype scores.

## Methodology

The framework consists of five main stages:

1. **Preprocessing:** Normalize gene expression, identify highly variable genes, scale the data, and compute principal components.
2. **Diffusion geometry:** Construct a cosine-based k-nearest-neighbor graph with Gaussian connectivity and calculate diffusion-map coordinates.
3. **Neighborhood construction:** Build overlapping cellular neighborhoods in diffusion space using Milo.
4. **Phenotype association:** Test associations between neighborhood abundance and sample-level phenotypes, adjusting for relevant covariates.
5. **Cell-level scoring:** Project neighborhood association coefficients onto individual cells using distance-weighted aggregation.

## Installation

The implementation uses R and Python, with `reticulate` providing an interface to Scanpy.

### R dependencies

```r
install.packages(c("Seurat", "reticulate", "Matrix", "Rfast"))
```

Additional Bioconductor packages:

- `miloR`
- `SingleCellExperiment`
- `edgeR`

Optional preprocessing dependency:

- `lisi`

### Python dependencies

```bash
pip install numpy scanpy anndata
```

## Functions

### 1. `computeDiffusionCoordinates.R`

**Function:** `computeDiffusionCoordinates()`

Computes diffusion-map coordinates from an existing PCA representation.

The function:

- Constructs a cosine-based k-nearest-neighbor graph.
- Computes Gaussian neighborhood connectivities.
- Calculates diffusion-map eigenvectors and eigenvalues.
- Scales the remaining eigenvectors by their corresponding eigenvalues raised to diffusion time.

**Input:** AnnData object containing PCA coordinates.

**Output:** Cell-by-component matrix of diffusion coordinates.

Default parameters:

- Nearest neighbors: `k = 30`
- Diffusion components: `n_components = 12`
- Diffusion time: `t = 3`

### 2. `buildDiffusionNeighborhoods.R`

**Function:** `buildDiffusionNeighborhoods()`

Constructs overlapping cellular neighborhoods using diffusion coordinates instead of conventional PCA coordinates.

The function integrates the diffusion embedding into Milo and constructs a k-nearest-neighbor graph followed by neighborhood sampling.

**Input:** SingleCellExperiment object and diffusion coordinates.

**Output:** Milo object containing diffusion-based neighborhoods.

### 3. `calculateNeighborhoodCenterMatrix.R`

**Function:** `calculateNeighborhoodCenterMatrix()`

Calculates the center of each Milo neighborhood in diffusion space.

Each neighborhood center is obtained by averaging the diffusion coordinates of its constituent cells.

**Input:** Milo object and cell diffusion coordinates.

**Output:** Neighborhood-by-component matrix of diffusion-space centroids.

### 4. `calculateCellNeighborhoodDistanceMatrix.R`

**Function:** `calculateCellNeighborhoodDistanceMatrix()`

Calculates Euclidean distances between individual cells and their associated neighborhood centers in diffusion space.

Distances are computed in chunks for memory efficiency and restricted to cell-neighborhood memberships.

**Input:** Neighborhood centers, diffusion coordinates, and Milo object.

**Output:** Sparse matrix of cell-neighborhood distances.

### 5. `runNeighborhoodAssociationTest.R`

**Function:** `runNeighborhoodAssociationTest()`

Tests associations between neighborhood cell abundance and sample-level phenotypes using Milo's neighborhood differential abundance framework.

The function supports multiple phenotypes and optional covariate adjustment.

**Input:** Milo object, phenotype variables, covariates, and sample identifiers.

**Output:** Neighborhood-level association results, including log fold changes and statistical significance.

### 6. `calculateCellScores.R`

**Function:** `calculateCellScores()`

Computes continuous cell-level phenotype scores by aggregating neighborhood association coefficients.

For each cell, neighborhood log fold changes are combined using distance-dependent softmax weights, giving greater influence to nearby neighborhood centers.

**Input:** Milo object, neighborhood association results, cell-neighborhood distance matrix, and phenotype names.

**Output:** Cell-level phenotype score matrix.

### 7. `filterAndSubsampleByCellCount.R`

**Function:** `filterAndSubsampleByCellCount()`

Filters samples with insufficient cell counts and subsamples samples with unusually large cell counts to reduce donor-level imbalance.

**Input:** Seurat object containing sample identifiers.

**Output:** Filtered and subsampled Seurat object.

### 8. `removePoorlyMixedCells.R`

**Function:** `removePoorlyMixedCells()`

Uses the Local Inverse Simpson's Index (LISI) to identify cells with insufficient local mixing across a specified metadata variable.

Cells below a user-defined LISI threshold are excluded.

**Input:** Seurat object, cell embedding, metadata variable, and LISI threshold.

**Output:** Filtered Seurat object.

