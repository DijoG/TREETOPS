# TREETOPS

The **TREETOPS** package implements a Growing Tree Region (GTR) algorithm for locating treetops in Canopy Height Models (CHMs) derived from low point density LiDAR data.

Designed to work within the **lidR** package ecosystem, **TREETOPS** is intended to be coupled with CHM-based tree segmentation methods to complete an individual tree detection workflow. 

***Update (27.10.2025):*** A parallelized processing option for treetop location has been added, providing dramatically faster processing times.

## Original paper
*A new method for individual treetop detection with low-resolution aerial laser scanned data*

https://doi.org/10.1007/s40808-024-02060-w


## Required packages 

tidyverse, terra, sf, data.table, crayon, future.apply


## Installation

```r
devtools::install_github("DijoG/TREETOPS")
```

## Example
Demonstration of how to use **TREETOPS**. 

### Data preparation

```r
TREETOPS::check_PACKS()
# Either:
> No need, packages already installed.
# or:
> All packages installed.
require(lidR)

# Forest point cloud (low resolution example data)
LASfile <- system.file("extdata", "MixedConifer.laz", package = "lidR")
Alas <- readLAS(LAS, filter = "-drop_z_below 0") 

# Black, white and green color palette for visualizing CHM
bgcol <- function(x)
{
  col = grDevices::colorRampPalette(c("grey1", "white", "forestgreen"))
  return(col(x))
}
```

### Obtain CHM

```r
CHM <- rasterize_canopy(Alas, 0.5, pitfree(subcircle = 0.25))
```

### Compute treetops 

```r
# 1) Original function 
tictoc::tic()
treetops_original <- get_TREETOPS(CHM, min_H = 5, level_increment = 0.2)
tictoc::toc()
# 811.9/60 ~ 13.5 mins

# 2) Optimized sequential 
tictoc::tic()
treetops_optimized <- get_TREETOPS_optimized(CHM, min_H = 5, level_increment = 0.2)
tictoc::toc()
# 752.48/60 ~ 12.5 mins 

# 3) Optimized with parallel processing
# One tile per core (4 cores × 1 tile per core) 
tictoc::tic()
treetops_parallel <- get_TREETOPS_optimized(CHM, 
                                            min_H = 5, 
                                            level_increment = 0.2,
                                            use_parallel = TRUE, 
                                            n_cores = 4, 
                                            cores_tile = 1)
tictoc::toc()
# 180.68/60 ~ 3 mins
```

### Visualize

```r
plot(CHM, main = "CHM 0.5 pitfree ~ treetops", col = bgcol(50), axes = F)
plot(sf::st_geometry(treetops_parallel), add = T, pch = 1, col = "firebrick3")
```
<img align="bottom" src="https://raw.githubusercontent.com/DijoG/storage/main/README/TREETOPS_01.png">

### Finalize treetops (filter)

```r
?TREETOPS::finalize_TREETOPS
fin_treetops <- TREETOPS::finalize_TREETOPS(treetops_parallel, distance = 5, min_H = 5)
```

### Visualize

```r
plot(CHM, main = "CHM 0.5 pitfree ~ filtered treetops", col = bgcol(50), axes = F)
plot(sf::st_geometry(fin_treetops), add = T, pch = 16, col = "firebrick3")
```
<img align="bottom" src="https://raw.githubusercontent.com/DijoG/storage/main/README/TREETOPS_02.png">

