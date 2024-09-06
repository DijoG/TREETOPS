The **TREETOPS** package provides functions utilizing a growing tree region (GTR) treetop identification algorithm using low point density LiDAR data derived Canopy Height Model (CHM).The package is meant to be used within the framework of the lidR package and intends to be coupled with CHM-based tree segmentation methods. 

### Original paper
*A new method for individual treetop detection with low-resolution aerial laser scanned data*

https://doi.org/10.1007/s40808-024-02060-w


### Required packages (make sure they are installed)

tidyverse, terra, sf, data.table, crayon


### Installation

```r
devtools::install_github("DijoG/TREETOPS")
```

# Example
Demonstration of how to use **TREETOPS**. 

### Data preparation

```r
TREETOPS::check_PACKS()
>All packages installed.
require(lidR);require(tidyverse);require(terra)

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

### Compute treetops (takes longer depending on the performance of your computer)

```r
treetops <- TREETOPS::get_TREETOPS(CHM, min_H = 5)
```

### Visualize

```r
plot(CHM, main = "CHM 0.5 pitfree ~ treetops", col = bgcol(50), axes = F)
plot(sf::st_geometry(treetops), add = T, pch = 1, col = "firebrick3")
```
<img align="bottom" src="https://raw.githubusercontent.com/DijoG/storage/main/README/TREETOPS_01.png">

### Finalize treetops (filter)

```r
?TREETOPS::finalize_TREETOPS
fin_treetops <- TREETOPS::finalize_TREETOPS(treetops, distance = 5, min_H = 5)
```

### Visualize

```r
plot(CHM, main = "CHM 0.5 pitfree ~ filtered treetops", col = bgcol(50), axes = F)
plot(sf::st_geometry(fin_treetops), add = T, pch = 16, col = "firebrick3")
```
<img align="bottom" src="https://raw.githubusercontent.com/DijoG/storage/main/README/TREETOPS_02.png">

