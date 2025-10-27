usethis::use_pipe()

#' Check packages
#' @return installs missing packages
#' @export
check_PACKS <- function() {
  required = c("terra", "tidyverse", "sf", "crayon", "data.table", "future.apply", "stringr")
  installable = required[!(required %in% installed.packages()[,"Package"])]
  if(length(installable) > 0) {
    install.packages(installable)
    cat("All packages installed.")
  }
  else {
    cat("No need, packages already installed.\n")
  }
}

#' Getting height levels I
#'
#' Function for getting height bin (numeric vector).
#' @param CHM_g chm from las, SpatRaster
#' @param level_increment numeric, height bin 
#' @param min_H numeric, under this value no trees are reasonable
#' @return numeric vector 
#' @export
get_HB <- function(CHM_g, level_increment, min_H) {
  return(seq(min(terra::values(CHM_g), na.rm = T) + min_H, max(terra::values(CHM_g), na.rm = T), level_increment) %>%
           round(2) %>%
           rev)
}

#' Getting height levels II
#' 
#' Function for getting height levels (stacked SpatRaster), uses output of 'get_HB()'.
#' @param CHM_g chm from las, SpatRaster
#' @param height_bin numeric height bin vector returned by 'get_HB()' 
#' @return stacked SpatRaster
#' @export
get_LEVELz <- function(CHM_g, height_bin) {
  
  ra = list()
  
  for (h in seq_along(height_bin)) {
    
    empty_r = CHM_g
    empty_r = terra::ifel(empty_r >= height_bin[h], 1, 0)
    cell = terra::cells(empty_r, 1)
    if (length(cell$Z) < 4) {
      empty_r = c(1,2,3)
    }
    
    ra[[h]] = empty_r
    names(ra[[h]]) = height_bin[h]
  }
  
  RA = ra[sapply(ra, function(x) class(x)[1] == "SpatRaster")]
  RAout = terra::rast(RA)
  
  cat(crayon::silver("\n_____ Height bin ready\n"))
  return(RAout)
}

#' Getting markers
#' 
#' Function for getting markers, its ids and heights (Z), used in 'get_TREETOPS()'. 
#' @param CHM_g chm from las, SpatRaster
#' @param GTR_ccomponent /result of get_CCC/ GTR patches from layer2 to layer3, level raster is layer3, SpatRaster (values n patches)
#' @param GTR_marker /result of get_CCC/ GTR markers from layer1 to layer2, level raster is layer2, sf spatial points
#' @return sf object
#' @export
get_MARKER <- function(CHM_g, GTR_ccomponent, GTR_marker) {
  
  # Check for NULL inputs
  if (is.null(GTR_ccomponent) || is.null(GTR_marker) || nrow(GTR_marker) == 0) {
    return(NULL)
  }
  
  # Patches of second, marker of first
  marker = tryCatch({
    extract_result = terra::extract(GTR_ccomponent, GTR_marker)
    data.table::data.table(extract_result) %>%
      dplyr::rename(marker = ID) %>%
      dplyr::mutate_all(~replace(., is.nan(.), 0)) %>%
      dplyr::filter(patches != 0)
  }, error = function(e) {
    return(NULL)
  })
  
  if (is.null(marker) || nrow(marker) == 0) {
    return(NULL)
  }
  
  P_unique = marker$patches %>% unique
  
  # Treetops as markers
  TTL = list()
  for (i in seq_along(P_unique)) {
    
    # i-th patch - FIXED: Better error handling
    patchi = tryCatch({
      # Create mask for this patch
      patch_mask = terra::ifel(GTR_ccomponent == P_unique[i], 1, NA)
      # Apply mask to CHM
      terra::mask(CHM_g, patch_mask) %>% 
        terra::trim()
    }, error = function(e) {
      NULL
    })
    
    if (is.null(patchi) || terra::ncell(patchi) == 0) {
      next
    }
    
    # i-th marker
    mark = marker %>%
      dplyr::filter(patches == P_unique[i]) %>%
      dplyr::select(marker) %>%
      dplyr::pull()
    
    mark_pts = GTR_marker %>%
      dplyr::filter(marker %in% mark)
    
    if (nrow(mark_pts) == 0) {
      next
    }
    
    # add Z (height) and treeID infos
    ttops = tryCatch({
      terra::extract(patchi, mark_pts)
    }, error = function(e) {
      NULL
    })
    
    if (is.null(ttops)) {
      next
    }
    
    ttops_Z = mark_pts %>%
      dplyr::mutate(Z = ttops$Z,
                    treeID = ttops$ID)
    
    TTL[[i]] = ttops_Z
  }
  
  if (length(TTL) == 0) {
    return(NULL)
  }
  
  return(dplyr::bind_rows(TTL))
}

#' Getting GTR and FETR
#' 
#' Function for binarizing patches into GROWING TREE REGIONS (GTR ~ 1) and FIRST EMERGING TREE REGION (FETR ~ 0), used in 'get_CCC()'.
#' @param lmo_seg segregated patches defined in get_CCC, SpatRaster stack
#' @param new_tr new tree region, SpatRaster (value 1)
#' @return SpatRaster (value 1 GTR pixels, value 0 non-GTR pixels)
#' @export
get_GTR <- function(lmo_seg, new_tr) {
  
  for (i in 1:terra::nlyr(lmo_seg)) {
    
    # FIX: Handle different return types from terra::cells()
    cells_lmo = terra::cells(terra::ifel(lmo_seg[[i]] == 1, 1, NA))
    cell_n_lmo = if (is.list(cells_lmo)) {
      if (!is.null(cells_lmo$cell)) length(cells_lmo$cell) else 0
    } else {
      length(cells_lmo)
    }
    
    cells_newtr = terra::cells(new_tr)
    cell_n_lmo_newtr = if (is.list(cells_lmo) && is.list(cells_newtr)) {
      if (!is.null(cells_lmo$cell) && !is.null(cells_newtr$cell)) {
        sum(cells_lmo$cell %in% cells_newtr$cell)
      } else {
        0
      }
    } else if (!is.list(cells_lmo) && !is.list(cells_newtr)) {
      sum(cells_lmo %in% cells_newtr)
    } else {
      0
    }
    
    if (cell_n_lmo_newtr != cell_n_lmo & cell_n_lmo_newtr > 0) {
      # 1 for GTR
      lmo_seg[[i]] = terra::ifel(lmo_seg[[i]] == 1, 1, NA)
    } else {
      # 0 for FETR
      lmo_seg[[i]] = terra::ifel(lmo_seg[[i]] == 1, 0, NA)
    }
  }
  return(max(lmo_seg, na.rm = T))
}

#' Getting patches
#'
#' Function for getting the patches (clumps) and centroids of GTR, used in 'get_MARKER()' that is embedded in 'get_TREETOPS()'.
#' @param level_raster /result of get_LEVELz()/ i-th layer of SpatRaster stack
#' @param level_rasterMINone /result of get_LEVELz()/ i - 1 layer of SpatRaster stack
#' @param level level layer, always level_rasterMINone, however the algorithm needs this extra parameter
#' @return list (sf object and terra SpatRaster)
#' @export
get_CCC <- function(level_raster, level_rasterMINone, level = NULL) {
  
  # First emerged tree region 
  m = terra::mask(level_rasterMINone, level_raster,  maskvalues = 1, updatevalue = 2)
  if (!1 %in% unique(terra::values(m))) {
    stop(stringr::str_c("\n meaning:\n There is no growing tree from ", names(level_raster), " to ", names(level_rasterMINone), "\n Increase level increment or tile (i.e. chm) size!\n OR ~ delete ", names(level_raster), " from height bin!:)"))
  }
  
  m = terra::ifel(m != 1, NA, m)
  
  # Clumps labeling
  clumps_lmo = terra::patches(level_rasterMINone, zeroAsNA = TRUE, allowGaps = FALSE)
  clumps_lmo_seg = terra::segregate(clumps_lmo)
  
  # 3b) GTR
  gtr = get_GTR(clumps_lmo_seg, m)
  clumps_gtr = terra::patches(terra::ifel(gtr != 1, NA, gtr), zeroAsNA = TRUE, allowGaps = FALSE)
  
  # Get valid patches - FIXED: Use proper method to get unique patch values
  patch_values = terra::values(clumps_gtr)
  valid_patches = unique(patch_values)
  valid_patches = valid_patches[!is.na(valid_patches) & valid_patches != 0]
  
  if (length(valid_patches) == 0) {
    return(list(GTR_marker = NULL, GTR_ccomponent = clumps_gtr))
  }
  
  # Centroid for each of the patches of GTR (clumps_gtr)
  out_l = matrix(NA, nrow = length(valid_patches), ncol = 3)
  valid_count = 0
  
  for (i in seq_along(valid_patches)) {
    patch_id = valid_patches[i]
    
    # FIXED: Use proper method to get cells for each patch
    patch_cells = which(terra::values(clumps_gtr) == patch_id)
    
    if (length(patch_cells) > 0) {
      valid_count = valid_count + 1
      coords = terra::xyFromCell(clumps_gtr, patch_cells)
      cmean = colMeans(coords)
      out_l[valid_count, 1] = cmean[[1]]
      out_l[valid_count, 2] = cmean[[2]]
      out_l[valid_count, 3] = as.double(level)
    }
  }
  
  # Remove unused rows
  if (valid_count == 0) {
    return(list(GTR_marker = NULL, GTR_ccomponent = clumps_gtr))
  }
  
  out_l = out_l[1:valid_count, , drop = FALSE]
  
  # Create vector points
  pts_l = terra::vect(cbind(out_l[,1], out_l[,2]), crs = terra::crs(level_raster))
  terra::values(pts_l) = out_l[,3]
  
  pts_l = sf::st_as_sf(pts_l)
  names(pts_l)[1] = "Z_level"
  pts_l$marker = 1:nrow(pts_l)
  
  return(list(GTR_marker = pts_l, GTR_ccomponent = clumps_gtr))
}

#' Parallel processing with controlled tile distribution
#' @param CHM_g SpatRaster object, raster chm derived from las
#' @param min_H numeric, minimum (Z) value of chm
#' @param level_increment numeric, level cutting increment
#' @param n_cores number of parallel workers
#' @param cores_tile number of tiles per core (default = 1)
#' @return sf object
#' @export
get_TREETOPS_parallel_controlled <- function(CHM_g, min_H, level_increment = 0.2, 
                                             n_cores = NULL, cores_tile = 1) {
  
  if (!requireNamespace("future.apply", quietly = TRUE)) {
    stop("Please install future.apply: install.packages('future.apply')")
  }
  
  if (is.null(n_cores)) {
    n_cores = parallel::detectCores() - 1
  }
  
  # Calculate total tiles
  total_tiles = n_cores * cores_tile
  
  cat(crayon::cyan("_____ Creating", total_tiles, "tiles (",n_cores, "cores ×", cores_tile, "tiles per core)\n"))
  
  # Write CHM to temporary file to avoid serialization issues
  temp_dir = tempdir()
  chm_file = file.path(temp_dir, "chm_temp.tif")
  terra::writeRaster(CHM_g, chm_file, overwrite = TRUE)
  
  # Get CHM info for tiling
  nrows = terra::nrow(CHM_g)
  ncols = terra::ncol(CHM_g)
  xmin = terra::xmin(CHM_g)
  xmax = terra::xmax(CHM_g)
  ymin = terra::ymin(CHM_g)
  ymax = terra::ymax(CHM_g)
  xres = terra::xres(CHM_g)
  yres = terra::yres(CHM_g)
  crs_info = terra::crs(CHM_g)
  
  # Simple row-based splitting
  rows_per_tile = ceiling(nrows / total_tiles)
  
  # Set up parallel plan
  future::plan(future::multisession, workers = n_cores)
  
  # Process each tile in parallel
  tile_results = future.apply::future_lapply(1:total_tiles, function(tile_idx) {
    
    tryCatch({
      # Calculate row range for this tile
      start_row = (tile_idx - 1) * rows_per_tile + 1
      end_row = min(tile_idx * rows_per_tile, nrows)
      
      # Skip if tile is too small
      if ((end_row - start_row) < 10) {
        return(NULL)
      }
      
      # Calculate Y coordinates for this tile (note: Y decreases with row increase)
      tile_ymax = ymax - (start_row - 1) * yres
      tile_ymin = ymax - end_row * yres
      
      # Create spatial extent for this tile
      tile_ext = terra::ext(xmin, xmax, tile_ymin, tile_ymax)
      
      # Load CHM from file and crop to tile (each worker loads independently)
      chm_loaded = terra::rast(chm_file)
      tile_CHM = terra::crop(chm_loaded, tile_ext)
      
      # Skip if tile has no data
      if (all(is.na(terra::values(tile_CHM)))) {
        return(NULL)
      }
      
      cat(crayon::blue("Processing tile", tile_idx, "of", total_tiles, 
                       paste0("(rows ", start_row, "-", end_row, ")\n")))
      
      # Use ORIGINAL get_TREETOPS for each tile
      tile_treetops = get_TREETOPS(
        CHM_g = tile_CHM,
        min_H = min_H,
        level_increment = level_increment
      )
      
      return(tile_treetops)
      
    }, error = function(e) {
      return(NULL)
    })
    
  }, future.seed = TRUE)
  
  # Reset to sequential
  future::plan(future::sequential)
  
  # Clean up temp file
  unlink(chm_file)
  
  # Combine results from all tiles with better error handling
  valid_results = list()
  for (i in seq_along(tile_results)) {
    if (!is.null(tile_results[[i]]) && inherits(tile_results[[i]], "sf") && nrow(tile_results[[i]]) > 0) {
      valid_results[[length(valid_results) + 1]] = tile_results[[i]]
    }
  }
  
  if (length(valid_results) == 0) {
    warning("No treetops found in any tile.")
    return(NULL)
  }
  
  # Combine all treetops
  all_treetops = dplyr::bind_rows(valid_results)
  
  cat(crayon::silver("_____ Combined", length(valid_results), "tiles, found", nrow(all_treetops), "treetops\n"))
  
  return(all_treetops)
}

#' Getting treetops
#' 
#' MAIN FUNCTION for obtaining treetops.
#' @param CHM_g SpatRaster object, raster chm derived from las
#' @param min_H numeric, minimum (Z) value of chm, higher values can be treetops
#' @param level_increment numeric, level cutting increment (default = 0.2m)
#' @return sf object
#' @export
get_TREETOPS <- function(CHM_g, min_H, level_increment = 0.2) {
  
  # 1)
  height_bin = get_HB(CHM_g, level_increment, min_H)
  
  # 2)
  RAS = get_LEVELz(CHM_g, height_bin)
  
  cat(crayon::cyan("_____ Getting markers at height level:\n"))
  
  LOUT = list()
  for (ini in seq_along(names(RAS))) {
    
    iniplus = ini+1
    
    # Check if we have enough layers
    if (iniplus + 1 > length(names(RAS))) {
      break
    }
    
    cat(crayon::cyan("_____", names(RAS)[iniplus], "\r"))
    
    first = names(RAS)[ini]
    second = names(RAS)[iniplus]
    third = names(RAS)[iniplus+1]
    
    # 3a) Get markers with error handling
    marks = tryCatch({
      get_MARKER(CHM_g,
                 get_CCC(RAS[[second]], RAS[[third]], third)$GTR_ccomponent,
                 get_CCC(RAS[[first]], RAS[[second]], second)$GTR_marker)
    }, error = function(e) {
      NULL
    })
    
    # FIXED: Handle NULL marks properly
    if (!is.null(marks) && nrow(marks) > 0) {
      LOUT[[ini]] = marks
    }
    
    if (ini == length(names(RAS))-2) {
      break
    }
  }
  
  # Filter out NULL results
  LOUT = LOUT[!sapply(LOUT, is.null)]
  
  if (length(LOUT) == 0) {
    warning("No treetops found at any level.")
    return(NULL)
  }
  
  LOUT = dplyr::bind_rows(LOUT)
  
  cat(crayon::silver("\n_____ Done\n"))
  return(LOUT)
}

#' Memory-optimized version of get_TREETOPS
#' 
#' MAIN FUNCTION for obtaining treetops with memory optimization.
#' @param CHM_g SpatRaster object, raster chm derived from las
#' @param min_H numeric, minimum (Z) value of chm, higher values can be treetops
#' @param level_increment numeric, level cutting increment (default = 0.2m)
#' @param use_parallel logical, whether to use parallel processing (default = FALSE)
#' @param n_cores number of cores for parallel processing (default = NULL)
#' @param cores_tile number of tiles per core (default = 1)
#' @return sf object
#' @export
get_TREETOPS_optimized <- function(CHM_g, min_H, 
                                   level_increment = 0.2, 
                                   use_parallel = FALSE, 
                                   n_cores = NULL, 
                                   cores_tile = 1) {
  
  # If parallel is requested, use the controlled tiled parallel version
  if (use_parallel) {
    cat(crayon::cyan("_____ Using controlled tiled parallel processing\n"))
    return(get_TREETOPS_parallel_controlled(
      CHM_g = CHM_g,
      min_H = min_H,
      level_increment = level_increment,
      n_cores = n_cores,
      cores_tile = cores_tile
    ))
  }
  
  # Convert to memory if not already
  if (!terra::inMemory(CHM_g)) {
    CHM_g = terra::rast(CHM_g)  # Ensure it's loaded
  }
  
  # 1) Get height bins
  height_bin = get_HB(CHM_g, level_increment, min_H)
  
  cat(crayon::cyan("_____ Processing height levels (optimized):\n"))
  
  LOUT = list()
  level_count = length(height_bin) - 2
  
  for (ini in seq_len(level_count)) {
    iniplus = ini + 1
    cat(crayon::cyan("_____ Level", ini, "of", level_count, "- Height:", height_bin[iniplus], "\r"))
    
    # Create level rasters on demand instead of storing all
    level_second = terra::ifel(CHM_g >= height_bin[iniplus], 1, 0)
    level_third = terra::ifel(CHM_g >= height_bin[iniplus + 1], 1, 0)
    level_first = terra::ifel(CHM_g >= height_bin[ini], 1, 0)
    
    # Use the original CCC function (proven to work)
    ccc_result_second_third = get_CCC(level_second, level_third, height_bin[iniplus + 1])
    ccc_result_first_second = get_CCC(level_first, level_second, height_bin[iniplus])
    
    if (!is.null(ccc_result_second_third$GTR_ccomponent) && !is.null(ccc_result_first_second$GTR_marker)) {
      marks = get_MARKER(CHM_g, 
                         ccc_result_second_third$GTR_ccomponent, 
                         ccc_result_first_second$GTR_marker)
      
      if (!is.null(marks) && nrow(marks) > 0) {
        LOUT[[ini]] = marks
      }
    }
    
    # Clean up
    rm(level_second, level_third, level_first, ccc_result_second_third, ccc_result_first_second)
    gc()
  }
  
  # Combine results
  LOUT = LOUT[!sapply(LOUT, is.null)]
  
  if (length(LOUT) == 0) {
    warning("No treetops found. Consider adjusting parameters.")
    return(NULL)
  }
  
  result = dplyr::bind_rows(LOUT)
  cat(crayon::silver("\n_____ Done (optimized)\n"))
  return(result)
}

#' Getting final treetops
#' 
#' MAIN FUNCTION for distance based treetop filtering using the output of 'get_TREETOPS()'.
#' @param sf_TREETOPS sf object (output of 'get_TREETOPS()')
#' @param distance numeric, distance between treetops in meters
#' @param min_H numeric, minimum height of trees
#' @param max_H numeric, maximum height of trees, if not set (default ~ NULL)
#' the maximum height value of the input sf object (sf_TREETOPS)
#' @return sf object
#' @export
finalize_TREETOPS <- function(sf_TREETOPS, distance, min_H, max_H = NULL) {
  
  sf_TREETOPS =
    sf_TREETOPS %>%
    dplyr::arrange(desc(Z_level), desc(Z)) %>%
    dplyr::mutate(treeID = 1:nrow(.))
  ter_TREETOPS = terra::vect(sf_TREETOPS)
  
  # Removing treetops falling in distance ~ radius (m) 
  nearID = data.table::data.table(terra::nearby(ter_TREETOPS, distance = distance))
  
  # Filtering height
  max_H = ifelse(is.numeric(max_H), max_H, max(sf_TREETOPS$Z))
  TTout =
    sf_TREETOPS %>%
    dplyr::slice(-unique(nearID$to)) %>%
    dplyr::filter(Z >= min_H & Z <= max_H) %>%
    tibble::add_column(Z_range = stringr::str_c(as.character(min_H), " - ", as.character(max_H)), .before = "marker")
  
  return(TTout %>%
           dplyr::mutate(treeID = 1:nrow(.)))
}

#' Helper function for obtaining metrics
#' 
#' Used in 'M_match_TREETOPS()'. 
#' @param profile ~  height (Z) profile between reference tree and test treetops
#' @param Zstart height of reference tree
#' @param Zstop height of test treetop
#' @return data frame
#' @export
M_diff_PROFILE <- function(profile, Zstart, Zstop) {
  
  min_p = c()
  mean_p = c()
  max_p = c()
  diff_height = c()
  for (p in seq_along(profile)) {
    
    if (Zstart %in% profile[[p]][length(profile[[p]])] &
        Zstop[p] %in% profile[[p]][length(profile[[p]])-1]) {
      profile[[p]] = c(profile[[p]][length(profile[[p]])], profile[[p]])
      profile[[p]] = profile[[p]][-length(profile[[p]])]
    }
    
    if (Zstart %in% profile[[p]][length(profile[[p]])-1] &
        Zstop[p] %in% profile[[p]][length(profile[[p]])-2]) {
      profile[[p]] = c(profile[[p]][length(profile[[p]])-1], profile[[p]])
      profile[[p]] = profile[[p]][-c(length(profile[[p]]), length(profile[[p]])-1)]
    }
    
    min_p[p] = min(profile[[p]])
    mean_p[p] = mean(profile[[p]])
    max_p[p] = max(profile[[p]])
    diff_height[p] = Zstart - min_p[p]
  }
  return(data.frame(diff_h = diff_height,
                    min_h = min_p,
                    mean_h = mean_p,
                    max_h = max_p)) 
}

#' Matching treetops 
#'
#' The function finds the 3 nearest test tree points (knn), matches one (height difference and distance) and returns:
#' - height difference between observed and test tree,
#' - distance bewteen observed and test tree.
#' @param refTT sf object, reference trees 
#' @param testTT sf object, either GTR or VWF results, extracted treetops
#' @param CHM_g SpatRAster, Gaussian smoothed chm generated from las
#' @param eval_maxHD numeric or double, minimum height difference (m) btw. reference and test trees (default = 2)
#' @param eval_maxD numeric or double, minimum distance btw. reference (m) and test trees (default = 4)
#' @return data.table
#' @export
M_match_TREETOPS <- function(refTT, testTT, CHM_g, eval_maxHD = 2, eval_maxD = 4) {
  
  refTT =
    refTT %>%
    dplyr::mutate(ID = 1:nrow(.)) %>%
    dplyr::select(ID) %>%
    tibble::add_column(Z = terra::extract(CHM_g, terra::vect(refTT))$focal_mean, .after = 1)
  
  testTT = 
    testTT %>%
    tibble::add_column(ID = 1:nrow(testTT), .before = 1) %>%
    dplyr::select(ID, Z)
  
  treeID = refTT$ID
  
  cat("_____ There are", max(treeID), "reference trees.\n")
  
  dfout = list()
  for (tID in treeID) {
    
    cat("_____ treeID:", tID, "\r")
    
    # Getting the 3 nearest test trees (to the given reference tree)
    nearID = data.table::data.table(terra::nearby(terra::vect(refTT %>%
                                                                dplyr::filter(treeID == tID)), 
                                                  terra::vect(testTT), 
                                                  k = 3)) %>%
      tidyr::pivot_longer(cols = k1:k3, values_to = "id_test", names_to = "centroids") %>%
      dplyr::rename("id_ref" = id) %>%
      data.table::data.table()
    
    if (nrow(nearID) >= 1) {
      fromTT = 
        refTT %>%
        dplyr::filter(treeID == tID)
      toTT = 
        testTT %>%
        dplyr::slice(nearID$id_test)
    } else {
      next
    }
    
    # Getting distance between reference tree and test trees
    D = sf::st_distance(map_dfr(seq_len(nrow(toTT)), ~dplyr::bind_rows(fromTT)), 
                        toTT, 
                        by_element = T) %>%
      as.numeric()
    
    # Obtaining tree coordinates and Z (height) profile
    L = list()
    V = list()
    for (i in 1:nrow(toTT)) {
      L[[i]] = sf::st_linestring(rbind(as.numeric(sf::st_coordinates(fromTT)[-3]),
                                       as.numeric(sf::st_coordinates(toTT[i,]))))
      V[[i]] = terra::extract(CHM_g, terra::vect(L[[i]]))$focal_mean
    }
    
    Zstart = fromTT$Z
    Zstop = toTT$Z
    
    # Calculating Z (height) difference
    P = M_diff_PROFILE(V, Zstart, Zstop)
    
    # Saving results
    dfout[[tID]] = 
      data.table::data.table(distance = D,
                             diff_h = P$diff_h,
                             ref_ID = tID,
                             test_ID = toTT$ID,
                             min_h = P$min_h,
                             mean_h = P$mean_h,
                             max_h = P$max_h,
                             Z_start = Zstart,
                             Z_stop = Zstop) 
  }
  
  # Filtering on max height (Z) difference (eval_maxHD) combined with distance (eval_maxD)
  # Output
  cat("\n_____ Done")
  return(dplyr::bind_rows(dfout) %>%
           dplyr::mutate(diff_h_reftest = abs(Z_start - Z_stop)) %>%
           dplyr::group_by(ref_ID) %>%
           dplyr::filter(distance == min(distance)) %>%
           dplyr::ungroup() %>% 
           dplyr::mutate(match = ifelse(diff_h_reftest <= eval_maxHD & distance <= eval_maxD, 1, 0),
                         n_test = nrow(testTT)) %>%
           data.table::data.table())
}