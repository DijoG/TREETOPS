
#' Getting height levels I
#'
#' Function for getting height bin (numeric vector).
#' @param CHM_g chm from las, SpatRaster
#' @param level_increment numeric, height bin 
#' @param min_H numeric, under this value no trees are reasonable
#' @return numeric vector 
#' @export
get_HB <- function(CHM_g, level_increment, min_H) {
  return(seq(min(values(CHM_g), na.rm = T) + min_H, max(values(CHM_g), na.rm = T), level_increment) %>%
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
    #if (length(cell$focal_mean) < 4) {
    #  empty_r = c(1,2,3)
    #}
    
    ra[[h]] = empty_r
    names(ra[[h]]) = height_bin[h]
  }
  
  RA = ra[sapply(ra, function(x) class(x)[1] == "SpatRaster")]
  RAout = rast(RA)
  
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
  
  # Patches of second, marker of first
  marker =
    data.table::data.table(extract(GTR_ccomponent, GTR_marker)) %>%
    rename(marker = ID) %>%
    mutate_all(~replace(., is.nan(.), 0)) %>%
    filter(patches != 0)
  
  P_unique =
    marker$patches %>% 
    unique
  
  # Treetops as markers
  TTL = list()
  for (i in seq_along(P_unique)) {
    
    # i-th patch
    patchi = mask(CHM_g, ifel(GTR_ccomponent == P_unique[i], 1, NA)) %>% trim()
    
    # i-th marker
    mark =
      marker %>%
      filter(patches == P_unique[i]) %>%
      select(marker) %>%
      pull
    
    mark_pts =
      GTR_marker %>%
      filter(marker %in% mark)
    
    # add Z (height) and treeID infos
    ttops =
      terra::extract(patchi, mark_pts)
    ttops_Z =
      mark_pts %>%
      mutate(Z = ttops$focal_mean,
             treeID = ttops$ID)
    
    TTL[[i]] = ttops_Z
  }
  return(bind_rows(TTL))
}

#' Getting GTR and FETR
#' 
#' Function for binarizing patches into GROWING TREE REGIONS (GTR ~ 1) and FIRST EMERGING TREE REGION (FETR ~ 0), used in 'get_CCC()'.
#' @param lmo_seg segregated patches defined in get_CCC, SpatRaster stack
#' @param new_tr new tree region, SpatRaster (value 1)
#' @return SpatRaster (value 1 GTR pixels, value 0 non-GTR pixels)
#' @export
get_GTR <- function(lmo_seg, new_tr) {
  
  for (i in 1:nlyr(lmo_seg)) {
    
    cell_n_lmo = length(cells(ifel(lmo_seg[[i]] == 1, 1, NA)))
    cell_n_lmo_newtr = cells(ifel(lmo_seg[[i]] == 1, 1, NA)) %in% cells(new_tr) %>% sum
    
    if (cell_n_lmo_newtr != cell_n_lmo & cell_n_lmo_newtr > 0) {
      # 1 for GTR
      lmo_seg[[i]] = ifel(lmo_seg[[i]] == 1, 1, NA)
    } else {
      # 0 for FETR
      lmo_seg[[i]] = ifel(lmo_seg[[i]] == 1, 0, NA)
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
  m = mask(level_rasterMINone, level_raster,  maskvalues = 1, updatevalue = 2)
  if (1 %in% unique(values(m))) {
    m = ifel(m != 1, NA, m)
  } else {
    stop(str_c("\n meaning:\n There is no growing tree from ", names(level_raster), " to ", names(level_rasterMINone), "\n Incease level increment or tile (i.e. chm) size!\n OR ~ delete ", names(level_raster), " from height bin!:)"))
  }
  
  # Clumps labeling
  clumps_lmo = patches(level_rasterMINone, zeroAsNA = TRUE, allowGaps = FALSE)
  clumps_lmo_seg = segregate(clumps_lmo)
  
  # 3b) GTR
  gtr = get_GTR(clumps_lmo_seg, m)
  clumps_gtr = patches(ifel(gtr != 1, NA, gtr), zeroAsNA = TRUE, allowGaps = FALSE)
  
  # Centroid for each of the patches of GTR (clumps_gtr)
  out_l = matrix(NA, nrow = max(values(clumps_gtr), na.rm = T), ncol = 3)
  for (i in 1:max(values(clumps_gtr), na.rm = T)) {
    cmean = colMeans(xyFromCell(clumps_gtr, which(values(clumps_gtr) == i)))
    out_l[i, 1] = cmean[[1]]
    out_l[i, 2] = cmean[[2]]
    out_l[i, 3] = as.double(level)
  }
  
  pts_l = vect(cbind(out_l[,1],out_l[,2]), crs = crs(level_raster)) # projection
  values(pts_l) = out_l[,3]
  
  pts_l = sf::st_as_sf(pts_l)
  names(pts_l)[1] = "Z_level"
  pts_l$marker = 1:nrow(pts_l)
  
  return(list("GTR_marker" = pts_l,
              "GTR_ccomponent" = clumps_gtr))
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
  height_bin <- get_HB(CHM_g, level_increment, min_H)
  
  # 2)
  RAS = get_LEVELz(CHM_g, height_bin)
  
  cat(crayon::cyan("_____ Getting markers at height level:\n"))
  
  LOUT = list()
  for (ini in seq_along(names(RAS))) {
    
    iniplus = ini+1
    cat(crayon::cyan("_____", names(RAS)[iniplus], "\r"))
    
    first = names(RAS)[ini]
    second = names(RAS)[iniplus]
    third = names(RAS)[iniplus+1]
    
    # 3a)
    marks = get_MARKER(CHM_g,
                       # 3c)
                       get_CCC(RAS[[second]], RAS[[third]], third)$GTR_ccomponent,
                       get_CCC(RAS[[first]], RAS[[second]], second)$GTR_marker)
    
    
    if (nrow(marks) == 0) {
      marks = c(1,2,3)
    }
    
    LOUT[[ini]] = marks
    
    if (ini == length(names(RAS))-2) {
      break
    }
    
  }
  
  LOUT = LOUT[sapply(LOUT, function(x) class(x)[1] == "sf")]
  LOUT = bind_rows(LOUT)
  
  cat(crayon::silver("\n_____ Done\n"))
  return(LOUT)
  
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
    arrange(desc(Z_level), desc(Z)) %>%
    mutate(treeID = 1:nrow(.))
  ter_TREETOPS = vect(sf_TREETOPS)
  
  
  # Removing treetops falling in distance ~ radius (m) 
  nearID = data.table::data.table(nearby(ter_TREETOPS, distance = distance))
  
  # Filtering height
  max_H = ifelse(is.numeric(max_H), max_H, max(sf_TREETOPS$Z))
  TTout =
    sf_TREETOPS %>%
    slice(-unique(nearID$to)) %>%
    filter(Z >= min_H & Z <= max_H) %>%
    add_column(Z_range = str_c(as.character(min_H), " - ", as.character(max_H)), .before = "marker")
  
  return(TTout %>%
           mutate(treeID = 1:nrow(.)))
  
  
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
    mutate(ID = 1:nrow(.)) %>%
    select(ID) %>%
    add_column(Z = extract(CHM_g, vect(refTT))$focal_mean, .after = 1)
  
  testTT = 
    testTT %>%
    add_column(ID = 1:nrow(testTT), .before = 1) %>%
    select(ID, Z)
  
  treeID = refTT$ID
  
  cat("_____ There are", max(treeID), "reference trees.\n")
  
  dfout = list()
  for (tID in treeID) {
    
    cat("_____ treeID:", tID, "\r")
    
    # Getting the 3 nearest test trees (to the given reference tree)
    nearID = data.table::data.table(nearby(vect(refTT %>%
                                                  filter(treeID == tID)), 
                                           vect(testTT), 
                                           k = 3)) %>%
      pivot_longer(cols = k1:k3, values_to = "id_test", names_to = "centroids") %>%
      rename("id_ref" = id) %>%
      data.table::data.table()
    
    if (nrow(nearID) >= 1) {
      fromTT = 
        refTT %>%
        filter(treeID == tID)
      toTT = 
        testTT %>%
        slice(nearID$id_test)
    } else {
      next
    }
    
    # Getting distance between reference tree and test trees
    D = sf::st_distance(map_dfr(seq_len(nrow(toTT)), ~bind_rows(fromTT)), 
                        toTT, 
                        by_element = T) %>%
      as.numeric()
    
    # Obtaining tree coordinates and Z (height) profile
    L = list()
    V = list()
    for (i in 1:nrow(toTT)) {
      L[[i]] = sf::st_linestring(rbind(as.numeric(sf::st_coordinates(fromTT)[-3]),
                                       as.numeric(sf::st_coordinates(toTT[i,]))))
      V[[i]] = terra::extract(CHM_g, vect(L[[i]]))$focal_mean
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
  return(bind_rows(dfout) %>%
           mutate(diff_h_reftest = abs(Z_start - Z_stop)) %>%
           group_by(ref_ID) %>%
           filter(distance == min(distance)) %>%
           ungroup() %>% 
           mutate(match = ifelse(diff_h_reftest <= eval_maxHD & distance <= eval_maxD, 1, 0),
                  n_test = nrow(testTT)) %>%
           data.table::data.table())
}
