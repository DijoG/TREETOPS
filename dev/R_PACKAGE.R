

setwd("D:/TREETOPS")
devtools::load_all(".")

# Make a new R project 
# Check wd
getwd()

# Roxygenise (man:))
roxygen2::roxygenise()

require(TREETOPS)

?get_HB
?get_MARKER
?get_LEVELz

# install TREETOPS
Sys.setenv(R_REMOTES_STANDALONE = "true")
devtools::install_github("DijoG/TREETOPS", build = F)


