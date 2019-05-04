## ---- echo=FALSE---------------------------------------------------------
knitr::opts_chunk$set(warning=FALSE, 
                      message=FALSE, 
                      width=500)

## ------------------------------------------------------------------------
files <- c(system.file('extdata', 'pathways.txt', package='ActivePathways'),
           system.file('extdata', 'subgroups.txt', package='ActivePathways'),
           system.file('extdata', 'pathways.gmt', package='ActivePathways'),
           system.file('extdata', 'legend.pdf', package ='ActivePathways'))

## ---- eval=FALSE---------------------------------------------------------
#  # Not run.
#  res <- ActivePathways(dat, gmt, cytoscape.file.dir='path/to/results/directory')

## ------------------------------------------------------------------------
cat(paste(readLines(files[1])[1:3], collapse='\n'))
cat(paste(readLines(files[2])[1:5], collapse='\n'))
cat(paste(readLines(files[3])[1], collapse='\n'))

