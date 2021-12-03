## ---- echo=FALSE--------------------------------------------------------------
knitr::opts_chunk$set(warning=FALSE, 
                      message=FALSE, 
                      width=500)
options(max.print=35)

## -----------------------------------------------------------------------------
scores <- read.table(
system.file('extdata', 'Adenocarcinoma_scores_subset.tsv', package = 'ActivePathways'), 
header = TRUE, sep = '\t', row.names = 'Gene')
scores <- as.matrix(scores)
scores

## -----------------------------------------------------------------------------
scores[is.na(scores)] <- 1

## -----------------------------------------------------------------------------
library(ActivePathways)
gmt.file <- system.file('extdata', 'hsapiens_REAC_subset.gmt', package = 'ActivePathways')
ActivePathways(scores, gmt.file)

## -----------------------------------------------------------------------------
nrow(ActivePathways(scores, gmt.file, significant = 0.05))
nrow(ActivePathways(scores, gmt.file, significant = 0.1))

## -----------------------------------------------------------------------------
gmt <- read.GMT(gmt.file)
names(gmt[[1]])

# Pretty-print the GMT
gmt[1:3]

# Look at the genes annotated to the first term
gmt[[1]]$genes

# Get the full name of Reactome pathway 2424491
gmt$`REAC:2424491`$name

## -----------------------------------------------------------------------------
gmt <- Filter(function(term) length(term$genes) >= 10, gmt)
gmt <- Filter(function(term) length(term$genes) <= 500, gmt)

## -----------------------------------------------------------------------------
ActivePathways(scores, gmt)

## -----------------------------------------------------------------------------
ActivePathways(scores, gmt.file, geneset.filter = c(10, 500))

## ---- eval=FALSE--------------------------------------------------------------
#  write.GMT(gmt, 'hsapiens_REAC_subset_filtered.gmt')

## -----------------------------------------------------------------------------
background <- makeBackground(gmt)
background <- background[background != 'TP53']
ActivePathways(scores, gmt.file, background = background)

## -----------------------------------------------------------------------------
sort(merge_p_values(scores, 'Fisher'))[1:5]
sort(merge_p_values(scores, 'Brown'))[1:5]

## -----------------------------------------------------------------------------
scores2 <- cbind(scores[, 'CDS'], merge_p_values(scores[, c('X3UTR', 'X5UTR', 'promCore')], 'Brown'))
colnames(scores2) <- c('CDS', 'non_coding')
scores[c(2179, 1760),]
scores2[c(2179, 1760),]

ActivePathways(scores, gmt.file)
ActivePathways(scores2, gmt.file)

## -----------------------------------------------------------------------------

nrow(ActivePathways(scores, gmt.file))
nrow(ActivePathways(scores, gmt.file, cutoff = 0.01))


## -----------------------------------------------------------------------------

nrow(ActivePathways(scores, gmt.file))
nrow(ActivePathways(scores, gmt.file, correction.method = 'none'))


## -----------------------------------------------------------------------------
res <- ActivePathways(scores, gmt.file)
res

## -----------------------------------------------------------------------------
res$overlap[1:3]

## -----------------------------------------------------------------------------
res[res$term.id == "REAC:422475","evidence"]

## ---- eval = FALSE------------------------------------------------------------
#  result.file <- paste('ActivePathways_results.csv', sep = '/')
#  export_as_CSV (res, result.file) # remove comment to run
#  read.csv(result.file, stringsAsFactors = F)[1:3,]

## ---- eval=FALSE--------------------------------------------------------------
#  result.file <- paste('ActivePathways_results2.txt', sep = '/')
#  data.table::fwrite(res, result.file, sep = '\t', sep2 = c('', ',', ''))
#  cat(paste(readLines(result.file)[1:2], collapse = '\n'))

## ---- eval=FALSE--------------------------------------------------------------
#  res <- ActivePathways(scores, gmt.file, cytoscape.file.tag = "enrichmentMap__")

## -----------------------------------------------------------------------------
files <- c(system.file('extdata', 'enrichmentMap__pathways.txt', package='ActivePathways'),
           system.file('extdata', 'enrichmentMap__subgroups.txt', package='ActivePathways'),
           system.file('extdata', 'enrichmentMap__pathways.gmt', package='ActivePathways'),
           system.file('extdata', 'enrichmentMap__legend.pdf', package='ActivePathways'))

## ---- eval=FALSE--------------------------------------------------------------
#  gmt.file <- system.file('extdata', 'hsapiens_REAC_subset.gmt', package = 'ActivePathways')
#  scores.file <- system.file('extdata', 'Adenocarcinoma_scores_subset.tsv', package = 'ActivePathways')
#  
#  scores <- read.table(scores.file, header = TRUE, sep = '\t', row.names = 'Gene')
#  scores <- as.matrix(scores)
#  scores[is.na(scores)] <- 1
#  
#  res <- ActivePathways(scores, gmt.file, cytoscape.file.tag = "enrichmentMap__")

## -----------------------------------------------------------------------------
cat(paste(readLines(files[1])[1:5], collapse='\n'))
cat(paste(readLines(files[2])[1:5], collapse='\n'))
cat(paste(readLines(files[3])[18:19], collapse='\n'))

