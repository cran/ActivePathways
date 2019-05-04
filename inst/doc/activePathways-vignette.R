## ---- echo=FALSE---------------------------------------------------------
knitr::opts_chunk$set(warning=FALSE, 
                      message=FALSE, 
                      width=500)
options(max.print=35)

## ------------------------------------------------------------------------
scores <- read.table(system.file('extdata', 'Adenocarcinoma_scores_subset.tsv', package='ActivePathways'), header=TRUE, sep='\t', row.names='Gene')
scores <- as.matrix(scores)
scores

## ------------------------------------------------------------------------
scores[is.na(scores)] <- 1

## ------------------------------------------------------------------------
library(ActivePathways)
gmt.file <- system.file('extdata', 'hsapiens_REAC_subset.gmt', package='ActivePathways')
ActivePathways(scores, gmt.file)

## ------------------------------------------------------------------------
nrow(ActivePathways(scores, gmt.file, return.all=FALSE, significant=0.05))
nrow(ActivePathways(scores, gmt.file, return.all=FALSE, significant=0.1))
nrow(ActivePathways(scores, gmt.file, return.all=TRUE, significant=0.05))

## ------------------------------------------------------------------------
gmt <- read.GMT(gmt.file)
names(gmt[[1]])

# Pretty-print the GMT
gmt[1:3]

# Look at the genes annotated to the first term
gmt[[1]]$genes

# Get the full name of Reactome pathway 2424491
gmt$`REAC:2424491`$name

## ------------------------------------------------------------------------
gmt <- Filter(function(term) length(term$genes) >= 10, gmt)
gmt <- Filter(function(term) length(term$genes) <= 500, gmt)

## ---- eval=FALSE---------------------------------------------------------
#  write.GMT(gmt, 'hsapiens_REAC_subset_filtered.gmt', path = ".")

## ------------------------------------------------------------------------
ActivePathways(scores, gmt)

## ------------------------------------------------------------------------
background <- makeBackground(gmt)
background <- background[background != 'TP53']
ActivePathways(scores, gmt.file, background=background)

## ------------------------------------------------------------------------
scores
merge_p_values(scores, 'Fisher')
merge_p_values(scores, 'Brown')

## ------------------------------------------------------------------------
scores2 <- cbind(scores[, 'CDS'], merge_p_values(scores[, c('X3UTR', 'X5UTR', 'promCore')], 'Brown'))
colnames(scores2) <- c('CDS', 'non_coding')

scores # with separate columns for 3' UTR, 5'UTR & Promoter

scores2 # with a single 'noncoding' column for the merged p-value of 3' UTR, 5'UTR & Promoter

ActivePathways(scores, gmt.file) # results on the dataset with separate columns for 3' UTR, 5'UTR & Promoter

ActivePathways(scores2, gmt.file) # results on the dataset with a single 'noncoding' column for the merged p-value of 3' UTR, 5'UTR & Promoter, reducing the impact of noncoding mutations

## ------------------------------------------------------------------------
res <- ActivePathways(scores, gmt.file)
res

## ------------------------------------------------------------------------
res$overlap[1:3]

## ------------------------------------------------------------------------
result.file <- paste(system.file('extdata', package='ActivePathways'), 'example_write.csv', sep='/')

## ---- eval = FALSE-------------------------------------------------------
#  data.table::fwrite(res, result.file)

## ---- echo=FALSE---------------------------------------------------------
knitr::kable(read.csv(result.file)[1:3, 1:4], caption = "Write results as csv file using data.table::fwrite - first 4 columns shown here")

## ---- eval=FALSE---------------------------------------------------------
#  res <- ActivePathways(scores, gmt.file, cytoscape.file.dir = system.file('extdata', package='ActivePathways'))

## ------------------------------------------------------------------------
files <- paste(system.file('extdata', package='ActivePathways'), 
               c('pathways.txt', 
                 'subgroups.txt',
                 'pathways.gmt', 
                 'legend.pdf'), 
               sep='/')

cat(paste(readLines(files[1])[1:5], collapse='\n'))
cat(paste(readLines(files[2])[1:5], collapse='\n'))
cat(paste(readLines(files[3])[1], collapse='\n'))

