## ----echo=FALSE---------------------------------------------------------------
knitr::opts_chunk$set(warning=FALSE, 
                      message=FALSE, 
                      width=500)
options(max.print=35)
library("ggplot2")
library("data.table")

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
gmt_file <- system.file('extdata', 'hsapiens_REAC_subset.gmt', package = 'ActivePathways')
ActivePathways(scores, gmt_file)

## -----------------------------------------------------------------------------
fname_data_matrix <- system.file('extdata', 
		'Differential_expression_rna_protein.tsv',
		package = 'ActivePathways')

pvals_FCs <- read.table(fname_data_matrix, header = TRUE, sep = '\t')

example_genes <- c('ACTN4','PIK3R4','PPIL1','NELFE','LUZP1','ITGB2')
pvals_FCs[pvals_FCs$gene %in% example_genes,]

## -----------------------------------------------------------------------------
pval_matrix <- data.frame(
		row.names = pvals_FCs$gene, 
		rna = pvals_FCs$rna_pval, 
		protein = pvals_FCs$protein_pval)

pval_matrix <- as.matrix(pval_matrix)
pval_matrix[is.na(pval_matrix)] <- 1

pval_matrix[example_genes,]

## -----------------------------------------------------------------------------
dir_matrix <- data.frame(
		row.names = pvals_FCs$gene, 
		rna = pvals_FCs$rna_log2fc, 
		protein = pvals_FCs$protein_log2fc)

dir_matrix <- as.matrix(dir_matrix)
dir_matrix <- sign(dir_matrix)
dir_matrix[is.na(dir_matrix)] <- 0

dir_matrix[example_genes,]

## -----------------------------------------------------------------------------
constraints_vector <- c(1,1)
constraints_vector

# constraints_vector <- c(1,-1)

## -----------------------------------------------------------------------------
directional_merged_pvals <- merge_p_values(pval_matrix, 
		method = "DPM", dir_matrix, constraints_vector)

merged_pvals <- merge_p_values(pval_matrix, method = "Brown")

sort(merged_pvals)[1:5]
sort(directional_merged_pvals)[1:5]

## -----------------------------------------------------------------------------
pvals_FCs[pvals_FCs$gene == "PIK3R4",]

pval_matrix["PIK3R4",]
dir_matrix["PIK3R4",]

merged_pvals["PIK3R4"]
directional_merged_pvals["PIK3R4"]

## -----------------------------------------------------------------------------
lineplot_df <- data.frame(original = -log10(merged_pvals),
			  modified = -log10(directional_merged_pvals))

ggplot(lineplot_df) +
	geom_point(size = 2.4, shape = 19,
		aes(original, modified,
		    color = ifelse(original <= -log10(0.05),"gray",
                                    ifelse(modified > -log10(0.05),"#1F449C","#F05039")))) +
	labs(title = "",
		 x ="Merged -log10(P)",
		 y = "Directional Merged -log10(P)") + 
            geom_hline(yintercept = 1.301, linetype = "dashed",
		       col = 'black', size = 0.5) +
            geom_vline(xintercept = 1.301, linetype = "dashed",
		       col = "black", size = 0.5) + 
            geom_abline(size = 0.5, slope = 1,intercept = 0) +
	    scale_color_identity()

## -----------------------------------------------------------------------------
constraints_vector <- c(-1,1)

constraints_vector <- c(1,-1)

## -----------------------------------------------------------------------------
constraints_vector <- c(0,0)

## -----------------------------------------------------------------------------
constraints_vector <- c(1,-1)

## -----------------------------------------------------------------------------
constraints_vector <- c(1,1,-1)

## -----------------------------------------------------------------------------
fname_GMT2 <- system.file("extdata", "hsapiens_REAC_subset2.gmt", 
		package = "ActivePathways")

## -----------------------------------------------------------------------------
enriched_pathways <- ActivePathways(
		pval_matrix, gmt = fname_GMT2, cytoscape_file_tag = "Original_")

constraints_vector <- c(1,1)
constraints_vector
dir_matrix[example_genes,]

enriched_pathways_directional <- ActivePathways(
		pval_matrix, gmt = fname_GMT2, cytoscape_file_tag = "Directional_", merge_method = "DPM",
		scores_direction = dir_matrix, constraints_vector = constraints_vector)

## -----------------------------------------------------------------------------
pathways_lost_in_directional_integration = 
		setdiff(enriched_pathways$term_id, enriched_pathways_directional$term_id)

pathways_lost_in_directional_integration
enriched_pathways[enriched_pathways$term_id %in% pathways_lost_in_directional_integration,] 

## -----------------------------------------------------------------------------
wnt_pathway_id <- "REAC:R-HSA-3858494"

enriched_pathway_genes <- unlist(
		enriched_pathways[enriched_pathways$term_id == wnt_pathway_id,]$overlap)
enriched_pathway_genes

## -----------------------------------------------------------------------------
pathway_gene_pvals = pval_matrix[enriched_pathway_genes,]
pathway_gene_directions = dir_matrix[enriched_pathway_genes,]

directional_conflict_genes = names(which(
		pathway_gene_directions[,1] != pathway_gene_directions[,2] &
		pathway_gene_directions[,1] != 0 & pathway_gene_directions[,2] != 0))

pathway_gene_pvals[directional_conflict_genes,]
pathway_gene_directions[directional_conflict_genes,]
length(directional_conflict_genes)

## -----------------------------------------------------------------------------
fname_GMT2 <- system.file("extdata", "hsapiens_REAC_subset2.gmt", 
        package = "ActivePathways")

# Perform integrative pathway enrichment analysis with no directionality
enriched_pathways <- ActivePathways(
        pval_matrix, gmt = fname_GMT2, cytoscape_file_tag = "Original_")

# Perform directional integration and pathway enrichment analysis
constraints_vector <- c(1,1)

enriched_pathways_directional <- ActivePathways(
        pval_matrix, gmt = fname_GMT2, cytoscape_file_tag = "Directional_", merge_method = "DPM",
        scores_direction = dir_matrix, constraints_vector = constraints_vector)

## -----------------------------------------------------------------------------
# Merge the results from both analyses
merged_results <- merge_results(
  enriched_pathways, enriched_pathways_directional,
  gmt_file = fname_GMT2,
  output_prefix = "Aggregated",
  col_colors = c("#FF0000", "#00FF00", "#FFFFF0")
)

## -----------------------------------------------------------------------------
pathways_lost_in_directional_integration = 
        setdiff(enriched_pathways$term_id, enriched_pathways_directional$term_id)
pathways_lost_in_directional_integration

enriched_pathways[enriched_pathways$term_id %in% pathways_lost_in_directional_integration,]

## -----------------------------------------------------------------------------
nrow(ActivePathways(scores, gmt_file, significant = 0.05))
nrow(ActivePathways(scores, gmt_file, significant = 0.1))

## -----------------------------------------------------------------------------
gmt <- read.GMT(gmt_file)
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
ActivePathways(scores, gmt_file, geneset_filter = c(10, 500))

## ----eval=FALSE---------------------------------------------------------------
#  write.GMT(gmt, 'hsapiens_REAC_subset_filtered.gmt')

## -----------------------------------------------------------------------------
background <- makeBackground(gmt)
background <- background[background != 'TP53']
ActivePathways(scores, gmt_file, background = background)

## -----------------------------------------------------------------------------
sort(merge_p_values(scores, 'Fisher'))[1:5]
sort(merge_p_values(scores, 'Brown'))[1:5]

sort(merge_p_values(scores, 'Stouffer'))[1:5]
sort(merge_p_values(scores, 'Strube'))[1:5]

## -----------------------------------------------------------------------------
scores2 <- cbind(scores[, 'CDS'], merge_p_values(scores[, c('X3UTR', 'X5UTR', 'promCore')], 'Brown'))
colnames(scores2) <- c('CDS', 'non_coding')
scores[c(2179, 1760),]
scores2[c(2179, 1760),]

ActivePathways(scores, gmt_file)
ActivePathways(scores2, gmt_file)

## -----------------------------------------------------------------------------

nrow(ActivePathways(scores, gmt_file))
nrow(ActivePathways(scores, gmt_file, cutoff = 0.01))


## -----------------------------------------------------------------------------

nrow(ActivePathways(scores, gmt_file))
nrow(ActivePathways(scores, gmt_file, correction_method = 'none'))


## -----------------------------------------------------------------------------
res <- ActivePathways(scores, gmt_file)
res

## -----------------------------------------------------------------------------
res$overlap[1:3]

## -----------------------------------------------------------------------------
unlist(res[res$term_id == "REAC:422475","evidence"])

## ----eval = FALSE-------------------------------------------------------------
#  result_file <- paste('ActivePathways_results.csv', sep = '/')
#  export_as_CSV (res, result_file) # remove comment to run
#  read.csv(result_file, stringsAsFactors = F)[1:3,]

## ----eval=FALSE---------------------------------------------------------------
#  result_file <- paste('ActivePathways_results2.txt', sep = '/')
#  data.table::fwrite(res, result_file, sep = '\t', sep2 = c('', ',', ''))
#  cat(paste(readLines(result_file)[1:2], collapse = '\n'))

## ----eval=FALSE---------------------------------------------------------------
#  res <- ActivePathways(scores, gmt_file, cytoscape_file_tag = "enrichmentMap__")

## -----------------------------------------------------------------------------
files <- c(system.file('extdata', 'enrichmentMap__pathways.txt', package='ActivePathways'),
           system.file('extdata', 'enrichmentMap__subgroups.txt', package='ActivePathways'),
           system.file('extdata', 'enrichmentMap__pathways.gmt', package='ActivePathways'),
           system.file('extdata', 'enrichmentMap__legend.pdf', package='ActivePathways'))

## ----eval=FALSE---------------------------------------------------------------
#  gmt_file <- system.file('extdata', 'hsapiens_REAC_subset.gmt', package = 'ActivePathways')
#  scores_file <- system.file('extdata', 'Adenocarcinoma_scores_subset.tsv', package = 'ActivePathways')
#  
#  scores <- read.table(scores_file, header = TRUE, sep = '\t', row.names = 'Gene')
#  scores <- as.matrix(scores)
#  scores[is.na(scores)] <- 1
#  
#  res <- ActivePathways(scores, gmt_file, cytoscape_file_tag = "enrichmentMap__")

## -----------------------------------------------------------------------------
cat(paste(readLines(files[1])[1:5], collapse='\n'))
cat(paste(readLines(files[2])[1:5], collapse='\n'))
cat(paste(readLines(files[3])[18:19], collapse='\n'))

## -----------------------------------------------------------------------------
res <- ActivePathways(scores, gmt_file, cytoscape_file_tag = "enrichmentMap__", color_palette = "Pastel1")

## -----------------------------------------------------------------------------
res <- ActivePathways(scores, gmt_file, cytoscape_file_tag = "enrichmentMap__", custom_colors = c("violet","green","orange","red"))

