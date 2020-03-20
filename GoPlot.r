library(GOplot)
library(tidyverse)
library(splitstackshape)

df <- read.table("degs_david_GO_chart.txt", sep='\t', header=TRUE, fill=TRUE)
df <- cSplit(df, "Term", "~")
df <- cSplit(df, "Category", "_")
df <- df %>% select(Category_2, Term_1, Term_2, Genes, Benjamini)
colnames(df) <- c('Category', 'ID', 'Term', 'Genes', 'adj_pval')
df <- df[df$adj_pval < 0.17,]
df

#df2 <- read.delim("GO_bp(benja).txt", sep='\t', header=TRUE, fill=TRUE) 
#chart_data <- read.table("degs_david_GO_table.txt", sep='\t', header=TRUE, fill=TRUE)

df2 <- read.delim("mrnas.degs.aml (2).csv", sep=',', header=TRUE) 
df2 <- df2 %>% select(X, log2FoldChange)
colnames(df2) <- c('ID', 'logFC')
df2

circ <- circle_dat(df, df2)

# Generate a simple barplot
GOBar(subset(circ, category == 'BP'))

# Facet the barplot according to the categories of the terms 
GOBar(circ, display = 'multiple')

# Facet the barplot, add a title and change the colour scale for the z-score
GOBar(circ, display = 'multiple', title = 'Z-score coloured barplot', zsc.col = c('yellow', 'black', 'cyan'))

# Generate the bubble plot with a label threshold of 3
GOBubble(circ, labels = 3)

# Add a title, change the colour of the circles, facet the plot according to the categories and change the label threshold
GOBubble(circ, title = 'Bubble plot', colour = c('orange', 'darkred', 'gold'), display = 'multiple', labels = 3)  

# Colour the background according to the category
GOBubble(circ, title = 'Bubble plot with background colour', display = 'multiple', bg.col = T, labels = 3)  

# Reduce redundant terms with a gene overlap >= 0.75...
reduced_circ <- reduce_overlap(circ, overlap = 0.75)
# ...and plot it
GOBubble(reduced_circ, labels = 2.8)

#Generate a circular visualization of the results of gene-annotation enrichment analysis
GOCircle(circ)

df_term <- as.data.frame(df$Term)
colnames(df_term) <- c('process')
df_term

# Now it is time to generate the binary matrix
chord <- chord_dat(circ, df2, df_term$process)
head(chord)

# Generate the matrix with a list of selected genes
chord <- chord_dat(data = circ, genes = df2)
# Generate the matrix with selected processes
chord <- chord_dat(data = circ, process = df_term$process)

# Create the plot
chord <- chord_dat(data = circ, genes = df2, process = df_term$process)
GOChord(chord, space = 0.05, gene.order = 'logFC', gene.space = 0.3, gene.size = 2)

# Display only genes which are assigned to at least three processes
GOChord(chord, limit = c(3, 0), gene.order = 'logFC')

# First, we use the chord object without logFC column to create the heatmap
GOHeat(chord[,-8], nlfc = 0)

# Now we create the heatmap with logFC values and user-defined colour scale
GOHeat(chord, nlfc = 1, fill.col = c('red', 'yellow', 'green'))

GOCluster(circ, df_term$process, clust.by = 'logFC', term.width = 2)

GOCluster(circ, df_term$process, clust.by = 'term', lfc.col = c('darkgoldenrod1', 'black', 'cyan1'))

l1 <- subset(circ, term == 'synapse assembly', c(genes,logFC))
l2 <- subset(circ, term == 'axon guidance', c(genes,logFC))
l3 <- subset(circ, term == 'integral component of plasma membrane', c(genes,logFC))
GOVenn(l1,l2,l3, label = c('synapse assembly', 'axon guidance', 'integral component of plasma membrane'))