# Paisley Samuel
# BMME8053 Intro. to Bioinformatics
# Fall 2023
# Weeks 9 & 10 Assignment
# ____________________________________________ #




#### 16S Functional analyses

## I couldn't find any specific packages that included everything in only one package but here I include 2 packages: microbiomeMarker and ggpicrust2
## MicrobiomeMarker includes a lot of different datasets including 16S data and it includes functions to import picrust2 data into R as a phyloseq object and statistical test functions
## ggpicrust2 is geared more towards pathway annotation and visualization of the picrust2 output




# Packages used

# installing microbiomeMarker
# if (!requireNamespace("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")
# BiocManager::install("microbiomeMarker")

library(microbiomeMarker) #uses primarily phyloseq objects to navigate data
library(phyloseq)
library(ggpicrust2) #install.packages("ggpicrust2")
library(tidyverse)
library(ggh4x)

## microbiomeMarker

# The microbiomeMarker package function import_picrust2 uses the PICRUst2 output 
# and creates a phyloseq object using certain files of the output:
# predicted function abundance profile = otu_table
# functional traits = tax_table (only 1 column if you don't have descriptions, 2 columns if you do have descriptions)
# metadata of the samples = sample_data

# Example from the package
sam_tab <- system.file(
  "extdata", "picrust2_metadata.tsv",
  package = "microbiomeMarker") #retrieving the file path of the metadata
feature_tab <- system.file(
  "extdata", "path_abun_unstrat_descrip.tsv.gz", #retrieving the file path of the function abundance profile
  package = "microbiomeMarker") 
# creating the phyloseq object
ps <- import_picrust2(feature_tab, sam_tab, trait = "PATHWAY") # for trait= , you should choose function trait that you used in PICRUst2 (options = "PATHWAY", "COG", "EC", "KO", "PFAM", "TIGRFAM", "PHENO"), this dataset used "PATHWAY"

# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 195 taxa and 24 samples ]           -> predicted function abundance profile
# sample_data() Sample Data:       [ 24 samples by 3 sample variables ]  -> metadata
# tax_table()   Taxonomy Table:    [ 195 taxa by 2 taxonomic ranks ]     -> functional traits (this dataset included the description so there are 2 "taxonomic ranks" Picrust_trait & Picrust_description)

#Melting phyloseq object into a dataframe
ps_melted <- psmelt(ps)



#___________________________________________________________________________#


## ggpicrust2

#Using ggpicrust2 to make a heatmap of example pathways
#loading in data
data("metacyc_abundance")
data("metadata")
#creating a differential abundance df of the pathway
metacyc_daa_results_df <- pathway_daa(abundance = metacyc_abundance %>%
                                        column_to_rownames("pathway"),
                                      metadata = metadata, group = "Environment", daa_method = "LinDA")
#annotating the da df
annotated_metacyc_daa_results_df <- pathway_annotation(pathway = "MetaCyc",
                                                       daa_results_df = metacyc_daa_results_df, ko_to_kegg = FALSE)
#filtering out significant results
feature_with_p_0.05 <- metacyc_daa_results_df %>% filter(p_adjust < 0.05)
#creating the heatmap of the annotated df and only showing significant results
pathway_heatmap(abundance = metacyc_abundance %>%
                  right_join(annotated_metacyc_daa_results_df %>%
                               select(all_of(c("feature","description"))), by = c("pathway" = "feature")) %>%
                  filter(pathway %in% feature_with_p_0.05$feature) %>%
                  select(-"pathway") %>%
                  column_to_rownames("description"), metadata = metadata, group = "Environment")

#creating a PCA plot of the pathway abundance df
pathway_pca(metacyc_abundance %>% column_to_rownames("pathway"), metadata, "Environment")

