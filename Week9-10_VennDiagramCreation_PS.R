# Paisley Samuel
# November 2023
# ____________________________________________ #


### Creating Venn Diagrams ###

## Code I used in my thesis but modified it with some new features

#Setting Working directory and Seed
setwd("insert path to folder containing the data files")
set.seed(1998) # choose a number to ensure you get the same results each time you run the code

#Packages used throughout entire script
library(phyloseq)
library(microbiome)
library(vegan)
library(microbiomeutilities)
library(eulerr)
library(gplots)
library(tidyverse)



## Creating abundance table by loading in feature and metadata tables
dat<-read.csv("feature_Y123_ADJUSTED.csv", header=TRUE, row.names = 1)
dat <- t(dat)  # rows should be samples so we have to transpose the dataframe (df); transposing converts dataframe to matrix
row.names(dat) # row names should now be the sample names after transposing
metadata <- read.csv("Metadata-Diversity_BATCH.csv", header = TRUE, row.names = 1)
dat <- as.data.frame(dat)  # need to convert back into a df to complete the rest of the steps

#this chunk of code matches the samples in the feature table to the samples in
#the metadata (since there were some samples that were filtered out after data cleaning)
common.rownames <- intersect(rownames(dat), rownames(metadata)) #finds the common row names (sample names in this case) between the two dfs
dat <- dat[common.rownames,]              #subsets the data to only include the samples from the common.rownames list
metadata <- metadata[common.rownames,]
all.equal(rownames(dat),rownames(metadata)) #run this code to ensure that both dfs have the same samples (should read TRUE if all samples match)

#finding the singletons and doubletons (ASVs that appear only 1 or 2 times)
otu.abund<-which(colSums(dat)>2)

library(vegan)
dat.dom<-dat[,otu.abund] #subsets the data to only include the ASVs that were NOT singletons or doubletons
dat.pa<-decostand(dat.dom, method ="pa")   #converting dominant taxa to presence/absence abundance
dat.otus.01per<-which(colSums(dat.pa) > (0.01*nrow(dat.pa)))  #finding the taxa that appear at least 1% of the time
dat.01per<-dat.dom[,dat.otus.01per]  #subsetting the dominant taxa abundance table to include only the taxa that appear at least 1% of the time


## Making phyloseq objects (using dat.01per abundance)
library(phyloseq)
library(microbiome)
library(microbiomeutilities)
asvdat <- as.data.frame(t(dat.01per)) #taxa has to be rows now so the df was transposed and converted back to a df
taxdat <- read.csv("taxonomy_Y123_edited&cleaned.csv", header = TRUE, row.names = 1) #loading in taxonomy table
meta <- read.csv("Metadata-Diversity_BATCH.csv", header = TRUE, row.names = 1) #loading in metadata (or you can just use the metadata variable that is loaded)
asvmat <- data.matrix(asvdat)
taxmat <- as.matrix(taxdat) # use as.matrix NOT as.data.matrix as the latter will convert the data into numbers
ASV <- otu_table(asvmat, taxa_are_rows = TRUE)
TAX <- tax_table(taxmat)
META <- sample_data(meta)
pseq <- phyloseq(ASV,TAX,META) #creating phyloseq object that contains the metadata, taxonomy, and abundances



# simple way to count number of samples in each group
table(meta(pseq)$Year, useNA = "ifany") #useNA will show you the # of NAs you have in the group if you have any ("ifany")

##   1       2        3    
##  157     210     174      

#convert to relative abundance using transform function in microbiome package 
pseq_rel <- microbiome::transform(pseq, "compositional")

#Make a list of Years (or variable of choice)
years <- unique(as.character(meta(pseq_rel)$Year)) #since the years in the given metadata file are loaded in as numeric, you have to convert them to characters first in order to make a list of the unique levels of the variable
print(years)
# [1] "1" "2" "3"

#format names to show taxa and checking that it worked
pseq_rel_f <- format_to_besthit(pseq_rel)
taxa_names(pseq_rel_f)[1:5]

#Write a for loop to go through each of the years one by one and combine identified core taxa into a list.
list_core <- c() # an empty object to store information

for (n in years){ # for each variable n in Year (or variable of choice)
  ps_sub <- subset_samples(pseq_rel_f, Year == n) # Chooses samples from Year = n
  core_m <- core_members(ps_sub, # ps.sub is phyloseq selected with only samples from g 
                         detection = 0.001, # 0.001 in at least 75% samples (change these values to suit your question)
                         prevalence = 0.75)
  print(paste0("No. of core taxa in ", n, " : ", length(core_m))) # prints the # of  core taxa identified in each year.
  list_core[[n]] <- core_m # add to a list core taxa for each group.
}

print(list_core) #shows ASV id with taxa information (format = ASV ID:taxa info.)

#changing the names of the years (they only have the numbers, ex. year 1 only says 1) for the example data
names(list_core)
# "1" "2" "3"
names(list_core) <- c("Year 1", "Year 2", "Year 3")
names(list_core)
# "Year 1" "Year 2" "Year 3"

##Plotting venn diagram using eulerr
library(eulerr)
plot(venn(list_core),fills = c("tomato3", "steelblue3", "springgreen3"))

## using gplots to list the taxa in each of the venn diagram interactions
library(gplots)
gplots_ven <- venn(list_core, show.plot = F)
venn_intersections <- attr(gplots_ven, "intersections")  #venn_interactions contains a list of the taxa within each interaction on the diagram

#OPTIONAL
#accessing an element in the venn_interaction list and separating ASV id and taxa info.
##this can be streamlined better, especially if you have more variable combinations
library(tidyverse)
Y13 <- as.data.frame(venn_intersections$`Year 1:Year 3`)
Y13 <- separate(Y13, "venn_intersections$`Year 1:Year 3`",c("FeatureID","Taxonomy") , sep=":") 

#make sure to rename the ASV column the same as your original taxonomy file so 
#you can see the entire taxonomic information for that ASV in the original taxonomy file

#OPTIONAL
#Convert the list of the venn diagram interactions into a dataframe that can be exported
vennint<- as.data.frame(unlist(venn_intersections, recursive = TRUE, use.names = TRUE)) #unlisting and turning the resulting vector into a dataframe
#Renaming column and exporting as a csv 
colnames(vennint) <- "ASV-ID:Taxonomy"
write.csv(vennint, "VennDiagramInteractions_ExampleCoreYear.csv")
