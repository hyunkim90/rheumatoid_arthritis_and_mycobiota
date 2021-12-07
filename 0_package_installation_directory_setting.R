## import libraries : phyloseq and microbiome
library(dplyr)
library(forcats) 
library(metagenomeSeq)
library(vegan)
library(phyloseq)
library(microbiome)
library(ggplot2)
library(scales)
library(grid)
library(reshape2)
library(seqtime)
library(agricolae)
library(RColorBrewer)
library(xlsx)
library(magrittr)
library(indicspecies)
library(Hmisc)
library(igraph)
library(qgraph)
library(randomForest)
library(multifunc)
library(FSA)
library(rcompanion)
library(seqinr)


# Set plotting theme
theme_set(theme_bw())

### setting input and output path
### We can then load the biom file with phyloseq function import_biom. We extract the OTU table with OTU abundances and the taxonomy table from the resulting phyloseq object.
bac_phylo=import_biom("OTU_table_final_16s.biom")

### merge with metadata
# Import sample metadata

## in metadata erase # (This step is essential)
map <- read.table(file = 'sample_metadata_16s.tsv', sep = '\t', header = TRUE)
class(map)

map <- sample_data(map)
class(map)

head(map)
dim(map)
summary(map)
str(map)

summary(map)
colnames(map)
rownames(map)
nrow(map)

# Assign rownames to be Sample ID's
map$SampleID

rownames(map) <- map$SampleID
rownames(map)
dim(map)
# Merge biom data object with sample metadata + tree data(this is rooted!)
phy_tree = read_tree("rooted_tree_16s.nwk")
phy <- merge_phyloseq(bac_phylo, map, phy_tree)

class(phy)

phy   ## 1346 otus

## changing rank names
colnames(tax_table(phy)) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

phy



## Fungal community

### setting input and output path
### We can then load the biom file with phyloseq function import_biom. We extract the OTU table with OTU abundances and the taxonomy table from the resulting phyloseq object.
fun_phylo=import_biom("otu_table_final_its2_old.biom")

sample_names(fun_phylo) <- paste("F",sample_names(fun_phylo), sep = "")

### merge with metadata
# Import sample metadata
## maybe Gyeryueng data didn't work because it wasn't transformed to json file

## in metadata erase # (This step is essential)
f.map <- read.table(file = 'sample_metadata_its2.tsv', sep = '\t', header = TRUE)
f.map <- sample_data(f.map)

head(f.map)
dim(f.map)
summary(f.map)
str(f.map)

summary(f.map)
colnames(f.map)
rownames(f.map)
nrow(f.map)

# Assign rownames to be Sample ID's
f.map$SampleID
rownames(f.map) <- f.map$SampleID
rownames(f.map)
dim(f.map)
# Merge biom data object with sample metadata + tree data(this is rooted!)
fun_tree = read_tree("rooted_tree_its2_old.nwk")
fun <- merge_phyloseq(fun_phylo, f.map, fun_tree)

class(fun)
fun   ## 1629 otus (lastest DB 2020)/ 1641 OTUs (old DB 2017)

## changing rank names
colnames(tax_table(fun)) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

fun



### fungal ITS1
fun_phylo=import_biom("otu_table_final_its1.biom")

sample_names(fun_phylo) <- paste("F",sample_names(fun_phylo), sep = "")

### merge with metadata
# Import sample metadata
## maybe Gyeryueng data didn't work because it wasn't transformed to json file

## in metadata erase # (This step is essential)
f.map <- read.table(file = 'sample_metadata_its1.tsv', sep = '\t', header = TRUE)
f.map <- sample_data(f.map)

head(f.map)
dim(f.map)
summary(f.map)
str(f.map)

summary(f.map)
colnames(f.map)
rownames(f.map)
nrow(f.map)

# Assign rownames to be Sample ID's
f.map$SampleID
rownames(f.map) <- f.map$SampleID
rownames(f.map)
dim(f.map)
# Merge biom data object with sample metadata + tree data(this is rooted!)
fun_tree = read_tree("rooted_tree_its1.nwk")
fun <- merge_phyloseq(fun_phylo, f.map, fun_tree)

class(fun)
fun   ## 1412 OTUs (old DB 2017)

## changing rank names
colnames(tax_table(fun)) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

fun
