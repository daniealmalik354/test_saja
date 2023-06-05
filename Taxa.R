```{r fig.width=10, fig.height=8, message=FALSE, fig.path='Fig_1b/', fig.align='center'}

# Set working directory
setwd("~/Desktop/Coba_Coba")

library(ampvis2)
library(ggplot2)
library(ape)
library(vegan)
library(lattice)
library(permute)
library(dplyr)
library(scales)
library(grid)
library(reshape2)
library(phyloseq)
library(biomformat)

# Read in OTU table
otu_table=read.csv("otu_table_r.csv",sep=",",row.names=1)
otu_table=as.matrix(otu_table)
otu_table[otu_table>0] <- 1 # Present absent

#   Read in taxonomy
# Seperate by kingdom phylum, class, order, family, genus, species
taxonomy <- read.csv("Taxonomy.csv",header=T, sep=",", row.names=1)
taxonomy=as.matrix(taxonomy)

#   Read in metadata

metadata=read.csv("map.csv",sep=",",row.names=1)

#   Read in tree
#   phy_tree=read_tree("tree.nwk")

#   Import as phyloseq objects
OTU = otu_table(otu_table, taxa_are_rows=TRUE)
TAX = tax_table(taxonomy)
META = sample_data(metadata)
# (tree was already imported as a phyloseq object)

##################
UNDIP = phyloseq(OTU, TAX, META); UNDIP 

# Present-absent ASVs dengan jumlah total ASVs 2228
## Jadi, walaupun sudah kita buat dalam bentuk present-absent, jumlah ASVs akan tetap sama utk keseluruhan sanple

##################
## Jumlah proteobacteria taxa/ASVs di semua sample (4 samples)
## Yang disubset adalah yang keseluruhan (UNDIP_3)

UNDIP_prot <- UNDIP %>% subset_taxa(
  Phylum   == "Proteobacteria"); UNDIP_prot

## Hasil:
## phyloseq-class experiment-level object
## otu_table()   OTU Table:         [ 613 taxa and 4 samples ]
## sample_data() Sample Data:       [ 4 samples by 9 sample variables ]
## tax_table()   Taxonomy Table:    [ 613 taxa by 7 taxonomic ranks ]
## ARTINYA: Total ASVs di phylum Proteobaceria di keempat samples adalah: 613

##################
### Jumlah total taxa/ASVs di semua sample (1 samples; S2406)

# Mulai degan subset sample
UNDIP_S2406 <- subset_samples(UNDIP,InputFileName == "S2406"); UNDIP_S2406

## Hasilnya:
## phyloseq-class experiment-level object
## otu_table()   OTU Table:         [ 2228 taxa and 1 samples ]
## sample_data() Sample Data:       [ 1 samples by 9 sample variables ]
## tax_table()   Taxonomy Table:    [ 2228 taxa by 7 taxonomic ranks ]
## ARTINYA: Total ASVs sample S2406 adalah: 2228; Ini sama dengan ASV total, maka dari itu kita ingin lihat hanya ASVs di sample S2406 saja

# Caranya:
ssum <- sum(sample_sums(UNDIP_S2406)); ssum
## Hasilnya 1000, berarti ada 1000 ASVs total di sample S2406

### Jumlah proteobacteria taxa/ASVs di semua sample (1 samples; S2406)

UNDIP_prot1 <- UNDIP_S2406 %>% subset_taxa(
  Phylum   == "Proteobacteria"); UNDIP_prot1

## Hasil:
## phyloseq-class experiment-level object
## otu_table()   OTU Table:         [ 613 taxa and 1 samples ]
## sample_data() Sample Data:       [ 1 samples by 9 sample variables ]
## tax_table()   Taxonomy Table:    [ 613 taxa by 7 taxonomic ranks ]

ssum <- sum(sample_sums(UNDIP_prot1)); ssum
## Hasil: [1] 330
## Ini menytaakan ada 330 ASVs di Phylum Proteobacteria pada sample S2406
## Sedangkan total ASVs di phylum Proteobaceria adalah: 613
## Total ASVs di sample S2406 adalah 1000
## Sedangkan total ASVs di phylum Proteobaceria adalah: 613
## Total ASVs di sample S2406 adalah 1000
```
