#Loading all the required packages and dealing with conflicts

library(rentrez)
library(seqinr)
library(stringr)
library(phangorn)
library(Biostrings)
library(muscle)
library(DECIPHER)
library(BiocManager)
library(ape)
library(RSQLite)
library(httr)
library(jsonlite)
library(tidyverse)
install.packages(ggtree)
conflicted::conflict_prefer("filter", "dplyr")
conflicted::conflicts_prefer(dplyr::rename)
library(viridis)
library(rglobi)
library(httr)
library(jsonlite)

####PART 1 - Download data from NCBI cleaning data up for further use----

#investigate the databases on NCBI 

entrez_dbs()


#investigate the data by searching on the NCBI website

Chromadorea18S_search <- entrez_search(db = "nuccore", term = "(Chromadorea[ORGN] AND 18S rRNA AND 200:1500[SLEN])", use_history = T)

Chromadorea18S_search
Chromadorea18S_search$web_history 

Nematoda_taxonomy_Search <- entrez_search(db = "taxonomy", term = "Nematoda[LING]", use_history = T)

Nematoda_taxonomy_Search
Nematoda_taxonomy_Search$web_history

Searchterms <- entrez_db_searchable(db = "taxonomy")
entrez_db_summary(db = "taxonomy")


#Load source function to download larger data sets from NCBI

source("Sourcefiles/Entrez_Functions.R")

#Fetching the sequences in FASTA format and writing to file, multiple files are made each with 100 sequences

FetchFastaFiles(searchTerm = "Chromadorea[ORGN] AND 18S rRNA AND 200:1500[SLEN]", seqsPerFile = 100, fastaFileName = "Chromadorea_18S")

#Combing the many 100 sequence nematode files into one data frame on R

Chromadorea_18S <- MergeFastaFiles(filePattern = "Chromadorea_18S*")

#Cleaning up sequence data to only have the sequence and the species name and removing uncomplete species 

df_Chromadorea_18S <- Chromadorea_18S %>%
  separate(Title, into = c("Code", "Genus", "Species", "rest"), 
           sep = " ", 
           extra = "merge", 
           fill = "right", 
           convert = FALSE)%>%
  mutate(Species = paste(Genus, Species))%>%
  mutate(spaces = str_count(Species, "[\\s\\.\\d]")) %>%
  select(Code, Species, Sequence, spaces)%>%
  filter(spaces == 1, !is.na(spaces), !is.na(Species))

length(unique(df_Chromadorea_18S$Species))

#looking at the sequence length and removing sequences longer than 1000 base pairs 

df_Chromadorea_18S <- df_Chromadorea_18S %>%
  mutate(seqlength = nchar(df_Chromadorea_18S$Sequence))

hist(x = df_Chromadorea_18S$seqlength, xlab = "Sequence Legnth", ylab = "Frequency (No. Sequences")

df_Chromadorea_18S <- df_Chromadorea_18S %>%
  filter(seqlength <= 1000) %>%
  filter(seqlength >=500)

hist(x = df_Chromadorea_18S$seqlength, xlab = "Sequence Legnth", ylab = "Frequency (No. Sequences")

length(unique(df_Chromadorea_18S$Species))





#### PART 3 -Downloading trait data on parasitism from GloBI using pagination----



#defining the limit and skip for pagination 

limit <- 100  
skip <- 0  

# Define the source taxon and interaction type
source_taxon <- "Chromadorea"
interaction_type <- "parasiteOf"

# Initialize empty data frame to store all results
Chromadorea_trait <- data.frame()

# Loop to fetch all data through the API
repeat {
  # Build the URL with query parameters
  url <- paste0(
    "https://api.globalbioticinteractions.org/interaction?",
    "sourceTaxon=", source_taxon,
    "&interactionType=", interaction_type,
    "&limit=", limit,
    "&skip=", skip
  )

  # Send a GET request to the API
  response <- GET(url)
  
  # Parse the response as JSON
  interactions_chunk <- fromJSON(content(response, "text"), flatten = TRUE)
  
  # If no results are retrieved, break the loop
  if (length(interactions_chunk$data) == 0) {
    break
  }
  
  # Append to the complete results
  Chromadorea_trait <- rbind(Chromadorea_trait, as.data.frame(interactions_chunk$data))
  
  # Update the skip value for pagination
  skip <- skip + limit
}
  
#cleaning up the data to only have complete species 

names(Chromadorea_trait)

length(unique(Chromadorea_trait$V2))

df.Chromadorea_trait <- Chromadorea_trait %>%
  mutate(spaces_source = str_count(V2, "[\\s\\.\\d]"),spaces_target = str_count(V12, "[\\s\\.\\d]")) %>%
  filter(spaces_source == 1, !is.na(spaces_source), spaces_target == 1, !is.na(spaces_target)) %>%
  select(Species = V2, Interaction = V10, Target_species = V12, Target_taxonomy = V13) %>%
  distinct() 

length(unique(df.Chromadorea_trait$Species))



#### PART 4 - Organizing vertebrate host species by class ---- 

# Define a vector of possible classes of vertebrates
vertebrate_classes <- c("Mammalia", "Aves", "Reptilia", "Amphibia", "Actinopterygii", "Agnatha", "Chondrichthyes", "Dipnomorpha")

# Function to extract the class of vertebrates
extract_class <- function(taxonomy) {
  # Split the taxonomy string by the delimiter
  taxonomy_levels <- unlist(str_split(taxonomy, " \\| "))
  
  # Find the first match for vertebrate classes
  class_found <- base::intersect(taxonomy_levels, vertebrate_classes)
  
  # Return the first found class or NA if none is found
  return(ifelse(length(class_found) > 0, class_found[1], NA))
}

# Create a new column with the class of vertebrates
df.Chromadorea_traitVert <- df.Chromadorea_trait %>%
  mutate(Class = sapply(Target_taxonomy, extract_class)) %>%
  filter(!is.na(Class)) %>%
  distinct(Species, .keep_all = TRUE)


length(unique(df.Chromadorea_traitVert$Species))


#combing trait data with the sequence data to reduce the number of species in the centroid calculation and therefore reduce the computation time 
?semi_join
df.Chromadorea_all <- full_join(df.Chromadorea_traitVert, df_Chromadorea_18S, join_by("Species" == "Species"))

df.Chromadorea_all <- df.Chromadorea_all %>%
  filter(!is.na(Sequence), !is.na(Target_species))

length(unique(df.Chromadorea_all$Species))
unique(df.Chromadorea_all$Species)


  
####PART 6 - finding a centroid sequence for each species---- 

#I think that I need to run that muscle function at the bottom before I run this -because we need an aligned data set before we do this -- but maybe its also a filtering problem 
#I also need to do some filtering before I do this!!!

#Making the sequences into a DNAstringset

class(df.Chromadorea_all$Sequence)
df.Chromadorea_all$Sequence2 <- DNAStringSet(df.Chromadorea_all$Sequence)
class(df.Chromadorea_all$Sequence2)

#Looking at the sequences on an online browser

BrowseSeqs(df.Chromadorea_all$Sequence2)

#Grouping sequences by species 

grouped_sequences <- split(df.Chromadorea_all$Sequence2, df.Chromadorea_all$Species)

#Creating the function to find the centroid 

calculate_centroid <- function(seqs) {
  # Perform multiple sequence alignment
  alignment <- DNAStringSet(muscle::muscle(seqs, gapOpening = -3000))
  
  # Convert the alignment to a DNAbin object to calculate distances
  alignment_dnabin <- as.DNAbin(alignment)
  
  # Calculate pairwise distance matrix using TN93 model
  dist_matrix <- dist.dna(as.DNAbin(alignment), model = "TN93")
  
  # Calculate centroid (sequence with the lowest sum of pairwise distances)
  centroid_index <- which.min(rowSums(as.matrix(dist_matrix)))
  
  # Extract the sequence of the centroid
  centroid_sequence <- as.character(alignment[centroid_index]) 
  
  return(centroid_sequence)
}


#Apply the centroid calculation for each specimens and save it as a data frame 


centroids <- lapply(grouped_sequences, calculate_centroid)
class(centroids) 
class(df.centroids)

#### PART 5 - Phylogeny tree construction ----


#Centroid conversion to a phyDat format

df.centroids <- centroids %>%
   as.data.frame() %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column() %>%
  mutate(Sequence_split = strsplit(as.character(V1), "")) %>%
  unnest_wider(Sequence_split, names_sep= "V1_") %>%
  mutate(rowname = gsub("\\.", " ", rowname)) %>%
  column_to_rownames(var = "rowname") %>%
  select(-V1) %>%
  as.matrix()

class(df.centroids)
DNA.Chromadorea_centroid <- as.phyDat(df.centroids)
class(DNA.Chromadorea_centroid)

#Parisomonious tree construction 

dist.Chromadorea <- dist.hamming(DNA.Chromadorea_centroid)

tree.Chromadorea <- NJ(dist.Chromadorea)

parsimony(tree.Chromadorea, DNA.Chromadorea_centroid)

treeRatchet <- pratchet(DNA.Chromadorea_centroid, start=tree.Chromadorea, maxit=100,
                        minit=5, k=5, trace=0)

tree.best <- acctran(treeRatchet, DNA.Chromadorea_centroid)

plotTree(tree.best, fsize = 0.6)

fit.tree <- pml(tree.Chromadorea, data = DNA.Chromadorea_centroid)

fit.model <-optim.pml(fit.tree, model = "GTR")
ml_branch_lengths <- fit.tree$tree.Chromadorea$edge.length

#compute the branch lengths 
plot(tree.best)
nodelabels()
?nodelabels

df.calibration <- makeChronosCalib(tree.best)
?pml
chronos(treeRatchet, method = "correlated ", power = 1)


#Trait correlation 



#converting the centroid sequences to a data frame

df.centroids.trait <- centroids %>%
  as.data.frame() %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column(var = "Species") %>%
  dplyr::rename("Centroid" = V1) %>%
  mutate(Species = gsub("\\.", " ", Species))

class(df.centroids.trait)
class(df.centroids.trait$Species)

# combing the centroid data with the trait data to make sure trait data corresponds to the species in the analysis  
?semi_join
df.Chromadorea_centroid <- full_join(df.Chromadorea_traitVert, df.centroids.trait, join_by("Species" == "Species"))

Chromadorea_centroid <- df.Chromadorea_centroid %>%
  filter(!is.na(Centroid), !is.na(Target_species)) %>%
  select(Class, Species) 
  
vector_Chromadorea <- setNames(Chromadorea_centroid$Class,Chromadorea_centroid$Species)
  
class(df.Chromadorea_centroid$Class)

#Trait correlation 

tratitree <- phylosig(tree.best, vector_Chromadorea, method="lambda", test=FALSE, nsim=1000, se=NULL, start=NULL)

?phylosig

#Converting the centroids to a DNAstringset

df.Chromadorea_centroid$Centroid <- as.vector(df.Chromadorea_centroid$Centroid)
class(df.Chromadorea_centroid$Centroid)


set.seed(3)
data(Laurasiatherian)
dm <- dist.hamming(Laurasiatherian)
tree <- NJ(dm)
parsimony(tree, Laurasiatherian)
treeRA <- random.addition(Laurasiatherian)
treeSPR <- optim.parsimony(tree, Laurasiatherian)
# lower number of iterations for the example (to run less than 5 seconds),
# keep default values (maxit, minit, k) or increase them for real life
# analyses.
treeRatchet <- pratchet(Laurasiatherian, start=tree, maxit=100,
                        minit=5, k=5, trace=0)
# assign edge length (number of substitutions)
treeRatchet <- acctran(treeRatchet, Laurasiatherian)
# remove edges of length 0
treeRatchet <- di2multi(treeRatchet)
plot(midpoint(treeRatchet))
add.scale.bar(0,0, length=100)
parsimony(c(tree,treeSPR, treeRatchet), Laurasiatherian


