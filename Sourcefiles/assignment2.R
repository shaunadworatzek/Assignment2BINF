#Loading the needed packages to open data from NCBI and GloBI into R 

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

#Loading needed packages to clean up the data 

library(tidyverse)
conflicted::conflict_prefer("filter", "dplyr")
library(viridis)

####PART 1 - Download data from NCBI cleaning data up for further use----

#look to see if the data I need is present on the NCBI website

gene18S_search <- entrez_search(db = "nuccore", term = "(Onchocercidae[ORGN] AND 18S rRNA AND 200:1500[SLEN])", use_history = T)

gene18S_search
gene18S_search$web_history ## This stores the ID on the NCBI server.

#Load source function
source("Sourcefiles/Entrez_Functions.R")


# Fetching the sequences in FASTA format and writing to file.
FetchFastaFiles(searchTerm = "Onchocercidae[ORGN] AND 18S rRNA AND 200:1500[SLEN]", seqsPerFile = 100, fastaFileName = "Onchocercidae_18S")

# Merging the sequences into a dataframe.
Onchocercidae_18S <- MergeFastaFiles(filePattern = "Onchocercidae_18S*")

#Cleaning up sequence data to only have the sequence and the species name and removing uncomplete species 

df_Onchocercidae_18S <- Onchocercidae_18S %>%
  separate(Title, into = c("Code", "Genus", "Species", "rest"), 
           sep = " ", 
           extra = "merge", 
           fill = "right", 
           convert = FALSE)%>%
  mutate(Species = paste(Genus, Species))%>%
  mutate(spaces = str_count(Species, "[\\s\\.\\d]")) %>%
  select(Code, Species, Sequence, spaces)%>%
  filter(spaces == 1, !is.na(spaces), !is.na(Species))

length(unique(df_Onchocercidae_18S$Species))

#looking at the sequence length 

df_Onchocercidae_18S.explore <- df_Onchocercidae_18S %>%
  mutate(seqlength = nchar(df_Onchocercidae_18S$Sequence))

hist(x = df_Onchocercidae_18S.explore$seqlength, xlab = "Sequence Legnth", ylab = "Frequency (No. Sequences")

#Making the sequences into a DNAstringset

class(df_Onchocercidae_18S$Sequence)
df_Onchocercidae_18S$Sequence2 <- DNAStringSet(df_Onchocercidae_18S$Sequence)
class(df_Onchocercidae_18S$Sequence2)

names(df_Onchocercidae_18S$Sequence2) <- df_Onchocercidae_18S$code
names(df_Onchocercidae_18S$Sequence2)

#Looking at the seqiences on an online broswer

BrowseSeqs(df_Onchocercidae_18S$Sequence2)
  
####PART 2 - finding a centroid sequence for each species---- 

#I think that I need to run that muscle function at the bottom before I run this -because we need an aligned data set before we do this -- but maybe its also a filtering problem 
#I also need to do some filtering before I do this!!!

#Grouping sequences by species 

grouped_sequences2 <- split(df_Onchocercidae_18S$Sequence2, df_Onchocercidae_18S$Species)

#Creating the function to find the centroid 

calculate_centroid <- function(seqs) {
  # Perform multiple sequence alignment
  alignment <- muscle::muscle(seqs, gapOpening = -3000)
  
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


#Apply the centroid calculation for each species

centroids2 <- lapply(grouped_sequences2, calculate_centroid)



#### PART 3 -Downloading trait data on parasitism from GloBI----

install.packages("rglobi")
library(rglobi)

??get_interactions_by_taxa
get_interaction_types()
df.parasite_Onc2 <- get_interactions_by_taxa(sourcetaxon='Onchocercidae', interactiontype='parasiteOf', returnobservations = TRUE, opts = list(),
                                     showfield = c("source_taxon_external_id", "source_taxon_name",
                                                   "interaction_type", "target_taxon_external_id",
                                                   "target_taxon_name", "study_external_id"))


#cleaning up the data to only have complete species 

names(df.parasite_Onc)

length(unique(df.parasite_Onc$source_taxon_name))

df.parasite.oncsub <- df.parasite_Onc %>%
  mutate(spaces_source = str_count(source_taxon_name, "[\\s\\.\\d]")) %>%
  mutate(spaces_target = str_count(target_taxon_name, "[\\s\\.\\d]")) %>%
  select(source_taxon_name, interaction_type, target_taxon_name, latitude, longitude, spaces_source, spaces_target) %>%
  filter(spaces_source == 1, !is.na(spaces_source), spaces_target == 1, !is.na(spaces_target), !is.na(source_taxon_name))
  
length(unique(df.parasite.oncsub$source_taxon_name))

#combing the two data sets 
?semi_join
df.onc1 <- anti_join(df.parasite.oncsub, df_Onchocercidae_18S, join_by("source_taxon_name" == "Species"))
df.onc2 <- semi_join(df_Onchocercidae_18S, df.parasite.oncsub, join_by("Species" == "source_taxon_name"))

df_Onchocercidae_traitgene <- full_join(df_Onchocercidae_18S, df.parasite.oncsub, join_by("Species" == "source_taxon_name"), keep = TRUE, relationship = "many-to-many")
?full_join

length(unique(df_Onchocercidae_traitgene$Species))
unique(df_Onchocercidae_traitgene$Species)
length(unique(df_Onchocercidae_traitgene$source_taxon_name))
unique(df_Onchocercidae_traitgene$source_taxon_name)

df.parasite <- df_Onchocercidae_traitgene %>%
  filter(!is.na(Species), !is.na(source_taxon_name))

length(unique(df.parasite$Species))
length(unique(df.parasite$source_taxon_name))
unique(df.parasite$Species)



#Starting the alignment - I might want to remove this

df.parasite.alignment <- DNAStringSet(muscle::muscle(df.parasite$Sequence2, maxiters = 2), use.names = TRUE)

df.parasite.alignment

BrowseSeqs(df.parasite.alignment)


length(all.equal(df.parasite.oncsub$source_taxon_name, df_Onchocercidae_18S$Species))
