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
library(httr)
library(jsonlite)
#Loading needed packages to clean up the data 

library(tidyverse)
conflicted::conflict_prefer("filter", "dplyr")
conflicted::conflicts_prefer(dplyr::rename)
library(viridis)

####PART 1 - Download data from NCBI cleaning data up for further use----

#look to see if the data I need is present on the NCBI website

Nematoda18S_search <- entrez_search(db = "nuccore", term = "(Nematoda[ORGN] AND 18S rRNA AND 200:1500[SLEN])", use_history = T)

Nematoda18S_search
Nematoda18S_search$web_history ## This stores the ID on the NCBI server.

#Load source function
source("Sourcefiles/Entrez_Functions.R")


# Fetching the sequences in FASTA format and writing to file.
FetchFastaFiles(searchTerm = "Nematoda[ORGN] AND 18S rRNA AND 200:1500[SLEN]", seqsPerFile = 100, fastaFileName = "Nematoda_18S")

# Merging the sequences into a dataframe.
Nematoda_18S <- MergeFastaFiles(filePattern = "Nematoda_18S*")

#Cleaning up sequence data to only have the sequence and the species name and removing uncomplete species 

df_Nematoda_18S <- Nematoda_18S %>%
  separate(Title, into = c("Code", "Genus", "Species", "rest"), 
           sep = " ", 
           extra = "merge", 
           fill = "right", 
           convert = FALSE)%>%
  mutate(Species = paste(Genus, Species))%>%
  mutate(spaces = str_count(Species, "[\\s\\.\\d]")) %>%
  select(Code, Species, Sequence, spaces)%>%
  filter(spaces == 1, !is.na(spaces), !is.na(Species))

length(unique(df_Nematoda_18S$Species))

#looking at the sequence length 

df_Nematoda_18S.explore <- df_Nematoda_18S %>%
  mutate(seqlength = nchar(df_Nematoda_18S$Sequence))

hist(x = df_Nematoda_18S.explore$seqlength, xlab = "Sequence Legnth", ylab = "Frequency (No. Sequences")


#Looking at the sequences on an online browser

BrowseSeqs(df_Onchocercidae_18S$Sequence2)


####PART 2 - finding a centroid sequence for each species---- 

#I think that I need to run that muscle function at the bottom before I run this -because we need an aligned data set before we do this -- but maybe its also a filtering problem 
#I also need to do some filtering before I do this!!!

#Making the sequences into a DNAstringset

class(df_Nematoda_18S$Sequence)
df_Nematoda_18S$Sequence2 <- DNAStringSet(df_Nematoda_18S$Sequence)
class(df_Nematoda_18S$Sequence2)

names(df_Nematoda_18S$Sequence2) <- df_Nematoda_18S$code
names(df_Nematoda_18S$Sequence2)


#Grouping sequences by species 

grouped_sequences <- split(df_Nematoda_18S$Sequence2, df_Nematoda_18S$Species)

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


#Apply the centroid calculation for each species

centroids <- lapply(grouped_sequences, calculate_centroid)



#### PART 3 -Downloading trait data on parasitism from GloBI using pagination----

install.packages("rglobi")
library(rglobi)
library(httr)
library(jsonlite)

#defining the limit and skip for pagination 

limit <- 100  
skip <- 0  

# Define the source taxon and interaction type
source_taxon <- "Nematoda"
interaction_type <- "parasiteOf"

# Initialize empty data frame to store all results
Nematoda_trait <- data.frame()

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
  Nematoda_trait <- rbind(Nematoda_trait, as.data.frame(interactions_chunk$data))
  
  # Update the skip value for pagination
  skip <- skip + limit
}
  
#cleaning up the data to only have complete species 

names(Nematoda_trait)

length(unique(Nematoda_trait$V2))

df.Nematoda_trait <- Nematoda_trait %>%
  mutate(spaces_source = str_count(V2, "[\\s\\.\\d]"),spaces_target = str_count(V12, "[\\s\\.\\d]")) %>%
  filter(spaces_source == 1, !is.na(spaces_source), spaces_target == 1, !is.na(spaces_target)) %>%
  select(Species = V2, Interaction = V10, Target_species = V12, Target_taxonomy = V13) %>%
  distinct() 

length(unique(df.Nematoda_trait$Species))


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
df.Nematoda_traitVert <- df.Nematoda_trait %>%
  mutate(Class = sapply(Target_taxonomy, extract_class)) %>%
  filter(!is.na(Class)) %>%
  distinct(Species, .keep_all = TRUE)
  



#### PART 4 - combing the the centroid sequences for each species with the trait data for each species 
?semi_join
df.Nematoda_all <- full_join(df.Nematoda_trait, df_Nematoda_18S, join_by("Species" == "Species"))

df.onc2 <- semi_join(df_Onchocercidae_18S, df.parasite.oncsub, join_by("Species" == "Species"))

df_Onchocercidae_traitgene <- full_join(df_Onchocercidae_18S, df.parasite.oncsub, join_by("Species" == "source_taxon_name"), keep = TRUE, relationship = "many-to-many")
?full_join

length(unique(df_Onchocercidae_traitgene$Species))
unique(df_Onchocercidae_traitgene$Species)
length(unique(df_Onchocercidae_traitgene$source_taxon_name))
unique(df_Onchocercidae_traitgene$source_taxon_name)

df.Nematoda <- df.Nematoda_all %>%
  filter(!is.na(Sequence), !is.na(Target_species))

length(unique(df.Nematoda_traitVert$Species))
length(unique(df.Nematoda$Sequence))


  






#Starting the alignment - I might want to remove this

df.parasite.alignment <- DNAStringSet(muscle::muscle(df.parasite$Sequence2, maxiters = 2), use.names = TRUE)

df.parasite.alignment

BrowseSeqs(df.parasite.alignment)


length(all.equal(df.parasite.oncsub$source_taxon_name, df_Onchocercidae_18S$Species))
