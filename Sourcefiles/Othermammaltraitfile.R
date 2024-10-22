#opening Mammal trait data into R and seleting for body size, genus and species, then combing genus and species into one column

install.packages("remotes")
remotes::install_github("RS-eco/traitdata", build_vignettes = T, force=T, build = FALSE)
library(traitdata)
data(pantheria)
vignette("data_info")
data(trait_glossary)


mammal.trait <- pantheria %>%
  select(AdultBodyMass_g, Genus, Species) %>%
  filter(!is.na(AdultBodyMass_g)) %>%
  unite(Target_species, Genus, Species, sep = " ")

#combing mammal trait data with global interactions data

df.traitdata <- full_join(mammal.trait, df.Chromadorea_trait, join_by("Target_species" == "Target_species"))

df.Chromadorea_trait2 <- df.traitdata %>%
  filter(!is.na(Species), !is.na(AdultBodyMass_g)) %>%
  select(Species, AdultBodyMass_g) %>%
  distinct(Species, .keep_all = TRUE)
