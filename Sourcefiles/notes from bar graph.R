#filtering the trait and centroid data to select for the parasite species and mammillian order 

extra.info <- df.Trait_DNA %>%
  select(Parasite_species, Mammillian_order, AdultBodyMass_g)


length(extra.info$Parasite_species) #141 species confirming the length is correct 
length(unique(extra.info$Mammillian_order)) #11 mammalian orders are present 
class(extra.info) #data frame

#the next few steps are done to make sure that the trait data and the tree have the same order of parasitic species for the ggtree construction 

#getting the tip labels from the tree
tree_tips <- tree.best$tip.label

#getting the species names from the trait data
trait_species <- extra.info$Parasite_species  

#checking if all species names are present in the tree
all(trait_species %in% tree_tips) 

#reordering the trait data to match the order of the tree tip labels
extra.info <- extra.info[match(tree_tips, extra.info$Parasite_species), ]

node <- 1:Ntip(tree.best)

p2 <- ggtree(tree.best, branch.length = 'none') 
  
  


facet_plot(p2, panel = "Adult body mass (g)", geom = geom_barh, data = extra.info, mapping = aes(x = AdultBodyMass_g, fill = Mass_grouping), stat = 'identity')+
  scale_fill_manual(name = "Mass grouping", labels = c("2 to 135000", "135000 to 270000", "270000 to 406000", "406000 to 541000", "541000 to 677000"), 
                    values = c("#2F4F4F", "#CD3333", "#7AC5CD", "#D2691E", "#EE6A50"))
            
          
levels=c("2 to 135000", "135000 to 270000", "270000 to 406000", "406000 to 541000", "541000 to 677000")))
?geom_tippoint

?facet_plot
