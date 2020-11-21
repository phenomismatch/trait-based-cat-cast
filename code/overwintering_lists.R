
# Assign overwintering stage to cat count species list ("Caterpillars Count! genus_species list - Sheet1.csv") 
# from discover life table ("DiscLifeMothInventory_Traits.xlsx - Species inventory.csv")
library(tidyverse)
count <- read_csv("CatCountSpeciesList.csv")
count <- read_csv("/Users/marianaabarcazama/Desktop/Projects/trait-based-cat-cast/data/CatCountSpeciesList.csv")
discover <- read_csv("/Users/marianaabarcazama/Desktop/Projects/trait-based-cat-cast/data/DiscLifeTraits.csv")

count <- count %>% 
  mutate_if(is.character, factor) %>% 
  mutate(sp = scientific_name)
discover <- discover %>% 
  select(sp ="Species group", "Overwintering") %>% 
  mutate(winter = ifelse(is.na(Overwintering), "unknown",Overwintering)) %>% 
  mutate_if(is.character, factor) 

# complex <- tibble(sp = c("Desmia funeralis", "Desmia maculalis",
# "Herpetogramma pertextalis",
# "Herpetogramma thestealis",
# Orgyia definita--leucostigma                                                                                                                                                           
# Halysidota harrisii--tessellaris
# Apantesis nais--phalerata--vittata                                                                                                                                                     
# [44] Zanclognatha marcidilinea--jacchusalis      
# Renia adspergillus--fraternalis  
# Crambidia pallida--uniformis
# Dichomeris inversella--kimballi   
# Melanolophia candaria -- signataria
# Hypagyrtis esther--unipunctata  
# Eulithis diversilineata--gracilineata     
# Macaria aequiferaria--bicolorata--bisignata
# Macaria aemulataria--promiscuata 
# Pero ancetaria--honestaria     
# Hydriomena divisaria--pluviata--transfigurata                                                                                                                                          
# [113] Leptostales pannaria                                                                                                                                                                   
# [114] Idaea tacturata                                                                                                                                                                        
# [115] Operophtera brumata--bruceata                                                                                                                 
# Orthosia alurina--hibisci--rubescens                                                                                                             
# Acrobasis angusella--ostryella

newcount <- left_join(count, discover, by = "sp")
unique(newcount$winter)
shared_sp <- intersect(count$sp, discover$sp)

diapa <- discover %>% 
  filter(sp %in% shared_sp)
Ana_dr <- Ana %>% 
  filter(set %in% list_complete) %>% 