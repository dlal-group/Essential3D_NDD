#generate protein file for swiss models that were obatined from the swiss model repository 
library(tidyverse)
library(seqinr)
library(bio3d)


convert_aa <- function(Amino){
  paste0(substring(Amino,1,1),tolower(substring(Amino,2,3))) %>% 
    a() %>% 
    return(.)
  
}

fasta <- read.fasta("Example_Monomeric_Swiss_model/Example_input/SLC2A1.fasta")

paraz <- read_delim("Example_Monomeric_Swiss_model/Example_input/SLC2A1_Paraz.txt", delim = "\t",skip = 11)
mtr <- read_delim("Example_Monomeric_Swiss_model/Example_input/SLC2A1_MTR.txt", delim = "\t")

tibble(ref_aa = fasta$ali %>% as.vector(), index = 1:length(fasta$ali %>% as.vector())) %>% 
  left_join(read.pdb("Example_Monomeric_Swiss_model/Example_input/SLC2A1_model.pdb") %>% 
                .$atom %>% 
                filter(elety == "CA") %>% 
                mutate(ref_aa = convert_aa(resid)) %>% 
                select(ref_aa,resno,chain,x,y,z), by = c("ref_aa" = "ref_aa","index" = "resno")) %>% 
  rename(Uniprot_position = "index",
         aminoacid = "ref_aa") %>% 
  mutate(PARAZ = paraz$`(SCORE-MEAN)/STD`,
         MTR = mtr$mtr) %>% 
  write_delim("Example_Monomeric_Swiss_model/Example_output/SLC2A1_protein_structure_chain_A_score_input.txt", delim = "\t")
  

