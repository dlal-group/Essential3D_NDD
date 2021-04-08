###20.03.2020
###Add x,y,z coordinates to the processed xml file 
library(tidyverse)
library(bio3d)

convert_aa <- function(Amino){
  paste0(substring(Amino,1,1),tolower(substring(Amino,2,3))) %>% 
    a() %>% 
    return(.)
  
}

#select chain of interest 
chain_sel = "A"

read_delim("Example_Monomeric_PDB/Example_output/SLC2A1_protein_structure_chain_A.txt", delim = "\t") %>%
  mutate(chain = chain_sel) %>% 
  left_join(read.pdb("Example_Monomeric_PDB/Example_input/5eqg.pdb1") %>% 
              .$atom %>% 
              as_tibble() %>% 
              filter(elety == "CA") %>% 
              mutate(ref_aa = convert_aa(resid)) %>% 
              select(resno,ref_aa,,chain,x,y,z), by = c("Position_in_structure" ="resno", "fastaseq" = "ref_aa","chain" = "chain")) %>% 
  rename(Uniprot_position = "index",
         aminoacid = "fastaseq") %>% 
  write_delim("Example_Monomeric_PDB/Example_output/SLC2A1_protein_structure_chain_A_score_input.txt", delim = "\t")
  