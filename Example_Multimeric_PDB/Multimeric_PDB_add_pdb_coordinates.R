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
structure_segments.list <- list.files("Example_Multimeric_PDB/Example_output/Sequence_to_structure/")

#obtain available chains 
ow.df <- tibble(chains = str_split(structure_segments.list,"_",simplify = T) %>% .[,4] ,name = structure_segments.list,genes = str_split(structure_segments.list,"_",simplify = T) %>% .[,2])


complete_complex <- tibble()
for(chain_sel in ow.df$chains %>% unique(.)){
  
  ow_sel.df <- ow.df %>% filter(chains == chain_sel)
  
  for(i in 1:nrow(ow_sel.df)){
    
    gene = ow_sel.df$genes %>% unique()
    
    if(i == 1){
      
      combined.df <- read_delim(paste0("Example_Multimeric_PDB/Example_output/Sequence_to_structure/",ow_sel.df$name[i]), delim = "\t") %>%
        mutate(chain = chain_sel) %>%
        left_join(read.pdb("Example_Multimeric_PDB/Example_input/6d6t.pdb1") %>% 
                    .$atom %>% 
                    as_tibble() %>% 
                    filter(elety == "CA",
                           chain == chain_sel) %>% 
                    mutate(ref_aa = convert_aa(resid)) %>% 
                    select(resno,ref_aa,,chain,x,y,z), by = c("Position_in_structure" ="resno", "fastaseq" = "ref_aa","chain" = "chain")) %>% 
        rename(Uniprot_position = "index",
               aminoacid = "fastaseq") 
      
    }else{
      
      addition.df <- read_delim(paste0("Example_Multimeric_PDB/Example_output/Sequence_to_structure/",ow_sel.df$name[i]), delim = "\t") %>%
        mutate(chain = chain_sel) %>% 
        left_join(read.pdb("Example_Multimeric_PDB/Example_input/6d6t.pdb1") %>% 
                    .$atom %>% 
                    as_tibble() %>% 
                    filter(elety == "CA",
                           chain == chain_sel) %>% 
                    mutate(ref_aa = convert_aa(resid)) %>% 
                    select(resno,ref_aa,,chain,x,y,z), by = c("Position_in_structure" ="resno", "fastaseq" = "ref_aa","chain" = "chain")) %>% 
        rename(Uniprot_position = "index",
               aminoacid = "fastaseq") 
      
      combined.df[which(!is.na(addition.df$x)),c(2:9)] <- addition.df[which(!is.na(addition.df$x)),c(2:9)]
      
    }
    
    
    
  }
  
  combined.df$gene <- gene
  
  complete_complex <- rbind(complete_complex,combined.df)
  
}

write_delim(complete_complex,"Example_Multimeric_PDB/Example_output/6d6t_protein_structure_score_input.txt", delim = "\t")

  


