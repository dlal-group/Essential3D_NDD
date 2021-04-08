require(XML)
library(seqinr)
library(here)
library(tidyverse)
library(bio3d)


convert_aa <- function(Amino){
  paste0(substring(Amino,1,1),tolower(substring(Amino,2,3))) %>% 
    a() %>% 
    return(.)
  
}


##Download file that contains macthes the aminoacid positions from the canonical transcript to the protein positions from ftp://ftp.ebi.ac.uk/pub/databases/msd/sifts/xml/5eqg.xml

#Load xml file 
xmlfile <- xmlParse("Example_Monomeric_PDB/Example_input/SLC2A1.xml.gz") 

xmllist <- xmlToList(xmlfile)

#Select pdb ID and desired uniprot ID from protein of interest. Nessesary if the PDb file contains several proteins 
uniprot <- "P11166"

pdbid <- "5eqg"


#Obtain protein structure annotations for selected protein 
entity <- 3
chains <- c()
while (!is.null(xmllist[entity]$entity)){
  
  chains <- c(chains,xmllist[entity]$entity$.attrs[[2]])
  if (length(unique(chains)) != entity-2){
    print("Several entitys per Chain!")
  }
  entity = entity + 1
}
if(entity == 3){
  print("The position of entity doesnt start at 3")
}
entity <- c(3:(entity-1))



iterate_df <- data.frame(matrix(ncol = 3))
colnames(iterate_df) <- c("entity","segment","uniprot")
for (entitys in entity){
  
  for (segment in 1:(length(xmllist[entitys]$entity)-1)){
    entity_seg_reslist_uni <- c()
    for (reslist in 1:(length(xmllist[entitys]$entity[segment]$segment$listResidue))){
      
      if(!is.null(xmllist[entitys]$entity[segment]$segment$listResidue[reslist]$residue[2]$crossRefDb[[1]])){
        
        if(xmllist[entitys]$entity[segment]$segment$listResidue[reslist]$residue[2]$crossRefDb[[1]] == "UniProt"){
          entity_seg_reslist_uni <- c(entity_seg_reslist_uni,xmllist[entitys]$entity[segment]$segment$listResidue[reslist]$residue[2]$crossRefDb[[3]])
          
        }
      }
      
    }
    
    iterate_df[nrow(iterate_df)+1,] <- NA
    iterate_df$entity[nrow(iterate_df)] <- entitys
    iterate_df$segment[nrow(iterate_df)] <- segment
    
    iterate_df$uniprot[nrow(iterate_df)] <- ifelse(is.null(entity_seg_reslist_uni),"no uniprot",ifelse(length(unique(entity_seg_reslist_uni))==0 | length(unique(entity_seg_reslist_uni))>1,"to much uniprot",unique(entity_seg_reslist_uni)))
  }
}   
iterate_df <- iterate_df[2:nrow(iterate_df),]    

if(any(iterate_df$uniprot == "to much uniprot")){
  print("more than one uniprot in a segment!")
}
iterate_df <- iterate_df[which(iterate_df$uniprot != "no uniprot"),]  


entity_uniprot_matches <- sapply(1:nrow(iterate_df), function(x){ifelse(any(iterate_df$uniprot[x] == uniprot),which(iterate_df$uniprot[x] == uniprot),0)})


### function which iterates over each entity
for (match in 1:(length(entity_uniprot_matches))){
  entity <- iterate_df$entity[match]
  segment <- iterate_df$segment[match]
  
  
  #Create a Table based on the protein fasta sequence. Checks if the PDB file contains the selected gene
  if(entity_uniprot_matches[match] != 0 & length(which(is.na(entity_uniprot_matches))) == 0){
    
    fasta_file <- read.fasta("Example_Monomeric_PDB/Example_input/SLC2A1.fasta") 
    len = fasta_file$ali %>% length()
    
    data <- data.frame(matrix(1:len,ncol = 1))
    colnames(data) <- "index"
    
    data$fastaseq <- as.vector(fasta_file$ali)
    
    
    data[c("Position_in_structure","MTR","Paraz")] <- NA
    
    
    #Add Paraz-score 
    
    data$Paraz <- read_delim("Example_Monomeric_PDB/Example_input/SLC2A1_Paraz.txt",delim = "\t",skip = 11) %>% .$`(SCORE-MEAN)/STD`
      
    #Add MTR 
    data$MTR <- read_delim("Example_Monomeric_PDB/Example_input/SLC2A1_MTR.txt", delim = "\t") %>% .$mtr
    
    #Add corresponding protein positios 
    for (reslist in 1:length(xmllist[entity]$entity[segment]$segment$listResidue)){
      #iterate over each residuelist elemeent 
      if(xmllist[entity]$entity[segment]$segment$listResidue[reslist]$residue[2]$crossRefDb[[1]] == "UniProt" & xmllist[entity]$entity[segment]$segment$listResidue[reslist]$residue[2]$crossRefDb[[2]] == "UniProt" & xmllist[entity]$entity[segment]$segment$listResidue[reslist]$residue[2]$crossRefDb[[3]] == iterate_df$uniprot[match]){ #uniprot line matches as expected 
        
        if(xmllist[entity]$entity[segment]$segment$listResidue[reslist]$residue[1]$crossRefDb[[1]] == "PDB" & xmllist[entity]$entity[segment]$segment$listResidue[reslist]$residue[1]$crossRefDb[[3]] == pdbid){
          
          # check if line is indeed  pdb line 
          if(xmllist[entity]$entity[segment]$segment$listResidue[reslist]$residue[1]$crossRefDb[[4]] == "null" | is.na(convert_aa(xmllist[entity]$entity[segment]$segment$listResidue[reslist]$residue[1]$crossRefDb[[5]]))){
            
            print(paste0("no or unknown aminaocid ",reslist))
            
          }else{
            if(data$fastaseq[which(data$index == as.numeric(xmllist[entity]$entity[segment]$segment$listResidue[reslist]$residue[2]$crossRefDb[[4]]))] == convert_aa(xmllist[entity]$entity[segment]$segment$listResidue[reslist]$residue[1]$crossRefDb[[5]]))            
              
              ###add position to data 
              pos_add <- which(data$index == as.numeric(xmllist[entity]$entity[segment]$segment$listResidue[reslist]$residue[2]$crossRefDb[[4]])) ##poisition where all info is added
              
              data$Position_in_structure[pos_add] <- as.numeric(xmllist[entity]$entity[segment]$segment$listResidue[reslist]$residue[1]$crossRefDb[[4]])
              
              
          }
        } else {
          
          print("PDB Line is not at position one")
        }
        
        
      }
      
    }

    chain <- xmllist[entity]$entity[segment]$segment$listResidue$residue$crossRefDb[[6]]
    
    write.table(data,paste0("Example_Monomeric_PDB/Example_output/SLC2A1_protein_structure_chain_",chain,".txt"), sep = "\t", dec =".",row.names = F)
    
  }
  
}




