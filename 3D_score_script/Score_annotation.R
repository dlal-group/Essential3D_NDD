
####15.03.2022
## Script to annotate 3D normalized scores & Essential3D sites on monomer and multimeric protein structures

###Load packages
library(tidyverse)
library(Rfast)
library(bio3d)
library(seqinr)

##Load required functions
#Calculation of Paraz-3D and MTR-3D and preparation pvEnriched3D calculation
add_3d_scores <- function(input_file.df,radius){
  neighbours <- c()  #residues in 3D bubble
  neighbour_patient <- c() #number of patient variants inside 3D bubble
  neighbour_gnomad <- c() #number of control variants inside 3D bubble
  neighbour_p_outside <- c() #number of patient variants outside 3D bubble
  neighbour_g_outside <- c() #number of control variants outside 3D bubble
  paraz_mean <- c()
  mtr_mean <- c()
  
  neighbour.list <- list()
  
  dist.matrix <- Dist(input_file.df %>% select(x,y,z)) #pairwise distances of all residues 
  
  for(self.pos in 1:nrow(input_file.df)){ ### iterates each protein position to and to identify neighbouring residues and count variants
    distance <- dist.matrix[self.pos,]
    
    paraz_mean <- c(paraz_mean,mean(input_file.df$Paraz_score[which(distance<radius)], na.rm = T))
    mtr_mean <- c(mtr_mean,mean(input_file.df$MTR_score[which(distance<radius)], na.rm = T))
    
    neighbour_patient <- c(neighbour_patient,sum(input_file.df$N_Pathogenic[which(distance<radius)]))
    neighbour_p_outside <- c(neighbour_p_outside,sum(input_file.df$N_Pathogenic[which(distance>=radius)]))
    
    neighbour_gnomad<- c(neighbour_gnomad,sum(input_file.df$N_Control[which(distance<radius)]))
    neighbour_g_outside<- c(neighbour_g_outside,sum(input_file.df$N_Control[which(distance>=radius)]))
    
    neighbours <- c(neighbours,length(which(distance<radius)))
    
  }
  
  ##zscore calculation 
  paraz_mean <- ifelse(is.na(paraz_mean), Inf,paraz_mean)
  mtr_mean <- ifelse(is.na(mtr_mean), Inf, mtr_mean)
  sdparaz <- sd(paraz_mean[paraz_mean != Inf])
  sdmtr <- sd(mtr_mean[mtr_mean != Inf])
  mean_paraz <- mean(paraz_mean[paraz_mean != Inf])
  mean_mtr <- mean(mtr_mean[mtr_mean != Inf])
  
  
  neighbour.list$paraz_mean <- paraz_mean
  neighbour.list$mtr_mean <- mtr_mean
  neighbour.list$paraz_3d <- (paraz_mean-mean_paraz)/sdparaz
  neighbour.list$mtr_3d <- (mtr_mean-mean_mtr)/sdmtr
  neighbour.list$paraz_3d <- ifelse(is.na(neighbour.list$paraz_3d),Inf,neighbour.list$paraz_3d)
  neighbour.list$mtr_3d <- ifelse(is.na(neighbour.list$mtr_3d),Inf,neighbour.list$mtr_3d)
  
  neighbour.list$outps <- neighbour_p_outside
  neighbour.list$outgs <- neighbour_g_outside
  neighbour.list$n_patient <- neighbour_patient
  neighbour.list$n_gnomad <- neighbour_gnomad
  neighbour.list$n_neighbour <- neighbours
  neighbour.list$position <- input_file.df$Position_in_protein
  
  return(neighbour.list)
}
##pvEnriched3D annotation
pvEnriched3D <- function(neighbour.list,input_file.df){
  total_number_patient <- sum(input_file.df$N_Pathogenic)
  total_number_gnomad <- sum(input_file.df$N_Control)
  ftest_result <- c()
  ftest_pvalue <- c()
  ftest_lci <- c()
  ftest_uci <- c()
  ftest_tot_numg <- c()
  ftest_tot_nump <- c()
  
  for(pos in 1:length(neighbour.list[[1]])){
    ftest_input.matrix <- matrix(c(neighbour.list$n_patient[pos],neighbour.list$n_gnomad[pos],
                                   total_number_patient- neighbour.list$n_patient[pos],
                                   total_number_gnomad- neighbour.list$n_gnomad[pos]), byrow = T, ncol = 2) ##creates fisher's exact test input matrix 
    
    ftest <- fisher.test(ftest_input.matrix)  ##perfrom fisher's exact test 
    ftest_result <- c(ftest_result,as.numeric(ftest$estimate))
    ftest_pvalue <- c(ftest_pvalue,as.numeric(ftest$p.value))
    ftest_lci <- c(ftest_lci,as.numeric(ftest$conf.int[1]))
    ftest_uci <- c(ftest_uci,as.numeric(ftest$conf.int[2]))
    ftest_g_out <- c(ftest_tot_numg,total_number_gnomad- neighbour.list$n_gnomad[pos])
    ftest_p_out <- c(ftest_tot_nump,total_number_patient- neighbour.list$n_patient[pos])
    
    
    
    
  }
  
  print(paste("The lowest pvalue obtained is:", min(ftest_pvalue)))
  
  ftestlist <- list()
  ftestlist$inps <- neighbour.list$n_patient 
  ftestlist$ings <- neighbour.list$n_gnomad
  ftestlist$outps <- neighbour.list$outps
  ftestlist$outgs <- neighbour.list$outgs
  ftestlist$odd <- ftest_result
  ftestlist$pvalue <- ftest_pvalue
  ftestlist$lci <- ftest_lci
  ftestlist$uci <- ftest_uci
  ftestlist$tot_g <- ftest_g_out
  ftestlist$tot_p <- ftest_p_out
  ftestlist$total_g <- total_number_gnomad
  ftestlist$total_p <- total_number_patient
  
  return(ftestlist)
}

#Load your input file with the position of each residue in the structure, the scores which will be converted to 3D and the number of Patient and control variants at each residue position

input_name <- "SCN2A_6j8e_input"
save_name <- str_replace(input_name,"input","annotated") ##file name to save the output

input_file.df <- read_delim(paste0("Input_files/",input_name,".txt"), delim = "\t")

radius = 12 #selected radius in angstrom used do determine in which radius the 3D normalization should take place (size of 3D bubble)

#add structure coordinates of selected structure to input file
pdb_structure.df <- read.pdb("Input_files/6j8e.pdb1") %>% 
  .$atom %>% 
  as_tibble() %>% 
  filter(elety == "CA") %>% #C-alpha atom is set as reference for the position of each protein residue 
  select(resno, chain, resid,x,y,z) %>% 
  setNames(c("Structure_position","Structure_chain","Amino_acids","x","y","z")) %>% 
  mutate(Amino_acids = paste0(str_sub(Amino_acids,1,1),str_sub(Amino_acids,2,3) %>% tolower()) %>% a())

input_file_w_struc.df <- input_file.df %>% 
  left_join(pdb_structure.df,by = c("Amino_acids" = "Amino_acids","Structure_position" = "Structure_position","Structure_chain" = "Structure_chain")) #Only protein positions covered by the selected protein structure are considered



scores_3D.list <-  add_3d_scores(input_file_w_struc.df,radius) #Calculation of Paraz-3D and MTR-3D and prepration pvEnriched3D calculation

pvEnriched3D.list <- pvEnriched3D(scores_3D.list,input_file_w_struc.df) #Calculation of pvEnriched3D's


#combine all annotations and derive Essential3D/Non-Essential3D annotations
cbind(input_file_w_struc.df[,c(1:8)],scores_3D.list$n_neighbour,pvEnriched3D.list$lci,pvEnriched3D.list$uci,pvEnriched3D.list$odd,pvEnriched3D.list$pvalue,scores_3D.list$mtr_3d,scores_3D.list$paraz_3d) %>% 
  setNames(c("Position_in_protein","Amino_acids","Structure_position","Structure_chain","Paraz-score","MTR-score","N_control","N_pathogenic","Residues_in_bubble","pvEnriched3D_lCI","pvEnriched3D_uCI","pvEnriched3D_Odds","pvEnriched3D_Pvalue","mtr3dscore","paraz3dscore")) %>% 
  #annotate pvEnriched3D
  mutate(pvEnriched3D = ifelse(pvEnriched3D_Pvalue <0.05 & pvEnriched3D_Odds >1,"pvEnriched3D",
                            ifelse(pvEnriched3D_Pvalue <0.05 & pvEnriched3D_Odds <1,"None-pvEnriched3D","neutral")),
         Paraz_available = ifelse(any(paraz3dscore !=Inf),"yes","no"),
         mtr_available = ifelse(any(mtr3dscore !=Inf),"yes","no"),
         ##Annotate Essential-3D annotation
         Essential3D = ifelse(paraz3dscore >0 & 
                                    Paraz_available == "yes" & 
                                    mtr_available == "yes"&
                                    paraz3dscore != Inf & 
                                    mtr3dscore != Inf &
                                    mtr3dscore <0 & 
                                    pvEnriched3D == "pvEnriched3D","Essential-3D",
                                  ifelse(Paraz_available == "no"&
                                           mtr_available == "yes"&
                                           mtr3dscore <0 & 
                                           mtr3dscore != Inf &
                                           pvEnriched3D == "pvEnriched3D","Essential-3D",
                                         ifelse(Paraz_available == "yes"&
                                                  mtr_available == "no"&
                                                  paraz3dscore >0 & 
                                                  paraz3dscore != Inf &
                                                  pvEnriched3D == "pvEnriched3D","Essential-3D",
                                                ifelse(Paraz_available == "no"&
                                                         mtr_available == "no"&
                                                         pvEnriched3D == "pvEnriched3D","Essential-3D",
                                                       #None_essential3D
                                                       ifelse(paraz3dscore <0 & 
                                                                Paraz_available == "yes" & 
                                                                mtr_available == "yes"&
                                                                paraz3dscore != Inf & 
                                                                mtr3dscore != Inf &
                                                                mtr3dscore >0 & 
                                                                pvEnriched3D == "None-pvEnriched3D","None-Essential-3D",
                                                              ifelse(Paraz_available == "no"&
                                                                       mtr_available == "yes"&
                                                                       mtr3dscore >0 & 
                                                                       mtr3dscore != Inf &
                                                                       pvEnriched3D == "None-pvEnriched3D","None-Essential-3D",
                                                                     ifelse(Paraz_available == "yes"&
                                                                              mtr_available == "no"&
                                                                              paraz3dscore <0 & 
                                                                              paraz3dscore != Inf &
                                                                              pvEnriched3D == "None-pvEnriched3D","None-Essential-3D",
                                                                            ifelse(Paraz_available == "no"&
                                                                                     mtr_available == "no"&
                                                                                     pvEnriched3D == "None-pvEnriched3D","None-Essential-3D","neutral")))))))),
         pvEnriched3D = ifelse(pvEnriched3D == "pvEnriched3D" ,"yes","no"),
         #invert mtr3dscore to ensure higher z-score account for higher missense constraints in 3D
         mtr3dscore = mtr3dscore*-1,
         ##binary pConserved3D and mIntolerant annotation
         Conservation3D = ifelse(!is.na(paraz3dscore) & paraz3dscore >0,"yes","no"),
         mIntolerant3D = ifelse(!is.na(mtr3dscore) & mtr3dscore >0,"yes","no")) %>% 
  select(1,2,3,4,7,8,5,6,9,15,20,14,21,16,19) %>% 
  setNames(c("Position_in_protein","Amino_acids","Structure_position","Structure_chain","N_control","N_pathogenic","Paraz-score","MTR-score","Residues_in_bubble","pConservation3D (score)","pConservation3D","mIntolerance3D (score)","mIntolerant3D","pvEnriched3D","Essential3D")) %>% 
  write_delim(paste0("Output_file/",save_name,".txt"), delim = "\t")
  
