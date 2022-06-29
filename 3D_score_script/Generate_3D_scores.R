####11.04.2022

###Hotspot score via fishers exact test --> evaluation with odds ratios 
library(tidyverse)
library(Rfast)


#2. load mappable protein residues 

load_mapable <- function(file_type,file_path){
  suma_file.df <- read.table(paste0(file_path,file_type), sep = "\t", dec = ".", header = T, quote = "", stringsAsFactors = F)
  suma_file_filter.df <- suma_file.df[which(!is.na(suma_file.df$x)),] ## filtered to obtain only mapable positions 
  return(suma_file_filter.df)
}


##3. count variants either aminoacid position wise --> "pos" or by the number of different substitutions "sub"

count_sub_variants <- function(suma_variant.df){
  variant_sub <- sapply(suma_variant.df$gmiss, function(x){
    if (!is.na(x)){
      strsplit(x,";") %>%
        unlist() %>%
        unique() %>% 
        length()
    } else{
      0
    }
  }) %>%
    as.vector()
  
  return(variant_sub)
}


###4. access variants in neighbourhood (defined by sphere radius) 
number_of_neighbour <- function(variants.patient,variants.gnomad,patient.df,patient_general,radius){
  neighbours <- c()  ###remember rule of pythagoras in 3d = a?+b?+c? = distance of hypothemuse
  neighbour_patient <- c()
  neighbour_gnomad <- c()
  neighbour_p_outside <- c()
  neighbour_g_outside <- c()
  paraz_mean <- c()
  mtr_mean <- c()
  
  neighbour.list <- list()
  
  dist.matrix <- Dist(patient.df %>% select(x,y,z))
  
  for(self.pos in 1:nrow(patient.df)){ ### iterates each aminoacid position to calculate all distances 
    distance <- dist.matrix[self.pos,]
    
    paraz_mean <- c(paraz_mean,mean(patient_general$PARAZ[which(distance<radius)], na.rm = T))
    mtr_mean <- c(mtr_mean,mean(patient_general$MTR[which(distance<radius)], na.rm = T))
    
    neighbour_patient <- c(neighbour_patient,sum(variants.patient[which(distance<radius)]))
    neighbour_p_outside <- c(neighbour_p_outside,sum(variants.patient[which(distance>radius)]))
    
    neighbour_gnomad<- c(neighbour_gnomad,sum(variants.gnomad[which(distance<radius)]))
    neighbour_g_outside<- c(neighbour_g_outside,sum(variants.gnomad[which(distance>radius)]))
    
    neighbours <- c(neighbours,length(which(distance<radius)))

  }
  
  ##zscore brechnung 
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
  neighbour.list$position <- patient.df$Uniprot_position
  
  return(neighbour.list)
}


fisher_score <- function(neighbour.list,variants.patient,variants.gnomad){
  total_number_patient <- sum(variants.patient)
  total_number_gnomad <- sum(variants.gnomad)
  fisher_result <- c()
  fisher_pvalue <- c()
  fisher_lci <- c()
  fisher_uci <- c()
  fish_tot_numg <- c()
  fish_tot_nump <- c()

  for(pos in 1:length(neighbour.list[[1]])){
    fisher_input.matrix <- matrix(c(neighbour.list$n_patient[pos],neighbour.list$n_gnomad[pos],
                                    total_number_patient- neighbour.list$n_patient[pos],
                                    total_number_gnomad- neighbour.list$n_gnomad[pos]), byrow = T, ncol = 2) ##creates fisher matrix which looks like ordered by row: 
    #Number of patients in the sphere; Number of Gnomad variants in the sphere
    #Number of patients outside the sphere; Number of Gnomad variants outside the sphere

    fishtest <- fisher.test(fisher_input.matrix)  ##perfrom fisher test 
    fisher_result <- c(fisher_result,as.numeric(fishtest$estimate))  ##save odds ratios in a vector 
    fisher_pvalue <- c(fisher_pvalue,as.numeric(fishtest$p.value))
    fisher_lci <- c(fisher_lci,as.numeric(fishtest$conf.int[1]))
    fisher_uci <- c(fisher_uci,as.numeric(fishtest$conf.int[2]))
    fish_g_out <- c(fish_tot_numg,total_number_gnomad- neighbour.list$n_gnomad[pos])
    fish_p_out <- c(fish_tot_nump,total_number_patient- neighbour.list$n_patient[pos])
    
    
    
    
  }
  ### converts all Inf's to the maximum odd ratio which was observed and adds +1 --> Only patient variants without any control variant in the sphere 
  #fisher_result[which(fisher_result == Inf)] <-  max(fisher_result[which(fisher_result != Inf)])+1
  
  
  print(paste("The lowest pvalue obtained is:", min(fisher_pvalue)))
  
  fishlist <- list()
  fishlist$inps <- neighbour.list$n_patient 
  fishlist$ings <- neighbour.list$n_gnomad
  fishlist$outps <- neighbour.list$outps
  fishlist$outgs <- neighbour.list$outgs
  fishlist$odd <- fisher_result
  fishlist$pvalue <- fisher_pvalue
  fishlist$lci <- fisher_lci
  fishlist$uci <- fisher_uci
  fishlist$tot_g <- fish_g_out
  fishlist$tot_p <- fish_p_out
  fishlist$total_g <- total_number_gnomad
  fishlist$total_p <- total_number_patient

  return(fishlist)
}

count_subs <- function(df){ ##once again needed in 6. 
  sapply(df$gmiss, function(x){
    ifelse(is.na(x),0, #count subs 
           str_split(x,";") %>% 
             unlist() %>% 
             unique() %>% 
             length())
  })
}

#filter inf to max+1 
remove_inf <- function(vec){
  if(max(vec) == Inf){
    new_vec <- vec[vec != Inf]
    maxim <- max(new_vec)+1
  }
  out <- sapply(vec, function(x){
    ifelse(x == Inf,maxim,x)
  })
  return(out)
}


###6. combine everything to calculate the score 

density_test_combined <- function(file_patient, file_gnomad, file_path, radius){### calculate the scores from scratch, file_type and file_path define protein files which should be loaded. In this example the experimental input files i have sent you were taken (protein: DDX3X)
  
  suma_patient.df <- load_mapable(file_patient,file_path)
  suma_gnomad.df <- load_mapable(file_gnomad,file_path)
  
  read_delim(paste0(file_path,file_patient), delim = "\t", guess_max = 10000) %>% 
    mutate(pscount = count_subs(.)) %>% 
    mutate(pscount = ifelse(is.na(pscount),"none",pscount)) %>% 
    filter(!is.na(x)) %>% 
    dplyr::rename(pconsq = "gmiss") %>% 
    mutate(pconsq = ifelse(is.na(pconsq),"none",pconsq))-> patient_general
  
  read_delim(paste0(file_path,file_gnomad), delim = "\t", guess_max = 10000) %>% 
    mutate(gscount = count_subs(.)) %>% 
    filter(!is.na(x)) %>% 
    dplyr::rename(gconsq = "gmiss") %>% 
    mutate(gconsq = ifelse(is.na(gconsq),"none",gconsq))-> gnomad_general
  
  
  variants.patient <- count_sub_variants(suma_patient.df)
  variants.gnomad <- count_sub_variants(suma_gnomad.df)
  
  
  
  variant_neighbour.list <- number_of_neighbour(variants.patient,variants.gnomad,suma_patient.df,patient_general,radius)
  
  fisher_out <- fisher_score(variant_neighbour.list,variants.patient,variants.gnomad)

  
  result.df <- cbind(patient_general[,c(1:3,5:9)],gnomad_general$gconsq,gnomad_general$gscount,patient_general$pconsq,patient_general$pscount,variant_neighbour.list$n_neighbour,fisher_out$ings, fisher_out$outgs,fisher_out$inps, fisher_out$outps,fisher_out$lci,fisher_out$uci,remove_inf(fisher_out$odd),fisher_out$pvalue,variant_neighbour.list$mtr_mean,variant_neighbour.list$paraz_mean,variant_neighbour.list$mtr_3d,variant_neighbour.list$paraz_3d)


  colnames(result.df) <- c("Uniprot_position","aminoacid","Position_in_structure","x","y","z","paraz","mtr","gconsq","gscount","pconsq","pscount","aa12","gsIn","gsOut","psIn","psOut","lCI","uCI","OR","pvalue","mtr3dmean","paraz3dmean","mtr3dscore","paraz3dscore")

  
  return(result.df)
  
}

#experimentally single 

radii = c(6,8,10,12,14,16,18)


for(radius in radii){
  file_path_expi <- "Molecular_solved_structures/Score_input_files_one_per_gene/"
  file_output <- paste0("Molecular_solved_structures/Output_scores_alternatives/A",radius,"/")
  file_list <- list.files(file_path_expi)

  for( i in 1:(length(file_list)/2)){
    
    name = paste0(paste(str_split(file_list[i*2],"_",simplify = T)[1:2], collapse = "_"), "_","_output_hotcold3d_mtr3d_paraz3d_full_struc.txt")
    
    write_delim(density_test_combined(file_list[i*2],file_list[i*2-1],file_path_expi,radius), paste0(file_output,name), delim = "\t")
    
    
  }
}

