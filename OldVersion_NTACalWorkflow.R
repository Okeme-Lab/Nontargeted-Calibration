#-------Load Libraries-------
library(ggplot2)
library(ggsankey)
library(classyfireR)
library(reticulate)
library(RColorBrewer)
library(tidyverse)
library(rcompanion)
library(writexl)
library(rstatix)
library(tibble)
library(lmtest)
library(gridExtra)
library(grid)
library(dplyr)
library(ggpubr)
library(tidyr)
library(stringr)
library(broom)
library(Metrics)
library(devtools)
library(dplyr)
library(gridExtra)
library(cowplot)
library(plyr)
library(naniar)
library(ggh4x)
library(EnvStats)
library(readr)
library(ggvenn)
library(ggVennDiagram)
library(readxl)
library(pheatmap)
library(janitor)
library(ggrepel)
library(reticulate)
library(networkD3)
library(reticulate)
library(insight)
library(tree)
library(ape)
library(ctxR)
library(svglite)
library(ggtext)

#python
#load required packages
py_require("pandas") 
py_require("numpy")
py_require('matplotlib')
py_require("scipy")
py_require("khipu")
py_require('khipu-metabolomics')

#import packages
pd <- import("pandas")
np <- import("numpy")
plt <- import("matplotlib.pyplot")
ttest_ind <- import("scipy.stats")$ttest_ind
khipu <- import("khipu")
khipu_extend <- import("khipu.extended")

#-------set API key for ctxR-------
#set default API key using personal key

#-------Set working directory-------
wd <- "C://Users//scrol//OneDrive - McMaster University//Joseph Okeme's files - Alyss//New passives calibration//Alyss New Calculations//Updated Workflow2"
setwd(wd)

#-------Load functions-------
ms_dial_clean <- function(df){
  df <- as.data.frame(df)
  df <- df[-c(1:3),]  #remove unwanted rows#
  names(df) <- df[1,] #make first row the column names#
  df <- df[-c(1),] #remove duplicate first row#
  names(df) <- gsub(" ", "_", names(df)) #replace spaces in column names with underscore#
  names(df) <- gsub("/", ".", names(df)) #replace slashes in column names with period#
  names(df) <- gsub("-", "_", names(df)) #replace dashes in column names with underscore#
  names(df) <- gsub("%", "percent", names(df)) #replace percentages with the word percent#
  index <- match("1",names(df)) #create variable for column index for the first instance of "1"#
  colnames(df)[index] <- "Sample.Avg" #set the corresponding column for sample average from 1 to sample.avg#
  colnames(df)[index+1] <- "Sample.Std" #set the correponding column for sample std from 1 sto sample.std"#
  df$Sample.Avg = as.numeric(df$Sample.Avg)
  df$Sample.Std = as.numeric(df$Sample.Std)
  df <- df %>% as.data.frame(row.names = 1:nrow(.)) #reset row number#
  df <- subset(df, select = -c(Post_curation_result, Fill_percent, Isotope_tracking_parent_ID,
                               Isotope_tracking_weight_number)) #remove unneeded columns#
  df #so output is a data frame#
} #cleans a msdial 5 output, assuming samples haven't been labelled with different categories
clean_agilentquant <- function(df){
  df_clean <- df[-c(1),] #remove first row
  names(df_clean) <- gsub(" ", "_", names(df_clean)) #replace spaces in column names with underscore#
  names(df_clean) <- gsub("-", "!", names(df_clean)) #replace spaces in column names with underscore#
  
  df_clean <- df_clean %>%
    subset(select = -c(Sample, `...2`, `...3`, `...7`, `...6`, `...5`)) %>%
    dplyr::rename(Data_File = `...4`) %>%
    dplyr::rename(Sample = Data_File)
  
  start_col <- 3
  step <- 2
  
  # Loop through the columns starting from `start_col` and rename every other column
  for(i in seq(start_col, ncol(df_clean), by = step)) {
    if (i > 1) {
      colnames(df_clean)[i] <- paste0(colnames(df_clean)[i - 1], "_istd")
    }
  }
  
  
  df_clean_list <- list()
  n_cols <- ncol(df_clean)-1
  start <-2
  
  #for loop to get every other column and bind together
  
  rename_by_pos = function(df, index, new_name){
    colnames(df)[index] = new_name
    df
  } #function to rename columns based on their index
  
  for(i in seq(start, n_cols, by = 2)){
    if (i + 1 <= ncol(df_clean)) {
      subset_df <- df_clean[,c(1, i:(i + 1))]
    }
    subset_df2<-subset_df %>%
      rename_by_pos(2,"Retention_Time") %>%
      rename_by_pos(3,"Area")
    compound <- colnames(subset_df[,c(2)])
    subset_df2[,'Compound'] = compound
    subset_df2 <- subset_df2 %>%
      separate(col = Compound, into = c("Compound", "Discard"), sep = "_R") %>% #split word results into second column
      subset(select = -c(Discard))
    df_clean_list[[i]] <- subset_df2
  }
  
  df <- bind_rows(df_clean_list) %>%
    mutate(Compound = gsub("!", "-", Compound)) %>%# reverse earlier name change
    mutate(Compound = gsub("_", " ", Compound)) # reverse earlier name change
  
  df
} #cleans agilent quant output, exported in wide format with retention time and area column per compound
annotate <- function(df, formulas, structures, unreal_level2){
  df1 <- df %>% 
    mutate(Simple_dot_product = as.numeric(Simple_dot_product)) %>% #make filtering parameteres numerics - nulls will turn into NAs
    mutate(Matched_peaks_count = as.numeric(Matched_peaks_count)) %>%
    mutate(Matched_peaks_percentage = as.numeric(Matched_peaks_percentage)) %>%
    mutate(Alignment_ID = as.numeric(Alignment_ID)) %>%
    filter(MS.MS_assigned == 'True') %>% #filter for things with MSMS
    filter(!is.na(Simple_dot_product) & !is.na(Matched_peaks_percentage) & !is.na(Matched_peaks_count)) #filter out things without a potential spectral match
  
  level2 <- df1 %>% 
    mutate(Simple_dot_product = as.numeric(Simple_dot_product)) %>%  #ensure these are still numeric
    mutate(Matched_peaks_count = as.numeric(Matched_peaks_count)) %>% 
    mutate(Matched_peaks_percentage = as.numeric(Matched_peaks_percentage)) %>%
    filter(Matched_peaks_count >= 3 & Simple_dot_product >= sqrt(0.84) & Matched_peaks_percentage >= 0.5) %>% #filter for level 2 criteria being met 
    filter(!Alignment_ID %in% unreal_level2)  #filter out anything manually checked and ocnfirmed not to be level 2
  if(nrow(level2) != 0){
    level2 <- level2 %>% mutate(level = 2) #mark as level 2 - if not an if statement, will get error when nothing is level 2
  }
  
  highest_rank_structures <- structures %>% #get highest ranked sirius structures for level 3s
    group_by(mappingFeatureId) %>%
    mutate(structurePerIdRank = as.numeric(structurePerIdRank)) %>% 
    slice_min(order_by = structurePerIdRank, n = 1) %>% #get highest ranked (i.e., rank 1) for each structure with SIRIUS prediction
    ungroup()
  
  test2 <- length(unique(highest_rank_structures$structurePerIdRank))  #make sure every structure has something ranked 1 (otherwise data may have been lost when importing SIRIUS results)
  if(test2 != 1){
    warning('Formatting issue')
  }
  
  df3 <- df %>%
    filter(!Alignment_ID %in% level2$Alignment_ID) %>% #filter for features not level 2
    filter(Alignment_ID %in% highest_rank_structures$mappingFeatureId) #anything with a highest rank predicted sirius structure is now level 3

  structure_info <- highest_rank_structures %>% #these will be used to replace results in MS-DIAL for level 3 to match SIRIUS
    select(c(name, smiles, mappingFeatureId, InChIkey2D, molecularFormula)) %>%
    dplyr::rename(Alignment_ID = mappingFeatureId)
  
  level3_sirius <- merge(df3, structure_info, by = c('Alignment_ID')) %>% #merge level 3s with the sirius structures, transfer name, inchikey and smles
    mutate(Metabolite_name = name) %>% #transfer metabolite name
    mutate(INCHIKEY = InChIkey2D) %>% #transfer inchikey
    mutate(SMILES = smiles) %>% #transfer smiles
    mutate(Formula = molecularFormula) %>% #transfer formula
    mutate(m.z_similarity = NA) %>% #not used anymore
    mutate(Weighted_dot_product = NA) %>% #make all this information NAs because this is level 3
    mutate(Reverse_dot_product = NA) %>% #not used because no spectral match
    mutate(Matched_peaks_count = NA) %>% #not used because no spectral match
    mutate(Matched_peaks_percentage = NA) %>% #not used because no spectral match
    subset(select = -c(smiles, name, InChIkey2D, molecularFormula)) #remove the sirius columns (which are now duplicates)
  if(nrow(level3_sirius) != 0){
    level3_sirius <- level3_sirius %>% mutate(level = 3) 
  } #mark level 3 features as level 3
  
  all <- rbind(level2, level3_sirius) #bind level 2 and 3 into a df
  
  lvl5 <- df %>% #remaining features are level 5
    filter(!Alignment_ID %in% all$Alignment_ID) %>% #filter for non-level 2 or 3 features
    mutate(level = 5) %>% #mark as level 5
    mutate(Formula = NA) %>% #remove because level 5 formula isn't known
    mutate(Metabolite_name = NA) %>% #remove - no identity for level 5
    mutate(INCHIKEY = NA) %>% #remove - no identity for level 5
    mutate(SMILES = NA) %>% #remove - no identity for level 5
    mutate(m.z_similarity = NA) %>% #remove - no spectral match
    mutate(Weighted_dot_product = NA) %>% #remove - no spectral match
    mutate(Reverse_dot_product = NA) %>% #remove - no spectral match
    mutate(Matched_peaks_count = NA) %>% #remove - no spectral match
    mutate(Matched_peaks_percentage = NA)  #remove - no spectral match
  
  test3 <- nrow(df) - nrow(all) - nrow(lvl5) #check that the number of features is same as the start of this function (i.e., nothing lost)
  if(test3 != 0){
    warning('Lost data')
  }
  
  all2 <- rbind(all, lvl5) #combine all results, now annotated
  return(all2)
}
dedup <- function(df){
  annotated <- df %>% #first, filter for level 2 and 3 compounds that can be deduplicated
    mutate(INCHIKEY = ifelse(level > 3, NA, INCHIKEY)) %>% #make sure INCHIKEY doesn't exist for things less confident than level 3
    mutate(INCHIKEY = ifelse(INCHIKEY == 'null', NA, INCHIKEY)) %>%
    mutate(INCHIKEY = ifelse(INCHIKEY == 'Unknown', NA, INCHIKEY)) %>%
    mutate(`2D_INCHIKEY` = INCHIKEY) %>% #make the 2D inchikey inchikey to eztract 2D inchikeys for dedupicating 
    mutate(`2D_INCHIKEY` = gsub("\\-.*", "", `2D_INCHIKEY`)) %>%  #get 2D inchikeys by removing everything after the first dash
    mutate(Metabolite_name = ifelse(Metabolite_name == 'null', NA, Metabolite_name)) %>%
    mutate(Metabolite_name = ifelse(Metabolite_name == 'Unknown', NA, Metabolite_name)) #remove unknown metabolite names
  
  #first, isolate features to deduplicate (level 3+)
  features_to_dedup <- annotated %>% 
    filter(level < 4)
  
  #isolate features not to be deduplicated
  features_to_not_dedup<- annotated %>%
    filter(level >= 4)
  
  #deduplicate features with inchikey and with known name - this is to stope unecessary data loss for the NAs which would all be marked as duplicates of each other 
  with_inchikey_and_name <- features_to_dedup %>%
    filter(!is.na(`2D_INCHIKEY`) & !is.na(Metabolite_name)) %>%
    group_by(`2D_INCHIKEY`) %>% #deduplicate duplicate inchikeys
    mutate(level = as.numeric(level)) %>% #ensure these are numeric for deduplicating
    mutate(Mean = as.numeric(Mean)) %>%
    slice_min(order_by = level, n = 1) %>% #first keep most confident annotation
    slice_max(order_by = Mean, n = 1) %>% #then, keep highest mean abundance
    ungroup()
  
  dedup_across_methods2 <- features_to_dedup %>% #now, with things without inchikeys, deduplicate by metabolite name
    filter(is.na(`2D_INCHIKEY`) & !is.na(Metabolite_name)) %>%
    rbind(., with_inchikey_and_name) %>% #combine so everything with metabolite name is in one data frame
    group_by(Metabolite_name) %>%
    mutate(level = as.numeric(level)) %>%
    mutate(Mean = as.numeric(Mean)) %>%
    slice_min(order_by = level, n = 1) %>% #dedup across methods, keeping highest level then mean 
    slice_max(order_by = Mean, n = 1) %>% #mean because total score not comparable across EI and ESI
    ungroup() 
  
  dedup_across_methods3 <- dedup_across_methods2 %>%
    filter(!is.na(`2D_INCHIKEY`))
  dedup_across_methods4 <- dedup_across_methods2 %>%
    filter(is.na(`2D_INCHIKEY`)) 
  
  dedup_across_methods5 <- features_to_dedup %>% #now, with things without metabolite name, deduplicate by metabolite name
    filter(!is.na(`2D_INCHIKEY`) & is.na(Metabolite_name)) %>%
    rbind(., dedup_across_methods3) %>% #bind wiht things with 2D INCHIKEY previously deduplicated by metabolite name
    group_by(`2D_INCHIKEY`) %>%
    mutate(level = as.numeric(level)) %>%
    mutate(Mean = as.numeric(Mean)) %>%
    slice_min(order_by = level, n = 1) %>% #dedup across methods, keeping highest level then mean 
    slice_max(order_by = Mean, n = 1) %>% #mean because total score not comparable across EI and ESI
    ungroup() 
  
  annotated_dedup <- features_to_dedup %>%
    filter(is.na(Metabolite_name) & is.na(`2D_INCHIKEY`)) %>%
    rbind(., dedup_across_methods5, dedup_across_methods4, features_to_not_dedup)
  return(annotated_dedup)
}
get_empirical_features <- function(df, khipu_df){
  khipu_results_clean <- khipu_df
  
  #1. get a cleaner khipu output
  khipu_results_clean_subset <- khipu_results_clean %>% 
    dplyr::rename(Alignment_ID = id) %>% #rename to be able to merge with feature table
    subset(select = c(Alignment_ID, interim_id, neutral_formula_mass)) #get only columns neede from khipu output
  
  #2. filter for features in the khipu output (degenerate features), then merge with khipu results
  khipu_features1 <- df %>%  
    filter(Alignment_ID %in% khipu_results_clean$id) %>%
    merge(., khipu_results_clean_subset, by = c('Alignment_ID')) 
  
  #3. for degenerate features, filter for those with MSMS and deduplicate within MSMS
  khipu_features_MSMS <- khipu_features1 %>% #prioritize MSMS
    mutate(Total_score = as.numeric(Total_score)) %>%
    mutate(Mean = as.numeric(Mean)) %>%
    mutate(m.z_similarity = as.numeric(m.z_similarity)) %>%
    filter(MS.MS_assigned == 'True') %>% #filter for with MSMS
    group_by(interim_id) %>% #interim ID groups degenerate features of the same empirical one together
    slice_max(order_by = Total_score, n = 1) %>% #prioritize score
    ungroup() %>%
    group_by(interim_id) %>%
    slice_max(order_by = Mean, n = 1) %>% #keep highest abundance
    ungroup()
  
  #4. get the khipu features without MSMS
  khipu_features_noMSMS <- khipu_features1 %>% 
    filter(MS.MS_assigned != 'True') %>%
    mutate(Total_score = as.numeric(Total_score)) %>%
    mutate(m.z_similarity = as.numeric(m.z_similarity)) %>%
    filter(!interim_id %in% khipu_features_MSMS$interim_id) %>%
    group_by(interim_id) %>%
    slice_max(order_by = Mean, n = 1) %>% #keep highest abundance features
    ungroup()
  
  #5. combine deduplicated khipu features
  khipu_features <- rbind(khipu_features_MSMS, khipu_features_noMSMS)
  
  #6. cehck that there are no degenerate features remaining
  test <- khipu_features %>% group_by(interim_id) %>% filter(n()>1)
  if(nrow(test) !=0){
    warning('duplicates in khipu features remain')
  }
  
  #7. get the singletons
  not_khipu_features <- df %>%
    filter(!Alignment_ID %in% khipu_results_clean$id) %>%
    mutate(Total_score = as.numeric(Total_score)) %>% #mutate to same data type to allow binding with empirical featutres
    mutate(m.z_similarity = as.numeric(m.z_similarity))
  
  empirical_features_list <- bind_rows(khipu_features, not_khipu_features) %>% #now can use this to cross reference SIRIUS
    separate(col = Metabolite_name, into = c("discard", "Metabolite_name"), sep = ": ", extra = 'merge', fill = 'left') %>%
    subset(select = -c(discard))
  
  return(empirical_features_list)
}
convert_to_smi_tab_delimited <- function(df, label) {
  smi_content <- ""  # Initialize output
  
  #Loop through each row
  for (i in 1:nrow(df)) {
    smile_text <- df$SMILES[i]
    label_text <- as.character(df[[label]][i])  # Safely extract label column as string
    
    #Escape backslashes in SMILES strings
    smile_text <- gsub("\\\\", "\\\\\\\\", smile_text)
    
    #Append tab-delimited line (SMILES \t Label)
    smi_line <- paste(smile_text, label_text, sep = "\t")
    smi_content <- paste(smi_content, smi_line, sep = "\n")
  }
  
  return(smi_content)
}

#-------Set variables-------
#goal: define variables for automated code
Manual <-  read_csv("SurrRedo3 - BacktoCSV2.csv", col_types = cols(.default = "c"), col_names = F) 
#containing manual integrations for surrogate standard, 13C18-TPP

ISTD <-'MTPhP' #internal standard as it appears in the mass hunter output; note MTPhP refers to 13C18-TPP

Method_Blank_String <- '_MB_' #string that only method blanks have
Solvent_Blank_String <- '_SB_UP_' #string that only solvent blanks have
Sample_Blank_String <- '_F_' #string that both samples and field blanks have
Sample_String <- '_S_' #stirng that only samples have
Blank_String <- '_Bk_' #string that only field blanks have
Spiked_Solvent_String <- '_SB_SP_' #string that only spiked solvent blanks have
all_runs_string = 'JO_' #appears in names of all runs and distinguishes data file names from MS-DIAL columns
Mode = 'pos' #'pos' or 'neg' for khipu (ESI mode)
wd_khipu = "C://Users//scrol//OneDrive - McMaster University//Joseph Okeme's files - Alyss//New passives calibration//Alyss New Calculations//Updated Workflow2" #second working directory for khipu if needed (useful to separate feature tables for different modes as they are all named the same)
label_samples <- function(df){ #write function to decide how to label samples for blank matching later.
  df <- df  %>%
    mutate(Sampler = ifelse(grepl('_A_', Sample), 'Active', NA)) %>%
    mutate(Sampler = ifelse(grepl('_P_', Sample) & grepl('_D_', Sample), 'PDMS-Passive', Sampler)) %>%
    mutate(Sampler = ifelse(grepl('_P_', Sample) & grepl('_G_', Sample), 'GFF-Passive', Sampler)) %>%
    mutate(Location = ifelse(grepl('_F_', Sample), 'Field', NA)) %>%
    mutate(Location = ifelse(grepl('_L_', Sample), 'Lab', Location)) %>%
    mutate(Sample_Type = ifelse(grepl(Blank_String, Sample), 'Blank', NA)) %>%
    mutate(Sample_Type = ifelse(grepl(Sample_String, Sample), 'Sample', Sample_Type)) %>%
    mutate(Blank_Batch = ifelse(grepl('_A_', Sample), #adds the batch for each blank used in matching blanks to appropriate samples
                                'A1',
                                ifelse(grepl('_D_', Sample),
                                       'D1',
                                       'G1'))) %>%
    mutate(Batch = ifelse(grepl('1', Sample),
                          1,
                          ifelse(grepl('2', Sample),
                                 2,
                                 ifelse(grepl('3', Sample),
                                        3,
                                        ifelse(grepl('4', Sample),
                                               4,
                                               5))))) 
  df
} #rewrite this as needed

add_batch <- function(df){ #write function to add back the batch.
  df <- df  %>%
    mutate(Batch = ifelse(grepl('1', Sample),
                          1,
                          ifelse(grepl('2', Sample),
                                 2,
                                 ifelse(grepl('3', Sample),
                                        3,
                                        ifelse(grepl('4', Sample),
                                               4,
                                               5))))) 
  df
}

merge_blanks_by = c('Sampler', 'Location', 'Blank_Batch') #should match what's in the previous function

khipu_extend$adduct_search_patterns

adduct_search_patterns = list(khipu_extend$adduct_search_patterns[[1]], khipu_extend$adduct_search_patterns[[2]], list(17.02655, 'NH3')) #adducts to consider when deduplicating
khipu_mz = 5 #tolerance m/z for khipu
khipu_rt = 2 #retention time tolerance for khipu

Blank_Threshold = 0.35 #ratio value (e.g., 0.1 = sample 10 blank)
detection_freq_threshold = 50 #in percent

unreal_level2 <- c(13639, 12148, 19833, 8961, 10333, 11634, 7055, 10019, 12713, 12469, 17823, 9912,
                   7439, 15306,  7004, 8185, 4702, 10693, 2229, 6730, 8871, 19770, 9623, 14388,
                   9175, 14962, 20346, 20642, 1510, 1709, 2144, 4448, 5060, 5906, 5908,
                   6809, 7440, 9622, 10335, 12468, 14390, 15113, 18676, 19831, 4447, 5059, 8187, 5909) #manually checked spectra that are not level 2 (obtained from an original list of potential level 2 candidates )

#-------Pre-cleaning if necessary for istd/sstd-------
#goal 1: have a manual integrations file where one column is compound, one is sample and one is area
#sample names were odd -> rename

#sample names were odd -> rename
Manual_clean <- Manual %>% #clean data frame
  subset(select = -c(X1, X2, X5, X6, X7)) %>% #remove column not needed
  row_to_names(2, remove_rows_above = FALSE) #make the second row the row names

manual_list <- list() #make empty list

for(i in seq(3,(ncol(Manual_clean)-2), by = 3)){#loop through all columns by groups of three starting from column 3
  Manual_clean[1,i+1] = Manual_clean[1,i] #make first row of the second column equal to name of the compound 
  Manual_clean[1,i+2] = Manual_clean[1,i] #make first row of the third column equal to name of the compound
  df <-Manual_clean[,c(1:2, i:(i+2))] #make a new df with the first two columns (sample names), and the group of three columns
  df[,'Compound'] <- df[1, 3] #make a new column for the compound, equal to the 1st row of the third column
  df <- df %>%
    mutate(Compound = gsub( " Results", "", Compound)) #remove the Results part of the compound name#
  df <- df[-c(1),] #remove first row now that it's in its own column
  manual_list[[i]] <- df #add all dfs to a list
}

Manual_clean <- bind_rows(manual_list) %>% #turn cleaned list into a df
  dplyr::rename(Sample = `Data File`) %>%
  dplyr::rename(Area = Resp.) %>%
  mutate(Compound = gsub(" Method", "", Compound)) %>%
  subset(select = c(Sample, Area, Compound)) %>%
  mutate(Sample = gsub("30819", "", Sample)) %>% #remove 30819 from names to enable filtering later#
  mutate(Sample = gsub("JOQT", "JO", Sample)) %>%#remove QT from names to enable filtering later#
  filter(!grepl("p2e", Sample) | grepl("repeat", Sample)) %>%
  filter(!grepl("p3a", Sample) | grepl("repeat", Sample)) %>% #remove two Data_Files that had been rerun#
  filter(!grepl("HV", Sample)) %>% #remove extra Data_Files#
  filter(!grepl("308", Sample)) %>% #remove extra Data_Files#
  filter(!grepl("q", Sample)) %>% #remove extra Data_Files#
  filter(!grepl("Q", Sample)) %>% #remove extra Data_Files#
  filter(!grepl("S2", Sample)) %>% #remove extra Data_Files#
  filter(!grepl("ACN", Sample)) %>% #remove extra Data_Files#
  filter(!is.na(Area))%>% #remove blanks that had no spike and therefore no peak area
  mutate(Sample = gsub("_repeat", "", Sample)) %>% #drop the repeat from the name#
  mutate(Sample = ifelse(grepl("Q", Sample) & !grepl("p", Sample),
                         gsub("JO_A", "JO_L_MB_A_", Sample),
                         ifelse(grepl("q", Sample) & grepl("p", Sample),
                                gsub("JO_p", "JO_L_MB_P_D_", Sample),
                                ifelse(grepl("S2", Sample),
                                       gsub("JO_S2", "JO_L_SB_SP_S2_", Sample),
                                       ifelse(grepl("ACN2", Sample),
                                              gsub("JO_ACN2", "JO_L_SB_UP_ACN2_", Sample),
                                              ifelse(grepl("bl", Sample),
                                                     gsub("JO_A", "JO_F_Bk_A_", Sample),
                                                     ifelse(grepl("p0", Sample) & !grepl("d", Sample) & !grepl("e", Sample) & !grepl("f", Sample),
                                                            gsub("JO_p0", "JO_F_Bk_P_D_", Sample),
                                                            ifelse(grepl("p0", Sample) & !grepl("b", Sample) & !grepl("c", Sample),
                                                                   gsub("JO_p0", "JO_F_Bk_P_G_", Sample),
                                                                   ifelse(grepl("A", Sample),
                                                                          gsub("JO_A", "JO_F_S_A_", Sample),
                                                                          ifelse(grepl("p", Sample) & !grepl("d", Sample) & !grepl("e", Sample) & !grepl("f", Sample),
                                                                                 gsub("JO_p", "JO_F_S_P_D_", Sample),
                                                                                 gsub("JO_p", "JO_F_S_P_G_", Sample))))))))))) #renaming

#-------Get surrogate standard areas--------
ISTD_Only <- Manual_clean %>% 
  filter(Compound %in% ISTD) %>% # filter for internal standard in manual integrations
  subset(select = c(Compound, Sample, Area)) #only use needed columns to keep df clean

#-------Upload Feature Table-------
feature_table <- read_delim(file = "Post-Level 2/NTA_PAS_Cal_PeakAreas_PostLvl2.txt", col_types = cols(.default = "c"), delim = "\t") #upload alignment table (from MSDial 5, do NOT label samples as blank or sample in MSDial)

feature_table_clean <- ms_dial_clean(feature_table) %>% #clean alginment table 
  mutate(Average_Mz = as.numeric(Average_Mz)) %>% #was uploaded as character to prevent data loss -> turn back to numeric for calculations
  mutate(`Average_Rt(min)` = as.numeric(`Average_Rt(min)`)) %>%
  mutate(Alignment_ID  = as.numeric(Alignment_ID))

#-------additional cleaning of alignment table if needed--------

feature_table_clean <- ms_dial_clean(feature_table) %>%
  mutate(Average_Mz = as.numeric(Average_Mz)) %>%
  mutate(`Average_Rt(min)` = as.numeric(`Average_Rt(min)`)) %>%
  mutate(Alignment_ID  = as.numeric(Alignment_ID)) %>%
  pivot_longer(cols = starts_with("JOQT"), names_to = "sample",
               values_to = "Area") %>%
  mutate(sample = gsub("30819", "", sample)) %>% #remove 30819 from names to enable filtering later#
  mutate(sample = gsub("JOQT", "JO", sample)) %>%#remove QT from names to enable filtering later#
  filter(!grepl("p2e", sample) | grepl("repeat", sample)) %>%
  filter(!grepl("p3a", sample) | grepl("repeat", sample)) %>% #remove two samples that had been rerun#
  filter(!grepl("HV", sample)) %>% #remove extra samples#
  mutate(sample = gsub("_repeat", "", sample)) %>% #drop the repeat from the sample name#
  mutate(sample = ifelse(grepl("Q", sample) & !grepl("p", sample),
                         gsub("JO_A", "JO_L_MB_A_", sample),
                         ifelse(grepl("q", sample) & grepl("p", sample),
                                gsub("JO_p", "JO_L_MB_P_D_", sample),
                                ifelse(grepl("S2", sample),
                                       gsub("JO_S2", "JO_L_SB_SP_S2_", sample),
                                       ifelse(grepl("ACN2", sample),
                                              gsub("JO_ACN2", "JO_L_SB_UP_ACN2_", sample),
                                              ifelse(grepl("bl", sample),
                                                     gsub("JO_A", "JO_F_Bk_A_", sample),
                                                     ifelse(grepl("p0", sample) & !grepl("d", sample) & !grepl("e", sample) & !grepl("f", sample),
                                                            gsub("JO_p0", "JO_F_Bk_P_D_", sample),
                                                            ifelse(grepl("p0", sample) & !grepl("b", sample) & !grepl("c", sample),
                                                                   gsub("JO_p0", "JO_F_Bk_P_G_", sample),
                                                                   ifelse(grepl("A", sample),
                                                                          gsub("JO_A", "JO_F_S_A_", sample),
                                                                          ifelse(grepl("p", sample) & !grepl("d", sample) & !grepl("e", sample) & !grepl("f", sample),
                                                                                 gsub("JO_p", "JO_F_S_P_D_", sample),
                                                                                 gsub("JO_p", "JO_F_S_P_G_", sample))))))))))) %>%
  pivot_wider(names_from = sample, values_from = Area)

#-------upload MS2 SIRIUS results-------
#structures must be exported as excel files to avoid data loss (exporting as .tsv or .csv then importing to excel also doesn't work)
#classes excel doesn't work
classes <- read_tsv("SIRIUS_PAScopy//canopus_formula_summary_all.tsv",  col_types = cols(.default = "c"), comment = "") 
formulas <- read_excel("SIRIUS_PAScopy//formula_identifications_all.xlsx", col_types = "text")
structures <- read_excel("SIRIUS_PAScopy//structure_identifications_all.xlsx", col_types = "text")

#quality check
structures_test <- structures %>%
  group_by(mappingFeatureId) %>%
  dplyr::summarise(min = min(structurePerIdRank)) 

unique(structures_test$min) #should give 1 and only 1 (make sure no formulas were lost - everything should have something ranked 1)

#-------make feature table for khipu-------
all_samples <- colnames(feature_table_clean[,grep(Sample_Blank_String, colnames(feature_table_clean))]) #get list of all samples/blanks
all_runs <-  colnames(feature_table_clean[,grep(all_runs_string, colnames(feature_table_clean))]) #get list of all runs

feature_table_khipu <- feature_table_clean %>%
  select(c(Alignment_ID, Average_Mz, `Average_Rt(min)`, all_of(all_runs))) #khipu accepts columns for id, mass, retention time and intensities 

write_tsv(feature_table_khipu,'feature_table_khipu.tsv') #export khipu table to import later 

#-------Prepare data for khipu-------
print_color('Preparing data for khipu', color = 'green') #status message

intensity_column_x <- colnames(feature_table_khipu[,grep(all_runs_string, colnames(feature_table_khipu))])[1] #first run's index
intensity_column_y <- colnames(feature_table_khipu[,grep(all_runs_string, colnames(feature_table_khipu))])[length(colnames(feature_table_khipu[,grep(all_runs_string, colnames(feature_table_khipu))]))] #last run's index

intensity_column_x_index <- match(intensity_column_x, colnames(feature_table_khipu))-1 #-1 to match python indexing
intensity_column_y_index <- match(intensity_column_y, colnames(feature_table_khipu)) #doesn't work if -1 for some reason. (column doesn't get read in python)
mz_col_index_number <- match("Average_Mz", colnames(feature_table_khipu))-1 #get index for mass column (-1 because python indices start at 0)
rt_col_index_number <- match("Average_Rt(min)", colnames(feature_table_khipu))-1 #get index for retention time column 

peaklist <- khipu_extend$local_read_file(paste(wd_khipu, '/feature_table_khipu.tsv', sep = ''), start_col = as.integer(intensity_column_x_index), end_col = as.integer(intensity_column_y_index)) #import peak list to python
khipu_list <- khipu_extend$peaklist_to_khipu_list(peaklist, isotope_search_patterns=list(), #turn off isotope search 
                                     adduct_search_patterns=adduct_search_patterns, mz_tolerance_ppm=khipu_mz, 
                                     rt_tolerance=khipu_rt, mode= Mode) #find duplicates/adducts

khipu_feature_list = khipu_extend$export_empCpd_khipu_list(khipu_list[[1]]) #get list of empirical compounds (first df of khipu_list)
print_color('Cleaning khipu output', color = 'green')

khipu_feature_list <- pd$DataFrame(khipu_feature_list) #turninto dataframe using pandas
test_list <- list() #get blank list
for (i in 1:nrow(khipu_feature_list)){ #clean khipu output list
  df <- khipu_feature_list[i,] #isolate each row
  df2 <- df %>%
    subset(select = c(interim_id, neutral_formula_mass)) #subset non-list, useful columns
  MS1_pseudo_Spectra <- df$MS1_pseudo_Spectra #isolate column with intensities
  intensity_list <- list() #initiate list
  for(m in 1:length(MS1_pseudo_Spectra)){
    outer_list <- MS1_pseudo_Spectra[[m]] #get each most outer list
    temp_list <- list() #initiate list
    for(z in 1:length(outer_list)){
      sub_list <- outer_list[[z]] #for the sublist, get each sub list 
      cols_bound <- bind_cols(sub_list, .name_repair = 'universal') #bind lists into one df as columns
      temp_list[[z]] <- cols_bound
    }
    temp_list_df <- bind_rows(temp_list) #start binding rows back into a data frame
    intensity_list[[m]] <- temp_list_df
  }
  intensity_df <- bind_rows(intensity_list) %>%
    mutate(interim_id = df2$interim_id) %>%
    mutate(neutral_formula_mass = df2$neutral_formula_mass) #put neutral mass as a column
  test_list[[i]] <- intensity_df
} 
khipu_results_clean <- bind_rows(test_list) 
write.csv(khipu_results_clean, 'khipu_results_clean.csv', row.names = F) #export
khipu_results_clean <- read.csv('khipu_results_clean.csv') #import so you don't have to run khipu again

#-------Normalize to sstd-------
spiked_solvent <- colnames(feature_table_clean[,grep(Spiked_Solvent_String, colnames(feature_table_clean))]) #get list of all method blanks
spike_blank_sample <- c(spiked_solvent, all_samples) #ilst of everything that has surrogate standard

feature_table_long <- feature_table_clean %>%
  pivot_longer(cols = all_of(spike_blank_sample), names_to = 'Sample', values_to = 'Sample_Area') %>% #pivot longer to do calculations
  merge(., ISTD_Only, by = c('Sample')) %>% #merge with surrogate standard so sstd area is in a column
  subset(select = -c(Compound)) %>% #remove compound column (came from std data frame indicating what the sstd is)
  mutate(Sample_Area = as.numeric(Sample_Area)) %>%
  dplyr::rename(ISTD_Area = Area) %>%
  mutate(ISTD_Area = as.numeric(ISTD_Area)) %>%
  mutate(ISTD_Normalized_Area = Sample_Area/ISTD_Area) %>% #normalize
  label_samples() #mark what each sample type is for later merging with blanks

#-------Blank Subtract-------
Blanks <- feature_table_long %>%
  filter(grepl(Blank_String, Sample)) %>% #filter for blanks
  select(c(Sample, Alignment_ID, `Average_Rt(min)`, Average_Mz, Metabolite_name, MS.MS_assigned, INCHIKEY, Matched_peaks_percentage,Formula,
           SMILES, MS.MS_matched, m.z_similarity, Simple_dot_product, Weighted_dot_product, Reverse_dot_product, Matched_peaks_count, Total_score, MS1_isotopic_spectrum, MS.MS_spectrum, 
           all_of(merge_blanks_by), ISTD_Normalized_Area)) %>% #select only necessary columns
  dplyr::rename(Blank_Area = ISTD_Normalized_Area) %>% #rename
  dplyr::rename(Blank_Sample = Sample) %>% #rename
  group_by(Alignment_ID, Blank_Batch) %>% #group by blank type
  dplyr::mutate(Avg_PDMS_Blank = ifelse(grepl("_D_", Blank_Sample),
                                        mean(Blank_Area),
                                        NA)) %>% #average blank for PDMS-Mesh-PAS
  dplyr::mutate(Avg_GFF_Blank = ifelse(grepl("_G_", Blank_Sample),
                                        mean(Blank_Area),
                                        NA)) %>% #average blank for GFF-Mesh-PAS
  dplyr::mutate(Avg_Active_Blank = ifelse(grepl("_A_", Blank_Sample),
                                        mean(Blank_Area),
                                        NA)) %>% #average blank for AAS
  slice(1) %>% #only need average, so just keep one replcaite
  ungroup()

Samples <- feature_table_long %>%
  filter(grepl(Sample_String, Sample)) %>% #filter for samples
  select(c(Sample, Alignment_ID, `Average_Rt(min)`, Average_Mz, Metabolite_name, MS.MS_assigned, INCHIKEY, Matched_peaks_percentage, Formula, 
           SMILES, MS.MS_matched, m.z_similarity, Simple_dot_product, Weighted_dot_product, Reverse_dot_product, Matched_peaks_count, Total_score, MS1_isotopic_spectrum, MS.MS_spectrum,
           all_of(merge_blanks_by), ISTD_Normalized_Area)) %>% #select necessary columns
  dplyr::rename(Sample_Area = ISTD_Normalized_Area) #rename

Samples_Blank_Merge <- merge(Samples, Blanks, by = c('Alignment_ID', "Average_Rt(min)", 'Average_Mz', 'Matched_peaks_percentage', 'Metabolite_name', 'MS.MS_assigned', 'INCHIKEY', 'Formula',
                                                     'SMILES', 'MS.MS_matched', 'm.z_similarity', 'Simple_dot_product', 'Weighted_dot_product', 'Reverse_dot_product', 'Matched_peaks_count', 'Total_score', 'MS1_isotopic_spectrum', 
                                                     'MS.MS_spectrum', merge_blanks_by)) %>% #merge samples and blanks
  dplyr::mutate(Blank_Sub = ifelse(grepl("_D_", Sample),
                                        Sample_Area - Avg_PDMS_Blank,
                                        ifelse(grepl("_G_", Sample),
                                               Sample_Area - Avg_GFF_Blank,
                                               Sample_Area - Avg_Active_Blank))) %>% #subtract appropriate blank
  dplyr::mutate(Blank_Ratio = ifelse(grepl("_D_", Sample),
                                   Avg_PDMS_Blank/Sample_Area,
                                   ifelse(grepl("_G_", Sample),
                                          Avg_GFF_Blank/Sample_Area,
                                          Avg_Active_Blank/Sample_Area))) #determine blank-to-signal ratio

sanity_check <- nrow(Samples_Blank_Merge)-nrow(Samples) 
if(sanity_check != 0){
  warning('Check samples and blanks')
}

#-------filter out high blanks-------
#filter to retain blanks below the blank threshold
Samples_Blank_Merge_Filtered <- Samples_Blank_Merge %>%
  filter(Blank_Ratio < Blank_Threshold)

Percent_detection_filter <- Samples_Blank_Merge_Filtered %>%
  select(-c(all_of(merge_blanks_by), Blank_Area, Avg_PDMS_Blank, Avg_GFF_Blank,
            Avg_Active_Blank, Blank_Ratio, Sample_Area, Blank_Sample)) %>% #remove columns that differ for different samples per feature (interferes with pivot wider)
  mutate(Blank_Sub = as.numeric(Blank_Sub)) %>% #make sure values are numbers
  pivot_wider(names_from = Sample, values_from = Blank_Sub)  #pivot

pivot_check1 <- nrow(Percent_detection_filter) - length(unique(Percent_detection_filter$Alignment_ID))
if(pivot_check1 != 0){
  warning('Check pivot')
} #ensure pivot worked correctly; if failed, there are still columns that differ#

#-------calculate percent detection-------

all_samples_deployed <- colnames(Percent_detection_filter[,grep(Sample_String, colnames(Percent_detection_filter))]) #get list of all samples

#get LOD, detection count, detection frequency and mean
#note: mean excludes non-detects and will be used in later filtering
Percent_detection_filter[,'LOD'] <- (apply(Percent_detection_filter[,which(colnames(Percent_detection_filter) %in% all_samples_deployed)], 1, min, na.rm = TRUE))/5
Percent_detection_filter[,'Detection_Count'] = rowSums(!is.na(Percent_detection_filter[grepl(Sample_String, colnames(Percent_detection_filter))])) #detection  for all sa,[;es]
Percent_detection_filter[,'Detection_Frequency'] = Percent_detection_filter$Detection_Count*100/length(all_samples_deployed) #calculate detection frequency
Percent_detection_filter[,'Mean'] <- apply(Percent_detection_filter[,which(colnames(Percent_detection_filter) %in% all_samples_deployed)], 1, mean, na.rm = T) #mean abundnace (excluding non-detects)

Percent_detection_filter <- Percent_detection_filter %>%
  mutate(Mean = ifelse(is.na(Mean), 0, Mean)) #if no detections (i.e., mean = NA), replace w 0

#-------get level2 before detection frequency-------
level2_nodetfreq <- annotate(Percent_detection_filter, formulas, structures, unreal_level2 = unreal_level2) %>%
  filter(level == 2) %>%
  dedup(.) %>% #deuplicate
  ungroup() 
#write_xlsx(level2_nodetfreq,'level2_nodetfreq_nooff.xlsx')

smi_data_lvl2 <- convert_to_smi_tab_delimited(level2_nodetfreq, label = 'Alignment_ID')

#Write the tab-delimited .smi file to get opera values
writeLines(smi_data_lvl2, "OPERA Code Exports/lvl2_OPERA_SMILES.smi")

#see level 2s before deduplicating to consdier other peaks as possiblities
level2_nodetfreq_nodedup <- annotate(Percent_detection_filter, formulas, structures, unreal_level2 = c(0)) %>%
  filter(level == 2) %>%
  group_by(INCHIKEY) %>%
  filter(n()>1) %>%
  ungroup(.) 

#import level 2 opera to view and add to tbale
#level2_opera <- read_csv('lvl2_OPERA_SMILES-smi_OPERA2.9Pred.csv') %>% 
  #dplyr::rename(`2D_INCHIKEY` = MoleculeID) %>%
  #merge(., level2_nodetfreq, by = c('2D_INCHIKEY'))

#-------filter for detection frequency, replace non-detects--------
#filter for highly detected features
Starting_Features_Detfreq_50 <- Percent_detection_filter %>%
  filter(Detection_Frequency >= detection_freq_threshold) 

#replace nondetects with LOD
Starting_Features_Prelim <- Starting_Features_Detfreq_50 %>%
  pivot_longer(cols = all_of(all_samples_deployed), names_to = 'Sample', values_to = 'Area', values_drop_na = F) %>%
  mutate(Area = ifelse(is.na(Area), LOD, Area)) %>% 
  pivot_wider(names_from = Sample , values_from = Area) 

sanity_check2 <- nrow(Starting_Features_Prelim)-(nrow(Starting_Features_Detfreq_50)) #make sure no data was lost between pivoting
if(sanity_check2 != 0){
  warning('Error in pivot')
}

#-------deduplicate degenerate features-------
empirical_features_list <- get_empirical_features(Starting_Features_Prelim, khipu_results_clean) %>% #get list of empirical features
  separate(col = Metabolite_name, into = c("discard", "Metabolite_name"), sep = ": ", extra = 'merge', fill = 'left') %>% #clean up metabolite names to get rid of the "No MS2: ", etc. in front of names
  subset(select = -c(discard)) #remove column not needed

#-------assign annotation confidence levels-------
#goal: annotate features
annotated <-annotate(empirical_features_list, formulas, structures, unreal_level2 = unreal_level2) 

#-------deduplicate-------
#goal: deduplicate level 2 and 3 features with the same annotation
annotated_dedup <- dedup(annotated)

Level2_Together <- annotated_dedup %>%
  filter(level == 2) %>%
  select(c(Alignment_ID, Metabolite_name, Average_Mz, Formula))

#Exports a data frame as a .csv file#
DF_to_CSV <- function(df){
  write.csv(df, paste(deparse(substitute(df)), ".csv", sep = ""))
}

DF_to_CSV(Level2_Together)

#-------visualize annotations for detection frequency across all compounds-------
dedup_features_bar <- annotated_dedup %>%
  group_by(level) %>%
  mutate(level = as.character(level)) %>%
  mutate(level = ifelse(level ==2, '2a/2b', level)) %>%
  dplyr::summarise(n = n()) %>%
  ggplot(aes(x = level, y = n)) +
  geom_bar(fill = 'grey', stat = "identity", width = 0.7)+
  labs(x = "Level", y = "Number of Compounds") +  
  theme_classic()+
  theme(axis.title.x = element_text(size = 25), axis.title.y = element_text(size = 25), axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20), legend.title = element_text(size = 20), legend.text = element_text(size = 20),
        title = element_text(size = 25), strip.text = element_text(size = 25), strip.background = element_rect(colour="black", fill="white"))+ #adjust text size#
  geom_text(aes(label = n), position = position_dodge(width = 0.9), vjust = -0.25, size = 6) # add n on each bar
dedup_features_bar
save_plot(plot = dedup_features_bar, 'Figures/dedup_features_bar.svg', base_width = 8, base_height = 7)

#-------analysis separating by sampler type-------
Sampler_type_separate <- Percent_detection_filter

Active <- colnames(Sampler_type_separate[,grep('_A_', colnames(Sampler_type_separate))]) #get list of all samples
PDMS <- colnames(Sampler_type_separate[,grep('_D_', colnames(Sampler_type_separate))]) #get list of all samples
GFF <- colnames(Sampler_type_separate[,grep('_G_', colnames(Sampler_type_separate))]) #get list of all samples

Sampler_type_separate[,'Detection_Count_A'] = rowSums(!is.na(Sampler_type_separate[grepl('_A_', colnames(Sampler_type_separate))]))
Sampler_type_separate[,'Detection_Frequency_A'] = Sampler_type_separate$Detection_Count_A*100/length(Active)
Sampler_type_separate[,'Mean_A'] <- apply(Sampler_type_separate[,which(colnames(Sampler_type_separate) %in% Active)], 1, mean, na.rm = T)

Sampler_type_separate[,'Detection_Count_D'] = rowSums(!is.na(Sampler_type_separate[grepl('_D_', colnames(Sampler_type_separate))]))
Sampler_type_separate[,'Detection_Frequency_D'] = Sampler_type_separate$Detection_Count_D*100/length(PDMS)
Sampler_type_separate[,'Mean_D'] <- apply(Sampler_type_separate[,which(colnames(Sampler_type_separate) %in% PDMS)], 1, mean, na.rm = T)

Sampler_type_separate[,'Detection_Count_G'] = rowSums(!is.na(Sampler_type_separate[grepl('_G_', colnames(Sampler_type_separate))]))
Sampler_type_separate[,'Detection_Frequency_G'] = Sampler_type_separate$Detection_Count_G*100/length(GFF)
Sampler_type_separate[,'Mean_G'] <- apply(Sampler_type_separate[,which(colnames(Sampler_type_separate) %in% GFF)], 1, mean, na.rm = T)

#Identifies the number of days a chemical was detected for filtering#
Count_num_days <- function(df){
  df_A <- df%>% #start with active
    filter(Sampler == "Active") %>%
    filter(Area > LOD) %>%
    group_by(Alignment_ID) %>%
    dplyr::mutate(Days_Detected_A = length(unique(Batch))) %>% #counts the number of days a feature was detected on
    slice(1) %>% #value will be duplicated, remove duplicates
    ungroup() %>%
    select(c(Alignment_ID, Days_Detected_A))
  
  df_D <- df%>% #repeat for PDMS
    filter(Sampler == "PDMS-Passive") %>%
    filter(Area > LOD) %>%
    group_by(Alignment_ID) %>%
    dplyr::mutate(Days_Detected_D = length(unique(Batch))) %>%
    slice(1) %>%
    ungroup() %>%
    select(c(Alignment_ID, Days_Detected_D))
  
  df_G <- df%>% #repeat for GFF
    filter(Sampler == "GFF-Passive") %>%
    filter(Area > LOD) %>%
    group_by(Alignment_ID) %>%
    dplyr::mutate(Days_Detected_G = length(unique(Batch))) %>%
    slice(1) %>%
    ungroup() %>%
    select(c(Alignment_ID, Days_Detected_G))
  
  df2 <- merge(df, df_A, by = "Alignment_ID", all.x = TRUE) %>%
    merge(., df_D, by = "Alignment_ID", all.x = TRUE) %>%
    merge(., df_G, by = "Alignment_ID", all.x = TRUE) %>% #combine
    mutate(Days_Detected_A = ifelse(is.na(Days_Detected_A), 0, Days_Detected_A)) %>%
    mutate(Days_Detected_D = ifelse(is.na(Days_Detected_D), 0, Days_Detected_D)) %>%
    mutate(Days_Detected_G = ifelse(is.na(Days_Detected_G), 0, Days_Detected_G)) %>% #replace any NA means with 0
    select(-Batch) #Batch no longer required, will cause issues with pivoting if not removed

  return(df2)
}

Separate_filter <-  Sampler_type_separate %>%
  mutate(Mean_D = ifelse(is.na(Mean_D), 0, Mean_D)) %>% 
  mutate(Mean_A = ifelse(is.na(Mean_A), 0, Mean_A)) %>%
  mutate(Mean_G = ifelse(is.na(Mean_G), 0, Mean_G)) %>% #replace any NA means with 0
  filter(Detection_Frequency_D >= 50 | Detection_Frequency_A >= 50 |
           Detection_Frequency_G >= 50) %>% #filter for detection frequency
  pivot_longer(cols = all_of(all_samples_deployed), names_to = 'Sample', values_to = 'Area', values_drop_na = F) %>%
  mutate(Area = ifelse(is.na(Area), LOD, Area)) %>% #change to LOD (this is across all samples) for things that are non-detects
  label_samples() %>% #add identifying information
  Count_num_days() %>% #apply day counting function
  filter(Days_Detected_D >= 4 | Days_Detected_A >= 4 |
           Days_Detected_G >= 4) %>% #filter for at least 4 time points
  select(-c(all_of(merge_blanks_by), Sample_Type)) %>% #remove extraneous columns
  pivot_wider(names_from = Sample , values_from = Area) %>% #return to wide form
  separate(col = Metabolite_name, into = c("discard", "Metabolite_name"), sep = ": ", extra = 'merge', fill = 'left') #clean up metabolite names to get rid of the "No MS2: ", etc. in front of names

pivot_check2 <- nrow(Separate_filter) - length(unique(Separate_filter$Alignment_ID))
if(pivot_check2 != 0){
  warning('Check pivot')
} #ensure pivot worked correctly; if failed, there are still columns that differ#

#finish annotation and filtering like before
Separate_emp <- get_empirical_features(Separate_filter, khipu_results_clean)
Separate_annotated <- annotate(Separate_emp, formulas, structures, unreal_level2 = unreal_level2)
Separate_dedup <- dedup(Separate_annotated)

#tests; can be ignored
test <- Separate_dedup %>% filter(level >3 & MS.MS_assigned == 'True') %>% mutate(type = 'not annotated')
ksdgsdg <-Separate_dedup %>% filter(level < 4) %>% mutate(type = 'annotated')

tegesg <- rbind(test, ksdgsdg) %>%
  ggplot(aes(x = Average_Mz, fill = type))+
  geom_density(alpha = 0.5)

#visualize the annotation levels
separate_dedup_features_bar <- Separate_dedup %>%
  group_by(level) %>%
  mutate(level = as.character(level)) %>%
  mutate(level = ifelse(level ==2, '2a/2b', level)) %>%
  dplyr::summarise(n = n()) %>%
  ggplot(aes(x = level, y = n)) + #get the number of features annotated per level
  geom_bar(fill = 'grey', stat = "identity", width = 0.7)+
  labs(x = "Level", y = "Number of Compounds") +  
  theme_classic()+
  theme(axis.title.x = element_text(size = 25), axis.title.y = element_text(size = 25), axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20), legend.title = element_text(size = 20), legend.text = element_text(size = 20),
        title = element_text(size = 25), strip.text = element_text(size = 25), strip.background = element_rect(colour="black", fill="white"))+ #adjust text size#
  geom_text(aes(label = n), position = position_dodge(width = 0.9), vjust = -0.25, size = 6) # add n on each bar
separate_dedup_features_bar

#-------compare feature profiles between sampler type-------
#1. deduplicate within each sampler type separately
Dedup_D <- Separate_annotated %>%
  filter(Detection_Frequency_D >= 50 & Days_Detected_D >= 4) %>% #filter for 100% detection
  select(-c(Mean)) %>% #remove to rename
  dplyr::rename(Mean = Mean_D) %>% #this is renamed for deduplicate to deduplicate by mean for this specific location/lens
  dedup(.) %>% #deduplicate
  mutate(Alignment_ID_INCHIKEY = Alignment_ID) %>% #preserve original alignment id before replacing it
  mutate(Alignment_ID_INCHIKEY = ifelse(!is.na(`2D_INCHIKEY`), `2D_INCHIKEY`, Alignment_ID)) %>% #so can align level 4/5 between groups by alignment id, but level 2 and 3 by inchikey to account for deduplication differences
  group_by(Alignment_ID_INCHIKEY) %>%
  slice(1) %>% #to address duplicates remaining after filtering when mean is the same
  mutate(Sampler = 'PDMS') %>%
  ungroup() %>%
  dplyr::rename(Mean_D = Mean)

Dedup_G <- Separate_annotated %>%
  filter(Detection_Frequency_G >= 50 & Days_Detected_G >= 4) %>%
  select(-c(Mean)) %>%
  dplyr::rename(Mean = Mean_G) %>%
  dedup(.) %>%
  mutate(Alignment_ID_INCHIKEY = Alignment_ID) %>%
  mutate(Alignment_ID_INCHIKEY = ifelse(!is.na(`2D_INCHIKEY`), `2D_INCHIKEY`, Alignment_ID)) %>%
  group_by(Alignment_ID_INCHIKEY) %>%
  slice(1) %>%  #to address duplicates remaining after filtering when mean is the same
  mutate(Sampler = 'GFF') %>%
  ungroup() %>%
  dplyr::rename(Mean_G = Mean)
  
Dedup_A <- Separate_annotated %>%
  filter(Detection_Frequency_A >= 50 & Days_Detected_A >= 4) %>%
  select(-c(Mean)) %>%
  dplyr::rename(Mean = Mean_A) %>%
  dedup(.) %>%
  mutate(Alignment_ID_INCHIKEY = Alignment_ID) %>%
  mutate(Alignment_ID_INCHIKEY = ifelse(!is.na(`2D_INCHIKEY`), `2D_INCHIKEY`, Alignment_ID)) %>% #so can align level 4/5 by alignment id, but level 2 and 3 by inchikey to account for deduplication differences
  group_by(Alignment_ID_INCHIKEY) %>%
  slice(1) %>% #to address duplicates remaining after filtering when mean is the same
  mutate(Sampler = 'Active') %>%
  ungroup() %>%
  dplyr::rename(Mean_A = Mean)

#2. bind back into one data frame
Seperate_dedup_individ_bound <- rbind(Dedup_D, Dedup_G, Dedup_A) #combine into one

Total_Number_of_Features <- length(unique(Seperate_dedup_individ_bound$Alignment_ID))

#3. visualize as bar graphs
features_per_replicate <- Seperate_dedup_individ_bound %>%
  group_by(Sampler) %>%
  dplyr::summarise(n = n()) %>%
  ggplot(aes(x = reorder(Sampler, -n), y = n)) +
  geom_bar(fill = '#D95F02', stat = "identity", width = 0.7, alpha = 0.7)+
  labs(x = "", y = "Number of Features") +  
  theme_classic()+
  # geom_point(size = 5) + 
  #  geom_segment( aes(x=Sampler, xend=Sampler, y=0, yend=n), linewidth = 2)+
  theme(axis.title.x = element_text(size = 30),axis.text.x=element_text(size=25),
        axis.text.y=element_text(size=20), axis.title.y=element_text(size=25), legend.text=element_text(size=20),
        legend.title=element_text(size=20))+ 
  # theme(plot.margin = unit(c(1, 0, 0, 0.75),"inches"))+
  geom_text(aes(label = n), position = position_dodge(width = 0.9), colour = 'black', vjust = -0.5, size = 8)+ # add n on each bar
  theme(axis.text.x = element_text(angle = 270, hjust=-0.05)) #rotate x axis labels
features_per_replicate

#-------venn diagram for coated and uncoated-------
#we are comparing features detected between coating types so use the one with deduplication within sampler types
Seperate_dedup_individ_bound$Sampler <- factor(Seperate_dedup_individ_bound$Sampler,
                                levels = c("PDMS", "GFF", "Active"),
                                labels = c("PDMS", "GFF", "Active")) #reorder layers#

venn_data <- Seperate_dedup_individ_bound %>%
  unstack(Alignment_ID ~ Sampler) #separate alignment ids/inchikeys present by sampler type#

Accent_Palette = brewer.pal(8, "Accent") #colours for venn diagram
Sampler_Colours = Accent_Palette[c(1,2,4)]
UR_Colours = Accent_Palette[c(1,2,3)]

#visualize as venn diagram
sampler_type_ven <- ggvenn(venn_data, stroke_size = 0.7,fill_alpha = 0.7,count_column = TRUE,
                          fill_color = Sampler_Colours, set_name_size = 10, text_size = 6, show_percentage = TRUE)+ #venn diagram between sampler types
  theme(legend.position="none")+
  theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"))+
  theme(legend.title = element_text(size = 20),
        legend.text = element_text(size = 18),
        title = element_text(size = 25),
        plot.background = element_rect(fill = "white", colour = "white"),
        panel.border = element_blank(),
        strip.text = element_text(size = 18)) # Customize theme
sampler_type_ven

save_plot(plot = sampler_type_ven, 'Figures/Figure_S2.svg', base_height = 10, base_width = 13)

Feature_Levels_Passive <- Seperate_dedup_individ_bound %>%
  filter(Sampler != "Active") %>%
  group_by(level, Sampler, Alignment_ID) %>%
  slice(1) %>%
  ungroup() %>%
  group_by(level, Sampler) %>%
  mutate(level = as.character(level)) %>%
  mutate(level = ifelse(level == 2, '2a/2b', level)) %>%
  mutate(level = ifelse(level > 3, '4/5', level)) %>%
  dplyr::mutate(n = n()) %>%
  ggplot(aes(x = level, fill = Sampler)) +
  geom_bar(width = 0.7, position = position_dodge(width = 0.8)) +
  labs(x = "Confidence Level", y = "Number of Features", fill = "PAS Sorbent") +
  scale_y_continuous(expand = c(0, 0, 0.1, 0.1))+ #make the bars flush to the x axis#
  theme_classic()+
  theme(axis.title.x = element_text(colour = "black", size = 25), 
        axis.title.y = element_text(colour = "black", size = 25), #adjust text size for axis titles#
        axis.text.x = element_text(colour = "black", size = 20),
        axis.text.y = element_text(colour = "black", size = 20), #adjust axis text size and colour#
        plot.title = element_text(size = 15), #adjust plot title text size#
        axis.ticks = element_line(colour="black", linewidth=1), #adjust tick mark size and colour#
        axis.line = element_line(colour="black", linewidth=1), #adjust axis size and colour#
        panel.border = element_blank(), #remove border#
        panel.grid.major=element_blank(), #remove gridlines#
        panel.grid.minor=element_blank(), #remove gridlines#
        #legend.title = element_text(colour = "black", size = 15),
        #legend.text = element_text(colour = "black", size = 15), #adjust legend title and text#
        strip.text = element_text(colour = "black", size = 20), #adjust strip text if faceting#
        strip.background = element_rect(colour="black", fill="white")) +
  geom_text(aes(y = n + 200, label = n), position = position_dodge(width = 0.9), vjust = 0.5, size = 6) # add n on each bar
Feature_Levels_Passive

Feature_Levels_Active <- Seperate_dedup_individ_bound %>%
  filter(Sampler == "Active") %>%
  group_by(level, Sampler, Alignment_ID) %>%
  slice(1) %>%
  ungroup() %>%
  group_by(level, Sampler) %>%
  mutate(level = as.character(level)) %>%
  mutate(level = ifelse(level == 2, '2a/2b', level)) %>%
  mutate(level = ifelse(level > 3, '4/5', level)) %>%
  dplyr::mutate(n = n()) %>%
  ggplot(aes(x = level)) +
  geom_bar(width = 0.7, show.legend = FALSE, fill = "green") +
  labs(x = "Confidence Level", y = "Number of Features") +
  scale_y_continuous(expand = c(0, 0, 0.1, 0.1))+ #make the bars flush to the x axis#
  theme_classic()+
  theme(axis.title.x = element_text(colour = "black", size = 25), 
        axis.title.y = element_text(colour = "black", size = 25), #adjust text size for axis titles#
        axis.text.x = element_text(colour = "black", size = 20),
        axis.text.y = element_text(colour = "black", size = 20), #adjust axis text size and colour#
        plot.title = element_text(size = 15), #adjust plot title text size#
        axis.ticks = element_line(colour="black", linewidth=1), #adjust tick mark size and colour#
        axis.line = element_line(colour="black", linewidth=1), #adjust axis size and colour#
        panel.border = element_blank(), #remove border#
        panel.grid.major=element_blank(), #remove gridlines#
        panel.grid.minor=element_blank(), #remove gridlines#
        #legend.title = element_text(colour = "black", size = 15),
        #legend.text = element_text(colour = "black", size = 15), #adjust legend title and text#
        strip.text = element_text(colour = "black", size = 20), #adjust strip text if faceting#
        strip.background = element_rect(colour="black", fill="white")) +
  geom_text(aes(y = n + 200, label = n), position = position_dodge(width = 0.9), vjust = 0.5, size = 6) # add n on each bar
Feature_Levels_Active

Feature_Levels_All <- Seperate_dedup_individ_bound %>%
  group_by(level, Sampler, Alignment_ID) %>%
  slice(1) %>%
  ungroup() %>%
  group_by(level, Sampler) %>%
  mutate(level = as.character(level)) %>%
  mutate(level = ifelse(level == 2, '2a/2b', level)) %>%
  mutate(level = ifelse(level > 3, '4/5', level)) %>%
  dplyr::mutate(n = n()) %>%
  ggplot(aes(x = level, fill = Sampler)) +
  geom_bar(width = 0.7, position = position_dodge(width = 0.8)) +
  labs(x = "Confidence Level", y = "Number of Features", fill = "Sampler Type") +
  scale_fill_manual(values = Sampler_Colours)+
  scale_y_continuous(expand = c(0, 0, 0.1, 0.1))+ #make the bars flush to the x axis#
  theme_classic()+
  theme(axis.title.x = element_text(colour = "black", size = 25), 
        axis.title.y = element_text(colour = "black", size = 25), #adjust text size for axis titles#
        axis.text.x = element_text(colour = "black", size = 20),
        axis.text.y = element_text(colour = "black", size = 20), #adjust axis text size and colour#
        plot.title = element_text(size = 15), #adjust plot title text size#
        axis.ticks = element_line(colour="black", linewidth=1), #adjust tick mark size and colour#
        axis.line = element_line(colour="black", linewidth=1), #adjust axis size and colour#
        panel.border = element_blank(), #remove border#
        panel.grid.major=element_blank(), #remove gridlines#
        panel.grid.minor=element_blank(), #remove gridlines#
        legend.title = element_text(colour = "black", size = 15),
        legend.text = element_text(colour = "black", size = 15), #adjust legend title and text#
        legend.position = "bottom",
        strip.text = element_text(colour = "black", size = 20), #adjust strip text if faceting#
        strip.background = element_rect(colour="black", fill="white")) +
  geom_text(aes(y = n + 200, label = n), position = position_dodge(width = 0.8), vjust = 0.5, size = 6) # add n on each bar
Feature_Levels_All

save_plot(plot = Feature_Levels_All, 'Figures/Figure_S3.svg', base_height = 9, base_width = 8)

####-----Calculating equivalent air volumes-----####

Volumes <- read_csv("C://Users//scrol//OneDrive - McMaster University//Joseph Okeme's files - Alyss//New passives calibration//Celine_Ellie_Sample_filled_2.csv", col_types = cols(.default = "c"))

#Reformats the file containing air volumes#
clean_sample_list <- function(df){
  df <- as.data.frame(df)
  names(df) <- gsub(" ", "_", names(df)) #replace spaces in column names with underscore#
  names(df) <- gsub("/", ".", names(df)) #replace slashes in column names with period#
  names(df) <- gsub("-", "_", names(df)) #replace dashes in column names with underscore#
  df <- df %>% as.data.frame(row.names = 1:nrow(.)) #reset row number#
  return(df) #so output is a data frame#
}

#Create a data frame containing the air volumes collected by the active samplers#
AirVolumes <- clean_sample_list(Volumes)%>%
  filter(grepl("_A", Sample_ID) & grepl("sample", Sample_Type))%>% #select only active samplers#
  select(-c(Volume_of_Air_for_Active, Injeciton_Vol, Matrix_Type,
            Sample_Type, Deployment_period))%>% #remove irrelevant information#
  dplyr::rename ("Active_Air_Volume" = Volume_of_Air_for_AAS)%>% #rename for simplicity#
  mutate(Deployment_Time = as.numeric(Deployment_Time),
         Batch = as.numeric(Batch),
         Active_Air_Volume = as.numeric(Active_Air_Volume))%>% #convert to numeric#
  mutate("L_per_Day" = Active_Air_Volume/Deployment_Time)%>% #create column for average daily air intake#
  mutate(Sample_ID = gsub("30819", "", Sample_ID))%>%
  mutate(Sample_ID = gsub("JOQT", "JO", Sample_ID))%>% #rename samples to prevent filtering problems#
  dplyr::rename("Sample" = Sample_ID) %>%
  mutate(Sample = ifelse(grepl("Q", Sample) & !grepl("p", Sample),
                         gsub("JO_A", "JO_L_MB_A_", Sample),
                         ifelse(grepl("q", Sample) & grepl("p", Sample),
                                gsub("JO_p", "JO_L_MB_P_D_", Sample),
                                ifelse(grepl("S2", Sample),
                                       gsub("JO_S2", "JO_L_SB_SP_S2_", Sample),
                                       ifelse(grepl("ACN2", Sample),
                                              gsub("JO_ACN2", "JO_L_SB_UP_ACN2_", Sample),
                                              ifelse(grepl("bl", Sample),
                                                     gsub("JO_A", "JO_F_Bk_A_", Sample),
                                                     ifelse(grepl("p0", Sample) & !grepl("d", Sample) & !grepl("e", Sample) & !grepl("f", Sample),
                                                            gsub("JO_p0", "JO_F_Bk_P_D_", Sample),
                                                            ifelse(grepl("p0", Sample) & !grepl("b", Sample) & !grepl("c", Sample),
                                                                   gsub("JO_p0", "JO_F_Bk_P_G_", Sample),
                                                                   ifelse(grepl("A", Sample),
                                                                          gsub("JO_A", "JO_F_S_A_", Sample),
                                                                          ifelse(grepl("p", Sample) & !grepl("d", Sample) & !grepl("e", Sample) & !grepl("f", Sample),
                                                                                 gsub("JO_p", "JO_F_S_P_D_", Sample),
                                                                                 gsub("JO_p", "JO_F_S_P_G_", Sample)))))))))))
#Simplified version of previous data frame
AirVolumes2 <- AirVolumes %>%
  group_by(Batch) %>%
  dplyr::mutate(AvgM3 = mean(Active_Air_Volume)/1000) %>%
  ungroup() %>%
  dplyr::mutate(TotalM3 = sum(unique(AvgM3)))

#Create a data frame for the deployment duration and batch of the PDMS PAS#
PDMS_Time <- clean_sample_list(Volumes)%>%
  mutate(Sample_ID = gsub("30819", "", Sample_ID))%>%
  mutate(Sample_ID = gsub("JOQT", "JO", Sample_ID))%>% #rename samples to prevent filtering problems#
  dplyr::rename("Sample" = Sample_ID) %>%
  mutate(Sample = ifelse(grepl("Q", Sample) & !grepl("p", Sample),
                         gsub("JO_A", "JO_L_MB_A_", Sample),
                         ifelse(grepl("q", Sample) & grepl("p", Sample),
                                gsub("JO_p", "JO_L_MB_P_D_", Sample),
                                ifelse(grepl("S2", Sample),
                                       gsub("JO_S2", "JO_L_SB_SP_S2_", Sample),
                                       ifelse(grepl("ACN2", Sample),
                                              gsub("JO_ACN2", "JO_L_SB_UP_ACN2_", Sample),
                                              ifelse(grepl("bl", Sample),
                                                     gsub("JO_A", "JO_F_Bk_A_", Sample),
                                                     ifelse(grepl("p0", Sample) & !grepl("d", Sample) & !grepl("e", Sample) & !grepl("f", Sample),
                                                            gsub("JO_p0", "JO_F_Bk_P_D_", Sample),
                                                            ifelse(grepl("p0", Sample) & !grepl("b", Sample) & !grepl("c", Sample),
                                                                   gsub("JO_p0", "JO_F_Bk_P_G_", Sample),
                                                                   ifelse(grepl("A", Sample),
                                                                          gsub("JO_A", "JO_F_S_A_", Sample),
                                                                          ifelse(grepl("p", Sample) & !grepl("d", Sample) & !grepl("e", Sample) & !grepl("f", Sample),
                                                                                 gsub("JO_p", "JO_F_S_P_D_", Sample),
                                                                                 gsub("JO_p", "JO_F_S_P_G_", Sample))))))))))) %>%
  filter(Sample_Type == "sample")%>% #select only samples#
  filter(grepl("_p", Sample) & !grepl("repeat", Sample)
         & !grepl("HV", Sample))%>% #select only PAS samples and remove the repeats#
  filter(grepl("a", Sample) | grepl("b", Sample)
         | grepl("c", Sample))%>% #select only the PDMS replicates#
  mutate(Deployment_Time = as.numeric(Deployment_Time),
         Batch = as.numeric(Batch))%>% #convert to numeric#
  select(c(Sample, Deployment_Time, Batch)) #retain only important columns; adjust as necessary#

#Do the same for GFF#
GFF_Time <- clean_sample_list(Volumes)%>%
  mutate(Sample_ID = gsub("30819", "", Sample_ID))%>%
  mutate(Sample_ID = gsub("JOQT", "JO", Sample_ID))%>% #rename samples to prevent filtering problems#
  dplyr::rename("Sample" = Sample_ID) %>%
  mutate(Sample = ifelse(grepl("Q", Sample) & !grepl("p", Sample),
                         gsub("JO_A", "JO_L_MB_A_", Sample),
                         ifelse(grepl("q", Sample) & grepl("p", Sample),
                                gsub("JO_p", "JO_L_MB_P_D_", Sample),
                                ifelse(grepl("S2", Sample),
                                       gsub("JO_S2", "JO_L_SB_SP_S2_", Sample),
                                       ifelse(grepl("ACN2", Sample),
                                              gsub("JO_ACN2", "JO_L_SB_UP_ACN2_", Sample),
                                              ifelse(grepl("bl", Sample),
                                                     gsub("JO_A", "JO_F_Bk_A_", Sample),
                                                     ifelse(grepl("p0", Sample) & !grepl("d", Sample) & !grepl("e", Sample) & !grepl("f", Sample),
                                                            gsub("JO_p0", "JO_F_Bk_P_D_", Sample),
                                                            ifelse(grepl("p0", Sample) & !grepl("b", Sample) & !grepl("c", Sample),
                                                                   gsub("JO_p0", "JO_F_Bk_P_G_", Sample),
                                                                   ifelse(grepl("A", Sample),
                                                                          gsub("JO_A", "JO_F_S_A_", Sample),
                                                                          ifelse(grepl("p", Sample) & !grepl("d", Sample) & !grepl("e", Sample) & !grepl("f", Sample),
                                                                                 gsub("JO_p", "JO_F_S_P_D_", Sample),
                                                                                 gsub("JO_p", "JO_F_S_P_G_", Sample))))))))))) %>%
  filter(Sample_Type == "sample")%>% #select only samples#
  filter(grepl("_p", Sample) & !grepl("repeat", Sample)
         & !grepl("HV", Sample))%>% #select only PAS samples and remove the repeats#
  filter(!grepl("a", Sample) | !grepl("b", Sample)
         | !grepl("c", Sample))%>% #remove the PDMS replicates#
  mutate(Deployment_Time = as.numeric(Deployment_Time),
         Batch = as.numeric(Batch))%>% #convert to numeric#
  select(c(Sample, Deployment_Time, Batch)) #retain only important columns; adjust as necessary#

#From the earlier de-duplicated data, prepare for further processing#

#Create long versions of the dedup tables#
Dedup_A_long <- Dedup_A %>%
  pivot_longer(cols = all_of(all_samples_deployed), names_to = 'Sample', values_to = 'Area', values_drop_na = F) %>%
  filter(grepl("_A_", Sample)) %>%
  filter(Area > LOD)

Dedup_D_long <- Dedup_D %>%
  pivot_longer(cols = all_of(all_samples_deployed), names_to = 'Sample', values_to = 'Area', values_drop_na = F) %>%
  filter(grepl("_D_", Sample)) %>%
  filter(Area > LOD)

Dedup_G_long <- Dedup_G %>%
  pivot_longer(cols = all_of(all_samples_deployed), names_to = 'Sample', values_to = 'Area', values_drop_na = F) %>%
  filter(grepl("_G_", Sample)) %>%
  filter(Area > LOD)

#Reformatting#
uActive <- Dedup_A_long%>%
  mutate(Alignment_ID = as.numeric(Alignment_ID),
         Average_Mz = as.numeric(Average_Mz))%>% #convert to numeric#
  add_batch() %>%
  select(-c(Total_score, Detection_Count_A, Detection_Count_D, Detection_Count_G, 
            Detection_Frequency_A, Detection_Frequency_D, Detection_Frequency_G,
            Mean_A, Mean_D, Mean_G)) #remove extraneous information#

#Repeat for PDMS#
uPDMS <- Dedup_D_long%>%
  mutate(Alignment_ID = as.numeric(Alignment_ID),
         Average_Mz = as.numeric(Average_Mz))%>%
  add_batch() %>%
  select(-c(Total_score, Detection_Count_A, Detection_Count_D, Detection_Count_G, 
            Detection_Frequency_A, Detection_Frequency_D, Detection_Frequency_G,
            Mean_D, Mean_A, Mean_G))

#Repeat for GFF#
uGFF <- Dedup_G_long%>%
  mutate(Alignment_ID = as.numeric(Alignment_ID),
         Average_Mz = as.numeric(Average_Mz))%>%
  add_batch() %>%
  select(-c(Total_score, Detection_Count_A, Detection_Count_D, Detection_Count_G, 
            Detection_Frequency_A, Detection_Frequency_D, Detection_Frequency_G,
            Mean_G, Mean_A, Mean_D))

##Because AAS were deployed in sets throughout the PAS deployment period,
##only AAS samples that were collected while a particular PAS was also
##deployed should be used to calculate air volume for that PAS; the following
##code calculates a concentration in U/L, where U represents arbitrary units for
##peak area, from each AAS and averages them into groups based on the deployment times##

#For the first deployment window (11/25/22 to 11/29/22)#
uActive_B1 <- merge(uActive, AirVolumes, by = c("Sample", "Batch"))%>% #add air volumes to the AAS data#
  mutate("Concentration" = Area/Active_Air_Volume)%>% #create new column for Conc by dividing peak area by air volume#
  filter(Batch == 1)%>% #retain only AAS samples collected during first deployment window#
  dplyr::group_by(Alignment_ID)%>% #group by Alignment_ID#
  dplyr::mutate("Avg_Conc_B1" = mean(Concentration))%>% #create a new column that calculates average Conc for batch 1 samples#
  select(c(Alignment_ID, Avg_Conc_B1))%>% #retain only Alignment_ID and average Conc#
  slice(1)%>% #remove repeats#
  ungroup()

#Repeat for the next window (11/25/22 to 11/29/22 and 11/29/22 to 12/03/22)#
uActive_B2 <- merge(uActive, AirVolumes, by = c("Sample", "Batch"))%>%
  mutate("Concentration" = Area/Active_Air_Volume)%>%
  filter(Batch == 1 | Batch == 2)%>% #selecting now the samplers from the first and second window#
  dplyr::group_by(Alignment_ID)%>%
  dplyr::mutate("Avg_Conc_B2" = mean(Concentration))%>%
  select(c(Alignment_ID, Avg_Conc_B2))%>%
  slice(1)%>%
  ungroup()

#Repeat for the next window#
uActive_B3 <- merge(uActive, AirVolumes, by = c("Sample", "Batch"))%>%
  mutate("Concentration" = Area/Active_Air_Volume)%>%
  filter(Batch == 1 | Batch == 2 | Batch == 3)%>%
  dplyr::group_by(Alignment_ID)%>%
  dplyr::mutate("Avg_Conc_B3" = mean(Concentration))%>%
  select(c(Alignment_ID, Avg_Conc_B3))%>%
  slice(1)%>%
  ungroup()

#Repeat for the next window#
uActive_B4 <- merge(uActive, AirVolumes, by = c("Sample", "Batch"))%>%
  mutate("Concentration" = Area/Active_Air_Volume)%>%
  filter(Batch != 5)%>%
  dplyr::group_by(Alignment_ID)%>%
  dplyr::mutate("Avg_Conc_B4" = mean(Concentration))%>%
  select(c(Alignment_ID, Avg_Conc_B4))%>%
  slice(1)%>%
  ungroup()

#Repeat for the next window#
uActive_B5 <- merge(uActive, AirVolumes, by = c("Sample", "Batch"))%>%
  mutate("Concentration" = Area/Active_Air_Volume)%>%
  dplyr::group_by(Alignment_ID)%>%
  dplyr::mutate("Avg_Conc_B5" = mean(Concentration))%>%
  select(c(Alignment_ID, Avg_Conc_B5))%>%
  slice(1)%>%
  ungroup()

#Combine the active concentrations and eliminate those with poor RSD#
AllActiveConc <- join_all(list(uActive_B1, uActive_B2, uActive_B3, uActive_B4,
                               uActive_B5), by = "Alignment_ID", type = "full")%>%
  mutate(Avg_Conc = rowMeans(.[,2:6]),
         Conc_SD = apply(.[2:6], 1, sd, na.rm = TRUE))%>%
  mutate(Conc_RSD = Conc_SD*100/Avg_Conc)%>%
  filter(Conc_RSD <= 20)

#Perform equivalent air volume calculations#

#Exposed surface area of the sampler#
RSA_PDMS <- (0.508*0.508*2)

#Calculate EAV for the PDMS PAS samples#
PDMS_Uptake <- join_all(list(AllActiveConc, uPDMS), by = "Alignment_ID", type = "full") %>% #combine the PDMS data with the AAS concentrations
  mutate("Days" = ifelse(Batch == 1,
                         4,
                         ifelse(Batch == 2,
                                8,
                                ifelse(Batch == 3, 
                                       12,
                                       ifelse(Batch == 4,
                                              17,
                                              20))))) %>% #use batch to create a column denoting how long each sample was exposed#
  mutate("Vol" = (Area / ifelse(Batch == 1,
                                     Avg_Conc_B1,
                                     ifelse(Batch == 2,
                                            Avg_Conc_B2,
                                            ifelse(Batch == 3,
                                                   Avg_Conc_B3,
                                                   ifelse(Batch == 4,
                                                          Avg_Conc_B4,
                                                          Avg_Conc_B5)))))/1000) %>% #calculate PDMS EAV using the appropriate AAS conc#
  mutate(Vol = Vol/RSA_PDMS) %>% #normalize to surface area#
  subset(select = c(Alignment_ID, Average_Mz, Metabolite_name, INCHIKEY, `2D_INCHIKEY`,
                    Alignment_ID_INCHIKEY, Formula, SMILES, level, Sampler, Batch,
                    Days, Vol)) %>%
  filter(!is.na(Vol)) #remove NAs#

#Filter out and amend variability between replicates
EqAirVol_Filtering <- function(df){
  df2 <- df %>%
    group_by(Alignment_ID, Batch) %>% #group replicates
    dplyr::mutate(Avg_Vol = mean(Vol)) %>% #calculate average Veq between replicates
    dplyr::mutate(Num = n()) %>% #number of replicates
    dplyr::mutate(Max = max(Vol)) %>% #highest Veq among replicates
    dplyr::mutate(Min = min(Vol)) %>% #lowest Veq among replicates
    dplyr::mutate(RelRange = (Max - Min)/max(Avg_Vol)*100) %>% #range relative to average Veq
    dplyr::mutate(MaxDiff = Max - Avg_Vol) %>% #difference between the highest Veq and the mean
    dplyr::mutate(MinDiff = Avg_Vol - Min) %>% #difference between the lowest Veq and the mean
    dplyr::mutate(Range_Filt = ifelse((RelRange >= 25) & (MaxDiff > MinDiff) & (Vol == Max) & (Num == 3),
                                      "Remove", #if the range is >=25%, there are three replicates, and the highest Veq is furthest, tag it for removal
                                      ifelse((RelRange >= 25) & (MaxDiff < MinDiff) & (Vol == Min) & (Num == 3),
                                             "Remove", #if the range is >=25%, there are three replicates, and the lowest Veq is furthest, tag it for removal
                                             "Keep"))) %>% #retain everything else
    ungroup() %>%
    filter(Range_Filt == "Keep") %>% #eliminate replicates marked for removal
    group_by(Alignment_ID, Batch) %>% #repeat all previous calculations after removing outlier replicates
    dplyr::mutate(Avg_Vol = mean(Vol)) %>%
    dplyr::mutate(Max = max(Vol)) %>%
    dplyr::mutate(Min = min(Vol)) %>%
    dplyr::mutate(RelRange = (Max - Min)/max(Avg_Vol)*100) %>%
    dplyr::mutate(MaxDiff = Max - Avg_Vol) %>%
    dplyr::mutate(MinDiff = Avg_Vol - Min) %>%
    dplyr::mutate(Range_Filt = ifelse((RelRange >= 25),
                                      "Remove",
                                      "Keep")) %>% #now mark all replicates with range >=25% for removal; this will remove entire time points that are too variable
    ungroup() %>%
    filter(Range_Filt == "Keep") %>% #filter out the tagged replicates
    group_by(Alignment_ID, Batch) %>%
    dplyr::mutate(Avg_Vol = mean(Vol)) %>% #calculate mean Veq again
    ungroup()
    
  df3 <- df2 %>% #filter for signs of equilibration/non-linearity
    select(Alignment_ID, Batch, Avg_Vol) %>%  #only need ID, batch, and average Veq
    group_by(Alignment_ID, Batch, Avg_Vol) %>%
    slice(1) %>% #there will be unnecessary duplicates, keep only one
    ungroup() %>%
    pivot_wider(names_from = Batch, values_from = Avg_Vol) %>%
    filter(!is.na(`5`)) %>% #eliminate any features without a day 20 Veq
    filter(`5` > `4` & `5` > `3`) %>% #eliminate features with day 20 Veq lower than Day 17 or Day 12; indicative of equilibrating/non-linearity
    filter(`1` < `3`) #also eliminate features with Day 12 lower than Day 4 to help select for increasing Veq
  
  df4 <- df2 %>%
    filter(Alignment_ID %in% df3$Alignment_ID) %>% #df3 acts as a list of acceptable features, filter based on it
    group_by(Alignment_ID) %>%
    dplyr::mutate(N_Days = length(unique(Batch))) %>% #number of time points
    dplyr::mutate(N_Samples = n()) %>% #number of individual data points
    ungroup() %>%
    filter((N_Days == 5) & (N_Samples > 7)) #filter for all five time points and at least 7 data points (which amounts to 50% detection)
  
  return(df4)
    
}

#Apply Veq filtering#
PDMS_Uptake_stats <- EqAirVol_Filtering(PDMS_Uptake)

#number of features
length(unique(PDMS_Uptake$Alignment_ID))
length(unique(PDMS_Uptake_stats$Alignment_ID))

#Repeat for GFF#

RSA_GFF <- (0.238*0.238*3.14*2) #exposed sampler surface area#

#Calculate EAV for the GFF PAS samples#
GFF_Uptake <- join_all(list(AllActiveConc, uGFF), by = "Alignment_ID", type = "full")%>% #combine the GFF data with the AAS Concs
  mutate("Days" = ifelse(Batch == 1,
                         4,
                         ifelse(Batch == 2,
                                8,
                                ifelse(Batch == 3, 
                                       12,
                                       ifelse(Batch == 4,
                                              17,
                                              20)))))%>% #use batch to create a column denoting how long each sample was exposed#
  mutate("Vol" = Area / (ifelse(Batch == 1,
                                    Avg_Conc_B1,
                                    ifelse(Batch == 2,
                                           Avg_Conc_B2,
                                           ifelse(Batch == 3,
                                                  Avg_Conc_B3,
                                                  ifelse(Batch == 4,
                                                         Avg_Conc_B4,
                                                         Avg_Conc_B5)))))/1000)%>% #calculate GFF EAV using the appropriate AAS Conc#
  mutate(Vol = Vol/RSA_GFF)%>% #normalize to surface area#
  subset(select = c(Alignment_ID, Average_Mz, Metabolite_name, INCHIKEY, `2D_INCHIKEY`,
                    Alignment_ID_INCHIKEY, Formula, SMILES, level, Sampler, Batch,
                    Days, Vol)) %>%
  filter(!is.na(Vol)) #remove NAs#

#Add mean and sd#
GFF_Uptake_stats <- EqAirVol_Filtering(GFF_Uptake)#%>%
  #filter(Batch < 5 | Alignment_ID %in% Poor_GFF_Curve_List == FALSE)

####-----Linear regression of data-----####

#Perform regression#
Reg_PDMS <- PDMS_Uptake_stats%>%
  group_by(Alignment_ID)%>% #cluster by Alignment_ID#
  do(., tidy(lm(Vol ~ 0 + Days, data = .))) #for each group, calculate regression statistics forcing through the origin#
summary(Reg_PDMS)

#Calculate R-squared#
PDMS_Rsq <- PDMS_Uptake_stats%>%
  split(.$Alignment_ID)%>% #split by Alignment_ID; alternatively, use PDMS_split made previously#
  map(~lm(Vol ~ 0 + Days, data = .))%>% #perform regression#
  map(summary)%>% 
  map_dbl("r.squared") #calculate R-squared#
PDMS_Rsq_DF <- as.data.frame(PDMS_Rsq) #convert to data frame so it's usable#

PDMS_Rsq_DF <- tibble::rownames_to_column(PDMS_Rsq_DF, "Alignment_ID") #change row headers into a column so it can be merged#

#Add R-squared to regression table#
Reg_PDMS_Rsq <- merge(Reg_PDMS, PDMS_Rsq_DF, by = "Alignment_ID")%>%
  dplyr::rename("Rsquared" = PDMS_Rsq)

#Perform R-squared test#
Reg_PDMS_Rsq[,'test'] = Reg_PDMS_Rsq$statistic*Reg_PDMS_Rsq$Rsquared

#Filter out R-squared less than 0.8#
PDMS_Rsq_filtered <- Reg_PDMS_Rsq%>%
  filter(Rsquared >= 0.8)%>%
  mutate("Sorbent" = "PDMS")

#Reformat for use in plots#
PDMS_stats <- PDMS_Rsq_filtered%>%
  select(c(Alignment_ID, term, estimate, std.error, Rsquared, Sorbent))%>%
  pivot_wider(names_from = term, values_from = estimate)%>%
  dplyr::rename("Uptake_Rate" = Days)%>%
  select(Alignment_ID, Uptake_Rate, std.error, Rsquared, Sorbent)

DF_to_CSV(PDMS_stats) #export regression#

#Repeat for GFF#

#Perform regression#
Reg_GFF <- GFF_Uptake_stats%>%
  group_by(Alignment_ID)%>% #cluster by Alignment_ID#
  do(., tidy(lm(Vol ~ 0+ Days, data = .))) #for each group, calculate regression statistics forcing through the origin#
summary(Reg_GFF)

#Calculate R-squared#
GFF_Rsq <- GFF_Uptake_stats%>%
  split(.$Alignment_ID)%>% #split by Alignment_ID; alternatively, use PDMS_split made previously#
  map(~lm(Vol ~ 0 + Days, data = .))%>% #perform regression#
  map(summary)%>% 
  map_dbl("r.squared") #calculate R-squared#
GFF_Rsq_DF <- as.data.frame(GFF_Rsq) #convert to data frame so it's usable#

GFF_Rsq_DF <- tibble::rownames_to_column(GFF_Rsq_DF, "Alignment_ID") #change row headers into a column so it can be merged#

#Add R-squared to regression table#
Reg_GFF_Rsq <- merge(Reg_GFF, GFF_Rsq_DF, by = "Alignment_ID")%>%
  dplyr::rename("Rsquared" = GFF_Rsq)

#Perform R-squared test#
Reg_GFF_Rsq[,'test'] = Reg_GFF_Rsq$statistic*Reg_GFF_Rsq$Rsquared

GFF_Rsq_filtered <- Reg_GFF_Rsq%>%
  filter(Rsquared >= 0.8)%>%
  mutate("Sorbent" = "GFF")

GFF_stats <- GFF_Rsq_filtered%>%
  select(c(Alignment_ID, term, estimate, std.error, Rsquared, Sorbent))%>%
  pivot_wider(names_from = term, values_from = estimate)%>%
  dplyr::rename("Uptake_Rate" = Days)%>%
  select(Alignment_ID, Uptake_Rate, std.error, Rsquared, Sorbent)

DF_to_CSV(GFF_stats) #export regression#

#Display both sorbent types in the same file#

PDMS_stats_comb <- PDMS_stats%>%
  mutate("Sorbent" = "PDMS") #add sorbent type identifier#

GFF_stats_comb <- GFF_stats%>%
  mutate("Sorbent" = "GFF") #add sorbent type identifier#

#Combine#
Both_PAS_stats <- bind_rows(PDMS_stats_comb, GFF_stats_comb)
Both_PAS_stats$Sorbent <- factor(Both_PAS_stats$Sorbent,
                                 levels = c("PDMS", "GFF"),
                                 labels = c("PDMS", "GFF")) #reorder layers#

DF_to_CSV(Both_PAS_stats) #export regression#

####------Creating uptake curves------####

#Function for creating a scatterplot with regression; df is your data frame,
#x is the x value, y is the y value, x_name and y_name are axis titles.
#Intended to be mapped onto a large data frame that has been split.#
Plot_Maker <- function(df, x, y, x_name, y_name, Max_X = 20,
                       Tick_Spacing = 4){
  df_plot <- df%>%
    ggplot(aes({{x}}, {{y}}))+ #call ggplot2 and assign x and y#
    geom_point(size = 4)+ #use scatterplot#
    theme_bw()+ #aesthetics
    labs(x = x_name, y = y_name)+ #axis titles#
    scale_x_continuous(breaks = seq(0, Max_X, by = Tick_Spacing), expand = c(0, 0, 0.1, 0.1), limits = c(0, Max_X))+
    scale_y_continuous(expand = c(0, 0, 0.1, 0.1), limits = c(0, NA))+
    stat_regline_equation(aes(label = after_stat(eq.label)),
                          formula = {{y}} ~ 0 + {{x}}, hjust = 0, vjust = 0.25)+ #add line equation#
    stat_regline_equation(aes(label = after_stat(rr.label)),
                          formula = {{y}} ~ 0 + {{x}}, hjust = 0, vjust = 1.75)+ #add r-squared#
    geom_smooth(method=lm, formula = {{y}} ~ 0 + {{x}}, se=FALSE)+ #add linear regression#
    theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank())+ #remove gridlines#
    ggtitle(df[1, 1])+ #add title; remove if not using a split data frame or other list#
    theme(axis.title.x = element_text(size = 20), 
          axis.title.y = element_text(size = 20), #adjust text size for axis titles#
          axis.text.x = element_text(colour="black", size = 15),
          axis.text.y = element_text(colour="black", size = 15), #adjust axis text size and colour#
          legend.title = element_text(size = 20), 
          legend.text = element_text(size = 20), #adjust legend title and text#
          plot.title = element_text(size = 15), #adjust plot title text size#
          strip.text = element_text(size = 12), #adjust strip text if faceting#
          axis.ticks = element_line(colour="black", linewidth=1), #adjust tick mark size and colour#
          axis.line = element_line(colour="black", linewidth=1), #adjust axis size and colour#
          panel.border = element_blank(), #remove border#
          strip.background = element_rect(colour="black", fill="white")) #black and white background#
  return(df_plot)
}

plotPDMS_filter <- PDMS_Uptake_stats%>%
  filter(Alignment_ID %in% PDMS_stats$Alignment_ID == TRUE)

#Use Plot_Maker to create plots from data frame#
PlotPDMS_RSD <- plotPDMS_filter%>%
  split(.$Alignment_ID)%>% #split by Alignment_ID#
  map(~Plot_Maker(., Days, Vol, "Days", "Equiv Air Volume")) #map function#

#Export the plots; this exports all plots as separate pages in a single PDF#
#pdf("Figures/All_PDMS_Curves.pdf",
    #width = 7,
    #height = 5) #call function and name the PDF#
#PlotPDMS_RSD #list containing the plots#
#dev.off() #close the PDF

#Repeat for GFF#

#Filter out low R-squared#
plotGFF_filter <- GFF_Uptake_stats%>%
  filter(Alignment_ID %in% GFF_stats$Alignment_ID == TRUE)

#Create the plots#
PlotGFF_RSD <- plotGFF_filter%>%
  split(.$Alignment_ID)%>% #split by Alignment_ID#
  map(~Plot_Maker(., Days, Vol, "Days", "Equiv. Air Volume")) #map function#

#Export the plots; this exports all plots as separate pages in a single PDF#
#pdf("Figures/All_GFF_Curves.pdf",
   # width = 7,
    #height = 5) #call function and name the PDF#
#PlotGFF_RSD #list containing the plots#
#dev.off() #close the PDF

#Plotting GFF and PDMS on the same graph#
PDMSComb <- plotPDMS_filter%>%
  select(c(Alignment_ID, Days, Vol))%>% #from RSD-filtered data, retain only columns needed to graph#
  mutate("PAS_Type" = "PDMS") #add sorbent identifier

#Repeat for GFF#
GFFComb <- plotGFF_filter%>%
  select(c(Alignment_ID, Days, Vol))%>%
  mutate("PAS_Type" = "GFF")

#Combine#
Both_PAS <- bind_rows(PDMSComb, GFFComb)%>%
  arrange(Alignment_ID)
Both_PAS$PAS_Type <- factor(Both_PAS$PAS_Type,
                            levels = c("PDMS", "GFF"),
                            labels = c("PDMS", "GFF")) #reorder layers#

####Box and whisker plot####

#Makes a box plot#
Box_Plot <- function(df, x, y, fill_by, x_title, y_title, max_val_from, showlegend = FALSE,
                     Colour_Palette = "Accent"){
  stats.summ <- function(x){
    upper_limit <- max(max_val_from)
    return(data.frame(y = 1.2 * upper_limit,
                      label = paste("n = ", length(x), "\n",
                                    "mean = ", round(mean(x), 2), "\n",
                                    "median = ", round(median(x), 2), "\n")))
  }
  dfplot <- df%>%
    ggplot(aes(x = {{x}}, y = {{y}}, fill = {{fill_by}}))+
    geom_boxplot(show.legend = showlegend)+ #box plot; legend is unnecessary#
    #geom_dotplot(binaxis='y', stackdir='center', dotsize=1) +
    theme_bw()+ #aesthetics#
    scale_fill_brewer(palette = Colour_Palette) +
    labs(x = x_title, y = y_title)+ #axis titles#
    stat_summary(fun.data = stats.summ, geom = "text", hjust = 0.5, vjust = 0.9, show.legend = FALSE)+ #adds the stats from get.n#
    theme(axis.title.x = element_text(size = 20), 
          axis.title.y = element_text(size = 20), #adjust text size for axis titles#
          axis.text.x = element_text(colour = "black", size = 15),
          axis.text.y = element_text(colour = "black", size = 15), #adjust axis text size and colour#
          plot.title = element_text(size = 15), #adjust plot title text size#
          strip.text = element_text(size = 12), #adjust strip text if faceting#
          axis.ticks = element_line(colour="black", linewidth=1), #adjust tick mark size and colour#
          axis.line = element_line(colour="black", linewidth=1), #adjust axis size and colour#
          panel.border = element_blank(), #remove border#
          panel.grid.major=element_blank(), #remove gridlines#
          panel.grid.minor=element_blank(), #remove gridlines#
          strip.background = element_rect(colour="black", fill="white"))
  
}

#Makes a box plot with facet wrapping#
Box_Plot_Facet <- function(df, x, y, fill_by, x_title, y_title, max_val_from, facet,
                           Colour_Palette = "Accent"){
  stats.summ <- function(x){
    upper_limit <- max(max_val_from)
    return(data.frame(y = 1.2 * upper_limit,
                      label = paste("n = ", length(x), "\n",
                                    "mean = ", round(mean(x), 2), "\n",
                                    "median = ", round(median(x), 2), "\n")))
  }
  dfplot <- df%>%
    ggplot(aes(x = {{x}}, y = {{y}}, fill = {{fill_by}}))+
    geom_boxplot(show.legend = FALSE)+ #box plot; legend is unnecessary#
    theme_bw()+ #aesthetics#
    labs(x = x_title, y = y_title)+ #axis titles#
    scale_fill_brewer(palette = Colour_Palette) +
    facet_wrap(enquo(facet))+
    stat_summary(geom = "text", fun.data = stats.summ, size = 500, hjust = 0.5, vjust = 0.9, show.legend = FALSE)+ #adds the stats from get.n#
    theme(axis.title.x = element_text(size = 20), 
          axis.title.y = element_text(size = 20), #adjust text size for axis titles#
          axis.text.x = element_text(colour = "black", size = 15),
          axis.text.y = element_text(colour = "black", size = 15), #adjust axis text size and colour#
          plot.title = element_text(size = 15), #adjust plot title text size#
          strip.text = element_text(size = 12), #adjust strip text if faceting#
          axis.ticks = element_line(colour="black", linewidth=1), #adjust tick mark size and colour#
          axis.line = element_line(colour="black", linewidth=1), #adjust axis size and colour#
          panel.border = element_blank(), #remove border#
          panel.grid.major=element_blank(), #remove gridlines#
          panel.grid.minor=element_blank(), #remove gridlines#
          strip.background = element_rect(colour="black", fill="white"))
  
}

PDMS_Curve_List <- 23719 #this feature had a non-linear uptake curve after filtering

PDMS_UR_fullfilt <- PDMS_stats%>%
  filter(Alignment_ID %in% PDMS_Curve_List == FALSE) %>% #remove non-linear feature
  mutate(Uptake_Rate = round(Uptake_Rate, digits = 2),
         Rsquared = round(Rsquared, digits = 2),
         std.error = round(std.error, digits = 2),
         Mean = round(mean(Uptake_Rate), digits = 2),
         STDEV = round(sd(Uptake_Rate), digits = 2),
         GeoMean = round(exp(mean(log(Uptake_Rate))), digits = 2),
         GeoSTDEV = round(exp(sd(log(Uptake_Rate))), digits = 2)) #formatting for easier use in tables

GFF_UR_fullfilt <- GFF_stats%>%
  mutate(Uptake_Rate = round(Uptake_Rate, digits = 2),
         Rsquared = round(Rsquared, digits = 2),
         std.error = round(std.error, digits = 2),
         Mean = round(mean(Uptake_Rate), digits = 2),
         STDEV = round(sd(Uptake_Rate), digits = 2),
         GeoMean = round(exp(mean(log(Uptake_Rate))), digits = 2),
         GeoSTDEV = round(exp(sd(log(Uptake_Rate))), digits = 2)) #formatting for easier use in tables

#Remake uptake curve files post filtering#

plotprepPDMS_UR_fullfilt <- PDMS_Uptake_stats%>%
  filter(Alignment_ID %in% PDMS_UR_fullfilt$Alignment_ID == TRUE)

#Use Plot_Maker to create plots from newly filtered data#
plotPDMS_UR_fullfilt <- plotprepPDMS_UR_fullfilt%>%
  split(.$Alignment_ID)%>% #split by INCHIKEY#
  map(~Plot_Maker(., Days, Vol, "Days", "Equiv Air Volume")) #map function#

#Export the plots; this exports all plots as separate pages in a single PDF#
#pdf("Figures/All_PDMS_RSD_Curves_fullfilt.pdf",
    #width = 7,
    #height = 5) #call function and name the PDF#
#plotPDMS_UR_fullfilt #list containing the plots#
#dev.off() 

#Repeat for GFF#
plotprepGFF_UR_fullfilt <- GFF_Uptake_stats%>%
  filter(Alignment_ID %in% GFF_UR_fullfilt$Alignment_ID == TRUE)

#Use Plot_Maker to create plots from newly filtered data#
plotGFF_UR_fullfilt <- plotprepGFF_UR_fullfilt%>%
  split(.$Alignment_ID)%>% #split by INCHIKEY#
  map(~Plot_Maker(., Days, Vol, "Days", "Equiv Air Volume")) #map function#

#Export the plots; this exports all plots as separate pages in a single PDF#
#pdf("Figures/All_GFF_RSD_Curves_fullfilt.pdf",
    #width = 7,
    #height = 5) #call function and name the PDF#
#plotGFF_UR_fullfilt #list containing the plots#
#dev.off()

#Export tables of the uptake rates#
DF_to_CSV(PDMS_UR_fullfilt)
DF_to_CSV(GFF_UR_fullfilt)

#Combine the filtered data#
NewFilt_Both_PAS <- bind_rows(PDMS_UR_fullfilt, GFF_UR_fullfilt)
NewFilt_Both_PAS$Sorbent <- factor(NewFilt_Both_PAS$Sorbent,
                                   levels = c("PDMS", "GFF"),
                                   labels = c("PDMS", "GFF")) #reorder layers#

#Export the combined table#
DF_to_CSV(NewFilt_Both_PAS)

####Comparing uptake rates####

NewFilt_BP_Rates <- Box_Plot(NewFilt_Both_PAS, Sorbent, Uptake_Rate, Sorbent, element_blank(),
                             expression(paste("Uptake Rate (m"^"3"*" day"^"-1"*" dm"^"-2"*")")),
                             NewFilt_Both_PAS$Uptake_Rate)
NewFilt_BP_Rates

save_plot(plot = NewFilt_BP_Rates, 'Figures/NewFilt_BP_Rates.svg', base_width = 8, base_height = 7)

#Remove outliers by IQR#

IQR_PDMS <- PDMS_UR_fullfilt$Uptake_Rate #isolate PDMS uptake rates#
quart_PDMS <- quantile(IQR_PDMS) #determine quartiles#
IQR1.5_PDMS <- 1.5*IQR(IQR_PDMS) #create variable for 1.5*IQR#

PDMS_ufence <- as.numeric(quart_PDMS[4]) + IQR1.5_PDMS #determine upper limit from Q3 + 1.5IQR#
PDMS_lfence <- as.numeric(quart_PDMS[2]) - IQR1.5_PDMS #determine lower limit from Q1 - 1.5IQR#

#Repeat for GFF#
IQR_GFF <- GFF_UR_fullfilt$Uptake_Rate
quart_GFF <- quantile(IQR_GFF)
IQR1.5_GFF <- 1.5*IQR(IQR_GFF)

GFF_ufence <- as.numeric(quart_GFF[4]) + IQR1.5_GFF
GFF_lfence <- as.numeric(quart_GFF[2]) - IQR1.5_GFF

#Remove outliers from PDMS data#
PDMS_UR_fullfilt_IQR <- PDMS_UR_fullfilt%>%
  filter(Uptake_Rate < PDMS_ufence & Uptake_Rate > PDMS_lfence)%>% #filters outliers#
  mutate(Mean = round(mean(Uptake_Rate), digits = 2),
         STDEV = round(sd(Uptake_Rate), digits = 2),
         GeoMean = round(exp(mean(log(Uptake_Rate))), digits = 2),
         GeoSTDEV = round(exp(sd(log(Uptake_Rate))), digits = 2)) %>%#redo mean and sd calculations#
  mutate(RSD = 100*STDEV/Mean)

#Repeat for GFF#
GFF_UR_fullfilt_IQR <- GFF_UR_fullfilt%>%
  filter(Uptake_Rate < GFF_ufence & Uptake_Rate > GFF_lfence)%>% #filters outliers#
  mutate(Mean = round(mean(Uptake_Rate), digits = 2),
         STDEV = round(sd(Uptake_Rate), digits = 2),
         GeoMean = round(exp(mean(log(Uptake_Rate))), digits = 2),
         GeoSTDEV = round(exp(sd(log(Uptake_Rate))), digits = 2)) %>%#redo mean and sd calculations#
  mutate(RSD = 100*STDEV/Mean)

#Combine#
NewFilt_Both_PAS_IQR <- bind_rows(PDMS_UR_fullfilt_IQR, GFF_UR_fullfilt_IQR)
NewFilt_Both_PAS_IQR$Sorbent <- factor(NewFilt_Both_PAS_IQR$Sorbent,
                                       levels = c("PDMS", "GFF"),
                                       labels = c("PDMS", "GFF")) #reorder layers#

#Plot IQR corrections#
NewFilt_IQR_BP_Rates <- Box_Plot(NewFilt_Both_PAS_IQR, Sorbent, Uptake_Rate, Sorbent, "Sorbent Type",
                                 expression(paste("Uptake Rate (m"^"3"*" day"^"-1"*" dm"^"-2"*")")),
                                 NewFilt_Both_PAS_IQR$Uptake_Rate)
NewFilt_IQR_BP_Rates

save_plot(plot = NewFilt_IQR_BP_Rates, 'Figures/NewFilt_IQR_BP_Rates.svg', base_width = 8, base_height = 7)

####Checking for significance####

#Before IQR#

#Check for normality#
PDMS_shapiro <- shapiro.test(PDMS_UR_fullfilt$Uptake_Rate)
GFF_shapiro <- shapiro.test(GFF_UR_fullfilt$Uptake_Rate)

Shapiro_Test = c(PDMS_shapiro, GFF_shapiro)

for(i in Shapiro_Test[c(2,6)]){
  print(i)
  if (i <= 0.05) {
    print("Normally distributed")
  } else {
    print("Not normally distributed")
  }
}

#Check for significance between the sorbent types without normality#
Sorbents_Wilcox <- wilcox.test(PDMS_UR_fullfilt$Uptake_Rate, GFF_UR_fullfilt$Uptake_Rate, alternative = "two.sided")
Sorbents_Wilcox$p.value

#Check for significance between the sorbent types with normality#
Sorbents_ttest <- t.test(PDMS_UR_fullfilt$Uptake_Rate, GFF_UR_fullfilt$Uptake_Rate, alternative = "two.sided", var.equal = FALSE)
Sorbents_ttest$p.value

#After IQR#

#Check for normality#
PDMS_shapiro_IQR <- shapiro.test(PDMS_UR_fullfilt_IQR$Uptake_Rate) #p < 0.05; not normally distributed#
GFF_shapiro_IQR <- shapiro.test(GFF_UR_fullfilt_IQR$Uptake_Rate) #p > 0.05; normally distributed#

Shapiro_Test_IQR = c(PDMS_shapiro_IQR, GFF_shapiro_IQR)

for(i in Shapiro_Test_IQR[c(2,6)]){
  print(i)
  if (i < 0.05) {
    print("Not normally distributed")
  } else {
    print("Normally distributed")
  }
}

#Check for significance between the sorbent types without normality#
Sorbents_Wilcox_IQR <- wilcox.test(PDMS_UR_fullfilt_IQR$Uptake_Rate, GFF_UR_fullfilt_IQR$Uptake_Rate, alternative = "two.sided")
Sorbents_Wilcox_IQR$p.value

#Check for significance between the sorbent types with normality#
Sorbents_ttest_IQR <- t.test(PDMS_UR_fullfilt_IQR$Uptake_Rate, GFF_UR_fullfilt_IQR$Uptake_Rate, alternative = "two.sided", var.equal = FALSE)
Sorbents_ttest_IQR$p.value

####-----list of features-----####

#Long bound table#
Separate_dedup_individual_long <- bind_rows(Dedup_A_long, Dedup_D_long, Dedup_G_long)

All_Features <- Seperate_dedup_individ_bound %>%
  subset(select = c(Alignment_ID, Metabolite_name, INCHIKEY, SMILES, Formula, Sampler, level, Average_Mz, m.z_similarity)) %>%
  mutate(Found_On = ifelse(Alignment_ID %in% Dedup_A_long$Alignment_ID == TRUE
                           & Alignment_ID %in% Dedup_D_long$Alignment_ID == TRUE
                           & Alignment_ID %in% Dedup_G_long$Alignment_ID == TRUE,
                           "All",
                           ifelse(Alignment_ID %in% Dedup_A_long$Alignment_ID == FALSE
                                  & Alignment_ID %in% Dedup_D_long$Alignment_ID == TRUE
                                  & Alignment_ID %in% Dedup_G_long$Alignment_ID == TRUE,
                                  "PDMS and GFF",
                                  ifelse(Alignment_ID %in% Dedup_A_long$Alignment_ID == TRUE
                                         & Alignment_ID %in% Dedup_D_long$Alignment_ID == FALSE
                                         & Alignment_ID %in% Dedup_G_long$Alignment_ID == TRUE,
                                         "Active and GFF",
                                         ifelse(Alignment_ID %in% Dedup_A_long$Alignment_ID == TRUE
                                                & Alignment_ID %in% Dedup_D_long$Alignment_ID == TRUE
                                                & Alignment_ID %in% Dedup_G_long$Alignment_ID == FALSE,
                                                "Active and PDMS",
                                                ifelse(Alignment_ID %in% Dedup_A_long$Alignment_ID == TRUE
                                                       & Alignment_ID %in% Dedup_D_long$Alignment_ID == FALSE
                                                       & Alignment_ID %in% Dedup_G_long$Alignment_ID == FALSE,
                                                       "Active Only",
                                                       ifelse(Alignment_ID %in% Dedup_A_long$Alignment_ID == FALSE
                                                              & Alignment_ID %in% Dedup_D_long$Alignment_ID == TRUE
                                                              & Alignment_ID %in% Dedup_G_long$Alignment_ID == FALSE,
                                                              "PDMS Only",
                                                              "GFF Only"))))))) #identify what samplers they are found on

#-------OPERA export------- 

#import list of SMILES known to cause error#
bad_smiles1 <- read.delim(file = "OPERA 09/New_Opera_Fails.txt") %>%
  mutate(Alignment_ID = as.character(Alignment_ID))
bad_smiles2 <- read.delim(file = "OPERA 09/OPERA_newest_fails.txt") %>%
  mutate(Alignment_ID = as.character(Alignment_ID))

bad_smiles <- bind_rows(bad_smiles1, bad_smiles2)

opera_export <- Seperate_dedup_individ_bound %>%
  filter(SMILES != 'null' & !is.na(SMILES)) %>%#filter for features with SMILES since OPERA requires tis
  filter(SMILES %in% bad_smiles$SMILES == FALSE) #filter out error SMILES#

smi_data <- convert_to_smi_tab_delimited(opera_export, label = 'Alignment_ID') #using alignment ID in case something has a smiles but not alignment ID

#Write the tab-delimited .smi file
writeLines(smi_data, "OPERA 09/OPERA_SMILES.smi")

#-------import OPERA-------

OPERA_all <- read_csv("OPERA 09/09 OPERA Results/OPERA_All.csv", col_names = T) %>%
  dplyr::rename(Alignment_ID = MoleculeID) %>% #id was 2D inchikey
  merge(., opera_export, by = c('Alignment_ID')) %>% #merge back with rest of information for each feature
  subset(select = c(`2D_INCHIKEY`, LogKOA_pred, Alignment_ID)) %>% #get only infromation needed
  group_by(Alignment_ID) %>% #prevent duplicates when later merging (use alignment id in case some don't have 2d inchikey)
  slice(1) %>% #prevent duplicates when later merging
  ungroup()

#see if there's any molecules missing (except those OPERA failed for)
missing_opera <- merge(OPERA_all, opera_export, by = c('Alignment_ID', '2D_INCHIKEY'), all = T) %>% #had failed smiles
  filter(is.na(LogKOA_pred))
smi_data2 <- convert_to_smi_tab_delimited(missing_opera, label = 'Alignment_ID') #export molecules missing OPERA predictions
writeLines(smi_data2, "OPERA_SMILES_Additional.smi")

#merge back with things detected in different sampler types
OPERA_A <- Dedup_A %>%
  filter(SMILES != 'null'& !is.na(SMILES)) %>%
  filter(level < 4) %>% #KOA only for level 2 and 3
  merge(., OPERA_all, by = c('Alignment_ID', '2D_INCHIKEY')) #merge with other information

OPERA_D <- Dedup_D %>%
  filter(SMILES != 'null'& !is.na(SMILES)) %>%
  filter(level < 4) %>% #KOA only for level 2 and 3
  merge(., OPERA_all, by = c('Alignment_ID', '2D_INCHIKEY')) #merge with other information

OPERA_G <- Dedup_G %>%
  filter(SMILES != 'null'& !is.na(SMILES)) %>%
  filter(level < 4) %>% #KOA only for level 2 and 3
  merge(., OPERA_all, by = c('Alignment_ID', '2D_INCHIKEY')) #merge with other information

all_koas_filtered <- bind_rows(OPERA_A, OPERA_D, OPERA_G) #bind back together with duplicates across the ther types remaining

OPERA_Sampler_n <- all_koas_filtered %>% 
  group_by(Alignment_ID) %>%
  filter(n()==1)%>% #filter for features uniquely found in each lens
  group_by(Sampler) %>% 
  dplyr::summarise(n= n(), min = min(LogKOA_pred), max = max(LogKOA_pred)) %>%
  ungroup()

Level2KOAs <- all_koas_filtered %>%
  filter(level == 2)

All_Features_KOAadded <- all_koas_filtered %>%
  select(c(Alignment_ID, Sampler, LogKOA_pred)) %>%
  merge(., All_Features, by = c("Alignment_ID", "Sampler"), all.y = TRUE)

DF_to_CSV(Level2KOAs)

####Denisty plots####

#Makes density plots#
Density_Plot <- function(df, x, fill.by, x_name, y_name, legend_title, legend_pos = "left",
                         colour_palette = Sampler_Colours, low_x = 0, high_x = 14, x_tick = 1){
  dplot <- df%>%
    ggplot(aes(x = {{x}}, fill = {{fill.by}}))+
    geom_density(alpha = 0.5)+
    labs(x = x_name, y = y_name, fill = legend_title)+
    theme_bw()+
    scale_fill_manual(values = colour_palette) +
    scale_y_continuous(labels = function(x) format(x, nsmall = 3), expand = c(0, 0))+
    scale_x_continuous(breaks = seq(low_x, high_x, by = x_tick), expand = c(0, 0), limits = c(low_x, high_x))+
    theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank())+
    theme(axis.title.x = element_text(size = 20), 
          axis.title.y = element_text(size = 20), #adjust text size for axis titles#
          axis.text.x = element_text(colour="black", size = 15),
          axis.text.y = element_text(colour="black", size = 15), #adjust axis text size and colour#
          legend.title = element_text(size = 20), 
          legend.text = element_text(size = 20), #adjust legend title and text#
          legend.position = legend_pos,
          axis.ticks = element_line(colour="black", linewidth=1), #adjust tick mark size and colour#
          axis.line = element_line(colour="black", linewidth=1), #adjust axis size and colour#
          panel.border = element_blank(), #remove border#
          strip.background = element_rect(colour="black", fill="white")) #black and white background#
  return(dplot)
}

#Makes density plots#
Density_Plot_Facet <- function(df, x, x_name, y_name,
                               facet_by, outline_colour,
                               fill_colour){
  dplot <- df%>%
    ggplot(aes(x = {{x}}))+
    geom_density(alpha = 0.5, colour = outline_colour, fill = fill_colour)+
    labs(x = x_name, y = y_name)+
    theme_bw()+
    facet_wrap(enquo(facet_by))+
    scale_y_continuous(expand = c(0, 0))+
    scale_x_continuous(expand = c(0, 0))+
    theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank())+
    theme(axis.title.x = element_text(size = 20), 
          axis.title.y = element_text(size = 20), #adjust text size for axis titles#
          axis.text.x = element_text(colour="black", size = 15),
          axis.text.y = element_text(colour="black", size = 15), #adjust axis text size and colour#
          legend.title = element_text(size = 20), 
          legend.text = element_text(size = 20), #adjust legend title and text#
          axis.ticks = element_line(colour="black", linewidth=1), #adjust tick mark size and colour#
          axis.line = element_line(colour="black", linewidth=1), #adjust axis size and colour#
          panel.border = element_blank(), #remove border#
          strip.background = element_rect(colour="black", fill="white")) #black and white background#
  return(dplot)
}

All_Features_KOAadded$Sampler <- factor(All_Features_KOAadded$Sampler,
                                    levels = c("PDMS", "GFF", "Active"),
                                    labels = c("PDMS", "GFF", "Active")) #reorder layers#

KOAs_for_Plots <- all_koas_filtered %>%
  select(c(Alignment_ID, LogKOA_pred, Average_Mz, Sampler))

Detection_Comparison <- All_Features_KOAadded %>%
  filter(!is.na(LogKOA_pred) & Sampler != "Active") %>% #don't need active
  group_by(Sampler) %>%
  dplyr::mutate(n = n()) %>%
  ungroup() %>%
  Density_Plot(., LogKOA_pred, Sampler,
               expression(paste("Predicted log"[10]*" K"[OA])),
               "Density", element_blank(), "bottom", Sampler_Colours)
Detection_Comparison

save_plot(plot = Detection_Comparison, 'Figures/Detection_Comparison.svg', base_width = 8, base_height = 7)

Detection_Comparison_MZ <- All_Features_KOAadded %>%
  filter(Sampler != "Active") %>%
  group_by(Sampler) %>%
  dplyr::mutate(n = n()) %>%
  ungroup() %>%
  Density_Plot(., Average_Mz, Sampler,
               "m/z", "Density", element_blank(), "bottom",
               low_x = 0, high_x = 1750, x_tick = 150)
Detection_Comparison_MZ

save_plot(plot = Detection_Comparison_MZ, 'Figures/Detection_Comparison_MZ.svg', base_width = 8, base_height = 7)

All_Features_KOAadded[, "Has_Good_UR"] = ifelse(All_Features_KOAadded$Alignment_ID %in% PDMS_UR_fullfilt$Alignment_ID == TRUE & All_Features_KOAadded$Sampler == "PDMS",
                                                "Rate-derived PDMS",
                                                ifelse(All_Features_KOAadded$Alignment_ID %in% GFF_UR_fullfilt$Alignment_ID == TRUE & All_Features_KOAadded$Sampler == "GFF",
                                                       "Rate-derived GFF",
                                                       "No rate derived"))

All_Features_KOAadded$Has_Good_UR <- factor(All_Features_KOAadded$Has_Good_UR,
                                        levels = c("Rate-derived PDMS", "Rate-derived GFF", "No rate derived"),
                                        labels = c("Rate-derived PDMS", "Rate-derived GFF", "No rate derived")) #reorder layers#

####-----Comparing uptake rates across studies-----####
StudyURs_imp <- read.csv("C://Users//scrol//OneDrive - McMaster University//Joseph Okeme's files - Alyss//New passives calibration//Manuscript for PAS//PaperURs.csv")%>%
  pivot_longer(cols = 1:11, names_to = "Study", values_to = "Uptake_Rate")%>%
  filter(Uptake_Rate != 0)%>%
  mutate(Sampler = ifelse(grepl("_up", Study),
                          "Unsheltered PUF",
                          ifelse(grepl("_db", Study),
                                 "Double bowl PUF",
                                 ifelse(grepl("_sb", Study) | grepl("PUF", Study),
                                        "Single bowl PUF",
                                        ifelse(grepl("_PDMS", Study),
                                               "Unsheltered PDMS",
                                               "XAD-Mesh-PAS")))))%>%
  mutate(Study = gsub("_", " ", Study))%>%
  mutate(Study = gsub("al", "al.", Study))%>%
  mutate(Study = gsub(" up", "", Study))%>%
  mutate(Study = gsub(" db", "", Study))%>%
  mutate(Study = gsub(" sb", "", Study))%>%
  mutate(Study = gsub(" PDMS", "", Study))%>%
  mutate(Study = gsub(" XAD", "", Study))%>%
  mutate(Study = gsub(" PUF", "", Study))%>% #reformatting
  group_by(Study, Sampler)%>%
  dplyr::mutate(n = n())%>% #number of uptake rates
  ungroup()

OurStudyURs <- NewFilt_Both_PAS_IQR %>%
  dplyr::rename("Sampler" = Sorbent)%>%
  mutate(Sampler = ifelse(Sampler == "PDMS",
                          "PDMS-Mesh-PAS",
                          "GFF-Mesh-PAS"))%>%
  mutate(Study = "Present Work")%>%
  select(c(Study, Uptake_Rate, Sampler)) %>%
  group_by(Study, Sampler)%>% #reformatting
  dplyr::mutate(n = n())%>% #number of uptake rates
  ungroup()

Comparing_Study_URs <- bind_rows(StudyURs_imp, OurStudyURs) %>%
  arrange(Uptake_Rate) #combine

Comparing_Study_URs$Study = factor(Comparing_Study_URs$Study,
                                   levels = c("Present Work", Comparing_Study_URs$Study),
                                   labels = c("Present Work", Comparing_Study_URs$Study)) #reorder layers#

Comparing_Study_URs$Sampler = factor(Comparing_Study_URs$Sampler,
                                     levels = c("PDMS-Mesh-PAS", "GFF-Mesh-PAS",
                                                "Unsheltered PDMS", "Double bowl PUF",
                                                "Single bowl PUF", "Unsheltered PUF",
                                                "XAD-Mesh-PAS"),
                                     labels = c("PDMS-Mesh-PAS", "GFF-Mesh-PAS",
                                                "Unsheltered PDMS", "Double bowl PUF",
                                                "Single bowl PUF", "Unsheltered PUF",
                                                "XAD-Mesh-PAS")) #reorder layers#

Comparing_Study_URs <- Comparing_Study_URs %>%
  group_by(Study, Sampler) %>%
  dplyr::mutate(Mean = mean(Uptake_Rate)) %>%
  dplyr::mutate(SD = sd(Uptake_Rate)) %>%
  dplyr::mutate(RSD = 100*SD/Mean) %>%
  dplyr::mutate(ID = paste(Study, Sampler)) %>% #calculate stats
  ungroup()

BP_Comparing_Study_URs <- Comparing_Study_URs %>%
  ggplot(aes(x = Study, y = Uptake_Rate, fill = Sampler))+
  geom_boxplot()+ #box plot; legend is unnecessary#
  scale_fill_brewer(palette = "Accent")+
  theme_bw()+ #aesthetics#
  labs(x = element_blank(), y = expression(paste("Uptake Rate (m"^"3"*" day"^"-1"*" dm"^"-2"*")")))+ #axis titles#
  scale_x_discrete(labels = function(x){str_wrap(x, 10)}) +  
  theme(axis.title.x = element_text(size = 20), 
        axis.title.y = element_text(size = 20), #adjust text size for axis titles#
        axis.text.x = element_text(colour = "black", size = 15),
        axis.text.y = element_text(colour = "black", size = 15), #adjust axis text size and colour#
        plot.title = element_text(size = 15), #adjust plot title text size#
        strip.text = element_text(size = 12), #adjust strip text if faceting#
        axis.ticks = element_line(colour="black", linewidth=1), #adjust tick mark size and colour#
        axis.line = element_line(colour="black", linewidth=1), #adjust axis size and colour#
        panel.border = element_blank(), #remove border#
        panel.grid.major=element_blank(), #remove gridlines#
        panel.grid.minor=element_blank(), #remove gridlines#
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 15),
        legend.background = element_rect(colour = "black", linewidth = 1),
        #legend.position = "bottom",
        strip.background = element_rect(colour="black", fill="white"))
BP_Comparing_Study_URs

save_plot(plot = BP_Comparing_Study_URs, 'Figures/Figure_3.svg', base_width = 15, base_height = 8)

#-------visualize classification-------

#1. filter for level 2 and 3 compounds
Seperate_dedup_individ_bound_level32 <- Seperate_dedup_individ_bound %>% 
  filter(level == 2 | level == 3) %>%
  mutate(ClassIDing = paste(Alignment_ID, " ", Formula))

#2. filter for confident classifications
classes1 <- classes %>%
  #filter(mappingFeatureId %in% Seperate_dedup_individ_bound_level32$Alignment_ID) %>% #filter for features surviving filtering (using Alignment ID original so it's the number, and no inchikeys to match mapping feature ID)
  #filter(molecularFormula %in% Seperate_dedup_individ_bound_level32$Formula) %>%
  mutate(ClassIDing = paste(mappingFeatureId, " ", molecularFormula)) %>%
  filter(ClassIDing %in% Seperate_dedup_individ_bound_level32$ClassIDing) %>%
  mutate(`ClassyFire#most specific class Probability` = as.numeric(`ClassyFire#most specific class Probability`)) %>%
  mutate(`ClassyFire#superclass probability` = as.numeric(`ClassyFire#superclass probability`)) %>%
  mutate(`ClassyFire#class Probability` = as.numeric(`ClassyFire#class Probability`)) %>%
  mutate(`ClassyFire#subclass Probability` = as.numeric(`ClassyFire#subclass Probability`)) %>%
  mutate(`ClassyFire#level 5 Probability` = as.numeric(`ClassyFire#level 5 Probability`)) %>%
  pivot_longer(cols = c(`ClassyFire#superclass probability`, `ClassyFire#class Probability`, `ClassyFire#subclass Probability`,
                        `ClassyFire#level 5 Probability`)) %>% 
  mutate(class_level = ifelse(name == 'ClassyFire#superclass probability', 1, NA)) %>% #these levels will be used later when keeping everything before a certain classification
  mutate(class_level = ifelse(name == 'ClassyFire#class Probability', 2, class_level)) %>%
  mutate(class_level = ifelse(name == 'ClassyFire#subclass Probability', 3, class_level)) %>%
  mutate(class_level = ifelse(name == 'ClassyFire#level 5 Probability', 4, class_level)) %>%
  filter(value >= 0.6) %>% #filter for only high confidence classifciations
  group_by(mappingFeatureId) %>% 
  dplyr::mutate(max_level = max(class_level)) %>% #this is to be able to filter for everything that comes before and including the most specific classification at adequate confidence
  ungroup() %>%
  filter(class_level < max_level | class_level == max_level) %>% #filter for everythig more broad than the most specific confident prediction
  pivot_wider(names_from = name, values_from = value) %>%
  dplyr::rename(Superclass = `ClassyFire#superclass`) %>% 
  dplyr::rename(Class = `ClassyFire#class`) %>%
  dplyr::rename(Subclass = `ClassyFire#subclass`) %>%
  dplyr::rename(`Level 5` =  `ClassyFire#level 5`) %>%
  dplyr::rename(Alignment_ID =  mappingFeatureId) %>%
  filter(!is.na(Superclass)) #these should not appear in plot

#3. get number of each classification remaining after filtering to put in manuscript text
unique_superclass <- unique(classes1$Superclass)
unique_class <- unique(classes1$Class)
unique_subclass <- unique(classes1$Subclass)
unique_level5 <- unique(classes1$`Level 5`)
number_compounds_classified <- unique(classes1$Alignment_ID)
lost_compounds_duuring_classification <- length(number_compounds_classified)-length(unique(Seperate_dedup_individ_bound_level32$Alignment_ID))

#Plots#

#List of classifications for each sampler#
PDMS_AllClasses <- classes1 %>%
  filter(Alignment_ID %in% Dedup_D$Alignment_ID) %>%
  mutate(Sampler = "PDMS")

GFF_AllClasses <- classes1 %>%
  filter(Alignment_ID %in% Dedup_G$Alignment_ID) %>%
  mutate(Sampler = "GFF")

Active_AllClasses <- classes1 %>%
  filter(Alignment_ID %in% Dedup_A$Alignment_ID) %>%
  mutate(Sampler = "Active")

#Recombined#

AllClasses <- bind_rows(PDMS_AllClasses, GFF_AllClasses, Active_AllClasses)

#Just superclasses# 

AllSuperclasses <- AllClasses %>%
  filter(class_level == 1)

#Compare peak area from each superclass#

#Create data frames with peak areas for each sorbent#
PDMS_PeakAreas <- uPDMS %>%
  filter(Batch == 5) %>%
  group_by(Alignment_ID) %>%
  dplyr::mutate(Avg_Area = mean(Area)) %>%
  ungroup() %>%
  mutate(Area = Area/RSA_PDMS) %>% #normalize to surface area#
  pivot_wider(names_from = Sample, values_from = Area) %>%
  mutate(Sampler = "PDMS") %>%
  select(c(Alignment_ID, Avg_Area, Sampler, Metabolite_name))

#Add classes#
PDMS_PeakAreas_Classes <- merge(PDMS_AllClasses, PDMS_PeakAreas, by = c("Alignment_ID", "Sampler")) %>%
  filter(class_level == 1) %>%
  group_by(Superclass) %>%
  dplyr::mutate(Total_Area = sum(Avg_Area),
                N_Superclass = n()) %>%
  ungroup() %>%
  group_by(Sampler) %>%
  dplyr::mutate(N_Sorbent = n()) %>%
  ungroup()

#Repeat for GFF#
GFF_PeakAreas <- uGFF %>%
  filter(Batch == 5) %>%
  #filter(Metabolite_name != "Erucamide" & Metabolite_name != "Oleamide") %>%
  group_by(Alignment_ID) %>%
  dplyr::mutate(Avg_Area = mean(Area)) %>%
  ungroup() %>%
  mutate(Area = Area/RSA_GFF) %>% #normalize to surface area#
  pivot_wider(names_from = Sample, values_from = Area) %>%
  mutate(Sampler = "GFF") %>%
  select(c(Alignment_ID, Avg_Area, Sampler, Metabolite_name))

GFF_PeakAreas_Classes <- merge(GFF_AllClasses, GFF_PeakAreas, by = c("Alignment_ID", "Sampler")) %>%
  filter(class_level == 1) %>%  
  group_by(Superclass) %>%
  dplyr::mutate(Total_Area = sum(Avg_Area),
                N_Superclass = n()) %>%
  ungroup() %>%
  group_by(Sampler) %>%
  dplyr::mutate(N_Sorbent = n()) %>%
  ungroup()

#Combine#
PAS_PeakAreas <- bind_rows(PDMS_PeakAreas_Classes, GFF_PeakAreas_Classes) %>%
  dplyr::rename(Sorbent = Sampler)

PAS_PeakAreas$Sorbent <- factor(PAS_PeakAreas$Sorbent,
                                levels = c("PDMS", "GFF"),
                                labels = c("PDMS", "GFF")) #reorder layers#

Superclass_With_PeakArea_Both <- PAS_PeakAreas%>%
  filter(N_Superclass > 10) %>%
  ggplot(aes(x = reorder(str_wrap(Superclass, 14), N_Superclass), fill = Total_Area))+
  geom_hline(aes(yintercept = y), data.frame(y = c(0:5) * 100), color = "black") + 
  geom_bar(position = "dodge", show.legend = TRUE)+
  # Lollipop shaft for mean gain per region
  geom_segment(aes(x = reorder(str_wrap(Superclass, 14), Total_Area), y = 0, 
                   xend = reorder(str_wrap(Superclass, 14), Total_Area), 
                   yend = 500), linetype = "dashed", color = "black") +
  #scale_y_log10()+
  scale_y_continuous(limits = c(-50, 550), expand = c(0, 0))+
  ggplot2::annotate("text", x = 0.45, y = 120, label = "100", size = 4.5) +
  ggplot2::annotate("text", x = 0.45, y = 220, label = "200", size = 4.5) +
  ggplot2::annotate("text", x = 0.45, y = 320, label = "300", size = 4.5) +
  ggplot2::annotate("text", x = 0.45, y = 420, label = "400", size = 4.5) +
  ggplot2::annotate("text", x = 0.45, y = 520, label = "500", size = 4.5) +
  theme_bw()+
  scale_fill_distiller(palette = "YlOrBr", direction = 1)+
  coord_radial(start = -0.04, expand = TRUE, clip = "off")+
  facet_wrap(~Sorbent, ncol = 1, nrow = 2)+
  labs(x = element_blank(), fill  = str_wrap("Total Normalized Peak Area (PAU)", 20))+ #axis titles#
  theme(#axis.title.x = element_text(size = 20), 
    #axis.title.y = element_text(size = 20), #adjust text size for axis titles#
    axis.text.x = element_text(colour = "black", size = 15),
    #axis.text.y = element_text(colour = "black", size = 7), #adjust axis text size and colour#
    legend.title = element_text(size = 15),
    #legend.margin = margin(2,2,2,2, "cm"),
    legend.background = element_rect(colour = "black", linewidth = 1),
    legend.frame = element_rect(colour = "black", linewidth = 1),
    legend.ticks = element_line(colour = "black", linewidth = 1),
    axis.text.y = element_blank(),
    legend.text = element_text(size = 15), #adjust legend title and text#
    legend.key.width = unit(2, "cm"),
    legend.key.height = unit(1, "cm"),
    legend.position = "bottom",
    plot.title = element_text(size = 15), #adjust plot title text size#
    plot.background = element_rect(fill = "white"),
    strip.text = element_text(size = 15), #adjust strip text if faceting#
    #axis.ticks = element_line(colour="black", linewidth = 1), #adjust tick mark size and colour#
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    #axis.line = element_line(colour = "black", linewidth = 1), #adjust axis size and colour#
    axis.line = element_blank(),
    panel.border = element_blank(), #remove border#
    #panel.spacing.x = unit(8, "cm"),
    panel.grid.major=element_blank(), #remove gridlines#
    panel.grid.minor=element_blank(), #remove gridlines#
    strip.background = element_rect(colour = "black", fill = "white", linewidth = 1))
Superclass_With_PeakArea_Both

Superclass_and_PChem <- merge(AllSuperclasses, All_Features_KOAadded, by = c("Alignment_ID", "Sampler"), all.x = TRUE) %>%
  select(c(Alignment_ID, Superclass, Class, Subclass, Sampler, LogKOA_pred, Average_Mz)) #table containing classes and physchem properties

#Add classes to uptake rates#
URs_SuperclassesPchem <- NewFilt_Both_PAS %>%
  mutate(Sampler = Sorbent) %>%
  merge(., Superclass_and_PChem, by = c("Alignment_ID", "Sampler")) %>%
  #select(c(Alignment_ID, Uptake_Rate, Class, Sampler, LogKOA_pred, Average_Mz)) %>%
  filter(Superclass != "Nucleosides, nucleotides, and analogues")

NumRDF_Superclass <- URs_SuperclassesPchem %>%
  group_by(Sorbent) %>%
  dplyr::mutate(N_RDF_Sorbent = n()) %>%
  ungroup() %>%
  group_by(Sorbent, N_RDF_Sorbent, Superclass) %>%
  dplyr::summarize(n = n()) #number of features with uptake rates in each superclass

#simplify by sorbent
NumRDF_Superclass_PDMS <- NumRDF_Superclass %>%
  filter(Sorbent == "PDMS")

NumRDF_Superclass_GFF <- NumRDF_Superclass %>%
  filter(Sorbent == "GFF")

#table for export
PDMS_ClassTable <- PAS_PeakAreas %>%
  filter(Sorbent == "PDMS") %>%
  group_by(Superclass) %>%
  dplyr::summarize(N_features = n(),
                   TotalPA = sum(Avg_Area)) %>%
  merge(., NumRDF_Superclass_PDMS, by = "Superclass", all = TRUE) %>%
  select(-Sorbent) %>%
  mutate(proportion = 100*n/N_RDF_Sorbent) %>%
  mutate(n = ifelse(is.na(n),
                    0,
                    n),
         N_RDF_Sorbent = ifelse(is.na(N_RDF_Sorbent),
                    0,
                    N_RDF_Sorbent))

DF_to_CSV(PDMS_ClassTable)

GFF_ClassTable <- PAS_PeakAreas %>%
  filter(Sorbent == "GFF") %>%
  group_by(Superclass) %>%
  dplyr::summarize(N_features = n(),
                   TotalPA = sum(Avg_Area)) %>%
  merge(., NumRDF_Superclass_GFF, by = "Superclass", all = TRUE) %>%
  select(-Sorbent) %>%
  mutate(proportion = 100*n/N_RDF_Sorbent) %>%
  mutate(n = ifelse(is.na(n),
                    0,
                    n),
         N_RDF_Sorbent = ifelse(is.na(N_RDF_Sorbent),
                                0,
                                N_RDF_Sorbent))

DF_to_CSV(GFF_ClassTable)

#plot of superclasses with proportion of features with uptake rates
Superclass_With_RDF_Both <- PAS_PeakAreas%>%
  merge(., NumRDF_Superclass, by = c("Sorbent", "Superclass"), all = TRUE) %>%
  filter(N_Superclass > 10 | n > 0) %>%
  mutate(proportion = 100*n/N_RDF_Sorbent) %>%
  mutate(total_proportion = 100*N_Superclass/N_Sorbent) %>%
  ggplot(aes(x = reorder(str_wrap(Superclass, 14), N_Superclass), fill = proportion))+
  geom_hline(aes(yintercept = y), data.frame(y = c(0:5) * 100), color = "black") + 
  geom_bar(position = "dodge", show.legend = TRUE)+
  # Lollipop shaft for mean gain per region
  geom_segment(aes(x = reorder(str_wrap(Superclass, 14), Total_Area), y = 0, 
                   xend = reorder(str_wrap(Superclass, 14), Total_Area), 
                   yend = 500), linetype = "dashed", color = "black") +
  #scale_y_log10()+
  scale_y_continuous(limits = c(-50, 550), expand = c(0, 0))+
  ggplot2::annotate("text", x = 0.45, y = 120, label = "100", size = 4.5) +
  ggplot2::annotate("text", x = 0.45, y = 220, label = "200", size = 4.5) +
  ggplot2::annotate("text", x = 0.45, y = 320, label = "300", size = 4.5) +
  ggplot2::annotate("text", x = 0.45, y = 420, label = "400", size = 4.5) +
  ggplot2::annotate("text", x = 0.45, y = 520, label = "500", size = 4.5) +
  theme_bw()+
  scale_fill_distiller(palette = "Blues", direction = 1)+
  coord_radial(start = -0.04, expand = TRUE, clip = "off")+
  facet_wrap(~Sorbent, ncol = 1, nrow = 2)+
  labs(x = element_blank(), fill  = str_wrap("Proportion of Features with Uptake Rates (%)", 20))+ #axis titles#
  theme(#axis.title.x = element_text(size = 20), 
    #axis.title.y = element_text(size = 20), #adjust text size for axis titles#
    axis.text.x = element_text(colour = "black", size = 15),
    #axis.text.y = element_text(colour = "black", size = 7), #adjust axis text size and colour#
    legend.title = element_text(size = 15),
    #legend.margin = margin(2,2,2,2, "cm"),
    legend.background = element_rect(colour = "black", linewidth = 1),
    legend.frame = element_rect(colour = "black", linewidth = 1),
    legend.ticks = element_line(colour = "black", linewidth = 1),
    axis.text.y = element_blank(),
    legend.text = element_text(size = 15), #adjust legend title and text#
    legend.key.width = unit(2, "cm"),
    legend.key.height = unit(1, "cm"),
    legend.position = "bottom",
    plot.title = element_text(size = 15), #adjust plot title text size#
    plot.background = element_rect(fill = "white"),
    strip.text = element_text(size = 15), #adjust strip text if faceting#
    #axis.ticks = element_line(colour="black", linewidth = 1), #adjust tick mark size and colour#
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    #axis.line = element_line(colour = "black", linewidth = 1), #adjust axis size and colour#
    axis.line = element_blank(),
    panel.border = element_blank(), #remove border#
    #panel.spacing.x = unit(8, "cm"),
    panel.grid.major=element_blank(), #remove gridlines#
    panel.grid.minor=element_blank(), #remove gridlines#
    strip.background = element_rect(colour = "black", fill = "white", linewidth = 1))
Superclass_With_RDF_Both

####Representative uptake curves####

Rep_Curve_List <- c(1367, 1814, 4516, 4841) #overlapping features between sorbents

Rep_Curves <- Both_PAS %>%
  filter(Alignment_ID %in% Rep_Curve_List)

Rep_Curves$Alignment_ID <- factor(Rep_Curves$Alignment_ID,
                                  levels = c(1367, 1814, 4516, 4841),
                                  labels = c("Feature 1367", "Feature 1814", "Feature 4516", "Feature 4841")) #formatting

Rep_Curves_Plot <- Rep_Curves %>%
  ggplot(aes(x = Days, y = Vol, colour = PAS_Type, shape = PAS_Type)) +
  geom_point(size = 4)+ #use scatterplot#
  scale_color_brewer(palette = "Accent") +
  scale_y_continuous(breaks = seq(0, 25, by = 5), expand = c(0, 0), limits = c(0, 25)) +
  scale_x_continuous(breaks = seq(0, 20, by = 4), expand = c(0, 0, 0.1, 0), limits = c(0, 20)) +
  theme_bw()+ #aesthetics#
  facet_wrap(~Alignment_ID, nrow = 1, ncol = 4, scales = "free")+
  labs(x = "Deployment time (days)", y = expression(paste("Equivalent air volume (m"^3*")")),
       colour = "PAS Sorbent", shape = "PAS Sorbent") +
  geom_smooth(method = lm, formula = y ~ 0 + x, se = FALSE)+ #add linear regression#
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank())+ #remove gridlines#
  theme(axis.title.x = element_text(size = 20), 
        axis.title.y = element_text(size = 20), #adjust text size for axis titles#
        axis.text.x = element_text(colour="black", size = 15),
        axis.text.y = element_text(colour="black", size = 15), #adjust axis text size and colour#
        legend.title = element_text(size = 15), 
        legend.text = element_text(size = 15), #adjust legend title and text#
        plot.title = element_text(size = 15), #adjust plot title text size#
        strip.text = element_text(size = 12), #adjust strip text if faceting#
        axis.ticks = element_line(colour="black", linewidth=1), #adjust tick mark size and colour#
        axis.line = element_line(colour="black", linewidth=1), #adjust axis size and colour#
        panel.border = element_blank(), #remove border#
        strip.background = element_rect(colour="black", fill="white"),
        legend.position = "top") #black and white background#
Rep_Curves_Plot

#generic uptake rates to be in manuscript
Generic_PDMS = mean(PDMS_UR_fullfilt_IQR$Uptake_Rate)
Generic_GFF = mean(GFF_UR_fullfilt_IQR$Uptake_Rate)

Generic_PDMS_SD = sd(PDMS_UR_fullfilt_IQR$Uptake_Rate)
Generic_GFF_SD = sd(GFF_UR_fullfilt_IQR$Uptake_Rate)

#features with uptake rates on both samplers#

PDMS_Overlap <- NewFilt_Both_PAS %>%
  filter(Alignment_ID %in% GFF_UR_fullfilt$Alignment_ID == TRUE & Sorbent == "PDMS")

GFF_Overlap <- NewFilt_Both_PAS %>%
  filter(Alignment_ID %in% PDMS_UR_fullfilt$Alignment_ID == TRUE & Sorbent == "GFF")

URs_on_both <- bind_rows(PDMS_Overlap, GFF_Overlap)

#plot uptake rates from each sorbent against one another
PDMSvGFF <- URs_on_both %>%
  select(c(Alignment_ID, Sorbent, Uptake_Rate)) %>%
  pivot_wider(names_from = Sorbent, values_from = Uptake_Rate) %>%
  ggplot(aes(x = GFF, y = PDMS)) +
  geom_point(size = 4)+ #use scatterplot#
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", linewidth = 1) +
  #geom_abline(intercept = 0, slope = 0.5, linetype = "dashed", linewidth = 1) +
  #geom_abline(intercept = 0, slope = 2, linetype = "dashed", linewidth = 1) +
  #scale_color_brewer(palette = "Accent") +
  theme_bw()+ #aesthetics#
  scale_y_continuous(breaks = seq(0, 2, by = 0.5), expand = c(0, 0, 0.1, 0), limits = c(0, 2)) +
  scale_x_continuous(breaks = seq(0, 2, by = 0.5), expand = c(0, 0, 0.1, 0), limits = c(0, 2)) +
  labs(x = expression(paste("Rate from GFF(m"^"3"*" day"^"-1"*" dm"^"-2"*")")),
       y = expression(paste("Rate from PDMS (m"^"3"*" day"^"-1"*" dm"^"-2"*")"))) +
  #geom_smooth(method = lm, formula = 1 ~ 1, se = FALSE) + #add linear regression#
  #stat_regline_equation(aes(label = after_stat(eq.label)),
    #                    formula = y ~ x, hjust = 0, vjust = 0.25)+ #add line equation#
  #stat_regline_equation(aes(label = after_stat(rr.label)),
    #                    formula = y ~ x, hjust = 0, vjust = 1.75)+ #add r-squared#
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank())+ #remove gridlines#
  theme(axis.title.x = element_text(size = 20), 
        axis.title.y = element_text(size = 20), #adjust text size for axis titles#
        axis.text.x = element_text(colour="black", size = 15),
        axis.text.y = element_text(colour="black", size = 15), #adjust axis text size and colour#
        legend.title = element_text(size = 15), 
        legend.text = element_text(size = 15), #adjust legend title and text#
        plot.title = element_text(size = 15), #adjust plot title text size#
        strip.text = element_text(size = 12), #adjust strip text if faceting#
        axis.ticks = element_line(colour="black", linewidth=1), #adjust tick mark size and colour#
        axis.line = element_line(colour="black", linewidth=1), #adjust axis size and colour#
        panel.border = element_blank(), #remove border#
        strip.background = element_rect(colour="black", fill="white")) #black and white background#
PDMSvGFF

#uptake rates vs KOA
KOAvUR <- NewFilt_Both_PAS %>%
  dplyr::rename(Sampler = Sorbent) %>%
  merge(., All_Features_KOAadded, by = c("Alignment_ID", "Sampler")) %>%
  filter(!is.na(LogKOA_pred)) %>%
  ggplot(aes(x = LogKOA_pred, y = Uptake_Rate, colour = Sampler)) +
  scale_colour_brewer(palette = "Accent") +
  #geom_hline(yintercept = Generic_GFF, linetype = "dashed", linewidth = 1) +
  #geom_hline(yintercept = (Generic_GFF + 2*Generic_GFF_SD), linetype = "dashed", colour = "red", linewidth = 1) +
  #geom_hline(yintercept = (Generic_GFF - 2*Generic_GFF_SD), linetype = "dashed", colour = "red", linewidth = 1) +
  geom_point(size = 5, alpha = 0.8)+ #use scatterplot#
  scale_y_continuous(breaks = seq(0, 10, by = 1), expand = c(0.01, 0, 0, 0), limits = c(0, 10)) +
  scale_x_continuous(breaks = seq(3, 11, by = 1), expand = c(0, 0), limits = c(3, 11)) +
  theme_bw()+ #aesthetics#
  labs(x = expression(paste("Predicted log"[10]*" K"[OA]*"")),
       y = expression(paste("Uptake Rate (m"^"3"*" day"^"-1"*" dm"^"-2"*")")),
       colour = "Sorbent Type") +
  #geom_smooth(method = lm, formula = y ~ x, se = FALSE) + #add linear regression#
  #stat_regline_equation(aes(label = after_stat(eq.label)), show.legend = FALSE,
  #                      formula = y ~ x, hjust = 0, vjust = 0.25)+ #add line equation#
  #stat_regline_equation(aes(label = after_stat(rr.label)), show.legend = FALSE,
  #                      formula = y ~ x, hjust = 0, vjust = 1.75)+ #add r-squared#
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank())+ #remove gridlines#
  theme(axis.title.x = element_text(size = 20), 
        axis.title.y = element_text(size = 20), #adjust text size for axis titles#
        axis.text.x = element_text(colour="black", size = 15),
        axis.text.y = element_text(colour="black", size = 15), #adjust axis text size and colour#
        legend.title = element_text(colour="black", size = 20), 
        legend.text = element_text(colour="black", size = 20), #adjust legend title and text#
        plot.title = element_text(size = 15), #adjust plot title text size#
        strip.text = element_text(size = 12), #adjust strip text if faceting#
        axis.ticks = element_line(colour="black", linewidth=1), #adjust tick mark size and colour#
        axis.line = element_line(colour="black", linewidth=1), #adjust axis size and colour#
        panel.border = element_blank(), #remove border#
        strip.background = element_rect(colour="black", fill="white")) #black and white background#
KOAvUR

#uptake rate vs m/z
MZvUR <- NewFilt_Both_PAS %>%
  dplyr::rename(Sampler = Sorbent) %>%
  merge(., All_Features_KOAadded, by = c("Alignment_ID", "Sampler")) %>%
  #filter(Uptake_Rate <=2) %>%
  ggplot(aes(x = Average_Mz, y = Uptake_Rate, colour = Sampler)) +
  #geom_hline(yintercept = Generic_PDMS, linetype = "dashed", linewidth = 1) +
  #geom_hline(yintercept = (Generic_PDMS + 2*Generic_PDMS_SD), linetype = "dashed", colour = "red", linewidth = 1) +
  #geom_hline(yintercept = (Generic_PDMS - 2*Generic_PDMS_SD), linetype = "dashed", colour = "red", linewidth = 1) +
  geom_point(size = 5, alpha = 0.8)+ #use scatterplot#
  scale_colour_brewer(palette = "Accent") +
  scale_y_continuous(breaks = seq(0, 10, by = 1), expand = c(0.01, 0, 0, 0), limits = c(0, 10)) +
  scale_x_continuous(breaks = seq(100, 650, by = 50), expand = c(0, 0), limits = c(90, 650)) +
  theme_bw()+ #aesthetics#
  labs(x = "m/z", y = expression(paste("Uptake Rate (m"^"3"*" day"^"-1"*" dm"^"-2"*")")),
       colour = "Sorbent Type") +
  #geom_smooth(method = lm, formula = y ~ x, se = FALSE) + #add linear regression#
  #stat_regline_equation(aes(label = after_stat(eq.label)), show.legend = FALSE,
  #                      formula = y ~ x, hjust = 0, vjust = 0.25)+ #add line equation#
  #stat_regline_equation(aes(label = after_stat(rr.label)), show.legend = FALSE,
  #                      formula = y ~ x, hjust = 0, vjust = 1.75)+ #add r-squared#
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank())+ #remove gridlines#
  theme(axis.title.x = element_text(size = 20), 
        axis.title.y = element_text(size = 20), #adjust text size for axis titles#
        axis.text.x = element_text(colour="black", size = 15),
        axis.text.y = element_text(colour="black", size = 15), #adjust axis text size and colour#
        legend.title = element_text(colour="black", size = 20), 
        legend.text = element_text(colour="black", size = 20), #adjust legend title and text#
        plot.title = element_text(size = 15), #adjust plot title text size#
        strip.text = element_text(size = 12), #adjust strip text if faceting#
        axis.ticks = element_line(colour="black", linewidth=1), #adjust tick mark size and colour#
        axis.line = element_line(colour="black", linewidth=1), #adjust axis size and colour#
        panel.border = element_blank(), #remove border#
        strip.background = element_rect(colour="black", fill="white")) #black and white background#
MZvUR

#used to determine new generic uptake for GFF based on pchem properties
URs_and_PChem <- NewFilt_Both_PAS %>%
  dplyr::rename(Sampler = Sorbent) %>%
  filter(Sampler == "GFF") %>%
  merge(., All_Features_KOAadded, by = c("Alignment_ID", "Sampler")) %>%
  mutate(LogKOA_pred = ifelse(is.na(LogKOA_pred),
                              0,
                              LogKOA_pred)) %>%
  filter(Average_Mz <= 250 & Uptake_Rate <=3 & LogKOA_pred <8.5) %>%
  group_by(Sampler) %>%
  dplyr::mutate(Mean = mean(Uptake_Rate),
                STDEV = sd(Uptake_Rate)) %>%
  dplyr::mutate(RSD = 100*STDEV/Mean) %>%
  ungroup()

####-----Regression of peak areas-----####

PDMS_AllPeakAreas <- uPDMS %>%
  filter(Days_Detected_D >= 4) %>%
  group_by(Alignment_ID) %>%
  dplyr::mutate(n = n()) %>%
  dplyr::mutate(PercentDet = n/13) %>%
  ungroup() %>%
  filter(PercentDet >= 0.5) %>%
  select(c(Alignment_ID, Area, Batch, Sampler)) %>%
  mutate(Days = ifelse(Batch == 1,
                       4,
                       ifelse(Batch == 2,
                              8,
                              ifelse(Batch == 3,
                                     12,
                                     ifelse(Batch == 4,
                                            17,
                                            20))))) %>%
  
  dplyr::rename(Vol = Area) %>%
  EqAirVol_Filtering() %>%
  dplyr::rename(Area = Vol)

Reg_PDMS_Areas <- PDMS_AllPeakAreas %>%
  filter(Alignment_ID %in% Dedup_A$Alignment_ID == FALSE | Alignment_ID %in% PDMS_UR_fullfilt$Alignment_ID) %>%
  group_by(Alignment_ID)%>% #cluster by Alignment_ID#
  do(., tidy(lm(Area ~ 0 + Days, data = .))) #for each group, calculate regression statistics forcing through the origin#
summary(Reg_PDMS_Areas)

RD_PDMS_Area <- Reg_PDMS_Areas %>%
  filter(Alignment_ID %in% PDMS_UR_fullfilt$Alignment_ID)

#Calculate R-squared#
PDMS_Area_Rsq <- PDMS_AllPeakAreas%>%
  split(.$Alignment_ID)%>% #split by Alignment_ID; alternatively, use PDMS_split made previously#
  map(~lm(Area ~ 0 + Days, data = .))%>% #perform regression#
  map(summary)%>% 
  map_dbl("r.squared") #calculate R-squared#
PDMS_Area_Rsq_DF <- as.data.frame(PDMS_Area_Rsq) #convert to data frame so it's usable#

PDMS_Area_Rsq_DF <- tibble::rownames_to_column(PDMS_Area_Rsq_DF, "Alignment_ID") #change row headers into a column so it can be merged#

#Add R-squared to regression table#
Reg_PDMS_Areas_Rsq <- merge(Reg_PDMS_Areas, PDMS_Area_Rsq_DF, by = "Alignment_ID")%>%
  dplyr::rename("Rsquared" = PDMS_Area_Rsq)

#Perform R-squared test#
Reg_PDMS_Areas_Rsq[,'test'] = Reg_PDMS_Areas_Rsq$statistic*Reg_PDMS_Areas_Rsq$Rsquared

#Filter out R-squared less than 0.8#
PDMS_Area_Rsq_filtered <- Reg_PDMS_Areas_Rsq%>%
  filter(Rsquared >= 0.8)%>%
  mutate("Sorbent" = "PDMS")

#Reformat for use in plots#
PDMS_Areas_stats <- PDMS_Area_Rsq_filtered%>%
  select(c(Alignment_ID, term, estimate, std.error, Rsquared, Sorbent))%>%
  pivot_wider(names_from = term, values_from = estimate)%>%
  dplyr::rename("Uptake_Rate" = Days)%>%
  select(Alignment_ID, Uptake_Rate, std.error, Rsquared, Sorbent)

GFF_AllPeakAreas <- uGFF %>%
  filter(Days_Detected_D >= 4) %>%
  group_by(Alignment_ID) %>%
  dplyr::mutate(n = n()) %>%
  dplyr::mutate(PercentDet = n/13) %>%
  ungroup() %>%
  filter(PercentDet >= 0.5) %>%
  select(c(Alignment_ID, Area, Batch, Sampler)) %>%
  mutate(Days = ifelse(Batch == 1,
                       4,
                       ifelse(Batch == 2,
                              8,
                              ifelse(Batch == 3,
                                     12,
                                     ifelse(Batch == 4,
                                            17,
                                            20))))) %>%
  dplyr::rename(Vol = Area) %>%
  EqAirVol_Filtering() %>%
  dplyr::rename(Area = Vol)

Reg_GFF_Areas <- GFF_AllPeakAreas %>%
  filter(Alignment_ID %in% Dedup_A$Alignment_ID == FALSE | Alignment_ID %in% GFF_UR_fullfilt$Alignment_ID) %>%
  group_by(Alignment_ID) %>% #cluster by Alignment_ID#
  do(., tidy(lm(Area ~ 0 + Days, data = .))) #for each group, calculate regression statistics forcing through the origin#
summary(Reg_GFF_Areas)

RD_GFF_Area <- Reg_GFF_Areas %>%
  filter(Alignment_ID %in% GFF_UR_fullfilt$Alignment_ID)

RD_ActiveConc_GFF <- AllActiveConc %>%
  filter(Alignment_ID %in% NewFilt_Both_PAS$Alignment_ID) %>%
  select(c(Alignment_ID, Avg_Conc)) %>%
  merge(., RD_GFF_Area, by = "Alignment_ID")

#Calculate R-squared#
GFF_Area_Rsq <- GFF_AllPeakAreas%>%
  split(.$Alignment_ID)%>% #split by Alignment_ID; alternatively, use GFF_split made previously#
  map(~lm(Area ~ 0 + Days, data = .))%>% #perform regression#
  map(summary)%>% 
  map_dbl("r.squared") #calculate R-squared#
GFF_Area_Rsq_DF <- as.data.frame(GFF_Area_Rsq) #convert to data frame so it's usable#

GFF_Area_Rsq_DF <- tibble::rownames_to_column(GFF_Area_Rsq_DF, "Alignment_ID") #change row headers into a column so it can be merged#

#Add R-squared to regression table#
Reg_GFF_Areas_Rsq <- merge(Reg_GFF_Areas, GFF_Area_Rsq_DF, by = "Alignment_ID")%>%
  dplyr::rename("Rsquared" = GFF_Area_Rsq)

#Perform R-squared test#
Reg_GFF_Areas_Rsq[,'test'] = Reg_GFF_Areas_Rsq$statistic*Reg_GFF_Areas_Rsq$Rsquared

#Filter out R-squared less than 0.8#
GFF_Area_Rsq_filtered <- Reg_GFF_Areas_Rsq%>%
  filter(Rsquared >= 0.8)%>%
  mutate("Sorbent" = "GFF")

#Reformat for use in plots#
GFF_Areas_stats <- GFF_Area_Rsq_filtered%>%
  select(c(Alignment_ID, term, estimate, std.error, Rsquared, Sorbent))%>%
  pivot_wider(names_from = term, values_from = estimate)%>%
  dplyr::rename("Uptake_Rate" = Days)%>%
  select(Alignment_ID, Uptake_Rate, std.error, Rsquared, Sorbent)

DF_to_CSV(PDMS_Areas_stats)
DF_to_CSV(GFF_Areas_stats)

#overlapping features with linear peak areas
PeakAreaOverlap <- PDMS_Areas_stats %>%
  filter(Alignment_ID %in% GFF_Areas_stats$Alignment_ID)

NumPAFeatures <- length(unique(PDMS_Areas_stats$Alignment_ID)) + length(unique(GFF_Areas_stats$Alignment_ID))
NumPAFeaturesUnique <- bind_rows(PDMS_Areas_stats, GFF_Areas_stats) %>%
  dplyr::summarize(n = length(unique(.$Alignment_ID)))

#Representative uptake curves for peak areas#

Rep_PACurve_List <- c(630, 1814, 4266, 4345, 5828, 9555, 7370, 13719) #overlapping peak area uptake curves bewteen sorbents

Rep_PACurves <- bind_rows(PDMS_AllPeakAreas, GFF_AllPeakAreas) %>%
  filter(Alignment_ID %in% Rep_PACurve_List)

Rep_PACurves$Sampler <- factor(Rep_PACurves$Sampler,
                                    levels = c("PDMS", "GFF"),
                                    labels = c("PDMS", "GFF")) #reorder layers

Rep_PACurves$Alignment_ID <- factor(Rep_PACurves$Alignment_ID,
                                    levels = c(630, 1814, 4266, 4345, 5828, 9555, 7370, 13719),
                                    labels = c("Feature 630", "Feature 1814", "Feature 4266", "Feature 4345", 
                                               "Feature 5828", "Feature 9555", "Feature 7370", "Feature 13719")) #formatting

Rep_PACurves_Plot <- Rep_PACurves %>%
  ggplot(aes(x = Days, y = Area, colour = Sampler, shape = Sampler)) +
  geom_point(size = 4)+ #use scatterplot#
  scale_color_brewer(palette = "Accent") +
  #scale_y_continuous(breaks = seq(0, 25, by = 5), expand = c(0, 0), limits = c(0, 25)) +
  scale_x_continuous(breaks = seq(0, 20, by = 4), expand = c(0, 0, 0.1, 0), limits = c(0, 20)) +
  theme_bw()+ #aesthetics#
  facet_wrap(~Alignment_ID, scales = "free", nrow = 2, ncol = 4)+
  labs(x = "Deployment time (days)", y = "Peak Area (PAU)",
       colour = "PAS Sorbent", shape = "PAS Sorbent") +
  geom_smooth(method = lm, formula = y ~ 0 + x, se = FALSE)+ #add linear regression#
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank())+ #remove gridlines#
  theme(axis.title.x = element_text(size = 20), 
        axis.title.y = element_text(size = 20), #adjust text size for axis titles#
        axis.text.x = element_text(colour="black", size = 15),
        axis.text.y = element_text(colour="black", size = 15), #adjust axis text size and colour#
        legend.title = element_text(size = 15), 
        legend.text = element_text(size = 15), #adjust legend title and text#
        plot.title = element_text(size = 15), #adjust plot title text size#
        strip.text = element_text(size = 12), #adjust strip text if faceting#
        axis.ticks = element_line(colour="black", linewidth=1), #adjust tick mark size and colour#
        axis.line = element_line(colour="black", linewidth=1), #adjust axis size and colour#
        panel.border = element_blank(), #remove border#
        strip.background = element_rect(colour="black", fill="white")) #black and white background#
Rep_PACurves_Plot

save_plot(plot = Rep_PACurves_Plot, "Figures/Figure_4.svg", base_height = 8, base_width = 15)

####-----Finalized Figures for Manuscript-----####

preFigure2 <- ggarrange(NewFilt_IQR_BP_Rates, PDMSvGFF,
                        labels = c("B", "C"), nrow = 1, ncol = 2, hjust = c(-4, -4), align = "h",
                        font.label = list(colour = "black", size = 25))
preFigure2 #preparing part of figure 2

Figure2 <- ggarrange(Rep_Curves_Plot, preFigure2,
                     labels = c("A", ""), hjust = c(-3.7, 0), vjust = c(2.5, NA, 1.5),
                     nrow = 2,
                     font.label = list(colour = "black", size = 25))
Figure2 #figure 2

save_plot(plot = Figure2, "Figures/Figure_2.svg", base_height = 10, base_width = 10)

Figure5 <- ggarrange(KOAvUR, MZvUR, Detection_Comparison, Detection_Comparison_MZ,
                     labels = c("A", "B", "C", "D"),
                     nrow = 2, ncol = 2, hjust = -3,
                     font.label = list(colour = "black", size = 35),
                     common.legend = TRUE, legend = "bottom")
Figure5 #figure 5

save_plot(plot = Figure5, "Figures/Figure_5.svg", base_height = 10, base_width = 15)

#used to create figure 6 in PowerPoint due to formatting issues with export
save_plot(plot = Superclass_With_PeakArea_Both, "Figures/Superclass_With_PeakArea_Both.svg", base_height = 12, base_width = 8)
save_plot(plot = Superclass_With_RDF_Both, "Figures/Superclass_With_RDF_Both.svg", base_height = 12, base_width = 8)

#number of features with uptake rates and KOAs for each sorbent
PDMS_KOA <- All_Features_KOAadded %>%
  filter(Alignment_ID %in% PDMS_UR_fullfilt$Alignment_ID
         & !is.na(LogKOA_pred))

GFF_KOA <- All_Features_KOAadded %>%
  filter(Alignment_ID %in% GFF_UR_fullfilt$Alignment_ID
         & !is.na(LogKOA_pred))

#counts number of unique features
NumFeatures <- function(df){
  Num <- length(unique(df$Alignment_ID))
  return(Num)
}

BothKOA_PAS_URs <- bind_rows(PDMS_KOA,GFF_KOA)

NumFeatures(BothKOA_PAS_URs)
NumFeatures(PDMS_KOA)
NumFeatures(GFF_KOA)

#make a table with uptake rates and pchem; note some features will not have uptake rates carried into this, add them manually
RDF_Pchem_Table <- merge(BothKOA_PAS_URs, URs_SuperclassesPchem, by = c("Alignment_ID", "Sampler", "LogKOA_pred", "Average_Mz"), all.x = TRUE) %>%
  filter(Sampler != "Active") %>%
  select(c(Alignment_ID, Sampler, LogKOA_pred, Average_Mz, Superclass, Uptake_Rate)) %>%
  group_by(Alignment_ID, Sampler, LogKOA_pred, Average_Mz, Superclass, Uptake_Rate) %>%
  slice(1) %>%
  ungroup() %>%
  filter((Alignment_ID %in% PDMS_KOA$Alignment_ID & Sampler == "PDMS") | Alignment_ID %in% GFF_KOA$Alignment_ID & Sampler == "GFF")

DF_to_CSV(RDF_Pchem_Table)

#outlier uptake rates to be starred in the regression table
Outliers <- NewFilt_Both_PAS %>%
  filter(Sorbent == "GFF" | Alignment_ID %in% PDMS_UR_fullfilt_IQR$Alignment_ID == FALSE) %>%
  filter(Sorbent == "PDMS" | Alignment_ID %in% GFF_UR_fullfilt_IQR$Alignment_ID == FALSE)

#features with stable AAS measurements and overlapping with PAS
ActiveGood <- AllActiveConc %>%
  filter(Alignment_ID %in% uPDMS$Alignment_ID | Alignment_ID %in% uGFF$Alignment_ID)

#AAS measurements for PDMS features with uptake rates
PDMS_AConc <- AllActiveConc %>%
  filter(Alignment_ID %in% PDMS_UR_fullfilt$Alignment_ID) %>%
  mutate(Avg_Conc = signif(Avg_Conc, digits=3),
         Conc_SD = signif(Conc_SD, digits=3)) %>%
  select(c(Alignment_ID, Avg_Conc, Conc_SD))

#AAS measurements for GFF features with uptake rates
GFF_AConc <- AllActiveConc %>%
  filter(Alignment_ID %in% GFF_UR_fullfilt$Alignment_ID) %>%
  mutate(Avg_Conc = signif(Avg_Conc, digits=3),
         Conc_SD = signif(Conc_SD, digits=3)) %>%
  select(c(Alignment_ID, Avg_Conc, Conc_SD))

DF_to_CSV(PDMS_AConc)
DF_to_CSV(GFF_AConc)

####Veq Tables####

#Veq done again to retain Sample column this time#
PDMS_UptakeTable <- join_all(list(AllActiveConc, uPDMS), by = "Alignment_ID", type = "full") %>% #combine the PDMS data with the AAS concentrations
  mutate("Days" = ifelse(Batch == 1,
                         4,
                         ifelse(Batch == 2,
                                8,
                                ifelse(Batch == 3, 
                                       12,
                                       ifelse(Batch == 4,
                                              17,
                                              20))))) %>% #use batch to create a column denoting how long each sample was exposed#
  mutate("Vol" = (Area / ifelse(Batch == 1,
                                Avg_Conc_B1,
                                ifelse(Batch == 2,
                                       Avg_Conc_B2,
                                       ifelse(Batch == 3,
                                              Avg_Conc_B3,
                                              ifelse(Batch == 4,
                                                     Avg_Conc_B4,
                                                     Avg_Conc_B5)))))/1000) %>% #calculate PDMS EAV using the appropriate AAS conc#
  mutate(Vol = Vol/RSA_PDMS) %>% #normalize to surface area#
  subset(select = c(Alignment_ID, Average_Mz, Metabolite_name, INCHIKEY, `2D_INCHIKEY`,
                    Alignment_ID_INCHIKEY, Formula, SMILES, level, Sampler, Batch,
                    Days, Vol, Sample)) %>%
  filter(!is.na(Vol)) #remove NAs#

PDMS_Uptake_statsTable <- EqAirVol_Filtering(PDMS_UptakeTable)

PDMS_EAVs <- PDMS_Uptake_statsTable %>%
  filter(Alignment_ID %in% PDMS_UR_fullfilt$Alignment_ID) %>%
  mutate(Vol = signif(Vol, digits=3)) %>%
  mutate(Days = as.character(Days)) %>%
  mutate(Days = paste("Day",Days)) %>%
  mutate(Replicate = ifelse(grepl("a", Sample),
                            "a",
                            ifelse(grepl("b", Sample),
                                   "b",
                                   "c"))) %>% #add replicate
  select(c(Alignment_ID, Replicate, Days, Vol)) %>%
  pivot_wider(names_from = Days, values_from = Vol) %>%
  arrange(Alignment_ID, Replicate)

#repeat for GFF
GFF_UptakeTable <- join_all(list(AllActiveConc, uGFF), by = "Alignment_ID", type = "full")%>% #combine the GFF data with the AAS Concs
  mutate("Days" = ifelse(Batch == 1,
                         4,
                         ifelse(Batch == 2,
                                8,
                                ifelse(Batch == 3, 
                                       12,
                                       ifelse(Batch == 4,
                                              17,
                                              20)))))%>% #use batch to create a column denoting how long each sample was exposed#
  mutate("Vol" = Area / (ifelse(Batch == 1,
                                Avg_Conc_B1,
                                ifelse(Batch == 2,
                                       Avg_Conc_B2,
                                       ifelse(Batch == 3,
                                              Avg_Conc_B3,
                                              ifelse(Batch == 4,
                                                     Avg_Conc_B4,
                                                     Avg_Conc_B5)))))/1000)%>% #calculate GFF EAV using the appropriate AAS Conc#
  mutate(Vol = Vol/RSA_GFF)%>% #normalize to surface area#
  subset(select = c(Alignment_ID, Average_Mz, Metabolite_name, INCHIKEY, `2D_INCHIKEY`,
                    Alignment_ID_INCHIKEY, Formula, SMILES, level, Sampler, Batch,
                    Days, Sample, Vol)) %>%
  filter(!is.na(Vol)) #remove NAs#

#Add mean and sd#
GFF_Uptake_statsTable <- EqAirVol_Filtering(GFF_UptakeTable)

GFF_EAVs <- GFF_Uptake_statsTable %>%
  filter(Alignment_ID %in% GFF_UR_fullfilt$Alignment_ID) %>%
  mutate(Vol = signif(Vol, digits=3)) %>%
  mutate(Days = as.character(Days)) %>%
  mutate(Days = paste("Day",Days)) %>%
  mutate(Replicate = ifelse(grepl("d", Sample),
                            "d",
                            ifelse(grepl("e", Sample),
                                   "e",
                                   "f"))) %>%
  select(c(Alignment_ID, Replicate, Days, Vol)) %>%
  pivot_wider(names_from = Days, values_from = Vol) %>%
  arrange(Alignment_ID, Replicate)

DF_to_CSV(PDMS_EAVs)
DF_to_CSV(GFF_EAVs)

#important numbers for use in manuscript
Important_Numbers <- data.frame(Start = "start") %>%
  mutate(NumRDF = length(unique(NewFilt_Both_PAS$Alignment_ID)),
                            NumPDMS = length(unique(PDMS_UR_fullfilt$Alignment_ID)),
                            NumGFF = length(unique(GFF_UR_fullfilt$Alignment_ID)),
                            PDMS_Gen = unique(PDMS_UR_fullfilt_IQR$Mean),
                            PDMS_GenSD = unique(PDMS_UR_fullfilt_IQR$STDEV),
                            PDMS_GenRSD = 100*(PDMS_GenSD/PDMS_Gen),
                            GFF_Gen = unique(GFF_UR_fullfilt_IQR$Mean),
                            GFF_GenSD = unique(GFF_UR_fullfilt_IQR$STDEV),
                            GFF_GenRSD = 100*(GFF_GenSD/GFF_Gen),
                            NumFeatures = length(unique(Seperate_dedup_individ_bound$Alignment_ID)),
                            NumL2and3 = length(unique(Seperate_dedup_individ_bound_level32$Alignment_ID)),
                            NumLClasses = length(unique(classes1$Alignment_ID)),
                            NumKOA = length(unique(OPERA_all$Alignment_ID)),
                            PDMS_HighKOA = max(PDMS_KOA$LogKOA_pred),
                            PDMS_LowKOA = min(PDMS_KOA$LogKOA_pred),
                            GFF_HighKOA = max(GFF_KOA$LogKOA_pred),
                            GFF_LowKOA = min(GFF_KOA$LogKOA_pred),
                            PDMS_HighMZ = max(PDMS_KOA$Average_Mz),
                            PDMS_LowMZ = min(PDMS_KOA$Average_Mz),
                            GFF_HighMZ = max(GFF_KOA$Average_Mz),
                            GFF_LowMZ = min(GFF_KOA$Average_Mz))

KOAsinPDMS <- OPERA_all %>%
  filter(Alignment_ID %in% Dedup_D$Alignment_ID)

KOAsinGFF <- OPERA_all %>%
  filter(Alignment_ID %in% Dedup_G$Alignment_ID)

NumFeatures(Dedup_D)
NumFeatures(Dedup_G)

####Uptake Curves for SI####

PDMS_CurveList1 <- PDMS_UR_fullfilt %>%
  slice(1:16)

PDMS_CurveList2 <- PDMS_UR_fullfilt %>%
  slice(17:32)

PDMS_CurveList3 <- PDMS_UR_fullfilt %>%
  slice(33:48)

PDMS_CurveList4 <- PDMS_UR_fullfilt %>%
  slice(49:59)

PDMS_Colours = Accent_Palette[1]

PDMS_Curves1 <- plotPDMS_filter %>%
  filter(Alignment_ID %in% PDMS_CurveList1$Alignment_ID) %>%
  ggplot(aes(x = Days, y = Vol, colour = PDMS_Colours)) +
  geom_point(size = 4, show.legend = FALSE)+ #use scatterplot#
  scale_colour_manual(values = PDMS_Colours) +
  #scale_y_continuous(breaks = seq(0, NULL, by = 5), expand = c(0, 0), limits = c(0, NULL)) +
  scale_x_continuous(breaks = seq(0, 20, by = 4), expand = c(0, 0, 0.1, 0), limits = c(0, 20)) +
  theme_bw()+ #aesthetics#
  facet_wrap(~Alignment_ID, nrow = 4, ncol = 4, scales = "free")+
  labs(x = "Deployment time (days)", y = expression(paste("Equivalent air volume (m"^3*")"))) +
  geom_smooth(method = lm, formula = y ~ 0 + x, se = FALSE, show.legend = FALSE)+ #add linear regression#
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank())+ #remove gridlines#
  theme(axis.title.x = element_text(size = 20), 
        axis.title.y = element_text(size = 20), #adjust text size for axis titles#
        axis.text.x = element_text(colour="black", size = 15),
        axis.text.y = element_text(colour="black", size = 15), #adjust axis text size and colour#
        legend.title = element_text(size = 15), 
        legend.text = element_text(size = 15), #adjust legend title and text#
        plot.title = element_text(size = 15), #adjust plot title text size#
        strip.text = element_text(size = 12), #adjust strip text if faceting#
        axis.ticks = element_line(colour="black", linewidth=1), #adjust tick mark size and colour#
        axis.line = element_line(colour="black", linewidth=1), #adjust axis size and colour#
        panel.border = element_blank(), #remove border#
        strip.background = element_rect(colour="black", fill="white"),
        legend.position = "top") #black and white background#
PDMS_Curves1

PDMS_Curves2 <- plotPDMS_filter %>%
  filter(Alignment_ID %in% PDMS_CurveList2$Alignment_ID) %>%
  ggplot(aes(x = Days, y = Vol, colour = PDMS_Colours)) +
  geom_point(size = 4, show.legend = FALSE)+ #use scatterplot#
  scale_colour_manual(values = PDMS_Colours) +
  #scale_y_continuous(breaks = seq(0, NULL, by = 5), expand = c(0, 0), limits = c(0, NULL)) +
  scale_x_continuous(breaks = seq(0, 20, by = 4), expand = c(0, 0, 0.1, 0), limits = c(0, 20)) +
  theme_bw()+ #aesthetics#
  facet_wrap(~Alignment_ID, nrow = 4, ncol = 4, scales = "free")+
  labs(x = "Deployment time (days)", y = expression(paste("Equivalent air volume (m"^3*")"))) +
  geom_smooth(method = lm, formula = y ~ 0 + x, se = FALSE, show.legend = FALSE)+ #add linear regression#
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank())+ #remove gridlines#
  theme(axis.title.x = element_text(size = 20), 
        axis.title.y = element_text(size = 20), #adjust text size for axis titles#
        axis.text.x = element_text(colour="black", size = 15),
        axis.text.y = element_text(colour="black", size = 15), #adjust axis text size and colour#
        legend.title = element_text(size = 15), 
        legend.text = element_text(size = 15), #adjust legend title and text#
        plot.title = element_text(size = 15), #adjust plot title text size#
        strip.text = element_text(size = 12), #adjust strip text if faceting#
        axis.ticks = element_line(colour="black", linewidth=1), #adjust tick mark size and colour#
        axis.line = element_line(colour="black", linewidth=1), #adjust axis size and colour#
        panel.border = element_blank(), #remove border#
        strip.background = element_rect(colour="black", fill="white"),
        legend.position = "top") #black and white background#
PDMS_Curves2

PDMS_Curves3 <- plotPDMS_filter %>%
  filter(Alignment_ID %in% PDMS_CurveList3$Alignment_ID) %>%
  ggplot(aes(x = Days, y = Vol, colour = PDMS_Colours)) +
  geom_point(size = 4, show.legend = FALSE)+ #use scatterplot#
  scale_colour_manual(values = PDMS_Colours) +
  #scale_y_continuous(breaks = seq(0, NULL, by = 5), expand = c(0, 0), limits = c(0, NULL)) +
  scale_x_continuous(breaks = seq(0, 20, by = 4), expand = c(0, 0, 0.1, 0), limits = c(0, 20)) +
  theme_bw()+ #aesthetics#
  facet_wrap(~Alignment_ID, nrow = 4, ncol = 4, scales = "free")+
  labs(x = "Deployment time (days)", y = expression(paste("Equivalent air volume (m"^3*")"))) +
  geom_smooth(method = lm, formula = y ~ 0 + x, se = FALSE, show.legend = FALSE)+ #add linear regression#
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank())+ #remove gridlines#
  theme(axis.title.x = element_text(size = 20), 
        axis.title.y = element_text(size = 20), #adjust text size for axis titles#
        axis.text.x = element_text(colour="black", size = 15),
        axis.text.y = element_text(colour="black", size = 15), #adjust axis text size and colour#
        legend.title = element_text(size = 15), 
        legend.text = element_text(size = 15), #adjust legend title and text#
        plot.title = element_text(size = 15), #adjust plot title text size#
        strip.text = element_text(size = 12), #adjust strip text if faceting#
        axis.ticks = element_line(colour="black", linewidth=1), #adjust tick mark size and colour#
        axis.line = element_line(colour="black", linewidth=1), #adjust axis size and colour#
        panel.border = element_blank(), #remove border#
        strip.background = element_rect(colour="black", fill="white"),
        legend.position = "top") #black and white background#
PDMS_Curves3

PDMS_Curves4 <- plotPDMS_filter %>%
  filter(Alignment_ID %in% PDMS_CurveList4$Alignment_ID) %>%
  ggplot(aes(x = Days, y = Vol, colour = PDMS_Colours)) +
  geom_point(size = 4, show.legend = FALSE)+ #use scatterplot#
  scale_colour_manual(values = PDMS_Colours) +
  #scale_y_continuous(breaks = seq(0, NULL, by = 5), expand = c(0, 0), limits = c(0, NULL)) +
  scale_x_continuous(breaks = seq(0, 20, by = 4), expand = c(0, 0, 0.1, 0), limits = c(0, 20)) +
  theme_bw()+ #aesthetics#
  facet_wrap(~Alignment_ID, nrow = 3, ncol = 4, scales = "free")+
  labs(x = "Deployment time (days)", y = expression(paste("Equivalent air volume (m"^3*")"))) +
  geom_smooth(method = lm, formula = y ~ 0 + x, se = FALSE, show.legend = FALSE)+ #add linear regression#
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank())+ #remove gridlines#
  theme(axis.title.x = element_text(size = 20), 
        axis.title.y = element_text(size = 20), #adjust text size for axis titles#
        axis.text.x = element_text(colour="black", size = 15),
        axis.text.y = element_text(colour="black", size = 15), #adjust axis text size and colour#
        legend.title = element_text(size = 15), 
        legend.text = element_text(size = 15), #adjust legend title and text#
        plot.title = element_text(size = 15), #adjust plot title text size#
        strip.text = element_text(size = 12), #adjust strip text if faceting#
        axis.ticks = element_line(colour="black", linewidth=1), #adjust tick mark size and colour#
        axis.line = element_line(colour="black", linewidth=1), #adjust axis size and colour#
        panel.border = element_blank(), #remove border#
        strip.background = element_rect(colour="black", fill="white"),
        legend.position = "top") #black and white background#
PDMS_Curves4

GFF_CurveList1 <- GFF_UR_fullfilt %>%
  slice(1:12)

GFF_CurveList2 <- GFF_UR_fullfilt %>%
  slice(13:24)

GFF_CurveList3 <- GFF_UR_fullfilt %>%
  slice(25:34)

GFF_Colours = Accent_Palette[2]

GFF_Curves1 <- plotGFF_filter %>%
  filter(Alignment_ID %in% GFF_CurveList1$Alignment_ID) %>%
  ggplot(aes(x = Days, y = Vol, colour = GFF_Colours)) +
  geom_point(size = 4, show.legend = FALSE)+ #use scatterplot#
  scale_colour_manual(values = GFF_Colours) +
  #scale_y_continuous(breaks = seq(0, NULL, by = 5), expand = c(0, 0), limits = c(0, NULL)) +
  scale_x_continuous(breaks = seq(0, 20, by = 4), expand = c(0, 0, 0.1, 0), limits = c(0, 20)) +
  theme_bw()+ #aesthetics#
  facet_wrap(~Alignment_ID, nrow = 3, ncol = 4, scales = "free")+
  labs(x = "Deployment time (days)", y = expression(paste("Equivalent air volume (m"^3*")"))) +
  geom_smooth(method = lm, formula = y ~ 0 + x, se = FALSE, show.legend = FALSE)+ #add linear regression#
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank())+ #remove gridlines#
  theme(axis.title.x = element_text(size = 20), 
        axis.title.y = element_text(size = 20), #adjust text size for axis titles#
        axis.text.x = element_text(colour="black", size = 15),
        axis.text.y = element_text(colour="black", size = 15), #adjust axis text size and colour#
        legend.title = element_text(size = 15), 
        legend.text = element_text(size = 15), #adjust legend title and text#
        plot.title = element_text(size = 15), #adjust plot title text size#
        strip.text = element_text(size = 12), #adjust strip text if faceting#
        axis.ticks = element_line(colour="black", linewidth=1), #adjust tick mark size and colour#
        axis.line = element_line(colour="black", linewidth=1), #adjust axis size and colour#
        panel.border = element_blank(), #remove border#
        strip.background = element_rect(colour="black", fill="white"),
        legend.position = "top") #black and white background#
GFF_Curves1

GFF_Curves2 <- plotGFF_filter %>%
  filter(Alignment_ID %in% GFF_CurveList2$Alignment_ID) %>%
  ggplot(aes(x = Days, y = Vol, colour = GFF_Colours)) +
  geom_point(size = 4, show.legend = FALSE)+ #use scatterplot#
  scale_colour_manual(values = GFF_Colours) +
  #scale_y_continuous(breaks = seq(0, NULL, by = 5), expand = c(0, 0), limits = c(0, NULL)) +
  scale_x_continuous(breaks = seq(0, 20, by = 4), expand = c(0, 0, 0.1, 0), limits = c(0, 20)) +
  theme_bw()+ #aesthetics#
  facet_wrap(~Alignment_ID, nrow = 3, ncol = 4, scales = "free")+
  labs(x = "Deployment time (days)", y = expression(paste("Equivalent air volume (m"^3*")"))) +
  geom_smooth(method = lm, formula = y ~ 0 + x, se = FALSE, show.legend = FALSE)+ #add linear regression#
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank())+ #remove gridlines#
  theme(axis.title.x = element_text(size = 20), 
        axis.title.y = element_text(size = 20), #adjust text size for axis titles#
        axis.text.x = element_text(colour="black", size = 15),
        axis.text.y = element_text(colour="black", size = 15), #adjust axis text size and colour#
        legend.title = element_text(size = 15), 
        legend.text = element_text(size = 15), #adjust legend title and text#
        plot.title = element_text(size = 15), #adjust plot title text size#
        strip.text = element_text(size = 12), #adjust strip text if faceting#
        axis.ticks = element_line(colour="black", linewidth=1), #adjust tick mark size and colour#
        axis.line = element_line(colour="black", linewidth=1), #adjust axis size and colour#
        panel.border = element_blank(), #remove border#
        strip.background = element_rect(colour="black", fill="white"),
        legend.position = "top") #black and white background#
GFF_Curves2

GFF_Curves3 <- plotGFF_filter %>%
  filter(Alignment_ID %in% GFF_CurveList3$Alignment_ID) %>%
  ggplot(aes(x = Days, y = Vol, colour = GFF_Colours)) +
  geom_point(size = 4, show.legend = FALSE)+ #use scatterplot#
  scale_colour_manual(values = GFF_Colours) +
  #scale_y_continuous(breaks = seq(0, NULL, by = 5), expand = c(0, 0), limits = c(0, NULL)) +
  scale_x_continuous(breaks = seq(0, 20, by = 4), expand = c(0, 0, 0.1, 0), limits = c(0, 20)) +
  theme_bw()+ #aesthetics#
  facet_wrap(~Alignment_ID, nrow = 3, ncol = 4, scales = "free")+
  labs(x = "Deployment time (days)", y = expression(paste("Equivalent air volume (m"^3*")"))) +
  geom_smooth(method = lm, formula = y ~ 0 + x, se = FALSE, show.legend = FALSE)+ #add linear regression#
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank())+ #remove gridlines#
  theme(axis.title.x = element_text(size = 20), 
        axis.title.y = element_text(size = 20), #adjust text size for axis titles#
        axis.text.x = element_text(colour="black", size = 15),
        axis.text.y = element_text(colour="black", size = 15), #adjust axis text size and colour#
        legend.title = element_text(size = 15), 
        legend.text = element_text(size = 15), #adjust legend title and text#
        plot.title = element_text(size = 15), #adjust plot title text size#
        strip.text = element_text(size = 12), #adjust strip text if faceting#
        axis.ticks = element_line(colour="black", linewidth=1), #adjust tick mark size and colour#
        axis.line = element_line(colour="black", linewidth=1), #adjust axis size and colour#
        panel.border = element_blank(), #remove border#
        strip.background = element_rect(colour="black", fill="white"),
        legend.position = "top") #black and white background#
GFF_Curves3

save_plot(plot = PDMS_Curves1, "Figures/PDMS_Curves1.svg", base_height = 16, base_width = 16)
save_plot(plot = PDMS_Curves2, "Figures/PDMS_Curves2.svg", base_height = 16, base_width = 16)
save_plot(plot = PDMS_Curves3, "Figures/PDMS_Curves3.svg", base_height = 16, base_width = 16)
save_plot(plot = PDMS_Curves4, "Figures/PDMS_Curves4.svg", base_height = 12, base_width = 16)

save_plot(plot = GFF_Curves1, "Figures/GFF_Curves1.svg", base_height = 12, base_width = 16)
save_plot(plot = GFF_Curves2, "Figures/GFF_Curves2.svg", base_height = 12, base_width = 16)
save_plot(plot = GFF_Curves3, "Figures/GFF_Curves3.svg", base_height = 12, base_width = 16)
