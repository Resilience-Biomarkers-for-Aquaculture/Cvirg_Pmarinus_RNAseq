### Converting fetchNGS samplesheet to rnaseq input samplesheet 
## https://nf-co.re/rnaseq/3.16.1/docs/usage 
## samplesheet information above 

library(tidyverse)

fetchNGS_samplehseet <- 
  read.csv("C:/Users/EmmaStrand/MyProjects/Resilience_Biomarkers_Aquaculture/samplesheet_rnaseq_Cvir_disease_set1.csv") 

rnaseq_samplesheet <- fetchNGS_samplehseet %>%
  dplyr::select(sample, fastq_1, fastq_2) %>%
  mutate(strandedness = "auto")

rnaseq_samplesheet %>% 
  write.csv("C:/Users/EmmaStrand/MyProjects/Resilience_Biomarkers_Aquaculture/samplesheet_rnaseq_Cvir_disease_set1.csv",
            row.names = FALSE)













