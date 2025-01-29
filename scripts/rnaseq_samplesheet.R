### Converting fetchNGS samplesheet to rnaseq input samplesheet 
## https://nf-co.re/rnaseq/3.16.1/docs/usage 
## samplesheet information above 

library(tidyverse)

fetchNGS_samplehseet <- 
  read.csv("C:/Users/EmmaStrand/MyProjects/Resilience_Biomarkers_Aquaculture/Cvirg_Pmarinus_RNAseq/data/rnaseq_samplesheets/samplesheet_fetchNGS_dataset5.csv") 

rnaseq_samplesheet <- fetchNGS_samplehseet %>%
  dplyr::select(sample, fastq_1, fastq_2) %>%
  mutate(strandedness = "auto")

rnaseq_samplesheet %>% 
  write.csv("C:/Users/EmmaStrand/MyProjects/Resilience_Biomarkers_Aquaculture/Cvirg_Pmarinus_RNAseq/data/rnaseq_samplesheets/samplesheet_rnaseq_dataset5.csv",
            row.names = FALSE)





