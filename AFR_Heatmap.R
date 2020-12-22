### HeatMaps for Manuscript ###
## December 22nd 2020 ##
library(dplyr)
library(ggplot2)
library(pheatmap)
library(cowplot)
library(RColorBrewer)
library(heatmap3)

filesAfrican <- c("~/Dropbox/Grad_School/FlexTyper/FinalData/CombinedAfrican_UniqReads_k50s10m100u.tsv")
CountDataA <- read.csv(file.path(filesAfrican), sep='\t',header=TRUE,row.names = "INfo")
names(CountDataA)[names(CountDataA) == "ERR1955404"] <- "EUR_Father"
names(CountDataA)[names(CountDataA) == "ERR1955491"] <- "EUR_Mother"
names(CountDataA)[names(CountDataA) == "ERR2304597"] <- "EUR_Child"
names(CountDataA)[names(CountDataA) == "ERR1955507"] <- "EAS_Father"
names(CountDataA)[names(CountDataA) == "ERR1955495"] <- "EAS_Mother"
names(CountDataA)[names(CountDataA) == "ERR2304569"] <- "EAS_Child"
names(CountDataA)[names(CountDataA) == "ERR1955420"] <- "AFR_Father"
names(CountDataA)[names(CountDataA) == "ERR1955443"] <- "AFR_Mother"
names(CountDataA)[names(CountDataA) == "ERR2304556"] <- "AFR_Child"

CountDataA$EUR_Father <- as.numeric(as.character(CountDataA$EUR_Father))
CountDataA$EUR_Mother <- as.numeric(as.character(CountDataA$EUR_Mother))
CountDataA$EUR_Child <- as.numeric(as.character(CountDataA$EUR_Child))
CountDataA$EAS_Father <- as.numeric(as.character(CountDataA$EAS_Father))
CountDataA$EAS_Mother <- as.numeric(as.character(CountDataA$EAS_Mother))
CountDataA$EAS_Child <- as.numeric(as.character(CountDataA$EAS_Child))
CountDataA$AFR_Father <- as.numeric(as.character(CountDataA$AFR_Father))
CountDataA$AFR_Mother <- as.numeric(as.character(CountDataA$AFR_Mother))
CountDataA$AFR_Child <- as.numeric(as.character(CountDataA$AFR_Child))

dim(CountDataA)

# Filter 
AFR_uniq <- CountDataA %>%
  filter((EUR_Father < 5) & (EUR_Mother <5 ) & (EUR_Child < 5) & 
           (EAS_Father < 5) & (EAS_Mother <5 ) & (EAS_Child < 5) & 
           (AFR_Father > 10) & (AFR_Mother > 10 ) & (AFR_Child > 10)) 

nrow(AFR_uniq)

EUR_uniq <- CountDataA %>%
  filter((AFR_Father > 5) & (AFR_Mother <5 ) & (AFR_Child < 5) & 
           (EAS_Father < 5) & (EAS_Mother <5 ) & (EAS_Child < 5) & 
           (EUR_Father > 10) & (EUR_Mother > 10 ) & (EUR_Child > 10)) 
nrow(EUR_uniq)


EAS_uniq <- CountDataA %>%
  filter((EUR_Father < 5) & (EUR_Mother <5 ) & (EUR_Child < 5) & 
           (AFR_Father < 5) & (AFR_Mother <5 ) & (AFR_Child < 5) & 
           (EAS_Father > 10) & (EAS_Mother > 10 ) & (EAS_Child > 10)) 
nrow(EAS_uniq)

CombinedUniq <- as.matrix(bind_rows(AFR_uniq, EUR_uniq, EAS_uniq))
CombinedUniq <- CombinedUniq[,c(2,3,7,1,4,9,6,5,8)]

# generate heatmap for African contigs 
afr_heatmap <- pheatmap(log(CombinedUniq+1), cluster_rows = FALSE, cluster_cols = FALSE, show_rownames = FALSE, show_colnames = TRUE, cellwidth = 18, cellheight = 0.4) 

## KIR Specific Contigs ## 

KIR_Sample_ethnicity <- read.csv("~/Dropbox/Grad_School/FlexTyper/FlexTyper_manuscript/KIR_Sample_ethnicity.csv")

#import all the files 
ERR1955404 <- read.delim("~/Dropbox/Grad_School/FlexTyper/FinalData/ERR1955404_KPIMarkers.tsv")
ERR1955420 <- read.delim("~/Dropbox/Grad_School/FlexTyper/FinalData/ERR1955420_KPIMarkers.tsv")
ERR1955443 <- read.delim("~/Dropbox/Grad_School/FlexTyper/FinalData/ERR1955443_KPIMarkers.tsv")
ERR1955491 <- read.delim("~/Dropbox/Grad_School/FlexTyper/FinalData/ERR1955491_KPIMarkers.tsv")
ERR1955495 <- read.delim("~/Dropbox/Grad_School/FlexTyper/FinalData/ERR1955495_KPIMarkers.tsv")
ERR1955507 <- read.delim("~/Dropbox/Grad_School/FlexTyper/FinalData/ERR1955507_KPIMarkers.tsv")
ERR2304556 <- read.delim("~/Dropbox/Grad_School/FlexTyper/FinalData/ERR2304556_KPIMarkers.tsv")
ERR2304569 <- read.delim("~/Dropbox/Grad_School/FlexTyper/FinalData/ERR2304569_KPIMarkers.tsv")
ERR2304597 <- read.delim("~/Dropbox/Grad_School/FlexTyper/FinalData/ERR2304597_KPIMarkers.tsv")

#remove hits less than 3
ERR1955404 <- filter(ERR1955404, count >= 3)
ERR1955420 <- filter(ERR1955420, count >= 3)
ERR1955443 <- filter(ERR1955443, count >= 3)
ERR1955491 <- filter(ERR1955491, count >= 3)
ERR1955495 <- filter(ERR1955495, count >= 3)
ERR1955507 <- filter(ERR1955507, count >= 3)
ERR2304556 <- filter(ERR2304556, count >= 3)
ERR2304569 <- filter(ERR2304569, count >= 3)
ERR2304597 <- filter(ERR2304597, count >= 3)

list_data <- list(ERR1955404,
                  ERR1955420,
                  ERR1955443,
                  ERR1955491,
                  ERR1955495,
                  ERR1955507,
                  ERR2304556,
                  ERR2304569,
                  ERR2304597)
names(list_data) <- c('ERR1955404',
                      'ERR1955420',
                      'ERR1955443',
                      'ERR1955491',
                      'ERR1955495',
                      'ERR1955507',
                      'ERR2304556',
                      'ERR2304569',
                      'ERR2304597')

for (r in 1:nrow(KIR_Sample_ethnicity)){
  s = KIR_Sample_ethnicity[r,"sample"]
  i = KIR_Sample_ethnicity[r,"info"] 
  s <- as.character(s)
  i <- as.character(i)
  h = list_data[[s]] %>% filter(info == i) %>% nrow
  a <- filter(list_data[[s]], info==i)
  c = sum(a$count)
  KIR_Sample_ethnicity[r,"hits"] = h
  KIR_Sample_ethnicity[r,"count"] = c
  KIR_Sample_ethnicity[r,"avghits"] = c/h
  
}

#change NA to 0
KIR_Sample_ethnicity[is.na(KIR_Sample_ethnicity)] <- 0

#KIR Heatmap
EUR_Father <- filter(KIR_Sample_ethnicity, sample=="ERR1955404")
EUR_Father <- EUR_Father[,c(4)]
AFR_Father <- filter(KIR_Sample_ethnicity, sample=="ERR1955420")
AFR_Father <- AFR_Father[,c(4)]
AFR_Mother <- filter(KIR_Sample_ethnicity, sample=="ERR1955443")
AFR_Mother <- AFR_Mother[,c(4)]
EUR_Mother <- filter(KIR_Sample_ethnicity, sample=="ERR1955491")
EUR_Mother <- EUR_Mother[,c(4)]
EAS_Mother <- filter(KIR_Sample_ethnicity, sample=="ERR1955495")
EAS_Mother <- EAS_Mother[,c(4)]
EAS_Father <- filter(KIR_Sample_ethnicity, sample=="ERR1955507")
EAS_Father <- EAS_Father[,c(4)]
AFR_Child <- filter(KIR_Sample_ethnicity, sample=="ERR2304556")
AFR_Child <- AFR_Child[,c(4)]
EAS_Child <- filter(KIR_Sample_ethnicity, sample=="ERR2304569")
EAS_Child <- EAS_Child[,c(4)]
EUR_Child <- filter(KIR_Sample_ethnicity, sample=="ERR2304597")
EUR_Child <- EUR_Child[,c(4)]
KIRSums = cbind(AFR_Father,
                AFR_Mother, 
                AFR_Child,
                EUR_Father,
                EUR_Mother,
                EUR_Child,
                EAS_Father, 
                EAS_Mother, 
                EAS_Child)

KIR_Heatmap <- pheatmap(log(KIRSums+1), cluster_cols = FALSE, cellwidth = 18, cellheight = 4.5, show_colnames = FALSE, 
                        labels_row = c("KIR3DS1","","","","","KIR3DL1","","","","","","","","","","","","","","","","","","","","","","","") ) 
plot_grid(afr_heatmap$gtable, KIR_Heatmap$gtable, ncol=1, align="v", axis=r) 
ggsave("heatmaps.png", width=3.5, height=4.5, dpi = 600)
