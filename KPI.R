#load modules
library(dplyr)
library(ggplot2)
library(pheatmap)

#gene set
list_genes <- list(KIR_3DS1, 
                  KIR_3DP1, 
                  KIR_3DP12DL4, 
                  KIR_3DL2, 
                  KIR_3DL1, 
                  KIR_2DS5, 
                  KIR_2DS4, 
                  KIR_2DS3, 
                  KIR_2DS2, 
                  KIR_2DS1, 
                  KIR_2DP1, 
                  KIR_2DL5, 
                  KIR_2DL4, 
                  KIR_2DL3, 
                  KIR_2DL2, 
                  KIR_2DL1, 
                  KIR_2DL2L32DP1, 
                  KIR_3DL32DS2, 
                  KIR_2DS3S5A2DS1, 
                  KIR_cA01tA02, 
                  KIR_cB03tA01, 
                  KIR_cA01tA01, 
                  KIR_cA01tB01, 
                  KIR_cB01tA01, 
                  KIR_cB01tB01, 
                  KIR_cB02tA01, 
                  KIR_cB02tB01)
    
names(list_genes) <- c("3DS1", 
                       "3DP1", 
                       "3DP12DL4", 
                       "3DL3", 
                       "3DL2", 
                       "3DL1", 
                       "2DS5", 
                       "2DS4", 
                       "2DS3", 
                       "2DS2", 
                       "2DS1", 
                       "2DP1", 
                       "2DL5", 
                       "2DL4", 
                       "2DL3", 
                       "2DL2", 
                       "2DL1", 
                       "2DL2L32DP1", 
                       "3DL32DS2", 
                       "2DS3S5A2DS1", 
                       "cA01tA02", 
                       "cB03tA01", 
                       "cA01tA01", 
                       "cA01tB01", 
                       "cB01tA01", 
                       "cB01tB01", 
                       "cB02tA01", 
                       "cB02tB01")

KIR_Sample_ethnicity <- read.csv("KIR_Sample_ethnicity.csv")

#import all the files 
ERR1955404 <- read.delim("ERR1955404_KPIMarkers.tsv")
ERR1955420 <- read.delim("ERR1955420_KPIMarkers.tsv")
ERR1955443 <- read.delim("ERR1955443_KPIMarkers.tsv")
ERR1955491 <- read.delim("ERR1955491_KPIMarkers.tsv")
ERR1955495 <- read.delim("ERR1955495_KPIMarkers.tsv")
ERR1955507 <- read.delim("ERR1955507_KPIMarkers.tsv")
ERR2304556 <- read.delim("ERR2304556_KPIMarkers.tsv")
ERR2304569 <- read.delim("ERR2304569_KPIMarkers.tsv")
ERR2304597 <- read.delim("ERR2304597_KPIMarkers.tsv")

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


#KIR Histogram
KIR_3DL1.ERR1955404 <- filter(ERR1955404, info=="3DL1")
KIR_3DL1.ERR1955420 <- filter(ERR1955420, info=="3DL1")
KIR_3DL1.ERR1955443 <- filter(ERR1955443, info=="3DL1")
KIR_3DL1.ERR1955491 <- filter(ERR1955491, info=="3DL1")
KIR_3DL1.ERR1955495 <- filter(ERR1955495, info=="3DL1")
KIR_3DL1.ERR1955507 <- filter(ERR1955507, info=="3DL1")
KIR_3DL1.ERR2304556 <- filter(ERR2304556, info=="3DL1")
KIR_3DL1.ERR2304569 <- filter(ERR2304569, info=="3DL1")
KIR_3DL1.ERR2304597 <- filter(ERR2304597, info=="3DL1")

KIR_3DL1.ERR1955404$id <- "ERR1955404"
KIR_3DL1.ERR1955420$id <- "ERR1955420"
KIR_3DL1.ERR1955443$id <- "ERR1955443"
KIR_3DL1.ERR1955491$id <- "ERR1955491"
KIR_3DL1.ERR1955495$id <- "ERR1955495"
KIR_3DL1.ERR1955507$id <- "ERR1955507"
KIR_3DL1.ERR2304556$id <- "ERR2304556"
KIR_3DL1.ERR2304569$id <- "ERR2304569"
KIR_3DL1.ERR2304597$id <- "ERR2304597"
KIR_3DL1 = rbind(KIR_3DL1.ERR1955404, 
                 KIR_3DL1.ERR1955420,
                 KIR_3DL1.ERR1955443,
                 KIR_3DL1.ERR1955491,
                 KIR_3DL1.ERR1955495, 
                 KIR_3DL1.ERR1955507,
                 KIR_3DL1.ERR2304556, 
                 KIR_3DL1.ERR2304569, 
                 KIR_3DL1.ERR2304597)
KIR_3DS1.ERR1955404 <- filter(ERR1955404, info=="3DS1")
KIR_3DS1.ERR1955420 <- filter(ERR1955420, info=="3DS1")
KIR_3DS1.ERR1955443 <- filter(ERR1955443, info=="3DS1")
KIR_3DS1.ERR1955491 <- filter(ERR1955491, info=="3DS1")
KIR_3DS1.ERR1955495 <- filter(ERR1955495, info=="3DS1")
KIR_3DS1.ERR1955507 <- filter(ERR1955507, info=="3DS1")
KIR_3DS1.ERR2304556 <- filter(ERR2304556, info=="3DS1")
KIR_3DS1.ERR2304569 <- filter(ERR2304569, info=="3DS1")
KIR_3DS1.ERR2304597 <- filter(ERR2304597, info=="3DS1")

KIR_3DS1.ERR1955404$id <- "ERR1955404"
KIR_3DS1.ERR1955420$id <- "ERR1955420"
KIR_3DS1.ERR1955443$id <- "ERR1955443"
KIR_3DS1.ERR1955491$id <- "ERR1955491"
KIR_3DS1.ERR1955495$id <- "ERR1955495"
KIR_3DS1.ERR1955507$id <- "ERR1955507"
KIR_3DS1.ERR2304556$id <- "ERR2304556"
KIR_3DS1.ERR2304569$id <- "ERR2304569"
KIR_3DS1.ERR2304597$id <- "ERR2304597"
KIR_3DS1 = rbind(KIR_3DS1.ERR1955404, 
                 KIR_3DS1.ERR1955420,
                 KIR_3DS1.ERR1955443,
                 KIR_3DS1.ERR1955491,
                 KIR_3DS1.ERR1955495, 
                 KIR_3DS1.ERR1955507,
                 KIR_3DS1.ERR2304556, 
                 KIR_3DS1.ERR2304569, 
                 KIR_3DS1.ERR2304597)
L <-ggplot(KIR_3DL1, aes(x = count, fill= id)) + theme(plot.title = element_text(hjust = 0.5), axis.ticks = element_blank(),  axis.line = element_line(size = 0.3, colour = "black"), axis.title.x = element_blank(),panel.background=element_blank(),  legend.position = c(0.8, 0.8),legend.title = element_blank()) + geom_histogram(aes(y = ..count..), bins=30,position="identity", alpha =0.8) + scale_x_continuous(expand = expansion(mult = c(0, 0))) + scale_y_continuous(expand = expansion(mult = c(0, 0)))
S <-ggplot(KIR_3DS1, aes(x = count, fill= id)) + theme(plot.title = element_text(hjust = 0.5), axis.ticks = element_blank(),  axis.line = element_line(size = 0.3, colour = "black"), axis.title.x = element_blank(),panel.background=element_blank(),  legend.position = c(0.8, 0.8),legend.title = element_blank()) + geom_histogram(aes(y = ..count..), bins=30,position="identity", alpha =0.8) + scale_x_continuous(expand = expansion(mult = c(0, 0))) + scale_y_continuous(expand = expansion(mult = c(0, 0)))
plot_grid(S, L, ncol=2)

#KIR Heatmap
KIRSums.ERR1955404 <- filter(KIR_Sample_ethnicity, sample=="ERR1955404")
KIRSums.ERR1955404 <- KIRSums.ERR1955404[,c(4)]
KIRSums.ERR1955420 <- filter(KIR_Sample_ethnicity, sample=="ERR1955420")
KIRSums.ERR1955420 <- KIRSums.ERR1955420[,c(4)]
KIRSums.ERR1955443 <- filter(KIR_Sample_ethnicity, sample=="ERR1955443")
KIRSums.ERR1955443 <- KIRSums.ERR1955443[,c(4)]
KIRSums.ERR1955491 <- filter(KIR_Sample_ethnicity, sample=="ERR1955491")
KIRSums.ERR1955491 <- KIRSums.ERR1955491[,c(4)]
KIRSums.ERR1955495 <- filter(KIR_Sample_ethnicity, sample=="ERR1955495")
KIRSums.ERR1955495 <- KIRSums.ERR1955495[,c(4)]
KIRSums.ERR1955507 <- filter(KIR_Sample_ethnicity, sample=="ERR1955507")
KIRSums.ERR1955507 <- KIRSums.ERR1955507[,c(4)]
KIRSums.ERR2304556 <- filter(KIR_Sample_ethnicity, sample=="ERR2304556")
KIRSums.ERR2304556 <- KIRSums.ERR2304556[,c(4)]
KIRSums.ERR2304569 <- filter(KIR_Sample_ethnicity, sample=="ERR2304569")
KIRSums.ERR2304569 <- KIRSums.ERR2304569[,c(4)]
KIRSums.ERR2304597 <- filter(KIR_Sample_ethnicity, sample=="ERR2304597")
KIRSums.ERR2304597 <- KIRSums.ERR2304597[,c(4)]
KIRSums = cbind(KIRSums.ERR1955404, 
                 KIRSums.ERR1955420,
                 KIRSums.ERR1955443,
                 KIRSums.ERR1955491,
                 KIRSums.ERR1955495, 
                 KIRSums.ERR1955507,
                 KIRSums.ERR2304556, 
                 KIRSums.ERR2304569, 
                 KIRSums.ERR2304597)
pheatmap(log(KIRSums+1), cluster_cols = FALSE)
