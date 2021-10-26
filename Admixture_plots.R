
setwd("~/Box Sync/aDNA/lions (deflami2@illinois.edu)/Cape_lions/Witt_pipeline/PCA_55")
populations <- c("Cape_lion", "Tanzania", "South_Africa", "Botswana", "Namibia", "Malawi", "Kenya", "India", "Gabon", "Ethiopia", "DRC", "Congo", "CAR", "Zambia", "Zimbabwe")
library(ggplot2)
library(zoom)

#cov = as.matrix(read.table("pca_46_all.cov"))
#cov = data.frame(cov)

pop = read.table("Sample_info_55.txt") #tab-separated file, each row is an individual, can be as simple as 1 column for individual ID, 1 column for population (or island, in my case)
##admixture bar plot:

tbl = read.table("pca_55_all.K2.a0.qopt")
tbl_df <- data.frame(tbl[,]) 
View(tbl_df)

tbl_df$pop = pop[,3] #add pop lable to column
tbl_df$ind = pop[,2] #add individual label to column
tbl_df$srt = pop[,5] #ad the order by which to sort (needs to be standardised across K values, so whichever order works best for K =2 )
#View(tbl_df)

ord <- tbl_df[order(tbl_df$srt),] #order by srt column 
#View(ord)

bp = barplot(t(as.matrix(ord[,1:2])), space = c(0.2), col = rainbow(3), border = NA, names.arg = ord$pop,las=2, font.axis=1, cex.names = 0.6)

#K=3:
tbl3 = read.table("pca_55_all.K3.a0.qopt")
tbl3_df <- data.frame(tbl3[,]) 
#View(tbl3_df)
tbl3_df$pop = pop[,3] #add pop lable to column
tbl3_df$ind = pop[,2] #add individual label to column
tbl3_df$srt = pop[,5]
#View(tbl3_df)
ord3 <- tbl3_df[order(tbl3_df$srt),]
#View(ord3)
bp3 = barplot(t(as.matrix(ord3[,1:3])), space = c(0.2), col = rainbow(4), border = NA, names.arg = ord3$pop,las=2, font.axis=1, cex.names = 0.6)

#K=4:
tbl4 = read.table("pca_55_all.K4.a0.qopt")
tbl4_df <- data.frame(tbl4[,]) 
#View(tbl4_df)
tbl4_df$pop = pop[,3] #add pop lable to column
tbl4_df$ind = pop[,2] #add individual label to column
tbl4_df$srt = pop[,5]
#View(tbl4_df)
ord4 <- tbl4_df[order(tbl4_df$srt),]
#View(ord4)
bp4 = barplot(t(as.matrix(ord4[,1:4])), space = c(0.2), col = rainbow(4), border = NA, names.arg = ord4$pop,las=2, font.axis=1, cex.names = 0.6)

#K=5:
tbl5 = read.table("pca_55_all.K5.a0.qopt")
tbl5_df <- data.frame(tbl5[,]) 
View(tbl5_df)
tbl5_df$pop = pop[,3] #add pop lable to column
tbl5_df$ind = pop[,2] #add individual label to column
tbl5_df$srt = pop[,5]
#View(tbl5_df)
ord5 <- tbl5_df[order(tbl5_df$srt),]
#View(ord5)
bp5 = barplot(t(as.matrix(ord5[,1:5])), space = c(0.2), col = rainbow(5), border = NA, names.arg = ord5$pop,las=2, font.axis=1, cex.names = 0.6)

#K=6:
tbl6 = read.table("pca_55_all.K6.a0.qopt")
tbl6_df <- data.frame(tbl6[,]) 
View(tbl6_df)
tbl6_df$pop = pop[,3] #add pop lable to column
tbl6_df$ind = pop[,2] #add individual label to column
tbl6_df$srt = pop[,5]
#View(tbl6_df)
ord6 <- tbl6_df[order(tbl6_df$srt),]
#View(ord6)
bp6 = barplot(t(as.matrix(ord6[,1:6])), space = c(0.2), col = rainbow(6), border = NA, names.arg = ord6$pop,las=2, font.axis=1, cex.names = 0.6)

