# Rx_identifier for ancient elephants
# This script was modified from Mittnik et al 2016 (original author: Chuan-Chao Wang)
# Based on the ratio of X chromosome-derived shotgun sequencing data to the autosomal coverage to establish the probability of an XX or XY karyotype for ancient samples.
# example:samtools view -q 30 -b sampleID.fa.merged.bam > SampleID_q30.bam
        # samtools index SampleID_q30.bam
        # samtools idxstats SampleID_q30.bam > SampleID.idxstats
        # Rscript Rx_identifier.r SampleID > SampleID.Rx
		

args=(commandArgs(TRUE))
PREFIX=as.character(args[1])

idxstats<-read.table(paste(PREFIX,'.idxstats',sep=''),header=F,nrows=1000,row.names=1) #Armstrong lion genome has 19 resolved chromosomes (18 autosomal and 1 X) - modify this value according to the chr number of the taxon of interest
c1 <- c(as.numeric(idxstats[,1])) #this column in the idxstats file indicates the number of reference reads
c2 <- c(as.numeric(idxstats[,2])) #this column in the idxstats file indicates the number of mapped reads
total_ref <- sum(c1)
total_map <- sum(c2)
  
LM <- lm(c1~c2)
summary(LM)  

#calculate the normalized ratio of each chr:
Rt1 <- (idxstats[177,2]/total_map)/(idxstats[177,1]/total_ref) #chr name A1 
Rt2 <- (idxstats[4,2]/total_map)/(idxstats[4,1]/total_ref) #chr name A1
Rt3 <- (idxstats[221,2]/total_map)/(idxstats[221,1]/total_ref) #chr name A2
Rt4 <- (idxstats[223,2]/total_map)/(idxstats[233,1]/total_ref) #chr name A3
Rt5 <- (idxstats[124,2]/total_map)/(idxstats[124,1]/total_ref) #chr name B1
Rt6 <- (idxstats[267,2]/total_map)/(idxstats[267,1]/total_ref) #chr name BB2
Rt7 <- (idxstats[148,2]/total_map)/(idxstats[148,1]/total_ref) #chr name B3
Rt8 <- (idxstats[261,2]/total_map)/(idxstats[261,1]/total_ref) #chr name B4
Rt9 <- (idxstats[13,2]/total_map)/(idxstats[13,1]/total_ref) #chr name C1
Rt10 <- (idxstats[141,2]/total_map)/(idxstats[141,1]/total_ref) #chr name C2
Rt11 <- (idxstats[99,2]/total_map)/(idxstats[99,1]/total_ref) #chr name D1
Rt12 <- (idxstats[147,2]/total_map)/(idxstats[147,1]/total_ref) #chr name D2
Rt13 <- (idxstats[142,2]/total_map)/(idxstats[142,1]/total_ref) #chr name D3
Rt14 <- (idxstats[380,2]/total_map)/(idxstats[380,1]/total_ref) #chr name D4
Rt15 <- (idxstats[898,2]/total_map)/(idxstats[898,1]/total_ref) #chr name E1
Rt16 <- (idxstats[308,2]/total_map)/(idxstats[308,1]/total_ref) #chr name E2
Rt17 <- (idxstats[176,2]/total_map)/(idxstats[176,1]/total_ref) #chr name E3
Rt18 <- (idxstats[150,2]/total_map)/(idxstats[150,1]/total_ref) #chr name F1
Rt19 <- (idxstats[135,2]/total_map)/(idxstats[135,1]/total_ref) #chr name F2
Rt20 <- (idxstats[391,2]/total_map)/(idxstats[391,1]/total_ref) #chr name X

#calculate averaged normalized ratio of the X chr
tot <- c(Rt20/Rt1,Rt20/Rt2,Rt20/Rt3,Rt20/Rt4,Rt20/Rt5,Rt20/Rt6,Rt20/Rt7,Rt20/Rt8,Rt20/Rt9,Rt20/Rt10,Rt20/Rt11,Rt20/Rt12,Rt20/Rt13,Rt20/Rt14,Rt20/Rt15,Rt20/Rt16,Rt20/Rt17,Rt20/Rt18,Rt20/Rt19)
Rx <- mean(tot)
cat("Rx :",Rx,"\n")
confinterval <- 1.96*(sd(tot)/sqrt(18))
CI1 <- Rx-confinterval
CI2 <- Rx+confinterval
cat("95% CI :",CI1, CI2,"\n")

if (CI1 > 0.8) {print ("Sex assignment:The sample should be assigned as Female")
} else if (CI2 < 0.6) {print ("Sex assignment:The sample should be assigned as Male")
} else if (CI1 > 0.6 & CI2 > 0.8) {print ("Sex assignment:The sample is consistent with XX but not XY")
} else if (CI1 < 0.6 & CI2 < 0.8) {print ("Sex assignment:The sample is consistent with XY but not XX")
} else print ("Sex assignment:The sample could not be assigned")

print ("***It is important to realize that the assignment is invalid, if there is no correlation between the number of reference reads and that of the mapped reads***")


