#setwd("~/Box Sync/aDNA/lions (deflami2@illinois.edu)/Cape_lions/2021_2022/Data&analyses/gen_distribution")

library(ggplot2)
library(readr)
library(patchwork)

#read per choromosome site files from ANGSD (filtered to -doCounts 1 -minQ 30 -minInd 28 -dumpCounts 2 -doDepth 1)
cov_177 <- read_tsv("ScFovFW_177.tsv", col_names = FALSE)
cov_4 <- read_tsv("ScFovFW_4.tsv", col_names = FALSE)
cov_221 <- read_tsv("ScFovFW_221.tsv", col_names = FALSE) 
cov_223 <- read_tsv("ScFovFW_223.tsv", col_names = FALSE)
cov_124 <- read_tsv("ScFovFW_124.tsv", col_names = FALSE)
cov_267 <- read_tsv("ScFovFW_267.tsv", col_names = FALSE)
cov_148 <- read_tsv("ScFovFW_148.tsv", col_names = FALSE)
cov_261 <- read_tsv("ScFovFW_261.tsv", col_names = FALSE)
cov_13 <- read_tsv("ScFovFW_13.tsv", col_names = FALSE)
cov_141 <- read_tsv("ScFovFW_141.tsv", col_names = FALSE)
cov_99 <- read_tsv("ScFovFW_99.tsv", col_names = FALSE)
cov_147 <- read_tsv("ScFovFW_147.tsv", col_names = FALSE)
cov_142 <- read_tsv("ScFovFW_142.tsv", col_names = FALSE)
cov_380 <- read_tsv("ScFovFW_380.tsv", col_names = FALSE)
cov_898 <- read_tsv("ScFovFW_898.tsv", col_names = FALSE)
cov_308 <- read_tsv("ScFovFW_308.tsv", col_names = FALSE)
cov_176 <- read_tsv("ScFovFW_176.tsv", col_names = FALSE)
cov_150 <- read_tsv("ScFovFW_150.tsv", col_names = FALSE)
cov_135 <- read_tsv("ScFovFW_135.tsv", col_names = FALSE)
cov_391 <- read_tsv("ScFovFW_391.tsv", col_names = FALSE)

#calculate historgrams per chromosome

his_177 <- ggplot(cov_177, aes(X2))+
    geom_histogram(binwidth=100000, fill="white", color="black") +
    labs(x="Position on chromosome", y = "Count per 100kb bin")+
    theme_classic()

his_4 <- ggplot(cov_4, aes(X2))+
  geom_histogram(binwidth=100000, fill="white", color="black") +
  labs(x="Position on chromosome", y = "Count per 100kb bin")+
  theme_classic()

his_221 <- ggplot(cov_221, aes(X2))+
  geom_histogram(binwidth=100000, fill="white", color="black") +
  labs(x="Position on chromosome", y = "Count per 100kb bin")+
  theme_classic()

his_223 <- ggplot(cov_223, aes(X2))+
  geom_histogram(binwidth=100000, fill="white", color="black") +
  labs(x="Position on chromosome", y = "Count per 100kb bin")+
  theme_classic()

his_124 <- ggplot(cov_124, aes(X2))+
  geom_histogram(binwidth=100000, fill="white", color="black") +
  labs(x="Position on chromosome", y = "Count per 100kb bin")+
  theme_classic()

his_267 <- ggplot(cov_267, aes(X2))+
  geom_histogram(binwidth=100000, fill="white", color="black") +
  labs(x="Position on chromosome", y = "Count per 100kb bin")+
  theme_classic()

his_148 <- ggplot(cov_148, aes(X2))+
  geom_histogram(binwidth=100000, fill="white", color="black") +
  labs(x="Position on chromosome", y = "Count per 100kb bin")+
  theme_classic()

his_261 <- ggplot(cov_261, aes(X2))+
  geom_histogram(binwidth=100000, fill="white", color="black") +
  labs(x="Position on chromosome", y = "Count per 100kb bin")+
  theme_classic()

his_13 <- ggplot(cov_13, aes(X2))+
  geom_histogram(binwidth=100000, fill="white", color="black") +
  labs(x="Position on chromosome", y = "Count per 100kb bin")+
  theme_classic()

his_141 <- ggplot(cov_141, aes(X2))+
  geom_histogram(binwidth=100000, fill="white", color="black") +
  labs(x="Position on chromosome", y = "Count per 100kb bin")+
  theme_classic()

his_99 <- ggplot(cov_99, aes(X2))+
  geom_histogram(binwidth=100000, fill="white", color="black") +
  labs(x="Position on chromosome", y = "Count per 100kb bin")+
  theme_classic()

his_147 <- ggplot(cov_147, aes(X2))+
  geom_histogram(binwidth=100000, fill="white", color="black") +
  labs(x="Position on chromosome", y = "Count per 100kb bin")+
  theme_classic()

his_142 <- ggplot(cov_142, aes(X2))+
  geom_histogram(binwidth=100000, fill="white", color="black") +
  labs(x="Position on chromosome", y = "Count per 100kb bin")+
  theme_classic()

his_380 <- ggplot(cov_380, aes(X2))+
  geom_histogram(binwidth=100000, fill="white", color="black") +
  labs(x="Position on chromosome", y = "Count per 100kb bin")+
  theme_classic()

his_898 <- ggplot(cov_898, aes(X2))+
  geom_histogram(binwidth=100000, fill="white", color="black") +
  labs(x="Position on chromosome", y = "Count per 100kb bin")+
  theme_classic()

his_308 <- ggplot(cov_308, aes(X2))+
  geom_histogram(binwidth=100000, fill="white", color="black") +
  labs(x="Position on chromosome", y = "Count per 100kb bin")+
  theme_classic()

his_176 <- ggplot(cov_176, aes(X2))+
  geom_histogram(binwidth=100000, fill="white", color="black") +
  labs(x="Position on chromosome", y = "Count per 100kb bin")+
  theme_classic()

his_150 <- ggplot(cov_150, aes(X2))+
  geom_histogram(binwidth=100000, fill="white", color="black") +
  labs(x="Position on chromosome", y = "Count per 100kb bin")+
  theme_classic()

his_135 <- ggplot(cov_135, aes(X2))+
  geom_histogram(binwidth=100000, fill="white", color="black") +
  labs(x="Position on chromosome", y = "Count per 100kb bin")+
  theme_classic()

his_391 <- ggplot(cov_391, aes(X2))+
  geom_histogram(binwidth=100000, fill="white", color="black") +
  labs(x="Position on chromosome", y = "Count per 100kb bin")+
  theme_classic()

pdf(file = "test.pdf"

(his_177 + his_4) / (his_221 + his_223) / (his_124 + his_267)/(his_148 + his_261)

dev.off()
#ggsave ('test.png', p)

#names of chromosomes (associated with domestic cat)
#ScFovFW_177
#ScFovFW_4
#ScFovFW_221
#ScFovFW_223
#ScFovFW_124
#ScFovFW_267
#ScFovFW_148
#ScFovFW_261
#ScFovFW_13
#ScFovFW_141
#ScFovFW_99
#ScFovFW_147
#ScFovFW_142
#ScFovFW_380
#ScFovFW_898
#ScFovFW_308
#ScFovFW_176
#ScFovFW_150
#ScFovFW_135
#ScFovFW_391


