#make principal component eigenvalue plot as seen in HCHS paper, fig. 1
#by Angela Andaleon (aandaleon@luc.edu)
library(data.table)
library(dplyr)
library(ggplot2)
library(reshape2)
pcs <- fread("../kingpc.ped")
region <- fread("../region.txt", header = T)
region <- region %>% dplyr::select("IID", "region")
pcs <- left_join(region, pcs, c("IID" = "V2"))
pcs$V1 <- NULL
pcs$V3 <- NULL
pcs$V4 <- NULL
pcs$V5 <- NULL
pcs$V6 <- NULL
pcs$FID <- NULL #Columns leftover from PLINK FID, irrelevant rn

colnames(pcs) <- c("IID", "REGION", "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10", "PC11", "PC12", "PC13", "PC14", "PC15", "PC16", "PC17", "PC18", "PC19", "PC20")
pcs <- pcs[complete.cases(pcs),]
pcs <- pcs %>% dplyr::select(IID, REGION, PC1, PC2, PC3, PC4, PC5)
pcs_t <- data.table::transpose(pcs) #why did this suddenly break
pcs_t <- pcs_t[-c(1:2),]
IDs <- pcs$IID
PCs <- as.data.frame(c("PC1", "PC2", "PC3", "PC4", "PC5"))

colnames(PCs) <- "PCs"
pcs_t <- cbind(PCs, pcs_t)
colnames(pcs_t) <- c("PCs", IDs)
melted <- melt(pcs_t, id.vars = "PCs")
melted$PCs <- factor(melted$PCs, levels = c("PC1", "PC2", "PC3", "PC4", "PC5"))
region <- pcs %>% dplyr::select(IID, REGION)
melted <- left_join(melted, region, by = c("variable" = "IID"))

#rotation plot
fig1 <- ggplot(data = melted, aes(x = PCs, y = as.numeric(value), group = variable)) + 
  geom_line(aes(color = REGION)) + 
  scale_color_brewer(palette = "Set1") + 
  labs(color = "Region", y = "eigenvalue") + 
  theme_bw() + 
  theme(text = element_text(size = 15))

pdf("rotation_plot.pdf", width = 6, height = 3)
print(fig1)
dev.off()